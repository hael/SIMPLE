# Learned pose and state assignment: a student distilled from the probabilistic search

Proposal of 2026-11-05, revised on 2026-11-06 with Hans and after a methodology review. SIMPLE's
matching machinery supplies per-dataset soft labels over 3D direction and state, and a small model
learns to propose candidate poses for particles streaming off the microscope. Whether this saves
time against the actual stream route is an acceptance measurement. Nothing is built. Sections 1 and
2 are the idea and what exists, sections 3 to 9 the design, section 10 the validation, section 11
the order of work, section 12 the open decisions for Hans, sections 13 to 15 the risks, the rules
and the sources.

## 1. What is proposed

A pose student: a small network, trained per dataset on the probabilistic posterior of the dense
3D table over a fixed grid of (state, direction) classes. It reads a particle's CTF-multiplied
polar Fourier transform in a low-resolution band and proposes a distribution over the grid.
Circular convolutions and pooling provide approximate in-plane rotation invariance; shifts are
learned as tolerance from the calculator's own phase ramps and found exactly by the shift map in
the candidate search that follows. Two students, one per half-set, each assign the new
particles of their own half, against their own half's references and spectrum. Their designs
and gates are selected independently within that half, with micrograph-level development,
calibration and untouched audit cohorts. A bounded candidate search recovers the full rotation and
shift; an empirically calibrated confidence measures whether that seed has the right state and lies
inside the local refinement's capture range. A row whose candidate coverage or confidence is
insufficient takes the existing route. The stream baseline is sub-class assignment followed by
continuous refinement, so an end-to-end speedup must include training, fallback and refinement. No
new command-line key: the model file travels through `infile` and `fname` as the class-average
quality model's does, and every threshold and constant that decided lives in the model file.

What sets this apart from the published attempts, which learn in real space, is the polar
implementation: the representation the search already uses, with the in-plane rotation a cyclic
shift of the input, the translation handled exactly by the Fourier-Bessel shift map, the CTF
applied as the search applies it, the sigma2 whitening of the objective, and the table's own
posterior as the teacher. Everything the student sees and learns from is what the probabilistic
search computes.

The design's constraints: deterministic teacher scores, preserved target uncertainty, same-half
production with declared upstream provenance, calibrated candidate coverage and full-seed
confidence, and measured complete-stage savings. Sections 4, 7, 8 and 10 define the implementation
and acceptance requirements.

## 2. Where the idea stands

### 2.1 Prior work

- Unsupervised amortized pose inference (cryoPoseNet, cryoAI, cryoFIRE): an encoder learns the poses
  jointly with the map, with no labels. Encoders converge slowly and get stuck; the semi-amortized
  method of Shekarforoush et al. (NeurIPS 2024) and DRGN-AI both retreat to per-particle optimization
  after an amortized warm start. That is evidence against learning poses from scratch, not against the
  supervised setting.
- Supervised pose assignment, which is this setting: the CESPED benchmark (Sanchez-Garcia et al.,
  arXiv 2311.06194; the package is checked out at `~/src/cesped`, section 2.4) and cryoPARES
  (Sanchez-Garcia et al., bioRxiv 10.1101/2025.03.04.641536). cryoPARES trains on a refined apo
  dataset and infers on ligand-bound datasets of the same protein; a U-Net feeds an Image2Sphere head
  that outputs a distribution over an SO(3) grid of 294,912 rotations; the CTF enters as a
  phase-flipped input channel; in-plane rotation and shift are learned by augmentation. Its numbers:
  network inference above 30,000 particles per minute per GPU; the raw top-1 pose within 5° of the
  RELION pose for about 30% of the particles, above 60% after a local refinement; resolutions 0.1 to
  0.6 Å worse than RELION on six datasets; a direction-normalised confidence rejects junk as well as 2D
  classification does; two RELION runs disagree by more than 5° on over 20% of the particles. Trained
  per protein, up to 100 epochs on GPUs, checkpoints of 40 GB.
- CESPED's own numbers: eight EMPIAR entries of at least 200,000 particles at 1.5 Å per pixel with
  poses re-estimated by RELION auto-refine, each split at random into two halves, a model trained on
  one half and judged on the other; its Image2Sphere baseline (an ImageNet ResNet152 on 256-pixel
  real-space images, an orthographic projection onto a hemisphere, S2 and SO(3) convolutions at
  lmax 8, a softmax over 294,912 rotations, phase-flipped input, no shift prediction) reaches mean
  angular errors of 13° to 45° per entry (28.7° over the benchmark) and maps from its poses 1.3 Å
  worse than the truth at FSC 0.143 and 3.4 Å at 0.5; a local refinement within 6° of its poses
  brings one entry's map correlation from 0.969 to 0.997 at a threefold speed-up over auto-refine.
  The same entry simulated gives 6.0° against 13.3° real, and consensus labels 8.3°: the labels and
  the noise set the error as much as the model.
- Lessons taken: the raw pose is coarse and a local refinement supplies the precision; the confidence
  is worth as much as the pose; the student is judged on soft targets and downstream recovery; the
  model is dataset-specific, so it is trained per dataset and retrained when the references change.

### 2.2 What SIMPLE has

- The dense 3D table (`src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90`): its grid is
  `build%eulspace` at `params%nspace` times the states; `fill_tab_range` fills a particle's row over
  every (state, direction) reference, drawing an in-plane rotation per reference from the row
  posterior at the particle's temperature (`score_ref` calls the calculator's `gen_posterior_val`,
  which requires `objfun=euclid` and takes a keyed variate) and polishing a selected subset in shift
  and rotation through the shift map. `posterior_row_stats` gives the normalised row's mode, `pmax`,
  `neff` and `k99`; `top_posterior_extras` supplies the weighted top-k. Skipping `assign_batch`
  therefore does not make the fill a deterministic teacher. The deterministic argmin over the
  in-plane grid at a given shift exists beside it, as the diagnostic posterior's path
  (`deterministic_best` over `gen_best_objfun_val` and `gen_best_euclid_val_cross`) under euclid. The registration pass of solve3D_addon (`set_refine3D_register_policy` in
  `src/main/solve/simple_solve3D_controller.f90`) is one full-cohort `prob` iteration at a fixed
  `NSPACE_REGISTER` (1000) against frozen references: a table fill at a chosen grid, as a step.
- The polar calculator (`src/main/pftc/simple_polarft_calc.f90`): the particle sections over the
  Friedel half of the circle (`pfts_ptcls`), the references per parity (`pfts_refs(:,:,:,parity)`:
  a particle is matched against its own half's reference, the gold standard), the memoised
  full-ring angular spectra of the CTF-multiplied particle (`ft_ptcl_ctf`) and of the references
  (`ft_ref`), the CTF matrices (`ctfmats`), the sigma2 shell weights as two global spectra, one per
  half (`sigma2_noise(shell, half)`, `memoize_sigma2_weights`), the shifted particle
  (`memoize_shifted_ptcl`), and the correlation over all in-plane rotations with the shift search per
  reference.
- The polar shift map (`src/main/pftc/simple_polar_shiftmap.f90`, the calculator's
  `simple_polarft_shiftmap.f90`; `doc/implementation_notes/completed/polar_shift_search.md`): the
  cross term of a reference/particle pair over every rotation and a disk of shifts from the
  Fourier-Bessel expansion, exact in the shift direction, gridded in magnitude (`gen_shiftmap`,
  `shiftmap_mode`, `shiftmap_shift`). The dynamic centring
  (`polar_shift_search_dynamic_centring.md`) keeps per particle the last move `shincr` and the width
  `shsig`, numbered slots of the ori module.
- The stream's route into 3D (`doc/policies/stream/stream_3D_ingestion_policy.md` section 10;
  `src/main/stream/stages/simple_stream_stage_solve3D.f90`, `src/main/solve/simple_state_assignment.f90`):
  no probabilistic search on particles; sub-class averages matched against reprojections of the
  merged state references (`LP_ASSIGN` 15 Å, `NSPACE_ASSIGN` 1000 directions), the particles inherit
  the sub-class state and pose, then one iteration of `refine3D refine=cont` from the
  seeds `set_refine_seeds` sets (the halves of the partition, one per micrograph by `eo_transaction`;
  the nearest grid direction), with the halves' starting volumes reconstructed by refine3D from the
  mapped poses. `refine=cont` keeps the state fixed and bounds the rotation by `athres_cont` (10° on
  this route) and the shifts by `trs`. Its recoverable full pose, including the state, is the bar;
  viewing-direction agreement on a 6° grid alone is insufficient.
- The logistic core (`src/main/quality/simple_quality_logistic.f90`), binary and at most sixteen
  features, is not the student's shape but is the right tool for the confidence calibration of
  section 8; the learn/evaluate/promote lifecycle and the versioned key=value model file of
  `src/main/cavg_quality/simple_cavg_quality_model.f90` are the pattern.
- OpenBLAS is linked (`cmake/Dependencies.cmake`) and `sgemm` is called from Fortran
  (`src/utils/math/simple_linalg.f90`, flex); the OpenMP BLAS policy keeps the library on the calling
  thread inside a team (`doc/implementation_notes/completed/blas_openmp_policy.md`). Nothing in the
  repo does backpropagation: the network is new code.
- Private programs: `add_private_prg` in `src/main/ui/simple_private_prgs.f90`, the `case` in
  `production/simple_private_exec_driver.f90`, a thin commander (`assign_cavgs_states` in
  `src/main/commanders/simple/simple_commanders_cavgs.f90` is the model), started by the stage as a
  job with the 3D jobs' resources.
- Validation data (section 10): the Ribosembly three-state benchmark on dell2
  (`/data2/multistate_bench`: 80,000 simulated particles, box 128 at 3 Å per pixel, SNR 0.01, true
  poses, shifts uniform to plus or minus 10 pixels, CTFs drawn per particle, labels at 30/60/10
  percent, half-sets, truth maps with their FSC at 12 Å, a 20,000 subset); `/data/simulated_vcrnapol`
  on the Dell (30,000 particles, two states, SNR 0.1, truth); the real multi-state `/data2/trpv1_c2`
  (its `6_refine3D_states`) and `/data/bgal`; the small membrane protein `/data/exp_gate_FlhB` and
  streptavidin `/data/10335-streptav`; two CESPED entries from Zenodo (section 2.4);
  `/data2/stream_holjunk` for the stream route with junk; `multistate_workflows`
  (`src/main/commanders/test/simple_commanders_test_highlevel_msworkflows.f90`) as the pattern for a
  nightly route on simulated beta-galactosidase states.

### 2.3 What XEMA's two-fold method contributes

XEMA (`~/src/XEMA`, Hans's cross-linking mass-spectrometry engine) motivates label-free partitions,
held-out choices, published agreement and replayable model records (the sources in section 15).
Those are useful principles, but its cross-scoring arrangement is not a cryo-EM gold-standard
guarantee. Here each production half owns its teacher references, spectrum, student, model selection and
confidence gate. Inner micrograph splits serve model development; a separate cohort
calibrates the gate and an untouched cohort audits it. The two halves may choose different rungs.
Their agreement is diagnostic and never selects a common winner or gate. Counts are compared as
paired micrograph-clustered outcomes, not independent Poisson observations. One route supplies each
row, fixed-chunk reductions make a fit reproducible across thread counts, and the model records
every data-dependent choice and all constants. Sections 7 and 8 define these boundaries.

### 2.4 What the CESPED checkout teaches

`~/src/cesped` (oxpig/cesped) holds the benchmark's loader, the Image2Sphere baseline, its training
configuration and its evaluation. It settles several things.

- **Held-out halves, and an evaluation through reconstruction.** An entry's half 0 trains (with a
  70/30 split inside it for validation and hyperparameters) and half 1 is judged, then the reverse.
  The evaluation reconstructs from the predicted poses and reports the FSC against the truth map
  beside the truth halves' own FSC, so that the method's loss is separated from the map's quality.
  That cross-half judgement is an offline benchmark here (section 8), and tier 4 adopts the same two
  numbers on real data (section 10).
- **Means hide the tail.** CESPED reports mean angular errors, which a bimodal distribution (most
  particles within a few degrees, a fraction wildly wrong) turns into tens of degrees; cryoPARES
  later reported the fraction within 5°. The reads here are fractions of recoverable seeds and the
  size of the set that holds the posterior, never a mean angle.
- **Seed, then refine.** The only result of theirs that reaches the truth's quality is the local
  refinement within 6° of the predicted pose, at a threefold speed-up; the student alone is not
  the product. That is this design's consumer: a seed, a search at the predicted directions, the
  local refinement.
- **Labels and confidence.** Their label is the nearest grid rotation to one RELION pose, weighted
  in the loss by RELION's `rlnMaxValueProbDistribution`; a cleaner label on the same particles cut
  the error by a third. The teacher's soft score row at the student's band (section 4) is the same
  idea carried through: the whole row, not one pose and one weight.
- **Translations were the stated gap.** Image2Sphere predicts no shift, the evaluation used the
  truth shifts, and joint rotation and translation estimation heads the authors' list of
  improvements. The shift map in the candidate search at the predicted directions (sections 5 and 8) is this
  design's answer.
- **CTF.** Their loader offers `ctf_multiply` beside `phase_flip`, and the paper used the flip; no
  one measured the difference. The CTF-multiplied input of section 3 is a declared ablation read
  on the Ribosembly bench (tier 2): the same rung on the multiplied and on the flipped input.
- **Capacity.** Their encoder is an ImageNet-pretrained ResNet152 on 65,000 real-space pixels,
  trained for up to 400 epochs on four A100s; this design's network has about a million weights
  on 10,000 band-limited polar values with the rotation built in and the shift found by the map.
  The lower capacity is plausible on the lower-dimensional input and is a risk the ladder's wider
  rungs bound (section 6).
- **The head.** Their distribution over SO(3) is band-limited at lmax 8, so its peaks cannot be
  much sharper than 20° whatever the 3.75° readout grid; and it classifies the in-plane angle as a
  class dimension, which the authors name as not separable from the in-plane translation. The flat
  softmax over directions with the in-plane angle left to the search is the cheaper factoring.
- **The data.** Eight public entries with RELION poses, two membrane proteins in a micelle or a
  nanodisc among them (the GPCR complexes 10786 and 11120 are the hardest, at 33° and 45°), two
  symmetric ones (C2, D2), and published baselines. Tier 4 takes two CESPED entries as its real
  data beside the Dell's own, so that the numbers compare with Image2Sphere's and cryoPARES's on the
  same particles. No entry is on this box yet; they come from Zenodo, with a small TEST entry for
  smoke tests.

## 3. The student's input

`simple_polar_feats` (`src/main/learn/`): one routine builds the input of a particle from the polar
calculator, called by training and prediction alike so that the two never drift apart.

1. **The band.** Rings `kfromto(1)` to `box*smpd/lp` with `lp` the student's low-pass. The band is a
   rung of the ladder (section 6), 15 Å first, then 10 and 8 Å; the number of rings depends on
   `box*smpd/lp` alone, about 25 at 15 Å for the benchmark's box 128 at 3 Å per pixel.
2. **The CTF-multiplied particle.** The tensor is built in the angular-sample domain from
   `pfts_ptcls(:,k,i)*ctfmats(:,k,i)`, with the CTF itself as a third channel. The calculator's
   `ft_ptcl_ctf` is the angular FFT of this Friedel-extended ring: its index is the harmonic order,
   not the angle, so it is the shift map's input (section 5) and never an angle-axis CNN's. The product
   follows the search's cross term `Q_k = REF CTF conj(X)`. With `X = CTF P + N` it has the signal
   `CTF^2 P` and the noise `CTF N`: a shell near a zero carries no signal and enters at no weight,
   where a phase flip would feed its noise at full amplitude. The multiplication is invertible away
   from the zeros, so it changes the conditioning of what the first layer has to learn, not the
   information; it does not by itself prove an optimal representation, and the phase-flipped input is
   the ablation of tier 2. The CTF channel supplies the varying envelope and the astigmatism; whether
   the model uses it correctly is tested by defocus and acquisition holdouts. Nothing is restored:
   no phase flip, Wiener filter or deconvolution in any rung; the output is a class, not an image.
3. **The shell weights.** The shells are weighted as the search weights them under objfun=euclid,
   the objective of the biomolecular workflows this is for: `k/sigma2(k)`, where sigma2 is the half's
   global spectrum, never a per-particle one: the calculator's weights are two spectra, one per half, and a part's per-particle residual
   records exist only to be summed into them. An unseen particle therefore needs nothing beyond its
   half, assigned at import, and its own weighted power (`wsqsum_ptcl_band`), which scales the band.
   The weights belong to the model, not to the project: each student's record holds the weights it was
   trained with, and prediction applies the weights of the student it applies, never the project's
   current sigma2, so a later refinement that moves the spectrum does not move the student's input
   distribution and a retrain picks the new spectrum up. A project without a spectrum is given one by the
   sigma2 bootstrap (`calc_pspec`) before training or prediction, never scored under another
   objective; the model records the spectrum it was trained with. Tensor
   amplitudes use `sqrt(w_k)` so that a squared norm applies the shell weight once. The logarithm of the weighted power is one scalar side input whose usefulness for
   the confidence is measured; it does not authorise a particle-rejection decision.
4. **The full circle.** The Friedel half (`pftsz` angles) extended by its conjugate to `nrots` angles,
   then angularly low-pass filtered and resampled to `NANG` (128). Circular convolution is
   equivariant to integer cyclic shifts on this grid; arbitrary physical in-plane angles require
   interpolation and are approximate. The resampling and filter convention is recorded, and the error
   is tested against angular increments finer than the direction grid, with astigmatic CTFs.
5. **The tensor.** Three channels (the real and imaginary parts of the CTF-multiplied particle, the
   CTF) by the rings by `NANG`, about 9,600 values at 15 Å for the benchmark's box. The phases are
   kept throughout; nothing is reduced to a power, which at 15 Å would see a particle's size and mass
   and separate no conformational states of one particle. Augmentation lives in the same module: the
   calculator's own phase ramp for shifts (`memoize_shifted_ptcl`), cyclic shifts for integer
   angular-bin rotations, the declared interpolation for other angles, with the CTF channel rotated
   consistently with the observation.

## 4. The class grid and the teacher

- **The grid.** `oris%spiral` at `NSPACE_STUDENT` (1000) directions, the class-average matching's grid,
  about 6° apart, over the whole sphere in c1 and the asymmetric unit as `build%eulspace` is built for
  a symmetric point group, times the states: `nstates*NSPACE_STUDENT` classes, 3,000 for three states.
  The in-plane angle is not a class: the student is invariant to it by construction and the search at
  the predicted directions recovers it (section 8). The student's output is a distribution over
  (state, direction) modulo the point group; no rotation or shift is regressed, the candidate search
  supplies the in-plane angle and the shift, and Euler angles stay the project's interchange.
- **The references are the half's own.** The teacher of half A scores half A's particles against the
  per-state reconstructions of half A, as the matcher's parity references are prepared
  (`prepare_refs_sigmas_and_pftc`; refine3D reconstructs them from the project's halves at its first
  iteration). The stream's merged state volumes serve the sub-class matching and are never a
  teacher's reference. The model records the references' provenance (the project, the rows, the
  band, a hash of each volume).
- **The teacher is the table's tempered posterior, made deterministic, at the student's grid and
  band.** A new scoring entry beside the table fill, reusing the matcher's preparation, the calculator's objective grids and the shift
  map's kernels, with no keyed draw anywhere. For every class c the deterministic best over the
  complete in-plane grid at the particle's prepared shift (`gen_best_objfun_val`), then, for the
  `NTEACH_SHIFT` (64) classes of lowest loss, the shift map over the declared disk and the polish,
  `L_c = min_(rotation, shift) loss_c`; the other classes keep their unshifted loss and the row
  records which classes were polished. The rule is the same for every particle and every class, the
  budget declared, the disk, prior, band, objective, discretisation and polish policy recorded. The existing
  sampled rows do not meet this contract; the diagnostic posterior's deterministic path
  (`deterministic_best`) is the step to extend.
- **The soft-label definition.** `t_c = exp[-tau_i (L_c-L_min)] / sum_j exp[-tau_i (L_j-L_min)]`.
  These are normalised profile scores, not a likelihood marginalised over nuisance poses. `tau_i` is the recorded `prob_beta` times the particle's Gaussian-noise temperature
  (`posterior_temperature`), the factor the table draws at. Objective, temperature and frozen shell weights are explicit
  in the model. Integrating likelihoods over rotations and shifts would be a different teacher,
  requiring the prior and the cell measure.
- **What the earlier search contributes.** The halves, the per-state half reconstructions, which
  particles have a converged pose (the candidate training rows), and the initial shift support. The
  teacher and the prediction use the same image preparation and initial shift support; converged
  shifts or centring widths that a new particle will not have cannot restrict the training searches.
- **The targets retain uncertainty.** The loss uses the full soft rows in bounded batches. There is
  no uniform label smoothing, no fixed top-k renormalisation and no inverse-frequency class weighting
  in the probability head: with 3,000 classes, 5% of uniform smoothing would cap the mass in 64
  classes at about 95% even for a perfect fit, so the 99% coverage of section 8 could never be met;
  and reweighting the classes changes the meaning of the predicted probabilities. A later sparse
  encoding must preserve the omitted mass and its distribution, or bound its approximation error; a
  broad row stays broad. A micrograph-stratified draw may improve coverage, but a nonuniform draw is
  corrected by its known inclusion weights for the intended particle population. Rare-view and
  per-state performance are reported separately rather than changing the probabilities' meaning.
- **Training and evaluation targets have different roles.** Training may use the band-specific soft
  scores. Every rung is evaluated against one frozen, same-half full-search reference pose (the
  converged run's, at its band) and the actual downstream refinement, with the truth where it
  exists. A coarser band does not win by agreeing with its own easier teacher.
- **The teacher's resolvability is read before anything is trained.** The deterministic rows give,
  per row and band, `pmax`, `neff` and `k99`; their distributions over the training particles say
  whether the band carries a concentrated assignment under that score model, before the ladder runs.
  These are conditional on the references, the objective and the temperature, not an
  information-theoretic bound. A band at which the median `k99` spans a ring of directions is
  reported as such, and the student's claim for that dataset is assessed against the common reference
  and the truth (section 10).
- **The consumer's grid is decoupled** from the student's by the local search: the student's
  directions seed an in-plane and shift search and then `refine=cont` within `athres_cont`, and
  refine3D's own grid (`set_projs` picks the nearest direction in the job's grid) can be 2,500 or
  5,000 as the stage likes.
- **Cost.** Per particle, the in-plane best for every class is one correlation over the rotations
  per reference from the memoised spectra, and the shift map runs for 64 classes, so the teacher is
  of the order of a billion floating-point operations per particle and band, on at most `NTRAIN_MAX`
  (200,000) particles drawn by whole micrographs, plus the common-reference evaluations. It is an
  estimate: the complete nuisance searches, the per-half reference preparation, the calibration and
  the audits are charged to training and measured in tier 5 before any time saving is claimed.


## 5. Shifts and the Fourier-Bessel shift map

On the half ring the shift matrix of s = rho (cos a, sin a) is exp(i z_k sin(theta_p - a)) with
z_k = 2 pi k rho / box, and the Jacobi-Anger expansion turns the cross term of a reference/particle pair
over every rotation and a whole disk of shifts into one Bessel-weighted sum over angular orders n:

    a_n(irot, rho) = (1/nrots) sum_m exp(-2 pi i m (irot-1)/nrots) sum_k w_k J_n(z_k(rho)) R^_k[m-n] conj(X^_k[m]),
    C(irot, a, rho) = Re sum_n exp(-i n a) a_n(irot, rho),

R^ and X^ the memoised full-ring spectra of the reference and of the CTF-multiplied particle. It is
exact, the direction a is exact, only rho is gridded, and the orders kept follow the truncation rule:
about 20 at a 30 Å radius in the 15 Å band, about 40 at 8 Å. Three things in it matter here.

- The memoised spectra are the Fourier-Bessel domain (ring k, order m). A rotation is the phase
  exp(i m alpha) per order; a shift is a convolution along m with the kernel J_n(z_k) exp(-i n a),
  whose width grows with k and rho. No basis is diagonal in both, so an architecture absorbs one
  nuisance exactly and must search or learn the other. The circular convolution along the angle
  absorbs the rotation; a convolution along m wide enough to absorb a 30 Å shift at the band's outer
  ring would span about 40 orders, which is the map's work, not a network's.
- The map is linear in the reference, so a reference may be a reprojection or a class average in the band.
- The CTF enters as the search applies it, in the cross spectrum, per particle and exact.

The shift is handled in two places, and in neither by the network's architecture.

- **Training: tolerance by augmentation.** The network of section 6 is trained with shifts drawn
  uniformly in the disk through the calculator's shift matrices (`memoize_shifted_ptcl`), the
  exact phase ramps, and with rotations as cyclic shifts: it learns to put a shifted particle in
  the right class. Its accuracy against the shift magnitude is a declared read (section 10).
- **Prediction: the map and the polish at the predicted directions.** For every candidate
  (state, direction) of section 8, the reprojection is searched over all in-plane rotations and the
  declared disk of shifts by the map (`gen_shiftmap`, `shiftmap_mode`) and the polish, from the
  particle's prepared shift: the pool's 2D shift in the stream, the centring record where one
  exists, zero otherwise. The search keeps the full declared support, with a disk retry or the
  fallback when the result hits a boundary.

**Only the accepted search pose becomes the centring record.** An accepted full search result
supplies the shift and `shincr` in the native-pixel frame. `shsig` stays at its unknown sentinel
unless its value comes from an existing objective-based uncertainty that passes coverage tests.

**Cost.** Per candidate direction the map costs (N+1) nrho (K nrots multiply-adds and one FFT) plus
the readout: at 15 Å with 20 orders, 8 magnitudes, 25 shells and about 400 rotations, 10 to 15
million floating-point operations, so a candidate set of a few directions costs a few times the
network and the set of 64 of section 8 about a billion operations per particle. These are kernel
estimates, not evidence of a speedup over p07's sub-class route; candidate counts, retries and
fallback are part of section 10's timing.

## 6. The ladder and the network

**The ladder.** The designs, simplest first: `cnn16`, `cnn32` and `cnn64` (16, 32 and 64 channels
in the first block, doubling per block), each at the bands 15, 10 and 8 Å, visited in cost order:
every design at 15 Å, then at 10 Å, then at 8 Å. Each half chooses its rung independently under section 7's recorded development
rule and common evaluation target. Measured costs may change the traversal order, but the order is
frozen before the selection run. A small membrane protein in a micelle is, at 15 Å, dominated by the
micelle, whose outline gives the tilt of the membrane normal and little of the azimuth about it; the
azimuth may sit at 8 to 12 Å, in the transmembrane helices, which is why the search's late stages run
there, and why the band is a rung.

**The network** (`simple_circ_cnn`, `src/main/learn/`): an encapsulated class with private components,
`new`/`kill`, forward, backward, one Adam step, read and write.

- On the tensor of section 3: three blocks of (a 3 by 3 convolution with circular padding along the
  angle and zero padding along the rings, bias, ReLU, average pooling along the rings only), with
  16, 32 and 64 channels; angular stride one throughout, then the mean over the angle for invariance
  to integer cyclic shifts (interpolation still limits the physical rotation invariance). Angular
  decimation is not used without an explicit error and cost comparison: three stride-two pools would
  leave a 22.5° angular step at `NANG` 128, coarser than the direction grid. Then a dense layer to
  256 with ReLU, the side scalar appended, and a dense layer to the classes with a softmax. About one
  million weights, most in the last layer; the per-particle cost is counted on the retained layout
  before it is quoted.
- The convolution as im2col plus `sgemm`, the backward pass as the two transposed products; weights,
  activations and gradients `real(sp)` and allocatable; one Adam state per weight array.
- **Bit-identical on any thread count.** The mini-batch is cut into fixed chunks; the chunks are
  computed in parallel, each into its own gradient buffer indexed by chunk, not by thread, with
  per-thread activation buffers sized by `omp_get_max_threads()` read before the region and the BLAS
  on the calling thread inside it; the chunks are reduced serially in chunk order after the region.
  The same sum at one thread and at any team, bit for bit, and the tester asks for identity, not a
  tolerance.
- The loss: the cross-entropy against the complete soft targets, with the sampling corrections of
  section 4; Adam at 1e-3, batches of 256, weight decay 1e-5; at most 30 epochs. Early stopping
  uses a fixed micrograph holdout inside the fitting cohort, distinct from the model-selection
  folds, the calibration and the audit. No refit absorbs calibration or audit rows without new
  holdouts.
- **The model file.** A key=value text header (`model_version=1`) and the weights of the two students
  as two binary files beside it (decision 7). The header records each half's design and band, the
  grid with `nspace`, `nstates` and `pgrp`, `NANG`, the frozen weights, the reference and spectrum
  hashes and source cohorts, the teacher's definition and temperatures, the training
  inclusion weights, every split and micrograph fingerprint, the calibration and shortlist rules,
  the uncertainty limits, the selection statistics, the timings, and each half's promotion or
  decline reason. The record distinguishes model scores, empirical coverage and calibrated seed
  confidence; a reader refuses incompatible grids, objectives, references or partition history.

## 7. Two halves: development, calibration, audit and decline

- **The halves** are the project's even-odd partition, one half per micrograph in the stream
  (`eo_transaction`). Each student is trained and applied within its own half. The record includes
  particle, state, direction and micrograph counts and stable stack/particle identities. Prediction
  verifies the recorded half of every overlapping training row with `ensure_field_halves`;
  additional rows are permitted when they follow the same partition rule, so a growing project need
  not have the identical row fingerprint of its training snapshot. Rows used for fitting, early
  stopping, design selection or calibration keep the incumbent route: the first delivery applies no
  out-of-sample calibration to an in-sample prediction. Batch reuse on those rows needs cross-fitted
  students inside the same reconstruction half.
- **The held-out roles.** Whole micrographs are assigned before any fit to development, calibration
  and audit cohorts. In the stream the audit cohort is the next publication's new rows, which no map,
  spectrum or fit has seen; in batch validation it is a held-out acquisition block. Development
  is split into two folds by a recorded hash bit, with early-stopping rows reserved inside each
  training fold. Calibration and audit particles must not have contributed to the teacher's maps
  or spectra: earlier frozen artifacts, or artifacts rebuilt without those cohorts; where that
  cannot be met the validation is reported as conditional and promotion waits for independent audit
  data. No calibration or audit statistic selects a rung or changes its weights, and no final refit
  absorbs these cohorts.
- **A common development objective.** Each rung predicts both held-out development folds. For every
  row: shortlist inclusion of the common-reference solution, the correct state, and full-pose
  recovery after the candidate search and the actual local refinement, plus runtime and fallback.
  Symmetry-aware full rotations and native-pixel shifts are used, with the capture criteria fixed
  before the comparison (section 8). Band-specific teacher-mode agreement is a diagnostic, not the
  selection target.
- **The selection rule.** Each half traverses the frozen cost order on its own. Candidates are
  compared on the same rows by paired differences in successful seeds, with whole micrographs
  resampled under fixed keys, never individual particles. Before the first benchmark the confidence
  level, the minimum useful improvement, the multiplicity adjustment for the finite ladder, the
  resampling count and the cost tie-break are fixed in the model policy (decision 16). The cheaper
  incumbent stays unless the adjusted interval supports the declared improvement; the incumbent route
  is an available outcome. The square root of a count is not a paired comparison's standard error.
  There is no joint winner across the reconstruction halves.
- **Calibration and audit.** After choosing and fitting one rung, the calibration cohort sets the
  distribution calibration, the shortlist coverage and the full-pose confidence gate. The model is
  frozen and evaluated once on the audit cohort. Promotion requires the predeclared seed-risk,
  coverage, per-state, reconstruction and timing criteria, with micrograph-clustered intervals. An
  audit failure declines that half's model; tuning after inspecting it requires a fresh audit cohort,
  and repeated retraining cannot reuse the same audit observations for an eventual pass.
- **Agreement is diagnostic.** An offline comparison may score untouched data with both students and
  report symmetry-aware direction and state agreement and the symmetric Kullback-Leibler divergence.
  This neither chooses a rung nor substitutes for the half-map FSC or the truth, and never feeds
  training, thresholds, production poses or retraining. The selection table
  (`<fname>_selection.txt`) records each half's cohort counts, common-reference scores, paired
  intervals, costs, choice and audit.
- **Decline by counts, before any fit.** Fewer than `MIN_HALF_PTCLS` (5,000) fitting particles or
  `MIN_HALF_STATE` (500) fitting modes in an active state declines that half; the development folds
  must satisfy the same floors; the calibration and audit cohorts need enough independent
  micrographs and per-state observations for the declared risk limits (decision 16). An unavailable
  comparison or an unsupported gate declines to the incumbent route; it never promotes the first
  rung by default. A populated half cannot certify the other half's model.

## 8. The two programs, the confidence and the arrangement

- **`learn_pose_student`**, a private program: the strategy routine `run_learn_pose_student(cline)` on
  the caller's `(params, build)` with a thin commander. Inputs: `projfile`, `nstates`, `lp`, `nspace`
  (the student grid), `mskdiam`, `pgrp`, `fname`, `nthr`; the references are the project's per-state
  half reconstructions, read where the project registers them or reconstructed as refine3D does at
  its first iteration (charged to training). Steps: prepare the parity references and the calculator
  as refine3D's matcher does; assign and record the cohorts; decline unsupported halves; build the
  deterministic score rows of section 4 in batches at each band, with their row statistics; build the
  features; select and fit its rung under the common evaluation target;
  calibrate coverage and confidence; freeze and audit; write the model, the provenance, the timings
  and the selection table.
- **`predict_poses`**, a private program: `run_predict_poses(cline)` with a thin commander. Inputs:
  `projfile`, `infile` (the model), `mskdiam`, `trs`, `nthr`, the same references. Per selected row:
  the identity and eligibility check, the incumbent kept for development and calibration rows; the
  features under its own half's band and weights; the forward pass of that half's promoted student;
  the distribution calibration; the smallest set holding at least 99% of the calibrated mass, at
  least `topk` (3). A set that needs more than `KMAX_STUDENT` (64) classes, or a half or stratum that
  has not passed the empirical coverage audit, is recorded as a fallback without truncating or
  renormalising the set. For an eligible set, every candidate state's reprojection is searched over
  the in-plane rotations and the shifts with the existing objective, the shift map and the polish,
  under the full shift support of section 5. The confidence is evaluated for the returned full pose;
  the row is accepted only if both the coverage and the seed-risk gates pass. The row's `corr`, the
  calibrated seed confidence (`pconf`), the route, the candidate count and the retained mass are
  recorded; `pconf` and the route are numbered ori slots, the counts may live in the part's
  diagnostics. The accepted state is hard; an abstention leaves the row to the incumbent. The
  program writes a part file the caller merges, never the caller's project.
- **The student as a pruner.** A broad prediction may still allow a bounded candidate search, but
  searching 40 directions is useful only when that set has demonstrated inclusion of the relevant
  solution; high normalised mass alone does not prove inclusion. Inclusion is reported against the
  common full search and the truth, per state, view and acquisition condition. A broad or
  overconfident prediction that fails this test takes the fallback, and no speedup claim counts
  omitted correct directions as saved work.
- **The confidence is for a recoverable full seed.** Success is the correct state and a symmetry-aware
  SO(3) error and native-pixel shift error inside predeclared margins of the consumer's actual
  `athres_cont` and `trs`, followed by recovery in the paired local-refinement test. For p07 the
  rotation bound is 10°, so a direction-only criterion on the 6° grid cannot certify it, and
  `refine=cont` cannot repair a wrong state. On the calibration cohort a small binary calibration is
  fitted on the returned class's score, the candidate search's margin, the retained mass and the
  shift-boundary and ambiguity indicators, through the existing binary logistic core unchanged.
  `pconf_min` is chosen to bound the accepted-seed error and meet the declared useful accepted
  fraction, then tested on the untouched audit micrographs, each active state included. The model's
  mode mass and a per-cone median and MAD score are diagnostics, not calibrated probabilities of
  success. Insufficient calibration support, unresolved competing states or a failed audit declines
  the gate. Every limit, coefficient, threshold, coverage interval and decline reason is stored in the
  model, and a different consumer or capture range requires its own calibration and audit.
- **One source per particle.** Every particle's pose comes from exactly one route, the student's or
  the incumbent's, never a blend; the route is recorded per row and counted per publication in the
  log with the counts that decided.
- **The arrangement, and cryo-EM's gold standard** (decision 11). XEMA scores each half by the other
  half's model to avoid in-sample scoring. Cryo-EM additionally needs independent reconstruction
  noise: training on half B's map-derived labels and applying to half A introduces a cross-half
  dependency that a low-pass cutoff and a later local refinement do not provably remove. The
  production arrangement is therefore same-half throughout: maps, spectra, the data-dependent
  design choice and the gate. Any common initial model, frame or state alignment or shared
  low-resolution reference upstream is declared with its cutoff; same-half training cannot repair a
  dependence inherited there. Cross-scoring is an offline benchmark only, and its map comparisons and
  agreement statistics are not gold-standard half-to-half resolution.

## 9. Where it fits

- `src/main/learn/`, picked up by the recursive glob of `src/CMakeLists.txt` with no CMake change:
  `learn.inf` (one `!@descr:` line), `simple_polar_feats.f90`, `simple_circ_cnn.f90`,
  `simple_pose_student.f90` (the model: the two students, the grid, the calibration, the file I/O), their testers, and a `CLAUDE.md` importing a `simple-main-learn` skill proposed with the
  delivery and applied on approval (`scripts/check_guide_imports.py` checks the import).
- The search layer: the deterministic teacher entry beside `simple_eul_prob_tab.f90`, extending the
  diagnostic posterior's deterministic path and reusing the calculator's grids, the shift map's
  kernels and the row statistics; the production sampler's semantics unchanged. A sparse getter
  comes only with its tail-preservation contract.
- The strategy routines in `src/main/strategies/parallelization/simple_strategy_pose_student.f90`
  (shared memory: one job on the 3D jobs' threads), the commanders in
  `src/main/commanders/simple/simple_commanders_pose_student.f90`, registered with `add_private_prg`
  and the driver's `case`; the required keys of a private program as
  `.github/skills/simple-main-ui/references/cli_change_chain.md` lays out.
- **Consumer 1, the stream** (p07, section 10 of the ingestion policy):
  - After the first assignment's local refinement has been taken (`take_refine_result`), the stage
    starts a `learn_pose_student` job in `assign/it_<n>/student/`, with the 3D jobs' resources, on
    the assignment's project; the refinement's per-half state volumes are the teacher's references
    and are retained with their provenance. The teacher runs under euclid: the assignment's
    project gets its spectrum from the sigma2 bootstrap (`calc_pspec`), as the pool does for its own
    2D, before the fill. The model waits for the next publication, whose new
    rows are the audit, before promotion. A promoted model becomes the stage's `student.model`,
    persisted across restarts with its artifacts and cohort identities; an unavailable or declined
    half keeps its incumbent route.
  - On every later publication, before the expansion: a `predict_poses` job on the merged
    publication's selected rows, each judged by its own half. Initially the existing expansion and
    sub-class matching run on the complete selection, preserving the fallback averages and giving a
    paired control; accepted student poses are merged only after that route's normal selection and
    rejection, and other rows keep its results; all selected rows go to the local refinement with
    the existing seeds. Removing accepted particles from the expansion changes the fallback averages
    and needs a separate validation before it counts as saved work. The student does not bypass the
    rejection by fit and is not a particle-quality selector in this delivery.
  - The log counts accepted and fallback rows, coverage failures, the full-search audit results, the
    per-state populations, the candidate counts and the times of every step, with the acceptance
    over the eligible new rows and over the whole publication reported separately. The initial
    full-control integration may cost more than the incumbent; promotion as an accelerator requires
    measured net savings after a validated execution optimisation. `STUDENT_RETRAIN` (5)
    assignments is a proposed cadence, subject to the measured break-even point; a retrain uses the
    newest per-half references and fresh held-out cohorts, and the student's confidence is never a
    training-label filter. The GUI reports the route actually applied and whether a model is still
    under audit.
- **Consumer 2, batch**, not in the first delivery: the registration pass of solve3D_addon and
  refine3D's neighbourhood search seed a particle never searched from the student's directions.
- Documents updated with the delivery: `stream_3D_ingestion_policy.md` section 10 and
  `p07_overview.md`; a page `doc/algorithms/pose_student.md`; this note moves to `completed/` with
  its measurements.

## 10. Validation: the reads, declared

Six tiers, with the policy and endpoints frozen before the first benchmark. Methods are compared on
the same rows, with uncertainty from micrograph-level paired resampling and multiplicity handled as
in section 7. The numerical floors below are proposals, not measurements; decision 16 fixes the
remaining margins, risk limits, useful accepted fractions and cohort minima before a promotion run.
Changing a limit after observing an audit requires a new policy version and new audit data.

1. **Unit** (the fast gate, a `pose_student` suite registered in a suite table of
   `simple_commanders_test_class.f90`): a central-difference gradient check of every layer on a tiny
   network (relative error below 1e-2 in single precision; the dense layers checked exactly in
   double); a toy problem (projections of two Gaussian-blob volumes at SNR 1 over 20 directions and 2
   states, accuracy at least 0.95 on held-out images); the same model bit for bit at one thread and
   at a team of 4; a model-file round trip exact; the shift augmentation checked before the angular
   resampling (a shifted image agrees with the phase-ramped ring to single precision); angular
   samples distinguished from their FFT coefficients, integer-bin equivariance and the measured
   arbitrary-angle interpolation error; small independent oracles for the deterministic teacher scores
   under euclid, repeated under different RNG scopes, batches and thread teams to show that
   no draw enters a label; uniform and broad target rows keeping their tails and triggering the
   fallback when 99% needs more than 64 classes; no uncorrected class weighting; an incomplete
   calibration declining; an unaccepted search never setting `shsig` or shrinking the search support;
   half B's data and maps perturbed without moving half A's production decisions, changed
   halves refused and legitimately appended rows allowed; repeated fitting, development or
   calibration identities taking the fallback in a new publication; split disjointness and paired
   resampling at micrograph granularity; fixtures for the grid's symmetry handling (group-equivalent directions, distinct
   antipodal views) and for the shift's sign, frame and sampling round trips.
2. **Truth, the Ribosembly 20,000 subset, then all 80,000.** The deterministic teacher at each band
   with 1,000 directions and a student per half, evaluated on its own disjoint held-out cohorts.
   Truth-map training is a labelled oracle experiment, not proof that finite-data half references
   are independent; the run is repeated with independently reconstructed half references to test the
   production provenance. Reported first: the teacher's resolvability and its error against the
   truth. Then every rung against the same truth poses, the complete rotation and shift included:
   the accepted and fallback fractions, the actual candidate counts and omitted mass, shortlist
   inclusion of the truth and of the common full-search solution, and the joint recoverable-seed
   accuracy before and after the actual one-iteration local refinement, with a proposed accepted
   joint floor of 0.9 and separate per-state limits (separate 0.9 direction and state accuracies do
   not imply a 0.9 joint rate); shortlist inclusion audited against the 99% target with its
   uncertainty, not claimed from softmax mass; the medium and small states' recall, full confusion
   matrices and per-view coverage; the accuracy against the true shift in 2-pixel bins and the candidate search's
   shift recovery against the truth; the calibration reliability and the
   risk against accepted fraction; the excess error, shortlist inclusion and calibration within
   defocus bins against the teacher and the full search in the same bin, with held-out defocus
   ranges and acquisition blocks; the CTF ablation, the chosen rung on the multiplied against the
   flipped input, compared by paired differences; after refinement, the FSC at 12 Å against the
   truth maps within a proposed 0.03 of the truth-pose reconstructions (0.934, 0.985 and 0.964), with
   per-state maps and populations checked; and the same with the truth labels as the teacher, to
   separate the student's error from the teacher's.
3. **Truth, easy:** `/data/simulated_vcrnapol`, the same measurements and cohort separation, as a
   favourable control, not a substitute for the hard assembly and rare-state tests.
4. **Real:** `/data2/trpv1_c2` and `/data/bgal`; the hard cases where the band carries the least,
   `/data/exp_gate_FlhB` and `/data/10335-streptav`; and two CESPED entries, one soluble (10648,
   PKM2, D2) and one membrane protein (11120, the TSHR complex, the benchmark's hardest), each half
   trained and judged within itself, with, as CESPED reports them, the FSC of the map from the
   predicted poses against the reference map beside the reference halves' own FSC, so that the
   numbers compare with Image2Sphere's and cryoPARES's on the same particles. Eligible frozen half
   references and untouched acquisition blocks; every band against the same full-search evaluation
   and the incumbent route; full-pose recovery, candidate inclusion, accepted fraction, per-state and
   per-defocus excess errors, calibration, fallbacks and artifact provenance; paired one-iteration
   refinements from student and sub-class seeds, with a proposed limit of no more than 0.1 Å worse
   and per-state map and population checks. Independent half-map FSC requires the provenance
   contract; a shared-reference run is an engineering comparison, so labelled.
5. **Stream:** `/data2/stream_holjunk` through p07 with the student: the confidence's AUC against the
   sieve's accept/reject decision as a diagnostic, not a truth label for pose correctness; agreement
   and paired outcomes with the sub-class route; the chain test's success unchanged; restart,
   declined-model and insufficient-coverage fallback; new micrographs and repeated model versions
   without recycling an audit cohort; the fallback averages and the normal rejection keeping their
   support. The actual p07 route is timed on matched hardware and publications, including the
   reference preparation, every teacher band and fold, all fits, calibration, audits, inference,
   candidate searches, retries, expansion, fallback, local refinement and retraining, as wall time,
   core time, peak memory, accepted fraction and candidate-count distribution. An accelerator passes
   only when the training and retraining time plus the student-route publication time is below the
   incumbent's total over the declared acquisition horizon at the required quality; a nonpositive
   per-publication saving has no break-even point, and then no retrain cadence is justified.
6. **Nightly:** a `pose_student` route beside `multistate_workflows` on the simulated
   beta-galactosidase states: each half fitted and predicting its own disjoint holdout, asserting
   candidate coverage, accepted joint seed accuracy, useful accepted fraction, smallest per-state
   recall and fallback, within the existing entry and area suites rather than a new CTest process; an
   optional cross-half benchmark separate from the production test; the actual runtime recorded,
   as `doc/policies/test_environment_policy.md` requires.

## 11. Order of work

- Phase 0, freeze the teacher, the coordinate conventions, the provenance, the
  cohorts and the acceptance policy of decision 16. Build the deterministic teacher entry and
  `simple_polar_feats` with their testers, and run the teacher on the Ribosembly subset at each
  band. Report the teacher's error against the truth, its resolvability per band, the preserved
  target mass and the teacher's cost against p07's actual route. These data decide nothing about
  a model; they say what the labels carry before a network exists.
- Phase 1, `simple_circ_cnn` with its tester in the fast gate (tier 1), then the `cnn16` rung at
  15 Å on the subset, fitted and evaluated within each half's disjoint cohorts: candidate
  inclusion, full-pose recovery and the accuracy against the shift, as exploratory data that
  certify no tuned model on the same audit cohort.
- Phase 2, `simple_pose_student`, `learn_pose_student` and `predict_poses` with the independently
  selected ladder, the calibration, the untouched audit and the selection table; the Ribosembly and
  vcrnapol measurements (tiers 2 and 3). Gate: the predeclared tier-2 criteria on new audit data.
- Phase 3, consumer 1 in p07, the stream_holjunk run and the real-data measurements (tiers 4 and 5),
  the small particles and the CESPED entries included; first the reference provenance, then the run
  with the complete incumbent control. Gate: tiers 4 and 5, the net timing included; a later
  optimisation that changes the expansion's membership repeats the fallback and quality checks
  before claiming its saving.
- Phase 4, the nightly route (tier 6), the policy and overview updates, the algorithms page, the
  skill proposal.
- Phase 5, later: consumer 2.

Each phase ends with the fast gate, the touched unit suites and `git status`; the development stays
uncommitted for review as one diff.

## 12. Open decisions for Hans

The values below are proposals. Teacher determinism, preserved uncertainty, same-half production,
the untouched audit and capture-aware acceptance are design constraints. Decision 16 is the
mandatory acceptance policy; unresolved acceptance values block promotion, not exploratory
measurements.

1. The folder and the names: `src/main/learn/`, `simple_polar_feats`, `simple_circ_cnn`,
   `simple_pose_student`; the programs `learn_pose_student` and `predict_poses`.
2. The student grid: 1,000 directions, the class-average matching's. The alternative is 500 (a
   cheaper last layer, coarser seeds).
3. The bands of the ladder: 15, 10 and 8 Å in that order, or 15 Å alone.
4. The labels: deterministic profile-score soft targets at the student band, with the shift map and
   polish for the `NTEACH_SHIFT` (64) best classes per particle under one declared rule. A
   marginalised teacher or refined one-hot labels are separate measured alternatives; the sampled production fill is never a draw-free teacher.
5. The confidence's ori key `pconf` and the route key, written to `ptcl3D`.
6. The state from the student: a hard label (the frame of `multistate_posterior_state_draw.md`: hard
   labels, no weights), or a soft weight.
7. The model file: one text header with two binary weight files, or one file.
8. `predict_poses` runs the in-plane and shift search at the predicted directions itself, or hands
   the directions to the neighbourhood search of refine3D.
9. The retraining cadence: 5 assignments provisionally, contingent on drift and the measured
   break-even.
10. The shell weights: the search's own `k/sigma2(k)` with the half's global spectrum frozen in the
    model, or the particle's own ring powers.
11. Production scoring by each half's own student against its own parity references, the matcher's
    normal mode; the stream retains the local refinement's per-half volumes as the teacher's
    references. Cross-scoring is an offline benchmark with no gold-standard resolution claim.
12. The fitting floors: `MIN_HALF_PTCLS` 5,000 and `MIN_HALF_STATE` 500; the micrograph and per-state
    support floors of the calibration and audit cohorts under decision 16.
13. The accepted physical search supplies the shift and `shincr`; `shsig` stays unknown unless a
    validated objective-based uncertainty is available.
14. The CTF as the third input channel, or the CTF-multiplied particle alone.
15. `KMAX_STUDENT` (64): the largest set of classes searched per particle; insufficient 99% coverage
    within it always falls back.
16. Before any promotion benchmark: the confidence level and the micrograph resampling count; the
    adjustment for the finite ladder; the useful-improvement and cost rule; the full-rotation and
    shift capture margins; the maximum accepted joint seed error; the minimum useful accepted
    fraction; the per-state, rare-view and shortlist-inclusion limits; the calibration and audit
    micrograph minima; the reconstruction noninferiority and the acquisition horizon of the
    amortised timing. The frozen values are recorded in this note and in the model policy before
    use, with fresh audit data for any later revision.

## 13. Risks and limits

- **Low SNR.** The student cannot create missing pose information. It must supply a full seed inside
  the local refinement's capture range or a candidate set with demonstrated inclusion; otherwise the
  incumbent route supplies the row, and view-only grid agreement does not qualify a seed.
- **Small particles and small membrane proteins.** At 15 Å a micelle gives the tilt and not the
  azimuth; the finer bands are rungs and the teacher's resolvability is a conditional diagnostic.
  Pruning is allowed only when its candidate-inclusion audit passes. The small-particle datasets and
  the CESPED membrane protein of tier 4 test whether any useful accepted fraction remains, rather
  than assuming that a broad ring fits in 64 classes.
- **States.** A single image may not distinguish conformations. Soft training targets and competing
  state and direction hypotheses are kept until the physical candidate evaluation; a hard wrong state
  is not repaired by `refine=cont`, and the medium and small Ribosembly states' recall and the
  accepted-seed error are gates.
- **Preferred orientation.** The probabilities keep their meaning rather than being reweighted by
  inverse frequency; rare-view support and fallback are reported explicitly, since a high global
  accepted accuracy can conceal missing directions and poor maps.
- **Defocus as a confounder.** A classifier over CTF-multiplied inputs can learn the defocus where
  the training set confounds it with the class: states that sit in different micrographs, views that
  cluster by acquisition. Stratification and the CTF channel do not prove separation; the excess
  error and the calibration are compared within defocus bins, and later acquisition blocks are
  audited.
- **Drift.** Frozen maps and spectra stabilise replay, not future-data validity. Retraining requires
  fresh eligible cohorts and an audit; a new version is never certified on its own accepted
  predictions. p07's fixed session references and the per-half references of the student's
  retraining need compatible provenance before integration.
- **Cost on the Dell** (24 cores, no GPU): every timing is an estimate until tier 5. The
  deterministic nuisance searches, the complete targets, the independent development, calibration
  and audit, the physical candidate search and the preserved fallback support may dominate a small
  CNN, and a faster replacement for a dense particle search need not accelerate the actual
  sub-class stream route.
- **Compile time.** The three modules are small; no new CMake target; heavy arrays allocatable.

## 14. Rules the work follows

- `doc/policies/layering_and_project_ownership_policy.md`: the strategy routines on `(params, build)`
  with `run_<x>(cline)` entries and thin commanders; the programs read the project and write part
  files, the stage merges and writes its project once; the search for an existing step before writing
  one (the matcher's reference preparation, the deterministic objective kernels, the shift map, the
  seed contract of `simple_state_assignment`).
- `doc/policies/compile_time_policy.md` and the `simple-modern-fortran` skill: encapsulated classes
  with `new`/`kill`, heavy arrays allocatable, `use ..., only:` imports, the OpenMP rules, testers
  registered in a suite table, no new CMake target, no per-file flags.
- No new command-line key and no internal switch: the model travels through `infile` and `fname`,
  the thresholds, coordinate conventions and constants that decided live in the model file, with
  versioned policy defaults and explicit provenance. Required reference and calibration metadata is
  never replaced by an undocumented key or a silent fallback to merged inputs.
- The logistic core and the quality model are not changed; the `learn`, `evaluate` and `analyze`
  outputs of `model_cavgs_rejection` stay byte-identical.
- Nothing is committed; the whole development is one reviewable diff.

## 15. Sources

- cryoPARES: Sanchez-Garcia et al., bioRxiv 10.1101/2025.03.04.641536 (2025); the CESPED benchmark,
  arXiv 2311.06194 and its package `~/src/cesped` (`cesped/network/image2sphere.py`,
  `cesped/datamanager/relionStarDataset.py`, `cesped/datamanager/ctf.py`, `cesped/evaluateEntry.py`,
  `cesped/configs/`); cryoFIRE, arXiv 2210.07387 (NeurIPS 2022); semi-amortized pose inference,
  Shekarforoush et al., NeurIPS 2024.
- [CESPED evaluation protocol, section 2.5](https://arxiv.org/html/2311.06194v5#S2.SS5): cross-half
  supervised benchmarking and its map-to-reference comparison do not establish gold-standard
  half-to-half resolution.
- XEMA: `docs/algorithms/13_cross_validation.md`, `docs/algorithms/11_consensus_scoring.md`,
  `docs/design_notes/gold_standard_selection.md`, `docs/design_notes/small_dataset_scoring_model.md`,
  `src/consensus_scoring/xema_consensus_twofold.f90` and `xema_consensus_fitter.f90`.
- SIMPLE: `doc/implementation_notes/completed/polar_shift_search.md` and
  `polar_shift_search_dynamic_centring.md`; `doc/implementation_notes/planned/polar_shift_laplace_width.md`;
  `doc/implementation_notes/planned/multistate_posterior_state_draw.md`;
  `doc/policies/stream/stream_3D_ingestion_policy.md`; the sources read for the teacher's contract:
  `simple_eul_prob_tab.f90` (`score_ref`, `deterministic_best`), `simple_polarft_corr.f90`
  (`gen_posterior_val`, `gen_best_objfun_val`), and `simple_strategy3D_cont.f90` (the fixed-state search).
