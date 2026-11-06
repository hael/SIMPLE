# Learned pose and state assignment: a student distilled from the probabilistic search

Proposal of 2026-11-05, revised with Hans on 2026-11-06. SIMPLE's probabilistic machinery makes
per-dataset training labels, the posteriors over 3D direction and state of a subset of particles, and
a low-complexity machine trained on them assigns the particles streaming off the microscope faster than
the search can, as a first pass with the search behind it. Nothing is built. Sections 1 and 2 are the
idea and what exists, sections 3 to 9 the design, section 10 the validation, section 11 the order of
work, section 12 the open decisions for Hans, sections 13 to 15 the risks, the rules and the sources.

## 1. What is proposed, in one paragraph

A pose student: a small network, trained per dataset on the tempered posteriors of the dense 3D
probability table over a fixed grid of (state, direction) classes, that reads a particle's
CTF-multiplied polar Fourier transform in a low-resolution band and returns a posterior over the grid
and a confidence. Its input is shift-marginalised through the polar shift map against a small bank of
templates, and rotation-invariant through the cyclic structure of the polar transform. The teacher's
labels come from the search's own table fill at the student's grid and band. Two students, one per
half-set, each judging the particles of the other half, chosen from a ladder of designs by a two-fold
rule with no threshold parameter and published with their agreement. A prediction program gives each
new particle its state, the set of directions holding its posterior and, from a search over the
in-plane rotations and shifts at those directions alone, a full pose and a confidence; the confident
particles become seeds of the existing local refinement, the rest take the route they take today, and
where the student cannot assign it prunes, reducing a search over thousands of directions to tens. The
probabilistic search keeps two roles: it makes the labels and it judges the hard tail. No new
command-line key: the model file travels through `infile` and `fname` as the class-average quality
model's does, and every threshold and constant that decided lives in the model file.

## 2. Where the idea stands

### 2.1 Prior work

- Unsupervised amortized pose inference (cryoPoseNet, cryoAI, cryoFIRE): an encoder learns the poses
  jointly with the map, with no labels. Encoders converge slowly and get stuck; the semi-amortized
  method of Shekarforoush et al. (NeurIPS 2024) and DRGN-AI both retreat to per-particle optimization
  after an amortized warm start. That is evidence against learning poses from scratch, not against the
  supervised setting.
- Supervised pose assignment, which is this setting: the CESPED benchmark (arXiv 2311.06194) and
  cryoPARES (Sanchez-Garcia et al., bioRxiv 10.1101/2025.03.04.641536). cryoPARES trains on a refined
  apo dataset and infers on ligand-bound datasets of the same protein; a U-Net feeds an Image2Sphere
  head that outputs a distribution over an SO(3) grid of 294,912 rotations; the CTF enters as a
  phase-flipped input channel; in-plane rotation and shift are learned by augmentation. Its numbers:
  network inference above 30,000 particles per minute per GPU; the raw top-1 pose within 5° of the
  RELION pose for about 30% of the particles, above 60% after a local refinement; resolutions 0.1 to
  0.6 Å worse than RELION on six datasets; a direction-normalised confidence rejects junk as well as 2D
  classification does; two RELION runs disagree by more than 5° on over 20% of the particles. Trained
  per protein, up to 100 epochs on GPUs, checkpoints of 40 GB.
- Lessons taken: the raw pose is coarse and a local refinement supplies the precision; the confidence
  is worth as much as the pose; the student is judged against posteriors, not hard labels; the model is
  dataset-specific, so it is trained per dataset and retrained when the references change.

### 2.2 What SIMPLE has

- The dense 3D table (`src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90`): its grid is
  `build%eulspace` at `params%nspace` times the states; `fill_tab_range` fills a particle's row over
  every (state, direction) reference with the in-plane and shift search done; `posterior_row_stats`
  gives the row's tempered posterior (the mode, its mass `pmax`, `neff`, `k99`), `top_posterior_extras`
  the weighted top-k. The keyed draws are stochastic at `prob_beta` times the particle temperature by
  design. The registration pass of solve3D_addon (`set_refine3D_register_policy` in
  `src/main/solve/simple_solve3D_controller.f90`) is one full-cohort `prob` iteration at a fixed
  `NSPACE_REGISTER` (1000) against frozen references: a table fill at a chosen grid, as a step.
- The polar calculator (`src/main/pftc/simple_polarft_calc.f90`): the particle sections over the
  Friedel half of the circle, the memoised full-ring spectra of the CTF-multiplied particle
  (`ft_ptcl_ctf`) and of the references (`ft_ref`), the CTF matrices (`ctfmats`), the sigma2 shell
  weights as two global spectra, one per half (`sigma2_noise(shell, half)`, `memoize_sigma2_weights`),
  the shifted particle (`memoize_shifted_ptcl`), and the correlation over all in-plane rotations with
  the shift search per reference.
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
  state references (`LP_ASSIGN` 15 Å, `NSPACE_ASSIGN` 1000 directions), the particles inherit the
  sub-class state and pose, then one iteration of `refine3D refine=cont` under objfun=cc from the
  seeds `set_refine_seeds` sets (the halves of the partition, one per micrograph by `eo_transaction`;
  the nearest grid direction). The sub-class pose on a 6° grid is the bar a student seed has to meet.
- The logistic core (`src/main/quality/simple_quality_logistic.f90`), binary and at most sixteen
  features, is not the student's shape; the learn/evaluate/promote lifecycle and the versioned
  key=value model file of `src/main/cavg_quality/simple_cavg_quality_model.f90` are the pattern.
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
  streptavidin `/data/10335-streptav`; `/data2/stream_holjunk` for the stream route with junk;
  `multistate_workflows` (`src/main/commanders/test/simple_commanders_test_highlevel_msworkflows.f90`)
  as the pattern for a nightly route on simulated beta-galactosidase states.

### 2.3 What XEMA's two-fold method contributes

XEMA (`~/src/XEMA`, Hans's cross-linking mass-spectrometry engine) scores crosslinks with a learned
model under one rule: every choice the engine makes about a dataset is made inside one half of the
data and judged on the other, and the agreement of the halves is published beside the result
(`docs/algorithms/13_cross_validation.md`, `docs/design_notes/gold_standard_selection.md`,
`src/consensus_scoring/`). It took the rule from cryo-EM's gold standard, and the student takes it
back: label-free halves with a fingerprint checked on replay; choices made on inner quarters split by
the next hash bit; a ladder of designs visited simplest first under an Occam rule with no threshold
parameter and a joint choice of the two halves; an agreement statistic in the FSC's seat, reported and
never selected on; gates chosen per half and applied to the other, with the result declined when the
halves do not agree the gain is real; decline by counts before any fit; one source per row, never a
blend; fits bit-identical on any thread count by per-chunk partial sums reduced in chunk order; reads
declared before the run and a versioned record of every constant that decided. Each appears below at
the place it governs.

## 3. The student's input

`simple_polar_feats` (`src/main/learn/`): one routine builds the input of a particle from the polar
calculator, called by training and prediction alike so that the two never drift apart.

1. **The band.** Rings `kfromto(1)` to `box*smpd/lp` with `lp` the student's low-pass. The band is a
   rung of the ladder (section 6), 15 Å first, then 10 and 8 Å; the number of rings depends on
   `box*smpd/lp` alone, about 25 at 15 Å for the benchmark's box 128 at 3 Å per pixel.
2. **The CTF-multiplied particle.** The input is the particle times its CTF, the calculator's memoised
   `ft_ptcl_ctf`, with the CTF itself (`ctfmats`) as a third channel. It is the search's own
   convention: the calculator compares the CTF-multiplied particle with the CTF-free reference in the
   cross spectrum Q_k = REF CTF conj(X), and the shift map sees nothing else. It is also the matched
   filter: with X = CTF P + N the product CTF X = CTF^2 P + CTF N weights every shell by the CTF before
   any template is applied, the optimal linear detector of CTF P in white noise and, with the shell
   weights below, in coloured noise. A shell near a zero carries no signal and enters at no weight,
   where a phase flip would feed its noise at full amplitude. The information per shell is unchanged,
   the multiplication being invertible away from the zeros; what changes is the conditioning of what
   the first layer has to learn. The product's envelope is CTF^2, which moves with the defocus and the
   astigmatism, so two images of one view at two defoci differ in a way the network must not read as
   two views: the CTF channel explains the envelope to the network rather than leaving it to be
   learned as noise. Nothing is restored, by phase flip, Wiener filter or deconvolution, in any rung:
   the student's output is a class, not an image, and Wiener theory has no claim on the
   representation.
3. **The shell weights.** The shells are weighted as the search weights them, `k` under objfun=cc and
   `k/sigma2(k)` under objfun=euclid, where sigma2 is the half's global spectrum, never a per-particle
   one: the calculator's weights are two spectra, one per half, and a part's per-particle residual
   records exist only to be summed into them. An unseen particle therefore needs nothing beyond its
   half, assigned at import, and its own weighted power (`wsqsum_ptcl_band`), which scales the band.
   The weights belong to the model, not to the project: each student's record holds the weights it was
   trained with, and prediction applies the weights of the student it applies, never the project's
   current sigma2, so a later refinement that moves the spectrum does not move the student's input
   distribution and a retrain picks the new spectrum up. The stream, whose local refinement runs
   objfun=cc, trains and applies with `k`; the batch consumer under euclid with the spectrum; a consumer
   whose objective differs from the model's is refused, not silently re-weighted. The logarithm of the
   weighted power is one scalar side input: it carries the junk signal.
4. **The full circle.** The Friedel half (`pftsz` angles) extended by its conjugate to `nrots` angles,
   then averaged along the angle to `NANG` (128) samples, so that an in-plane rotation is a cyclic
   shift of the tensor and a circular convolution along the angle is exact.
5. **The tensor.** Three channels (the real and imaginary parts of the CTF-multiplied particle, the CTF)
   by the rings by `NANG`, about 9,600 values at 15 Å for the benchmark's box. The phases are kept
   throughout; nothing is reduced to a power, which at 15 Å would see a particle's size and mass and
   separate no conformational states of one particle. Augmentation lives in the same module and is
   exact in this domain: a shift is the calculator's own phase ramp (`memoize_shifted_ptcl`), an
   in-plane rotation a cyclic shift.

Under a bank rung (section 5) the network's input is not this tensor but the bank's rotation profiles
at the particle's shift cell, built from the same memoised spectra under the same weights.

## 4. The class grid and the teacher

- **The grid.** `oris%spiral` at `NSPACE_STUDENT` (1000) directions, the class-average matching's grid,
  about 6° apart, over the whole sphere in c1 and the asymmetric unit as `build%eulspace` is built for
  a symmetric point group, times the states: `nstates*NSPACE_STUDENT` classes, 3,000 for three states.
  The in-plane angle is not a class: the student is invariant to it by construction and the search at
  the predicted directions recovers it (section 8).
- **The teacher's fill is at the student's grid.** The particles of a converged run were searched at
  the run's grid, 2,500 or 5,000 directions, and at the run's band. Those rows are not the labels.
  `learn_pose_student` runs its own dense fill at the student grid (`nspace` = `NSPACE_STUDENT`) and the
  student band against the state maps, as the registration pass of solve3D_addon fills at a fixed
  grid, and the posterior rows of that fill are the targets. The grid of the labels is then the grid of
  the classes, the band of the labels the band of the input, the shift cell of each row the one the
  bank's map finds, under the weights of section 3: one consistent object, recorded in the model.
  Coarse-graining a finer row onto the student grid by summing each reference's mass into its nearest
  student direction is a correct marginalisation over the grid, but the row was made at the search's
  band, finer than the student's, and a target sharper than the input can resolve trains an
  overconfident student: the teacher must see what the student sees.
- **What the earlier search contributes.** The state maps, the halves, which particles have a
  converged pose (the training sample), each particle's shift and centring record (the fill's shift
  prior, `shift_prior_sigma2`, and the disk or the window) and the state labels of the balanced draw.
- **The targets.** The fill's posterior row kept sparse (the `NTARGET` (32) entries of largest mass,
  renormalised) and smoothed with 0.05 of the mass spread uniformly over the grid. Sample weight 1;
  class weights inverse to the direction frequency in the training set, capped at 10, so that the rare
  views count. Never the keyed draws: a student trained on draws learns the draw's noise.
- **In the stream** the same fill, on the assignment's project against the stage's state volumes
  (decision 4 of section 12); the refined poses of the local refinement as one-hot targets on the
  nearest grid class are the alternative.
- **The teacher's resolvability is read before anything is trained.** The fill gives, per row and band,
  `pmax`, `neff` and `k99`; their distributions over the training particles say whether the band
  carries the direction at all, for free and before the ladder runs. A band at which the median `k99`
  spans a ring of directions is reported as such, and the student's claim for that dataset is bounded
  by it (section 13).
- **The consumer's grid is decoupled** from the student's by the local search: the student's directions
  seed an in-plane and shift search and then `refine=cont` within `athres_cont`, and refine3D's own grid
  (`set_projs` picks the nearest direction in the job's grid) can be 2,500 or 5,000 as the stage likes.
- **Cost.** One dense fill at 1,000 directions and the band on at most `NTRAIN_MAX` (200,000)
  particles, drawn by whole micrographs as the solve3D sample draw balances, with the shift map. The
  band keeps the rings few, so a fraction of a second per particle per core: under an hour on the Dell,
  once per training and once per retrain. It is the expensive step the whole design amortises.

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
- The map is linear in the reference, so a reference may be anything in the band: a reprojection, a
  class average, a principal component of the reprojections, or a learned template. Its gradient with
  respect to the reference is the adjoint of the same transform.
- The CTF enters as the search applies it, in the cross spectrum, per particle and exact.

Three levels of shift handling, each a rung of the ladder (section 6):

1. **Augmentation only (`cnn16`, `cnn32`).** The network on the tensor of section 3, trained with shifts
   drawn uniformly in the disk through the calculator's shift matrices and rotations as cyclic shifts:
   it learns tolerance. The shift is found afterwards, in the search at the predicted directions, by
   the map and the polish. The cheapest rung, and the one whose accuracy falls with the shift; measured
   first.
2. **A fixed bank through the map (`pca8:linear`, `pca8+cnn16`).** A bank of `NBANK` (8) templates in
   the band, the leading eigenvectors of the Gram matrix of the state references' reprojections over
   the student grid (LAPACK, linked), each stored as a reference of the calculator (`set_ref_pft`).
   For one particle `gen_shiftmap` gives every template's cross term over all rotations and the disk.
   The bank's energy, the sum over templates of the squared cross term, is read per cell, and its best
   cell is the particle's shift (`shiftmap_mode`, `shiftmap_shift`): the subspace centring of the
   particle without a direction hypothesis. The templates' rotation profiles at that cell, `NBANK` by
   `nrots` real values, are the student's input: shift-marginalised exactly, rotation-equivariant
   through the FFT over rotations, the phases of the image kept in the correlations. The `linear` head
   takes rotation-invariant summaries of the profiles (each profile's maximum, mean and the magnitudes
   of its three lowest circular harmonics, about 40 numbers) into a multinomial softmax: the simplest
   phase-aware rung. The `cnn16` head takes the profiles whole into the circular convolution along the
   rotation axis with the templates as channels. No gradient enters the map in these rungs.
3. **A learned bank (`bank8+cnn16`).** The same with the templates trained: the map's adjoint carries
   the gradient to the templates, and the chosen cell is held fixed within a step. Only if the fixed
   templates leave accuracy on the table.

A two-pass scheme (a tolerant first pass, the map at the predicted reprojections, re-centring, a second
pass) is rung 1 followed by what prediction does anyway, and is not a rung of its own.

**The student's shift is the particle's first centring.** Under a bank rung `predict_poses` writes the
bank's cell as the row's shift, with `shincr` the move and `shsig` the width of the bank's energy over
the cells at the physical temperature, as the dynamic centring keeps them (the centring belongs to the
particle, not the stage), so that the local refinement and every later window search start centred. A
particle whose energy map has no clear cell, a width above the window, is unconfident whatever its
direction posterior says.

**Cost.** Per template the map costs (N+1) nrho (K nrots multiply-adds and one FFT) plus the readout:
at 15 Å with 20 orders, 8 magnitudes, 25 shells and about 400 rotations, 10 to 15 million
floating-point operations, so a bank of 8 is a few times the network and ten to a hundred times below
the search over thousands of reprojections.

## 6. The ladder and the network

**The ladder.** The designs, simplest first: `pca8:linear`, `cnn16`, `pca8+cnn16`, `bank8+cnn16`,
`cnn32`, each at the bands 15, 10 and 8 Å, visited in order of cost: every design at 15 Å, then at
10 Å, then at 8 Å. Which rung a dataset gets is chosen by the two-fold rule of section 7, never by a
constant. A small membrane protein in a micelle is, at 15 Å, the micelle, whose outline gives the tilt
of the membrane normal and nothing of the azimuth about it; the azimuth sits at 8 to 12 Å, in the
transmembrane helices, which is why the search's late stages run there, and why the band is a rung.

**The network** (`simple_circ_cnn`, `src/main/learn/`): an encapsulated class with private components,
`new`/`kill`, forward, backward, one Adam step, read and write.

- On the tensor of section 3: three blocks of (a 3 by 3 convolution with circular padding along the
  angle and zero padding along the rings, bias, ReLU, 2 by 2 average pooling), with 16, 32 and 64
  channels; then the mean over the angle (the rotation invariance; the pooling along the angle makes
  the equivariance approximate below the pooled step, which is finer than the grid); then a dense
  layer to 256 with ReLU, the side scalar appended, and a dense layer to the classes with a softmax.
  About one million weights, most of them in the last layer; about 25 million floating-point
  operations per particle at 15 Å, more at the finer bands with the input. On the bank's profiles: the
  same blocks as one-dimensional circular convolutions along the rotation axis with the templates as
  channels.
- The convolution as im2col plus `sgemm`, the backward pass as the two transposed products; weights,
  activations and gradients `real(sp)` and allocatable; one Adam state per weight array.
- **Bit-identical on any thread count.** The mini-batch is cut into fixed chunks; the chunks are
  computed in parallel, each into its own gradient buffer indexed by chunk, not by thread, with
  per-thread activation buffers sized by `omp_get_max_threads()` read before the region and the BLAS
  on the calling thread inside it; the chunks are reduced serially in chunk order after the region.
  The same sum at one thread and at any team, bit for bit, and the tester asks for identity, not a
  tolerance.
- The loss: the class-weighted cross-entropy against the soft targets; Adam at 1e-3, batches of 256,
  weight decay 1e-5; early stop on a held-out tenth of the training half; at most 30 epochs.
- **The model file.** A key=value text header (`model_version=1`) and the weights of the two students
  as two binary files beside it (decision 7). The header records everything that decided: the design
  and the band, the grid with `nspace`, `nstates` and `pgrp`, `NANG`, the shell weights with the
  spectrum they came from, the bank's templates, the floors and the margin rule, the halves'
  fingerprint and balance, both confidence thresholds and the decline, and the selection table of
  section 7, so that a prediction is a replay of a stated method and the reader is fail-closed.

## 7. Two halves: the selection, the agreement and the decline

- **The halves** are the project's even-odd partition, one half per micrograph in the stream
  (`eo_transaction`): label-free and from the data alone. The model file records a fingerprint of the
  training rows' halves and the balance per half (particles, teacher modes per state, the direction
  histogram), and `predict_poses` refuses a project whose halves differ (`ensure_field_halves`).
- **The quarters.** A choice made on the data is made on an inner two-fold split inside each half, by
  the next bit of a hash of the micrograph key: whole micrographs, as label-free as the halves. The
  labels are the half's teacher posteriors, never re-made on a quarter.
- **The Occam rule and the joint choice.** Per half, each rung is fitted on each quarter and scores the
  other; the half's count S_H(d) is the number of held-out particles the rung puts within 2 grid steps
  (about 12°) of the teacher's mode with the state right. A later rung replaces the incumbent only if
  its count exceeds the incumbent's by more than the square root of the incumbent's count, one Poisson
  standard error: no threshold parameter, the noise of the count sets the margin. The choice is joint:
  a rung is chosen only when both halves prefer it to the joint incumbent. The chosen rung is then
  fitted on all of each half; every other rung's outer pair is fitted too, for its agreement alone.
- **The agreement statistic**, in the FSC's seat, is on outputs, since network weights are not
  comparable: both students score each half's held-out tenth (student A has not trained on A's tenth,
  student B never saw half A), and the agreement is the fraction of those particles the two put within
  2 grid steps of each other with the same state, and the mean symmetric Kullback-Leibler divergence
  of their posteriors. Printed for every rung in the selection table (`<fname>_selection.txt`, one row
  per rung: S_A, S_B, the quarter counts, the agreement, each half's preference, the joint choice).
  Selecting on it would favour the dumbest model: two constant students agree perfectly.
- **Decline by counts, before any fit.** A half with fewer than `MIN_HALF_PTCLS` (5,000) training
  particles, or fewer than `MIN_HALF_STATE` (500) teacher modes in an active state, declines, and the
  model file says so; a quarter below the same floors makes every rung unavailable for its half, and
  the first rung stands. A deterministic function of counts the method computes anyway, never a
  comparison of models on the data.

## 8. The two programs, the confidence and the arrangement

- **`learn_pose_student`**, a private program: the strategy routine `run_learn_pose_student(cline)` on
  the caller's `(params, build)` with a thin commander. Inputs: `projfile`, `vol1..volN`, `nstates`,
  `lp`, `nspace` (the student grid), `mskdiam`, `pgrp`, `fname`, `nthr`. Steps: prepare the references
  and the calculator as refine3D's matcher does (`prepare_refs_sigmas_and_pftc` in
  `simple_strategy3D_matcher.f90`); read the teacher's resolvability and the floors; fill the rows of
  the training particles in batches (`fill_tab_range`, no draw) at each band of the ladder; take each
  particle's sparse posterior through one new getter of the table beside `posterior_row_stats`; build
  the bank and the features; walk the ladder on the quarters, choose jointly, fit the chosen rung on
  each half; calibrate the confidence; write the model and the selection table.
- **`predict_poses`**, a private program: `run_predict_poses(cline)` with a thin commander. Inputs:
  `projfile`, `infile` (the model), `vol1..volN`, `mskdiam`, `trs`, `nthr`. Per selected row: the
  features under the model's band and weights; the forward pass of the student of the other half; the
  classes holding 99% of its posterior, at least `topk` (3) and at most `KMAX_STUDENT` (64); for each,
  the correlation of the particle against that state's reprojection at that direction over all
  in-plane rotations with the calculator's shift map and polish, started at the bank's cell; the best
  (state, direction, in-plane angle, shift) becomes the row's pose with `corr` its correlation, the
  shift the row's centring record; the confidence goes to a new ori key (`pconf`) and the route to
  another, both numbered slots declared as the ori module requires. The state stays a hard label. The
  program writes its rows as a part file the caller merges, never the caller's project.
- **The student as a pruner where it cannot be an assigner.** A posterior spread over a ring of 40
  directions, the small membrane protein at a coarse band, costs a search over 40 reprojections instead
  of thousands, still the two orders of magnitude the design is for, and the pose it finds is the
  search's own. The confidence then applies to the pose the search found, not to the student's mode
  alone.
- **The confidence** is the mass of the posterior mode, normalised per direction cone by the median
  and MAD measured on the held-out tenth at training (cryoPARES's direction-normalised score), so that
  a preferred-orientation cone does not look confident by population alone. The threshold `pconf_min`
  is a gate chosen per half and applied to the other: optimised on half A's held-out tenth and on half
  B's, cross-applied, their distance published. The companion read: on each half's held-out tenth, the
  gain of the student route over the incumbent route (in the stream the sub-class inheritance, whose
  pose exists for the same rows) in particles within 2 grid steps of the teacher. A gain no greater
  than the square root of that half's incumbent count, or of opposite sign to the in-sample gain, means
  the halves have not agreed that the gain is real: the model is written with `declined=yes` and the
  reason, the consumer keeps its route for every row, and the log says so. The in-sample winner is
  reported and never applied. Both thresholds and the decision live in the model file, never a
  command-line key.
- **One source per particle.** Every particle's pose comes from exactly one route, the student's or the
  incumbent's, never a blend; the route is recorded per row and counted per publication in the log
  with the counts that decided.
- **The arrangement, and cryo-EM's gold standard** (decision 11). XEMA scores each half by the other
  half's model so that no row is scored by a model that saw it. Cryo-EM's gold standard asks something
  else: that half A's poses not depend on half A's map noise. A student trained on half B's labels and
  applied to half A carries B's labels into A's poses, but the labels are directions on a 6° grid from
  a band-limited search, where the two half-maps are one map, and the local refinement per half
  follows; the stream's sub-class route already crosses the halves at that band. Cross-scoring is the
  proposal because it makes every prediction out-of-sample for the training particles as well, which
  the batch consumer needs; the same-half arrangement is the alternative, as separate as refinement
  itself and out-of-sample for new particles anyway.

## 9. Where it fits

- `src/main/learn/`, picked up by the recursive glob of `src/CMakeLists.txt` with no CMake change:
  `learn.inf` (one `!@descr:` line), `simple_polar_feats.f90`, `simple_circ_cnn.f90`,
  `simple_pose_student.f90` (the model: the two students, the grid, the bank, the calibration, the
  file I/O), their testers, and a `CLAUDE.md` importing a `simple-main-learn` skill proposed with the
  delivery and applied on approval (`scripts/check_guide_imports.py` checks the import).
- The table: the getter of the sparse posterior row in `simple_eul_prob_tab.f90`.
- The strategy routines in `src/main/strategies/parallelization/simple_strategy_pose_student.f90`
  (shared memory: one job on the 3D jobs' threads), the commanders in
  `src/main/commanders/simple/simple_commanders_pose_student.f90`, registered with `add_private_prg`
  and the driver's `case`; the required keys of a private program as
  `.github/skills/simple-main-ui/references/cli_change_chain.md` lays out.
- **Consumer 1, the stream** (p07, section 10 of the ingestion policy):
  - After the first assignment's local refinement has been taken (`take_refine_result`), the stage
    starts a `learn_pose_student` job on the assignment's project against the stage's state volumes,
    in `assign/it_<n>/student/`, with the 3D jobs' resources. The model becomes the stage's
    (`student.model` in the stage's folder), kept across restarts like the references; a declined
    model leaves the route as it is.
  - On every later publication, before the expansion: a `predict_poses` job on the merged
    publication's selected rows. A row at or above `pconf_min` takes the predicted pose, state and
    centring and skips the sub-class route; the rest go through the expansion and the matching as
    today; all go to the local refinement together, whose seeds `set_refine_seeds` sets as before.
  - The log counts the confident rows, their agreement with the sub-class route where both exist
    (the free control: the sub-class assignment runs on every class anyway), and each route's
    per-state populations. The student is retrained every `STUDENT_RETRAIN` (5) assignments from the
    newest maps and poses. The GUI shows "assigning particles with the pose student".
- **Consumer 2, batch**, not in the first delivery: the registration pass of solve3D_addon and
  refine3D's neighbourhood search seed a particle never searched from the student's directions.
- Documents updated with the delivery: `stream_3D_ingestion_policy.md` section 10 and
  `p07_overview.md`; a page `doc/algorithms/pose_student.md`; this note moves to `completed/` with
  its measurements.

## 10. Validation: the reads, declared

Six tiers, declared before the first run. Where a count is compared the margin is one Poisson
standard error, as in section 7; fixed fractions appear only where the truth exists. The numbers are
provisional until Hans rules on the first measurement, as the floors of `multistate_workflows` were
set from its profiles.

1. **Unit** (the fast gate, a `pose_student` suite registered in a suite table of
   `simple_commanders_test_class.f90`): a central-difference gradient check of every layer on a tiny
   network (relative error below 1e-2 in single precision; the dense layers checked exactly in
   double); a toy problem (projections of two Gaussian-blob volumes at SNR 1 over 20 directions and 2
   states, accuracy at least 0.95 on held-out images); the same model bit for bit at one thread and
   at a team of 4; a model-file round trip exact; the feature module's augmentation exact (a shifted
   image's features equal the phase-ramped features to single precision); the bank's map against the
   calculator's shift search on a reference stored as a template.
2. **Truth, the Ribosembly 20,000 subset, then all 80,000.** The teacher is the fill against the truth
   maps at each band with 1,000 directions (the first `ms_posterior` run found the fill's argmin
   perfect at 10 Å with the truth maps); a student per half. Reported first: the teacher's
   resolvability per band. Then, on the other half: the accepted fraction at `pconf_min`; among the
   accepted, the direction within 2 grid steps of the truth in at least 0.9 of the cases and the state
   right in at least 0.9 (the medium assembly's recall on its own: it is the hard one); the accuracy
   of the augmentation-only rung against the true shift magnitude in bins of 2 pixels; the bank's
   shift recovery, the distance between the bank's cell and the true shift against the cell size; the
   accuracy and the state confusion by defocus bin, which must not vary beyond their Poisson resolution
   (the bench's CTFs are drawn per particle independently of the state, so a dependence is the
   student's own failure); a reliability curve of the confidence in ten bins; the selection table.
   Then `predict_poses` with one local refinement, and the FSC at 12 Å against the truth maps within
   0.03 of the reference FSC of the truth-pose reconstructions (0.934, 0.985 and 0.964). The same with
   the truth poses as the teacher, which separates the student's error from the teacher's.
3. **Truth, easy:** `/data/simulated_vcrnapol`, the same measurement, as the ceiling.
4. **Real:** `/data2/trpv1_c2` and `/data/bgal`, and the hard cases where the band carries the least,
   `/data/exp_gate_FlhB` and `/data/10335-streptav`. The teacher is the fill against the project's
   converged state maps; per band the teacher's resolvability, the accepted fraction and the size of
   the searched set per particle, the agreement of the search over that set with the full search; the
   agreement with the teacher's posterior mode on the other half; the accuracy by defocus bin; and the
   resolution after one local refinement from student seeds against the sub-class seeds of a p07-like
   assignment, not worse by more than 0.1 Å.
5. **Stream:** `/data2/stream_holjunk` through p07 with the student: the confidence's AUC against the
   sieve's accept/reject decision; the agreement with the sub-class route; the stream chain test's
   outcome unchanged.
6. **Nightly:** a `pose_student` route beside `multistate_workflows` on the simulated beta-galactosidase
   states, a few minutes: train on half 1 of the truth, predict half 2, floors on the agreement and the
   smallest per-state recall, registered as `doc/policies/test_environment_policy.md` section 1.5
   requires.

## 11. Order of work

- Phase 0, the feature module and the bank rung: `simple_polar_feats` with its tester (the band, the
  weights, the CTF channel, the full ring, the exact augmentation), then the fixed bank through the
  shift map with the `linear` head, on the Ribosembly subset at 15 Å, trained on one half and scored on
  the other. Gate: the teacher's resolvability, the bank's shift recovery and the linear head's state
  and direction accuracy reported. They decide nothing; they say how much the map's centring gives and
  how much the images carry before any network exists.
- Phase 1, `simple_circ_cnn` with its tester in the fast gate (tier 1).
- Phase 2, `simple_pose_student`, the table getter, `learn_pose_student` and `predict_poses` with the
  ladder, the halves and the selection table; the Ribosembly and vcrnapol measurements (tiers 2 and 3).
  Gate: the tier-2 criteria, or Hans's ruling on the numbers.
- Phase 3, consumer 1 in p07, the stream_holjunk run and the real-data measurements (tiers 4 and 5),
  the small particles included. Gate: tiers 4 and 5.
- Phase 4, the nightly route (tier 6), the policy and overview updates, the algorithms page, the
  skill proposal.
- Phase 5, later: consumer 2.

Each phase ends with the fast gate, the touched unit suites and `git status`; the development stays
uncommitted for review as one diff.

## 12. Open decisions for Hans (the proposals are in force until he rules)

1. The folder and the names: `src/main/learn/`, `simple_polar_feats`, `simple_circ_cnn`,
   `simple_pose_student`; the programs `learn_pose_student` and `predict_poses`.
2. The student grid: 1,000 directions, the class-average matching's. The alternative is 500 (a cheaper
   last layer, coarser seeds).
3. The bands of the ladder: 15, 10 and 8 Å in that order, or 15 Å alone.
4. The labels in the stream: the trainer's own fill at the student grid and band (one table fill per
   training), or the refined poses as one-hot targets.
5. The confidence's ori key `pconf` and the route key, written to `ptcl3D`.
6. The state from the student: a hard label (the frame of `multistate_posterior_state_draw.md`: hard
   labels, no weights), or a soft weight.
7. The model file: one text header with two binary weight files, or one file.
8. `predict_poses` runs the in-plane and shift search at the predicted directions itself, or hands
   the directions to the neighbourhood search of refine3D.
9. The retraining cadence in the stream: every 5 assignments.
10. The shell weights: the search's own, `k` or `k/sigma2(k)` with the half's global spectrum frozen in
    the model, or the particle's own ring powers.
11. The scoring arrangement: each half's new particles assigned by the other half's student
    (cross-scoring), or by its own half's student (section 8).
12. The floors: `MIN_HALF_PTCLS` 5,000 and `MIN_HALF_STATE` 500.
13. The bank: `NBANK` (8) principal components of the reprojections, or learned templates.
14. The student's shift as the particle's centring record (`shincr`, `shsig`), or the shift left to the
    local refinement.
15. The CTF as the third input channel, or the CTF-multiplied particle alone.
16. `KMAX_STUDENT` (64): the largest set of directions the consumer searches for one particle.

## 13. Risks and limits

- **Low SNR.** The student's ceiling is the information in one image. The bar is the sub-class seed on
  a 6° grid, not the search, and the search stays behind the student for the unconfident.
- **Small particles and small membrane proteins.** At 15 Å a micelle gives the tilt and not the
  azimuth; the teacher's resolvability says so before training, the finer bands are rungs, and where
  no band resolves the azimuth from one image the student prunes rather than assigns. The design does
  not claim more than the band carries, and the small-particle datasets of tier 4 are the read.
- **States.** A single image carries little state evidence, which is why the posterior over states
  stays soft even for the full search. The medium Ribosembly assembly is the test; the sub-class
  route's class-level vote stays for the unconfident rows.
- **Preferred orientation.** The class weights and the per-cone normalisation; the rare views may
  still fall to the fallback route, which is the right place for them.
- **Defocus as a confounder.** A classifier over CTF-multiplied inputs can learn the defocus where the
  training set confounds it with the class: states that sit in different micrographs, views that
  cluster by acquisition. The micrograph-balanced draw reduces it, the CTF channel lets the network
  separate it, and the accuracy by defocus bin is a read of every tier.
- **Drift.** The student is tied to its reference set and its spectrum. p07 fixes the references per
  session, the model freezes its weights, and the student is retrained from the newest maps.
- **Cost on the Dell** (24 cores, no GPU): the teacher's fill under an hour per training; training
  200,000 particles for 20 epochs about 3e14 floating-point operations, 25 minutes at 200 GFLOP/s;
  prediction a few milliseconds per particle per core under a bank rung plus the search over the
  predicted directions.
- **Compile time.** The three modules are small; no new CMake target; heavy arrays allocatable.

## 14. Rules the work follows

- `doc/policies/layering_and_project_ownership_policy.md`: the strategy routines on `(params, build)`
  with `run_<x>(cline)` entries and thin commanders; the programs read the project and write part
  files, the stage merges and writes its project once; the search for an existing step before writing
  one (the matcher's reference preparation, the table's fill and the registration pass, the shift map,
  the seed contract of `simple_state_assignment`).
- `doc/policies/compile_time_policy.md` and the `simple-modern-fortran` skill: encapsulated classes
  with `new`/`kill`, heavy arrays allocatable, `use ..., only:` imports, the OpenMP rules, testers
  registered in a suite table, no new CMake target, no per-file flags.
- No new command-line key and no internal switch: the model travels through `infile` and `fname`,
  the thresholds and the constants that decided live in the model file, the floors and margins are
  named parameters of the modules.
- The logistic core and the quality model are not changed; the `learn`, `evaluate` and `analyze`
  outputs of `model_cavgs_rejection` stay byte-identical.
- Nothing is committed; the whole development is one reviewable diff.

## 15. Sources

- cryoPARES: Sanchez-Garcia et al., bioRxiv 10.1101/2025.03.04.641536 (2025); the CESPED benchmark,
  arXiv 2311.06194; cryoFIRE, arXiv 2210.07387 (NeurIPS 2022); semi-amortized pose inference,
  Shekarforoush et al., NeurIPS 2024.
- XEMA: `docs/algorithms/13_cross_validation.md`, `docs/algorithms/11_consensus_scoring.md`,
  `docs/design_notes/gold_standard_selection.md`, `docs/design_notes/small_dataset_scoring_model.md`,
  `src/consensus_scoring/xema_consensus_twofold.f90` and `xema_consensus_fitter.f90`.
- SIMPLE: `doc/implementation_notes/completed/polar_shift_search.md` and
  `polar_shift_search_dynamic_centring.md`; `doc/implementation_notes/planned/polar_shift_laplace_width.md`;
  `doc/implementation_notes/planned/multistate_posterior_state_draw.md`;
  `doc/policies/stream/stream_3D_ingestion_policy.md`.
