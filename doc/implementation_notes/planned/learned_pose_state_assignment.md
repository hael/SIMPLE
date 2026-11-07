# Learned pose assignment: a student distilled from the probabilistic search

Proposal of 2026-11-05, revised on 2026-11-06 with Hans and after a methodology review. SIMPLE's
probabilistic search supplies per-dataset soft labels over the 3D direction, and a small model
trained on them proposes directions for particles the search has not seen, so that a local search
and the local refinement can finish the alignment. The first delivery is single-state: a pose
student, validated for poses before anything else. Whether it saves time is a measurement. Nothing is
built. Sections 1 and 2 are the idea and what exists, sections 3 to 9 the design, section 10 the
validation, section 11 the order of work, section 12 the open decisions for Hans, sections 13 to 15
the risks, the rules, the outlook and the sources.

## 1. What is proposed

A pose student: a small network, trained per dataset on the probabilistic posterior of the dense 3D
table over a fixed grid of directions. It reads a particle's CTF-multiplied polar Fourier transform
in a low-resolution band and proposes a distribution over the grid. Hermitian circular convolutions on the
half ring and pooling provide in-plane rotation invariance on the calculator's own rotation grid; shifts are learned as tolerance from the
calculator's own phase ramps and found exactly by the shift map in the candidate search that
follows. The student's direction and the width of its posterior are handed to refine3D's neighbourhood
search, which finds the in-plane rotation, the shift and its own posterior with the machinery it
has; a particle whose posterior is too broad for a neighbourhood takes the dense search. No new command-line key: the model file travels through `infile` and `fname` as the
class-average quality model's does, and every threshold and constant that decided lives in the
model file.

What sets this apart from the published attempts, which learn in real space, is the polar
implementation: the representation the search already uses, with the in-plane rotation a cyclic
shift of the input, the translation handled exactly by the Fourier-Bessel shift map, the CTF applied
as the search applies it, the sigma2 whitening of the objective, and the table's own posterior as
the teacher. Everything the student sees and learns from is what the probabilistic search computes.

Two questions decide the design's worth, and section 10 asks them first. On real data: split the
dataset into two view-balanced halves, reconstruct with the standard machinery from one half, train
the student on it, use the student on the other half, and compare the result with the standard
machinery's alignment of that half. And how few training images suffice: if the cycle of making
training rows, training, predicting and refining is fast enough on a small set sampled from classes
as the balanced sampling samples, the student can run inside the iterations and stages of refine3D
and solve3D, extending each iteration's alignment of its sampled fraction to the particles it did
not search.

Hans's rulings of 2026-11-06 are recorded in section 12; the open questions are the training-set
size and, later, the production arrangement of the halves.

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
- Related architectures, none of which combines this design's ingredients:
  - Polar images with cyclic convolution: CyCNN (Kim et al., 2020) maps images to polar
    coordinates and wraps its convolutions along the angle, so that a rotation is a translation,
    and beats plain CNNs on rotated benchmarks without augmentation; Polar Transformer Networks
    (Esteves et al., ICLR 2018) do the same in log-polar coordinates and first predict the polar
    origin with a small network, their answer to the translation. In cryo-EM, the Polar
    Transformers of Andén and Sagemüller (arXiv 2506.11283, a NeurIPS oral) use polar CNN feature
    extractors with attention over all relative angular shifts by FFT correlation, equivariant to
    the independent in-plane rotations of every image of a set, for joint denoising; up to a
    twofold MSE reduction at SNR 0.02 and about 20% better ab initio resolution on simulated data.
  - Equivariant pose encoders in cryo-EM: de Ruijter and Cesa (ICLR 2024 GEM workshop) build a
    D4-equivariant encoder on the fact that an in-plane rotation or mirror of the image transforms
    the pose predictably; it converges faster and more often than cryoAI's encoder, in real space
    and on simulated data.
  - Fourier-Bessel, rotation-invariant representations without a network: Zhao and Singer (2014)
    expand images in the Fourier-Bessel basis, where a rotation is a phase shift of the
    coefficients, and build a complete invariant from the bispectrum for viewing-direction
    classification; PolarGMM (2022) clusters noisy projections with a polar-coordinate Gaussian
    mixture fitted by EM. These are the ancestors of the representation here and of flex's polar
    E-step.
  - Complex-valued networks: Deep Complex Networks (Trabelsi et al., ICLR 2018) supply complex
    convolutions, complex batch normalisation and the modReLU activation of section 6, used in MRI
    reconstruction and radar; Harmonic Networks (Worrall et al., 2017) restrict filters to circular
    harmonics so that complex feature maps rotate by a phase under input rotation, the nearest
    thing to features that transform like Fourier coefficients.
  - Not found: a network that imposes the Hermitian symmetry of a real image's transform by
    conjugation-equivariant layers on the half plane (Fourier-domain operators keep the
    non-redundant half of a real FFT, without the equivariance), or one that keeps the CTF in the
    operator and learns from the search's own posterior. With the polar implementation itself,
    those are where this design stands apart (section 1).

### 2.2 What SIMPLE has

- The dense 3D table (`src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90`): its grid is
  `build%eulspace` at `params%nspace`; `fill_tab_range` fills a particle's row over every direction,
  drawing an in-plane rotation per reference from the row posterior at the particle's temperature
  (`score_ref` calls the calculator's `gen_posterior_val`, which requires `objfun=euclid` and takes a
  keyed variate) and polishing a selected subset in shift and rotation through the shift map.
  `posterior_row_stats` gives the normalised row's mode, `pmax`, `neff` and `k99`;
  `top_posterior_extras` supplies the weighted top-k. Skipping `assign_batch` therefore does not make
  the fill a deterministic teacher. The deterministic argmin over the in-plane grid at a given shift
  exists beside it, as the diagnostic posterior's path (`deterministic_best` over
  `gen_best_objfun_val` and `gen_best_euclid_val_cross`) under euclid. The registration pass of
  solve3D_addon (`set_refine3D_register_policy` in `src/main/solve/simple_solve3D_controller.f90`) is
  one full-cohort `prob` iteration at a fixed `NSPACE_REGISTER` (1000) against frozen references: a
  table fill at a chosen grid, as a step.
- The fractional update and the balanced sampling (`src/main/ori/simple_oris_sampling.f90`;
  `doc/policies/importance_sampling_fractional_update_policy.md`): each iteration of refine3D and
  solve3D searches a sampled fraction of the particles (`update_frac`, `nsample`,
  `sample4update_cnt`), drawn by class or class average when `balance` asks for it
  (`sample_balanced_mode`); the rows of the sampled particles are what the iteration computes anyway.
  `eo_rebalance=yes` reassigns the selected particles' halves by direction and score: a
  view-balanced split of a dataset exists.
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
- The neighbourhood search: `refine3D refine=prob_neigh` (`simple_eul_prob_tab_neigh.f90`) searches
  a sparse neighbourhood of each particle's current orientation, geometric under
  `prob_neigh_mode=geom` with the cone `prob_athres`, with its own in-plane and shift search, its
  own posterior draw and the centring record it writes (`shincr`, `shsig`); the dense `refine=prob`
  is the search without a neighbourhood. The local refinement: `refine3D refine=cont` continues
  from a discrete seed, bounding the rotation
  by `athres_cont` and the shifts by `trs`, and reconstructs the halves' starting volumes itself; the
  stream's 3D stage runs one iteration of it from the seeds `set_refine_seeds` sets
  (`src/main/solve/simple_state_assignment.f90`: the halves of the partition, the nearest grid
  direction). Its recoverable full pose is the bar; direction agreement on a 6° grid alone is
  insufficient.
- The logistic core (`src/main/quality/simple_quality_logistic.f90`), binary and at most sixteen
  features, is not the student's shape; the learn/evaluate/promote lifecycle and the versioned key=value model file of
  `src/main/cavg_quality/simple_cavg_quality_model.f90` are the pattern.
- OpenBLAS is linked (`cmake/Dependencies.cmake`) and `sgemm` is called from Fortran
  (`src/utils/math/simple_linalg.f90`, flex); the OpenMP BLAS policy keeps the library on the calling
  thread inside a team (`doc/implementation_notes/completed/blas_openmp_policy.md`). Nothing in the
  repo does backpropagation: the network is new code.
- Private programs: `add_private_prg` in `src/main/ui/simple_private_prgs.f90`, the `case` in
  `production/simple_private_exec_driver.f90`, a thin commander (`assign_cavgs_states` in
  `src/main/commanders/simple/simple_commanders_cavgs.f90` is the model).
- Validation data (section 10): the Ribosembly benchmark on dell2 (`/data2/multistate_bench`: 80,000
  simulated particles, box 128 at 3 Å per pixel, SNR 0.01, true poses, shifts uniform to plus or
  minus 10 pixels, CTFs drawn per particle, truth maps; its 48,000 particles of the medium assembly
  are a single-state truth set, its 24,000 of the small one another); `/data/simulated_vcrnapol` on
  the Dell (two states with truth, 17,000 and 13,000 particles, each a single-state set); the real
  `/data/trpv1` and `/data2/trpv1_c2`, `/data/bgal`; the small membrane protein `/data/exp_gate_FlhB`
  and streptavidin `/data/10335-streptav`; two CESPED entries from Zenodo (section 2.4), homogeneous
  refinements by construction; the single-state high-level tests on simulated 1jxy and
  beta-galactosidase (`src/main/commanders/test/simple_commanders_test_highlevel_cont.f90`) as the
  pattern for a nightly route.

### 2.3 What XEMA's two-fold method contributes

XEMA (`~/src/XEMA`, Hans's cross-linking mass-spectrometry engine) motivates label-free partitions,
held-out choices, published agreement and replayable model records (the sources in section 15).
Here each student owns its teacher references, spectrum, model selection and confidence gate;
micrograph-level inner splits serve model development, a separate cohort calibrates the gate and an
untouched cohort audits it. Counts are compared as paired micrograph-clustered outcomes, not
independent Poisson observations. One route supplies each row, fixed-chunk reductions make a fit
reproducible across thread counts, and the model records every data-dependent choice and all
constants. Sections 7 and 8 define these boundaries.

### 2.4 What the CESPED checkout teaches

`~/src/cesped` (oxpig/cesped) holds the benchmark's loader, the Image2Sphere baseline, its training
configuration and its evaluation. It settles several things.

- **Held-out halves, and an evaluation through reconstruction.** An entry's half 0 trains (with a
  70/30 split inside it for validation and hyperparameters) and half 1 is judged, then the reverse.
  The evaluation reconstructs from the predicted poses and reports the FSC against the reference map
  beside the reference halves' own FSC, so that the method's loss is separated from the map's
  quality. That is the split-half protocol of section 10, with their halves random where ours are
  view-balanced, and tier 4 adopts the same two numbers.
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
  the error by a third. The teacher's soft row at the student's band (section 4) is the same idea
  carried through: the whole row, not one pose and one weight.
- **Translations were the stated gap.** Image2Sphere predicts no shift, the evaluation used the
  truth shifts, and joint rotation and translation estimation heads the authors' list of
  improvements. The shift map in the candidate search at the predicted directions (sections 5 and 8)
  is this design's answer.
- **CTF.** Their loader offers `ctf_multiply` beside `phase_flip`, and the paper used the flip; no
  one measured the difference, and neither keeps the CTF where flex keeps it, in the operator. The
  student takes flex's formalism (sections 2.5 and 3) and reads its effect by defocus bin.
- **Capacity.** Their encoder is an ImageNet-pretrained ResNet152 on 65,000 real-space pixels,
  trained for up to 400 epochs on four A100s; this design's network has about a million weights
  on 10,000 band-limited polar values with the rotation built in and the shift found by the map.
  The lower capacity is plausible on the lower-dimensional input and is a risk the wider networks of
  section 6 bound.
- **The head.** Their distribution over SO(3) is band-limited at lmax 8, so its peaks cannot be
  much sharper than 20° whatever the 3.75° readout grid; and it classifies the in-plane angle as a
  class dimension, which the authors name as not separable from the in-plane translation. The flat
  softmax over directions with the in-plane angle left to the search is the cheaper factoring.
- **The data.** Eight public entries with RELION poses, two membrane proteins in a micelle or a
  nanodisc among them (the GPCR complexes 10786 and 11120 are the hardest, at 33° and 45°), two
  symmetric ones (C2, D2), and published baselines, all homogeneous refinements. Tier 4 takes two
  CESPED entries as its real data beside the Dell's own, so that the numbers compare with
  Image2Sphere's and cryoPARES's on the same particles. No entry is on this box yet; they come from
  Zenodo, with a small TEST entry for smoke tests.

### 2.5 What the flex machinery teaches

`src/main/flex` (the `simple-main-flex` skill; `doc/algorithms/heterogeneity_analysis/flex_pca.md`)
models Fourier data with a learned low-rank basis, in flex_pca over projections and in the PPCA
kernels of `cls_expansion` and `trajectory_denoise` over class members and frames. It has solved,
in this code base, the questions a network over Fourier transforms raises.

- **Hermitian by the half plane, never by redundancy.** Flex stores and computes on the native
  k <= 0 half plane; the k = 0 line is its own Friedel mate and is counted once (`cov_herm_inner`).
  Every inner product is the real half-plane product `<u,v> = sum_j wq_j Re(conj(u_j) v_j)` with the
  quadrature weight `wq` 2 for a coefficient whose mate lies off the stored half and 1 for a
  self-mate, which equals the full-plane real product by Parseval, and its polar sampler packs each
  sample as `sqrt(wq)` times (Re, Im) into a real vector (`polar_sample_particle_fused`), so that a
  Euclidean dot product on the packed reals is the Hermitian product. Taps that fall on k > 0 read
  the conjugate of the mirrored sample: the Friedel mate per tap. The student takes this whole:
  the half ring as the data, the shell weight packed as its square root, the conjugate wrap at the
  ring's ends, and layers that commute with conjugation, so that the network is a model of the
  Fourier transform of a real image and the full ring is never materialised (sections 3 and 6).
- **The CTF formalism.** Flex's observation is the whitened, CTF-modulated central section,
  `y_i = a_i C_i S_i P_i V + n_i`, and the CTF never leaves the operator: the data enter through
  the cross term `b_q = <T u_q, y_i>` with `T = C S P`, which for a real CTF is the CTF-free basis
  correlated with the CTF-multiplied data, and the Gram `G_qr = <T u_q, T u_r>` carries the CTF
  squared, so that a latent `A^-1 b` is never biased by the envelope and nothing is deconvolved
  from an image. In the polar E-step the Gram factorises as `sum_k w_i(k) C_qr(k)` with `w_i(k)`
  the ring mean of the CTF squared and the azimuthal variation reported as `tazim`, the
  approximation's own diagnostic; the contrast `a_i` is fitted in closed form against the projected
  mean and clamped to [0.1, 5]. The student takes the formalism whole for its first layer (section
  3, item 2): the cross term of each filter with the CTF-multiplied particle under the shell
  weights, normalised by that filter's energy under the particle's CTF squared, which is the Gram's
  diagonal. No CTF channel, no flip, no Wiener filter.
- **Two fits on disjoint halves, and the agreement as the reproducibility.** Flex fits two bases
  on disjoint particle halves (the mod-4 rule of `flex_pca_half_of`), frame-aligns and merges
  them, trusts each component by its match cosine between the independent fits rather than by the
  requested rank, and feeds the cross-fit FSC back as the SSNR ridge of the joint solve. The
  student's two production halves and the diagnostic agreement of section 7 are the same
  discipline; whether an agreement statistic may regularise, as the ridge does, rather than only
  report, is not proposed here.
- **Coefficient-split cross-validation.** The PPCA kernels choose rank and bandwidth without
  holding out a single particle: the fit band is split into two interleaved sets of 4 by 4
  Fourier blocks in a checkerboard (`coeff_sets`), the model fitted on one set scores the other,
  and the pre-images use each set's weights from the other set's model, so there is no selection
  bias and no Stein term. It is noted as the way flex chooses what cannot be chosen on particles.
- **A noise spectrum from the data.** `flex_traj_noise_spectrum` takes the per-shell noise power
  from the median half-squared difference of consecutive frames; flex_pca takes `sigma2` from the
  high-frequency shells. The student keeps the calculator's sigma2 and its bootstrap, the same
  quantity the objective uses.
- **Data-free initialisation.** Flex starts its subspace from geometry, the lowest lattice points
  realised as cosine and sine pairs, so that two runs are comparable and the EM never seeds inside
  a noise subspace. The student's weights start from a fixed seed recorded in the model, and the
  teacher's rows, not a data moment, are the first thing read.
- **A plane cache with a validated contract.** Prepped Fourier planes are cached on disk under a
  run token with a provenance header (the project's path hash and mtime, the rows, box, crop,
  sampling, the selection), written last and renamed into place atomically, a quarter of the free
  space at most, fetched by copy, deleted once nothing after the embedding reads planes
  (`simple_flex_pca_plane_cache`, `simple_flex_pca_planes`). The student's input tensors are read
  once per epoch, width and band; they are cached under the same contract (section 9).
- **The OpenMP discipline.** Per-thread scratch sized by `omp_get_max_threads()` at entry and
  indexed by the thread number; only the E-step body inside the region, bank builds and inserts
  outside; no FFTW plan in a region; the BLAS probed for the calling thread (`probe_blas`);
  reductions in chunks. The network's trainer follows the same rules (section 6).
- **Gates that are declared and that stand down.** The two merge gates are a disattenuated
  FSC-type ratio at 0.98 under complete linkage and a chi-square view test at p = 0.001 with a
  robust effect-size cut, and the view gate stands down when more than a third of the states are
  flagged. The student's gates are declared the same way and decline the same way.
- **One project gateway.** Every project read and write of flex goes through
  `simple_flex_pca_project_gateway`; workers never write the project; a cached embedding is reused
  only when particle identity and model provenance match. The student's `pconf` and route writes
  and its model's refusal of a changed project follow it.

## 3. The student's input

`simple_polar_feats` (`src/main/learn/`): one routine builds the input of a particle from the polar
calculator, called by training and prediction alike so that the two never drift apart.

1. **The band.** Rings `kfromto(1)` to `box*smpd/lp` with `lp` the student's low-pass. The band is one
   of the two hyperparameters (section 6), 15 Å first, then 10 and 8 Å; the number of rings depends on
   `box*smpd/lp` alone, about 25 at 15 Å for the benchmark's box 128 at 3 Å per pixel.
2. **The CTF, in flex's formalism.** The tensor is built in the angular-sample domain from
   `pfts_ptcls(:,k,i)*ctfmats(:,k,i)`, the CTF-multiplied particle, because that is what the cross
   term of a CTF-free filter with the observation is: the search's `Q_k = REF CTF conj(X)` and
   flex's `<T u, y>` alike correlate the CTF-free model with the CTF-multiplied data. The CTF stays
   in the operator: the first layer's filters are CTF-free, each filter's response is normalised by
   its energy under the particle's CTF squared and the shell weights, the Gram's diagonal, so that
   the envelope `CTF^2` of the product never reaches the deeper layers as structure, and a shell
   near a zero carries no signal and enters at no weight. The normalisation uses the ring mean of
   the CTF squared, as flex's polar E-step factorises its Gram, with the azimuthal variation
   reported as flex reports `tazim` (decision 14 keeps the full azimuthal CTF squared as the
   alternative). The calculator's `ft_ptcl_ctf` is the angular FFT of this Friedel-extended ring:
   its index is the harmonic order, not the angle, so it is the shift map's input (section 5) and
   never an angle-axis CNN's. Nothing is restored: no phase flip, Wiener filter or deconvolution in
   any design, as nothing is in flex; the output is a class, not an image. The per-particle contrast
   flex fits against its model is here the band's scale by the particle's weighted power (item 3).
3. **The shell weights.** The shells are weighted as the search weights them under objfun=euclid,
   the objective of the biomolecular workflows this is for: `k/sigma2(k)`, where sigma2 is the
   global spectrum of the particle's half, never a per-particle one: the calculator's weights are two
   spectra, one per half, and a part's per-particle residual records exist only to be summed into
   them. An unseen particle therefore needs nothing beyond its half, assigned at import, and its own
   weighted power (`wsqsum_ptcl_band`), which scales the band. The weights belong to the model, not
   to the project: the student's record holds the weights it was trained with, and prediction applies
   them, never the project's current sigma2, so a later refinement that moves the spectrum does not
   move the student's input distribution and a retrain picks the new spectrum up. A project without
   a spectrum is given one by the sigma2 bootstrap (`calc_pspec`) before training or prediction,
   never scored under another objective. Tensor amplitudes use `sqrt(w_k)` so that a squared norm
   applies the shell weight once. The logarithm of the weighted power is one scalar side input whose
   usefulness for the confidence is measured; it does not authorise a particle-rejection decision.
4. **The half ring, with the conjugate wrap.** The calculator's `pftsz` angles over the half turn
   are the data; the other half is their Friedel mate and is never materialised. The angular axis
   is the calculator's own rotation grid, so an in-plane rotation by one grid step is one cyclic
   shift, and what passes the ring's end returns conjugated: the convolution's padding at both ends
   is the conjugate of the opposite end, as flex's polar sampler takes the Friedel mate per tap.
   Rotation equivariance is then exact on the search's rotation grid, with no resampling and no
   interpolation (a resampled grid only if the cost demands it).
5. **The tensor.** Two channels by the rings by the half ring's `pftsz` angles: the real and
   imaginary parts of the CTF-multiplied particle, each scaled by `sqrt(w_k)` so that a Euclidean
   product of the packed reals is the objective's Hermitian product on the full ring, as flex packs
   its polar samples by `sqrt(wq)`. The particle's CTF squared, ring-averaged, travels beside the
   tensor for the first layer's normalisation, not as a channel. About 8,500 values at 15 Å for
   the benchmark's box. The phases are
   kept throughout; nothing is reduced to a power. Augmentation lives in the same module: the
   calculator's own phase ramp for shifts (`memoize_shifted_ptcl`), cyclic shifts with the
   conjugate wrap for rotations on the grid.

## 4. The direction grid and the teacher

- **The grid.** `oris%spiral` at `NSPACE_STUDENT` (1000) directions, about 6° apart, over the whole
  sphere in c1 and the asymmetric unit as `build%eulspace` is built for a symmetric point group. The
  in-plane angle is not a class: the student is invariant to it by construction and the search at
  the predicted directions recovers it (section 8). The student's output is a distribution over the
  directions modulo the point group; no rotation or shift is regressed, the candidate search supplies
  the in-plane angle and the shift, and Euler angles stay the project's interchange.
- **The references.** The teacher scores a particle against the reconstruction of the particles it
  is trained with: in the split-half protocol the training half's map, in the stream and the batch
  consumers the particle's own half's map as the matcher's parity references are prepared
  (`prepare_refs_sigmas_and_pftc`; refine3D reconstructs them from the project's halves at its first
  iteration). The model records the references' provenance (the project, the rows, the band, a hash
  of each volume).
- **The teacher is the posterior of the probabilistic machinery** (Hans, decision 4): the table's
  tempered posterior rows at the student's grid and band, as `fill_tab_range` computes them against
  the training half's map, each reference scored at its drawn in-plane rotation and the best
  candidates polished in shift through the shift map, with no assignment draw taken. In the loop of
  refine3D and solve3D these rows exist for the sampled particles anyway. If the in-plane draw
  proves to leave too much noise in the labels, the same rows with the deterministic best over the
  in-plane grid (`deterministic_best`, the diagnostic posterior's path) are the alternative, at the
  cost of a fill.
- **The soft label** is the row's posterior `t_c = exp[-tau_i (L_c-L_min)] / sum_j exp[-tau_i
  (L_j-L_min)]` with `tau_i` the recorded `prob_beta` times the particle's temperature
  (`posterior_temperature`), the factor the table draws at; objective, temperature and shell
  weights are recorded in the model.
- **What the earlier search contributes.** The map of the training particles, the halves, which
  particles have a converged pose (the candidate training rows), and the initial shift support. The
  teacher and the prediction use the same image preparation and initial shift support; converged
  shifts or centring widths that a new particle will not have cannot restrict the training searches.
- **The targets retain uncertainty.** The loss uses the full soft rows in bounded batches. There is
  no uniform label smoothing, no fixed top-k renormalisation and no inverse-frequency class weighting
  in the probability head: with 1,000 classes, 5% of uniform smoothing would cap the mass in 64
  classes at about 95% even for a perfect fit, so the 99% coverage of section 8 could never be met;
  and reweighting the classes changes the meaning of the predicted probabilities. A later sparse
  encoding must preserve the omitted mass and its distribution, or bound its approximation error; a
  broad row stays broad. A class-balanced draw of the training rows may improve coverage, but a
  nonuniform draw is corrected by its known inclusion weights for the intended particle population.
  Rare-view performance is reported separately rather than changing the probabilities' meaning.
- **Discrete grid, continuous pose** (Hans, decision 2). The teacher's rotations and shifts are
  continuous: the polish is continuous in shift and the local refinement in both. The student's
  output is discrete because the posterior over directions is multimodal at low SNR and the
  references are discrete: a softmax over a 6° grid holds competing directions with their masses,
  which a regressed rotation cannot. The grid is a sampling of a continuous posterior on the
  sphere: a continuous direction is the posterior-weighted mean of the directions in the mode's
  cone, with the cone's spread as its uncertainty, and the continuous in-plane rotation and shift
  come from the neighbourhood search, which has the references and is the better continuous
  estimator. Learning a continuous density on SO(3) directly, as Image2Sphere does with spherical
  harmonics at lmax 8, gives a smoother but blurrier posterior than the grid at 1,000 directions.
- **Training and evaluation targets have different roles.** Training may use the band-specific soft
  scores. Every design is evaluated against one frozen full-search reference pose (the standard
  machinery's, at its band) and the actual downstream refinement, with the truth where it exists. A
  coarser band does not win by agreeing with its own easier teacher.
- **The teacher's resolvability is read before anything is trained.** The deterministic rows give,
  per row and band, `pmax`, `neff` and `k99`; their distributions over the training particles say
  whether the band carries a concentrated assignment under that score model, before any training.
  These are conditional on the references, the objective and the temperature, not an
  information-theoretic bound. A band at which the median `k99` spans a ring of directions is
  reported as such, and the student's claim for that dataset is assessed against the common reference
  and the truth (section 10).
- **The consumer's grid is decoupled** from the student's by the local search: the student's
  directions seed an in-plane and shift search and then `refine=cont` within `athres_cont`, and
  refine3D's own grid (`set_projs` picks the nearest direction in the job's grid) can be 2,500 or
  5,000 as the stage likes.
- **Cost.** The rows cost what a table fill costs at the student grid and band, nothing when they
  are the iteration's own; the training-set size study of section 10 says how many particles they
  are needed for.

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
- The map is linear in the reference, so a reference may be a reprojection or a class average in the
  band.
- The CTF enters as the search applies it, in the cross spectrum, per particle and exact.

The shift is handled in two places, and in neither by the network's architecture.

- **Training: tolerance by augmentation.** The network of section 6 is trained with shifts drawn
  uniformly in the disk through the calculator's shift matrices (`memoize_shifted_ptcl`), the exact
  phase ramps, and with rotations as cyclic shifts: it learns to put a shifted particle in the right
  class. Its accuracy against the shift magnitude is a declared read (section 10).
- **Prediction: refine3D's neighbourhood search** (Hans, decision 7). `predict_poses` writes the
  student's direction as the particle's orientation and the width of its posterior as the cone;
  refine3D under `refine=prob_neigh prob_neigh_mode=geom` searches that cone with its own in-plane
  and shift search, the map and the polish, from the particle's prepared shift, and draws its pose
  from its own posterior. The student runs no search of its own.

**The centring record is the search's** (Hans, decision 13): the neighbourhood search writes the
shift, `shincr` and `shsig` from its own shift posterior, as it does for every particle it searches;
the student writes no shift.

**Cost.** The neighbourhood search over a cone of a few directions costs the map and the polish
for those directions, ten to a hundred times below the dense search over thousands; it is measured
in section 10, never assumed.

## 6. The hyperparameters and the network

**Two hyperparameters.** The network's width (16, 32 or 64 channels in the first block, doubling
per block) and the band (15, 10 or 8 Å). The cheapest, 16 channels at 15 Å, is tried first; a wider
network or a finer band replaces it only when held-out micrographs say it recovers clearly more
seeds (section 7). A small membrane protein in a micelle is, at 15 Å,
dominated by the micelle, whose outline gives the tilt of the membrane normal and little of the
azimuth about it; the azimuth may sit at 8 to 12 Å, in the transmembrane helices, which is why the
search's late stages run there, and why the band is a hyperparameter.

**The network** (`simple_circ_cnn`, `src/main/learn/`): an encapsulated class with private components,
`new`/`kill`, forward, backward, one Adam step, read and write.

- **Hermitian layers on the half ring.** The network is a model of the Fourier transform of a
  real image, so every layer commutes with conjugation, and the half ring with the conjugate wrap
  computes exactly what the full ring would. Three blocks of (a 3 by 3 convolution along the ring
  and the angle with real kernels applied alike to the real and imaginary channels, the complex
  signal's linear map with real weights, which is the only linear map that commutes with
  conjugation; conjugate-wrap padding along the angle and zero padding along the rings; a bias on
  the modulus; a modulus-gated nonlinearity `z * relu(|z| + b) / |z|`, which mixes the real and
  imaginary parts only through the conjugation-invariant modulus; average pooling along the rings
  only), with 16, 32 and 64 channels and angular stride one throughout. The first block's filters are
  CTF-free and their responses are divided by the filter's energy under the particle's
  ring-averaged CTF squared and the shell weights before the modulus gate, flex's Gram diagonal
  (section 3, item 2); the deeper blocks see structure, not the envelope. The invariant pooling is the angular mean over the full ring, which for a Hermitian
  signal is the mean of the real part over the half ring, the imaginary part cancelling
  identically, beside the mean of the modulus; both are real. Then a dense layer to 256 with ReLU,
  the side scalar appended, and a dense layer to the classes with a softmax. Rotation
  equivariance is exact on the calculator's rotation grid, and the invariance exact after the
  pooling. About half a million weights at 1,000 classes; the per-particle cost is counted on the
  retained layout before it is quoted.
- The convolution as im2col plus `sgemm`, the backward pass as the two transposed products; weights,
  activations and gradients `real(sp)` and allocatable; one Adam state per weight array.
- **Bit-identical on any thread count.** The mini-batch is cut into fixed chunks; the chunks are
  computed in parallel, each into its own gradient buffer indexed by chunk, not by thread, with
  per-thread activation buffers sized by `omp_get_max_threads()` read before the region and the BLAS
  on the calling thread inside it; the chunks are reduced serially in chunk order after the region.
  The same sum at one thread and at any team, bit for bit, and the tester asks for identity, not a
  tolerance.
- The loss: the cross-entropy against the complete soft targets; Adam at 1e-3, batches of 256,
  weight decay 1e-5; at most 30 epochs, with early stopping on held-out micrographs of the training
  set.
- **The model file.** A key=value text header (`model_version=1`) and the weights as a binary file
  beside it (two when the stream's two students share a file, decision 7). The header records the
  design and the band, the grid with `nspace` and `pgrp`, the angular grid (`pftsz`, or the resampled count), the frozen weights, the reference
  and spectrum hashes and source cohorts, the teacher's definition and temperatures, the held-out micrographs, the confidence threshold, the timings and the
  reason a model was kept or declined; a reader refuses incompatible grids, objectives or
  references.

## 7. Held-out micrographs, the choice of the design and the confidence

- **The halves.** In the split-half protocol of section 10 they are the two view-balanced halves of
  the dataset as `eo_rebalance` assigns them by direction and score, one to train on and one to
  test on. In production each half of the even-odd partition would have its own student, trained
  and applied within it; whether that is required, or the cross arrangement of the protocol is
  admissible for the gold-standard FSC, is decided later on tier 4's evidence (section 12).
- **Held-out micrographs.** Whole micrographs of the training half are held out before any fit,
  for early stopping, for the choice of width and band, and for the confidence threshold. The test
  half is never touched until the model is frozen; in the stream the next publication's new rows
  play that part.
- **The choice.** The cheapest design is kept unless a wider or finer one recovers clearly more
  seeds on the held-out micrographs: a seed is recovered when the neighbourhood search from the
  student's direction reaches the pose the full search gives the same particle, and the comparison
  is the paired difference over micrographs, resampled by micrograph.
- **The confidence** is the student's posterior mass inside the cone the neighbourhood search will
  cover. Its threshold is set on the held-out micrographs so that, above it, the full search's
  direction lies inside the cone in a declared fraction of the particles; below it the particle
  takes the dense search. The threshold and the fraction are recorded in the model.
- **Agreement is diagnostic.** Where two halves have students, their agreement on untouched
  particles is reported and never used to choose or to gate.
- **Decline.** A student that does not recover more seeds than the incumbent route on the held-out
  micrographs, or whose test-half result fails section 10's reads, is not used; the search runs as
  it would have.

## 8. The two programs

- **`learn_pose_student`**, a private program: the strategy routine `run_learn_pose_student(cline)`
  on the caller's `(params, build)` with a thin commander. Inputs: `projfile`, `lp`, `nspace` (the
  student grid), `mskdiam`, `pgrp`, `fname`, `nthr`; the references are the training particles' map
  (`vol1` in the split-half protocol, the project's half reconstructions in production). Steps:
  prepare the references and the calculator as refine3D's matcher does; hold out micrographs; fill
  the rows of the training particles at the student grid and band, or take the iteration's own;
  build the features; train, choose width and band, set the confidence threshold on the held-out
  micrographs; write the model and its record.
- **`predict_poses`**, a private program: `run_predict_poses(cline)` with a thin commander. Inputs:
  `projfile`, `infile` (the model), `mskdiam`, `nthr`. Per selected particle: the features under the
  model's band and weights; the forward pass; the mode direction, the cone holding 99% of the
  posterior, at most `KMAX_STUDENT` (64) directions, and the confidence; the mode written as the
  particle's orientation with the cone's width, the confidence to `pconf` and the route to its key,
  both numbered ori slots. A particle whose cone exceeds the limit or whose confidence is below the
  threshold is marked for the dense search. The program writes a part file the caller merges, never
  the caller's project; refine3D does the rest.
- **One source per particle.** Every particle's pose comes from the search that ran for it, the
  neighbourhood one seeded by the student or the dense one; the route is recorded per particle.

## 9. Where it fits

- `src/main/learn/`, picked up by the recursive glob of `src/CMakeLists.txt` with no CMake change:
  `learn.inf` (one `!@descr:` line), `simple_polar_feats.f90`, `simple_circ_cnn.f90`,
  `simple_pose_student.f90` (the model: the student, the grid, the threshold, the file I/O), the
  input tensors' disk cache inside `simple_polar_feats` under the plane cache's contract (a
  provenance header, a run token, written last and renamed into place, a quarter of the free space
  at most, deleted after training), their testers, and a `CLAUDE.md` importing a
  `simple-main-learn` skill proposed with the delivery and applied on approval
  (`scripts/check_guide_imports.py` checks the import).
- The search layer: a getter of a particle's posterior row from the table, beside
  `posterior_row_stats`; and the neighbourhood search's cone taken per particle from the student's
  width where today it is the global `prob_athres`.
- The strategy routines in `src/main/strategies/parallelization/simple_strategy_pose_student.f90`
  (shared memory: one job on the caller's threads), the commanders in
  `src/main/commanders/simple/simple_commanders_pose_student.f90`, registered with `add_private_prg`
  and the driver's `case`; the required keys of a private program as
  `.github/skills/simple-main-ui/references/cli_change_chain.md` lays out.
- **Consumer 1, batch.** The split-half protocol as it stands, and refine3D's neighbourhood search
  seeded by the student for particles never searched: the registration pass of solve3D_addon.
- **Consumer 2, in the loop of refine3D and solve3D.** Each iteration computes the rows of its
  sampled fraction anyway. If a few hundred to a few thousand of them train a useful student, the
  iteration can seed the neighbourhood search of the particles it did not sample and end with every
  particle aligned instead of the fraction, at the cost of a neighbourhood search instead of none.
  The training-set size study of section 10 decides it; the integration is a step of the
  distributed master beside the sample draw.
- **Consumer 3, the stream**, after the multi-state extension of section 15.
- Documents updated with the delivery: a page `doc/algorithms/pose_student.md`; the sampling and
  fractional-update policy for consumer 2; this note moves to `completed/` with its measurements.

## 10. Validation

**The split-half protocol**, the real-data test, on every dataset of tier 4 and, with the truth
beside it, on the sets of tiers 2 and 3:

1. The dataset is split into two view-balanced halves as `eo_rebalance` splits them, by direction
   and score from the standard machinery's alignment of the whole set.
2. The standard machinery reconstructs from the training half alone: its map, its spectrum and its
   rows are the teacher's.
3. `learn_pose_student` trains on the training half, with micrographs held out inside it, and
   `predict_poses` seeds the test half, which no map, spectrum or fit has seen; refine3D's
   neighbourhood search aligns the seeded particles and the dense search the rest; one iteration
   of the local refinement follows and the test half is reconstructed.
4. The reads: the fraction of the test half seeded; among the seeded, the fraction whose pose
   agrees with the standard machinery's independent alignment of the same particle within the
   local refinement's capture range; the cone size per particle; the FSC of the test half's map
   against the training half's map beside the FSC of the two halves under the standard machinery,
   CESPED's two numbers; and the time of every step against aligning the test half the standard
   way.
5. The halves are swapped and the protocol repeated.

**The training-set size study** (Hans, decision 12), on the same split: the student is trained on
class-balanced subsets of the training half of 200, 500, 1,000, 2,000, 5,000, 20,000 particles
and the whole half, as `sample_balanced_mode` draws them, and every read above is reported against
the subset size with the cycle time of rows, training, prediction, neighbourhood search and
refinement. The curve says how few images make a useful student and whether consumer 2 is feasible
at an iteration's sampled fraction, since the soft rows teach many directions per particle and the
Hermitian network shares its weights over the whole ring. The same with the iteration's own drawn
rows as the teacher beside recomputed deterministic rows (decision 4).

The tiers:

1. **Unit** (the fast gate, a `pose_student` suite registered in a suite table of
   `simple_commanders_test_class.f90`): a central-difference gradient check of every layer on a
   tiny network; a toy problem (projections of a Gaussian-blob volume at SNR 1 over 20 directions,
   accuracy at least 0.95 on held-out images); the same model bit for bit at one thread and at a
   team of 4; a model-file round trip exact; the shift augmentation checked before the angular
   resampling; the Hermitian contract: every layer's output for the conjugated input is the
   conjugate of the output, the half-ring computation with the conjugate wrap equals an explicit
   full-ring computation bit for bit, a rotation by one grid step of the input rotates every feature
   map by one step and leaves the pooled features unchanged, and the pooled features are real; the
   packed Euclidean product equals the objective's Hermitian product; the first layer's normalised
   response to a filter's own CTF-modulated template is one at every defocus, and its response to
   the search's cross term matches the calculator's correlation; broad target rows keeping their
   tails and sending the particle to the dense search when 99% needs more than 64 directions; the
   feature cache's contract refusing a changed project, selection or geometry; fixtures for the
   grid's symmetry handling (group-equivalent directions, distinct antipodal views).
2. **Truth, the Ribosembly medium assembly (48,000 particles), then the small one (24,000).** The
   split-half protocol with the truth beside it: the teacher's resolvability per band and its error
   against the truth first; then, against the truth poses, the seeded fraction, the cone sizes,
   whether the truth lies in the cone, the recovered-seed fraction before and after the local
   refinement, per view and per defocus bin against the teacher in the same bin, the training-set
   size curve, and after refinement the FSC at 12 Å against the truth map beside the truth-pose
   reconstruction's; the same with the truth poses as the teacher, to separate the student's error
   from the teacher's.
3. **Truth, easy:** `/data/simulated_vcrnapol`, one state's particles, the same measurements.
4. **Real:** `/data/bgal`, `/data/trpv1`; the hard cases where the band carries the least,
   `/data/exp_gate_FlhB` and `/data/10335-streptav`; and two CESPED entries, one soluble (10648,
   PKM2, D2) and one membrane protein (11120, the TSHR complex, the benchmark's hardest), so that
   the numbers compare with Image2Sphere's and cryoPARES's on the same particles. The split-half
   protocol on each, the training-set size curve on at least bgal and one CESPED entry, every band
   against the same full-search alignment.
5. **Cost:** on the same runs, the time of every step on matched hardware against the standard
   machinery's alignment of the test half, and for consumer 2 the cycle time at each subset size
   against one iteration of refine3D on the same particles.
6. **Nightly:** a `pose_student` route on the simulated single-state sets of the high-level tests
   (1jxy, beta-galactosidase): the split-half protocol in miniature, within the existing entry and
   area suites, with the actual runtime recorded as `doc/policies/test_environment_policy.md`
   requires.

## 11. Order of work

- Phase 0, the feature module and the rows: `simple_polar_feats` with its tester, the row getter,
  and the teacher's rows on the Ribosembly medium assembly at each band, with their resolvability
  and their error against the truth.
- Phase 1, `simple_circ_cnn` with its tester in the fast gate (tier 1), then the 16-channel network
  at 15 Å on the medium assembly under the split-half protocol, and the first points of the
  training-set size curve from 200 particles.
- Phase 2, `simple_pose_student`, `learn_pose_student`, `predict_poses` and the per-particle cone of
  the neighbourhood search; the Ribosembly and vcrnapol measurements (tiers 2 and 3).
- Phase 3, the real-data split-half runs and the training-set size curves (tiers 4 and 5), the
  small particles and the CESPED entries included.
- Phase 4, the nightly route (tier 6), the algorithms page, the skill proposal; the decision on
  consumer 2 from the measured cycle.
- Phase 5, later: consumer 2 in the distributed master, then the multi-state extension and the
  stream (section 15).

Each phase ends with the fast gate, the touched unit suites and `git status`; the development stays
uncommitted for review as one diff.

## 12. Decisions

Ruled by Hans on 2026-11-06:

1. The folder and the names: `src/main/learn/`, `simple_polar_feats`, `simple_circ_cnn`,
   `simple_pose_student`; the programs `learn_pose_student` and `predict_poses`.
2. The grid stays discrete for the direction; the continuous pose comes from the search and the
   local refinement (section 4).
3. The width and the band are the hyperparameters, tried cheapest first (section 6).
4. The labels are the posterior rows of the probabilistic machinery; the deterministic in-plane best
   is the alternative if the drawn rows prove noisy (section 4).
5. The confidence's ori key `pconf` and the route key, written to `ptcl3D`.
6. The model file: one text header with the weights as a binary file beside it.
7. The student hands its direction and cone to refine3D's neighbourhood search; it runs no search
   of its own (sections 5 and 8).
8. A consumer retrains when its references change, if the training-set size study makes that
   affordable; no fixed cadence.
9. The shell weights are the teacher's: the spectrum the rows were whitened with, frozen in the
   model; the per-particle residual spectra of the refinement are summed into it, and a new particle
   has none until it has a pose.
12. No fitting floor is set; the training-set size study starts at 200 particles and decides
    whether the student can serve the stages and iterations of solve3D.
13. The shift and the centring record come from the search's own posterior; the student writes no
    shift.
14. The CTF as in flex (sections 2.5 and 3).
15. `KMAX_STUDENT` (64): a cone needing more directions sends the particle to the dense search.

Not to be decided now:

10. The split of the split-half protocol: `eo_rebalance` by direction and score is used; a
    class-stratified draw is the alternative if the halves prove unbalanced.
11. Same-half production or the cross arrangement for the gold-standard FSC: after tier 4.

## 13. Risks and limits

- **Low SNR.** The student cannot create missing pose information; a cone that does not hold the
  true direction sends the neighbourhood search astray, which is why the dense search stays behind
  it and the cone's coverage is read on the test half.
- **Small particles and small membrane proteins.** At 15 Å a micelle gives the tilt and not the
  azimuth; the finer bands exist for that, and the small-particle datasets and the CESPED membrane
  protein of tier 4 say whether a useful seeded fraction remains.
- **Preferred orientation.** Rare views are reported on their own, since a high overall seeded
  fraction can conceal missing directions; the view-balanced split and the class-balanced training
  draw are the protocol's answer.
- **Defocus as a confounder.** The first layer's normalisation removes the CTF's envelope as flex's
  Gram does, but the information lost at the zeros still differs by defocus, and a training set
  that confounds defocus with the view (views that cluster by acquisition) can still teach it.
  The excess error and the calibration are compared within defocus bins against the teacher in the
  same bin, and later acquisition blocks are audited.
- **Drift.** A student is tied to the map and the spectrum it was trained against; it is retrained
  when they change, never judged on its own seeds.
- **The in-loop consumer.** Its teacher is the iteration's own rows, its map changes every
  iteration, and the next iteration's search is its check; whether the cycle is fast enough on a few
  hundred particles is the training-set size study's question, answered before anything is
  integrated.
- **Cost on the Dell** (24 cores, no GPU): every timing is an estimate until tier 5, and a faster
  replacement for a dense particle search need not accelerate a whole stage.
- **Compile time.** The three modules are small; no new CMake target; heavy arrays allocatable.

## 14. Rules the work follows

- `doc/policies/layering_and_project_ownership_policy.md`: the strategy routines on `(params, build)`
  with `run_<x>(cline)` entries and thin commanders; the programs read the project and write part
  files, the master merges and writes its project once; the search for an existing step before
  writing one (the matcher's reference preparation, the deterministic objective kernels, the shift
  map, the seed contract of `simple_state_assignment`, the sample draw of the distributed master).
- `doc/policies/compile_time_policy.md` and the `simple-modern-fortran` skill: encapsulated classes
  with `new`/`kill`, heavy arrays allocatable, `use ..., only:` imports, the OpenMP rules, testers
  registered in a suite table (a tester that needs private entities is a submodule of the module it
  tests), no new CMake target, no per-file flags.
- No new command-line key and no internal switch: the model travels through `infile` and `fname`,
  the thresholds, coordinate conventions and constants that decided live in the model file, with
  versioned policy defaults and explicit provenance. Required reference and calibration metadata is
  never replaced by an undocumented key or a silent fallback.
- The logistic core and the quality model are not changed; the `learn`, `evaluate` and `analyze`
  outputs of `model_cavgs_rejection` stay byte-identical.
- Nothing is committed; the whole development is one reviewable diff.

## 15. Outlook and sources

Once the pose student is validated, the grid is extended by the state, the teacher's rows run over
the states' references, and the stream's multi-state publications become the third consumer; the
heterogeneous Ribosembly and vcrnapol sets, whose single-state subsets validate the poses here, are
the truth for that step.

- cryoPARES: Sanchez-Garcia et al., bioRxiv 10.1101/2025.03.04.641536 (2025); the CESPED benchmark,
  arXiv 2311.06194 and its package `~/src/cesped` (`cesped/network/image2sphere.py`,
  `cesped/datamanager/relionStarDataset.py`, `cesped/datamanager/ctf.py`, `cesped/evaluateEntry.py`,
  `cesped/configs/`); cryoFIRE, arXiv 2210.07387 (NeurIPS 2022); semi-amortized pose inference,
  Shekarforoush et al., NeurIPS 2024.
- [CESPED evaluation protocol, section 2.5](https://arxiv.org/html/2311.06194v5#S2.SS5): cross-half
  supervised benchmarking and its map-to-reference comparison, the two numbers tier 4 reports.
- Related architectures: [CyCNN](https://arxiv.org/abs/2007.10588); [Polar Transformer
  Networks](https://arxiv.org/pdf/1709.01889); [Polar Transformers for cryo-EM image
  sets](https://arxiv.org/abs/2506.11283) and its [NeurIPS
  listing](https://neurips.cc/virtual/2026/oral/156032); [Equivariant amortized inference of poses
  for cryo-EM](https://arxiv.org/abs/2406.01630); [Rotationally invariant image representation for
  viewing direction classification in cryo-EM](https://arxiv.org/abs/1309.7643); [Probabilistic
  PolarGMM](https://arxiv.org/abs/2206.12959); [Deep Complex Networks](https://arxiv.org/abs/1705.09792);
  [Harmonic Networks](https://arxiv.org/abs/1612.04642).
- XEMA: `docs/algorithms/13_cross_validation.md`, `docs/algorithms/11_consensus_scoring.md`,
  `docs/design_notes/gold_standard_selection.md`, `docs/design_notes/small_dataset_scoring_model.md`,
  `src/consensus_scoring/xema_consensus_twofold.f90` and `xema_consensus_fitter.f90`.
- Flex: `doc/algorithms/heterogeneity_analysis/flex_pca.md`,
  `doc/implementation_notes/completed/trajectory_denoise_flex.md`,
  `doc/implementation_notes/completed/cls_expansion_kernel_preimages.md`; `simple_flex_pca_basis.f90`
  (`cov_herm_inner`), `simple_flex_pca_polar.f90` (`polar_sample_particle_fused`),
  `simple_flex_traj_denoise.f90` (`coeff_sets`, `flex_traj_noise_spectrum`),
  `simple_flex_pca_plane_cache.f90`, `simple_flex_pca_stages.f90` (`flex_pca_half_of`).
- SIMPLE: `doc/implementation_notes/completed/polar_shift_search.md` and
  `polar_shift_search_dynamic_centring.md`; `doc/implementation_notes/planned/polar_shift_laplace_width.md`;
  `doc/policies/importance_sampling_fractional_update_policy.md`;
  `doc/policies/stream/stream_3D_ingestion_policy.md`; the sources read for the teacher's contract:
  `simple_eul_prob_tab.f90` (`score_ref`, `deterministic_best`), `simple_polarft_corr.f90`
  (`gen_posterior_val`, `gen_best_objfun_val`), `simple_oris_sampling.f90` (`sample4update_cnt`,
  `sample_balanced_mode`) and `simple_strategy3D_cont.f90` (the fixed-state search).
