# Learned pose assignment: a student distilled from the probabilistic search

Proposal of 2026-11-05, revised on 2026-11-06 with Hans and after a methodology review. SIMPLE's
probabilistic search supplies per-dataset soft labels over the 3D direction, and a small model
trained on them proposes directions for particles the search has not seen, so that a local search
and the local refinement can finish the alignment. The first delivery is single-state: a pose
student, validated for poses before anything else. Whether it saves time is a measurement. Nothing is
built. Sections 1 and 2 are the idea and what exists, sections 3 to 9 the design, section 10 the
validation, section 11 the order of work, section 12 the open decisions for Hans, sections 13 to 15
the risks, the rules, the outlook and the sources.

**Changed after the run, 2026-11-09 (in the maintainer's tree):** the student's classes are the 500
mirror pairs of the grid, not its 1,000 directions, and HHCNN is invariant to the mirror; the mirror
is resolved by the search through the Euler relation (`euler_mirror`), as the opt-in mirrored
reference extraction (`mirr_proj=yes`) pairs the spiral. Sections 4 and 6 and decision 16 describe it;
the Progress entry "Mirror pairs (2026-11-09)" holds the change and its measurements. The phases'
entries keep the design of their day.

## 1. What is proposed

A pose student: a small network, trained per dataset on the probabilistic posterior of the dense 3D
table over a fixed grid of directions. The network is HHCNN, the Harmonic Hermitian convolutional neural network: Hermitian because every layer commutes with conjugation on the half ring with the conjugate wrap, so the network is a model of the Fourier transform of a real image; harmonic because the readout's class templates live in the angular harmonic domain of the polar transform and the rotation search is the cross-correlation over the cyclic shifts through the angular Fourier transform; convolutional on the polar Fourier representation (the name is Hans's ruling R5 of
2026-11-08, applied in phase 5). It reads a particle's CTF-multiplied polar Fourier transform
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
  the predicted directions recovers it (section 8). No rotation or shift is regressed, the candidate
  search supplies the in-plane angle and the shift, and Euler angles stay the project's interchange.
- **The classes are the mirror pairs** (Hans, 2026-11-09; decision 16). The spiral of
  `build_refspiral` holds the unmirrored directions of the asymmetric unit in its first half and their
  mirrors in its second, paired with a frame angle (`oris%get_mirror_pairs`, the pairing of
  `mirr_proj=yes`). A direction and its mirror are related by the Euler angles (`euler_mirror`), and
  the projection of one is the reflected image of the other. The student classifies the 500 pairs;
  its output, a distribution over pairs, is placed on each pair's unmirrored direction of the grid,
  and the search tells the two members apart with a cross-correlation (the neighbourhood search scores
  the cone about the direction and about its mirror). The teacher still scores all 1,000 directions;
  its rows are folded onto the pairs before the loss. The fold of ruling R1, below, is now the class
  space itself rather than a step after the prediction.
- **The folded sphere** (Hans, ruling R1, 2026-11-07, after phase 2's "What the cones are"). The
  student is Friedel-symmetric by construction: a view and its mirror view (the antipodal direction,
  whose projection is the reflected image; the spiral's mirror partner) are one class to it, and the
  network is not changed to separate them. The posterior is folded over the mirror pairs before the
  cone is read: the mass of a direction and of its partner are summed, the mode and the 99% set are
  taken on the folded sphere, and the cone's radius is the set's spread about that mode (a member's
  distance the smaller of its own and its partner's). The two possibilities are resolved by the
  search (section 8). Ruling R2 (the same day) asks whether the student is, in effect, a network of
  spectral powers without phases, which would make it shift-invariant too: `predict_poses` reports a
  shift check on its first particles (the posterior and the pooled features of each particle shifted
  through the calculator's phase ramps against its own), and the training's shift augmentation stands
  or falls by it (section 5; the finding in Progress).
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
  class. Its accuracy against the shift magnitude is a declared read (section 10). As built, the disk
  is three times `shift_sigma` in pixels (`SHIFT_AUG_NSIG`), the keyed draw per row and epoch, and
  the measurements (phases 1 to 3) used a 2 px width: the plan's 10 px disk learns nothing. Ruling
  R2's shift check (section 4) found the student not shift-invariant, so the augmentation stays,
  and on uncentred particles the student seeds about one in a hundred (phase 3's Progress).
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

**The network** (`simple_hhcnn`, `src/main/learn/`; `simple_circ_cnn` until ruling R5 renamed it in phase 5):
HHCNN, the Harmonic Hermitian convolutional neural network: Hermitian because every layer commutes with conjugation on the half ring with the conjugate wrap, so the network is a model of the Fourier transform of a real image; harmonic because the readout's class templates live in the angular harmonic domain of the polar transform and the rotation search is the cross-correlation over the cyclic shifts through the angular Fourier transform; convolutional on the polar Fourier representation. An encapsulated class with private components, `new`/`kill`, forward, backward, one Adam step, read and
write.

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
  (section 3, item 2); the deeper blocks see structure, not the envelope. Since 2026-11-09 the
  kernels are symmetric along the angle (the two angle offsets of a kernel share one weight, their
  gradients summed so that Adam keeps the symmetry bit for bit), so a reflected input gives reflected
  feature maps; with the readout's reflected score (next item) the network is invariant to the
  mirror on the rotation grid, as the mirror-pair classes require. The head as first
  designed (phases 1 and 2) pooled rotation-invariant statistics of the last block (the angular
  mean over the full ring, which for a Hermitian signal is the mean of the real part over the half
  ring, the mean of the modulus, and the magnitudes of eight harmonics of each), then a dense layer
  to 256 with ReLU, the side scalar appended, and a dense layer to the classes with a softmax; it
  was replaced by the readout of the next item (ruling R3), because the pooling discarded the
  angular structure that tells views apart and made a view and its mirror one class. Rotation
  equivariance is exact on the calculator's rotation grid. The per-particle cost is counted on the
  retained layout before it is quoted: 0.5 to 0.9 ms for the features and the forward pass at
  1,000 classes (phase 3).
- **The reflected score** (2026-11-09). A class scores the readout map and its reflection, whose
  angular spectrum is the conjugate, each maximised over the rotation grid, and takes the larger; the
  gradient follows the winner, through the conjugate spectrum when the reflection won. A class is a
  mirror pair, so its template matches either member, and the search decides which.
- **The readout replaces the invariant head** (Hans, ruling R3, 2026-11-07, after phase 3's folded
  cones and modulus-only variant showed the head to be the defect: it discards the angular structure
  that tells views apart, makes a view and its mirror one class, and without its phase-sensitive
  terms learns nothing). The three blocks stay. After them a pointwise complex-linear map with real
  weights reduces the channels to `NREAD` (8) and an average over ring bands reduces the rings to
  `KREAD` (8) bands (fewer when the last block has fewer rings); the result is Hermitian on the half
  ring and its full-ring angular spectrum is read to order `MREAD` (16), the even orders from the
  real part and the odd from the imaginary (the Hermitian extension's). One template per direction
  class lives in that domain (`NREAD` by `KREAD` by `MREAD` complex parameters, the order-0 term real),
  initialised from geometry as each class's reprojection through the network at the start and then
  trained with the rest. The class score is the maximum over the ring's rotation grid of the real
  Hermitian inner product of the feature map with the template (the spectra's product summed over
  channels and bands, the inverse transform of the ring's length, the real part), times a learned
  scale; the softmax over the classes is the posterior and the gradient flows through the maximising
  rotation. A template and its mirror differ, so the readout separates a view from its mirror by
  itself; ruling R1's fold and pair search stay as the backstop. The side scalar is no longer read
  (the interface keeps it). Shifts are unchanged (ruling R2: the features carry the phases). The
  model file's version is 2 and the weights file's 2: a model of the invariant head is refused.
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
- **The model file.** A key=value text header (`model_version=2` since the readout) and the weights as a binary file
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
  is the paired difference over micrographs, resampled by micrograph. As built, recovery is read
  on the student's own output without running the search: the teacher's direction, or its mirror
  partner, lies inside a folded cone within the cap and the route's bound (ruling R1); the paired
  difference is resampled 1,000 times by micrograph with a Lehmer generator seeded from the run
  seed and the wider design is kept when the 5% quantile of its mean gain is above zero; the
  held-out micrographs are every tenth stack, or every tenth selected particle of a single stack.
- **The confidence** is the student's posterior mass inside the cone the neighbourhood search will
  cover. Its threshold is set on the held-out micrographs so that, above it, the full search's
  direction lies inside the cone in a declared fraction of the particles; below it the particle
  takes the dense search. The threshold and the fraction are recorded in the model. As built, the
  declared fraction is 95% and the threshold the confidence of the longest confidence-sorted
  prefix of the eligible held-out rows that keeps it (above every confidence when none does), with
  the bound it was calibrated against; a student that seeds no held-out row therefore seeds nobody.
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
- **The folded cone's route and its search** (Hans, ruling R1, 2026-11-07). The mode, the cone and
  the confidence are read on the folded sphere (section 4); the mode written as the orientation is
  the heavier of the pair. The route rule bounds the folded cone's radius against the dense search's
  cell, `prob_athres` or the grid's spacing, whichever is larger, recorded in the model as the bound
  the threshold was calibrated with (`cone_bound_deg`), beside the count's cap: a folded set wider
  than the bound goes to the dense search. The confidence threshold is calibrated with recovery on
  the folded sphere: the teacher's direction, or its mirror partner, inside a cone within the cap and
  the bound. The neighbourhood search (`refine=prob_neigh prob_neigh_mode=geom`), given a particle's
  cone, scores the references within the cone about the stored direction and within the cone about
  its antipode (`euler_mirror`), both at the same radius, and takes the better pose as it takes any;
  a project without the slot runs exactly as before.
- **One source per particle.** Every particle's pose comes from the search that ran for it, the
  neighbourhood one seeded by the student or the dense one; the route is recorded per particle.

## 9. Where it fits

- `src/main/learn/`, picked up by the recursive glob of `src/CMakeLists.txt` with no CMake change:
  `learn.inf` (one `!@descr:` line), `simple_polar_feats.f90`, `simple_hhcnn.f90`,
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
- **Consumer 2 against the measured cycle** (phase 3, the proposal for Hans). On trpv1 the student's
  cycle at the curve's sizes (the rows, the training of the designs by the rule, the prediction, the
  search of the test half, `cont` and the reconstruction) costs 80 s of training at 180 rows, 2.5
  minutes at 900, 6.8 at 4,500 and 16 at the whole half of 17,752, plus 17 s of prediction and the
  searches, against 173 to 205 s for one dense `refine3D` iteration on the same test half; and it
  seeds nobody there, as on every real set at 15 Å (0 to 4% seeded), so the iteration would gain
  nothing. On the centred medium assembly the student seeds 46% of the particles from the whole
  half (826 s of training for three designs at 4,500 rows, 10 s of prediction, the cone search 35 s
  for the seeded rows against the dense search's 127 s for the rest), which would make a sampled
  iteration end with every particle aligned at the cost of a third of a dense search per seeded
  particle, but the whole-half training's early stopping and the real sets' broad teacher stand in
  the way. The proposal: consumer 2 is not warranted by the measured cycle as the programs stand;
  it becomes a candidate only where the particles are centred by an earlier stage, the teacher is
  sharp (a map the search already resolves) and the training is held to the cheapest design at a
  few thousand rows, which costs a few minutes against an iteration's three; the decision is the
  maintainer's. The sampling and fractional-update policy is not changed by this run.
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

Under ruling R3 (2026-11-07) the network's suite tests the readout: the gradient check through it,
the rotation test (a rotated input gives the same class scores and moves the maximising rotation by
the step), the mirror test (a class's tensor and its mirror image score differently against its
chiral template), the thread identity and the toy problem at its floor; the identity of the existing
routes is unaffected and not rerun. The order: build and test; the medium assembly at 15 Å, blank
and at the truth shifts, with the size-curve points at 5,000 rows and the whole half, the held-out
mode agreement and seed recovery reported against the invariant head's (5 to 7% and 0%); the real
sets of phase 3 run only if the new readout seeds at least a tenth of the held-out rows with the
teacher's direction inside the cone for 95% of them on the medium assembly, otherwise phase 3 ends
after the medium assembly with the result recorded. The prediction's cost per particle is reported
beside the coarse fill's own (2.2 ms per particle at the student's grid and band).

**The student-only read** (Hans, ruling R4, 2026-11-08). The reads above judge the student's direction
against one dense search at 15 Å from nothing, which on the real sets lands within 6° of the reference
alignment for only 24 to 69% of the particles, so disagreement with it does not say the student is wrong;
the student's orientations, drawn from the rows of the whole training set, may be better than that
alignment, and the maps judge. For every test half with a trained student, blank and at the shifts:
every selected row takes the student's folded mode as its direction regardless of the threshold, the
cap and the bound; the in-plane rotation and the shift come from the neighbourhood search restricted to
the cone about the mode and about its mirror at the larger of the folded 99% set's radius and the
bound, with no dense fallback; the test half is reconstructed from those poses with no local refinement
and again after a single `refine=cont` iteration; the reads are both maps' FSC against the training
half's map and, where it exists, the truth map (the value at 12 Å, the 0.143 and 0.5 shells) beside the
reference's reconstruction from the dense search's poses, the direction agreement of the modes, of the
cone search's poses and of the refined poses with the reference alignment (the truth on the truth sets)
for all rows beside the dense search's, and the cost per particle of the prediction and the cone search
beside the dense search's. On the real sets the agreement with the reference alignment is a diagnostic,
not the judge; the maps are. The read decides nothing by itself; the maintainer reads the maps. Its
table is in Progress (`scratch/keep/studentonly.sh`).

Under ruling R1 (2026-11-07) the reads of step 4 are taken with the folded cones, beside the
sphere-wide ones phase 2 recorded, on the medium assembly first (blank and at the truth shifts, 15
and 10 Å) and then on the real sets; the agreement of a direction with the reference is read on the
folded sphere where the student's class is meant (the cone), and on the full sphere for the final
pose the search returns. Ruling R2's shift check (section 4) is read on every prediction; if the
student is shift-invariant within noise, the blank and the truth-shift variants are one and the shift
augmentation is dropped from the training.

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
6. **Nightly:** a `pose_student` route on a simulated single-state set of the high-level tests: the
   split-half protocol in miniature, within the existing entry and area suites, with the actual
   runtime recorded as `doc/policies/test_environment_policy.md` requires. Built (phase 4) on the
   blob phantom of the multi-state fixture, not on 1jxy: beta-galactosidase is a D2 tetramer and at
   `c1` every posterior spreads over its four equivalent directions, so no cone is narrow and the
   student seeds nobody; the blob is asymmetric.

## 11. Order of work

- Phase 0, the feature module and the rows: `simple_polar_feats` with its tester, the row getter,
  and the teacher's rows on the Ribosembly medium assembly at each band, with their resolvability
  and their error against the truth.
- Phase 1, `simple_circ_cnn` (now `simple_hhcnn`, ruling R5) with its tester in the fast gate (tier 1), then the 16-channel network
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

1. The folder and the names: `src/main/learn/`, `simple_polar_feats`, `simple_circ_cnn` (renamed `simple_hhcnn` by ruling R5, 2026-11-08: the network is HHCNN, the Harmonic Hermitian convolutional neural network),
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
16. (2026-11-09) The mirror correspondence belongs to the search space, by the Euler angles, not to
    the network: the student classifies the grid's mirror pairs and is invariant to the mirror, the
    search resolves the pair (section 4). A tie of a pair's two templates in the polar Fourier domain
    was tried first and is not used (Progress, "Mirror pairs").

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

## 16. The run on the Dell

The development runs as the autonomous run `pose_student` on the Dell, with the engine of
`~/agent_runs/ms_regression/engine` and the spec in `~/tasks/pose_student_spec/` (run.conf,
prompt.md, phases.md with the launch procedure, rulings.md), the run root
`/data2/agent_runs/pose_student` behind the link `~/agent_runs/pose_student`. It is launched by
hand on the Dell after Hans has committed the synced tree and this note; the multi-state posterior
work runs on dell2 meanwhile and the two share no file. Sections 1 to 15 are the design; this
section is the run: the machine and the data, the phases with their exits, the rules, the file
table the fence enforces, and the Progress record the phases fill.

### 16.1 The machine and the data

**The Dell** (Oracle Linux 8.10, user elmlundho): Core Ultra 9 285, 24 cores (8 performance cores,
CPUs 0-7, and 16 efficiency cores; no hyperthreading), 62 GB; `/data2` 1.2 TB free, `/home` 16 GB
free, so nothing of the run goes under `/home`. GCC 15.2 from gcc-toolset-15, CMake 3.26,
`python3.12` for CMake (`python3` is 3.6, no numpy); `OPENBLAS_CORETYPE=Haswell` and
`TMPDIR=$HOME/tmp` are exported. The maintainer uses the machine during the run: timed
measurements record the load and are repeated under load; his checkout `~/src/SIMPLE` and its
`build/` are read-only for the run.

**Builds:** `./compile_clean.sh` (Release, tests on, the fast gate, install into `build/`) in about
4 minutes; increments by `make -C build -j20`, `scripts/run_fast_gate.sh build` (15 suites, 1-2
minutes), `make -C build install`. The source glob of `src/main` is recursive without
`CONFIGURE_DEPENDS`: reconfigure once after adding the `src/main/learn/` sources. The base HEAD is
exported with `git archive` into `scratch/keep/src_base` and built in `scratch/keep/build_base` for
the identity checks of phase 2; HEAD is never moved.

**Timings for planning** (the Dell, 2026-11-03): `refine3D` of trpv1, three iterations at 4 parts
by 6 threads, 17 minutes; `simulated_workflow_1jxy` 5 minutes at one thread; `test=nightly` about
2 hours; a table fill at 1,000 directions in the 15 Å band is unmeasured and is phase 0's first
measurement, on 2,000 particles before the whole half.

**Data** (read-only; projects copied into the run directory, `mkdir=yes`):

- The Ribosembly benchmark, `~/multistate_bench` (5.2 GB): `particles.mrcs` (80,000 particles, box
  128 at 3 Å per pixel), `bench.simple` with the true poses and shifts and the CTFs,
  `ground_truth_labels.txt` (one label per particle in stack order: 1 the small assembly, 24,000
  particles; 2 the medium, 48,000; 3 the large, 8,000), `reconstructed_volumes/` with the truth maps
  per label. The medium assembly is the run's first truth set, the small one the second. Its README
  gives the SNR (0.01), the shifts (uniform to plus or minus 10 pixels) and the CTF draw.
- `/data/simulated_vcrnapol`: two states with truth at SNR 0.1; `generate.sh` says how its
  `stack_0.1.mrcs`, `oris.txt` and `labels.txt` were made; state 1 (17,000 particles) is tier 3.
- trpv1: `/data2/trpv1_c2/5_solve3D/trpv1_c2.simple` (35,645 particles, c2, with its poses and
  `recvol_state01.mrc`), the stack `/data/trpv1/flipnew.mrc` (9.3 GB, box 256 at 1.2156 Å),
  `mskdiam=175`.
- bgal: `/data2/bgal_addon/15_abinitio3D/` (the solve3D of the selection; 5,513 particles, 2,756
  active, `mskdiam=180`), small enough only for the low end of the training-set size curve.
- streptavidin, `/data/10335-streptav/run_data.star` (RELION poses, `rlnAngleRot` present), and the
  small membrane protein FlhB, `/data/exp_gate_FlhB/BestFlhBparticles.star` with
  `BestFlhBparticles.mrcs`: imported with `import_starproject`; whether the FlhB file carries
  angles is checked in phase 3, and without them the set is aligned once with the standard machinery
  before the protocol or recorded as skipped.
- CESPED from Zenodo into `/data2/cesped/<entry>/`: the TEST entry (records 10998049 and 10998591,
  1 GB per half) for the smoke test of the import, then 10648 (records 8387534 and 8390711, 23 GB per
  half; the record's title reads 10647) and 11120 (records 8352296 and 8352300, 42 GB per half). Each
  half is a `particles_<h>.star`, an `info_<h>.json` and the stack in `particles_<h>.mrcs.chunk_N`
  pieces, concatenated in chunk order into one `.mrcs`, as `cesped/zenodo/downloadFromZenodo.py` in
  `~/src/cesped` does; the halves are the benchmark's own and are kept as such. Downloads start in
  phase 0 in the background and are used in phase 3 as far as they have arrived.

**The split-half protocol as the run instantiates it.** The view-balanced halves come from
`select_particles` (the program of the module `simple_ptcl_subset`; phase 0 found the program
name) with `eo_rebalance=yes`, which reassigns the selected particles' halves by direction
and score from the poses the project holds: the truth poses on the truth sets, the existing
alignment on the real sets (trpv1's solve3D, bgal's, the RELION poses of the imports). It needs
a score and an update count on every particle and the stacks' half-set markers, which on a truth
set one `refine3D refine=greedy_inpl maxits=1 trs=0 center=no` pass against the truth map
supplies (phase 0; `refine=eval` writes no score). The training
half's map is `reconstruct3D` from the training half's rows alone; the test half's reference poses
are the same alignment's, and on the truth sets the truth. A de novo solve3D of the training half
is the stricter variant and is run on the medium assembly in phase 3 if time allows.

### 16.2 The phases of the run

**Phase 0: foundations and the teacher's rows.**

1. Build the tree and the base export.
2. `src/main/learn/`: `learn.inf`; `simple_polar_feats` (the tensor of section 3 from the
   calculator: the band, the CTF-multiplied half ring packed by `sqrt(w_k)`, the ring-averaged CTF
   squared beside it, the conjugate-wrap padding, the shift augmentation through the calculator's
   shift matrices and the rotation by cyclic shift, the disk cache with the plane cache's contract)
   and its tester (a shifted image against the calculator's phase ramp; the wrap; the packed product
   against the objective's Hermitian product; the cache refusing a changed project, selection or
   geometry), registered as a suite of `unit_numerics` per the test policy.
3. The rows: a getter of a particle's posterior row from the dense table after the fill (the
   entries of largest mass with the tail's mass preserved) beside `posterior_row_stats`; the
   matcher's reference and calculator preparation exposed as a routine where the trainer needs it;
   `learn_pose_student` in its first form, a strategy routine `run_learn_pose_student(cline)` on
   `(params, build)` with a thin commander, registered as a private program: it prepares the
   references and the calculator at the student grid (`nspace`) and band (`lp`), fills the rows of
   the project's selected particles in batches with the production fill, writes them to the model
   directory with their resolvability (`pmax`, `neff`, `k99` per row and their distributions) and
   the cost per particle, and stops there until phase 1 adds the training.
4. The medium assembly's split-half data: its project from `bench.simple` and the labels, the
   view-balanced halves by `ptcl_subset eo_rebalance=yes`, the training half's map from its truth
   poses; the rows of the training half at 15, 10 and 8 Å; the teacher's resolvability per band,
   the row mode's error against the truth, and the cost per particle per band.
5. The CESPED TEST entry downloaded and imported; the downloads of 10648 and 11120 started in the
   background under `/data2/cesped`.

Exit: both builds pass the fast gate, the new suite included; the rows program runs on the medium
assembly's training half; the resolvability, the error against the truth and the costs are in
Progress.

**Phase 1: the network.**

1. `simple_hhcnn` (`simple_circ_cnn` until phase 5): the Hermitian layers of section 6 (real kernels applied alike to the real and
   imaginary channels, the conjugate wrap, the modulus gate with the first block's CTF
   normalisation, pooling along the rings, the angular mean of the real part and of the modulus,
   the dense layers and the softmax), forward and backward as im2col plus `sgemm`, Adam, the
   chunk-ordered reduction that is bit-identical on any thread count, read and write; its tester
   (the gradient check of every layer, the toy problem, the Hermitian contract of section 10's
   tier 1, the same model bit for bit at one thread and at a team of 4, the round trip), registered
   as a suite of `unit_numerics`.
2. `learn_pose_student` trains: the soft-target cross-entropy on the rows, early stopping on
   held-out micrographs, the model file with its record; `predict_poses` in its first form writes
   per particle the mode direction, the cone holding 99% of the posterior and the confidence to a
   text table in the model directory, no project write yet.
3. On the medium assembly: the student trained on the training half's rows at 15 Å with 16
   channels, predicted on the test half; the mode within 2 grid steps of the truth, the cone's
   coverage of the truth, the cone sizes; the first points of the training-set size curve at 200,
   1,000, 5,000 and all particles, with the training and prediction times.

Exit: the fast gate passes with the new suite; the exploratory numbers and the first points of the
curve are in Progress.

**Phase 2: the programs and the hand-off.**

1. `simple_pose_student` (the model: the grid, the threshold, the file I/O); `learn_pose_student`
   complete: held-out micrographs for early stopping, for the choice of width and band (the paired
   difference over micrographs, resampled by micrograph, the cheaper design kept unless the other
   recovers clearly more seeds) and for the confidence threshold; the selection table and the
   record of section 6.
2. `predict_poses` complete: the mode written as the particle's orientation, the cone's width, the
   confidence `pconf` and the route into a part file the caller merges, through new numbered ori
   slots in `simple_defs_ori`; a particle whose cone exceeds `KMAX_STUDENT` (64) directions or whose
   confidence is below the threshold marked for the dense search.
3. The neighbourhood search's cone per particle: `refine=prob_neigh prob_neigh_mode=geom` takes a
   particle's cone from its record where present and `prob_athres` otherwise, so that a project
   without student records runs exactly as before. Identity against the base build at a fixed seed
   and one thread on `simulated_workflow_1jxy` and on three `refine3D prob_neigh` iterations of
   trpv1, byte for byte on the declared outputs.
4. The split-half protocol end to end on the medium assembly: prediction, one iteration of
   `refine3D refine=prob_neigh` with the students' cones (the dense search for the marked rows),
   one iteration of `refine=cont`, the test half reconstructed; the reads of section 10 against the
   truth, at 15, 10 and 8 Å and the widths the choice rule picks; the complete training-set size
   curve at 15 Å with the cycle times (rows, training, prediction, neighbourhood search,
   refinement); the same on the small assembly and on vcrnapol's state 1.

Exit: the suites and the fast gate pass; the identity holds; the protocol's reads on the three
truth sets and the full curve are in Progress.

**Phase 3: real data and cost.**

1. The protocol on trpv1, bgal, streptavidin, FlhB where it has poses, the CESPED TEST entry, and
   10648 and 11120 as far as downloaded; the training-set size curves on trpv1 and 10648; the
   stricter variant with a de novo solve3D of the medium assembly's training half if time allows.
2. The cost table: every step of the student route on the test half (rows on the training half,
   training, prediction, neighbourhood search, refinement) against the dense search of the test
   half, as wall time, core time and peak memory on the quiet machine, three runs where wall time
   is reported; and the in-loop cycle time per subset size against one `refine3D` iteration on the
   same particles.

Exit: the reads and the costs are in Progress; no suite regressed.

**Phase 4: the nightly route and the documents.**

1. The `pose_student` nightly route: the split-half protocol in miniature on the simulated
   single-state set of the high-level tests, asserting cone coverage, seed recovery and fallback,
   registered in the registry, the test UI's help, the router and the CTest entry so that
   `scripts/check_test_registry.py` passes; its runtime recorded.
2. `doc/algorithms/pose_student.md`; this note's sections brought to what was built, the
   measurements in Progress, and the proposal on consumer 2 from the measured cycle, taken by Hans,
   not the run; skill text proposed in the report, not edited.
3. `test=nightly` on the final tree.

Exit: the nightly green, base failures recorded as such; the documents written.

**Phase 5: the report.** The review report beside the completed note
(`doc/implementation_notes/completed/learned_pose_state_assignment_report.md`), leading with what
was learned in plain words, then the evidence, the measurements, the skill text and what is left to
the maintainer, ending with a commit message; the note moved to `completed/` with `mv` and its row
of `doc/implementation_notes/README.md` updated.

### 16.3 The rules of the run

- Nothing is committed; git is read-only; the base is a `git archive` export.
- No new command-line key, public or internal: the band is `lp`, the grid `nspace`, the model
  `infile` and `fname`, the search's cone a per-particle record; if the work seems to need a key,
  stop and write it up.
- The existing routes are untouched in behaviour: a project without student records runs byte for
  byte as the base, checked in phase 2 and kept through phases 3 and 4. No edits to the stream, the
  flex folder, the quality model or the sampling policy.
- Single-state only; the multi-state extension and consumer 2's integration are out of scope, and
  the run measures what they would need.
- No floor is set from the first run of the code that is judged; the run measures, Progress
  reports, and Hans rules on the numbers (decision 16 of section 12 is his). The thresholds the
  programs set on held-out micrographs are recorded with their evidence.
- New or changed OpenMP code gets a tester at 16 threads, repeated, compared bit for bit with one
  thread; no string temporaries inside parallel regions; work arrays sized by the team in use; the
  BLAS on the calling thread inside a region.
- No `goto`; THROW macros never continued with `//&`; no `.and.` short circuit on allocatable
  arrays; `use ..., only:`; no large type as a plain component; a tester that needs private
  entities is a submodule of the module it tests.
- The datasets are read-only; copies live under `scratch`; bulky outputs are deleted once their
  numbers are in Progress; `scratch/keep/` holds the base export and build, the models and tables
  of the measurements and their logs; the CESPED downloads stay under `/data2/cesped`.
- Timed work runs with the load recorded; a measurement under load is repeated.
- Write prose for people: Progress entries and the report read without other documents, with
  SIMPLE terms explained and acronyms spelled out on first use; no item IDs in the report.

### File table

| file | why | phase |
| --- | --- | --- |
| `doc/implementation_notes/planned/learned_pose_state_assignment.md` | this plan and its progress (moved to `completed/` at the end of phase 5) | each |
| `doc/implementation_notes/completed/learned_pose_state_assignment.md` | this plan, completed | 5 |
| `doc/implementation_notes/README.md` | the note's row when it moves to completed | 5 |
| `doc/implementation_notes/completed/learned_pose_state_assignment_report.md` | the review report | 5 |
| `doc/algorithms/pose_student.md` | the algorithm page (phase 5: the HHCNN name, ruling R5) | 4 5 |
| `doc/algorithms/README.md` | the chapter list's entry for the pose student | 4 |
| `doc/policies/test_environment_policy.md` | the nightly entry's table row, its card and the entry count | 4 |
| `src/main/learn/learn.inf` | the folder's description line | 0 |
| `src/main/learn/simple_polar_feats.f90` | the input tensor, the augmentation, the cache (phase 5: the plan's path in a comment) | 0 1 2 3 5 |
| `src/main/learn/simple_polar_feats_tester.f90` | its tests | 0 1 2 3 |
| `src/main/learn/simple_circ_cnn.f90` | the Hermitian network, its trainer and I/O (moved to `simple_hhcnn.f90` in phase 5, ruling R5) | 1 2 3 |
| `src/main/learn/simple_circ_cnn_tester.f90` | its tests (moved to `simple_hhcnn_tester.f90` in phase 5) | 1 2 3 |
| `src/main/learn/simple_hhcnn.f90` | HHCNN, the network, its trainer and I/O (ruling R5: the module, the type `hhcnn`) | 5 |
| `src/main/learn/simple_hhcnn_tester.f90` | its tests (the suite `HHCNN`) | 5 |
| `src/main/learn/simple_pose_student.f90` | the model: the grid, the threshold, the record (phase 5: the renamed network) | 1 2 3 4 5 |
| `src/main/learn/simple_pose_student_tester.f90` | its tests | 1 2 3 4 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90` | the posterior row getter | 0 2 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_utils.f90` | the row record | 0 2 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_neigh.f90` | the cone per particle; phase 3: the cone about the antipode too (ruling R1) | 2 3 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_tester.f90` | the getter's and the cone's tests | 0 2 |
| `src/main/strategies/search/simple_strategy3D_matcher.f90` | the reference and calculator preparation exposed to the trainer | 0 2 |
| `src/main/pftc/simple_polarft_calc.f90` | the getters of a particle's CTF matrix and of the shell weights, beside `get_ptcl_pft` (phase 0: the feature module reads the calculator through them) | 0 |
| `src/main/pftc/simple_polarft_access.f90` | the getters' bodies | 0 |
| `src/main/strategies/parallelization/simple_strategy_pose_student.f90` | the strategy routines of the two programs (phase 5: the plan's path in a comment) | 0 1 2 3 4 5 |
| `src/main/commanders/simple/simple_commanders_pose_student.f90` | the thin commanders | 0 1 2 3 4 |
| `src/main/apis/simple_private_exec_api.f90` | the commanders' registration (phase 1: `predict_poses` in its first form) | 0 1 2 |
| `production/simple_private_exec_driver.f90` | the programs' cases (phase 1: `predict_poses`) | 0 1 2 |
| `src/main/ui/simple_private_prgs.f90` | the private programs' required keys (phase 1: `predict_poses`) | 0 1 2 |
| `src/defs/simple_defs_ori.f90` | the ori slots of the confidence, the route and the cone | 2 |
| `src/main/ori/simple_ori.f90` | the slots' names | 2 |
| `src/main/ori/simple_ori_tester.f90` | the slots' tests | 2 |
| `src/main/project/simple_binoris_tester.f90` | the record width's pins (the named slots, the padding) | 2 |
| `src/main/commanders/test/simple_commanders_test_class.f90` | the suite tables (phase 5: the suite `HHCNN`) | 0 1 2 5 |
| `src/main/ui/simple_test/simple_test_ui_class.f90` | the suites' help (phase 5: the suite key `hhcnn`) | 0 1 2 5 |
| `src/main/commanders/test/simple_commanders_test_highlevel_pose.f90` | the nightly route (phase 5: the plan's path in a comment) | 4 5 |
| `src/main/commanders/test/simple_commanders_test_highlevel.f90` | its interface | 4 |
| `src/main/commanders/test/simple_test_registry.f90` | its registry entry | 4 |
| `src/main/ui/simple_test/simple_test_ui_highlevel.f90` | its help | 4 |
| `src/main/exec/simple_test_exec_highlevel.f90` | its router case | 4 |
| `src/main/params/simple_parameters_tester.f90` | the test names the parameters tester lists | 4 |
| `production/CMakeLists.txt` | its CTest entry | 4 |

## Progress

(filled by the run, one entry per phase)

### Phase 0 (2026-11-06): foundations and the teacher's rows

**Files changed.** `src/main/learn/learn.inf`, `src/main/learn/simple_polar_feats.f90`,
`src/main/learn/simple_polar_feats_tester.f90` (new); `src/main/pftc/simple_polarft_calc.f90` and
`simple_polarft_access.f90` (the getters `get_ptcl_ctfmat` and `get_sigma2_weights`, added to the file
table for this phase); `src/main/strategies/search/probabilistic/simple_eul_prob_tab_utils.f90` (the
row record `prob_row_record` and the rows file `prob_row_file`), `simple_eul_prob_tab.f90` (the getter
`get_posterior_row`); `src/main/strategies/parallelization/simple_strategy_pose_student.f90` and
`src/main/commanders/simple/simple_commanders_pose_student.f90` (new: `run_learn_pose_student`,
`learn_pose_student_exec`, the thin commander); `src/main/apis/simple_private_exec_api.f90`,
`production/simple_private_exec_driver.f90`, `src/main/ui/simple_private_prgs.f90` (the private
program `learn_pose_student`; `pgrp` and `fname` added to the private programs' key dictionary, which
lacked them); `src/main/commanders/test/simple_commanders_test_class.f90` and
`src/main/ui/simple_test/simple_test_ui_class.f90` (the suite `polar features` of `unit_numerics`).
The matcher (`simple_strategy3D_matcher.f90`) was not edited: the reference and calculator
preparation it runs (`read_reprojection_model`, `prep_sigmas_objfun`, `alloc_ptcl_imgs`,
`build_batch_particles3D`, `materialize_reprojection_model_from_volumes`) is public already and the
strategy sequences it as `refine3D_exec` does.

**What was built.** `simple_polar_feats`: the student's input tensor of a particle from the polar
calculator through its public getters (the CTF-multiplied half ring over the band, real and
imaginary channels each scaled by the square root of the Euclidean shell weight k/sigma2(k) of the
particle's half and by the particle's weighted power, the ring mean of the CTF squared and its
azimuthal variation beside it, the log of the weighted power as the side scalar), the conjugate-wrap
padding and the rotation by cyclic shift, the shift augmentation through the calculator's phase
ramp, and the feature cache under the plane cache's contract (one direct-access record per
particle, the contract written last and the file renamed into place, refused on a changed project,
selection, geometry, spectrum or CTF treatment, at most a quarter of the free space). Its suite
`polar features` (38 checks, 0.01 s): a shifted ring against the calculator's own shift, the wrap
(padding, half turn, full turn, one grid step against the Friedel-extended ring), the packed
Euclidean product as the objective's weighted Hermitian product on the half ring and half the full
ring's, the weighted power as the calculator's self-term, the unit norm of the normalised tensor,
the ring mean of the CTF squared, the calculator path equal to the array path bit for bit, and the
cache's acceptance and nine refusals. No OpenMP code was added (the module is serial; the strategy
runs the existing fill), so no thread-identity check applies.

The getter `get_posterior_row` of the dense table returns, between a batch's fill and its
assignment, the row's entries of largest mass (at most a cap) with their masses and the tail's
mass, the row statistics (`pmax`, `neff`, `k99`, the valid count), the temperature (the tempering
factor times the particle's own, the factor the assignment draws at) and the mode's drawn in-plane
rotation with its refined shift; `prob_row_record` and `prob_row_file` in the table utilities hold
and stream these records with the whole row's losses.

`learn_pose_student` (a private program; its strategy routine on the caller's `(params, build)`)
prepares the project as refine3D's setup does (the halves' one maker, the sigma2 bootstrap, the run
seed), sets the band from `lp`, materialises the matcher's reprojection model at the student grid
(`nspace`, 1000 by default) from `vol1` with `center=no` so that frozen references move no
particle, and fills the rows of the selected particles batch by batch with the production fill
(`refine=prob`, `objfun=euclid`, the shift map's window path, `trs` 5 pixels by default), writing
beside `fname` the rows file (`_rows.bin`), the per-particle statistics with the error against the
particle's stored direction (`_rows.txt`), the feature cache (`_feats.bin`), the grid (`_grid.txt`)
and a header with the provenance (the project and its stamp, the selection's hash, the grid, the
band, the shell weights, the references' hash, the temperatures, the fill's budget) and the cost.

**The medium assembly's split-half data** (`scratch/keep/medium/`, `prep3.log`; the base build's
binaries, which this phase does not change). The label-2 particles (48,000) selected from the
benchmark's project (its stack path remapped to this machine), scored at their stored pose by
`refine3D refine=greedy_inpl maxits=1 trs=0` against the truth map (the balanced selection needs a
score and an update count on every particle; `refine=eval` writes no score: its strategy never calls
the assignment, a base behaviour recorded, not changed), then split by `select_particles
eo_rebalance=yes nspace=64`: T (eo 0) and S (eo 1), 24,000 each; over 144 direction bins of about
333 particles the mean difference between the halves is 11.3 and the largest 35, within 1-3% per
15-degree tilt band (`balance_M48.txt`). The exact truth poses were restored by `import_particles`
from the selected stack with an oritab merged from the two printouts (the truth project's poses, the
selection's halves; rows matched by their defocus), because refine3D's reference centring, on in
the scoring pass that ran with the default shift range, had mapped the truth map's centre-of-mass
shift onto every particle's stored shift (5-9 pixels per particle; a `center=no` lesson for every
later refine3D on the truth sets, which `learn_pose_student`'s commander applies). Projects:
`M48/M48.simple` (truth poses, eo = T/S) whose `reconstruct3D` gives the T map as its even map, the
S map as its odd one and their FSC (0.5 and 0.143 at 6.10 and 6.00 Å; the benchmark's README quotes
the truth-pose reconstruction's FSC at 12 Å as 0.985); `T_project/M48noeo.simple` (T alone, its
internal even/odd per index at first use) whose `reconstruct3D` from the truth poses
(`T_project/1_reconstruct3D/recvol_state01.mrc`) is the references of the rows; `S_project/`;
`T2000_project/` (every 12th row of T).

**The teacher's rows** (`learn_pose_student`, 24 threads on the otherwise idle machine, nspace
1000, `trs` 5, the shift map's window path, 55 references polished per particle, the shell weights
of the project's sigma2 bootstrap; `rows_T_lp<band>/`, `rows_T_<band>.log`):

| band | shells | fill per particle (wall, 24 threads) | program on 24,000 | peak RSS | pmax median (q10) | neff median (q90) | k99 median (q90) | mode-truth angle median, q90 | within 6° | truth grid point = mode | in the 64-entry cone | in the 99% set |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 15 Å | 2-26 | 1.88 ms | 56 s | 0.88 GB | 0.974 (0.627) | 1.14 (2.07) | 2 (3) | 2.68°, 4.31° | 99.25% | 77.6% | 100% | 96.3% |
| 10 Å | 2-38 | 2.62 ms | 76 s | 1.00 GB | 0.998 (0.747) | 1.01 (1.83) | 1 (2) | 2.61°, 3.98° | 99.88% | 83.4% | 100% | 95.8% |
| 8 Å | 2-48 | 3.46 ms | 97 s | 1.09 GB | 0.9996 (0.787) | 1.00 (1.72) | 1 (2) | 2.59°, 3.89° | 99.90% | 85.1% | 100% | 95.7% |

The first measurement, 2,000 particles at 15 Å, gave 1.72 ms per particle (fill 3.4 s, program 5.0 s,
peak RSS 0.71 GB) with the same distributions. The residual angle between the mode and the truth
(median 2.6°) is the grid's own spacing: about 6° between neighbours at 1,000 directions, so a random
direction lies about 3° from its nearest grid point; every mode is within 12° (two grid steps). No
row needs more than the 256 kept entries at any band; the two defocus bins (1.0-1.5 and 1.5-2.0 µm)
do not differ. The draws land on the row mode for 89%, 94% and 95% of the particles at 15, 10 and
8 Å, with an effective number of candidates of 1.37, 1.20 and 1.16 at the tempering factor 0.3.

**The deterministic alternative** (decision 4; `rows_T2000_lp15_diag/`, `prob_post_diag=yes`, which
needs `shift_sigma=0`, so these production rows take the seed path without the map and match the
map path's statistics to three digits): the best grid rotation per reference, no draw, puts the
truth's grid point at the mode for 88.6% of the 2,000 particles against 78.1% for the drawn rows;
the two modes agree for 89.0%, and where they differ the drawn mode is the neighbouring grid point
(the angle between the modes 0 at the median and 6.35° at the 90th percentile; the drawn mode's rank
in the deterministic row is 2 at the 90th percentile). The in-plane draw costs about ten points of
exact mode agreement, all neighbour swaps on the 6° grid; both teachers hold the truth inside the
cone for every particle.

**Costs for planning.** The dense fill at the student grid in the 15 Å band costs 1.9 ms per particle
wall on 24 threads (about 45 ms of thread time, of which the polish of 55 references is 80%): the
training half's rows at one band are a minute, and the rows of a 48,000-particle set at the three
bands under 10 minutes. The table fill the plan feared would be slow is not the cost of this design.

**CESPED.** The TEST entry (records 10998049 and 10998591) is downloaded, md5-verified per chunk
and concatenated in `/data2/cesped/TEST`; both halves imported with `import_starproject` as
`scratch/keep/cesped_test/test_h<half>/` (5,000 particles each, box 232 at 1.513 Å, poses, shifts and
CTFs present) from STAR copies whose image names carry absolute paths, because the import resolves
no bare stack name and prefixes the parent of the working directory to `import_dir` before
testing it against the working directory (`simple_commanders_starproject.f90` lines 47 and 69, the
second reading an unset string): a base defect recorded for the maintainer, not fixed. The downloads
of 10648 (records 8387534 and 8390711) and 11120 (records 8352296 and 8352300) run on under
`/data2/cesped` (`scratch/logs/cesped_chain.sh`), at about 1 GB per 8 minutes per entry.

**Exit.** Both builds pass the fast gate: the base export (`scratch/logs/build_base.log`, 15/15)
and the final tree with the new suite (`scratch/logs/build_final.log`, 15/15; the suite alone 38
checks, the area in reverse order 1,878 checks, all passed). The rows program ran on the training
half at the three bands; the resolvability, the error against the truth and the costs are above.
Base behaviours met and recorded, not changed: `refine3D refine=eval` writes no score; the
reference centring of `refine3D` (`center=yes` by default) moves a truth set's stored shifts;
`import_starproject` resolves no bare stack name and mishandles `import_dir`; `fname` must carry a
known extension. Journal: `state/journal_phase_0.md`; kept artefacts indexed in
`scratch/keep/summary.tsv`.

### Phase 1 (2026-11-06): the network, the training, the first size-curve points

**Files changed.** `src/main/learn/simple_circ_cnn.f90` and `simple_circ_cnn_tester.f90` (new: the
Hermitian network and its suite), `src/main/learn/simple_pose_student.f90` and
`simple_pose_student_tester.f90` (new: the model and its suite), `src/main/learn/simple_polar_feats.f90`
(the shift on a packed tensor, `shift_feat`), `src/main/strategies/parallelization/simple_strategy_pose_student.f90`
(the training after the fill; `run_predict_poses` and `predict_poses_exec`),
`src/main/commanders/simple/simple_commanders_pose_student.f90` (the `predict_poses` commander),
`src/main/apis/simple_private_exec_api.f90`, `production/simple_private_exec_driver.f90` and
`src/main/ui/simple_private_prgs.f90` (the private program `predict_poses`: `projfile`, `infile`, `mskdiam`
required; `nthr`, `fname` optional), `src/main/commanders/test/simple_commanders_test_class.f90` and
`src/main/ui/simple_test/simple_test_ui_class.f90` (the suites `Hermitian network` and `pose student` of
`unit_numerics`). No new command-line key.

**What was built.** `simple_circ_cnn`: three blocks of 3 by 3 real-kernel convolutions on the half ring
with the conjugate wrap along the angle and zero padding along the rings, the modulus gate
`z relu(|z| + b) / |z|`, average pooling along the rings, the invariant head, a dense layer of 256 with a
rectifier, the side scalar appended, the dense layer to the classes and the softmax; im2col plus `sgemm`
(the interface declared locally, `simple_linalg` keeps it private; this build's OpenMP OpenBLAS is
team-safe); the mini-batch cut into chunks of 8 computed in parallel into per-chunk gradient buffers and
reduced serially in chunk order, in double; Adam with decoupled weight decay; a Lehmer generator for the
seeded initialisation; the weights file (magic, version, shape, seed, arrays). Its suite (23 checks,
5 s at one thread): the Hermitian contract per layer (conjugation, the rotation by one grid step with the
wrap, the half ring against the explicit full ring bit for bit, the BLAS convolution against reference
loops, the pooled features real); the first block's CTF normalisation (the response to a filter's own
CTF-modulated template is one at four defocus curves with zeros inside the band; the raw response is the
packed product with the footprint, the objective's cross term by phase 0's contract); the gradient check
of every layer in single precision at the step 1e-2 (halved up to three times where a modulus or a
harmonic magnitude is near zero within the step) with the dense head in double at 1e-3; the same weights
and loss at one thread and at a team of 4, bit for bit; the round trip; the toy (projections of a six-blob
volume over 20 directions at box 32, rotated by a grid step, shifted within 2 px and noised at SNR 1 by
`add_gauran`'s rule on the particle path: real-space projections and noise, the image polarizer; width 8,
20 epochs at the rate 1e-2 in batches of 32 over 150 images per direction: held-out accuracy 0.988 on 400
images, 3 s). The thread identity was also run by hand at 24 threads, repeated.

`simple_pose_student`: the model (the geometry, the grid, the frozen shell weights, the network, the side
scalar's standardisation by the training rows, the provenance), the training (every tenth selected
particle held out, the keyed shift augmentation in a disk, early stopping on the held-out loss, the best
epoch's weights kept), the cone (the smallest set holding 99%, capped at 64, its mass and angular radius),
the header + weights file. Its suite (28 checks): the soft targets keep a row's tails, the cone is capped
at the dense search's marker, the file round trip is exact. `learn_pose_student` trains after the fill
and records the training in the model and the rows header; `predict_poses` (first form) writes per
particle the mode, `pmax`, `k99`, the cone's size, mass and radius, and, where the project holds a
direction, the error against it, to `<fname>_pred.txt`; no project write. The training budget is 2,500
optimiser steps in batches of a sixteenth of the training rows within 16 and 256 (30 epochs of the
training half at 256), the patience 250 steps, so that every point of the size curve gets the same number
of steps.

**Three defects found and fixed on the real rows, which the toy could not show.**

1. The first block's normaliser was fed the shell weights times the ring-averaged CTF squared, but the
   tensor already carries `sqrt(w_k)` ("a squared norm applies the shell weight once", section 3), so the
   weights entered twice: with `w_k` of 1e5 to 2e6 and unit-norm tensors the first block's outputs were
   about 1e-7 and the logits flat; the student learned the marginal over directions and nothing else, at
   every size. In the packed domain the data-domain filter is `u = w_tap / sqrt(w_k)` and flex's Gram
   diagonal `sum_k w_k ctf_k^2 u_k^2` is `sum_tap w_tap^2 ctf_k^2`: the ring-averaged CTF squared alone is
   the normaliser's input (`simple_pose_student` passes `feat%ctf2`). The toy, with unit weights, was
   blind to it; the contract check that catches it is now in the suite.
2. The plan's invariant pooling, the angular mean of the real part and of the modulus, is the zeroth
   angular harmonic of each feature map. Three blocks of 3 by 3 taps at angular stride one see 7 of the
   320 angular steps (8 degrees) before that mean, so the head is a bag of local-pattern statistics,
   where a direction's signature at 15 Å is the arrangement around the whole ring (the toy's ring of 16
   angles hid this: 7 steps cover 79 degrees). On the whole training half without augmentation the plan's
   head reached a held-out loss of 5.70 (uniform 6.91) and 3.0% exact-mode agreement; the head extended
   to the magnitudes of the angular harmonics 1 to 8 of each feature map reached 4.39 and 10.6%, and on
   the centred test half 30% of modes within two grid steps against 3% (below). The harmonic magnitudes
   are exactly rotation-invariant on the grid and real; the Hermitian extension's even full-ring orders are
   the real half-ring map's integer harmonics and the odd ones the imaginary map's half-integer harmonics,
   so the conjugation and rotation contracts hold unchanged. `CNN_NHARM = 8` in `simple_circ_cnn` (0 is
   the plan's head exactly); 18 features per channel and ring, 3,456 in all; 1.17 million weights against
   380 thousand. **A departure from section 6 for the maintainer's ruling.** The transforms are explicit
   sums, serial per sample (no FFT), so the thread identity holds.
3. The model header wrote the shell weights and the side scalar's standardisation with six significant
   digits; nine now (exact for single precision), so a read model builds the training's features exactly
   (the round-trip test found it).

Tried and dropped, with the evidence: the rate 1e-2 (no learning at all on the real rows; the toy liked
it); a data-dependent negative gate bias at initialisation (the gate at zero bias is the identity, so the
blocks start linear; the init changed nothing: 4.47 against 4.39, 10.5% against 10.6%). The width (32,
64) was not run: the time went to the shift question below.

**The medium assembly at 15 Å, 16 channels, the harmonic head** (`scratch/keep/medium/`; the training half
T of 24,000 with its truth map as the references, the test half S; every cone capped at 64 directions is
6.4% of the grid, the chance of the truth inside it; a mode within two grid steps (12°) by chance 1.3%).

*The centring.* The benchmark's particles carry truth shifts of 5.2 px on average, uniform to 10 px, in
both halves; `prepimg4align` shifts every particle by its stored shift before the polar transform, so
the training features are of centred particles, and S with its shifts zeroed ("the picking error stands")
is decentred by 5 px on average. A student trained without augmentation is therefore read twice: on S
centred (its truth shifts: the unseen-particle result at the shift the neighbourhood search would start
from if it were known) and on S zeroed, and on S decentred by exactly r pixels (the truth shift moved by
r along a seeded per-particle direction; `S_decentred/`), which is section 10's read of the accuracy
against the shift magnitude.

| training rows | shift disk | program (fill + training, 24 threads) | best epoch | held-out loss (uniform 6.91) | held-out exact mode | S centred: truth in cone / mode within 12° | S zeroed: in cone / within 12° |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 180 (+20 held out) | 0 | 18 s | 1 | 6.89 | 0 | 6.8% / 1.1% | 6.7% / 1.1% |
| 900 (+100) | 0 | 26 s | 1 | 6.88 | 0 | 7.0% / 1.2% | 6.8% / 1.1% |
| 4,500 (+500) | 0 | 1:49 | 9 of 23 | 6.74 | 1.8% | 22.3% / 5.6% | 8.9% / 2.0% |
| 21,600 (+2,400) | 0 | 9:48 (training 547 s) | 22 of 25 | 4.39 | 10.6% | 84.6% / 30.4% | 12.2% / 2.6% |
| 21,600 (+2,400) | 2 px | 12:12 (training 662 s) | 30 of 30 | 4.80 | 7.1% | 85.6% / 30.0% | 15.5% / 3.5% |
| 21,600 (+2,400) | 5 px | 12:15 (training 663 s) | 30 of 30 | 6.13 | 1.1% | 47.9% / 12.1% | 19.5% / 3.9% |
| 21,600 (+2,400) | 10 px | 8:12 (training 419 s) | 16 of 19 | 6.71 | 0.5% | 17.0% / 3.2% | 14.4% / 2.7% |

The plan's head (the angular mean alone) at disk 0: 21,600 rows: 5:59 (training 315 s), best epoch 15
of 18, held-out loss 5.75, held-out exact mode 2.6%, S centred 51.2% / 11.9% (the truth's rank median 61),
S zeroed 11.0% / 2.1%; its shift read 46.2% / 10.7% at 1 px, 33.7% / 7.3% at 2 px, 21.3% / 4.8% at 3 px,
9.4% / 2.1% at 5 px; 200 and 1,000 rows: chance (6.4-6.5% in the cone, 1.2-1.3% within two grid steps,
best epoch 1); 5,000 rows: best epoch 7 of 21, held-out loss 6.76, S centred 15.0% / 3.2%, S zeroed
7.6% / 1.1% (the harmonic head at 5,000: 22.3% / 5.6%). With the 10 px disk the plan's head learned nothing at 200 and
1,000 rows (held-out loss above uniform after the first epoch); the whole half was not run at that disk
with the plan's head (the first whole-half run at 10 px was under defect 1).

*The accuracy against the shift magnitude* (the whole half; the truth in the cone / the mode within two
grid steps on S decentred by r):

| disk | r = 0 | 1 px | 2 px | 3 px | 5 px | 8 px |
| --- | --- | --- | --- | --- | --- | --- |
| 0 | 84.6% / 30.4% | 76.5% / 25.6% | 51.8% / 14.3% | 26.8% / 6.3% | 10.6% / 2.1% | 7.8% / 1.4% |
| 2 px | 85.6% / 30.0% | 80.9% / 27.8% | 66.5% / 19.9% | 44.3% / 11.2% | 16.3% / 3.1% | 9.0% / 1.4% |
| 5 px | 47.9% / 12.1% | 47.0% / 11.6% | 44.6% / 10.6% | 40.2% / 9.4% | 28.3% / 5.9% | 14.9% / 2.8% |

A student trained on centred particles tolerates about one pixel at 15 Å (a shift of one pixel is a
phase ramp of 1.3 rad at the band's outer ring, k = 26 at box 128); the 2 px disk buys about a pixel at
no cost on centred particles; the 5 px disk gives a flat tolerance across its disk at about half the
narrow students' accuracy; the 10 px disk, the picking error's range here, defeats the fit. The 2 and
5 px students had not converged at the budget's 30 epochs (the held-out loss still falling), so a wider
disk also needs a longer schedule. The dense search's own shift range in the rows (`trs` 5) is not the
student's problem: the student is applied at whatever shift the particle has.

*The cones.* On centred unseen particles at 21,600 rows the posterior's 99% set has a median of 106
directions (71% of the particles capped at 64), the confidence (the mass inside the cone) a median of
0.97 and the truth's rank in the posterior a median of 11 (q90 92): the cone of 64 holds the truth for
85%, which the neighbourhood search would reach at a tenth of the dense search's references. At 4,500
rows the same student holds the truth for 22%; the curve rises steeply between 5,000 and 24,000 rows and
has not flattened, so decision 12's question (how few particles make a useful student) is open upward.

*Costs.* The rows at 24,000: 51 s (2.1 ms per particle with the shift map; 35 s on the seed path);
the training at 21,600 rows 16 s per epoch at the plan's head and 20-22 s with the harmonic head, 547 s
to the best epoch; prediction on 24,000 particles 10 s (0.29 ms per particle for the features and the
forward pass, the rest the reads); peak memory 1.1 GB. The student's forward pass per particle is a
small part of a dense fill's 1.9 ms, and the training of a whole half is under ten minutes.

**Exit.** The fast gate passes with the two new suites (`scratch/logs/build_gate_final.log`, 15/15; the
six static checks; `Hermitian network` 23 checks, `pose student` 28). The exploratory numbers and the
first points of the curve are above; the two design findings (the normaliser's input and the head) and
the shift question are for the maintainer: the harmonic head is in the tree as `CNN_NHARM = 8` with the
plan's head one constant away.

### Phase 2 (2026-11-07): the programs and the hand-off

**Files changed.** `src/defs/simple_defs_ori.f90` (the slots `pconf` 57, `cone` 58, `route` 59, absent at
0; 59 named slots within the record's 64 reals, so the project file's record width and every project of
the base read unchanged), `src/main/ori/simple_ori_tester.f90` and
`src/main/project/simple_binoris_tester.f90` (the slots and the padding pinned; the latter added to the
file table for this phase), `src/main/strategies/search/probabilistic/simple_eul_prob_tab_neigh.f90` (the
cone per particle), `src/main/learn/simple_pose_student.f90` (the threshold's record, the tensor-level
prediction, the grid's spacing), `src/main/learn/simple_pose_student_tester.f90`,
`src/main/strategies/parallelization/simple_strategy_pose_student.f90` (both programs completed),
`src/main/ui/simple_private_prgs.f90` (`part`, `nparts`, `outfile` optional keys of `predict_poses`; `part`
in the private key dictionary), `doc/implementation_notes/planned/learned_pose_state_assignment.md` (the
file table, this entry). No new command-line key.

**The hand-off.** `predict_poses` writes per selected particle its route: 1, the neighbourhood search from
the student's cone, when the cone holds 99% of the posterior within 64 directions and the confidence (the
posterior mass inside the cone, slot `pconf`) is at or above the model's threshold; 2, the dense search,
otherwise. A seeded particle's `e1, e2` become the mode's direction (its in-plane angle and shift stay the
search's), its slot `cone` the cone's angular radius about the mode, at least the student grid's spacing
(4.3° at 1,000 directions: direction 1 to its nearest neighbour; the grid's typical spacing is 6.5°). Run as a part (`part=`, `outfile=`) it writes its range's document as refine3D's
workers do, for a distributed caller's merge; run on its own it writes its project copy's `ptcl3D` segment,
never the caller's project. The table `<fname>_pred.txt` carries the route.

`refine3D refine=prob_neigh prob_neigh_mode=geom` takes a particle's cone from its record: every active
state's references whose direction lies within the cone of the particle's orientation under the point
group are scored and refined as the pooled neighbourhood's members are; a particle without the slot takes
its coarse cell exactly as before (the geometric mode never read `prob_athres`: its cell comes from
`nspace_sub`), so a project without student records runs the base's code path. An empty cone takes the
existing empty-cell fallback.

**The training complete.** `learn_pose_student` holds out every tenth stack of a project with several
(the micrographs) and every tenth selected particle of a single stack, trains the designs cheapest first
(16, then 32, then 64 channels) and keeps a wider one only when it recovers clearly more held-out seeds: a
seed is recovered when the teacher's direction (the dense search's mode for the same particle) lies in a
cone within the cap; the gain is the paired difference per micrograph, resampled 1,000 times by micrograph
with a Lehmer generator seeded from the run seed, and the wider design is kept when the 5% quantile of
the resampled mean gain is above zero. The confidence threshold is set on the held-out rows of the kept
design: the rows with a cone within the cap sorted by confidence, the longest prefix with at least 95% of
its teacher directions inside the cone; above every confidence when no prefix qualifies. The model records
the threshold, the declared fraction, the held-out share seeded and the recovered fraction among them;
the rows header records every candidate's recovery, held-out loss, best epoch, time and gain quantile, and
`<fname>_heldout.txt` the held-out table (particle, micrograph, confidence, 99% set, the teacher's and the
student's modes, recovered).

**Identity against the base export** (`scratch/identity/`). `simulated_workflow_1jxy` at one thread (the test fixes its own seed, 20260923) from the base build and from the tree, both passing their gates (docking correlations 0.87): the identity tool reports the 216 declared files the same, 0 differing or missing (`1jxy_identity_early.log`). Three `refine3D refine=prob_neigh prob_neigh_mode=geom` iterations of trpv1 (35,645 particles, c2, `box_crop=160`, `nparts=4 nthr=6`, `SIMPLE_SEED=20261107`) from the same project copy with the base build and with the tree (22-24 minutes each under the other runs): every volume of the three iterations, the four parts' half maps, the FSC binary, the reprojection models and the sigma2 state are identical byte for byte (44 files), and the particle fields (35,645 rows of `ptcl3D` and `ptcl2D`) and the stack segment are identical by their printouts; the 21 files that differ carry the run directory's path (the project's out, projinfo and jobproc segments) or the moment of the run (process ids, PDF plots, provenance sidecars, timing tables) (`trpv1_identity.txt`).

**The split-half protocol on the medium assembly at 15 Å** (`scratch/keep/medium/protocol/`). The training half T (24,000 particles with their truth poses, its map from them as the references) and the
test half S (24,000) of phase 0; the student at 15 Å trained on T with the shift augmentation disk of 2 px
(phase 1's shift read: the plan's 10 px disk learns nothing); the teacher's rows 56.7 s. The width rule:
16 channels recovered 3.9% of the held-out seeds (best epoch 28 of 30, held-out loss 4.86, 627 s), 32
channels 9.4% (30 of 30, 4.70, 1,529 s; the gain's 5% resampled quantile 4.6 points above zero, so 32 is
kept), 64 channels 8.6% (30 of 30, 4.69, 4,087 s; no gain over 32: the quantile -1.7 points): width 32.
The threshold 0.99 seeds 9.5% of the held-out rows, with the teacher's direction inside the cone for 99.6%
of them. The whole training 6,324 s with the three designs, under the identity runs' load.

The test half enters twice. *S at its truth shifts* (the particles centred, as an unseen particle would be
once a prepared shift is known): the student seeds 9.1% (2,185 particles) at the threshold, the truth
inside the cone for 98.9% of them (for 86% of all particles, the mode within two grid steps for 31%);
prediction 20 s; the neighbourhood search on the seeded rows 166 s (76 ms per particle: every one of
their cones has the radius 180°, see *What the cones are* below, so it scored all 20,000 directions of
`prob_neigh`'s grid for each of them), the dense search on the other 21,815 rows 210 s; after one
iteration of `refine=cont` the seeded rows end within 6° of the truth for 98.9% (12°: 99.7%) and the
dense rows for 98.8% (99.8%): the exhaustive search over the fine grid reaches the pose the dense search
gives. The test half's map against the training half's: FSC 0.949 at
12 Å; against the truth map 0.978. *S blank* (no poses, no shifts: the benchmark's picking error of 5 px on
average, to 10 px, stands, which is the protocol as section 10 states it): the student seeds 0.35% (85
particles) at the same threshold; those 85 end within 6° of the truth for 78.8% against 66.9% for the
23,915 dense rows (the reference, the dense search on every row: 66.9%, 70.8% within 12°); the map's FSC
0.912 at 12 Å against the training half's, 0.938 against the truth, the reference's values. The cycle on
the blank half: prediction 20 s, the neighbourhood search 14 s, the dense search 248 s, `cont` 9.8 s, the
reconstruction 5.6 s; the reference 242 + 8.5 + 5.4 s.

So on this benchmark the threshold, calibrated on centred held-out rows, declines nearly every decentred
particle, and the few it seeds, as the tenth it seeds of the centred particles, are aligned at least as
well as the dense search aligns the rest; but what the seeded rows get is not a local search, and the
recovery the threshold is calibrated on is granted by a two-mode set (*What the cones are*, below). The
decentring is this benchmark's own difficulty (the dense search itself falls from 98.8% to 66.9% within
6° when the shifts are unknown), and the student's shift tolerance (phase 1: about a pixel) is the
question the training's augmentation has to answer.

Two defects of the protocol's first attempt, found and fixed: the base build's `selection` wrote the
projects between the steps and its project writer, knowing 56 slots, dropped the student's three, so the
first neighbourhood jobs searched coarse cells instead of cones (the reads above are the rerun with the
tree's programs); and two 24-thread jobs side by side oversubscribe the machine's spinning OpenMP teams
(a dense search of 210 s took 1,444 s, a `cont` of 8.5 s 966 s), so every chain runs one job at a time.

**The training-set size curve at 15 Å with the cycle times.** The protocol end to end on the blank test half with the student trained on 200, 1,000, 5,000 rows
and the whole half (every design tried by the rule, one job at a time on the 24 cores; the whole half's
training ran under the identity runs' load):

| training rows | rows fill | the designs (16 / 32 / 64) | kept | held-out seeds recovered | threshold | seeded of S | dense search | neighbourhood search | `cont` + reconstruction | within 6° of the truth |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 180 (+20 held out) | 0.1 s | 17 / 48 s | 16 | 0 / 0 | none | 0 | 226 s on 24,000 | none | 9.3 + 5.5 s | 66.9% (the reference's) |
| 900 (+100) | 0.6 s | 25 / 65 s | 16 | 0 / 0 | none | 0 | 242 s | none | 9.2 + 5.3 s | 66.9% |
| 4,500 (+500) | 2.7 s | 133 / 290 s | 16 | 0 / 0 | none | 0 | 241 s | none | 9.1 + 5.3 s | 66.9% |
| 21,600 (+2,400) | 56.7 s | 627 / 1,529 / 4,087 s | 32 | 3.9% / 9.4% / 8.6% | 0.99 | 85 (0.35%) | 248 s on 23,915 | 14 s on 85 | 9.8 + 5.6 s | 66.9% (the 85 seeded: 78.8%) |

The same student on the centred test half (above) seeds 2,185 (9.1%) with the neighbourhood search at
166 s for them and the dense search at 210 s for the rest. The curve says what phase 1's said: at 15 Å on
this benchmark a student that seeds at all needs the whole training half; below it the threshold admits
no particle, and the protocol then costs the dense search plus a minute of training. The cycle times of a
seeded particle: prediction 0.8 ms, the neighbourhood search 76 ms, which at the sphere-wide cones the
student handed over is one in-plane scoring of each of `prob_neigh`'s 20,000 directions (the dense search
10 ms per particle on its 2,500-direction grid at this band); the cost of a cone that is a cone is not
measured on this set (vcrnapol below has them).

**The other bands and sets.** *The 10 Å band on the medium assembly* (the whole training half, the blank test half): the rows 71 s
(2.97 ms per particle at 37 rings); the designs 16 / 32 / 64 recovered 0.9% / 3.2% / 6.6% of the held-out
seeds (977 / 2,437 / 6,300 s; 32 over 16 at the quantile +1.7 points, 64 over 32 at +2.6), so the rule
keeps 64 channels here; the threshold 0.99 seeds 6.6% of the held-out rows with 99.4% of their teacher
directions in the cone. On the blank test half 61 particles (0.25%) are seeded and end within 6° of the
truth for 82.0% against 65.5% for the 23,939 dense rows (12°: 68.7%); the neighbourhood search 14 s, the
dense search 295 s, `cont` 11 s, the reconstruction 5 s; FSC at 12 Å 0.962 against the training half's map
and 0.976 against the truth (the 15 Å band's 0.912 and 0.938: the finer band's dense search gives the
better map). *The 8 Å band on the medium assembly*: the rows 95 s (3.97 ms per particle at 47 rings); the designs
recovered 0.5% / 1.9% / 4.3% of the held-out seeds (1,273 / 3,156 / 7,978 s; 32 over 16 at the quantile
+0.9 points, 64 over 32 at +1.8), so the rule keeps 64 channels here too, at 12,434 s for the three
designs; the threshold 0.99 seeds 4.4% of the held-out rows with 99.1% of their teacher directions in the
cone. On the blank test half 36 particles (0.15%) are seeded (prediction 98 s at this width and band) and
end within 6° of the truth for 80.6% (12°: 86.1%) against 68.4% (71.6%) for the 23,964 dense rows; the
neighbourhood search 13 s, the dense search 429 s, `cont` 13 s, the reconstruction 5 s; FSC at 12 Å 0.960
and 0.978. Across the bands the finer band costs more in rows (1.9, 3.0, 4.0 ms per particle) and in
training (the kept design's time 1,529, 6,300 and 7,978 s) and seeds fewer of the decentred particles
(0.35%, 0.25%, 0.15%), while the dense search's own result changes little (66.9%, 65.5%, 68.4% within
6°); the held-out recovery of the kept design falls from 9.4% at 15 Å to 6.6% and 4.3%. The band that
pays at 15 Å is the cheapest one, which the plan's order of trial expects.

*The small assembly at 15 Å* (24,000 particles of one class, halves of 12,000, its own truth map, the blank
test half): the dense search on every row ends within 6° of the truth for 47.9% (12°: 56.0%), FSC at 12 Å
0.562 against the training half's map and 0.612 against the truth, a harder set than the medium. The rows
35 s; the designs recovered 0.08% and 0% of the held-out seeds (657 and 1,047 s; the held-out mode 2.0%
and 1.1%), no gain, so 16 channels are kept, and the threshold 0.99 seeds 0.08% of the held-out rows and
nobody on the test half; the protocol then runs the dense search on all 12,000 rows (219 s), `cont` 5 s,
the reconstruction 3 s, and reads the reference's numbers to the digit (47.9%, 0.562 and 0.612), as a
student that seeds nobody must give. With 10,800 training rows and a weaker teacher this set sits below the
training-set size at which the student finds a confident mode, which the curve on the medium assembly
places between 5,000 and 24,000 rows.

*vcrnapol's state 1 at 15 Å* (`scratch/keep/vcrnapol/`: the generated set at SNR 0.1, 17,000 of its 30,000
particles, box 180 at 2.055 Å, halves of 8,500 by the small assembly's recipe, the map EMD-51275 the
set was simulated from as the truth, the blank test half): the dense search on every row ends within 6° of
the truth for 99.99% (68 s), FSC 0.997 at 12 Å against the training half's map and 0.998 against the
truth: an easy set with a near-perfect teacher, the opposite end from the small assembly. The rows 17 s
(2.0 ms per particle); the designs recovered 95.1%, 96.5% and 97.5% of the held-out seeds (823, 1,741 and
3,972 s; the held-out mode agreement 64.9%, 64.7% and 60.8%; 32 over 16 at the quantile +0.2 points,
kept; 64 over 32 at -0.1, not kept), 6,542 s for the three; the threshold 0.99 admits every held-out row
(96.5% recovered). On the test half 8,489 particles (99.9%) are seeded (prediction 10.7 s, 1.0 ms per
particle) and end within 6° of the truth for 99.69% (12°: 99.76%), the 11 dense rows for 100%, the
reference's every-row dense search for 99.99%; the neighbourhood search on the seeded rows 127 s against
the reference's 68 s on all rows; `cont` 7.8 s, the reconstruction 6.5 s; FSC 0.997 and 0.998, the
reference's. Where the student is confident of nearly everything, then, the hand-off costs more than the
dense search and loses 0.3 points, for the reason that follows.

**What the cones are** (read from the prediction and held-out tables, with the 1,000-direction spiral
regenerated by `make_oris even=yes`). On every seeded row of the medium assembly's runs (85, 2,185, 61
and 36 rows at 15 Å blank, 15 Å at the truth shifts, 10 Å and 8 Å) and of the small assembly's two rows
within the cap, the cone's radius is 180°: the 99% set holds the mode's neighbourhood and the antipode's
(the spiral's mirror partner), so the slot handed `prob_neigh` the whole sphere and the neighbourhood
search scored all 20,000 directions of its grid for every seeded particle. The seeded rows' better reads
above (78.8% against 66.9% within 6° at 15 Å, 82.0% and 80.6% at the finer bands) are the exhaustive fine
search's on confident rows, not the student's localisation, and the seeded rows' cycle times are that
search's. The held-out tables give the cause: among the medium assembly's confident held-out rows (a 99%
set within the cap) the student's mode lies within 12° of the teacher's for 48% and beyond 150° of it,
at the mirror view (the antipodal direction, whose projection is the reflected image), for 46% (at 8 Å
41% and 55%); "recovered" asks only that the teacher's direction lie in the 99% set, which the antipodal
second mode satisfies, so the width rule and the threshold were calibrated on a recovery the sphere-wide
set grants. On vcrnapol the mode is the teacher's for 98.5% of the held-out rows and the antipode's for
1.5%; 44% of the test half's cones are within 20° (a 99% set of one to three directions) and 54% span
the sphere (four or more directions: the fourth is the antipode); the narrow-cone rows end within 6° for
99.82%, the sphere-wide ones for 99.63%. Three things follow for the plan's section 10, for the
maintainer: the student as built does not separate a view from its mirror (whether the harmonic head's
magnitudes or the pooling over the ring discards the handedness is for the design to settle); a single
radius about the mode is the wrong summary of a two-mode set (two cones, or the set's directions, would
carry it); and the route rule by count (64 directions) lets a sphere-wide set through, where a bound on
the radius against the dense search's cell would not. No code changed for this in this phase: the
hand-off is as the plan specifies it, and the identity, the suites and the gate stand.

**Exit.** The suites and the fast gate pass on the final tree (`scratch/logs/fast_gate_phase2.log`, 15 of 15
with the six static checks; `pose student` 28 checks, `Hermitian network` 23, the slot pins of `unit_ori`
and `binoris`). The identity holds on both fixtures. The protocol's reads are above for the medium assembly
at 15 Å (both test-half variants, the reference), the curve's four points with their cycle times, the bands
10 and 8 Å, the small assembly and vcrnapol's state 1 (prepared from the ms_regression run's generated
stack with the small assembly's recipe). The finding on the cones is the phase's result for the
maintainer. The plan's section 6 text on the normaliser (phase 1) and the plan's `prob_athres` reading of
the geometric mode (it never read it) are recorded here, not rewritten in the sections.

### Phase 3 (2026-11-07): real data and cost

**Files changed.** `src/main/learn/simple_pose_student.f90` (ruling R1: the grid's mirror partners, the cone on the
folded sphere, the route bound in the model's record; ruling R2: the pooled features' getter),
`src/main/learn/simple_circ_cnn.f90` (the pooled feature vector of the last forward pass; ruling R3: the template readout in place of the invariant head) and `simple_circ_cnn_tester.f90` (the readout's rotation, mirror and gradient tests),
`src/main/learn/simple_pose_student_tester.f90` (the cone tests on the folded sphere),
`src/main/strategies/parallelization/simple_strategy_pose_student.f90` (the recovery and the threshold on the
folded sphere within the bound, the route rule, the folded diagnostics, the shift check),
`src/main/strategies/search/probabilistic/simple_eul_prob_tab_neigh.f90` (the cone about the antipode too; its
row of the file table extended to this phase), the plan's sections 4, 8 and 10 (the two rulings recorded) and this
entry; the scripts under `scratch/keep/` (the preparation of each set, `protocol_real.sh`, `reference_real.sh`,
`reads.py`, `subsets.sh`, `star31.py`, `maskcheck.py`, `shiftcheck.sh`) are the measurements' tooling. No new
command-line key.

**Ruling R1 applied** (Hans, 2026-11-07, on phase 2's "What the cones are"). The posterior is folded over the
spiral's mirror pairs (the grid carries its partner per direction) before the cone is read: the mode, the 99% set
and the radius are taken on the folded sphere, the written direction is the heavier partner, and the slot carries
the folded radius. The route admits a particle when its folded set is within the cap, its radius within the bound
(`prob_athres` or the grid's spacing, whichever is larger: 10° at the defaults, recorded in the model as
`cone_bound_deg`) and its confidence at or above the threshold, which is recalibrated with recovery on the folded
sphere (the teacher's direction or its mirror partner inside a cone within the cap and the bound). The neighbourhood
search scores the references within the cone about the stored direction and about its antipode (`euler_mirror`) at
the same radius and takes the better pose as it takes any; a project without the slot runs the base's path: the
1jxy identity rerun from the base export and from this tree at one thread, 216 declared files the same, none
differing (`scratch/identity/1jxy_identity_p3.log`). The fast gate passes with the folded cone tests (`scratch/logs/build_p3_fold2.log`,
15 of 15; pose student 34 checks).

**Ruling R2: is the student a network of spectral powers?** No. `predict_poses` now predicts its first 256 selected
particles again with their tensors shifted through the calculator's phase ramps by 2, 5 and 10 pixels in four
directions and reports, per magnitude, the posterior's largest total variation, how often the folded mode is kept
and the unshifted mode stays inside the shifted cone, the cones' Jaccard index, and the relative change of each
pooled feature group (`<fname>_shiftcheck.txt`; `scratch/keep/shiftcheck/`). On every trained model the posterior
moves with the shift: on the medium assembly's 15 Å model, read on the centred test half, a 2 px shift changes the
posterior by a total variation of 0.68 on average (0.80 at the 90th percentile), keeps the folded mode for 25% of the
shifted copies and the unshifted mode inside the shifted cone for 96% (the cones' Jaccard 0.47); at 5 px the total
variation is 0.95 and the mode is kept for 1%; at 10 px 0.97 and 0.3%. The pooled features say which carries the
shift: the angular mean of the real part changes by 65% of its size at 2 px (118% at 5 px, 132% at 10 px) and the
harmonic magnitudes of the real and the imaginary maps by 42% (57%, 58%), while the mean modulus changes by 3.4%
(7.3%, 10.1%). The 10 and 8 Å models, the small assembly's, vcrnapol's and trpv1's read the same way (the modulus
mean 1.7-3.6% at 2 px, the real-part mean 27-74%, the harmonics 22-44%; trpv1's model is the least sensitive, a
total variation of 0.28 at 2 px with the mode kept for 51%). So the student is not shift-invariant by construction:
it uses the phases, the shift augmentation stays, the blank and the truth-shift variants stay two reads, and phase
1's "about a pixel of tolerance" was this dependence measured coarsely. R2's second branch was tried as an experiment on a copy of the tree
(`scratch/variant_src`, built in `scratch/variant_build`; the two edits of `simple_circ_cnn.f90` zero the real-part
means and the harmonic magnitudes in the forward pass and their gradient in the backward one, so the dense head reads
the modulus means alone): that head is shift-invariant to the modulus feature's own drift (on the medium assembly at
15 Å a 2 px shift changes the posterior by a total variation of 0.036 and keeps the folded mode for 85% of the
copies; 0.09 and 71% at 5 px; 0.16 and 59% at 10 px) but learns almost nothing: its held-out loss stays at 6.77 for
16 and for 32 channels against 4.86 and 4.70 for the full head (the uniform posterior over 1,000 directions scores
6.91), the held-out mode agreement is 0.3% against 5.4%, no held-out seed is recovered, and in the tester the
Hermitian network's toy problem (held-out accuracy of at least 0.95 on 20 directions) fails. So the shift-dependent
features are the ones that carry the pose; a head without them is invariant at the loss of the task, and the
student keeps the phases. The variant is not in the tree.

**The protocol on the real sets, as instantiated.** The reference alignment of a set is the existing one (trpv1's and
bgal's solve3D, the RELION poses of the STAR imports, the CESPED entries' own); one `refine3D refine=greedy_inpl
maxits=1 trs=0 center=no` pass on a copy supplies the score and the half-set markers `select_particles
eo_rebalance=yes` needs, the view-balanced halves T and S are its, and the reference poses are restored from the
source project's printout (the CESPED entries keep the benchmark's own halves: half 0 trains, half 1 tests). The
training half's map is `reconstruct3D` from its rows. The test half enters blank (no poses, no shifts) and at the
reference shifts without angles; the reads are the final direction's agreement with the reference alignment under
the point group (within 6°, 12°, 20°), per route, the FSC of the test half's map against the training half's beside
the FSC of the two halves under the reference alignment, the seeded fraction and the cones' radii, and every step's
wall time, core time (user plus system) and peak memory from `/usr/bin/time -v`. The reference run is the dense
search (`refine=prob`, one iteration) on every row of the same test half, then the same `cont` iteration and
reconstruction. The RELION 3.0 STARs of streptavidin and FlhB (no optics block) were rewritten in the 3.1 layout for
`import_starproject`; FlhB's per-micrograph stack names were mapped onto its consolidated stack, whose order the
maintainer's defocus table confirms row for row.

**The medium assembly with the folded cones** (ruling R1, item 5; `scratch/keep/medium/protocol/med_lp*_fold_*`).
At 15 Å, trained on the whole training half with the 2 px augmentation as before: 16 channels recover 0.00% of the
held-out seeds under the folded rule (held-out loss 4.861, mode agreement 5.4%, 816 s) and 32 channels 0.04% (one
row of 2,400; loss 4.699, mode agreement 6.8%, 1,922 s), so 16 is kept and the threshold admits no row: the student
seeds nobody on the blank test half nor on the centred one, and the protocol's reads are the reference's to the
digit (blank: 66.9% within 6° of the truth, FSC 0.912 at 12 Å against the training half's map and 0.938 against the
truth; at the truth shifts 98.85%, 0.949 and 0.978). Beside the sphere-wide reads of phase 2 (85 seeded at 78.8%
within 6° on the blank half, 2,185 at 98.9% on the centred one) the folded reads say that those seeds were the
whole-sphere search's: no held-out particle of this set has a confident 99% set within 10° of its mode that holds the
teacher's direction. At 10 Å the same: 16 and 32 channels each recover one held-out row of 2,400 (0.04%; held-out
loss 5.246 and 4.939, mode agreement 4.8% and 5.0%; 1,002 and 2,418 s), 16 is kept with a threshold that admits
that one row, and on the test half nobody is seeded: the reads are the reference's (blank 65.6% within 6°, FSC 0.962
and 0.976; at the truth shifts 99.5%, 0.966 and 0.990). The cost of the folded protocol on this set is therefore the
dense search's plus the training (2,742 s at 15 Å and 3,426 s at 10 Å for two designs, 10 s of prediction). 

**Ruling R3: the template readout.** The invariant head is replaced (section 6): after the three Hermitian blocks
a pointwise complex-linear map with real weights to 8 channels, the mean over 8 ring bands (the rings the last
block has when fewer), the full-ring angular spectrum to order 16 (even orders from the real part, odd from the
imaginary), one template per direction class in that domain initialised from the class's reprojection through the
network at the start, the class score the maximum over the ring's rotation grid of the real Hermitian inner product
times a learned scale, the softmax, the gradient through the maximising rotation; the dense layers and the side
scalar's use are gone, the model and weights files are version 2. The network's suite (26 checks) adds the gradient
check through the readout (the pointwise map, the templates, the scale), the rotation test (the class scores
unchanged by a one-step rotation and every maximising rotation moved by one step; by a conjugation, by half a
turn), the mirror test (a class's tensor and its mirror image score differently against its chiral template), and
the toy problem at its floor; the thread identity and the round trip stand. The existing routes are untouched by
it (the identity is not rerun, as ruled).

*The medium assembly at 15 Å with the readout* (`med_lp15_r3_*`): 16 channels reach a held-out loss of 1.781 and a
mode agreement of 46.0% at the third epoch (the invariant head's 4.861 and 5.4%), recovering 44.8% of the held-out
seeds under the folded rule (0% before), in 163 s; 32 channels 1.723, 45.8% and 41.2% (319 s; no gain, 16 kept);
the threshold 0.99 seeds 44.8% of the held-out rows with every teacher direction inside the cone, so ruling R3's
gate for the real sets (a tenth, 95%) is passed. On the centred test half 11,055 of 24,000 particles (46.1%) are
seeded, their cones 6.5°, 7.8° and 9.5° wide at the 10th, 50th and 90th percentile (none beyond 20°), the folded
mode within two grid steps of the truth for 99.1% of all rows; the prediction costs 0.56 ms per particle for the
features and the forward pass against the coarse fill's 2.2 ms, the neighbourhood search 35 s for the seeded rows
(3.2 ms per particle) against the dense search's 127 s for the other 12,945 (9.8 ms per particle); after `cont`
the seeded rows end within 6° of the truth for 99.0% and the dense rows for 98.7%, the map FSC 0.949 and 0.978 as
the reference's. On the blank test half (uncentred) the student seeds 230 particles (1.0%; 94% of the rows beyond
the cap, the shift check's total variation 0.81 at 2 px) which end within 6° for 98.7% against 66.6% for the dense
rows, the map the reference's: the readout's student is a seeder of centred particles. The 5,000-row point (4,500 training rows, 500 held out)
keeps 64 channels by the rule (16: held-out loss 1.771, mode agreement 47.0%, 45.4% of the held-out seeds
recovered, 117 s; 32: 1.644, 52.0%, 53.0%, 199 s; 64: 1.540, 54.4%, 63.0%, 483 s), its threshold seeds 63.2% of
the held-out rows with 99.7% recovered, and on the blank test half it seeds 409 particles (1.7%) which end within
6° for 90.2% against 66.4% for the dense rows; that a tenth of the training rows recovers more held-out seeds
than the whole half says the whole half's training stops too early (its best epoch is the third of six at 85
steps per epoch, the budget's early stopping) and is for the maintainer.

**The real sets with the readout** (ruling R3 item 7 admitted them; the invariant head's runs of the same sets,
made earlier in the phase, seeded nobody and are in the journal). Each set's reference run is the dense search
(`refine=prob`, one iteration at 15 Å) on every row of its test half from no poses and no shifts, then `cont` and
the reconstruction, read against the set's reference alignment; the student is trained on the training half as on
the medium assembly and the protocol run on the blank test half and, where the halves carry the reference shifts,
on the test half at those shifts with the model reused.

*trpv1* (`scratch/keep/trpv1/`: 35,504 state-1 particles of the solve3D project, box 256 at 1.2156 Å, c2, mask
175 Å, the phase-flipped stack; halves of 17,752 by `eo_rebalance` after a `greedy_inpl` pass, the solve3D poses
restored; the two halves' maps under the reference alignment agree to FSC 0.97 at 12 Å). The reference on the
blank half puts 34.3% of the particles within 6° of the solve3D direction (42.6% within 12°, median 23°; FSC 0.5
at 8.2 Å and 0.143 at 6.6 Å against the training half's map), at the reference shifts 42.5% (median 10.3°; 8.0
and 6.5 Å): one iteration at 15 Å from nothing recovers a third of the poses the set's solve3D settled over many.
The student (the rows 39 s, 2.2 ms per particle): 16 channels reach a held-out loss of 4.247 and a mode agreement
of 11.6% (the invariant head's 4.826 and 4.5%), 32 channels 4.217 and 12.5%, each recovering one held-out seed of
1,775 under the folded rule; 16 is kept and seeds nobody on either half (16% of the blank half's rows have a 99%
set within the cap, 53% of the centred half's, but every cone is wider than the 10° bound; the folded mode lies
within two grid steps of the solve3D direction for 6% and 9%); prediction 0.69 ms per particle. The size curve at
200, 1,000 and 5,000 rows seeds nobody either (held-out mode agreement 9 to 11%). The protocol's reads are the
reference's; the costs are in the table.

*bgal* (`/data2/bgal_addon/15_abinitio3D/`, 4,859 state-1 particles, box 256 at 1.275 Å, mask 180 Å; halves of
2,430 and 2,429): the reference 24.4% within 6° of the solve3D direction on the blank half (median 77°; FSC 0.143
at 4.35 Å against the training half's map, the set being small and half aligned) and 35.0% at the shifts; the
student (16: held-out loss 3.117, mode agreement 23.5%; 32: 2.963, 24.7%; one held-out seed recovered of 243) seeds
nobody; the dense search 38 s.

*streptavidin* (`/data/10335-streptav/run_data.star`, 11,402 particles in 405 per-micrograph stacks, box 128 at
1.08 Å, mask 80 Å, the RELION poses; halves of 5,701): the reference 6.5% within 6° of the RELION direction on
the blank half (median 91°; FSC 0.143 at 9.9 Å) and 13.2% at the shifts (median 80°): one dense iteration at 15 Å
does not align a 50 Å particle from nothing; the student's held-out loss stays at 5.8 (mode agreement 3%), 7 of
5,701 test rows have a set within the cap, nobody is seeded.

*CESPED TEST* (5,000 particles per half of the 11120 entry's optics, box 232 at 1.513 Å, mask 180 Å; half 0
trains, half 1 tests; the two halves' maps under the RELION alignment agree to FSC 0.143 at 6.1 Å): the student (16: held-out loss 4.715, mode agreement 10.0%; 32:
4.663, 10.6%; no held-out seed recovered of 500) seeds nobody; the dense search from nothing 36.5% within 6° of
the RELION direction (median 11.6°; FSC 0.143 at 9.5 Å), at the shifts 42.3% (median 8.4°); 52 s on 5,000 rows.

*FlhB* (`/data/exp_gate_FlhB/`, 162,408 particles, box 288 at 0.822 Å, the RELION poses, mask 170 Å for the
alignment; halves of 81,204, the student trained on a class-balanced 20,000 of the training half): the reference
29.4% within 6° of the RELION direction (39.1% within 12°, median 31°; FSC 0.143 at 5.9 Å) in a dense search of
40 minutes (55,600 core-seconds, 29.7 GB); the student (the rows 1.55 ms per particle; 16: held-out loss 4.864,
mode agreement 11.6%, 0 recovered of 2,000; 32: 2 recovered) seeds nobody; prediction 76 s on the test half (0.94
ms per particle).

*CESPED 10648* (PKM2, D2; 117,478 particles per half, box 222 at 1.497 Å, mask 150 Å; half 0 trains, half 1
tests): the reference's dense search from nothing puts 68.9% of the test half within 6° of the RELION direction
(80.8% within 12°, median 3.2°; FSC 0.143 at 3.96 Å against the training half's map) in 16.7 minutes (22,900
core-seconds, 24.7 GB): a sharp teacher at 15 Å, the opposite of the other real sets. The size curve of the readout's
student on class-balanced subsets of the training half, each point the whole protocol on the blank test half:

| training rows | held-out loss (16 / 32) | mode agreement | held-out seeds recovered | kept | seeded of the test half | the seeded rows within 6° | the dense rows | neighbourhood search | dense search |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 180 (+20 held out) | 3.87 / 3.83 | 0% | 5% (1 of 20) | 16 | 3,572 (3.0%) | 66.4% | 68.8% | 78 s | 16.6 min |
| 900 (+100) | 3.75 / 3.91 | 11% | 6% | 16 | 5,204 (4.4%) | 68.3% | 68.7% | 91 s | 16.2 min |
| 4,500 (+500) | 4.07 / 4.03 | 8.6% | 3.6% | 16 | 4,156 (3.5%) | 68.6% | 68.7% | 82 s | 16.5 min |
| 18,000 (+2,000) | 4.48 / 4.48 | 7.9% | 0.95% | 16 | 2,361 (2.0%) | 77.3% | 68.7% | 67 s | 16.7 min |

The maps are the reference's at every point (FSC 0.143 at 3.9 to 4.0 Å). The student seeds a few percent of this
entry from 180 training rows on, those rows aligned as the dense search aligns the rest (and better at 18,000
rows: 77% against 69%), but the neighbourhood job's fixed cost (the 20,000-direction references) makes its
seconds per seeded particle higher than the dense search's; and the held-out loss worsens with the training set
(3.75 at 900 rows, 4.48 at 18,000): the early stopping takes the second to fifth epoch at 70 steps per epoch,
which, with the medium assembly's whole-half training stopping at the third of six, says the training budget's
early stopping is tuned for small sets and is for the maintainer.

*CESPED 11120* (the TSHR complex, C1; 194,019 and 194,020 particles per half, box 232 at 1.513 Å, mask 220 Å for
the alignment; half 0 trains, half 1 tests; the student on a class-balanced 20,000 of half 0): the reference's
dense search from nothing puts 45.8% of the test half within 6° of the RELION direction (61.9% within 12°, median
7.1°; FSC 0.143 at 4.39 Å against the training half's map) in 38.6 minutes (53,400 core-seconds); the student (the
rows 0.79 ms per particle; 16: held-out loss 3.559, mode agreement 17.7%, no held-out seed recovered of 2,000; 32:
3.480, 19.1%, one) seeds nobody and the protocol's reads are the reference's.

**What the real sets say.** On every real set the dense search's rows at 15 Å are a broad teacher (the reference
alignment is far from one iteration's own answer: 6 to 42% within 6°), and the readout's student learns that
breadth: its held-out mode agreement rises from 3-5% to 10-25% but no held-out particle gets a confident 99% set
within 10°, so the threshold admits nobody and the protocol costs the dense search plus the training. The medium
assembly's teacher, the same search against a map from truth poses, is sharp enough for the student to seed 46%
of the centred particles; a real set's teacher at 15 Å after one iteration is not. The student can only hand
over what its teacher resolves.

**Ruling R4: the student-only alignment, the maps judging** (`scratch/keep/studentonly.sh`; the runs `so_*` under
each set's `protocol/`). Every selected row of the test half takes the student's folded mode as its direction with
no threshold, cap or bound; the neighbourhood search in the cone about the mode and its mirror at the larger of the
99% set's radius and the bound (most cones are wide, so this search is near-exhaustive and costs more than the
dense search per particle) gives the in-plane rotation and the shift; the half is reconstructed from those poses
and again after one `cont` iteration. The table gives, per test half, the share of all rows whose direction lies
within 6° of the reference alignment (the truth on the truth sets) for the modes, the cone search's poses and
the refined poses beside the dense search's, and both maps' FSC against the training half's map (the value at 12
Å; the 0.5 and 0.143 shells) beside the reference's reconstruction from the dense search's poses; `cent.` is the
test half at the reference shifts (centred), `blank` without poses or shifts.

| set, variant | modes | cone search | after `cont` | dense | map, no refinement: at 12 Å; 0.5 / 0.143 | map after `cont` | reference map |
| --- | --- | --- | --- | --- | --- | --- | --- |
| medium, blank | 9.2% | 27.8% | 40.8% | 66.9% | 0.639; 11.3 / 7.1 Å | 0.819; 8.7 / 6.0 Å | 0.912; vs truth 0.938 |
| medium, cent. | 88.6% | 97.1% | 98.9% | 98.9% | 0.944; 6.1 / 6.0 Å (vs truth 0.970) | 0.949; 6.1 / 6.0 Å (vs truth 0.978) | 0.949; vs truth 0.978 |
| trpv1, blank | 21.3% | 25.0% | 30.5% | 34.3% | 0.741; 9.7 / 8.0 Å | 0.829; 8.4 / 7.2 Å | 8.2 / 6.6 Å |
| trpv1, cent. | 42.8% | 41.5% | 44.9% | 42.5% | 0.860; 8.2 / 6.0 Å | 0.884; 8.0 / 5.4 Å | 8.0 / 6.5 Å |
| bgal, blank | 21.0% | 13.8% | 17.9% | 24.4% | 0.646; 9.9 / 6.7 Å | 0.734; 9.1 / 5.1 Å | 8.0 / 4.35 Å |
| bgal, cent. | 56.4% | 34.7% | 36.2% | 35.0% | 0.824; 7.8 / 4.3 Å | 0.830; 7.6 / 4.3 Å | 7.6 / 4.3 Å |
| streptavidin, blank | 5.0% | 4.8% | 5.5% | 6.5% | 0.154; 13.8 / 11.5 Å | 0.206; 12.6 / 11.5 Å | 12.6 / 9.9 Å |
| streptavidin, cent. | 12.5% | 11.2% | 13.4% | 13.2% | 0.773; 8.1 / 4.3 Å | 0.820; 7.3 / 4.3 Å | 4.3 Å |
| CESPED TEST, blank | 22.4% | 31.5% | 35.4% | 36.5% | 0.336; 14.0 / 10.3 Å | 0.388; 14.0 / 9.5 Å | 9.5 Å |
| CESPED TEST, cent. | 28.9% | 45.4% | 48.2% | 42.3% | 0.414; 14.0 / 8.8 Å | 0.380; 14.0 / 9.2 Å | 10.3 Å |
| FlhB, blank | 26.0% | 26.6% | 33.6% | 29.4% | 0.842; 7.6 / 6.2 Å | 0.886; 7.4 / 5.9 Å | 5.9 Å |
| 10648, blank (900-row student) | 35.7% | 49.1% | 60.9% | 68.9% | 0.853; 4.9 / 4.3 Å | 0.906; 4.4 / 3.8 Å | 4.0 Å |
| 10648, blank (18,000-row student) | 38.8% | 48.7% | 60.3% | 68.9% | 0.855; 5.3 / 4.2 Å | 0.909; 4.4 / 3.9 Å | 4.0 Å |
| 11120, blank (18,000-row student) | 32.7% | 38.4% | 46.5% | 45.8% | 0.877; 7.6 / 4.6 Å | 0.908; 7.0 / 4.4 Å | 4.4 Å |
| small, blank | 7.7% | 28.9% | 38.8% | 47.9% | 0.362; 14.2 / 10.4 Å (vs truth 0.386) | 0.492; 14.2 / 9.4 Å (vs truth 0.530) | 0.562; vs truth 0.612 |
| vcrnapol, blank | 89.5% | 97.4% | 99.85% | 99.99% | 0.996; 4.2 / 4.1 Å (vs truth 0.996) | 0.997; 4.2 / 4.1 Å (vs truth 0.998) | 0.997; vs truth 0.998 |

What the student-only maps say. On the centred medium assembly the student's modes alone put 88.6% of all
particles within 6° of the truth and the cone search 97.1%; after one `cont` iteration the poses and the map are
the dense search's (98.9%, FSC 0.949 and 0.978), at a cone search of 135 s against the dense search's 231 s on the
same rows. On the uncentred half the student-only alignment falls far short (40.8% after `cont`, FSC 0.819 against
0.912): the student needs centred particles. On the real sets the student's orientations, drawn from the whole
training set's rows, give maps as good as or better than one dense search's wherever the particles are centred or
the set is large: trpv1 at the reference shifts (0.143 at 5.4 Å against 6.5 Å, 44.9% within 6° of the solve3D
direction against 42.5%), CESPED TEST at the shifts (0.143 at 9.2 Å against 10.3 Å; 48.2% against 42.3%), FlhB
(5.9 Å both; 33.6% against 29.4%), 10648 from a 900-row student (0.143 at 3.82 Å against 3.96 Å; 60.9% against
68.9%, the modes alone 35.7%; from an 18,000-row student 3.86 Å and 60.3%) and 11120 from an 18,000-row student
(4.39 Å both; 46.5% against 45.8%); on bgal at the shifts the modes alone agree with the solve3D direction for 56.4% of
the particles against the dense search's 35.0%, and the cone search then loses them (34.7%: the search inside a
wide cone leaves the student's direction), the maps equal. On the blank halves of the small real sets (bgal,
streptavidin) and on streptavidin at the shifts the student-only maps are no better than the dense search's, and
on uncentred particles everywhere the modes are worse than one dense search. The cost of the student-only route
is the prediction (0.5 to 1 ms per particle) plus the cone search, which with the wide cones the real sets'
students give is a near-exhaustive search over the 20,000-direction grid and costs 2 to 4 times the dense search
(trpv1 12 min against 3.4; FlhB 2.5 h against 40 min; 10648 33 to 38 min against 17; 11120 76 min against 39); on the centred medium assembly,
whose cones are narrow, it costs 0.6 of the dense search. The read decides nothing by itself; the maps are the
maintainer's. The small assembly's readout student (held-out mode agreement 20.8% against the invariant head's 2.0%, 0.6%
of the held-out seeds recovered, nobody seeded) gives student-only maps below the dense search's (uncentred
particles), and vcrnapol's (held-out mode agreement 78.2% against 64.9%, 98.6% of the held-out seeds recovered,
86% of the blank test half seeded at the dense search's accuracy, the neighbourhood search 28 s for them against
the dense search's 68 s for all) gives the dense search's maps from the modes alone and after `cont` at half the
cost (the cone search on all rows 32 s against 68 s).

**The cost table** (trpv1 at 15 Å, the test half of 17,752 particles at box 256, the training half of 17,752 for
the rows; the machine otherwise idle, 24 threads, `/usr/bin/time -v`; wall from three runs of the whole protocol
and of the reference, core time as user plus system seconds, peak memory as the maximum resident set). The
student's route: the rows of the training half (the coarse fill at the student's grid and band), the training of
the designs by the rule (16 and 32 channels here), the prediction on the test half, the neighbourhood search on
the seeded rows (none on trpv1) and the dense search on the rest, `cont`, the reconstruction; the dense route: the
dense search on every row, `cont`, the reconstruction.

| step | wall (three runs) | core time | peak memory | per particle |
| --- | --- | --- | --- | --- |
| the rows (the fill at 1,000 directions, 17,752 particles) | 39.7 s | inside the training's | inside the training's | 2.24 ms |
| the training (two designs, early stopping) | 961, 964, 948 s | 19,679 s | 2.7 GB | 60 ms wall, 1.2 core-s per training row |
| the prediction (17,752 particles) | 17, 17, 17 s | 345 s | 1.6 GB | 0.69 ms (features and forward), 0.96 ms with the reads and the writes |
| the dense search, every row (the reference) | 173, 204, 204 s | 4,320 s | 7.2 GB | 9.7 to 11.5 ms wall, 243 ms core |
| the dense search inside the protocol (the same rows, nobody seeded) | 205, 204, 205 s | 4,578 s | 7.2 GB | 11.5 ms wall, 258 ms core |
| `cont`, one iteration | 30, 29, 29 s | 486 s | 7.1 GB | 1.7 ms |
| the reconstruction | 23, 23, 22 s | 374 s | 1.9 GB | 1.3 ms |

The student-only route of ruling R4 on the same test half: the prediction as above and the cone search on every
row 12:05 wall, 16,512 s core, 19.7 GB (40.8 ms wall per particle: the wide cones make it a near-exhaustive search
over the 20,000-direction grid), against the dense search's 3:25; at the reference shifts 8:38 and 11,619 s. The
in-loop cycle at the curve's sizes on trpv1 (the rows, the training, the prediction, the search of the test half,
`cont` and the reconstruction): 180 rows 0.1 s + 80 s + 17 s + 194 s + 52 s; 900 rows 0.6 s + 2.5 min + 17 s +
194 s + 52 s; 4,500 rows 2.7 s + 6.8 min + 17 s + 194 s + 52 s; the whole half 40 s + 16 min + 17 s + 205 s +
52 s; against one dense iteration of 173 s on the same test half: at every size the student's cycle costs more
than the dense iteration it would replace, the training being the cost, and on trpv1 it seeds nobody.

**Exit.** The reads and the costs are above: the protocol with the folded cones on the medium assembly at 15 and
10 Å (ruling R1), the shift check and the modulus-only variant (ruling R2), the template readout with its suite and
the medium assembly's reads and curve points under it (ruling R3), the real sets under the readout (trpv1 with its
curve, bgal, streptavidin, CESPED TEST, FlhB at 20,000 rows, 10648's curve and 11120 at 20,000 rows), the
student-only read of ruling R4 on every test half with a trained student, and the cost table. The student is
cheaper than the coarse search it is trained on: the prediction's features and forward pass cost 0.5 to 0.9 ms per
particle against the fill's 2.2 ms (ruling R3, item 5), and on centred particles its cone search costs 0.6 of the
dense search; on the real sets' wide cones the cone search costs 2 to 4 times the dense search. No suite regressed:
the fast gate passes on the final tree (`scratch/logs/build_p3_readout.log`, 15 of 15; the Hermitian network 26
checks, the pose student 34, the polar features 39), the 1jxy identity held on the folded-cone tree
(`scratch/identity/1jxy_identity_p3.log`) and the readout touches no existing route (not rerun, as ruled). The
stricter variant (a de novo solve3D of the medium assembly's training half) was not run: the rulings' reruns took
the time. Left to the maintainer: the maps of the student-only read; the whole-half training's early stopping
(a tenth of the rows recovers more held-out seeds than the whole half on the medium assembly); the side scalar,
no longer read by the readout; the cone search's cost on wide cones against the dense search; the real sets'
teacher at 15 Å.

### Phase 4 (2026-11-08): the nightly route and the documents

**Files changed.** `src/main/commanders/test/simple_commanders_test_highlevel_pose.f90` (new: the `pose_student`
nightly route, a submodule of the high-level test family beside `_cont`),
`src/main/commanders/test/simple_commanders_test_highlevel.f90` (its commander type and interface),
`src/main/exec/simple_test_exec_highlevel.f90` (the router's case), `src/main/ui/simple_test/simple_test_ui_highlevel.f90`
(the test UI's help), `src/main/commanders/test/simple_test_registry.f90` (its registry row),
`production/CMakeLists.txt` (its CTest entry; the budget 42 to 43 with the reason),
`doc/policies/test_environment_policy.md` (the budget sentence, the entry table's row, its card; the file added to
the table), `doc/algorithms/pose_student.md` (new, the algorithm page) and `doc/algorithms/README.md` (its chapter;
added to the table), and this note (sections 5, 6, 7, 9 and 10 brought to what was built, the consumer-2 proposal,
this entry). `src/main/params/simple_parameters_tester.f90`, listed in the table for "the test names the parameters
tester lists", was not edited: it lists no test names (it binds one high-level command line, `cont_refine3D_1jxy`,
to check the suite's parameter binding), so a new entry needs nothing there; the registry check
(`scripts/check_test_registry.py .`: 38 registered selectors, 48 test programs, 48 router cases, 22 suite tables,
20 registry entries, consistent) is the gate that ties the seven files together. No new command-line key: the
route's programs run with the keys of phase 2 (`projfile`, `vol1`, `lp`, `mskdiam`, `nspace`, `pgrp`, `nthr`,
`fname`, `infile` and the existing `shift_sigma`) and the searches with `refine3D`'s own.

**The route** (`test=pose_student nthr=8`; the plan's section 10, tier 6; ruling R3, item 7). The split-half
protocol in miniature on the blob phantom of the multi-state fixture (`simple_test_state_fixture`: 2,000 simulated
particles of one asymmetric state at box 64 and 2 Å, SNR 1, no shift, their truth poses in the project; the fixture
is drawn fresh under the gate's seed). The odd rows are the training half, the even rows the test half (the
even/odd markers set in both `ptcl2D` and `ptcl3D`, as `reconstruct3D` requires). The training half's map is the
truth map soft-masked at `mskdiam` 100 Å. Stage `train`: `learn_pose_student` at `lp` 8 Å, `nspace` 1,000,
`shift_sigma` 2 (the default 10 px makes a 23 px disk on this box and the student learns nothing), both designs by
the rule (16 and 32 channels); the model's record is read back (threshold, bound, the held-out seeded and recovered
fractions, the training time). Stage `test`: the test half with blank poses, `predict_poses`, its table read back
(route, cone radius, confidence per row). The machinery's checks: the model records a threshold and a bound; every
test row routed; every seeded row's cone within the bound at or above the threshold; every dense row without a cone;
every folded cone radius at most 90°. Stages `search`, `control` and `dense`: one `refine3D refine=prob_neigh
prob_neigh_mode=geom` iteration on every test row (the cone where there is one, the coarse cell where there is none),
the same search on a copy with the student's records zeroed, and the dense search (`refine=prob`) on every row as the
protocol's fallback; the classic generator is reset to the gate's seed before each in-process search (each
`refine3D` reseeds it from the call count, so without the reset two in-process runs draw differently; the keyed
draws share the run's seed). The checks: every searched pose finite; the rows without a cone bit for bit the
no-student search's (the base path's identity inside the route); the seeded rows changed by their cones (reported,
not checked, when nobody is seeded). Stage `merged`: the seeded rows' poses from the cone search into the dense
project, `reconstruct3D` of that merge and of the dense search alone, the maps' FSC against the truth (both volumes
Fourier-transformed first: `image%fsc` reads Fourier components). Everything else is reported with no floor, by
R3's item 7: the held-out record, the seeded fraction, the truth inside the cone and the mode within 12° for the
seeded rows, the agreement with the truth within 6° and 12° per route beside the dense search's on the same rows,
the maps' FSC, every step's wall time. Every row goes to `metrics.tsv`; a profile sets floors later.

The route was first written on 1,000 simulated 1JYX particles at box 144 (the `cont` entry's set). At `c1` every
posterior of that D2 tetramer spreads over its four equivalent directions, so no folded cone is narrow and the
student seeds nobody; the gate's seeded-row checks would then never run. The blob is asymmetric and the route moved
to it. Two more corrections on the way: the first draft's rows without a cone took the coarse cell from a blank
pose (1.2% within 6° of the truth), not the protocol's dense fallback, so the dense search, the merge and the two
reconstructions were added; and the FSC read all zeros until the volumes were transformed.

**The measured run** (`scratch/pose_route/route.log`, `metrics.tsv`, `route.time`; by hand at 8 threads under
`/usr/bin/time -v` with a disk watcher, the machine otherwise idle): PASS in 6:49 wall, 2,444 s of core time,
peak 1.6 GB, 0.72 GB of run directory. The student: 16 channels held-out loss 2.733, mode agreement 68%, 3% of the
held-out seeds recovered (158 s); 32 channels 2.691, 71%, 3% (230 s); 16 kept, threshold 0.991, bound 10°; the
training of the two designs 394 s. Prediction 1.2 s for 1,000 particles (1.2 ms each). 15 of the 1,000 test rows
seeded (1.5%), the truth inside the cone for 40% of them and the mode within 12° for 67%; after the cone search
26.7% of the seeded rows within 6° of the truth (66.7% within 12°) against the dense search's 33.3% on the same
rows; the dense search on every row from blank 3.8% within 6° (the blob's views are alike in the 8 Å band) while
the maps reach FSC 0.143 at 8 Å (the band's edge) and 0.78 at 12 Å against the truth (the dense map 0.77). The
cone search 3.2 s, the control 3.1 s, the dense 7.1 s, the reconstructions 0.2 s. All eight checks hold. The
registry row and the CTest entry carry the measurement: deadline 2,400 s, 8 threads, 8 cores, 2.0 GB, 1.0 GB of
disk; the policy's card records the claim, the checks, the path and control, no profile yet, and the cost. The
budget is 43 (15 fast, 7 library, 20 highlevel, 1 platform); `check_descr.py` passes; the fast gate 15 of 15 on the
tree (`scratch/logs/build_p4.log`, `scratch/logs/fast_gate_phase4.log`).

**The documents.** `doc/algorithms/pose_student.md`: the problem, the model (the input tensor, the Hermitian
blocks, the template readout), the training (the teacher's rows, the augmentation, the designs and the held-out
rule, the threshold), the prediction and the hand-off (the fold, the route, the cone of the neighbourhood search),
what the reads found (phases 1 to 3 in plain words, the shift check among them), the rationale and the
implementation;
listed as chapter 16 of `doc/algorithms/README.md`. This note: section 5 (the augmentation's disk as built and the
2 px width the measurements used; ruling R2's finding), section 6 (the readout in place of the invariant head, the
model file's version 2, `nhid` gone), section 7 (recovery and the threshold on the folded sphere), section 9 (the
consumer-2 proposal from the measured cycle: not warranted by the measured cycle as the programs stand; a candidate
only where the particles are centred by an earlier stage, the teacher is sharp and the training is held to the
cheapest design at a few thousand rows; the decision is Hans's; the sampling and fractional-update policy unchanged),
section 10 (tier 6 as built). The skill text (`simple-main-learn`, and a line each for the strategies and the
probabilistic search guides) is proposed in phase 5's report, not edited. The note's row of
`doc/implementation_notes/README.md` still says "nothing built"; it changes when the note moves to `completed/`
(phase 5, the file table).

**`test=nightly` on the final tree** (`test=nightly ncunits=24`, 12:58 to 13:44 on 2026-11-08, the machine
otherwise idle; `scratch/nightly/nightly.log`, the report `scratch/nightly/nightly_report.tsv`, 1,169 rows): every
one of the 20 entries passed, no row failed, no base failure. The walls: `solve3D_addon` 118 s,
`cont_refine3D_1jxy` 329, `flex_pca_blobs` 367, `simulate_particles` 12, `pose_student` 550, `state_fixture` 6,
`species_discovery` 147, the four `single_atoms_stats` 25 to 33, `mini_stream_1jxy` 57, `mini_stream_6vxx` 52,
`stream_chains` 720, `stream_preproc` 82, `single_workflow_fcc` 602, `simulated_workflow_1jxy` 933,
`multistate_workflows` 2,521, `simulated_workflow_6vxx` 847, `single_workflow_wurtzite` 686; the tier 2,790 s.
The `pose_student` entry's reports inside the nightly equal the by-hand run's to the digit (the gate's seed: 1.5%
seeded, 26.7% within 6° after the cone search, FSC 0.78 at 12 Å), only the times differ on the shared machine
(the training 529 s, the prediction 1.6 s, the searches 4.3, 4.2 and 9.6 s; 550 s in all against the 2,400 s
deadline; peak 8.0 cores of the 8 reserved, 1.49 GiB). After the nightly, three strings from the route's first
draft (the route file's description line and the test UI's two help lines, which still said 1JYX) were corrected,
the tree rebuilt and the fast gate rerun (`scratch/logs/build_p4_final.log`): no behaviour changed.

**Exit.** The nightly is green on the final tree, with no base failure to record; the documents are written (the algorithm page, the policy's card, this
note's sections). Left to the maintainer: a 20-seed profile of the route's reports before any floor (the policy's
rule; the seeded fraction on the blob, 1.5%, is the number a floor would guard); the decision on consumer 2.

### Phase 5 (2026-11-08): the report, and ruling R5

**Files changed.** `src/main/learn/simple_circ_cnn.f90` and `simple_circ_cnn_tester.f90` moved with `mv` to
`simple_hhcnn.f90` and `simple_hhcnn_tester.f90` (ruling R5: the module `simple_hhcnn`, the submodule
`simple_hhcnn_tester`, the type `hhcnn`, `run_all_hhcnn_tests`; the constants `CNN_*` and the weights file's magic
number stay, as the other names do); `src/main/learn/simple_pose_student.f90` (the renamed import and type);
`src/main/commanders/test/simple_commanders_test_class.f90` (the suite `HHCNN` of `unit_numerics`);
`src/main/ui/simple_test/simple_test_ui_class.f90` (the suite key `hhcnn` in the `suite=` list, which the registry
check caught); `doc/algorithms/pose_student.md` (the network named in one sentence in the Model section, with the
contrast to Worrall's Harmonic Networks; the module's name under Implementation); this note (sections 1 and 6 name
the network and define it, the phase-1 line of section 11, decision 1 of section 12 and the module list of section
16 note the rename, the file table carries the moved files and this phase's rows, this entry; then moved to
`doc/implementation_notes/completed/` with `mv`); `doc/implementation_notes/README.md` (the note's row moved from
the open work to the completed investigations, pointing at the report and the completed note);
`doc/implementation_notes/completed/learned_pose_state_assignment_report.md` (new: the review report). The earlier
Progress entries keep the names of their day (`simple_circ_cnn`, the suite `Hermitian network`).

**Ruling R5 applied** (Hans, 2026-11-08: the network is HHCNN, the Harmonic Hermitian convolutional neural network:
Hermitian because every layer commutes with conjugation on the half ring with the conjugate wrap, so the network is
a model of the Fourier transform of a real image; harmonic because the readout's class templates live in the
angular harmonic domain of the polar transform and the rotation search is the cross-correlation over the cyclic
shifts through the angular Fourier transform; convolutional on the polar Fourier representation). The rename is a
rename: no behaviour, constant or file format changed. After it the tree was reconfigured (the glob of
`src/CMakeLists.txt` picks the moved files up), built and installed, the registry check is consistent, the fast gate
passes (`scratch/logs/build_p5b.log`, 15 of 15) and the three learn suites alone at 16 threads pass
(`scratch/unit/p5/unit_*.log`: polar features 39 checks, HHCNN 26 checks in 5.8 s with its thread identity, pose
student 34 checks). After the plan's move the five comments that pointed at its old path (the algorithm page, the
network, the features, the strategy, the nightly route) were repointed at `completed/`, and the final tree was
rebuilt with the gate (`scratch/logs/build_p5c.log`, 15 of 15) and the suites (`scratch/unit/p5c/`) green again. The first build after the rename (`scratch/logs/build_p5.log`) failed the gate's registry check
on the test UI's stale suite key `hermitian_network`; the key is `hhcnn` now and that file's row of the table
carries this phase.

**The report.** `doc/implementation_notes/completed/learned_pose_state_assignment_report.md` (340 lines): what the
change is, what was learned in plain words (the routes untouched, the teacher not the cost, the two heads and the
readout, the real sets' broad teacher, how few particles, the centring, the student-only maps, the costs), what
changed phase by phase with the diffs' paths, the decisions taken under delegation and where each is recorded, the
deviations from the plan and the base behaviours recorded, the evidence against the baseline, the measurements
(the protocol's reads per set, the cost table, the in-loop cycle), the skill text proposed, the in-loop consumer's
proposal, what is left to the maintainer, and the commit message. The report was cross-checked against a digest of
the Progress entries, the status files and the rulings (`scratch/report_digest.md`).

**Exit.** The report is written beside this note; this note is moved to `completed/` and its row of the README
updated; ruling R5 is applied with the gate and the suites green; no job of the run is running; nothing is
committed. Journal: `state/journal_phase_5.md`.

### Mirror pairs (2026-11-09, in the maintainer's tree)

Done by hand on the Dell in `~/src/SIMPLE`, after the run's code was brought into the working tree on
top of HEAD b2154da6f (Cyril's opt-in mirrored reference extraction, `mirr_proj=yes`), uncommitted.

**Why.** Cyril found that generating the second half of the reference spiral by mirroring the first
accelerates reference preparation (`mirror_ref_bank`, the frame angles of `oris%get_mirror_pairs`).
The student was trained over all 1,000 directions, so the network had to tell a view from its
mirror, a distinction the search makes exactly with a cross-correlation of the two candidates. Hans
ruled that the mirror correspondence be used through the Euler relation in the search space: the
student classifies the mirror pairs and leaves the pair to the search.

**A first attempt, not kept.** Two classes per pair with the partner's template tied to its
partner's mirror in the polar Fourier domain, `T_m(n) = conj(T_i(n)) exp(-i n psi)` (the
angular-harmonic form of `mirror_ref_bank`), and symmetric kernels. Without the frame angle the tie
blurred each template against a rotated copy of itself (the spiral's mirror partner is the reflection
rotated by an off-grid angle), and the 200-particle student fell from 69% to 30% within 6°; with the
frame angle it was level with the untied student at 200 and 5,000 particles and worse at 1,000. It
was the wrong level for the relation and was replaced by the pair classes.

**What changed.**

- `simple_hhcnn.f90`: the kernels are symmetric along the angle from construction (`impose_symmetry`)
  and stay so through training (`tie_kernel_grads` sums the two offsets' gradients onto both); the
  readout scores each template against the map and its reflection (the conjugate spectrum), the
  larger, and records which won (`refl`) for the backward pass.
- `simple_pose_student.f90`: the network has `nspace/2` classes, one per mirror pair from
  `oris%get_mirror_pairs` (checked against the grid's `mirr` key); the teacher's row is folded onto
  the pairs in the batch assembly; a prediction places each pair's mass on its unmirrored direction,
  so the interface over the grid, the cone, the strategy and `prob_neigh` are unchanged; the
  templates start from the unmirrored directions' reprojections. Model version 3.
- `simple_strategy_pose_student.f90`: the prediction's diagnostic rank of the reference direction is
  taken over its pair.
- Tests: `test_mirror_invariance` in the HHCNN suite (a sample mirrored by `mirror_ref_bank` at 0, 180
  and 75 degrees, a general frame on the rotation grid, scores against every class as the sample, to
  1e-4 of the largest score; measured 7.7e-8 at 0 and 180; the kernels symmetric bit for bit after
  training); the rotation test expects a backward step for a class the reflection scored; the toy
  problem classifies the spiral's 10 mirror pairs (held-out accuracy 1.000); `test_pairs` in the
  student's suite (all of a prediction's mass on the unmirrored directions). A frame angle off the
  rotation grid (73.7 degrees) changes the scores by 4.5%: the network is rotation-invariant on the
  grid's steps only, as for any particle. Fast gate 15 of 15; `unit_numerics` 1,926 assertions pass.

**The size curve, pairs against the 1,000-class student** (the run's final build), the medium
assembly at 15 Å, the same subsets, seed, 2 px shift augmentation and test half, centred particles
(the test half at its truth shifts), one seed per point; folded reads in both (`pred_stats.awk`):

| training particles | within 6° of truth, 1,000 classes | pairs | truth in the cone, 1,000 | pairs | cone size q90, 1,000 | pairs | over the cap, 1,000 | pairs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 200 | 69.2% | 69.6% | 93.9% | 94.9% | 64 | 64 | 14.1% | 18.9% |
| 1,000 | 81.7% | 85.5% | 98.8% | 99.0% | 64 | 64 | 24.7% | 14.1% |
| 5,000 | 87.9% | 85.6% | 97.8% | 98.2% | 30 | 58 | 4.9% | 9.1% |
| whole half (24,000) | 91.5% | 94.0% | 99.6% | 99.6% | 42 | 14 | 7.3% | 3.0% |

On the whole half the pair student also puts 99.5% within 12° (98.6%), and 97.0% of the particles
clear the confidence 0.99 (92.7%) with the truth in the cone for 99.7% of them. The 5,000 point goes
the other way; each design chose its width there separately, and one seed per point does not
separate that from noise. The reads' "mode = truth grid point" line is not comparable between the
two (the pair student always names the unmirrored direction), nor the rank (pairs against grid
points). Logs and runs: `/data2/agent_runs/mirror_tie/study_pairs.log`; the polar-tie attempt's in
`study_polar_tie.log` and `study_first_pass.log`.

**Open:** repeat the curve over three or four seeds before reading a trend at the small sizes; the
real sets and the student-only reconstructions of phase 3 have not been rerun with the pair student.

### Validation workflow and benchmarks (2026-11-09, in the maintainer's tree)

Done by hand on the Dell in `~/src/SIMPLE`, uncommitted, at Hans's request: a solve3D-driven validation
program and detailed benchmarks, to judge where the student belongs in production code.

**What was added.**

- `validate_pose_student` (simple_exec, developer visibility): solve3D's inputs, its `nsample` the
  size of a class-balanced greedy training sample of the 2D solution. solve3D runs on the sample with
  `nsample` at its size (full sampling), the student is trained on solve3D's poses and map (15 Å, a
  2 Å shift width: the particles are centred), and every selected particle is predicted. The report is
  described in `doc/algorithms/pose_student.md`, section "Validation and benchmarks". Controller in
  `simple_commanders_pose_student.f90`, UI in `simple_ui_solve3D.f90`, route in `simple_exec_solve3D.f90`.
- Bench records `POSE_STUDENT_LEARN` (rows, then per design: batch assembly, gradients, steps,
  held-out evaluation; the fill profile) and `POSE_STUDENT_PREDICT` (phases; per-particle steps in
  thread-seconds; utilisation and load balance), and `VALIDATE_POSE_STUDENT` (the phases).
- `test=pose_student_bench` (simple_test_exec, no CTest entry, test policy section 1.3): the HHCNN
  kernels on synthetic tensors (`simple_hhcnn_bench.f90`, a submodule of `simple_hhcnn`, with the
  network's `hhcnn_profile` counters).

**First read: bgal** (`/data2/bgal_addon/3_selection`, one stack, 2,756 selected particles in 50
classes, D2, box 256 at 1.275 Å, mask 180 Å). 1,500 training particles, 1,256 unseen; 24 threads.
Two runs of the same command (`run2`, `run3` in `/data2/agent_runs/validate_ps/`); solve3D took a
different stage path in each (multi-threaded, unseeded), so the teachers differ.

| read | run 2 | run 3 |
| --- | --- | --- |
| held-out rows seeded, recovered | 37.3%, 98.2% | 45.3%, 97.1% |
| training particles: mode within 6° of the teacher | 87.9% | 89.9% |
| the same, seeded training particles | 98.9% | 97.5% |
| training particles seeded | 49.2% | 70.3% |
| unseen particles seeded | 3.7% | 6.1% |
| distance to the 2D class's teacher medoid, median: the teacher's own (baseline) | 6.86° | 4.44° |
| the same, student on training particles | 6.64° | 4.49° |
| the same, student on unseen particles | 6.67° | 4.74° |
| unseen within 12° of the medoid (baseline: training particles' teacher poses) | 73.4% (68.8%) | 77.1% (72.8%) |
| teacher's shift minus the 2D shift (training), median and q90 | 1.6 px, 16.8 px | 1.9 px, 17.6 px |
| training rows within the training's shift disk (4.7 px) | | 71.3% |

The unseen particles' directions agree with their classes as well as the teacher's own poses do,
but few are seeded: their cones are wider and their confidence lower than the held-out training
rows'. The likeliest cause is the shift: the unseen rows carry their 2D shifts, and the teacher moved
a tenth of the training rows by more than 16.8 px from theirs, against a training disk of 4.7 px. A
first run with the default shift width (10 Å, a 23.5 px disk) learned almost nothing (held-out mode
agreement 7%, nobody seeded): the width matters more than any other setting seen so far.

**Cost** (run 2; run 3's solve3D ran one more stage, 604 s): solve3D 344 s (229 ms per training particle), training 178 s (two designs,
119 ms per training particle; the parallel gradients 91% of it), prediction 4.2 s for 2,756 particles
(1.5 ms each; ruling R2's shift check on 256 particles is half the prediction's thread time; the forward
pass alone is 12.6 thread-ms per particle, 0.5 ms at 24 threads).

**Kernels** (`test=pose_student_bench`, two runs, 24 threads, 160-angle half ring, 500 classes; rows
in `/data2/agent_runs/validate_ps/bench2/BENCH/`): the forward pass takes 4.3 to 18 ms per sample at
one thread, the readout 2.6 to 5.5 ms of it, so the readout dominates at 16 channels up to 25 rings; a
training step of 256 samples takes 140 to 2,400 ms at the team, 8 to 13 times faster than at one
thread (the one-thread step is timed once and moved by up to 50% between the runs); the step is
bit-identical at one thread and at the team in every row.

**Open:** the shift gap between the 2D solution and the teacher (a training disk that covers it, or
centring the unseen rows in the teacher's frame) before the seeded fraction of unseen particles is
read as the student's; a larger project; the readout's cost.

**Poses straight from the student (the same day, Hans's question).** Can the unseen particles be
reconstructed from the student's orientations with no refinement at all? The readout already finds
the in-plane rotation and the mirror (it maximises each class's score over the rotation grid and the
reflection) and threw them away; `pose_student%pose_of` now reads them out for the mode's pair
(`doc/algorithms/pose_student.md`, "The readout's pose"), `predict_poses` writes them as the table's
last three columns, and the validation reconstructs four maps with reconstruct3D. The network does
not estimate shifts: the unseen particles keep their 2D shifts. Run 4 (`/data2/agent_runs/validate_ps/run4`,
the same command as runs 2 and 3; this time the 64-channel design was kept):

| readout pose against the teacher, training particles | median | q90 | within 6° |
| --- | --- | --- | --- |
| direction (mirror as read out, point group folded) | 2.7° | 21.1° | 80.5% |
| in-plane angle, direction within 12° | 2.0° | 9.4° | 84.1% |
| direction, seeded particles | 2.1° | 4.3° | 95.6% |

| map (reconstruct3D, D2, no search) | particles | FSC 0.5 / 0.143 against the teacher's map | halves 0.5 / 0.143 |
| --- | --- | --- | --- |
| unseen: student's poses, 2D shifts | 1,256 | 10.2 / 7.6 Å | 13.6 / 8.8 Å |
| training: student's poses, 2D shifts | 1,500 | 9.6 / 6.4 Å | 12.6 / 8.8 Å |
| training: student's poses, teacher's shifts | 1,500 | 6.8 / 4.1 Å | 9.9 / 6.9 Å |
| training: teacher's poses and shifts (control) | 1,500 | (same particles) | 7.6 / 4.4 Å |

The unseen map against the teacher's is cross-validated: the two share no particle, and they agree to
10.2 Å at 0.5. The unseen and the training particles give the same half-map resolution from the
student's poses, so the student generalises: what it loses against the teacher is not the unseen
particles' novelty. The shifts cost more than the orientations: with the teacher's shifts the
student's poses reach 9.9 Å at 0.5 between halves, against 7.6 Å for the teacher's poses; with the 2D
shifts, 12.6 Å. The control's FSC against the teacher's map (2.55 Å, Nyquist) only says the two are
the same reconstruction.

**Open (added):** a shift read-out or a cheap shift estimate for the unseen particles; the readout's
in-plane angle is on the rotation grid (0.8° steps at this 220-angle half ring), finer than its measured error.

**Shifts registered at the read-out pose (the same day, Hans's option A).** The network tolerates
shifts but does not estimate them. Given the map, `predict_poses` now fits each particle's shift
against the map's reprojection at the student's pose, direction and rotation fixed, with the search's
fixed-rotation shift search (a 5 x 5 grid over 25 Å, then L-BFGS-B). Run 5
(`/data2/agent_runs/validate_ps/run5`, the same command as runs 2-4):

| shift against the teacher's, training particles | median | q90 | within 3 px |
| --- | --- | --- | --- |
| 2D shift | 1.59 px | 17.87 px | 64.9% |
| registered at the read-out pose | 0.59 px | 1.62 px | 96.6% |
| registered, direction within 6° | 0.54 px | 1.10 px | 99.5% |

| map (reconstruct3D, D2, no orientation search) | particles | FSC 0.5 / 0.143 against the teacher's map | halves 0.5 / 0.143 |
| --- | --- | --- | --- |
| unseen: student's poses, 2D shifts | 1,256 | 9.1 / 5.5 Å | 10.5 / 7.8 Å |
| unseen: student's poses, registered shifts | 1,256 | 7.8 / 4.7 Å | 10.2 / 7.6 Å |
| training: student's poses, 2D shifts | 1,500 | 7.8 / 5.1 Å | 10.2 / 7.8 Å |
| training: student's poses, registered shifts | 1,500 | 5.5 / 3.9 Å | 9.1 / 6.0 Å |
| training: student's poses, teacher's shifts | 1,500 | 4.7 / 3.9 Å | 8.6 / 5.1 Å |
| training: teacher's poses and shifts (control) | 1,500 | (same particles) | 6.5 / 4.1 Å |

The registration recovers the teacher's shifts to about half a pixel, and closes most of the gap to
the teacher's-shift ceiling on the training particles (9.1 against 8.6 Å at 0.5). On the unseen
particles the cross-validated agreement with the teacher's map improves from 9.1 to 7.8 Å, their own
halves only from 10.5 to 10.2 Å: their orientations, which no teacher measures, are now the limit. It
costs 0.5 thread-ms per particle against 22 for the forward pass. This run's maps are better than run
4's throughout (its teacher differs; solve3D is unseeded), so compare rows within a run, not across.

**Open (added):** the shift fit at the 15 Å band; a finer band for it alone; the unseen particles'
orientation error, which only a teacher pass over them would measure.
