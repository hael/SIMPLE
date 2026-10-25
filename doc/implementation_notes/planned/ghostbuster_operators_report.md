# Ghostbuster's operators in SIMPLE: report of the run

Run "wave_operators", 2026-10-23 to 2026-10-24, unattended on the Dell (Intel Core Ultra 9 285, 24 cores, no
GPU), against base `5111068bf`. The plan is `ghostbuster_operators.md` beside this report; its section 8,
"Progress", holds the detail and the evidence of every phase. Paths of the evidence are under the run directory
`/home/elmlundho/agent_runs/wave_operators`; each phase's diff is `review/phase_N.diff` there. Nothing is
committed: the change is the working tree against the base.

## What the run delivered

SIMPLE can now predict a particle image from a 3D electrostatic potential with wave optics, under four models
(projection, first Born and Rytov on the Ewald paraboloid, and multislice), and has the exact, hand-written
derivative of that prediction. It is built on SIMPLE's own machinery: the Kaiser-Bessel (KB) gather and scatter
of the PCG (preconditioned conjugate gradients) reconstruction's padded lattice, SIMPLE's CTF, its particle read
path, its PCG start solve, density and maximum-likelihood (ML) prior. specter, the package that contains
Ghostbuster, served only as the oracle that wrote reference files; nothing in the repository calls it.

On top of the operators:

- **`reconstruct3D_wave`**, a new shared-memory program, reconstructs half maps from a project with known poses
  and CTF under any of the four models, by nonlinear conjugate gradients from the flat PCG map.
- **`simulate_scattering`** now draws its images with the same operators, with corrected sign conventions and
  the defocus at the specimen midplane (section "simulate_scattering" below).

What it buys, on simulated data with a known answer:

- On a 300 A ribosome at 100 kV, where the Ewald curvature limits a flat map beyond 4.7 A, Rytov with the
  data's curvature sign raises the agreement with the truth between 3 and 2.4 A (the Fourier shell correlation,
  FSC, of the map against the truth) from 0.66 to 0.81; the wrong sign drops it to 0.21. specter's converged
  reconstruction reaches 0.83.
- On an ordinary protein at 300 kV all models agree, Rytov slightly ahead (0.762 against 0.733 between 3 and
  2 A), level with specter.
- On a 469-atom Pt nanoparticle every map finds all 469 atoms; the wave models place them more accurately
  (0.029 to 0.030 A against 0.035 A) and repair part of the band near the particle's own curvature limit (about
  1 A). specter's own multislice matches SIMPLE's over the last band (0.801 both).

Cost: a Rytov gradient pass costs 1.7 to 3 times a PCG accumulation; the ribosome's wave reconstruction took
about 9 min against 1.5 min for PCG at 24 threads. Multislice costs about 3 s a particle on one core at box 160,
which suits nanoparticle sets, not 100 000-particle cryo-EM sets.

## Phase by phase

- **Phase 0, oracle and conventions** (`review/phase_0.diff`: the plan only). specter installed in scratch and
  its tests for the parts SIMPLE reproduces passed (one oversampling-3 rounding failure, explained). The
  SIMPLE-to-specter convention table measured through RELION star files (blob positions within 0.02 pixel, CTF
  zeros within 0.0001 1/A). Reference cases and three ground-truth datasets (protein, ribosome, Pt
  nanoparticle) simulated by specter's multislice with their truth. Both findings about `simulate_scattering`
  (inverted phase, exit-face defocus) confirmed with numbers.
- **Phase 1, foundations** (`review/phase_1.diff`). The plane gather and scatter moved from
  `reconstructor_pcg` into `pcg_lattice` with the curvature; the reconstruction's observation-batch frame
  extracted from the PCG worker into `simple_matcher_ptcl_io`; the wavelength and interaction constant given one
  home (`simple_math_ctf`); the CTF's phase made public (`eval_phase`); the Poisson sampler moved to
  `simple_rnd`; the new `wave_field` (2D complex fields) and `nlcg_solver` (nonlinear conjugate gradients). PCG
  results byte-identical to the base.
- **Phase 2, the single-scatter operator** (`review/phase_2.diff`). `wave_imaging` (the 2D chain of one particle)
  and `wave_projector` (the operator over a particle set) for projection, first Born and Rytov, forward,
  linearisation and adjoint, threaded and reproducible; the fold of the scatter accumulator moved into the
  lattice. The parity campaign against specter, the oversampling measurement and the cost table.
- **Phase 3, the program** (`review/phase_3.diff`). `reconstruct3D_wave` end to end (keys, UI, router,
  commander, strategy and its objective); the PCG worker's accumulation of a state's half extracted for sharing;
  `simulate_scattering`'s single-scatter models on the shared operator, its keys renamed and retired; the
  library test; the ground-truth comparison on the protein and the ribosome.
- **Phase 4, multislice** (`review/phase_4.diff`). `wave_multislice` and the cube gather and scatter of the
  lattice; `wave_model=multislice` in the program; `simulate_scattering`'s multislice on the shared operator and
  `simple_scattering_simulator.f90` deleted; the ground truth on the nanoparticle with SINGLE's atom detection.
- **Phase 5, documentation and this report** (`review/phase_5.diff`). The wave-model policy
  (`doc/policies/3D/wave_model_policy.md`); the user-guide page `doc/user_guide/reconstruct3D_wave.md` and the
  updated `simulate_scattering.md` and index; the new sub-suites in `doc/policies/test_environment_policy.md`;
  `reconstruct3D_pcg_policy.md` pointed at the moved routines; the note of `simulate_scattering`'s author, Ruben, updated; a particle count corrected
  in the plan's Phase 4 entry (2 000 Pt particles, 1 000 per half, not 1 000).

## simulate_scattering: what changed, for its author

Ruben wrote `simulate_scattering` on 2026-10-23 (commit `f4cac4773`); this section is for him. Ruben, your
simulator was the starting point of the 2D wave engine; the run moved it onto the shared operator.
Your note `ghostbuster_scattering_simulator.md` carries the full account in a dated section; in short:

- **The images were those of `-V`.** The propagator `exp(+i pi lambda d k^2)` and the lens `exp(+i chi)` were
  both conjugated relative to Kirkland while the transmission was `exp(+i sigma V)`, so a positive potential came
  out bright at underfocus and the multiple-scattering terms were those of the negated potential. Against specter
  on the same potential and pose (Phase 0), the simulator's image of `V` correlated with specter's image of `-V`
  at 0.98 to 1.000 and with that of `V` at -0.90 to -0.93. The absorption had the right sign. Kirkland's signs
  are now used throughout.
- **Defocus referred to the exit face**, half a box below the midplane that SIMPLE's CTF fits and specter use
  (the best match was 48.6 to 51.4 A off for a 48 A half depth). It now refers to the midplane.
- **Keys.** `scattering_model` became `wave_model` (with `projection` and `rytov` added), `ewald` selects the
  curvature sign or, for multislice, the beam direction, and `absorption` is retired for `fraca`: with the object
  potential `(sqrt(1 - alpha^2) + i alpha) V` the weak-phase limit of every model is SIMPLE's CTF with amplitude
  contrast `alpha`, tested to `1e-5`, so they are one quantity. No aliases (release 4).
- **No real-space rotation.** Every model samples the potential's transform on the PCG lattice; multislice
  gathers the rotated spectrum cube from it. The data are made by the exact operator the reconstruction inverts.
  `rotvol` had not corrected the KB window's roll-off (8 to 10% weaker off the centre).
- **Lens and geometry.** The lens is SIMPLE's CTF phase, with a phase plate acting on the scattered waves only;
  shifts and phase shifts of the orientation table now apply.
- **Code and tests.** `simple_scattering_simulator.f90` and its tester are gone: the FFT buffers became
  `wave_field`, the lens and detector `wave_imaging` (`detect_intensity` holds your Poisson and Gaussian noise),
  multislice `wave_multislice`, the Poisson sampler `poisson_rnd`. Your multislice limits (zero, constant real
  and constant absorptive potentials) are kept; the first-Born constant-potential check gave way to the
  weak-phase identity and to first Born being Rytov's linearisation, and the vacuum is tested for every model.
  New tests can see the sign and the defocus reference (a weak sheet 16 A downstream images at `df + 16 A`, not `df - 16 A`).
  The FFTW plans now use SIMPLE's one planner lock (the unnamed critical section).

The PDB route, the potential route, the dose and noise, the defaults and the outputs are as you made them. The
program and its tests ran on Linux (Oracle Linux 8, gfortran 15) throughout this run.

## Decisions taken under delegation

Each is recorded, with its reason, in the plan section named.

- The FFTW planner lock is the unnamed critical section, not a named one (section 4.2, Phase 1).
- The window helpers, `plane_loc` and `rim_free` are module procedures of `pcg_lattice`, not type-bound, to
  leave `flex_pcg_t` unchanged (Phase 1).
- `ewald_curv` refuses a paraboloid more than 0.49 padded voxel past the Nyquist sphere (Phase 1).
- Rytov exponentiates at the exit face, as specter and the approximation's definition do (section 3.2,
  Phase 2); the midplane differed by 2 to 17%.
- A phase plate's shift leaves the unscattered beam alone, `T(0) = 1` (section 3.2, Phase 2).
- The wavelength stays SIMPLE's CTF value, shared by design (section 3.2, Phase 2).
- Oversampling stays 2 with no key; the wave operator keeps the PCG lattice's three-tap window, and the existing
  key `winsz` serves a caller who wants seven taps (Phases 2 and 3).
- specter's positive curvature sign is one voxel of defocus off the symmetric reference; SIMPLE keeps the
  symmetric one (section 3.1, Phase 2).
- The program's design details (section 3.6, Phase 3): the prior in volts carries `boxpd^3`; the start is the
  base PCG solve shrunk by the ML prior in closed form, written beside the wave maps; the Wiener filter is undone
  per shell; the step subset is 5% of a half within 256 to 1000 particles; the repeated passes use the in-memory
  pass buffer; `pcg_mskfile` is the support key; halves from `ensure_field_halves`.
- `ewald` defaults to `no` for projection and `neg` otherwise; projection with a curvature and multislice without
  one are refused (Phases 3 and 4).
- `wave_multislice` is a per-thread type run by `wave_projector` as its fourth model, not a second extension of
  the lattice with its own batch frame (section 4.1, Phase 4); its constant is `WAVE_MSLICE`, since Fortran is
  case-insensitive and `WAVE_MULTISLICE` would name the type.
- The ground-truth sets render each particle from its rotated atoms, are normalised by the dose so the gain is
  exact, and ship SIMPLE's orientation table beside their star files; the ribosome set propagates unpadded
  (Phase 0).
- specter reconstructs Rytov with its surface sampling (its real-space rotation would have taken 8 to 33 h);
  for multislice its real-space rotation is the only mode (Phases 3 and 4).
- The ML prior is left as the plan fixes it, from the flat start pair, so that Phase 3's results stand (Phase 4);
  see "Open".
- This phase added `doc/policies/3D/reconstruct3D_pcg_policy.md` to the file table, since the routines it names
  moved in Phase 1 (the item Phase 1 left for Phase 5).

## Deviations from the plan

- **The lattice's accuracy is ten times the plan's estimate.** Section 6 assumed about `1e-3`; the three-tap
  window interpolates a section to 0.7 to 1% of the image contrast for content off the box centre (Phases 1 and
  2). The parity tolerance was measured at zero curvature, where both packages are exact by another route
  (8.2e-3 of the image contrast, 1.3e-2 of a gradient), and the golden test's bound set at 1.25 times that; a
  seven-tap window brings the operator to 4e-5 to 9e-4. Far below any particle's noise.
- **specter at general poses under multislice.** specter's multislice exists only with its real-space trilinear
  rotation, whose own error (12.6 to 19% against its surface sampling) bounds parity at general poses; parity
  was established at quarter turns (to 7e-5 and 2e-3) and the gradient's 30% difference traced to specter's
  real-space gradient (Phase 4).
- **The orientation test's criterion** is a factor of four between the wrong and the right sign, not ten: the
  two signs differ by only 2% at box 24 (Phase 4, recorded before relied on).
- **The particle cache is not used** (it is two-dimensional only); the pass buffer serves the repeated passes
  (section 4.3).
- **The nanoparticle set's expectation** ("multislice recovers atoms the linear map misplaces or misses") did not
  arise: at this dose every map finds every atom; the difference is in position error and in the band near the
  curvature limit (Phase 4).
- **`simulate_scattering` kept `rotvol` for multislice between Phases 3 and 4**, then lost it as planned.

## Evidence against the baseline

- **PCG unchanged.** A fixed case (400 myoglobin particles, `reconstruct3D rec_backend=pcg`, point groups c1 and
  d2, objectives `cc` and `euclid`) gives byte-identical maps, FSC and sigma2 files against the base build after
  every phase (`scratch/identity_p4.log`, earlier `identity_*.log`); the 6 000-particle protein set likewise at
  24 threads (Phase 1). `pcg_recon` reports all 42 gate metrics identical; PCG's accumulation on the protein set took
  6.85 s after the move against 6.72 s at the base, within the run-to-run spread of 0.5 s
  (`scratch/timing_pcg/times.txt`).
- **Suites.** `unit_image` (925 assertions), `unit_reconstruction` (440), `unit_numerics` (1365),
  `lib_reconstruction` (86), `pcg_recon` pass (`scratch/suites_p4/`). The fast gate passes 14 of 15 entries; the
  failure, `unit_pftc_align2D3D`'s "a class tie draws one of the tied classes", is the base's own, identical on
  the base build (Phase 0). `check_test_registry.py` and `check_descr.py` pass. No CTest entry added.
- **Thread identity.** Every new parallel path is bitwise identical at 1, 11, 16 and 24 threads, each twice, in
  the testers and in the program (`scratch/thread_identity_p3.log`, `scratch/thread_identity_p4.log`).
- **Parity with specter** (`scratch/keep/parity/`): within the lattice's interpolation error for every model and
  sign (table in the plan's Phase 2 entry), and for multislice where specter's rotation is exact.
- **Ground truth** (`scratch/keep/gt3/scores_*.json`, `scratch/keep/gt4/table_pt.md`, `atoms.txt`): the
  numbers above; full tables in the plan's Phase 3 and Phase 4 entries.

## Open for the maintainer

- **The prior near the Ewald limit.** The ML prior comes from the flat start pair's FSC, which collapses exactly
  where the flat model fails, so it is strongest there. On the nanoparticle a ten times weaker prior (`tau=10`)
  recovers most of that band (truth FSC 0.47 to 0.50 at 1.12 A against 0.18 to 0.20); rebuilding the prior from
  the wave pair under the present prior would not, since that pair's FSC is held near zero there. A prior from a
  weakly regularised first wave solve is the candidate.
- **The ribosome's last band** stays 0.02 below specter's (0.810 against 0.833) while the half maps agree: an
  error common to both halves, in the observation preparation shared with PCG (edge taper at large
  delocalisation, background-based noise normalisation), outside this run's files.
- **Noise normalisation and the scale link.** The particles' background spread was 1.035 to 1.098 of the true
  noise on these sets, because the background carries delocalised signal; the nonlinear models' absolute scale
  depends on it. The fitted gain shows it.
- **SIMPLE's wavelength** is 3.4e-4 short of the exact relativistic value: about 0.1 rad of lens phase at
  Nyquist for 2 um defocus at 300 kV. It concerns the CTF as a whole.
- **Star interchange** (Phase 0): `export_relion` writes no poses or origins, and the star import and exports
  read RELION's origins in angstroms as pixels. Not touched by this run.
- **`objfun` defaults to `euclid`** (SIMPLE's global default), so `reconstruct3D_wave` stops on a project that
  has no registered noise spectra unless `objfun=cc` is given; every run of this work passed `cc`. The user guide
  says so; whether the program should default to `cc`, as a reconstruction from imported poses usually needs, is
  the maintainer's choice (found in Phase 5, which changes no code).
- **The base's fast-gate failure** in `unit_pftc_align2D3D` (code of base commit `0147279c4`).
- **Integration** (section 4.7 of the plan): `rec_backend=wave` in `reconstruct3D` and `refine3D`, distributed
  execution, Ewald correction in the PCG and gridding backends, CTF and pose refinement.
