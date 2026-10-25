# Wave mode in simulate_particles, a faster wave operator, and clean-up

Status: plan, 2026-10-24, decisions of the same day applied. Nothing
implemented. It follows the `wave_operators` run
(`ghostbuster_operators.md` and its report) and section 3 of
`wave_model_summary.md`, which explains the terms used here (wave model, flat
model, scale link).

It starts from a tree that holds the `wave_operators` changes.

## Decisions (Hans, 2026-10-24)

1. `simulate_particles` gets a wave option; `simulate_scattering` stays as it
   is (the route from a PDB model).
2. Wave mode takes an ordinary SIMPLE map and converts it to a potential with
   the scale link, from the dose, the voltage, the pixel size and the box.
3. Noise in wave mode: Poisson shot noise at the dose on the expected electron
   counts, then Gaussian detector noise in real space as `simulate_particles`
   adds it today.

## 1. simulate_particles in wave mode

- Keys `wave_model` (`projection|first_born|rytov|multislice`) and `ewald`
  (`neg|pos|no`), the keys `reconstruct3D_wave` uses, with the same defaults
  and refusals. Without `wave_model` the program does exactly what it does
  today, byte for byte.
- Wave mode needs `total_dose` (it sets both the scale and the shot noise) and
  `ctf=yes` (a phase object in focus has almost no contrast).
- **Input**: each state's map, as SIMPLE writes maps, divided by the scale link
  `kappa = 2 g sigma a N` (`g` the square root of `total_dose` times the pixel
  area, `sigma` the interaction constant at `kv`, `a` the pixel size, `N` the
  box). The scale link gets one home, a function next to the electron
  wavelength and interaction constant in `simple_math_ctf`, used by
  `reconstruct3D_wave` (which computes it inline today) and by
  `simulate_particles`.
- **Image**: the wave operator's prediction of each particle at its pose,
  shift and CTF (`wave_projector`), as expected electron counts per pixel,
  `E = total_dose a^2 I` with `I` the relative intensity.
- **Noise**: shot noise first, each pixel's count drawn from a Poisson
  distribution of mean `E`; then detector noise, Gaussian in real space, its
  standard deviation `sqrt(var(E) / snr)` from the noiseless image's variance,
  as `add_gauran` sets it today, and none for `snr` of 5 or more.
- **Output**: particles in SIMPLE's convention, `x = -(counts - D a^2) / sqrt(D a^2)`
  (contrast inverted, normalised as the scale link assumes), with the
  orientation table as today.
- Per-particle B-factors (`bfac`, `bfacerr`) act on the lens envelope of each
  particle, so the operator takes a B-factor per particle. Multiple states,
  `recthres`, `even`, parallel parts (`part`, `fromp`, `top`) work as today:
  only the per-particle image changes. Multislice is far slower than the
  single-scatter models; the guide says so, with the measured cost.
- The detector step (shot noise, then detector noise) is one routine shared
  with `simulate_scattering`'s detector, which keeps its own keys.

Tests: the weak-signal limit through the program's path (a weak map gives the
flat model's particle, up to second order, in the same units); the noise
statistics in an empty region (mean `D a^2`, variance `D a^2 + sdev^2` before
the normalisation); the scale link round trip (a wave-mode stack of a map,
reconstructed by `reconstruct3D` with PCG, returns the map's scale at low and
middle resolution within the second-order terms; with shot noise only, `snr` of
5 or more, since detector noise raises the background the scale link assumes to
be shot noise); the default path unchanged
(`simulate_particles` without `wave_model` byte-identical to the base).

## 2. A faster wave operator

Measured first: per-stage timers for one gradient pass (gather, lens set-up,
FFTs, image formation, scatter) at box 128 and 256, for Rytov, first Born,
projection and multislice, on the Release build, nothing else running. Then:

- **Lens set-up.** `wave_imaging%set_optics` rebuilds the lens on the whole
  fine grid for every particle on every pass, with several transcendental
  calls per point through a type-bound call. The parts that do not change
  between particles (the angle terms of the astigmatism, the B-factor
  envelope at a common B, the Rytov propagator and its inverse, which depend
  only on the wavelength) are computed once; per particle, one phase
  evaluation inline and one sine and cosine per point remain.
- **Real fields.** The intensity, its derivative and the real part taken in
  the adjoint are real fields transformed with complex FFTs; real-to-complex
  and complex-to-real transforms take half the work.
- **Work arrays.** Arrays allocated per particle inside the parallel loops
  (six planes per particle, two cubes per particle for multislice) become
  per-thread buffers allocated once, sized at the point of use.
- **Shift phase.** Computed per pixel with a sine and cosine; separable phase
  vectors along `h` and `k` give the same factor with a product.

Each change keeps the operator's tests passing at their current tolerances
(adjoints, finite differences, weak-signal identity, specter's golden images)
and the results identical at 1, 11, 16 and 24 threads; outputs may move at
float rounding only. The gain is reported per stage and on the ribosome case
of the previous run (`reconstruct3D_wave` against `reconstruct3D` with PCG,
24 threads, twice each).

## 3. Clean-up to SIMPLE's rules

- Comments to the `simple-modern-fortran` rules: file headers of at most five
  lines after `!@descr:`, procedure headers of at most three, no restating of
  the policy; the type headers of the new modules (up to 20 lines today)
  shortened, pointing at `doc/policies/3D/wave_model_policy.md` instead of
  the plan.
- No umbrella imports in the new modules (`simple_wave_field`,
  `simple_nlcg_solver` import `simple_core_module_api` today).
- The loop in `simulate_scattering`'s commander that writes predicted planes
  into an image moves into the operator, shared with `simulate_particles`.

## Execution

Common rules: build in the run's build tree, the fast gate and
`check_test_registry.py` at every phase end; new or changed OpenMP code tested
at 1, 11, 16 and 24 threads, twice each; no tolerance relaxed; timings on the
Release build, nothing else heavy running, twice each way; bugs found are fixed
in the phase that finds them; prose readable on its own.

0. **simulate_particles in wave mode** (section 1) with its tests. Exit: the
   tests pass; the default path is byte-identical to the base; a stack made in
   wave mode reconstructs with `reconstruct3D` to its map's scale.
1. **Speed** (section 2). Exit: the per-stage profile before and after
   recorded; the operator's tests pass unchanged; thread identity; the
   ribosome timing recorded.
2. **Clean-up, documents, report** (section 3; the user guide of
   `simulate_particles`, the wave-model policy, the test-environment rows, the
   summary's status). Exit: `check_descr.py` and `check_test_registry.py`
   pass; this plan moved to `completed/` with the report beside it.

### File table

| File | Purpose | Phases |
| --- | --- | --- |
| `doc/implementation_notes/planned/wave_simulate_speed.md` | this plan: progress | each |
| `src/utils/math/simple_math_ctf.f90` | the scale-link function | 0 |
| `src/main/commanders/simple/simple_commanders_sim.f90` | wave mode in `simulate_particles`; the shared image writing | 0, 2 |
| `src/main/simple_simulator.f90` | the detector step, if it lives here | 0 |
| `src/main/image/simple_wave_imaging.f90` | detector step, lens set-up, real transforms, buffers, comments | 0, 1, 2 |
| `src/main/image/simple_wave_field.f90` | real-to-complex transforms, imports, comments | 1, 2 |
| `src/main/volume/simple_wave_projector.f90` | per-particle B-factor, buffers, shift phase, image writing, comments | 0, 1, 2 |
| `src/main/volume/simple_wave_multislice.f90` | buffers, real transforms, comments | 1, 2 |
| `src/main/volume/simple_pcg_lattice.f90` | only if the profile points at the gather or scatter | 1 |
| `src/main/opt/simple_nlcg_solver.f90` | imports, comments | 2 |
| `src/main/strategies/parallelization/simple_rec3D_wave_strategy.f90` | the scale-link function; comments | 0, 2 |
| `src/main/params/simple_parameters_phases.f90` | checks of wave mode, if they belong there | 0 |
| `src/main/ui/simple/simple_ui_sim.f90` | `simulate_particles` keys | 0 |
| `src/main/commanders/test/simple_commanders_test_class.f90` | sub-suite registration | 0, 1 |
| `src/main/ui/simple_test/simple_test_ui_class.f90` | `suite=` help | 0, 1 |
| `doc/user_guide/simulate_particles.md` | wave mode | 2 |
| `doc/policies/3D/wave_model_policy.md` | scale-link function, simulation, speed | 2 |
| `doc/policies/test_environment_policy.md` | rows for new sub-suites | 2 |
| `doc/implementation_notes/planned/wave_model_summary.md` | status of its section 3 | 2 |
| `doc/implementation_notes/completed/wave_simulate_speed.md` | this plan, moved when done | 2 |
| `doc/implementation_notes/completed/wave_simulate_speed_report.md` | the run's report | 2 |

## Progress

### Phase 0 (run star_wave, Phase 1), 2026-10-24: wave mode in simulate_particles done

**In plain words.** `simulate_particles` can now draw particles with the wave
model: `wave_model=projection|first_born|rytov|multislice` (and `ewald`) take
an ordinary SIMPLE map, convert it to a potential in volts with the scale link
of `total_dose`, predict each particle's electron counts with the wave
operator, add Poisson shot noise at the dose and then Gaussian detector noise
set by `snr`, and write particles in SIMPLE's units. Reconstructed with
`reconstruct3D`, such a stack returns its map's scale to within a few percent
at low and middle resolution. Without `wave_model` the program writes exactly
the same files as before.

**Files changed.** `src/utils/math/simple_math_ctf.f90` (the scale-link
function `electron_scale_link`), `src/main/strategies/parallelization/simple_rec3D_wave_strategy.f90`
(uses it instead of its inline formula), `src/main/volume/simple_wave_projector.f90`
(lens B-factor per particle; `plane_to_image`, the loop that turns a predicted
plane into an image, moved here from `simulate_scattering`'s commander),
`src/main/simple_simulator.f90` (`simimg_wave`, the detector of wave mode),
`src/main/commanders/simple/simple_commanders_sim.f90` (the wave branch of
`simulate_particles`; `simulate_scattering` calls `plane_to_image`),
`src/main/ui/simple/simple_ui_sim.f90` (the keys `wave_model`, `ewald`,
`total_dose` and the help text), the new tester `src/main/simple_simulator_tester.f90`
with its registration as the "wave simulation" sub-suite of
`unit_reconstruction` (`simple_commanders_test_class.f90`,
`simple_test_ui_class.f90`), this plan and the run plan (file table).

**How it works.** Wave mode is chosen by the presence of `wave_model`; it
needs `total_dose` and `ctf=yes`, takes `ewald`'s default and refusals from
`reconstruct3D_wave`, and refuses `ewald` or `total_dose` without
`wave_model`. Each state's map is masked as in the flat mode and divided by
the scale link `kappa = 2 g sigma a N` (`g = sqrt(total_dose) a`, `sigma` the
interaction constant at `kv`, `a` the pixel size, `N` the box); one wave
operator per state predicts its particles in batches, in the output order of
the flat mode (the particles of a state that fall in a batch are consecutive
in that state's selection). Each predicted plane becomes the relative
intensity `I` of the native box; `simimg_wave` draws the counts from a
Poisson distribution of mean `D a^2 I` (`D` the dose), adds Gaussian noise of
standard deviation `sqrt(var(E) / snr)` for `snr` below 5 (the rule of
`add_gauran`, `E` the noiseless counts), and writes
`x = -(counts - D a^2) / sqrt(D a^2)`. The detector step is
`detect_intensity`, shared with `simulate_scattering`. `bfac` and `bfacerr`
are drawn per particle as in the flat mode and act on each particle's lens
envelope. Multiple states, `recthres`, `even` and parallel parts work as in
the flat mode.

**Decisions taken during the phase.**

- *The flat mode's particles are at half the scale the wave model refers to.*
  The flat chain projects the volume padded to twice the box and normalises
  over the padded box, without Kaiser-Bessel deapodisation: its particle is
  the sum of the map along the beam divided by `2N` (measured: the image sum
  is 1/64 x 0.959 of the map sum at box 32). The wave model's linear part is
  the PCG (preconditioned conjugate gradients) reconstruction's forward model,
  the sum divided by `N`, which the scale link and `reconstruct3D` use. The
  flat mode is not changed (its files stay byte-identical, and its absolute
  scale never matters: its noise is set relative to the image and
  reconstruction normalises particles to unit noise). The weak-signal test
  therefore compares the wave particle's linear part with twice the flat
  particle of the deapodised map.
- *Second order is checked by how it scales.* The even part's share of the
  odd part doubles when the map doubles; at phases below 0.01 radian the
  single-precision rounding of the unit intensity (about 1e-5) hides the
  second-order term, so the check runs at phases 0.015 and 0.030. A first
  bound, "the even part is at most twice the phase of the odd part", was
  wrong: the CTF suppresses the linear part at low frequency but not the
  second-order part (measured share: 0.13 to 0.14 times the phase).
- *The round trip needs a background without signal.* `reconstruct3D`
  normalises each particle by its background outside the mask. With the
  myoglobin map in a box of 64 at 1.5 micrometres of defocus the CTF spreads
  signal into that background (standard deviation 1.108 instead of 1), and
  the reconstruction came back at 0.84 of the map's scale at every resolution.
  With the map in a box of 128 at 0.6 micrometres (background 1.003) it comes
  back at the map's scale (below). This is the same effect the policy notes
  for the fitted gain of real data.
- The program-level bookkeeping is checked by hand rather than in a unit test:
  a commander ends the process, which a unit test may not.

**Tests and evidence** (paths under `/home/elmlundho/agent_runs/star_wave`).

- *Weak-signal limit through the program's chain* ("wave simulation":
  `test_weak_limit`, box 64 at 1 Angstrom, 0.45/0.35 micrometres astigmatic
  defocus, three smooth Gaussian blobs): the part of the wave particle that is
  odd in the map's sign equals twice the flat particle of the deapodised map
  to 5.7e-3 (projection) and 5.9e-3 (first Born), within the lattice's
  interpolation tolerance 1.1e-2 used by the operator's tests, at a largest
  phase of 3.7e-4 radian; uniform over all frequency shells (scale 2.00 to
  2.02, correlation 0.99998). The even part doubles with the map: ratio 1.975
  (projection) and 2.013 (first Born), tolerance 2 x 0.030.
- *Per-particle B-factors* (`test_bfac_per_particle`): a particle with its own
  B-factor predicts bitwise as an operator with that B-factor for all.
- *Noise in an empty region* (`test_detector_noise`, seeded): counts of mean
  39.98 (expected 40.00) and variance 59.66 (expected 59.85 = `D a^2` plus the
  detector variance), within four standard errors.
- `unit_reconstruction`: all 452 checks pass on the debug build
  (`scratch/t_unit1/unit_rec.log`); the library suite "wave reconstruction"
  passes on the Release build (`scratch/lib_wave.log`), so `reconstruct3D_wave`
  with the shared scale-link function still meets its checks.
- *Scale-link round trip* (`scratch/scale`, scripts `run_scale2.sh` and
  `compare_scale.py`): the previous run's myoglobin potential (from its
  identity case, box 64 at 1 Angstrom, padded to box 128), 3000 wave-mode
  particles at 40 electrons per square Angstrom, 0.6 micrometres defocus,
  shot noise only (`snr=10`), reconstructed with `reconstruct3D
  rec_backend=pcg`; the mean of the two unfiltered half maps against the input
  map, shell by shell: projection 1.00 to 1.06 from 128 to 32 Angstroms and
  0.95 to 1.00 from 26 to 4 Angstroms; Rytov 1.01 to 1.05 and 0.96 to 1.00.
  The map is not weak (largest phase 0.29 to 0.37 radian, root mean square
  0.08), so deviations of a few percent are the size of the second-order
  terms, positive at low resolution as for a phase object.
- *Default path unchanged* (`scratch/default_identity/run_cases.sh`, seed
  fixed with `SIMPLE_SEED`): five cases (noise, B-factors, two states, no CTF
  with even projections, `recthres`) give byte-identical stacks and
  orientation files on the base and the run's Release builds.
  `simulate_scattering` is byte-identical too (`scratch/scatter_identity`).
- *Bookkeeping at program level* (`scratch/bookkeeping`): with the same seed,
  a flat run (no noise) and a wave run (two different maps, `recthres`,
  B-factors, one part of the particle range) write identical orientation
  files and 66 images each; for every wave image the best-matching flat image
  (correlation at low resolution) is the one with the same index, 66 of 66.
  A review found and this phase fixed a size mismatch of the B-factor
  selection when `recthres` drops particles.
- *Threads* (`scratch/threads/threads.log`): wave-mode stacks (Rytov and
  multislice, two states, B-factors, detector noise) are bitwise identical at
  1, 11, 16 and 24 threads, twice each.
- *Cost* (`scratch/cost`, Release, 24 threads, nothing else running, 2000
  particles at box 128): flat 6.9 to 7.4 s, projection 2.7 s, first Born 2.5
  to 2.7 s, Rytov 2.4 s, multislice 91 to 97 s, about 1.2 to 1.4 ms per
  particle for the single-scatter models and 47 ms for multislice, 35 to 40
  times more.
- `check_test_registry.py` and `check_descr.py` pass. The fast gate
  (`scratch/fast_gate_p1.log`): 14 of 15 suites pass; the one failure,
  `unit_pftc_align2D3D` ("a class tie draws one of the tied classes"), is the
  base's own (Phase 0 of this run). `unit_reconstruction` takes 21.0 s in the
  gate (16.9 s before), total 21.0 s of real time against a budget of 60.

**Exit items.** The tests pass; the default path is byte-identical to the
base; a stack made in wave mode reconstructs with `reconstruct3D` to its map's
scale (above).

**Not covered.** The weak-signal test runs the models without curvature
(`ewald=no`), the case in which the flat model is their limit; the curved
models' limit is the flat model on the Ewald paraboloid, which the operator's
own tests cover.

### Phase 1 (run star_wave, Phase 2), 2026-10-24: a faster wave operator done

**In plain words.** One gradient pass of the single-scatter wave models
(projection, first Born, Rytov) now takes 15 to 20% less time on this machine.
The lens is set up about six times faster, the Fourier transforms of the real
fields cost a quarter less, and nothing else changed measurably. Multislice is
as fast as before: its time goes into the transforms of its slice-by-slice
propagation, which none of the planned changes touch. On the ribosome case of
the previous run the whole `reconstruct3D_wave` is about as fast as before
(within the run-to-run spread), because the operator's passes are only a fifth
of that program's time; the largest single cost there, the preparation of the
particle images at every pass, is outside this plan's files and is left as a
finding below.

**Files changed.** `src/main/image/simple_wave_imaging.f90` (lens set-up;
real-field transforms of the intensity, its derivative and the adjoint's real
part), `src/main/image/simple_wave_field.f90` (an optional real field with
real-to-complex and complex-to-real transforms), `src/main/volume/simple_wave_projector.f90`
(per-thread work buffers, separable shift phase), `src/main/volume/simple_pcg_lattice.f90`
(the rim passes of the scatter visit only the outer ring), the tester
`src/main/image/simple_wave_field_tester.f90` (direct checks of the real-field
transforms; added to the run plan's file table), this plan and the run plan.

**How it was measured.** No profiler is installed (and unprivileged
performance counters are off), so the profile comes from stage timers added
to scratch copies of the source only: a timer module and a script that
inserts the calls (`scratch/prof/simple_wave_prof.f90`, `instrument.py`), and
a driver that runs one gradient pass of the operator, 24 threads, at box 128
(1024 particles; multislice 96) and box 256 (512; multislice 48)
(`scratch/prof/profile_pass.f90`). Nothing timing-only is in the repository.
The stages: gather (the volume's section or cube), lens set-up, the 2D Fourier
transforms, the forward and adjoint image formation (which contain those
transforms), and the scatter of the gradient, which runs between particle
batches while the other threads wait. Other users' load changed during the
work (load averages 1 to 19), so the final numbers come from base and run
alternated back to back, twice each (`scratch/prof/final_profile.sh`, log
`scratch/prof/profile_final.log`).

**Where the time went, before.** For the single-scatter models at box 128,
per pass (thread-seconds, of about 19 in all): Fourier transforms 6.4, gather
3.1, lens set-up 1.4 to 1.6, the rest of the forward and adjoint chains about
1.3, and the scatter 0.20 s of wall time (a quarter of the pass) while 23 of 24
threads wait. Box 256 was similar. Multislice: its transforms and adjoint sweep
dominate, the scatter 15 to 19% of the wall time.

**What was changed, and what each change bought** (box 128, projection, unless
said otherwise; thread-seconds per pass, base then run):

- *Lens set-up* (`wave_imaging%set_optics`): the angle terms of the
  astigmatism are precomputed on the fine grid, the B-factor envelope is kept
  until the B-factor or pixel size changes, the Rytov propagator and its
  inverse until the wavelength or depth changes; per particle one phase
  (`ctf%eval_phase` with the angle terms expanded) and one sine and cosine per
  point remain. 1.41 to 0.22; box 256 3.50 to 0.42. The inverse propagator
  stays a separate factor: folded into the phase it no longer cancelled the
  propagator to single precision, and the test that first Born is Rytov's
  linearisation caught it.
- *Real fields* (`wave_field`, `wave_imaging`): the intensity, its derivative
  and the real part taken in the adjoint use real-to-complex and
  complex-to-real transforms out of place, the native disk read through the
  Hermitian symmetry; the transform of the adjoint's real part adds each disk
  point and its conjugate mirror into the half spectrum, which also covers the
  grid without oversampling, where the disk's two Nyquist points share a bin
  (a first version missed that, and the operator's transpose test without
  oversampling caught it). Transforms 6.39 to 4.76; box 256 17.7 to 13.3.
- *Work arrays* (`wave_projector`): the six planes per particle (and the cubes
  of multislice) allocated inside the parallel loops became per-thread buffers
  allocated once, sized by the team before each parallel region; the cubes are
  sized only when a pass needs them. *Shift phase*: separable vectors along h
  and k replace a sine and cosine per point. Together with the above, the
  forward chain 4.17 to 3.43 and the adjoint chain 3.56 to 3.03.
- *Scatter* (`pcg_lattice`, which the plan allows when the profile points at
  the gather or scatter): its serial pass over the wrapping rim scanned every
  plane point, and now visits only each row's outer segments, the same points
  in the same order (bitwise identical results, for the PCG reconstruction as
  well); 0.20 to 0.18 s. Two larger designs were tried and removed: all planes
  of a batch in one parallel region, and an owner-computes scatter (each
  thread adds, in the original order, only the parts of the windows in the
  lattice blocks it owns, so no barriers and identical results at any team).
  Neither was faster: one thread needs 2.5 ms per plane at box 128 (350 000
  scattered additions into a 135 MB accumulator), so the scatter is limited by
  memory, and both designs reach about 13 times on 24 threads, as the existing
  one does.

**Per pass, base then run** (two runs each, alternated):

| Model, box (particles) | Base (s) | Run (s) | Per particle (ms) |
| --- | --- | --- | --- |
| projection, 128 (1024) | 0.790, 0.793 | 0.660, 0.669 | 0.77 to 0.65 |
| first Born, 128 (1024) | 0.776, 0.777 | 0.644, 0.649 | 0.76 to 0.63 |
| Rytov, 128 (1024) | 0.835, 0.807 | 0.683, 0.680 | 0.80 to 0.67 |
| multislice, 128 (96) | 10.87, 11.92 | 11.60, 12.21 | 119 to 124 |
| projection, 256 (512) | 2.581, 2.582 | 2.187, 2.192 | 5.04 to 4.28 |
| first Born, 256 (512) | 2.506, 2.338 | 2.091, 1.798 | 4.73 to 3.80 |
| Rytov, 256 (512) | 2.234, 2.211 | 1.861, 1.858 | 4.34 to 3.63 |
| multislice, 256 (48) | 49.68, 50.27 | 49.89, 49.76 | 1041 to 1038 |

Gather and scatter are unchanged in every row; multislice's differences are
within the spread of its two runs.

**The ribosome case** (the previous run's set: 4000 particles, box 320,
1.2 Angstrom, 100 kV; `scratch/ribo`, project from its recipe,
`timing.sh` with 24 threads, base and run alternated, twice each):
`reconstruct3D_wave` (Rytov, `ewald=neg`, 4 iterations) 541.7 and 540.5 s on
the base, 535.7 and 563.8 s on the run; `reconstruct3D` with PCG 88.1 and 89.9
s against 86.0 and 81.7 s. No gain beyond the spread. Instrumented once each
(`prof_run.sh`, logs `scratch/ribo/gt3/prof_base.log` and `prof_run.log`), the
549 s of the base split as: preparing the observations 113.8 s, the gradient
passes 109.6 s, the flat PCG start 63.7 s, line search and curvature 24 s,
volume transforms 22 s, and about 190 s in the remaining steps; the run's
gradient passes took 100.1 s (9% less), its lens set-up 71 instead of 305
thread-seconds, its transforms 544 instead of 656.

**Finding for the maintainer.** At every objective evaluation
`read_rec_obs_batch` (`simple_matcher_ptcl_io`, shared with the PCG
reconstruction, not in this plan's files) prepares each particle again in a
serial loop (normalisation, crop, taper, a Fourier transform), although the
observations do not change between evaluations: 114 s of the ribosome's 549.
Keeping the prepared planes of a half in memory (about 1.6 GB at box 320) or
threading the loop would remove most of it; either is a change to the shared
observation frame, or a memory trade-off, for the maintainer to decide.

**Tests and evidence** (paths under `/home/elmlundho/agent_runs/star_wave`).

- The operator's tests pass unchanged, at their tolerances (adjoints, finite
  differences, the weak-phase identity, specter's golden images):
  `unit_reconstruction` 456 of 456 checks on the debug build
  (`scratch/t_unit2/unit_rec_p2.log`), including the new direct checks that
  the real-field transforms equal their complex counterparts to 1e-5 on the
  oversampled and the native grid; the library suite "wave reconstruction"
  passes on the Release build (`scratch/t_lib2/lib_wave_p2.log`).
- The results move at single-precision rounding only: against the base, the
  gradient of one pass differs by 1.4e-6 (box 128) and 3.6e-6 (box 256)
  relative, the loss by 1e-8 to 3e-8; the scatter change alone was bitwise
  identical.
- Threads: one gradient pass of every model is bitwise identical at 1, 11, 16
  and 24 threads, twice each (`scratch/prof/thread_identity.log`), and the
  operator's own thread-identity test passes.
- `check_test_registry.py` and `check_descr.py` pass; the fast gate
  (`scratch/fast_gate_p2.log`): 14 of 15, the one failure
  `unit_pftc_align2D3D` the base's own. No file mode changes.

**Exit items.** The per-stage profile before and after is recorded above; the
operator's tests pass unchanged; the results are identical across thread
counts; the ribosome timing is recorded.

### Phase 2 (run star_wave, Phase 3), 2026-10-24: clean-up, documents and report done

**In plain words.** The wave code now follows SIMPLE's comment and import
rules, and the documents describe what the run built.

**Files changed.** `src/main/image/simple_wave_field.f90`,
`src/main/image/simple_wave_imaging.f90`, `src/main/volume/simple_wave_multislice.f90`,
`src/main/volume/simple_wave_projector.f90`, `src/main/opt/simple_nlcg_solver.f90`,
`src/main/strategies/parallelization/simple_rec3D_wave_strategy.f90` (comments
and imports only); `doc/user_guide/simulate_particles.md`,
`doc/user_guide/user_guide.md`, `doc/policies/3D/wave_model_policy.md`,
`doc/policies/test_environment_policy.md`,
`doc/implementation_notes/planned/wave_model_summary.md`; the three plans,
moved to `doc/implementation_notes/completed/`, and the run's report beside
them.

**What was done.**

- *Comments.* The type headers of the six wave modules (8 to 21 lines) are
  three lines each, pointing at `doc/policies/3D/wave_model_policy.md`; the
  file headers point at the policy instead of the plan; three procedure
  headers of 4 to 10 lines are three lines. What only those comments held
  moved into the policy: the field's transform conventions (new section 3.1)
  and the solver's client contract and algorithm (new section 6.1).
- *Imports.* `simple_wave_field` and `simple_nlcg_solver` import from
  `simple_defs` and `simple_error` instead of the umbrella
  `simple_core_module_api`.
- *Image writing.* The loop of `simulate_scattering`'s commander that turned a
  predicted plane into an image moved into the operator
  (`wave_projector%plane_to_image`) in Phase 1, when `simulate_particles`
  needed it; `simulate_scattering`'s output stayed byte-identical.
- *Documents.* The `simulate_particles` guide has an example, the three keys,
  a "Wave mode" section (scale link, detector, units, the refusals, the two
  modes' scales, the measured cost); the policy names the scale-link
  function, describes simulation (new section 11), the real-field transforms,
  the per-thread buffers, the new tests and the measured costs; the test
  policy lists the "wave simulation" sub-suite; the summary marks its sections
  2 and 3 done.

**Not changed.** The testers of the previous run keep their long headers (up
to 30 lines); they are not in this plan's files. The older long comments of
`simulate_nanoparticle` in `simple_commanders_sim.f90` are outside the wave
code.

**Evidence** (`scratch/final_checks_p3.txt`, paths under
`/home/elmlundho/agent_runs/star_wave`). After the clean-up, on the debug
build: `check_descr.py` and `check_test_registry.py` pass;
`unit_reconstruction` 456 of 456 checks, `unit_numerics` 1365 of 1365,
`unit_project` 2914 of 2914; the fast gate 14 of 15, the one failure
`unit_pftc_align2D3D` the base's own (it fails identically on the unmodified
base, see Phase 0 of the run).

**Exit items.** `check_descr.py` and `check_test_registry.py` pass; this plan
moved to `completed/` with the report beside it.
