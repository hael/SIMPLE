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

(none yet)
