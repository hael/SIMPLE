# Report of the run star_wave: STAR origins, wave mode in simulate_particles, speed and clean-up

Run on the Dell (Intel Core Ultra 9 285, 24 cores, gfortran 15.2.1), 2026-10-24,
unattended, from the clean tree at `be6b1e510`. Nothing is committed: the change
is the working tree, and each phase's diff is in
`/home/elmlundho/agent_runs/star_wave/review/phase_N.diff` (N = 0 to 3). The
plans with their full progress sections sit beside this report:
`star_origins_angstrom.md` (Phase 0 and the STAR documents of Phase 3),
`wave_simulate_speed.md` (Phases 1 to 3) and the run plan
`star_origins_and_wave_simulate_run.md`. Paths to evidence below are under
`/home/elmlundho/agent_runs/star_wave`. No ruling was issued during the run.

## In one paragraph

SIMPLE now converts particle shifts between pixels and RELION's Angstrom
origins on every STAR path, `export_relion` writes poses, and a RELION 3.0
STAR file is refused with a clear message; an independent reader (specter)
agrees with SIMPLE's particle positions to 0.019 pixel in both directions,
where the base was 1.50 pixels off. `simulate_particles` has a wave mode that
images an ordinary SIMPLE map with the wave model, adds shot and detector
noise, and gives stacks that reconstruct to the map's scale; without the new
key the program writes the same bytes as before. One gradient pass of the
single-scatter wave models is 15 to 20% faster; multislice and the whole
`reconstruct3D_wave` are not measurably faster, for reasons found and recorded
below. The wave code follows SIMPLE's comment and import rules, and the user
guides and the wave-model policy describe what was built.

## Phase 0: STAR origins (`review/phase_0.diff`)

- One conversion per particle with the pixel size of the particle's stack:
  imports divide (in `finalise_ptcl_import`, formerly `check_stk_params`,
  after the stack's pixel size is final, guarded so it runs once), exports
  multiply (both writers of `write_ptcl2D_table`, the stream's export,
  `export_relion`). A nonzero shift without a stack pixel size stops the run
  with a message naming particle and stack.
- `export_relion` writes rotation, tilt, psi and origins from the 3D field when
  `has_valid_ptcl3D_alignment` says it holds poses, otherwise psi and origins
  from the 2D field.
- A file without `data_optics` is refused ("... has no data_optics block:
  RELION 3.1 or later STAR files are required"); the reader no longer treats an
  unnamed block as particles.
- Tests in the three STAR sub-suites of `unit_project`, with two stacks at
  different pixel sizes; deliberately broken conversions fail them (9 and 5
  failures). The external check through specter passes on `import_particles`,
  `import_starproject` (3D and 2D), `export_relion` and `export_starproject`
  (worst 0.019 pixel against 0.05). Base against run on the same project: only
  the two origin columns differ, by the pixel size 1.3.

## Phase 1: wave mode in simulate_particles (`review/phase_1.diff`)

- Keys `wave_model`, `ewald`, `total_dose`. Each state's map, masked as
  before, divided by the scale link `kappa = 2 g sigma a N`, now one function
  `electron_scale_link` used by `reconstruct3D_wave` too; one wave operator per
  state predicts the expected counts, with a lens B-factor per particle;
  Poisson shot noise, then Gaussian detector noise by `snr` as `add_gauran`
  sets it; particles in SIMPLE's units.
- Evidence: the default path byte-identical to the base in five cases (fixed
  seed, Release against Release); a wave-mode stack reconstructed by
  `reconstruct3D` (PCG) returns its map's scale shell by shell, 0.95 to 1.06
  down to 4 Angstroms (the map's phase reaches 0.3 rad, so a few percent is
  the second-order size); a new "wave simulation" sub-suite (weak limit,
  second order, per-particle B-factors, noise statistics); output bitwise
  identical at 1, 11, 16 and 24 threads; output order verified at program
  level, 66 of 66 images.
- Cost at box 128 (2000 particles, 24 threads): 2.4 to 2.7 s single-scatter,
  91 to 97 s multislice, 7 s for the Fourier-projection mode.

## Phase 2: a faster wave operator (`review/phase_2.diff`)

- Profiled first, with stage timers in scratch copies only (no profiler on the
  machine). Before: transforms a third of the thread time, gather a sixth,
  lens set-up 7 to 8%, and the gradient scatter a quarter of the wall time.
- Changes: lens set-up with precomputed angle terms and cached envelope and
  propagators (5 to 8 times faster); real-to-complex transforms for the real
  fields (transforms 25% cheaper); per-thread work buffers; separable shift
  phase; the scatter's rim passes restricted to the outer ring (bitwise
  identical, PCG included).
- Per pass, base then run, alternated: 0.77 to 0.65 ms per particle at box
  128 (projection), 5.04 to 4.28 ms at box 256; first Born and Rytov alike;
  multislice unchanged (its slice-by-slice propagation transforms dominate).
  Results move at single-precision rounding (gradients about 1e-6 relative)
  and are bitwise identical across thread counts.
- The ribosome case (4000 particles, box 320): `reconstruct3D_wave` 541.7 and
  540.5 s on the base, 535.7 and 563.8 s now: no gain beyond the spread.
  Instrumented, the operator's passes are a fifth of the run (109.6 to 100.1 s)
  and the preparation of the observations costs as much (see "Open").

## Phase 3: clean-up, documents, report (this phase's diff)

- Type headers of the six wave modules (8 to 21 lines) cut to three, file and
  procedure headers to the rules, pointing at the policy; what only those
  comments held moved into the policy (sections 3.1 and 6.1). The umbrella
  imports of `simple_wave_field` and `simple_nlcg_solver` replaced.
- Documents: a STAR section in `user_guide.md`; wave mode in
  `simulate_particles.md`; the policy (scale-link function, simulation as new
  section 11, real-field transforms, buffers, tests, costs); the test policy's
  sub-suite list; the summary's status.
- Checks after the clean-up: `check_descr.py`, `check_test_registry.py`,
  `unit_reconstruction` 456/456, `unit_numerics` 1365/1365, `unit_project`
  2914/2914 (`scratch/final_checks_p3.txt`).

## simulate_scattering: its state, for its author

`simulate_scattering` behaves as before: its output is byte-identical to the
base (`scratch/scatter_identity`). What changed around it:

- The loop that turned a predicted plane into an image moved into the operator
  as `wave_projector%plane_to_image` (it divides by the operator's gain, 1 in
  `simulate_scattering`), shared with `simulate_particles`.
- Its detector, `detect_intensity`, is now shared with the wave mode of
  `simulate_particles`, which calls it with a detector noise derived from
  `snr`; the routine itself is unchanged (including the clamp at zero after
  the Gaussian noise).
- The operator can now take a lens B-factor per particle; `simulate_scattering`
  still uses one `bfac` for all.
- The two programs now differ in what they write: `simulate_scattering`
  writes relative intensities (vacuum 1) from a potential in volts or a PDB
  model, with its own keys `total_dose` and `noise_sdev`; `simulate_particles`
  with `wave_model` writes SIMPLE's contrast-inverted particles normalised to
  unit shot noise from a SIMPLE map. Whether `simulate_scattering` should also
  write SIMPLE's particle units is the author's call.

## Decisions taken under delegation

Each is recorded, with its reason, in the plan named.

- Dead flag-table particle exports (`export_ptcls2D`, `export_ptcls3D`) and
  other uncalled STAR exports deleted rather than fixed (`star_origins_angstrom.md`).
- The pixel-size check stops only a nonzero shift; particle records always
  report a shift (same plan).
- The refusal is tested through `read_layout`, since unit tests cannot observe
  a hard stop (same plan).
- `initialise`/`kill` of `starproject` fixed, and an out-of-bounds scan in
  `export_relion` fixed (it crashed the base on a project without coordinates)
  (same plan).
- The external check adds the columns specter's reader requires (defocus;
  rotation and tilt for a 2D file) with RELION's zero default to a copy, and
  gives `export_starproject` a project with a micrograph (same plan).
- Wave mode's weak-limit test compares with twice the Fourier-projection
  particle of the deapodised map: the flat mode normalises over the
  twice-padded box (`wave_simulate_speed.md`, Phase 1).
- Second order tested by its scaling, at phases where single precision does
  not hide it (same).
- The round trip measured in a geometry without signal in the particles'
  background (box 128, 0.6 micrometres) after the first geometry showed the
  normalisation effect (same).
- Two scatter redesigns tried and removed for lack of gain; the profiling kept
  outside the repository (same, Phase 1 of that plan's numbering).
- The previous run's testers keep their long headers: not in the plan's files
  (same, Phase 2 of that plan's numbering).

## Deviations from the plan

- The STAR plan asked for the flag-table particle export to be fixed; it had
  no caller and was deleted. Its own report became this one report.
- Files added to the run plan's file table during the run: the three STAR
  testers (Phase 0), `simple_simulator_tester.f90` (Phase 1),
  `simple_wave_field_tester.f90` (Phase 2).
- The plan's weak-limit test spoke of "the flat model's particle in the same
  units"; the Fourier-projection mode is at half the wave model's scale, so the
  test applies the factor 2 and the deapodisation explicitly.
- Speed: the plan expected the four changes to speed up the operator; they did
  for the single-scatter models, not for multislice, and not visibly for the
  ribosome program as a whole.

## Evidence against the baseline

The base build (`scratch/keep/build_base`, Release of `be6b1e510`) was kept
for the whole run: STAR output compared column by column; specter's check run
on both (1.50 against 0.019 pixel); `simulate_particles` and
`simulate_scattering` byte-identical with a fixed seed; every timing alternated
base and run back to back. The fast gate fails one check,
`unit_pftc_align2D3D` ("a class tie draws one of the tied classes"), in every
phase; it fails identically on the unmodified base (debug and Release builds),
as the previous run also found.

## Open for the maintainer

- **Particle preparation in `reconstruct3D_wave`.** `read_rec_obs_batch`
  (`simple_matcher_ptcl_io`, shared with the PCG reconstruction) re-prepares
  every particle at every pass in a serial loop: 114 s of 549 on the ribosome,
  as much as the operator. Keeping the prepared planes of a half in memory
  (about 1.6 GB at box 320) or threading the loop would remove most of it.
- **The base's failing check** `unit_pftc_align2D3D`, untouched.
- **`scripts/relion2simplegui.pl`** still reads the RELION 3.0 `_rlnOriginX`
  (pixels); consistent in units, but a reader of the refused layout.
- **`import_starproject` automatic mode** looks for `run_data.star` and then
  reads `run_it<NNN>run_data.star`, a name that cannot exist (before this run).
- **`export_starproject`** writes no particle file for a project without
  micrographs (`write_ptcl2D_star`).
- **Testers' headers** of the previous run's wave testers exceed the comment
  rule.
- **`simulate_scattering`'s units** (above).
