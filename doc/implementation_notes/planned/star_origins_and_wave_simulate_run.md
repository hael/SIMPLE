# One run: STAR origins, then wave mode in simulate_particles, speed and clean-up

Status: run plan, 2026-10-24. It carries out two plans in one unattended run,
in this order:

- `star_origins_angstrom.md`: particle origins converted between pixels and
  Angstrom on every STAR import and export path, poses in `export_relion`, a
  clear refusal of STAR files without an optics block;
- `wave_simulate_speed.md`: wave mode in `simulate_particles`, a faster wave
  operator, and the clean-up of the wave code to SIMPLE's rules.

The two plans hold the design, the decisions and the tests; this plan holds
the phases and the file table of the run. Its phase numbers replace the ones in
the two plans. Each plan's own "Progress" section records its phases.

## Phases

0. **STAR origins**: the fix and tests of `star_origins_angstrom.md` (items 1
   to 5, its section "Tests", including the external check through specter).
   Exit: as that plan's Phase 0.
1. **Wave mode in simulate_particles**: section 1 of `wave_simulate_speed.md`
   with its tests. Exit: as that plan's Phase 0.
2. **Speed**: section 2 of `wave_simulate_speed.md`, profiled first. Exit: as
   that plan's Phase 1.
3. **Clean-up, documents, reports**: section 3 of `wave_simulate_speed.md`;
   the documents of both plans (`user_guide.md` for the STAR programs, the
   `simulate_particles` guide, the wave-model policy, the test-environment rows,
   the status in `wave_model_summary.md`); both plans moved to `completed/`
   with their progress filled in; this plan moved to `completed/` with the
   run's report beside it, one report covering both plans. Exit: as the two
   plans' last phases.

The common rules of `wave_simulate_speed.md` (section "Execution") apply to
every phase.

### File table

| File | Purpose | Phases |
| --- | --- | --- |
| `doc/implementation_notes/planned/star_origins_and_wave_simulate_run.md` | this plan | each |
| `doc/implementation_notes/planned/star_origins_angstrom.md` | its progress | 0, 3 |
| `doc/implementation_notes/planned/wave_simulate_speed.md` | its progress | 1, 2, 3 |
| `doc/implementation_notes/completed/star_origins_angstrom.md` | moved when done | 3 |
| `doc/implementation_notes/completed/wave_simulate_speed.md` | moved when done | 3 |
| `doc/implementation_notes/completed/star_origins_and_wave_simulate_run.md` | this plan, moved when done | 3 |
| `doc/implementation_notes/completed/star_origins_and_wave_simulate_run_report.md` | the run's report | 3 |
| `src/main/star/simple_starproject_utils.f90` | the unit attribute of `star_flag` | 0 |
| `src/main/star/simple_starproject.f90` | import and export through the flag table | 0 |
| `src/main/star/simple_starfile.f90` | `write_ptcl2D_table`, both writers | 0 |
| `src/main/star/simple_starproject_stream.f90` | the stream's particle export | 0 |
| `src/main/star/simple_relion.f90` | `export_relion` poses and origins | 0 |
| `src/main/commanders/simple/simple_commanders_starproject.f90` | only if the conversion point needs it | 0 |
| `src/main/commanders/simple/simple_commanders_relion.f90` | only if `export_relion` needs it | 0 |
| `src/utils/math/simple_math_ctf.f90` | the scale-link function | 1 |
| `src/main/commanders/simple/simple_commanders_sim.f90` | wave mode in `simulate_particles`; the shared image writing | 1, 3 |
| `src/main/simple_simulator.f90` | the detector step, if it lives here | 1 |
| `src/main/image/simple_wave_imaging.f90` | detector step, lens set-up, real transforms, buffers, comments | 1, 2, 3 |
| `src/main/image/simple_wave_field.f90` | real-to-complex transforms, imports, comments | 2, 3 |
| `src/main/volume/simple_wave_projector.f90` | per-particle B-factor, buffers, shift phase, image writing, comments | 1, 2, 3 |
| `src/main/volume/simple_wave_multislice.f90` | buffers, real transforms, comments | 2, 3 |
| `src/main/volume/simple_pcg_lattice.f90` | only if the profile points at the gather or scatter | 2 |
| `src/main/opt/simple_nlcg_solver.f90` | imports, comments | 3 |
| `src/main/strategies/parallelization/simple_rec3D_wave_strategy.f90` | the scale-link function; comments | 1, 3 |
| `src/main/params/simple_parameters_phases.f90` | checks of wave mode, if they belong there | 1 |
| `src/main/ui/simple/simple_ui_sim.f90` | `simulate_particles` keys | 1 |
| `src/main/commanders/test/simple_commanders_test_class.f90` | sub-suite registration | 0, 1, 2 |
| `src/main/ui/simple_test/simple_test_ui_class.f90` | `suite=` help | 0, 1, 2 |
| `doc/user_guide/user_guide.md` | units and exported fields of the STAR programs | 3 |
| `doc/user_guide/simulate_particles.md` | wave mode | 3 |
| `doc/policies/3D/wave_model_policy.md` | scale-link function, simulation, speed | 3 |
| `doc/policies/test_environment_policy.md` | rows for new sub-suites | 3 |
| `doc/implementation_notes/planned/wave_model_summary.md` | status | 3 |

## Progress

(recorded in the two plans; this section lists the phases as they end)
