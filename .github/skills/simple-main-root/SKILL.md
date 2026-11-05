---
name: simple-main-root
description: Use when working with the root modules directly under SIMPLE's src/main (builder and toolboxes, cmdline and job descriptions, convergence, the shared final reconstruction, refine3D stage plan and startup handoff, refine3D_auto/local and refine3D_states controllers, external-reference pose initialization, class-average work project, simulation) or src/main/solve (solve2D/solve3D stage controllers, run manifest, chunk plans, add-on report, multistate class-average assignment).
---

# `src/main` root modules and `src/main/solve`: builder, command line, convergence, controllers

Unrelated modules most programs share, and the solve2D/solve3D stage controllers in `solve/`;
controllers may run whole programs through commanders (policy section 4). Not here: parameters
(`params/`), a refine iteration and the master's steps (`strategies/`), project I/O (`project/`).

- `simple_builder.f90`: `builder` (`build`), the project and the general, 2D and 3D toolboxes.
- `simple_cmdline.f90`: `cmdline`, key/value parsing, job descriptions.
- `simple_convergence.f90`: per-iteration statistics and verdicts of the refine strategies.
- `simple_final_rec.f90`: `calc_final_rec(_addon)`, the native-sampling final map of solve3D,
  refine3D_auto, refine3D_states and classify3D_refs.
- `simple_refine3D_*`: `stage_plan` (frequency blocks), `startup` (`trail_seed`, full first
  update), `workflow_controller` (refine3D_auto/local), `states_controller` (refine3D_states).
- `simple_external_reference_pose_initialization.f90`; `simple_cavgs_work_project.f90` (class
  averages as rows, even 1..N, odd N+1..2N; a second builder in `exec_solve3D_cavgs` is ruled to
  go, reference item 8); `simple_simulator.f90`.
- `solve/`: stage controllers (`simple_solve3D_controller.f90`, a submodule of
  `simple_solve3D_utils.f90`; `simple_solve2D_controller.f90`), run manifest, chunk plans, add-on
  report, `simple_state_assignment.f90` (stream p07 and solve3D_pyramid): `references/solve_controllers.md`.

## Rules

Silent ones first: nothing reports these.

- `kill_general_tbox` kills `build%spproj` and nullifies `spproj_field`, and `build_general_tbox`
  begins with it: rebuilding the general toolbox drops the project. Remake volumes with `new_vols`.
- Building writes `params`: `build_general_tbox` sets `nsym`, `eullims`; `build_strategy3D_tbox`
  sets `nrots` and adds `proj` to the field. Copy `params` fields after building.
- Order: `init_params_and_build_*` (calls `params%new`) or `build_spproj`, `build_general_tbox`,
  `build_strategy*_tbox`. The builder reads `params%projfile` whole and takes no in-memory project:
  pass the holder's `build`, never re-read what it holds.
- `obs_mask` is per process: a preparation records its mask (`note_obs_mask`); a re-preparation
  (3D partial reconstruction, the buffer's reread, class restoration) asks `obs_on_box`, or the
  observation disagrees with its record.
- `cmdline`: `set` takes any key (`params%new` ignores the unbound; a parse stops on it); an
  undefined key reads as 0 or empty. A value is a real or a string: `set` with the other kind keeps
  the old one, which `get_*` and `gen_job_descr` may read. The `INPROCESS_KEYS` (left out of job
  descriptions) are ruled to be retired: add none.
- `check_conv2D`/`check_conv3D` widen `params%trs` (unless the cline pins `trs`) and write `STATS_FILE`
  in the cwd, which `check_conv3D` reads first: once per iteration, in the run directory.
  Shared-memory refine3D rebuilds `params` and builder each iteration: what must last goes on the
  cline (as `trs`) or in the one `convergence` object (annealing-tail start, two-back history).
- `refine3D_startup%begin_iteration` sets a full update in `params` and `cline`; every path of the
  iteration must reach `finish_iteration`, or the rest of the run stays at full update.
- `calc_final_rec` reads the whole project into the caller's `spproj` (unsaved edits are lost),
  writes back `out`, and builds its child cline from nothing (`prep_final_rec_cline`): a refinement
  key the final map must honour (`m_estimator`, `automsk`, `pcg_solvent`) is passed there or lost.
- `solve/`: run state is module-level in `simple_solve3D_utils` (`lpinfo`, `cline_refine3D`,
  `update_frac`, override flags) and `solve3D_chunks` runs solve3D then solve3D_addon in one
  process: set new state at every entry. Ruled to be converted to a type owned by each run: add no
  module state. A stage's cline is composed by `set_cline_refine3D` /
  `set_cline_refine2D_stage` and numbered from `endit`. A new key that shapes the solution joins
  `MANIFEST_INPUT_KEYS` (and `_REPLAY_KEYS`), or solve3D_addon runs it at its default.

## Tests

- Seconds each: `test=unit_core suite=command_line`; `test=unit_cart_align3D` `suite=pose_statistics`
  (convergence), `suite=startup_state`; `unit_project suite=cavgs_work_project`;
  `unit_reconstruction suite=wave_simulation` (simulator).
- `solve/`: `test=unit_project` (`solve3D_manifest`, `solve3D_addon_report`, `solve3D_chunks_plan`),
  `test=unit_reconstruction` (`solve3D_cavgs_stage_rules`, `state_assignment_rules`),
  `test=unit_pftc_align2D3D suite=solve2D_stage_rules`, `lib_reconstruction suite=addon_report_docking`.
- Builder, cmdline, convergence: `test=units`, then `test=simulated_workflow suite=1jxy`; a
  controller: its entry (skill `simple-main-commanders`); refactors: `test=identity`.

## See also

`doc/policies/layering_and_project_ownership_policy.md`, `doc/policies/3D/solve3D_policy.md`,
`doc/policies/2D/solve2D_policy.md`; skills `simple-main-params`, `simple-main-commanders`,
`simple-main-strategies-parallelization`, `simple-refine3d`.
