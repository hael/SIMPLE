# `src/main/solve`: the stage controllers and what a stage change must keep

The workflow rules are in `doc/policies/3D/solve3D_policy.md` (section 3: the stage controller),
`doc/policies/2D/solve2D_policy.md`, `doc/policies/3D/solve3D_addon_policy.md`,
`doc/policies/3D/solve3D_chunks_policy.md` and `doc/policies/3D/solve3D_pyramid_policy.md`.

## Files

- `simple_solve3D_utils.f90`: the run state (module variables), the entry helpers
  (`prep_class_command_lines`, `set_symmetry_class_vars`, the `lpinfo` set-up) and helpers used
  outside solve3D: `gen_ortho_reprojs4viz`, `write_solve3D_lowpass_snapshot`,
  `configure_final_pcg_solve_budget`, `write_final_rec_outputs` (by `simple_final_rec`, the refine3D
  commanders, `simple_refine3D_states_controller`, SINGLE's `autorefine3D_nano`).
- `simple_solve3D_controller.f90` (a submodule of the above): the ladder (`NSTAGES`, `MAXITS`,
  `NSPACE`, the stage constants) and `set_cline_refine3D`.
- `simple_solve2D_controller.f90`: `stage_params`, the stage count, `set_cline_refine2D_stage`, the
  seed pass, the class FRC=0.5 stage omission, the mask estimate. `NSAMPLE_DEFAULT_2D` is read by
  `strategies/search/simple_matcher_smpl_and_lplims.f90`.
- `simple_solve2D_chunk.f90`: one subset of `solve2D_chunks`.
- `simple_solve3D_manifest.f90`: the run manifest; `simple_solve3D_addon_report.f90`: the add-on's
  validation against the base run; `simple_solve3D_chunks_plan.f90`: the `solve3D_chunks` split.
- `simple_state_assignment.f90`: multi-state 3D by class-average assignment, shared by stream p07
  (`simple_stream_stage_solve3D`) and `solve3D_pyramid`.

The drivers are commanders: `exec_solve3D`, `exec_solve3D_cavgs`, `exec_solve3D_addon`
(`simple_commanders_solve3D.f90`), `exec_solve2D` (`simple_commanders_solve2D.f90`).

## The contract

1. **A stage is a whole program.** Per stage the driver composes the stage's cline
   (`set_cline_refine3D(params, istage, start_stage, ...)`, `set_cline_refine2D_stage`), runs
   `refine3D`/`refine2D` through its commander in process, then reads back `cls3D`, `ptcl3D`, `out`
   (3D) or `cls2D`, `out` (2D) from the file. The callee writes `endit` on the stage cline; the next
   stage's first iteration is `endit + 1` (`init_refine3D_iteration`,
   `set_refine2D_stage_iteration_policy`).
2. **Policies compose in order.** `build_refine3D_stage_cfg` runs update, symmetry, mode, backend,
   balance, gauref, trailing reconstruction, filtering, automask, envfsc, stage controls and search
   adjustments, a later routine may overwrite an earlier one's field, and `emit_refine3D_stage_cfg`
   writes the cline. Keys the user gave at entry (`refine`, `lp`, `lpstop`, `trs`, `prob_beta_end`)
   are recorded once as the `l_refine3D_*_override` flags and `l_prob_beta_end_override`
   (`prep_class_command_lines`) and honoured
   in the mode policy, the `lpinfo` set-up and the emit, not re-read from the user's cline.
3. **Stage numbers are coupled.** `ML_REG_START_STAGE` must match the case split of
   `set_refine3D_stage_controls`; `AGGRESSIVE_MARCH_STAGE`, `VOL_SUPPLIED_STAGE`, `PCG_SOLVENT_START_STAGE`,
   `AUTOMSK_STAGE`, `ENVFSC_STAGE` are defined from other stage constants. `MAXITS` feeds
   `solve3D_remaining_niters` and the manifest's ladder.
4. **Sampling is preserved.** Every particle-to-volume handoff keeps `box * smpd == box_crop *
   smpd_crop`; each stage's crop comes from `lpinfo(istage)`. The `refine3D` and `reconstruct3D`
   stage clines carry no native `box`, `smpd` or `smpd_crop` (`prep_class_command_lines` deletes
   them; the project is the authority).
5. **Stage omission is bounded.** 3D (`solve3D_stage_behind_fsc05`): never the first or last stage,
   the axis-search stage, or any stage of a multi-state run (`solve3D_stage_skipping`). 2D
   (`solve2D_stage_behind_frc05`): never the entry stage, the first `SKIP_MIN_STAGE` (2) stages, the
   last, or a stage without a planned limit; an omitted phase-1 stage's iterations go to the next
   phase-1 stage.
6. **Run state is per entry.** `lpinfo`, the stage clines (`cline_refine3D`, `cline_symmap`,
   `cline_reconstruct3D`, `cline_reproject`), `update_frac`, `nstates_glob`, `nstages_refine3D` and
   the mode and override flags are module variables of `simple_solve3D_utils`, set by the entry
   (`exec_solve3D`, `exec_solve3D_cavgs`, `prep_class_command_lines`, `set_symmetry_class_vars`).
   `solve3D_chunks` runs `solve3D` then `solve3D_addon` in one process, and `exec_solve3D_addon`
   runs `solve3D` in process: a value an entry does not set is the previous run's. Ruled to be
   converted to a type owned by each run: add no module variable.
7. **The manifest is the add-on's only input.** Written last, atomically, never fatally; `read`
   refuses unknown fields, a bad checksum, truncation and another `MANIFEST_VERSION`, so a new field
   bumps the version. A new key that shapes the solution goes in `MANIFEST_INPUT_KEYS`, and in
   `MANIFEST_REPLAY_KEYS` when the add-on must reuse it; otherwise the add-on runs at the default.
   `balance`, `nclust`, `mskdiam` are replayed but overridable.
8. **The class-average row layout is built twice.** Even averages are rows 1..N (`eo=0`), odd ones
   N+1..2N (`eo=1`): in `simple_cavgs_work_project` (refine3D_states_cavgs, which also honours
   `state_even`/`state_odd` and mirrors the selection to `ptcl2D`) and inline in
   `exec_solve3D_cavgs`, which does neither. Known defect, ruled for fixing: one builder for both;
   until then change the layout in both. The reseed rules and `simple_state_assignment` read it.
9. **The ending.** Full schedules, every multi-state schedule and a single-state run stopped at
   the axis-search stage or later end with `calc_final_rec` (a multi-state run first assigns the
   particles no stage updated), then the manifest.

## Tests

`test=unit_project` (`solve3D_manifest`, `solve3D_addon_report`, `solve3D_chunks_plan`,
`cavgs_work_project`), `test=unit_reconstruction` (`solve3D_cavgs_stage_rules`,
`state_assignment_rules`), `test=unit_pftc_align2D3D suite=solve2D_stage_rules`,
`test=lib_reconstruction suite=addon_report_docking`. High level: `test=simulated_workflow`
(solve2D and solve3D; `suite=1jxy` also runs solve3D distributed), `test=solve3D_addon`,
`test=multistate_workflows` (the `cavgs` route is solve3D_cavgs), `test=stream_chains` and
`test=unit_stream suite=solve_3D` (state assignment through p07). No entry runs `solve3D_chunks`,
`solve3D_pyramid` or `solve2D_chunks` end to end.
