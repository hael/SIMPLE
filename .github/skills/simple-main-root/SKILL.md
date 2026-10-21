---
name: simple-main-root
description: Use when working with the root modules directly under SIMPLE's src/main directory (builder, command-line parsing, convergence, the shared final reconstruction, refine3D stage planning and the refine3D_states controller, external-reference pose initialization, the class-average work project, simulation) or with the solve2D/solve3D controllers in src/main/solve.
---

# SIMPLE Main Root

This skill covers files directly in `src/main/`, not the subfolders, plus the
`src/main/solve/` controllers, which have no skill of their own.

## Key Files

- `simple_builder.f90`
- `simple_cmdline.f90`
- `simple_convergence.f90`
- `simple_final_rec.f90`
- `simple_refine3D_stage_plan.f90`
- `simple_refine3D_states_controller.f90`
- `simple_external_reference_pose_initialization.f90`
- `simple_cavgs_work_project.f90`
- `solve/simple_solve3D_controller.f90`, `solve/simple_solve2D_controller.f90`, `solve/simple_solve3D_utils.f90`

Probability tables moved to `strategies/search/probabilistic/`
(`.github/skills/simple-main-strategies/SKILL.md`).

## What Lives Here

- The builder: one `sp_project` plus the general, strategy2D and strategy3D
  toolboxes, built through `init_params_and_build_*` (which call
  `params%new(cline)`) and released with the matching `kill_*_tbox`;
  `new_vols` remakes the general toolbox's volumes at `box_crop` after a
  reprojection model released them
- CLI parsing and validation (`cmdline`, a large fixed-size type: never a plain
  component, see `doc/policies/compile_time_policy.md`)
- Command-line object handling before typed parameter materialization
- 2D/3D convergence reports and rules
- Shared controllers and plans for major workflows; `solve/` also holds the
  solve3D manifest, chunk planning, add-on report and multistate class-average
  assignment (shared by stream p07 and `solve3D_pyramid`)

## When Editing

- Check whether the behavior really belongs at the root or should live in a subsystem folder.
- For typed parameter fields, parsing, derived flags, or validation, use
  `.github/skills/simple-main-params/SKILL.md` instead.
- The builder reads the project from `params%projfile` only; it takes no
  in-memory project. A caller that owns a project passes its own `build` to a
  strategy routine instead of building a second one.
- Be careful with changes here because these modules often affect many executables and workflows.
