---
name: simple-refine3d
description: Use when working on SIMPLE's refine3D workflow, including CLI/UI entrypoints, the thin commanders and the strategy routines they wrap (refine3D strategy, volassemble_exec, prob_align_exec, commit_sigma2_exec), the distributed master as project owner, the per-iteration reprojection model, probabilistic pre-alignment and its dense/neighbourhood tables, 3D search and matcher strategies, the particle buffer and partial reconstruction, volume assembly, automasking, nonuniform filtering, multi-state behavior, distributed execution, and related solve3D refinement code paths.
---

# SIMPLE Refine3D

This is the single 3D refinement skill. Treat `refine3D` as a layered
fixed-point workflow: commanders choose the route, strategies control iteration,
execution mode and the reprojection model, matcher/search modules update
particle poses and write partial reconstructions, and the volume-assembly
strategy routine owns assembled-volume work.

## Start Here

1. `doc/policies/3D/refine3D_policy.md`
2. `src/main/ui/simple/simple_ui_refine3D.f90`
3. `src/main/exec/simple_exec_refine3D.f90`
4. `src/main/commanders/simple/simple_commanders_refine3D.f90`
5. `src/main/strategies/parallelization/simple_refine3D_strategy.f90`
6. `src/main/strategies/search/simple_strategy3D_matcher.f90`
7. `src/main/strategies/parallelization/simple_volassemble_strategy.f90`
8. `src/main/strategies/parallelization/simple_prob_align_strategy.f90`

For deeper navigation, read:

- `references/refine3d-code-map.md`
- `references/bayesian-model.md`
- `references/volume-postprocessing.md`

## Ownership Map

- UI/CLI contract: `ui`, `exec`, `simple_cmdline`, `src/main/params/simple_parameters.f90`
- Top-level command: `simple_commanders_refine3D.f90` (thin; `create_refine3D_strategy`)
- Iteration, execution mode, reprojection model: `simple_refine3D_strategy.f90`
- Particle-domain update/search: `refine3D_exec` in `simple_strategy3D_matcher.f90` and neighboring `search/` modules
- Probabilistic pre-alignment: `prob_align_exec` / `prob_align_neigh_exec` (`simple_prob_align_strategy.f90`), table workers `prob_tab` / `prob_tab_neigh` (`simple_commanders_prob.f90`), `search/probabilistic/simple_eul_prob_tab*.f90`
- Volume assembly: `volassemble_exec` (`simple_volassemble_strategy.f90`); PCG master: `simple_rec3D_pcg_strategy.f90`
- Sigma2 commit: `commit_sigma2_exec` (`simple_strategy_sigma2.f90`)
- Final reconstruction ending shared by all 3D workflows: `simple_final_rec.f90` (`calc_final_rec`)

## Strategy Routines And The Project Owner

- `volassemble_exec`, `prob_align_exec`, `prob_align_neigh_exec` and
  `commit_sigma2_exec` run on a caller's `(params, build)` with explicit
  arguments (thread count, `which_iter`, a `write` flag). They do not run
  `params%new`, set their own thread count (the caller restores its own) and
  draw random numbers without reseeding. Volume assembly returns the volume
  names it wrote (`vols`, `vols_even`, `vols_odd`), not `params%vols`.
- `commander_volassemble` (`simple_commanders_rec_distr.f90`),
  `commander_prob_align*` (`simple_commanders_prob.f90`) and
  `commander_commit_sigma2` (`simple_commanders_euclid.f90`) are thin
  `execute(cline)` wrappers (`run_*`) for shared-memory `refine3D`,
  `reconstruct3D` and other command-line callers. Put new logic in the
  routine, not the wrapper.
- The distributed master owns its project. It builds the reprojection model on
  its own builder (`materialize_reprojection_model_master`: general toolbox and
  segment reads plus `update_projinfo` at the first iteration only) and runs
  the pre-alignment, the fractional-update draw (`master_draws_sample3D`),
  gridding assembly (`write=.false.`) and the sigma commit on it. It writes the
  particle field once before the workers, once after the merge only when a
  polish pass follows, and once at the end of the iteration;
  `merge_algndocs(..., write=.false.)` merges in memory; `out` is written once,
  before postprocessing. The PCG master (`rec_backend=pcg`) additionally writes
  the particle field once at the end of its call (`execute_rec3D_pcg_distributed_master`
  has no `write` flag). Workers read only the file. Do not add project
  rereads, a second kept builder, or project writes from gridding assembly on
  this path.
- Shared-memory gridding assembly (`run_volassemble`) reads the project file,
  so the strategy writes the particle field before it
  (`publish_before_file_assembly`).

## Conceptual Boundary

- Particle-domain work: sampling, candidate generation, matching, pose/state updates, partial reconstructions
- Volume-domain work: assembly, even/odd handling, gridding correction, automask, nonuniform filtering, merged outputs

Preserve that split when refactoring.

## Working Rules

- Keep commanders thin, strategies iteration-focused, matcher/search modules
  search-focused, and numerical modules implementation-focused.
- Keep shared-memory and distributed refine3D behavior aligned.
- Use `parameters` as parsed workflow state. After `params%new(cline)` or a
  builder initializer, read parsed values from `params`; use `cmdline` mainly
  for pre-parse defaults and child command construction.
- Let project metadata own native sampling. Child command lines should not set
  native `box`/`smpd` or redundant `smpd_crop`; pass `box_crop` only for a
  deliberate crop.
- Keep command-line contracts registered in `src/main/params` before passing new
  keys on SIMPLE command lines.
- Use "probabilistic" precisely. SIMPLE uses probabilistic candidate generation
  and assignment with hard orientation/state updates in key paths; do not call a
  path full soft-assignment EM unless the code actually does that.
- The strategy materializes the iteration's reprojection model from the current
  Cartesian volumes (`reprojection_model_{even,odd}.bin`; prepared Cartesian
  references `cart_refvols_{even,odd}.bin` when a Cartesian pass runs). Its
  header is the authority for the matching band. Volume assembly never
  refreshes matcher references.
- Preserve the particle read contract (`doc/policies/particle_read_path_policy.md`):
  a 3D search pass reads its share once into the pass buffer
  (`fill_pass_buffer`, `simple_matcher_ptcl_io.f90`), and the reconstruction
  after the assignment barrier copies from it. With `box_crop < box` the search
  batches keep cropped observations as records (`open_record_buffer`,
  `put_record`, `seal_record_buffer`) and both reconstruction backends take
  them (`take_records`, `get_rec_observation`). Over budget the buffer stays
  empty and each phase reads the stacks once; outputs are the same.
- Preserve artifact handoffs: `assignment_part.dat`, partition-local partial
  reconstructions, the reprojection model files, state volumes, even/odd
  volumes, FSC outputs, state automasks and trailing chains are workflow
  contracts.
- For `solve3D` stage outputs, distinguish planned stage LP from saved
  low-pass diagnostic volumes. `_stageNN_lp.mrc` snapshots use measured state
  FSC resolution when available, with planned stage LP only as fallback.
- For multi-state `solve3D` with `pgrp_start != pgrp`, keep the symmetry-axis
  search and symmetrization state-local.

## Adjacent Files Worth Reading

- `doc/policies/3D/automasking_policy.md`
- `doc/policies/3D/reconstruct3D_pcg_policy.md`
- `doc/policies/3D/separate_alignment_and_reconstruction_for_multistate_peak_mem_reduction.md`
- `doc/policies/NU/nonuniform_filtering_policy.md`
- `src/main/strategies/search/simple_matcher_refvol_utils.f90` (reprojection model)
- `src/main/volume/simple_halfmap_diagnostics.f90`
- `src/main/nu_filt/*`
- `src/main/pftc/*`
- `src/main/ori/*`
- `src/main/project/*`

## Focused Companion Skills

- For solve3D `nsample`/`update_frac`, `sampled`/`updatecnt`,
  `prob_align`/`prob_tab` sampling reuse, or trailing-reconstruction weighting,
  read `.github/skills/simple-solve3d-importance-sampling/SKILL.md`.
- For fractional-update contracts, current partial reconstruction handoffs,
  the trailing accumulator chain or previous even/odd/rho compatibility, read
  `.github/skills/simple-frac-update-trailing/SKILL.md`.
- For `filt_mode=nonuniform`, `_nu_filt` reference products, spherical NU
  support, independent automask consumers, or `simple_nu_filter` changes, read
  `.github/skills/simple-nonuniform-regularization/SKILL.md`.

## Common Traps

- `refine3D` file references in docs may use logical names; verify real paths under `strategies/parallelization`, `strategies/search` (`search/probabilistic` for the tables), and `volume/`.
- Shared-memory and distributed routes should differ in launch/orchestration, not scientific ownership boundaries.
- `volassemble_exec` is the execution site for expensive volume-domain work; keep it out of the refine3D strategy, the matcher and the commander wrapper.
- `execute_from_project`, `execute_hand_on` and the builder's `spproj_src` are gone; do not reintroduce a handed-on project.
- The distributed fractional-update sample of a non-probabilistic pass is drawn once by the master over the whole project; workers reproduce it with `sample4update_reprod`.
- The master's dense table is `dists`/`inpls(nptcls, nrefs)` (each reference's distances contiguous); part-file layouts are in `references/refine3d-code-map.md`.
- `refine=cont` is the Cartesian pass (`strategy3D_cont`, `src/main/cftc`): it refines the stored poses, `inpl_cont` is forced off, and `athres_cont` (not `athres`) bounds its rotation. `pose_cont=yes` in a discrete mode makes the iteration strategy follow each discrete pass with a `refine=cont` polish over the same sample; `refine=cont pose_cont=yes` on a command line is refused. `cont_route=shift_then_joint` is a test seam, not a production option. A `refine=cont` run converges by its own motion rule (`check_cont_conv`, `pose_cont_status` slot 42).
- `refine=prob_state`, `refine3D_states pose_policy=fixed`, the `flex` key and docked `solve3D`/`solve3D_cavgs` are retired with no aliases: `flex_pca` is the only state initializer of `refine3D_states` (state-0/1 input, `nstates>=3`), and multi-state `solve3D` (`nstates>1`) refines independent states.
