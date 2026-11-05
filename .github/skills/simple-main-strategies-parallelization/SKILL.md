---
name: simple-main-strategies-parallelization
description: Use when working in SIMPLE's src/main/strategies/parallelization: the create_*_strategy shared-memory, distributed-master and worker strategies (refine3D, refine2D, reconstruct3D and PCG, make_cavgs, cls_expansion, flex_pca, calc_pspec, preprocess, motion_correct, ctf_estimate, pick, extract, reextract, refine_motion_model), the distributed master's iteration, sample draw and project writes, and the master steps on its own (params, build): volassemble_exec, prob_align_exec/_neigh/2D, cavgassemble_exec, commit_sigma2_exec.
---

# `strategies/parallelization`: execution strategies and the master's steps

Owns how a program runs: `create_<x>_strategy(cline)` picks shared memory or the distributed master
(refine2D/refine3D: `nparts` given and `part` not), the iteration, the scheduling of workers
(`qenv%gen_scripts_and_schedule_jobs`) and the master's steps on its own `(params, build)`:
`volassemble_exec`, `prob_align_exec`/`prob_align_neigh_exec`/`prob_align2D_exec`, `cavgassemble_exec`,
`commit_sigma2_exec`, the PCG master. A worker's pass is the matcher (`search/`); the class sums are
`class/`, the sigma2 store `sigma2/`, the volumes `volume/`.

## Rules

Silent ones first.

- Draw, write, then schedule. The master draws the iteration's sample once over the whole project
  (`prob_align*_exec`, or `sample_ptcls4iter3D`/`sample_ptcls4update2D` under `master_draws_sample3D/2D`)
  and writes the particle field before the workers start. A worker's `sample4update_reprod` takes the latest
  sampling in the file: an unwritten draw makes it reproduce the previous iteration's sample without error.
  Master and workers decide with the same `master_draws_sample*` on the same command line: change it there.
- refine3D master (`distr_execute_iteration`): model (`materialize_reprojection_model_master`) -> sample and
  field write -> remove partials, `prepare_canonical_sigma_update` -> workers ->
  `merge_algndocs(write=.false.)` (`.true.` when a polish pass follows; the polish reproduces that sample) ->
  `volassemble_exec(write=.false.)` or the PCG master (writes the field) -> `out` written before
  postprocessing -> convergence -> `commit_sigma2_exec` -> one field write. Under `sigma_commit_deferred=yes`
  (an `internal:` key, ruled to be retired) the last iteration leaves its ranges uncommitted for the caller
  that set it (`commit_deferred_sigma_update` in solve3D). refine2D master: stale part files deleted,
  carry-over check (`carryover_needs_full_update`), sample, sigma2 prepare, workers, `merge_algndocs`,
  `run_distributed_cavg_assembly`, commit, convergence (it reads the `res` that `cavger_gen2Dclassdoc`
  writes), the field in `distr_finalize_iteration`.
- The blend, the sigma2 commit and convergence read the master's particle field in memory: merge first.
- Shared-memory refine3D rebuilds `params` and the builder from the file every iteration: a field change
  not written by `inmem_finalize_iteration` is lost. The commanders it runs read the file: write the field
  before them (`publish_before_file_assembly`), `read_segment` after (`prob_align`).
- One write of the particle field per iteration; the extra writes are listed in the policy's section 8: add
  none. Write segments (`write_segment_inside`) mid-run; after a full read call `update_projinfo`.
- A `<x>_exec` never calls `params%new` or reseeds; `run_<x>` and the strategies' `initialize` build them.
- Do not import a concrete commander module into a strategy module beyond those the policy's section 8 lists.
  `simple_refine_motion_model_strategy.f90` importing `commander_rec3D` is an unlisted deviation, ruled to be
  moved to the reconstruction strategy routine: do not copy.

## Where to start

- `simple_refine3D_strategy.f90`, `simple_refine2D_strategy.f90`: both modes; the distributed
  `execute_iteration` is the master's order. `simple_prob_align_strategy.f90`, `simple_volassemble_strategy.f90`,
  `simple_cavgassemble_strategy.f90`, `simple_strategy_sigma2.f90`: the master steps and their `run_<x>`.
- `simple_rec3D_*`: reconstruct3D (service, gridding, PCG, wave). Preprocessing strategies share one pattern:
  mic-row part documents merged by `merge_algndocs`; extract-type steps assemble part projects.

## Tests

- `build/bin/simple_test_exec test=unit_cart_align3D suite=local_staging` (refine3D_local staging),
  `unit_reconstruction suite=rec3d_backend`, `unit_motion suite=fraction_dose_metadata`;
  `lib_reconstruction suite=fractional_reconstruction`, `wave_reconstruction`, `fraction_dose_export` (nightly).
- No tester runs a master step on its own `(params, build)`: the workflow checks of `simple-main-strategies`
  cover them.

## See also

`doc/policies/layering_and_project_ownership_policy.md` (sections 3-5, 8), `doc/policies/3D/refine3D_policy.md`
sections 3-4, `doc/policies/2D/solve2D_policy.md` section 3; skills `simple-main-sigma2`, `simple-main-class`,
`simple-frac-update-trailing`.
