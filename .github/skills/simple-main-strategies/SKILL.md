---
name: simple-main-strategies
description: Use when working in SIMPLE's src/main/strategies subsystem, including the shared-memory and distributed-master parallelization strategies, the strategy routines the refine3D/refine2D masters run on their own project (volume and class-average assembly, prob_align, sigma2 commit), 2D and 3D search strategies, matchers, probability tables, batching, reconstruction helpers, and the execution-policy layer between commanders and low-level scientific objects.
---

# SIMPLE `src/main/strategies`

This folder is the execution-policy and search-engine layer.

## Read First

- `parallelization/simple_refine3D_strategy.f90`
- `parallelization/simple_refine2D_strategy.f90`
- `parallelization/simple_preprocess_strategy.f90`
- `search/simple_strategy2D_matcher.f90`
- `search/simple_strategy3D_matcher.f90`
- `search/simple_matcher_3Drec.f90`

## Subareas

- `parallelization/` handles shared-memory vs distributed policy and workflow iteration (`create_<x>_strategy(cline)` picks the shared-memory or distributed type)
- `search/` holds 2D/3D alignment strategy families, matchers, batching, and search helpers
- `search/probabilistic/` holds the probability tables (`simple_eul_prob_tab*.f90`) and the `prob` 2D/3D strategies

## Master-Owned Project (distributed refine3D/refine2D)

- The distributed master owns its project. Its steps are strategy routines on
  the master's `(params, build)`:
  `simple_volassemble_strategy.f90` (`volassemble_exec`),
  `simple_prob_align_strategy.f90` (`prob_align_exec`,
  `prob_align_neigh_exec`, `prob_align2D_exec`),
  `simple_cavgassemble_strategy.f90` (`cavgassemble_exec`),
  `simple_strategy_sigma2.f90` (`commit_sigma2_exec`). Each also has a
  `run_<x>(cline)` entry that builds its own `(params, build)`; the commanders
  are thin wrappers of it.
- A `<x>_exec` routine runs at the `nthr` it is given; the caller saves its
  thread count (`omp_get_max_threads()`, `nthr_glob`) and restores it after.
- These strategy modules never import the commander module that wraps them
  (circular); a body that needs an in-process commander takes it as a
  `commander_base` dummy (`run_prob_align(cline, xprob_tab)`).
- The refine3D master writes its particle field once before the workers and once
  per iteration (once more after the merge when a polish pass follows):
  `merge_algndocs(..., write=.false.)` merges in memory and gridding assembly
  passes `write=.false.`. The PCG master (`rec_backend=pcg`) writes the field once
  more, at the end of its call (it has no `write` flag). The distributed refine2D
  master writes `ptcl2D` in `merge_algndocs`, after assembly and at the end of the
  iteration (three times, four in probabilistic mode); 2D assembly updates only the
  segments it changes. These are deviations from one write per iteration, listed in
  `doc/policies/layering_and_project_ownership_policy.md` section 8.
  Keep `update_projinfo` when the refine3D master builds its toolbox
  (`materialize_reprojection_model_master`), so writes go to the run's own
  project. Workers read only the file.
- 2D mirrors 3D: the matcher (`refine2D_exec`) searches, accumulates class
  sums and writes orientations; the strategy restores class averages
  (`restore_cavgs_inmem`) and decides convergence in both modes
  (`check_refine2D_convergence`).
- Shared-memory refine2D/refine3D still run these steps through the commanders
  (`execute(cline)`) or `run_<x>(cline)`, which read the project file; write
  the particle field first.
- Contracts: `doc/policies/3D/refine3D_policy.md` (section 4) and
  `doc/policies/2D/solve2D_policy.md` (section 3).

## Continuous Pose Refinement

- `search/simple_strategy3D_cont.f90` is the Cartesian 3D strategy of `refine=cont` (one `cartft_pose_opt` transaction per particle; `athres_cont` bounds the rotation).
- `parallelization/simple_refine3D_strategy.f90` schedules the polish (`pose_cont=yes`, `polish_follows`): a `refine=cont` pass after each discrete pass over the same sample, with the particle field written before and after it. `cont_route` is an internal test seam.
- Before a gridding assembly the shared-memory strategy writes the particle field (`publish_before_file_assembly`): its `run_volassemble` builds from the project file.

## Working Rule

This folder is usually the right home for “how a workflow runs” but not “what the raw data structure means.”
