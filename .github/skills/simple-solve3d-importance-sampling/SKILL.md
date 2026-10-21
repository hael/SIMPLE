---
name: simple-solve3d-importance-sampling
description: Use when working on SIMPLE solve3D sampling policy, nsample/update_frac behavior and the full-sampling switch, stage-wise sampling (balance units, greedy/stochastic, frac_best, fillin), prob_align/prob_tab sampling reuse and the distributed master's sample draw, sampled/updatecnt bookkeeping, inner probability-table candidate sampling, or the realized update fraction consumed by trailing reconstruction.
---

# SIMPLE Solve3D Importance Sampling

`solve3D` is de novo map determination: it couples ab initio 3D reconstruction
with initial 3D refinement in one stage schedule.

Use this skill to reason about the full contract between solve3D stage policy
and per-iteration particle sampling. Treat the implementation as a pipeline:

1. `exec_solve3D` sets the run's fixed `update_frac` from `nsample` (or forces
   full sampling), writes the sampling units, and the stage controller emits
   the per-stage sampling controls.
2. `simple_oris` stores the chosen subset through `sampled` and `updatecnt`.
3. `prob_align_exec` (probabilistic passes) or the distributed master
   (non-probabilistic fractional passes) draws once; table workers and
   `refine3D_exec` reproduce that exact subset instead of resampling.
4. `simple_eul_prob_tab*` performs importance-style reference and in-plane selection only inside the sampled subset.
5. Volume assembly does not resample particles, but it consumes the realized per-state update fraction to blend the trailing chain.

## Quick Start

Read [references/coupling-map.md](./references/coupling-map.md) first.

Then answer the user request in this order:

1. Identify whether the question is about policy selection, subset bookkeeping,
   probabilistic assignment, or reconstruction consequences.
2. Name the owner file for each step rather than treating "importance sampling"
   as one monolith.
3. State whether the subset is being newly sampled, reproduced from prior
   sampling, or inferred from `updatecnt`.
4. Separate particle-domain search decisions from volume-domain trailing-reconstruction effects.
5. Use exact field names such as `nsample`, `update_frac`, `sampled`,
   `updatecnt`, `frac_best`, `greedy_sampling`, `fillin`, `balance`
   (`none|class|cavg`), `nclust`, `cohort_sampling`, and `ufrac_trec`.

## Working Rules

- Do not describe volume assembly as the module that decides which particles
  get updated. It only consumes the result later for trailing reconstruction
  weighting.
- Do not describe `prob_align` as introducing a second independent sampling
  policy. It samples once, writes the subset to `ptcl3D`, and downstream
  `prob_tab` plus `refine3D_exec` reproduce that subset.
- Do not conflate probabilistic candidate selection with the outer
  fractional-update schedule. The outer schedule decides which particles
  participate; `simple_eul_prob_tab` decides which references and in-plane
  angles are sampled for those particles.
- Treat `sampled` as the marker for the current sampling round and `updatecnt`
  as persistent history across rounds.
- `solve3D` takes `nsample` only: `nsample_start`/`nsample_stop` are refused
  and the update fraction does not change between stages.
- There is no docked solve3D split, sticky class sampling or `sampled_only`
  cohort restriction; multi-state solve3D is independent from the start, and
  the only cohort rule is `refine3D_states` `cohort_sampling`.
- The sampler never uses 3D maps, poses or projection directions; sampling
  units come from 2D classes only.
- When explaining behavior, call out whether the code path is initial solve3D
  seeding, normal iteration, last-stage fill-in, or trailing reconstruction.
- Preserve the reconstruction read contract: the reconstruction of the sampled
  subset copies the images (or cropped records) the search pass already read
  into the pass buffer. Do not add a second particle-stack pass unless the
  repository policy explicitly changes.

## File-Driven Workflow

### Stage policy and initial subset

Inspect:

- `src/main/commanders/simple/simple_commanders_solve3D.f90`
- `src/main/solve/simple_solve3D_utils.f90`
- `src/main/solve/simple_solve3D_controller.f90`

Focus on `nsample` (`NSAMPLE_SOLVE3D_DEFAULT`), `UPDATE_FRAC_MAX`,
`FULL_SAMPLE_SWITCH_FRAC`, `force_full_sampling_mode`, `sample4update_class`,
`prob_neigh` (and `prob_neigh_mode`), `prob`, `NSPACE_MULTI_MAX`,
`greedy_sampling`, `frac_best`, `fillin`, and `trail_rec`.

### Particle-subset bookkeeping

Inspect:

- `src/main/ori/simple_oris_sampling.f90`
- `src/main/ori/simple_oris_getters.f90`
- `src/main/ori/simple_oris.f90`

Focus on `sample4update_class`, `sample4update_cnt`,
`sample4update_fillin`, `sample4update_missing`, `sample4update_reprod`,
`sample4update_rescore`, `sample4update_updated`, `incr_sampled_updatecnt`,
`get_update_frac`, `get_state_update_fracs`, `get_group_update_counts`, and the
unit-table helpers `class_sample_quotas` and `class_sample_sweep`. The sampling
units (`balance=none|class|cavg`) are produced by
`src/main/strategies/search/simple_view_partition_sampling.f90`.

### Probabilistic sampling reuse

Inspect:

- `src/main/strategies/parallelization/simple_prob_align_strategy.f90`
- `src/main/commanders/simple/simple_commanders_prob.f90` (table workers)
- `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90`
- `src/main/strategies/search/probabilistic/simple_eul_prob_tab_neigh.f90`
- `src/main/strategies/search/simple_matcher_smpl_and_lplims.f90`
- `src/main/strategies/search/simple_strategy3D_matcher.f90`

Preserve the sample-once-and-reproduce contract between `prob_align`,
`prob_tab`, and `refine3D_exec`, and between the distributed master's draw
(`master_draws_sample3D`) and its workers.

### Trailing reconstruction

Inspect:

- `src/main/strategies/parallelization/simple_volassemble_strategy.f90`
  (`volassemble_exec`; trailing blend, chain and gridding pair I/O inside
  `restore_state_from_parts`)
- `src/main/volume/simple_halfmap_diagnostics.f90`

Remember that volume assembly blends the current partial accumulators with the
per-state trailing chain using the realized update fraction (population rule).
It should not choose a new particle subset.

## Language To Prefer

- "outer fractional-update selection"
- "inner probabilistic candidate sampling"
- "sample once, then reproduce the subset"
- "`sampled` current round and `updatecnt` persistent history"
- "realized update fraction consumed by trailing reconstruction"
