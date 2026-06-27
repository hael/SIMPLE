---
name: simple-abinitio3d-importance-sampling
description: Use when working on SIMPLE abinitio3D sampling policy, update_frac or nsample behavior, prob_align/prob_tab sampling reuse, sampled/updatecnt bookkeeping, online matcher reconstruction I/O, or volassemble trailing-update weighting.
---

# SIMPLE Abinitio3D Importance Sampling

Use this skill to reason about the full contract between ab initio stage policy
and per-iteration particle sampling. Treat the implementation as a pipeline:

1. `exec_abinitio3D` decides the initial and stage-wise fractional-update policy.
2. `simple_oris` stores the chosen subset through `sampled` and `updatecnt`.
3. `prob_align` or `refine3D_exec` reuses that exact subset instead of resampling independently.
4. `simple_eul_prob_tab` performs importance-style reference and in-plane selection only inside the sampled subset.
5. `volassemble` does not resample particles, but it does consume the realized update fraction to mix trailing reconstructions.

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
5. Use exact field names such as `update_frac`, `sampled`, `updatecnt`,
   `frac_best`, `fillin`, `balance`, and `ufrac_trec`.

## Working Rules

- Do not describe `volassemble` as the module that decides which particles get
  updated. It only consumes the result later for trailing reconstruction
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
- When explaining behavior, call out whether the code path is initial ab initio
  seeding, normal iteration, final-stage fill-in, or trailing reconstruction.
- Preserve the online matcher single-read contract. If reconstruction or
  restoration is active, the matcher should reuse the already-read batch images
  after assignment. Do not add a second particle-stack pass to reconstruct the
  sampled subset unless the repository policy explicitly changes.

## File-Driven Workflow

### Stage policy and initial subset

Inspect:

- `src/main/commanders/simple/simple_commanders_abinitio.f90`
- `src/main/simple_abinitio_utils.f90`
- `src/main/simple_abinitio_controller.f90`

Focus on `update_frac`, `nsample*`, `UPDATE_FRAC_MAX`, `sample4update_class`,
`shc_smpl`, `prob`, `prob_neigh`, `frac_best`, `fillin`, and `trail_rec`.

### Particle-subset bookkeeping

Inspect:

- `src/main/ori/simple_oris_sampling.f90`
- `src/main/ori/simple_oris_getters.f90`
- `src/main/ori/simple_oris.f90`

Focus on `sample4update_class`, `sample4update_cnt`,
`sample4update_fillin`, `sample4update_reprod`, `sample4update_updated`,
`incr_sampled_updatecnt`, and `get_update_frac`.

### Probabilistic sampling reuse

Inspect:

- `src/main/commanders/simple/simple_commanders_prob.f90`
- `src/main/simple_eul_prob_tab.f90`
- `src/main/simple_eul_prob_tab_neigh.f90`
- `src/main/strategies/search/simple_strategy3D_matcher.f90`

Preserve the sample-once-and-reproduce contract between `prob_align`,
`prob_tab`, and `refine3D_exec`.

### Trailing reconstruction

Inspect:

- `src/main/commanders/simple/simple_commanders_rec_distr.f90`
- `src/main/volume/simple_reconstructor_eo.f90`

Remember that `volassemble` blends current and previous even/odd reconstructions
using the realized update fraction. It should not choose a new particle subset.

## Language To Prefer

- "outer fractional-update selection"
- "inner probabilistic candidate sampling"
- "sample once, then reproduce the subset"
- "`sampled` current round and `updatecnt` persistent history"
- "realized update fraction consumed by trailing reconstruction"
