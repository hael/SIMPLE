# Importance Sampling and Fractional Update Policy

## Scope

This document defines the current policy for how SIMPLE couples fractional particle updates to importance-sampled 3D search in:

- `abinitio3D`
- `refine3D`
- probabilistic pre-alignment
- orientation bookkeeping in `simple_oris`
- trailing reconstruction in `volassemble`

This is a workflow and ownership document, not just an implementation note. It explains which layer decides the active particle subset, which layer samples candidates inside that subset, and which layer consumes the realized update fraction later in reconstruction.

## Core design rule

The 3D workflow has two distinct sampling layers that must not be conflated:

- outer fractional-update sampling chooses which particles participate in the current iteration
- inner importance sampling chooses which reference and in-plane candidates are explored for those participating particles

In other words:

- `update_frac` or `nsample*` decides the active particle subset
- probabilistic tables and search strategies decide candidate sampling within that subset

That split should be preserved in both code and documentation.

## Public policy

The current 3D policy is:

1. choose an active particle subset for the iteration
2. record that subset in `ptcl3D` through `sampled` and `updatecnt`
3. if probabilistic pre-alignment is enabled, generate assignments only for that recorded subset
4. run the main matcher or search pass on the same subset
5. write partial reconstructions from those updated particles
6. let `volassemble` assemble volumes and, if trailing reconstruction is enabled, weight old versus new volumes using the realized update fraction

The user-facing controls that shape this policy include:

- `update_frac`
- `nsample`
- `nsample_start`
- `nsample_stop`
- `fillin`
- `balance`
- `frac_best`
- `trail_rec`
- `ufrac_trec`

## Core vocabulary

The following terms should be used consistently:

- outer fractional-update subset: the particle set selected for the current iteration
- current sampling round: the most recent subset marked through `sampled`
- update history: the per-particle selection history tracked through `updatecnt`
- probabilistic pre-alignment: the phase that computes assignment tables before the main matcher
- realized update fraction: the fraction recovered from current `sampled` versus active updated particles, not merely the CLI target

Avoid describing the whole workflow as one generic "importance sampling" mechanism. There are two different sampling contracts in play.

## Ownership

### `simple_commanders_abinitio`

`src/main/commanders/simple/simple_commanders_abinitio.f90` owns:

- initial ab initio sampling policy
- conversion from `nsample` or `nsample*` inputs into `update_frac`
- the initial class-balanced subset used to seed iterative ab initio
- class-sampling-statistics generation through `CLASS_SAMPLING_FILE`

It is the top-level owner of ab initio sampling intent.

### `simple_abinitio_controller`

`src/main/simple_abinitio_controller.f90` owns:

- stage-specific `refine3D` command generation inside ab initio
- dynamic `update_frac` scheduling across stages
- transitions between `shc_smpl`, `prob`, and `prob_neigh`
- stage-specific control of `frac_best`, `fillin`, and `trail_rec`

This layer couples search regime and fractional-update policy at the workflow level.

### `simple_matcher_smpl_and_lplims`

`src/main/strategies/search/simple_matcher_smpl_and_lplims.f90` owns:

- per-iteration selection of the active particle subset during 3D refinement
- dispatch between class-balanced sampling, update-count-biased sampling, fill-in sampling, or all-particle sampling

This is the main execution helper for outer subset selection during refinement.

### `simple_oris`

`src/main/ori/simple_oris_sampling.f90` and `src/main/ori/simple_oris_getters.f90` own:

- representation of the current sampling round through `sampled`
- representation of update history through `updatecnt`
- exact reproduction of a previous sampled subset
- recovery of the realized update fraction

This is the bookkeeping contract that all higher workflow layers depend on.

### `simple_commanders_prob`

`src/main/commanders/simple/simple_commanders_prob.f90` owns:

- sampling the subset for probabilistic pre-alignment
- persisting that subset before probability-table generation
- ensuring partition-local probability-table work is performed only for the chosen subset

This layer should sample once, then reuse.

### `simple_eul_prob_tab` and `simple_eul_prob_tab_neigh`

`src/main/simple_eul_prob_tab.f90` and `src/main/simple_eul_prob_tab_neigh.f90` own:

- importance-style candidate sampling over references and in-plane angles
- transformation of correlation or distance values into a candidate-selection distribution
- hard assignment output for downstream consumption

These modules do not decide the outer particle subset.

### `simple_strategy3D_matcher`

`src/main/strategies/search/simple_strategy3D_matcher.f90` owns:

- matcher execution on the active particle subset
- reuse of probabilistic assignments when present
- writing partial reconstructions from the updated subset

It is the main particle-domain execution engine after subset selection has already been decided.

### `volassemble`

`src/main/commanders/simple/simple_commanders_rec_distr.f90` owns:

- volume assembly from partial reconstructions
- trailing-reconstruction weighting
- consumption of the realized update fraction through `get_update_frac` or explicit `ufrac_trec`

It does not own outer subset selection and should not be described that way.

## Fractional-update policy in `abinitio3D`

### Initial policy

At ab initio startup, `exec_abinitio3D` computes a base `update_frac` from one of four sources:

- explicit `nsample`
- explicit `update_frac`
- dynamic lower and upper bounds through `nsample_start` and optional `nsample_stop`
- default sample-count bounds

The resulting value is capped at `UPDATE_FRAC_MAX = 0.9`.

This cap is intentional. The ab initio workflow is designed to remain in fractional-update mode rather than silently drifting into full-update behavior.

### Initial subset construction

Before iterative refinement starts, ab initio may create a class-balanced subset through `sample4update_class`.

That subset is not just a local temporary choice. It seeds the particle set that reconstruction and later search stages build on, and it initializes the `simple_oris` bookkeeping that later stages depend on.

### Stage policy

Inside `simple_abinitio_controller.f90`, the stage controller changes:

- search mode
- dynamic update fraction
- greediness of class-balanced selection
- whether fill-in is enabled
- whether trailing reconstruction is enabled

Current high-level stage policy is:

- stages 1 to 4: `shc_smpl`
- stages 5 to 6: `prob`
- stages 7 to 8: `prob_neigh`

This is an important architectural point: ab initio does not bolt probabilistic importance sampling onto a fixed outer update schedule. It co-evolves the search mode and the update policy stage by stage.

## Outer subset-selection policy during 3D refinement

### Normal update path

`sample_ptcls4update3D` applies the following policy:

- if `l_update_frac` is false, select all active particles
- if `balance=yes`, use class-balanced sampling through `sample4update_class`
- otherwise use update-count-aware sampling through `sample4update_cnt`

This means the normal 3D path prefers either:

- class-balance preservation
- under-updated particles

depending on the workflow settings.

### Fill-in path

When `fillin=yes`, refinement may switch to `sample4update_fillin`.

This is a different policy from normal balanced or count-biased update. The goal is not to preserve the current exploration distribution, but to fill gaps in update coverage late in the workflow.

That distinction should be preserved when reviewing or modifying the code.

## `simple_oris` bookkeeping policy

### `sampled`

`sampled` identifies the current sampling round.

Policy meaning:

- all particles with the latest `sampled` marker belong to the current active subset
- reproducibility of a subset across workflow stages should happen by reproducing the latest `sampled` round

### `updatecnt`

`updatecnt` tracks cumulative update history.

Policy meaning:

- particles with lower `updatecnt` are eligible for preferential selection in update-count-biased paths
- fill-in logic uses update history rather than just the last sampled round
- realized update fraction is interpreted relative to the updated active pool

### Realized update fraction

`simple_oris_getters.f90` computes `get_update_frac` from:

- the latest `sampled` round
- particles with `updatecnt > 0`
- particles with active `state > 0`

This is a crucial policy point:

- the workflow target may be `update_frac = x`
- the reconstruction contract later uses the realized fraction recovered from bookkeeping

So the bookkeeping state is authoritative for downstream trailing-reconstruction behavior.

## Probabilistic pre-alignment policy

### Sample once, then reuse

When `l_prob_align_mode` is enabled, `simple_commanders_prob.f90` samples the active particle subset before probability-table generation.

The policy is:

1. choose the active subset through the same outer-selection helpers used elsewhere
2. persist that subset into the project
3. run partition-local probability-table generation only for that subset
4. aggregate tables and emit a single assignment map

This policy must remain sample-once-and-reuse, not sample-again-later.

### Probability-table generation is not a second outer sampler

`exec_prob_tab` and `exec_prob_tab_neigh` reproduce the existing subset through `sample4update_reprod`.

That is the correct design. These steps should not independently pick a new particle subset, because that would break the contract between:

- pre-alignment
- assignment
- matcher execution

### Matcher reuse

When probabilistic mode is active, `refine3D_exec` also reproduces the same sampled subset rather than choosing a new one.

This guarantees that:

- the particles scored in the probability tables
- the particles updated in the matcher
- the particles contributing to partial reconstructions

all refer to the same current sampling round.

This is one of the most important invariants in the entire workflow.

## Inner importance-sampling policy

Within the already chosen active particle subset, `simple_eul_prob_tab*.f90` performs importance-style candidate selection.

The current policy is:

- compute candidate scores over references and in-plane rotations
- transform those values into a positive candidate-selection space
- use `angle_sampling` and `greedy_sampling` to choose likely candidates rather than relying only on strict top-1 selection
- perform hard assignment after the candidate-search phase

This layer is allowed to be probabilistic or greedy in candidate exploration, but it must remain conceptually downstream of outer subset selection.

That distinction should be explicit in documentation:

- outer subset selection is particle-level participation policy
- inner importance sampling is per-particle candidate-exploration policy

## Reconstruction and trailing-update policy

### Partial reconstructions

`simple_strategy3D_matcher.f90` writes partial reconstructions only for the particles in the current active subset.

That means volume assembly is already downstream of outer subset choice before `volassemble` even begins.

### `volassemble`

When `trail_rec=yes`, `volassemble` mixes previous and current even and odd volumes using `update_frac_trail_rec`.

Policy order:

- use explicit `ufrac_trec` if provided
- otherwise recover the realized update fraction from `ptcl3D` via `get_update_frac`

The weighting then becomes:

- previous contribution: `1 - update_frac_trail_rec`
- current contribution: `update_frac_trail_rec`

This is not a search policy. It is a reconstruction-consumption policy.

### Implication

Changes to sampling or bookkeeping semantics can change:

- which particles are updated
- how many particles are updated
- how strongly the current reconstruction overrides the previous reconstruction

So sampling-policy changes are also reconstruction-policy changes, even if `volassemble` itself is untouched.

## Invariants to preserve

- The outer particle subset must be chosen before probabilistic table generation.
- `prob_tab` and matcher execution must reproduce the same sampled subset in probabilistic mode.
- `sampled` must continue to represent the current sampling round.
- `updatecnt` must continue to represent cumulative update history.
- `get_update_frac` must continue to reflect realized rather than nominal update behavior.
- `volassemble` must remain a consumer of update-fraction state, not the producer of particle-selection policy.
- Documentation must distinguish outer subset selection from inner candidate importance sampling.

## Recommended reading order

For policy work:

1. `src/main/commanders/simple/simple_commanders_abinitio.f90`
2. `src/main/simple_abinitio_utils.f90`
3. `src/main/simple_abinitio_controller.f90`
4. `src/main/strategies/search/simple_matcher_smpl_and_lplims.f90`
5. `src/main/ori/simple_oris_sampling.f90`
6. `src/main/ori/simple_oris_getters.f90`
7. `src/main/commanders/simple/simple_commanders_prob.f90`
8. `src/main/simple_eul_prob_tab.f90`
9. `src/main/simple_eul_prob_tab_neigh.f90`
10. `src/main/strategies/search/simple_strategy3D_matcher.f90`
11. `src/main/commanders/simple/simple_commanders_rec_distr.f90`

## Guidance for future modifications

Changes are consistent with current policy when they:

- make the outer subset policy more explicit without duplicating it across layers
- preserve sample-once-and-reuse behavior in probabilistic mode
- improve observability of realized update fraction
- separate candidate-sampling logic from particle-subset logic more clearly
- keep trailing-reconstruction weighting tied to bookkeeping rather than ad hoc estimates

Changes are inconsistent with current policy when they:

- let probability-table generation silently resample particles
- let matcher execution use a different subset than pre-alignment
- repurpose `sampled` or `updatecnt` without updating all dependent workflow stages
- describe `volassemble` as owning particle-selection policy
- blur the distinction between outer fractional update and inner candidate importance sampling
