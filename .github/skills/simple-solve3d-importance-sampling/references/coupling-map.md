# Solve3D Importance Sampling and Fractional Updates

## Scope

This note explains how SIMPLE couples:

- outer particle-subset selection driven by `nsample` / `update_frac`
- inner importance-style sampling used by probabilistic orientation assignment
- downstream reconstruction behavior that depends on the realized subset

Keep these layers separate. The code is easier to understand when split into
policy, bookkeeping, assignment, and reconstruction.

## Core Contract

### 1. Solve3D owns the run's update fraction

`src/main/commanders/simple/simple_commanders_solve3D.f90` computes the
`update_frac` of the whole run before the stage loop.

Important details:

- `nsample` defaults to `NSAMPLE_SOLVE3D_DEFAULT` (10,000) and is the only
  public sample-size control; `nsample_start`/`nsample_stop` are refused.
- `update_frac = nsample * nstates / active`, capped by
  `UPDATE_FRAC_MAX = 0.9` so fractional update stays on.
- When `nsample / active > FULL_SAMPLE_SWITCH_FRAC` (0.9) the run forces full
  active sampling: the stage lines carry no `update_frac`, `nsample`, `fillin`
  or trailing reconstruction (`force_full_sampling_mode`).
- On the ini3D routes the initial subset is a greedy draw over the class units
  (`sample4update_class`, or `sample4update_cnt` without a 2D solution) with
  `updatecnt` set to 1 for the drawn particles.

That same command writes the sampling units to `CLASS_SAMPLING_FILE`
(`clssmp.bin`, `make_class_samples` in `simple_view_partition_sampling`): one
unit per selected 2D class with its `group` (0 under `balance=class`, the
class-average group under `balance=cavg`; class units under `none`, used only
by the initial draw), and prints the unit table with the expected visits per
particle before the first stage. Stages only read the file. Nothing in this
layer uses 3D maps, poses or projection directions.

### 2. The stage controller changes search mode and sampling together

`src/main/solve/simple_solve3D_controller.f90` (a submodule of
`simple_solve3D_utils`) emits the per-stage `refine3D` command line.

Search mode by stage:

- single-state: stages 1-2 `prob_neigh` with `prob_neigh_mode=shc`, stages
  3-5 `prob`, stages 6-8 `prob_neigh` with `prob_neigh_mode=state`
- multi-state (`nstates>1`, independent states): stages 1-2 `shc` (stage 2
  also exits on the scanned fraction, `scan_exit`), `prob` from stage 3, and
  `prob_neigh` with `prob_neigh_mode=state` from the first stage whose
  `nspace` reaches `NSPACE_MULTI_MAX` (2,500; stage 5,
  `multistate_prob_neigh_stage`); no stage exceeds 2,500 directions; by default
  the run stops after stage 5 (`nstages=5`, `lpstop=6`)
- `nstates` alone selects single-state or independent multi-state stage
  rules; solve3D has no `multivol_mode` input and no docked split

Sampling by stage:

- every stage carries the run's `update_frac` and `balance` (plus `nclust`
  for `cavg`)
- greedy (top 2D score per unit) until `STOCH_SAMPL_STAGE` (5; 4 for
  multi-state); then single-state uses `frac_best=0.5`, rising to 0.85 in
  stages 7-8, while multi-state uses `greedy_sampling=no` with `frac_best=1.0`
- `fillin=yes` on the last active stage of single-state runs only
- `trail_rec=yes` from `TRAILREC_STAGE_SINGLE` (5) single-state and
  `TRAILREC_STAGE_MULTI` (8) multi-state

Stage omission: a stage after the first, before the last and other than the
axis-search stage is skipped when the measured FSC=0.5 already resolves past
its planned limit (`solve3D_stage_behind_fsc05`, `>>> SOLVE3D OMITTING STAGE`).

This is the main place where the code couples search strategy, sampling size,
and reconstruction behavior.

## Bookkeeping Layer: `simple_oris`

### 3. `sampled` and `updatecnt` are the handshake

`src/main/ori/simple_oris_sampling.f90` defines the subset-management routines.

Key routines:

- `sample4update_class`: nested equal quota over the units of `CLASS_SAMPLING_FILE`
  (equal over groups, then over the units of a group, then lowest `updatecnt`);
  `class_sample_quotas` and `class_sample_sweep` give the expected quotas and the
  draws one full sweep takes (`max ceil(pop/quota)`)
- `sample4update_rescore`: the latest round again, `sampled`/`updatecnt` advanced
  once (the `cohort_sampling` schedule of `refine3D_states` frequency blocks)
- `sample4update_cnt`: prefers particles with lower `updatecnt`
- `sample4update_fillin`: fills in particles with low update history
- `sample4update_missing`: the particles never updated (`update_missing`)
- `sample4update_reprod`: reproduces the most recent sampled subset exactly
- `sample4update_updated`: retrieves all particles updated in a previous step

Interpretation:

- `sampled` identifies the current sampling round
- `updatecnt` counts how often a particle has been selected for update
- `incr_sampled_updatecnt` advances both the round marker and per-particle
  history when requested

### 4. Realized update fraction is recovered from bookkeeping

`src/main/ori/simple_oris_getters.f90` implements `get_update_frac`,
`get_group_update_counts` and `get_state_update_fracs`.

They do not return the command-line target. Per state (or globally):

- numerator `n`: active particles in the latest `sampled` round
- denominator `N`: active particles with `updatecnt > 0`

So the realized fraction is that of the currently updated subset within the
active updated pool. Under `m_estimator=flex` both are applied weight masses.

## Update Selection During Refinement

### 5. Non-probabilistic refine samples directly

`src/main/strategies/search/simple_matcher_smpl_and_lplims.f90` provides
`sample_ptcls4iter3D` and `sample_ptcls4update3D`.

Behavior:

- `update_missing`: `sample4update_missing`
- `fillin=yes`: every fifth iteration (`mod(which_iter,5) == 0`)
  `sample_ptcls4fillin` (`sample4update_fillin`)
- otherwise `sample_ptcls4update3D`: all active particles when `l_update_frac`
  is false; `sample4update_rescore` with `cohort_sampling=yes` and
  `which_iter /= startit`; `sample4update_class` for `balance=class|cavg`;
  `sample4update_cnt` for `balance=none`
- distributed: when `master_draws_sample3D` holds (fractional,
  non-probabilistic, not polish/missing/eval/sigma) the master draws once over
  the whole project and writes it; each worker reproduces its range with
  `sample4update_reprod`

### 6. Probabilistic refine samples once, then reproduces

`src/main/strategies/parallelization/simple_prob_align_strategy.f90` is the
key coupling point.

In `prob_align_exec` and `prob_align_neigh_exec` (run on the caller's
`(params, build)`; the distributed master's own, or via the
`commander_prob_align*` wrappers):

- `updatecnt` and `sampled` are cleared only at `startit == 1` on the first
  iteration
- the active subset is chosen by `sample_ptcls4update3D` or, on fill-in
  iterations, `sample_ptcls4fillin`
- the chosen subset is written to the project file
- partition jobs then compute probability tables only for that subset

In `exec_prob_tab` and `exec_prob_tab_neigh` (`simple_commanders_prob.f90`):

- no new sampling happens
- `sample4update_reprod` reproduces the subset created by `prob_align`
  (the worker stops if the field was never sampled)

In `src/main/strategies/search/simple_strategy3D_matcher.f90`
(`sample_particles_for_update`):

- probabilistic passes, the `pose_cont=yes` polish pass and distributed
  workers of a master-drawn sample use `sample4update_reprod`
- this ensures search and assignment operate on the exact same particles
- the reconstruction of the subset copies the images (or cropped records) the
  search pass read into the pass buffer instead of re-reading stacks

This is the most important contract in the whole workflow.

## Inner Importance Sampling

### 7. Probabilistic tables do their own candidate sampling

`src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90` and
`simple_eul_prob_tab_neigh.f90` implement the inner importance-style sampling
logic, with the samplers in `simple_eul_prob_tab_utils.f90`.

Important points:

- `angle_sampling` and `sample_bounded_dists` draw among the best candidates
  within the `prob_athres` angular bound
- `ref_assign` draws each state's projection with likelihood weights
  `exp(-dist)` over a top-K support (`sample_likelihood_dist`); the state label
  is the deterministic argmin over states
- `prob_neigh` restricts this to neighbourhood candidate sets
  (`prob_neigh_mode=state|geom|shc|snhc`)

This sampling is not the same as the outer `update_frac` selection:

- outer selection chooses which particles participate this iteration
- inner sampling chooses which candidate orientations or references are explored
  for each participating particle

## Reconstruction Consequences

### 8. Volume assembly consumes the realized fraction

`volassemble_exec` (`src/main/strategies/parallelization/simple_volassemble_strategy.f90`)
does not choose particles.

Its coupling to sampling appears when `trail_rec=yes`:

- the realized per-state fraction `f = n/N` is always computed
- the applied weight `u` is `ufrac_trec` only for a single-state run with a
  parsed override (`l_ufrac_trec_defined`), otherwise `f`
- current partial accumulators are scaled by `u/f` and the trailing chain by
  `(1-u)*N/M` (`population_blend_weights`); the chain records `M = N`

So the realized particle subset affects how strongly the new reconstruction
replaces the previous one. The details are in
`.github/skills/simple-frac-update-trailing/`.

## Practical Reading Order

For explanations:

1. `simple_commanders_solve3D.f90`
2. `simple_solve3D_controller.f90`
3. `simple_matcher_smpl_and_lplims.f90`
4. `simple_oris_sampling.f90`
5. `simple_prob_align_strategy.f90`
6. `simple_strategy3D_matcher.f90`
7. `search/probabilistic/simple_eul_prob_tab*.f90`
8. `simple_volassemble_strategy.f90`

For modifications:

1. change policy in solve3D or stage-controller code
2. confirm `simple_oris` bookkeeping still matches the intended subset semantics
3. confirm probabilistic and master-drawn paths reproduce rather than resample
4. confirm trailing reconstruction still uses the intended realized fraction
5. confirm reconstruction still takes the pass buffer's images or records
