# Refine3D Policy

## Scope

This document defines the current architectural policy for `refine3D` in SIMPLE.

It covers:

- command and strategy ownership
- probabilistic pre-alignment and particle-update flow
- the boundary between particle-domain work and volume-domain work
- the role of `volassemble` in shared-memory and distributed execution

This document is intentionally broader than a single implementation note. It is the top-level policy for how the `refine3D` workflow is divided conceptually and in code.

## Core design rule

`refine3D` is a layered fixed-point workflow with a strict boundary:

- particle-domain work stays in probabilistic alignment, search strategies, matcher preparation, orientation update, and partial reconstruction
- assembled-volume work stays in `volassemble` and its postprocessing policy/helpers

In other words:

- `refine3D` decides and updates particle poses
- `volassemble` builds and postprocesses state volumes from those decisions

That split should be preserved.

## Public policy

The current public workflow for `refine3D` is:

1. initialize execution mode and iteration state
2. optionally run probabilistic pre-alignment
3. run the particle-update matcher/search pass
4. write partial reconstructions
5. assemble and postprocess state volumes in `volassemble`
6. persist updated orientation and volume artifacts for the next iteration or downstream programs

The same policy applies to:

- shared-memory `refine3D`
- distributed-master `refine3D`

Only the process-launch mechanism differs between those routes. The scientific workflow and ownership boundaries should remain aligned.

## Execution ownership

### `simple_commanders_refine3D`

`simple_commanders_refine3D.f90` owns:

- the entrypoint for `refine3D`
- top-level defaults
- selection of execution strategy through `create_refine3D_strategy`

It should remain thin. It is not the place for low-level search logic or assembled-volume policy.

### `simple_refine3D_strategy`

`simple_refine3D_strategy.f90` owns:

- iteration control
- shared-memory versus distributed-master execution policy
- scheduler interaction
- iteration counters and run-finalization bookkeeping
- orchestration of pre-alignment, matcher execution, and `volassemble`

This layer may thread command-line state and workflow state across steps, but it should not absorb numerical postprocessing logic.

### `simple_commanders_prob` and probability-table modules

`simple_commanders_prob.f90` together with:

- `simple_eul_prob_tab.f90`
- `simple_eul_prob_tab_neigh.f90`

own the probabilistic pre-alignment phase.

That phase is responsible for:

- sampling particles for update
- generating partition-local probability tables
- aggregating those tables across partitions
- writing a single assignment artifact consumed by the matcher

This phase is particle-domain work, not volume-domain work.

### `simple_strategy3D_matcher`

`simple_strategy3D_matcher.f90` owns the core particle-update pass through `refine3D_exec`.

That includes:

- reference preparation
- search-strategy dispatch
- candidate evaluation
- orientation/state/in-plane/shift update
- Euclidean sigma update during search when applicable
- writing partial reconstructions for downstream assembly

This is the execution center of particle-domain refinement.

### `volassemble`

`simple_commanders_rec_distr.f90`, specifically `commander_volassemble`, is the single execution point for shared assembled-volume work used by `refine3D`.

That includes:

- reduction of partial reconstructions
- even/odd handling
- gridding correction
- merged-volume creation
- automask production
- nonuniform filtering

This policy now applies to both shared-memory and distributed `refine3D`.

## Particle-domain versus volume-domain boundary

The most important architectural boundary in `refine3D` is the distinction between particle-domain and volume-domain work.

### Particle-domain work

Particle-domain work includes:

- particle sampling for update
- probabilistic candidate generation
- state/projection/in-plane/shift assignment
- reference matching
- sigma updates during search
- writing partial reconstructions

The key code owners are:

- `simple_commanders_prob.f90`
- `simple_eul_prob_tab*.f90`
- `simple_strategy3D_matcher.f90`
- `simple_strategy3D_*.f90`
- `simple_matcher_3Drec.f90`

### Volume-domain work

Volume-domain work includes:

- summing partial reconstructions into state volumes
- even/odd averaging and merged-volume generation
- gridding correction
- automask generation and reuse
- nonuniform filtering
- artifact writing that downstream FSC and reference preparation consume

The key code owners are:

- `simple_commanders_rec_distr.f90`
- `simple_volume_postprocess_policy.f90`
- `simple_reconstructor_eo.f90`

Mixing these responsibilities makes the workflow harder to reason about and risks divergence between execution modes.

## Probabilistic refinement policy

`refine3D` supports probabilistic search-related modes, but the workflow is not a monolithic soft-assignment EM implementation.

The current policy is:

- probabilistic pre-alignment may build candidate tables and select assignments before the main matcher pass
- the main matcher consumes that information and applies hard orientation/state updates to the working orientation model
- downstream reconstruction and volume assembly operate on those assigned particle updates

This distinction should be preserved in documentation and code review:

- "probabilistic pre-alignment" is accurate
- "hard-assignment particle update" is accurate
- "fully soft volume-integrated EM" is not an accurate description of the current `refine3D` implementation

## Artifact policy

The workflow depends on conventional artifacts that act as handoff points between layers.

Examples include:

- assignment maps written by probabilistic alignment
- partition-local partial reconstructions
- state volumes and even/odd volumes
- state-specific automasks such as `automask3D_stateNN.mrc`

These artifacts are not incidental. They are part of the current contract between workflow phases and should be treated as stable unless a coordinated policy change is made.

## Why `volassemble` is the right execution point

Keeping assembled-volume work in `volassemble` is good policy because it:

- removes duplicated post-assembly logic from separate `refine3D` routes
- keeps the expensive shared-memory volume work in one place
- improves behavioral parity between shared-memory and distributed execution
- keeps `refine3D_strategy` focused on orchestration rather than low-level volume handling
- makes benchmarking more meaningful because the same step owns the same costs

This is the right direction and should not be rolled back.

## Current implementation status

The current implementation has already moved in the right direction.

In particular:

- `simple_refine3D_strategy.f90` now orchestrates probabilistic pre-alignment, matcher execution, and `volassemble`
- `simple_strategy3D_matcher.f90` is the core particle-update path
- `simple_volume_postprocess_policy.f90` has started to factor automask and nonuniform-filter decisions out of `volassemble`

That newer policy helper is important. It shows the intended refinement of the architecture:

- `volassemble` remains the execution site
- dedicated helpers own more of the policy branching

## Remaining rough edges

The main rough edge is not the direction of the architecture, but the remaining amount of policy detail still carried in `volassemble`.

Without helper objects, `volassemble` tends to accumulate decisions about:

- when automask regenerates
- whether an existing mask is compatible
- whether nonuniform filtering should use a state mask or spherical fallback
- which artifacts should be reused versus regenerated

That is too much policy density for one commander.

## Recommended boundary going forward

The recommended split is:

- `simple_commanders_refine3D.f90`
  Entry point and defaults only.
- `simple_refine3D_strategy.f90`
  Iteration control, execution-mode orchestration, scheduler interaction, and workflow threading.
- probability-table and matcher/search modules
  Particle-domain candidate generation, assignment, and partial reconstruction.
- `simple_volume_postprocess_policy.f90`
  Decide automask cadence, compatibility, nonuniform-mask precedence, and related artifact policy.
- `volassemble`
  Execute the plan for one iteration and state set.
- numerical modules
  Implement automasking, FSC masking, and nonuniform filtering details.

## Near-term improvements

The following improvements are consistent with current policy:

- preserve exact policy values end-to-end, including `automsk=tight`
- keep one compatibility rule for both state-mask production and state-mask consumption
- make automask cadence explicit and testable
- continue moving postprocessing branching into `simple_volume_postprocess_policy.f90`
- add end-to-end regression coverage that compares shared-memory and distributed `refine3D` artifact sets

## Rules to preserve during refactors

- Do not move assembled-volume postprocessing back into `refine3D_strategy`.
- Do not treat `volassemble` as a distributed-only helper.
- Do not merge probabilistic particle-update logic with volume postprocessing logic.
- Do not describe the current `refine3D` implementation as if it were a single undifferentiated Bayesian engine.
- Do preserve parity between shared-memory and distributed workflows.

## Workflow summary

1. `commander_refine3D` selects the execution strategy.
2. `refine3D_strategy` manages iteration state and execution mode.
3. probabilistic pre-alignment may sample particles and write an assignment map.
4. `refine3D_exec` performs particle-domain search, update, and partial reconstruction.
5. `volassemble` assembles state volumes and performs shared volume-domain postprocessing.
6. output artifacts are persisted for the next iteration and for downstream FSC/reference consumers.

That is the current `refine3D` policy and the architectural model future changes should follow.
