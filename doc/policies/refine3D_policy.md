# Refine3D Volume-Assembly Policy

## Scope

This document describes the current `refine3D` policy for volume-domain work and the intended evolution of that policy.

The key design shift is:

- particle-orientation work stays in the matcher/refinement path
- expensive shared-memory work on assembled volumes moves into `volassemble`

## Current policy

`refine3D` should use `volassemble` as the single execution point for volume-domain operations that are shared between execution modes.

That includes:

- assembly of partial reconstructions
- even/odd handling
- gridding correction
- merged-volume creation
- automask generation
- nonuniform filtering

This policy now applies to both:

- shared-memory `refine3D`
- distributed master `refine3D`

## Why this evolution is good

- It removes duplicated post-assembly logic from different refine3D paths.
- It keeps the expensive volume work in one place, which improves parity between execution modes.
- It makes workflow behavior easier to benchmark because the same step owns the same costs.
- It moves refine3D strategy code closer to orchestration and away from low-level volume manipulation.

## Where the current implementation is still rough

The direction is right, but `volassemble` is starting to carry too much policy detail itself.

Today it decides:

- whether automasking runs
- whether an existing mask is compatible
- how nonuniform filtering gets its mask
- when to reuse or regenerate artifacts

That is execution logic plus policy logic plus artifact-management logic in one commander.

## Policy boundaries

Recommended boundary split:

- `refine3D_strategy`: iteration control, scheduler interaction, state bookkeeping, command-line threading
- `volassemble`: execute volume-domain work for one iteration/state set
- policy/helper layer: decide automask cadence, mask compatibility, nonuniform-mask precedence, artifact naming
- numerical modules: implement mask generation, FSC masking, and nonuniform filtering

## Recommended next step

Keep the execution centralized in `volassemble`, but narrow the commander surface by introducing a dedicated volume postprocessing policy/service.

That helper should return an explicit plan such as:

- generate automask for states `{...}`
- reuse compatible masks for states `{...}`
- run nonuniform filtering with mask source `{state_mask|sphere}`

Then `volassemble` simply executes that plan.

## Practical improvements

- Preserve policy values exactly, including `automsk=tight`.
- Make automask cadence explicit and testable.
- Use one compatibility routine for both mask production and mask consumption.
- Prefer caller-owned state over module-global caches in volume postprocessing code.
- Add end-to-end regression coverage for shared-memory and distributed refine3D using the same expected artifact set.

## Workflow summary

1. refine/matcher code updates orientations and partial reconstructions
2. `volassemble` builds state volumes
3. `volassemble` performs volume-domain postprocessing
4. postprocessed artifacts are written back to the project/output model
5. downstream FSC/reference preparation consumes those artifacts by convention
