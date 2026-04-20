# Nonuniform Filtering Policy

## Scope

This document defines how nonuniform filtering is selected, where it runs, what mask it uses, and how it should evolve inside the current `refine3D` and `abinitio3D` architecture.

## Public policy

The user-facing control is `filt_mode`.

Supported values:

- `none`
- `uniform`
- `fsc`
- `nonuniform`

Only `filt_mode=nonuniform` activates the nonuniform volume filter path.

## Execution point

Nonuniform filtering is executed in `volassemble`, after:

- even/odd assembly
- gridding correction
- eo low-resolution reinsertion

This is intentional. The expensive shared-memory volume work now happens in one place for both shared-memory and distributed refine3D paths.

## Inputs

The nonuniform filter consumes:

- the current unfiltered even volume
- the current unfiltered odd volume
- optionally, auxiliary pre-filtered even/odd candidate pairs supplied by `volassemble`
- a logical support mask

Mask precedence is:

1. state-specific automask from `volassemble`, when `automsk != 'no'` and a compatible mask is available
2. spherical support mask derived from `mskdiam`

## Outputs

For each state, the filter produces:

- `vol_state_even_nu.mrc`
- `vol_state_odd_nu.mrc`
- `vol_state_nu.mrc`

These are derived products. The base even/odd and merged volumes remain the primary reconstruction outputs.

## Current implementation strategy

The filter currently:

1. builds a bank of low-pass filtered even/odd volumes from the unfiltered pair
2. optionally appends auxiliary pre-filtered even/odd pairs to that candidate bank
3. caches the low-pass-filtered base bank on disk
4. computes voxelwise objective maps across all candidates
5. chooses the best candidate per voxel
6. synthesizes filtered even/odd outputs from the selected candidate map

This design is scientifically reasonable and easy to debug, but it pays a large I/O and memory-traffic cost.

## Strengths of the current design

- Shared-memory volume work is centralized in one execution step.
- The mask used for automasking can also guide nonuniform filtering.
- The implementation is explicit and testable.
- Cached filtered volumes reduce repeated FFT/filter work when the same cutoff bank is reused.

## Current limitations

- Disk-backed cache files add avoidable filesystem traffic inside an already expensive volume step.
- The cutoff search still performs repeated full-volume passes.
- Policy and execution are mixed together in `volassemble`, which makes the commander harder to reason about.
- The nonuniform filter depends on module-level cached state, which limits composability and future concurrency.

## Recommended direction

Medium-term improvements:

- Option to dynamically change the filter bank?
- Return both the selected cutoff map and summary statistics as explicit outputs.

Architectural target:

- `volassemble` should orchestrate nonuniform filtering.
- A dedicated volume postprocessing service should own the filter-bank setup, objective evaluation, and output synthesis.
- When `ml_reg=yes`, `volassemble` should feed the `_unfil` pair as the base nonuniform input and may add the ML-regularized pair as an auxiliary candidate source.
