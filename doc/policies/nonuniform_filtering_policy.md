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
- an auxiliary effective resolution, in Angstrom, for each auxiliary candidate pair
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
3. maps all candidates onto a shared filter-bank coordinate axis
4. caches the low-pass-filtered base bank on disk
5. computes voxelwise objective maps across all candidates
6. chooses the best candidate per voxel
7. optionally refines the candidate map with a Potts/ICM spatial prior
8. synthesizes filtered even/odd outputs from the selected candidate map

This design is scientifically reasonable and easy to debug, but it pays a large I/O and memory-traffic cost.

## Optional Potts/ICM label regularization

The default nonuniform filter is still the unary voxelwise selector described above.
For testing, the implementation can optionally apply a Potts/ICM refinement to the candidate-label map before writing `_nu_filt` volumes.

The Potts/ICM stage is intended to reduce abrupt local jumps in the selected filter-bank label. It initializes from the ordinary voxelwise argmin, then runs a small number of ICM passes over the label map using the library 6-connected 3D pixel neighborhood (`neigh_4_3D`). The neighborhood penalty is evaluated on a candidate-coordinate axis rather than on raw label indices:

- base low-pass candidates use coordinates `1..n_base`
- auxiliary candidates are projected onto that same axis from their required effective resolutions
- one-step differences are tolerated by the current penalty setting
- jumps larger than one filter-bank step are penalized

Auxiliary candidate resolutions are mandatory whenever auxiliary candidate volumes are supplied. In the current `volassemble` ML-regularized path, the auxiliary candidate is the ML-regularized even/odd pair and its coordinate is derived from the state FSC(0.143) resolution, `res0143s(state)`.

### Activating Potts/ICM for iterative 3D refinement tests

There is no public command-line switch yet. For iterative `refine3D` testing, enable the Potts/ICM stage in one of two code-level ways:

1. Prefer the local call-site override in `src/main/commanders/simple/simple_commanders_rec_distr.f90`:

   ```fortran
   call optimize_nu_cutoff_finds(l_potts_prior=.true.)
   ```

   This activates Potts/ICM only for the Cartesian `volassemble` path used by iterative 3D refinement.

2. For broader developer testing, temporarily set the module default in `src/utils/filter/simple_nu_filter.f90`:

   ```fortran
   logical, parameter :: L_NU_POTTS_PRIOR_DEFAULT = .true.
   ```

   This affects every caller that does not pass an explicit `l_potts_prior` argument.

Keep the public default off until the regularization strength, convergence behavior, and reference-map continuity metrics have been evaluated on representative datasets.

## Strengths of the current design

- Shared-memory volume work is centralized in one execution step.
- The mask used for automasking can also guide nonuniform filtering.
- The implementation is explicit and testable.
- Cached filtered volumes reduce repeated FFT/filter work when the same cutoff bank is reused.
- The optional Potts/ICM stage is orthogonal to candidate-bank construction and can be enabled without changing output synthesis.

## Current limitations

- Disk-backed cache files add avoidable filesystem traffic inside an already expensive volume step.
- The cutoff search still performs repeated full-volume passes.
- Potts/ICM regularization currently has code-level activation only; it is not yet exposed as a public workflow parameter.
- Policy and execution are mixed together in `volassemble`, which makes the commander harder to reason about.
- The nonuniform filter depends on module-level cached state, which limits composability and future concurrency.

## Recommended direction

Medium-term improvements:

- Option to dynamically change the filter bank?
- Return both the selected cutoff map and summary statistics as explicit outputs.
- Promote Potts/ICM activation and strength to explicit workflow parameters if validation supports it.

Architectural target:

- `volassemble` should orchestrate nonuniform filtering.
- A dedicated volume postprocessing service should own the filter-bank setup, objective evaluation, and output synthesis.
- When `ml_reg=yes`, `volassemble` should feed the `_unfil` pair as the base nonuniform input and may add the ML-regularized pair as an auxiliary candidate source.
