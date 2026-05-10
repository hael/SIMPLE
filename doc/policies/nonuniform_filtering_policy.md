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

Nonuniform filtering is executed in Cartesian `volassemble`, after:

- even/odd assembly
- sampling-density correction and FSC estimation
- merged-volume restoration and write
- gridding correction
- optional trailing half-map blending
- per-state postprocessing planning and any required automask regeneration

This is intentional. The expensive shared-memory volume work happens in one place for both shared-memory and distributed `refine3D` paths, with `volassemble` owning the derived filtered references.

The filter must not consume even/odd maps after low-resolution insertion. Low-resolution even/odd blending is a registration-reference trick applied only after reference read/mask/filter and immediately before generating the reprojection model for 3D registration. `volassemble` must not write low-resolution-blended half-maps to disk, and those blended maps must not feed FSC, automasking, nonuniform filtering, or ordinary half-map handoffs.

If trailing reconstruction is active, the trailing blend is applied to the clean restored half-maps before automask generation and before `_nu_filt` products are generated. Because low-resolution-blended half-maps are never persisted, the ordinary previous even/odd half-maps remain valid trailing inputs.

## Inputs

The nonuniform filter consumes:

- the current unfiltered even volume
- the current unfiltered odd volume
- optionally, auxiliary pre-filtered even/odd candidate pairs supplied by `volassemble`
- an auxiliary effective resolution, in Angstrom, for each auxiliary candidate pair
- a logical support mask

Mask precedence is:

1. freshly regenerated state-specific automask from `volassemble`
2. existing compatible state-specific automask, when `automsk != 'no'`
3. spherical support mask derived from `mskdiam`

When `ml_reg=yes`, `volassemble` uses the `_unfil` even/odd pair as the base nonuniform input and supplies the ML-regularized even/odd pair as an auxiliary candidate source. Both pairs are clean of low-resolution insertion and receive any trailing blend before filtering. The auxiliary coordinate is derived from the state FSC(0.143) resolution, `res0143s(state)`.

## Outputs

For each state, the filter produces:

- `vol_state_even_nu_filt.mrc`
- `vol_state_odd_nu_filt.mrc`
- `vol_state_nu_filt.mrc`

Actual filenames are built by appending `NUFILT_SUFFIX`, currently `_nu_filt`, to the base even, odd, and merged state volume names. These are derived products. The base even/odd and merged volumes remain the primary reconstruction outputs.

## Current implementation strategy

The filter currently:

1. builds a bank of low-pass filtered even/odd volumes from the unfiltered pair
2. optionally appends auxiliary pre-filtered even/odd pairs to that candidate bank
3. maps all candidates onto a shared filter-bank coordinate axis
4. caches the low-pass-filtered base bank on disk
5. computes voxelwise objective maps across all candidates
6. chooses the best candidate per voxel
7. optionally refines the candidate map with an ordered-label spatial smoothing prior
8. logs selected-low-pass statistics and neighbor-continuity diagnostics
9. synthesizes filtered even/odd outputs from the selected candidate map
10. writes the merged `_nu_filt` volume as the average of the filtered even/odd outputs

This design is scientifically reasonable and easy to debug, but it pays a large I/O and memory-traffic cost.

In nonuniform mode, matcher reference loading tries `_nu_filt` even/odd references first, falls back to the regular even/odd references before filtered products exist, and avoids applying the ordinary low-pass filter on top of the nonuniform reference path.

## Optional ordered-label smoothing regularization

The default nonuniform filter is still the unary voxelwise selector described above.
For testing, the implementation can optionally apply an ordered-label smoothing refinement to the candidate-label map before writing `_nu_filt` volumes.

The smoothing stage is intended to reduce abrupt local jumps in the selected filter-bank label. It initializes from the ordinary voxelwise argmin, then runs a small number of ICM-style passes over the label map using the fully connected 26-neighbor 3D voxel neighborhood. Updates use an 8-color parity schedule so voxels updated within the same pass are not neighbors under the full 3x3x3 neighborhood. The neighborhood penalty is evaluated on a candidate-coordinate axis rather than on raw label indices:

- base low-pass candidates use coordinates `1..n_base`
- auxiliary candidates are projected onto that same axis from their required effective resolutions
- one-step differences are tolerated by the current penalty setting
- jumps larger than one filter-bank step are penalized with a convex hinge to discourage larger jumps more strongly
- neighbor penalties are normalized by the number of in-mask neighbors, so boundary and thin-mask voxels do not receive systematically weaker or stronger regularization
- ties preserve the current label within a small tolerance instead of drifting to the lowest candidate index

Auxiliary candidate resolutions are mandatory whenever auxiliary candidate volumes are supplied. In the current `volassemble` ML-regularized path, the auxiliary candidate is the ML-regularized even/odd pair and its coordinate is derived from the state FSC(0.143) resolution, `res0143s(state)`.

Diagnostics log the estimated smoothing beta, candidate and auxiliary counts, jump-penalty settings, candidate coordinates, changed voxels per iteration, and mean site energy.

### Activating ordered-label smoothing

For standalone testing, `nu_filt3D` exposes the smoothing stage through a typed command-line parameter:

```bash
simple_exec prg=nu_filt3D vol1=odd.mrc vol2=even.mrc smpd=1.0 potts_prior=yes
```

The public `potts_prior` name is retained for familiarity and for compatibility with the current `optimize_nu_cutoff_finds(l_potts_prior=...)` keyword, but the implementation is the ordered-label smoothing prior described above.

For iterative `refine3D` testing, enable the ordered-label smoothing stage in one of two code-level ways:

1. Prefer the local call-site override in `src/main/commanders/simple/simple_commanders_rec_distr.f90`:

   ```fortran
   call optimize_nu_cutoff_finds(l_potts_prior=.true.)
   ```

   The keyword name is retained for source compatibility. This activates ordered-label smoothing only for the Cartesian `volassemble` path used by iterative 3D refinement.

2. For broader developer testing, temporarily set the module default in `src/utils/filter/simple_nu_filter.f90`:

   ```fortran
   logical, parameter :: L_NU_LABEL_SMOOTH_DEFAULT = .true.
   ```

   This affects every caller that does not pass an explicit `l_potts_prior` argument.

Keep the default off until the regularization strength, convergence behavior, and reference-map continuity metrics have been evaluated on representative datasets.

## Strengths of the current design

- Shared-memory volume work is centralized in one execution step.
- The mask used for automasking can also guide nonuniform filtering.
- The implementation is explicit and testable.
- Cached filtered volumes reduce repeated FFT/filter work when the same cutoff bank is reused.
- The optional ordered-label smoothing stage is orthogonal to candidate-bank construction and can be enabled without changing output synthesis.

## Current limitations

- Disk-backed cache files add avoidable filesystem traffic inside an already expensive volume step.
- The cutoff search still performs repeated full-volume passes.
- Ordered-label smoothing regularization is exposed for standalone `nu_filt3D` testing, but not yet as a public iterative refinement workflow parameter.
- Policy and execution are mixed together in `volassemble`, which makes the commander harder to reason about.
- The nonuniform filter depends on module-level cached state, which limits composability and future concurrency.

## Recommended direction

Medium-term improvements:

- Option to dynamically change the filter bank?
- Return both the selected cutoff map and summary statistics as explicit outputs.
- Promote ordered-label smoothing activation and strength to explicit workflow parameters if validation supports it.

Architectural target:

- `volassemble` should orchestrate nonuniform filtering.
- A dedicated volume postprocessing service should own the filter-bank setup, objective evaluation, and output synthesis.
- When `ml_reg=yes`, `volassemble` should feed the `_unfil` pair as the base nonuniform input and may add the ML-regularized pair as an auxiliary candidate source.
