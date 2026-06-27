# Nonuniform Regularization Map

## Public Policy

The user-facing selector is `filt_mode`:

- `none`
- `uniform`
- `fsc`
- `nonuniform`
- `nonuniform_lpset`

`filt_mode=nonuniform` activates NU-filtered volume products while preserving
ordinary half-map matching unless another policy changes the bandwidth.
`filt_mode=nonuniform_lpset` activates the same NU products and promotes the
selected NU bandwidth into LP-set matching with merged registration-reference
topology. The public policy is documented in
`doc/policies/nonuniform_filtering_policy.md`.

`src/main/params/simple_parameters.f90` declares the value and
`simple_parameters_phases.f90` derives `l_nonuniform` and
`l_nonuniform_lpset`. Unlike `uniform` and `fsc`, plain `nonuniform` does not
set `l_lpauto`/`l_lpset`; filtering is handled when Cartesian state volumes are
assembled.

For staged `abinitio3D`, `simple_abinitio_controller.f90` enables NU filtering
only after the NU stage boundary, emits `nu_refine=no`, promotes staged
`nonuniform` to `nonuniform_lpset` before the disabled gold-standard stage, and
does not enable NU on the class-average route.

## Execution Path

Main owner:

- `src/main/commanders/simple/simple_commanders_rec_distr.f90`
- `commander_cartesian_volassemble%execute`

High-level order for each state:

1. Reduce partition-local even/odd reconstruction sums.
2. Sampling-density correct even/odd volumes.
3. Reinsert low-resolution merged content into the even/odd pair.
4. Apply gridding correction.
5. Write base even, odd, and merged state volumes.
6. Plan automask and nonuniform-mask behavior with `plan_state_postprocess`.
7. Generate or reuse the state automask if requested.
8. Build a logical nonuniform support mask from automask or spherical `mskdiam`.
9. Run `simple_nu_filter`.
10. Write `_nu_filt` even, odd, merged, and `_nu_locres` derived products.
11. Record the selected NU matching low-pass and refresh polar reference sections for the next matcher pass.

## Mask Policy

Policy owner:

- `src/main/volume/simple_vol_pproc_policy.f90`

`plan_state_postprocess(params, state, which_iter, l_nonuniform_mode, plan)` decides:

- whether automasking is enabled
- whether an existing state mask exists
- whether that mask is compatible with current `box_crop` and `smpd_crop`
- whether to regenerate or reuse automask
- which mask source nonuniform filtering should use

Nonuniform mask source precedence:

1. `NU_MASK_SOURCE_FRESH_AUTOMASK`
2. `NU_MASK_SOURCE_EXISTING_AUTOMASK`
3. `NU_MASK_SOURCE_SPHERICAL`

Automask file convention:

- `automask3D_stateNN.mrc`

Compatibility requires all three dimensions to match `box_crop` and sampling
distance to match `smpd_crop` within tolerance.

## Numerical Filter Module

Implementation owner:

- `src/main/nu_filt/simple_nu_filter.f90`
- `src/main/nu_filt/simple_nu_filter_*.f90`

Required lifecycle:

```fortran
call setup_nu_dmats(vol_even, vol_odd, l_mask, aux_resolutions, aux_even, aux_odd)
call optimize_nu_cutoff_finds()
call nu_filter_vols(vol_even_nu, vol_odd_nu)
call cleanup_nu_filter()
```

Current base low-pass bank:

- `20., 15., 12., 10., 8., 6., 5., 4.` Angstrom

Optional high-resolution extension support exists through:

- `extend_nu_filter_highres_shell_next`
- `extend_nu_filter_highres_shells`

As of this map, the normal `volassemble` path calls setup, optimization,
filtering, stats/continuity logging, output write, and cleanup. It calls
high-resolution extension only when `nu_refine=yes`.

## Candidate Selection Model

`setup_nu_dmats`:

- validates input even/odd dimensions and sampling
- validates mask shape and nonempty support
- builds Butterworth filters for the cutoff bank
- caches filtered base even/odd pairs on disk
- computes mask-packed objective costs for each candidate
- smooths each candidate objective over mask-normalized local support
- uses auxiliary replacement volumes only when their effective resolution is
  finer than the finest retained static-bank label

`optimize_nu_cutoff_finds`:

- picks the lowest objective candidate per voxel inside the mask
- runs ordered-label Potts smoothing before finalizing the cutoff map, except
  for degenerate exits such as a single label or numerical-zero prior strength
- stores base-bank winners in `filtmap`
- caches the finest base objective map for possible high-resolution extension
- releases the large objective storage

The ordered-label prior uses retained-bank coordinates rather than raw label
numbers. Adjacent retained-bank coordinate steps are tolerated; larger jumps are
penalized with a linear-quadratic hinge.

`nu_filter_vols`:

- reads filtered base candidate volumes from cache
- copies candidate voxels into output even/odd volumes according to `filtmap`
- copies auxiliary replacement voxels where the replacement label is selected

## ML-Regularized Coupling

When `params%l_ml_reg` is true and `params%l_nu_refine` is false,
`volassemble` uses:

- base input: `_unfil` even/odd pair
- auxiliary replacement: the ML-regularized even/odd pair currently held in
  `build%vol` and `build%vol2`

This lets voxelwise optimization replace the finest static-bank member with the
ML-regularized pair when it represents a finer effective resolution. When
`nu_refine=yes`, the auxiliary pair is not supplied because shell extension owns
the resolution-extension experiment.

## Outputs and Consumers

Output suffixes:

- `NUFILT_SUFFIX = '_nu_filt'`
- `NULOCRES_SUFFIX = '_nu_locres'`

For each state, `volassemble` writes:

- `vol_state_even_nu_filt.mrc`
- `vol_state_odd_nu_filt.mrc`
- `vol_state_nu_filt.mrc`
- `vol_state_nu_locres.mrc`

Filtered names are built via `add2fbody(..., params%ext, NUFILT_SUFFIX)`.
The local-resolution map uses `NULOCRES_SUFFIX`.

Matcher reference loading is in:

- `src/main/strategies/search/simple_matcher_refvol_utils.f90`

Plain `nonuniform` tries `_nu_filt` even/odd references first, then falls back
to regular even/odd files. `nonuniform_lpset` with active LP-set matching uses
merged registration-reference topology and prefers the merged `_nu_filt`
product. Fallback matters before the first assembled nonuniform products exist.

## Tests and Manual Checks

Standalone test driver:

- `production/tests/simple_test_nu_filter.f90`

Usage shape:

```text
simple_test_nu_filter even.mrc odd.mrc smpd mskdiam [out_even.mrc out_odd.mrc] [aux_even.mrc aux_odd.mrc aux_res]
```

For code changes, consider:

- focused build/test of `simple_test_nu_filter`
- checking cleanup removes `nu_filter_cache_even_k_*.mrc` and
  `nu_filter_cache_odd_k_*.mrc`
- checking nonuniform mode still falls back to regular references before
  `_nu_filt` files exist
- checking multi-state automask naming and compatibility behavior
- checking `ml_reg=yes` static-NU auxiliary replacement
- checking `nu_refine=yes` shell extension and `nu_highres_depth_stateNN.txt`
  persistence

## Current Limitations

- The base candidate bank is disk-backed, which creates avoidable filesystem traffic.
- The module uses module-level state, which makes concurrent use awkward.
- `volassemble` still contains substantial orchestration around mask creation,
  auxiliary setup, logging, output naming, and cleanup.

Preferred architectural direction:

- Keep `volassemble` as the orchestrator.
- Move filter-bank setup, objective evaluation, output synthesis, and summaries
  behind a dedicated volume postprocessing service when the next refactor is justified.
- Return selected-cutoff maps and summary statistics as explicit outputs rather
  than relying only on module-level state and logging.
