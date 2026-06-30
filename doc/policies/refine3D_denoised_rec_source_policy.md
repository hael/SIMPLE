# Particle Source Policy

This document defines the policy for projects that carry paired raw and
denoised particle representatives.

The current implementation exposes one public source selector:

```text
ptcl_src = raw | den
```

The selected source is used consistently for 3D matching, Euclidean sigma
estimation, and volume reconstruction. SIMPLE does not maintain separate
matching and reconstruction particle sources in this version.

Related policy documents:

- [refine3D_policy.md](refine3D_policy.md)
- [refine3D_auto_policy.md](refine3D_auto_policy.md)
- [refine3D_multi_policy.md](refine3D_multi_policy.md)
- [project_stack_indexing_policy.md](project_stack_indexing_policy.md)
- [sigma_calculation_policy.md](sigma_calculation_policy.md)
- [nonuniform_filtering_policy.md](nonuniform_filtering_policy.md)

## 1. Scientific Objective

The workflow supports projects with two stored image representations per
logical particle:

1. a raw representative image
2. a denoised representative image

At execution time, `ptcl_src` selects the representation used by the 3D
workflow:

```text
ptcl_src=raw  -> raw images supply matching, sigmas, and reconstruction
ptcl_src=den  -> denoised images supply matching, sigmas, and reconstruction
```

There is one Euclidean sigma stream:

```text
sigma2_noise_part*
sigma2_it_*.star
```

No `sigma2_match_*` channel is produced or consumed.

## 2. Required Project Metadata

Dual representative projects store the denoised stack path as metadata on the
same `os_stk` row as the raw stack.

Required fields:

```text
os_stk%stk       raw stack path
os_stk%stk_den   denoised stack path
```

`stk_den` is not derived from `stk` by suffix. Arbitrary denoised stack names
are valid.

The same particle row and physical image index identify the raw and denoised
representatives:

```text
ptcl3D row -> stkind -> os_stk row -> stk     -> raw image path
ptcl3D row -> stkind -> os_stk row -> stk_den -> denoised image path
ptcl3D row -> indstk -> physical image number in either paired stack
```

The denoised stack must have the same physical image order as the raw stack.

## 3. Import Contract

Import commands attach denoised representative paths to `os_stk%stk_den`.

Supported import parameters:

```text
stk_den       denoised stack paired with stk
stktab_den    denoised stack table paired row-by-row with stktab
```

For `stktab_den`, row `i` must be the denoised stack corresponding to raw
`stktab` row `i`.

## 4. Source Resolution

All 3D commands that read particle images through `ptcl_src` must resolve the
stack path through project metadata:

```text
ptcl_src=raw -> os_stk%stk
ptcl_src=den -> os_stk%stk_den
```

`ptcl_src=den` requires `oritype=ptcl3D`; 2D workflows continue to use the
historical raw particle path.

## 5. Command Behavior

`refine3D`, `refine3D_auto`, and `refine3D_multi` pass `ptcl_src` through to
their child 3D matching and sigma-calculation commands.

`calc_pspec` estimates the single sigma stream from the selected particle
source. `calc_group_sigmas` consolidates only `sigma2_noise_part*` into
`sigma2_it_*.star`.

Online Cartesian reconstruction in `refine3D` grids the selected particle
source. Offline `reconstruct3D` also honors `ptcl_src` and reads the selected
source.

## 6. Validation Expectations

Tests should verify:

- `ptcl_src=raw` preserves the historical raw-image behavior
- `ptcl_src=den` uses denoised images for matching
- `ptcl_src=den` estimates sigmas from denoised images
- `ptcl_src=den` reconstructs volumes from denoised images
- distributed worker command lines preserve `ptcl_src`
- `ptcl_src=den` fails clearly when `stk_den` is missing or the stack geometry
  is incompatible
