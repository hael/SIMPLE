# Particle Source Policy

This document defines the policy for projects that carry paired raw and
denoised particle representatives.

The current implementation exposes two public source selectors:

```text
ptcl_src = raw | den
rec_src  = match | raw | den
```

`ptcl_src` selects the particle images used for 3D matching/alignment and
Euclidean sigma estimation. `rec_src` selects the particle images used for 3D
volume reconstruction. The default `rec_src=match` follows `ptcl_src`.

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

At execution time, SIMPLE can select the representation independently for
matching and reconstruction:

```text
ptcl_src=raw -> raw images supply matching and matching sigmas
ptcl_src=den -> denoised images supply matching and matching sigmas

rec_src=match -> reconstruction follows ptcl_src
rec_src=raw   -> raw images supply 3D reconstruction
rec_src=den   -> denoised images supply 3D reconstruction
```

When `rec_src=den`, reconstruction is CC-style and unregularized:

```text
rec_src=den forces ptcl_src=raw
rec_src=den forces ml_reg=no
```

This keeps denoised-image reconstruction out of the Euclidean ML-regularized
path and avoids mixing denoised reconstruction images with denoised matching
residual sigmas.

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

All 3D commands that read particle images through `ptcl_src` or `rec_src` must
resolve stack paths through project metadata:

```text
raw -> os_stk%stk
den -> os_stk%stk_den
```

The `den` source requires `oritype=ptcl3D`; 2D workflows continue to use the
historical raw particle path.

## 5. Command Behavior

`refine3D`, `refine3D_auto`, `refine3D_multi`, and `abinitio3D` expose both
source selectors. Wrapper commands pass both selectors to child 3D matching and
reconstruction commands.

`calc_pspec` estimates the single sigma stream from `ptcl_src`. `calc_group_sigmas`
consolidates only `sigma2_noise_part*` into `sigma2_it_*.star`.

Online Cartesian reconstruction in `refine3D` grids the effective reconstruction
source. When `rec_src=match`, or when explicit `rec_src` equals `ptcl_src`, the
matcher reads each particle batch once and copies the pre-alignment images into
the reconstruction buffer. When the sources differ, the matcher reads the
matching source into the polar-alignment buffer and reads the reconstruction
source directly into the reconstruction buffer.

Offline `reconstruct3D` honors `rec_src` and reads that source directly.

## 6. Validation Expectations

Tests should verify:

- `ptcl_src=raw rec_src=match` preserves historical raw-image behavior
- `ptcl_src=den rec_src=match` uses denoised images for matching and
  reconstruction
- `ptcl_src=den rec_src=raw` uses denoised images for matching and raw images
  for reconstruction
- `ptcl_src=raw rec_src=den ml_reg=no` uses raw images for matching and
  denoised images for reconstruction
- `rec_src=den ptcl_src=den` resolves to raw matching and denoised reconstruction
- `rec_src=den ml_reg=yes` resolves to unregularized reconstruction
- distributed worker command lines preserve both source selectors
- `den` source selection fails clearly when `stk_den` is missing or the stack
  geometry is incompatible
