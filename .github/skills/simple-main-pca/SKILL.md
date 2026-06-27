---
name: simple-main-pca
description: Use when working in SIMPLE's src/main/pca subsystem, including PCA, PPCA, complex PPCA, kernel PCA, diffusion maps, and SVD-backed dimensionality-reduction code used in analysis and denoising workflows.
---

# SIMPLE `src/main/pca`

This folder owns PCA-family methods.

## Read First

- `simple_pca.f90`
- `simple_pca_svd.f90`
- `simple_ppca.f90`
- `simple_complex_ppca.f90`
- `simple_kpca_svd.f90`
- `simple_diffusion_maps.f90`
- `simple_diff_map_denoise.f90`

## Connections

- Often paired with image stacks, denoising, and analysis workflows
- Likely depends on `image/`, `project/`, and numeric utilities

## Working Rule

Different PCA variants may share concepts but not assumptions; avoid forcing a single abstraction unless the code already wants it.
