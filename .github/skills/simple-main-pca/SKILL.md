---
name: simple-main-pca
description: Use when working in SIMPLE's src/main/pca subsystem, including SVD PCA, probabilistic PCA, kernel PCA, and the sparse diffusion-map engines (CSR graphs, embeddings, graph denoising) used in denoising and analysis workflows.
---

# SIMPLE `src/main/pca`

This folder owns PCA-family methods.

## Read First

- `simple_pca.f90`: the abstract `pca` type (`new`, `get_feat`, `generate`, `master`, `kill`)
- `simple_pca_svd.f90`
- `simple_ppca.f90`
- `simple_kpca_svd.f90`
- `simple_diffusion_maps.f90`
- `simple_diff_map_graphs.f90`
- `simple_diff_map_denoise.f90`

## Structure

- `pca_svd`, `ppca` and `kpca_svd` extend `pca`; callers pick one by `pca_mode`
- The diffusion-map modules are separate engines (angularly gated kNN graphs in CSR form, embeddings, coefficient-projection denoising), not `pca` extensions

## Connections

- Callers: `commanders/simple/simple_commanders_imgops.f90` (`pca_mode`, including `diffusion_maps` and the PPCA plus kPCA residual hybrid), `simple_commanders_refine2D.f90` (`ppca_denoise_classes`), `motion/simple_denoise_movies.f90`
- Tests: `unit_numerics` ("PCA", "diffusion-map graphs")

## Working Rule

Different PCA variants share the `pca` interface but not their assumptions (noise model, kernel, pre-images); keep variant-specific logic in its own module.
