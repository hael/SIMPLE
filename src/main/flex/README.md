# Flex analysis subsystem

`src/main/flex` contains the workflow-specific implementation behind two
heterogeneity programs:

```text
simple_exec prg=flex_analysis ...   # diffusion-map pipeline
simple_exec prg=flex_pca ...        # projection-aware low-rank covariance (PCA)
```

They are independent pipelines that share the reconstruction and projected-latent
layers. `flex_analysis` is also consumed by `nano3D` and the trajectory chunker;
`flex_pca` is standalone.

## `flex_analysis` — diffusion map

The directory is organized by the data flow through the analysis:

1. `simple_flex_analysis_strategy.f90` owns the shared-memory, distributed
   master, and distributed worker strategies. It also coordinates registration,
   graph construction, embedding, pre-image selection, and reconstruction.
2. `simple_flex_diffmap_features.f90` prepares registered residual particle
   images and their diffusion-map feature vectors.
3. `simple_flex_diffmap_preimage.f90` selects representative latent states.
4. `simple_flex_diffmap_rec3D.f90` reconstructs the representative 3D
   pre-images and handles distributed reconstruction parts.
5. `simple_flex_projected_latent_model.f90` implements the projection-aware
   latent model used by the reconstruction update.
6. `simple_flex_reconstructor_latent_ops.f90` contains the flex-specific
   Fourier projection and backprojection operations used by that model.

`flex_embedding_result`, the compact embedding returned to `nano3D`, is part
of `simple_flex_analysis_strategy.f90`: it is the public result of that
strategy rather than a module of its own.

## `flex_pca` — projection-aware low-rank covariance

7. `simple_flex_pca_columns.f90` estimates regularized columns of the 3D
   conformational covariance from fixed-pose particles (matched-KB selected-column
   accumulation, regularized halfset merge, reduced projected-covariance solve) and
   derives the low-rank real eigenbasis, then the contrast-aware MAP embedding.
8. `simple_flex_pca_model.f90` is the driver: it owns the embedding cache and its
   resume path, kernel state weights and bandwidth selection, weighted state
   reconstruction, and the volume-space trajectory ordering.

`flex_pca` reuses `simple_flex_diffmap_rec3D.f90` for weighted state
reconstruction. That routine takes an optional `floor_rho` argument, off by
default, which applies a RELION-style shellwise density floor before the
gridding divide; only `flex_pca` opts in, because kernel weights in `[0,1]`
make `rho` small and erratic where occupancy is low.

Self-contained tests for both pipelines live in `../../../production/tests/`
(`simple_test_flex_pca.f90`, `simple_test_flex_projected_latent_model.f90`,
`simple_test_flex_diffmap_graph.f90`) and require no data.

The following integration points remain in their architectural layers:

- `../commanders/simple/simple_commanders_flex_analysis.f90` is the thin
  commander entrypoint.
- `../exec/simple_exec_denoise.f90`, `../apis/simple_private_exec_api.f90`, and
  `../ui/simple/simple_ui_denoise.f90` register the public and worker command.
- `../pca/simple_diff_map_graphs.f90` and
  `../pca/simple_diff_map_denoise.f90` remain shared diffusion-map engines; they
  are also used by `denoise_project`, `cls_split`, and other applications.
- `../volume/simple_reconstructor.f90` remains the shared reconstruction
  implementation used throughout SIMPLE.

The numerical contract and distributed-execution policy are documented in
`../../../doc/policies/flex_analysis_policy.md`. Focused tests live in
`../../../production/tests/simple_test_flex_projected_latent_model.f90` and
`../../../production/tests/simple_test_flex_diffmap_graph.f90`.

For isolating embedding/clustering from the pre-image estimator, run with
`preimage_mode=discrete`. This writes hard k-medoids labels into
`ptcl3D/state` in `outfile` (default `flex_cluster_states.simple`) and skips
pre-image reconstruction, so the clusters can be checked with ordinary
`reconstruct3D`.
