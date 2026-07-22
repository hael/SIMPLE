# Flex analysis subsystem

`src/main/flex` contains the workflow-specific implementation behind

```text
simple_exec prg=flex_analysis ...
```

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
