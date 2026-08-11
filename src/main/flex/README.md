# Flex heterogeneity subsystem

`src/main/flex` contains the workflow-specific implementation behind one
heterogeneity program:

```text
simple_exec prg=flex_pca ...        # projection-aware low-rank covariance (PCA)
```

The `flex_analysis` diffusion-map pipeline that used to live here was removed;
its embedding never produced usable states. The shared diffusion-map engines it
used (`../pca/simple_diff_map_graphs.f90`, `../pca/simple_diff_map_denoise.f90`,
`../pca/simple_diffusion_maps.f90`) remain, because `denoise_project`,
`cls_split`, and other applications still depend on them.

## `flex_pca` — projection-aware low-rank covariance

1. `simple_flex_pca_columns.f90` estimates regularized columns of the 3D
   conformational covariance from fixed-pose particles (matched-KB selected-column
   accumulation, regularized halfset merge, reduced projected-covariance solve) and
   derives the low-rank real eigenbasis, then the contrast-aware MAP embedding.
2. `simple_flex_pca_model.f90` is the driver: it owns the embedding cache and its
   resume path, kernel state weights and bandwidth selection, and weighted state
   reconstruction.
3. `simple_flex_pca_rec3D.f90` reconstructs the kernel-weighted states, combined
   and per halfset, in one pass through the gridding reconstructor. Its
   `floor_rho` argument applies a RELION-style shellwise density floor before the
   gridding divide, because kernel weights in `[0,1]` make `rho` small and erratic
   where occupancy is low.
4. `simple_flex_pca_distr.f90` and `simple_flex_pca_parts.f90` own the
   master/worker stage protocol and the part-file I/O.
5. `simple_flex_projected_latent_model.f90` implements the projection-aware
   latent model, and `simple_flex_reconstructor_latent_ops.f90` the flex-specific
   Fourier projection and backprojection operations it builds on.

Alongside the kernel state maps, the state stage writes `outfile` (default
`flex_hard_states.simple`): a copy of the project carrying the hard state label of
every embedded particle in `ptcl3D/state`, with unassigned particles left at
state 0. This lets the embedding and its state assignment be judged with a plain
`simple_exec prg=reconstruct3D projfile=<outfile> nstates=<n>`, independently of
the kernel-weighted backend that shares every upstream assumption with it.

Self-contained tests live in `../../../production/tests/`
(`simple_test_flex_pca.f90`, `simple_test_flex_projected_latent_model.f90`) and
require no data.

The following integration points remain in their architectural layers:

- `../commanders/simple/simple_commanders_flex_pca.f90` is the thin commander
  entrypoint.
- `../exec/simple_exec_denoise.f90`, `../apis/simple_private_exec_api.f90`, and
  `../ui/simple/simple_ui_denoise.f90` register the public and worker command.
- `../volume/simple_reconstructor.f90` remains the shared reconstruction
  implementation used throughout SIMPLE.
