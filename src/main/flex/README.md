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

## Ownership

The strategy owns roles, partitions and rounds; each producing module owns its
part I/O; the commander owns the defaults. The execution contract is
`doc/policies/flex_pca_policy.md`.

- `../strategies/parallelization/simple_flex_pca_strategy.f90`: shared-memory,
  distributed-master and worker strategies and the factory (`part=` is a
  worker, `nparts>1` a master). The master carries the qsys context and
  implements the rounds; the worker has one entry (`run_flex_pca_worker`).
- `simple_flex_pca_rounds.f90`: the distribution contract the domain receives
  (`flex_pca_rounds`: master/worker/nparts, `plan_partitions`, `run_stage`),
  the stage and fit identifiers, the mod-4 half rule.
- `../commanders/simple/simple_commanders_flex_pca.f90`: defaults, the
  factory, the four lifecycle calls.

## `flex_pca` modules

1. `simple_flex_pca_model.f90` is the driver (`run_flex_pca`, shared by the
   shared-memory and master strategies): validation and particle selection,
   sigma, the mean, the fit (single or paired), embedding, latent
   post-processing (ICA, UMAP readouts), the embedding cache and its resume
   path, outputs. Its helpers are split by ownership: `simple_flex_pca_util`
   (environment switches, chi-squared median, unimodality),
   `simple_flex_pca_gmm` (the tied-covariance GMM and GMM AUTO weights),
   `simple_flex_pca_weights` (kernel/equal-mass/on-axis placement, bandwidth
   selection, half masks), `simple_flex_pca_targets` (k-means, diffusion
   k-centre, FINCH, path and reliability-path targets, basis rotations).
2. `simple_flex_pca_em.f90` is the EM engine's parent (interfaces, the
   `probe_fit_t` per-fit state, environment helpers); its submodules split the
   work: `_env` (environment switches), `_fit` (eigenbasis entry, data-free
   initialiser, probe-state and basis I/O, the worker stage bodies), `_mean`
   (the projected mean and its scale), `_basis` (basis pooling), `_state`
   (`probe_fit_t` lifecycle), `_iter` (the EM iteration: E-step passes,
   the single and paired drivers), `_solve` (the coupled M-step solvers),
   `_embed` (the MAP embedding and its statistics), `_pose` (pose
   refinement, scoring and perturbation), `_polar` (the polar E-step
   accumulation and its ring/band helpers), `_estep` (the shared E-step pass
   over one or two fits, its part I/O), `_mstep` (the per-fit M-step),
   `_crossfsc` (the cross-fit-FSC driver context), `_pairmerge` (the half-set
   accumulator merge). `simple_flex_pca_crossfsc.f90` is the cross-fit FSC
   series and its ridge.
3. `simple_flex_pca_rec3D.f90` reconstructs the weighted states, combined
   and per halfset, in one pass through the gridding reconstructor
   (`floor_rho` applies a shellwise density floor before the gridding divide
   because kernel weights in `[0,1]` make `rho` small where occupancy is low);
   delivery is masked at `mskdiam` under the project-FSC low-pass.
4. `simple_flex_pca_merge.f90` is the two-gate state merge;
   `simple_flex_pca_ica.f90` the latent ICA rotation; `simple_umap.f90` (opt-in, `umap=yes`) the
   UMAP readout.
5. Part-file protocols live next to their producers: probe parts in
   `em_estep`, embedding statistics in `em_embed`, the mean scale in
   `em_mean`, the sigma decision in `model`, state-weight rounds in `rec3D`;
   part naming and the format identifiers are part of the contract in
   `simple_flex_pca_rounds.f90`.
6. `simple_flex_pca_polar.f90` (the polar E-step bank), `simple_flex_gpu.f90`
   (device kernels, `USE_FLEX_CUDA`) and
   `simple_flex_reconstructor_latent_ops.f90` (the projection-aware latent
   model: Fourier projection/backprojection, particle prep, the coupled
   M-step solve).

Alongside the state maps, the state stage writes the hard state label of every
embedded particle into `ptcl3D/state` of the run's own project, leaving
unassigned particles at state 0. `mkdir=yes` already gave the master a private
copy of the project, so this rewrites that copy inside the job directory and
never the project the user pointed at. Only `ptcl3D` is written, as refine3D does
for `nstates>1`. This lets the embedding and its state assignment be judged with a
plain `simple_exec prg=reconstruct3D projfile=<projfile> nstates=<n>`.

Self-contained tests live in `../../../production/tests/`
(`simple_test_flex_pca.f90`) and require no data.

Other integration points: `../exec/simple_exec_denoise.f90`,
`../apis/simple_private_exec_api.f90` and `../ui/simple/simple_ui_denoise.f90`
register the public and worker command; `../volume/simple_reconstructor.f90` is
the shared reconstruction implementation.
