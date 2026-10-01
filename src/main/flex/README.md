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

## Layout

Three subfolders, dependencies pointing downward (`run` -> `states` -> `fit`
-> the leaves); module names are unchanged, the source glob is recursive.

- `run/`: the application and its session (`application`, `run_types`,
  `records`), the distribution contract (`rounds`, `stages`, `artifacts`) and
  the services the application composes (`project_gateway`, `embedding_io`,
  `delivery_3d`, `state_service`).
- `fit/`: the probe fit (`simple_flex_probe_fit` and its four submodules),
  its state (`fit_types`), backends and kernels (`mstep`, `posterior`,
  `pcg`, `polar`, `reconstructor_latent_ops`, `planes`, `plane_cache`,
  `crossfsc`), and the services around it (`basis`, `embed`, `pairmerge`,
  `fit_driver`).
- `states/`: state inference and reconstruction (`weights`, `gmm`,
  `targets`, `deconv`, `merge`, `rec3D`, `states_backend`,
  `states_gridding`, `states_pcg`, `state_delivery`, `state_parts`,
  `weights_state`).
- the folder root: shared helpers (`util`), figures (`plot`, `umap`) and the
  fast-gate testers.

## Ownership

The strategy owns roles, partitions and rounds; each producing module owns its
part I/O; the commander owns the defaults. The architecture and its phases are
recorded in
`doc/refactoring_notes/flex_pca_architecture_audit_and_refactoring_plan_2026_09_18.md`
(there is no `doc/policies/flex_pca_policy.md`).

- `../strategies/parallelization/simple_flex_pca_strategy.f90`: shared-memory,
  distributed-master and worker strategies and the factory (`part=` is a
  worker, `nparts>1` a master). The master carries the qsys context and
  implements the rounds; every strategy hands the application one
  `flex_pca_rounds` and calls `app%run` (master, shared memory) or
  `app%run_worker` (worker).
- `simple_flex_pca_rounds.f90`: the distribution contract the domain receives
  (`flex_pca_rounds`: master/worker/nparts, `plan_partitions`, `run_stage`);
  the stage and fit identifiers and the mod-4 half rule are in
  `simple_flex_pca_stages.f90`.
- `../commanders/simple/simple_commanders_flex_pca.f90`: defaults, the
  factory, the four lifecycle calls.

## `flex_pca` modules

1. `simple_flex_pca_application.f90` is the application (`flex_pca_application`:
   `run` for the shared-memory and master strategies, `run_worker` for a
   worker): prepare -> obtain_embedding -> infer_states -> reconstruct_states
   -> publish over one `flex_run_session`, which owns every product of the
   run. `simple_flex_pca_run_types.f90` holds the run settings
   (`flex_run_settings`, the ONLY reader of the `SIMPLE_COV_*` /
   `SIMPLE_FLEX_*` environment keys, resolved once from `params` and the
   command line) and the session. Services the application composes:
   `simple_flex_pca_project_gateway` (input validation, sigma, the project
   FSC read, weight-store and project writes, out-segment publication),
   `simple_flex_pca_state_service` (latent deconvolution, auto settings,
   state placement with the population floor), `simple_flex_pca_embedding_io`
   (the embedding cache and deconvolution block), `simple_flex_pca_delivery_3d`
   (latent readouts, UMAP figures, covariance tables, eigenvolumes, manifest).
   Their helpers: `simple_flex_pca_util` (chi-squared median, unimodality,
   state-map naming), `simple_flex_pca_gmm` (the tied-covariance GMM and GMM
   AUTO weights), `simple_flex_pca_weights` (kernel/equal-mass/on-axis
   placement, bandwidth selection, half masks), `simple_flex_pca_targets`
   (k-means, diffusion k-centre, FINCH, path and reliability-path targets,
   basis rotations), `simple_flex_pca_deconv` (the latent measurement-error
   mixture).
2. The distribution contract: `simple_flex_pca_rounds.f90` (`flex_pca_rounds`:
   master/worker/nparts, `plan_partitions`, `run_stage(params, request)`),
   `simple_flex_pca_stages.f90` (the stage identifiers, the typed
   `flex_stage_request`, the fit identifiers and the mod-4 half rule) and
   `simple_flex_pca_artifacts.f90` (the part-file magic and part naming).
3. The fit: `simple_flex_pca_fit_types.f90` composes the per-fit state
   (`flex_fit` = spec, model, history, estep, mstep, diag, iter, each with its
   own `kill`) and the typed probe part payload (`flex_probe_part`, borrowed
   from and restored to the fit without copies); `simple_flex_pca_mstep.f90`
   owns the M-step backend (`flex_fit_mstep`: begin_iteration /
   accumulate_batch / reduce_local / apply_ridge / solve_halves /
   snapshot_for_pair_merge / kill_iteration, gridding or PCG storage);
   `simple_flex_pca_posterior.f90` is the posterior inference (SPD solves,
   the ECM contrast alternation, the MCFA mixture).
4. `simple_flex_probe_fit.f90` is the probe fit: `flex_probe_fit` (the fit
   state with its iteration procedures bound) and the fit engine. Like
   `polarft_calc`, one parent declares the type and its bindings and four
   submodules hold the bodies by topic: `_estep` (stage begin, polar bank,
   Cartesian/polar formers, per-particle solve, batch insert, thread reduce,
   the merged-list pass `fit_estep_pass` over one or two fits), `_update`
   (iteration begin, the M-step tail: coupled solve, FSC-Wiener merge,
   re-orthonormalisation, deflation, mixture update, the merge snapshot),
   `_engine` (`fit_engine_iterate`, one master loop for single, paired and
   worker fits; the distributed round; the probe part codec and its
   reduces) and `_crossfsc` (the cross-fit-FSC driver context).
   The services around it are plain modules in dependency order:
   `simple_flex_pca_basis.f90` (the projected mean and its scale, the
   data-free initialiser and its calibration, probe-state and basis I/O,
   pooling, deflation, cross-half angles), `simple_flex_pca_embed.f90` (the
   MAP embedding with per-particle contrast and statistics, basis composition
   across runs), `simple_flex_pca_pairmerge.f90` (the paired final stage:
   frame-align the half fits, merge their raw statistics, one joint solve)
   and `simple_flex_pca_fit_driver.f90` (the paired master run over the
   mod-4 halves and the worker stage bodies). `simple_flex_pca_crossfsc.f90`
   is the cross-fit FSC series and its ridge.
5. State reconstruction: `simple_flex_pca_rec3D.f90` is the service
   (`reconstruct_flex_weighted_states`: selects the backend, runs the state
   stage, hands every state's maps to the delivery); the backends extend
   `flex_states_backend` (`simple_flex_pca_states_backend.f90`: begin /
   accumulate_local_or_write_part / fold_parts / finalize_maps /
   delivery_policy / kill, and the `flex_state_maps` bundle) --
   `simple_flex_pca_states_gridding.f90` (kernel-weighted backprojection of
   every state in one pass; `floor_rho` applies a shellwise density floor
   before the gridding divide because kernel weights in `[0,1]` make `rho`
   small where occupancy is low) and `simple_flex_pca_states_pcg.f90` (one
   cold PCG solve per state and half on the spherical support);
   `simple_flex_pca_state_delivery.f90` is the common delivery (per-state
   eo-FSC, filtering and masking per the backend's declared policy, naming,
   publication, the project update); `simple_flex_pca_state_parts.f90` is the
   parts codec (the per-round weight table, the part names of both backends).
6. `simple_flex_pca_merge.f90` is the two-gate state merge;
   `simple_umap.f90` (`umap=yes`, the default) the UMAP readout.
7. Part-file protocols live next to their producers: probe parts in
   `simple_flex_probe_fit_engine` (the `flex_probe_part` codec), embedding
   statistics in `simple_flex_pca_embed`, the mean scale in
   `simple_flex_pca_basis`, the sigma decision in the project gateway,
   state-weight rounds and state parts in `simple_flex_pca_state_parts`; the
   magic and part naming are in `simple_flex_pca_artifacts.f90`.
8. `simple_flex_weights_state.f90` is the delivered state-weight store:
   the state stage publishes one `flex_weights_state_NNN.bin` per delivered
   state (that state's weight over every physical project row, zero outside
   the selection, a flag marking its hard-labelled particles, and its
   mass/neff/population/bandwidth/target scalars) and registers each in the
   out segment of the run's project copy as imgkind `flex_weights`, state
   `NNN`, beside `vol_flex` state `NNN`, so per-state selection and removal
   apply to weights and maps alike. The files follow the canonical sigma2
   store's rules (layout digest, generation-scoped candidate, atomic
   publish; bytes in `../../fileio/simple_flex_weights_file.f90`) and share
   generation, digest and `nstates` across one delivery. Consumers go
   through `flex_weights_consumable` / `flex_weights_load_state` (one
   state) or `flex_weights_load_all` (the set, cross-checked).
9. `simple_flex_pca_polar.f90` (the polar E-step bank),
   `simple_flex_pca_pcg.f90` (the coupled PCG operator of the M-step) and
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

Self-contained tests live in `simple_flex_pca_tester.f90` and
`simple_flex_pcg_tester.f90` (suites registered in
`../commanders/test/simple_commanders_test_class.f90`) and require no data.

Other integration points: `../exec/simple_exec_denoise.f90`,
`../apis/simple_private_exec_api.f90` and `../ui/simple/simple_ui_heterogeneity.f90`
register the public and worker command; `../volume/simple_reconstructor.f90` is
the shared reconstruction implementation.
