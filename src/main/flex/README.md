# Flex heterogeneity subsystem

`src/main/flex` contains the workflow-specific implementation behind one
heterogeneity program:

```text
simple_exec prg=flex_pca ...        # projection-aware low-rank covariance (PCA)
```

The `flex_analysis` diffusion-map pipeline that used to live here was removed;
its embedding never produced usable states, and so was the diffusion-map
`cls_expansion` (2026-10-01: it split on background and density, not structure).
The shared diffusion-map engines (`../pca/simple_diff_map_graphs.f90`,
`../pca/simple_diff_map_denoise.f90`, `../pca/simple_diffusion_maps.f90`)
remain, because `denoise_project` and other applications still depend on them.

## Layout

Three subfolders, dependencies pointing downward or acyclically within one layer
(`run` -> `states` -> `fit` -> the leaves). Every module and submodule has its
own source file; the source glob is recursive.

Shared cross-layer constants live in `src/defs/simple_defs_flex.f90`. The static
`scripts/check_flex_dag.py` check enforces one module/submodule per source and
the production-module layers, and runs before the fast test gate.

- `run/`: the application and its session (`application`, `run_types`,
  `records`), the distribution contract (`rounds`, `stages`, `artifacts`) and
  the services the application composes (`project_gateway`, `embedding_io`,
  `delivery_3d`, `state_service`).
- `fit/`: the probe fit (`simple_flex_probe_fit` and its four submodule files),
  its state (`fit_types`), backends and kernels (`mstep`, `posterior`,
  `pcg` with its self-test submodule `pcg_tester`, `polar`,
  `reconstructor_latent_ops`, `plane_cache`, `planes`, `crossfsc`), and the
  services around it (`basis`, `embed`, `pairmerge`, `fit_driver`).
- `states/`: state inference (`weights`, `gmm`, `targets`, `deconv`, `merge`).
  The maps are SIMPLE's: FLEX has no reconstruction backend (module 5).
- the folder root: shared helpers (`util`), figures (`plot`, `umap`) and the
  fast-gate testers.

There are 41 source files: 32 production modules, the probe fit's four
submodule files, the PCG self-test submodule, three testers and
`cls_expansion` (outside the layering). `scripts/check_flex_dag.py`
rejects multiple module units in one source as well as layering violations and
cycles. CMake discovers `src/main` through a recursive glob without
`CONFIGURE_DEPENDS`; reconfigure once before building after adding or restoring
source files.

## Ownership

The strategy owns roles, partitions and rounds; each producing module owns its
part I/O; the commander owns the defaults. The refactorings that shaped it are
recorded under `doc/refactoring_notes/`: `completed/flex_refactoring_plan_2026_10_06.md`,
the earlier audit beside it, and `completed/flex_on_simple_machinery.md` with its report
`completed/flex_on_simple_machinery_report.md` (FLEX on SIMPLE's reconstruction, weight set and
clustering). There is no `doc/policies/flex_pca_policy.md`.

- `../strategies/parallelization/simple_flex_pca_strategy.f90`: shared-memory,
  distributed-master and worker strategies and the factory (`part=` is a
  worker, `nparts>1` a master). The master carries the qsys context and
  implements the rounds; every strategy hands the application one
  `flex_pca_rounds` and calls `app%run` (master, shared memory) or
  `app%run_worker` (worker).
- `run/simple_flex_pca_rounds.f90`: the distribution contract the domain receives
  (`flex_pca_rounds`: master/worker/nparts, `plan_partitions`, `run_stage`);
  stage/fit identifiers and the mod-4 half rule are in
  `simple_flex_pca_stages.f90`, and the per-run artifact catalog is in
  `simple_flex_pca_artifacts.f90`.
- `../commanders/simple/simple_commanders_flex_pca.f90`: defaults, the
  factory, the four lifecycle calls.

## `flex_pca` modules

1. `simple_flex_pca_application.f90` is the application (`flex_pca_application`:
   `run` for the shared-memory and master strategies, `run_worker` for a
   worker): prepare -> obtain_embedding -> infer_states -> reconstruct_states
   -> publish over one `flex_run_session`, which owns every product of the
   run. `simple_flex_pca_run_types.f90` holds that session; supported run
   configuration is carried by the typed `parameters` object, while the
   application reads command-line definedness only where it changes workflow
   ownership (resume, worker partition, explicit state ceiling). Services the application composes:
   `simple_flex_pca_project_gateway` (input validation, sigma, the project
   FSC read, the state weight set and hard labels, the delivered maps,
   out-segment publication of the maps and the embedding artifact),
   `simple_flex_pca_state_service` (latent deconvolution, auto settings,
   state placement with the population floor, bandwidth-CV trial maps),
   `simple_flex_pca_embedding_io` (the embedding artifact: latents,
   precisions, noise scale and deconvolution block, written atomically), `simple_flex_pca_delivery_3d`
   (latent readouts, UMAP figures, covariance tables, eigenvolumes, manifest).
   Their helpers: `simple_flex_pca_util` (chi-squared median, unimodality,
   state-map naming), `simple_flex_pca_gmm` (state weights from a
   tied-covariance mixture, `simple_gmm`, and the GMM AUTO placement),
   `simple_flex_pca_weights` (kernel/equal-mass/on-axis placement and
   bandwidth-CV numerics, with no project or reconstruction ownership),
   `simple_flex_pca_targets` (k-means targets in the latent metric,
   `simple_kmeans`; the diffusion embedding whose k-center cells,
   `simple_kcenter`, give targets; path and reliability-path targets; basis
   rotations), `simple_flex_pca_deconv` (the noise calibration and the
   latent mixture prior by extreme deconvolution, `simple_xd_gmm`). The
   clustering itself is in `../../utils/clustering`.
2. The distribution contract has one source per module:
   `simple_flex_pca_rounds.f90` (`flex_pca_rounds`: master/worker/nparts,
   `plan_partitions`, `run_stage(params, request)`),
   `simple_flex_pca_stages.f90` (stage identifiers, typed
   `flex_stage_request`, fit identifiers and the mod-4 half rule) and
   `simple_flex_pca_artifacts.f90` (part-file magic and the artifact catalog
   owned by `flex_pca_rounds`). Every codec receives the rounds object; there
   is no process-global part directory.
3. The fit: `simple_flex_pca_fit_types.f90` composes the per-fit state
   (`flex_fit` = spec, model, history, estep, mstep, diag, iter, each with its
   own `kill`) and the typed probe part payload (`flex_probe_part`, borrowed
   from and restored to the fit without copies); `simple_flex_pca_mstep.f90`
   owns the M-step backend (`flex_fit_mstep`: begin_iteration /
   accumulate_batch / reduce_local / apply_ridge / solve_halves /
   snapshot_for_pair_merge / kill_iteration, gridding or PCG storage);
   `simple_flex_pca_posterior.f90` is the posterior inference (SPD solves,
   fixed-contrast latent moments and the MCFA mixture).
4. `simple_flex_probe_fit.f90` is the probe fit: `flex_probe_fit` (the fit
   state with its iteration procedures bound) and the fit engine. Like
   `polarft_calc`, one parent declares the type and its bindings and four
   submodule files hold the bodies by topic: `_estep` (stage begin, polar bank,
   Cartesian/polar formers, per-particle solve, batch insert, thread reduce,
   the merged-list pass `fit_estep_pass` over one or two fits), `_update`
   (iteration begin, the M-step tail: coupled solve, FSC-Wiener merge,
   re-orthonormalisation, deflation, mixture update and the merge snapshot),
   `_engine` (`fit_engine_iterate`, one master loop for single, paired and
   worker fits; the distributed round; the probe part codec and its
   reduces) and `_crossfsc` (the cross-fit-FSC driver context).
   The services around it are plain modules in dependency order:
   `simple_flex_pca_basis.f90` (the projected mean and its scale, the
   data-free initialiser and its calibration, probe-state and basis I/O,
   pooling, deflation, cross-half angles), `simple_flex_pca_embed.f90` (the
   MAP embedding with per-particle contrast and statistics),
   `simple_flex_pca_pairmerge.f90` (the paired final stage:
   frame-align the half fits, merge their raw statistics, one joint solve)
   and `simple_flex_pca_fit_driver.f90` (the paired master run over the
   mod-4 halves and the worker stage bodies). `simple_flex_pca_crossfsc.f90`
   is the cross-fit FSC series and its ridge. Its file, record and driver
   context are lifecycle-owned objects rather than free persistent records.
5. State maps come from SIMPLE's reconstruction service
   (`../strategies/parallelization/simple_rec3D_service.f90`). The state service's
   `reconstruct_state_halves` calls it in the FLEX process with the weight
   columns as a caller-owned table, at the covariance box, raw (shell density
   floor, no filter, mask, prior or trailing chain); the bandwidth
   cross-validation and the merge read these trial halves. The delivered
   weights are published as the project's state weight set
   (`../project/simple_state_weight_set.f90`: one file per state and a manifest
   published last, registered in the out segment as imgkind `state_weights`),
   and the project gateway's `deliver_state_maps` runs `reconstruct3D
   m_estimator=flex` from it at the native box and registers each state's map
   and FSC as ordinary `vol` and `fsc` entries.
6. `simple_flex_pca_merge.f90` is the two-gate state merge (the map gate by
   complete linkage, `simple_hac`); the UMAP readout
   (`umap=yes`, the default) is in `simple_umap.f90`.
7. Part-file protocols live next to their producers: probe parts in
   `simple_flex_probe_fit_engine` (the `flex_probe_part` codec), embedding
   statistics in `simple_flex_pca_embed`, the mean scale in
   `simple_flex_pca_basis`, the sigma decision in the project gateway; the
   magic and part naming are in `simple_flex_pca_artifacts.f90`.
8. The embedding artifact (`simple_flex_pca_embedding_io`) is written to a
   `.part` name and renamed into place, and registered in the out segment of
   the run's project copy as imgkind `flex_embedding`; a resumed run reads its
   noise scale from the artifact.
9. `simple_flex_pca_polar.f90` (the allocatable per-fit polar E-step bank),
   `simple_flex_pca_pcg.f90` (the coupled PCG operator of the M-step, on the
   padded-lattice geometry and support it shares with `simple_reconstructor_pcg`,
   `../volume/simple_pcg_lattice.f90`, and the application-owned, once-resampled
   support environment copied into resident fits) and
   `simple_flex_reconstructor_latent_ops.f90` (the projection-aware latent
   model: the projected mean and basis through `simple_reconstructor`'s
   projection and multi-volume gather, particle preparation through
   `prep_imgs4rec`, the coupled M-step solve; insertion is
   `simple_reconstructor`'s `insert_planes_multi`).
   `simple_flex_pca_plane_cache.f90` owns the disk cache (directory, run token
   and space budget shared with the particle cache; built under a temporary
   name and renamed into place); `simple_flex_pca_planes.f90` owns the
   application-owned `flex_plane_store`, which fetches and retains already
   prepared resident planes and owns one allocatable disk-cache object. The
   cache versions its payload against the absolute project path and
   modification stamp, geometry and master selection. Probe, polish and embed
   workers adopt it only when their partition lies within the completed cache;
   once the master has the embedding, `release` frees the resident planes and
   deletes the cache.
   The store is passed explicitly through particle-read paths and its `kill`
   resets both resident and disk-cache state; neither module retains run state
   at module scope. The application allocates its session, and all stateful
   FLEX objects have explicit lifecycle ownership; the scoped production tree
   retains no mutable run state at module scope.
   Both the coupled rank-4 solve and `simple_reconstructor_pcg` use the rank-1
   recurrence and outcome in `../opt/simple_pcg_solver.f90`; small client adapters
   remap contiguous storage and retain each domain's operator and preconditioner.

Alongside the state weight set, the application writes the hard state label of every
embedded particle into `ptcl3D/state` of the run's own project, leaving
unassigned particles at state 0. `mkdir=yes` already gave the master a private
copy of the project, so this rewrites that copy inside the job directory and
never the project the user pointed at. Only `ptcl3D` is written, as refine3D does
for `nstates>1`. This lets the embedding and its state assignment be judged with a
plain `simple_exec prg=reconstruct3D projfile=<projfile> nstates=<n>`.

Self-contained tests live in `simple_flex_pca_tester.f90`,
`simple_flex_pcg_tester.f90` (asserting the checks of the white-box self-test in
`fit/simple_flex_pca_pcg_tester.f90`) and `run/simple_flex_pca_application_tester.f90`
(suites registered in
`../commanders/test/simple_commanders_test_class.f90`) and require no data.

Other integration points: `../exec/simple_exec_denoise.f90`,
`../apis/simple_private_exec_api.f90` and `../ui/simple/simple_ui_heterogeneity.f90`
register the public and worker command; `../volume/simple_reconstructor.f90` is
the shared reconstruction implementation.

## `cls_expansion`

`simple_flex_cls_expansion.f90` is the per-class covariance model behind `cls_expansion`: a
CTF- and noise-weighted PPCA on the members of one 2D class in the class frame, cross-fitted,
divisive placement into `ncls` subclasses, posterior-precision kernel weights, CTF-corrected
weighted sub-class averages and a cross-half reproducibility per subclass. Arrays in, arrays out.
The strategy (`../strategies/parallelization/simple_cls_expansion_strategy.f90`) prepares the
planes through `transform_ptcls(keep_ft=.true.)`, owns the part files, the merge and the project
writes; the commander (`../commanders/simple/simple_commanders_denoise.f90`) runs the split, one
greedy in-plane `refine2D` round and the restoration. Method and benchmarks:
`doc/implementation_notes/completed/flex_cls_expansion_2026_10_01.md`.
