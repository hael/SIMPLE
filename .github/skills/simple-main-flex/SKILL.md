---
name: simple-main-flex
description: Use when working in SIMPLE's src/main/flex subsystem, including flex_pca heterogeneity analysis (the projection-aware low-rank latent volume model: paired probe EM over mod-4 halves, gridding or PCG M-step, MAP embedding, latent deconvolution, state placement, kernel weights, two-gate merge), its application, rounds, plane cache and project gateway (state weight set, hard ptcl3D labels, delivered state maps), and the PPCA kernels of cls_expansion and trajectory_denoise.
---

# `src/main/flex`: the flex_pca heterogeneity model and its PPCA kernels

Owns the numerics of `flex_pca` (`fit/`: probe EM, M-step, embedding; `states/`: deconvolution,
placement, weights, merge), its application and project services (`run/`), and the arrays-in,
arrays-out class models of `cls_expansion` and `trajectory_denoise`. Not owned: the role, partition
and qsys rounds (`strategies/parallelization/simple_flex_pca_strategy.f90`, passed down as a
`flex_pca_rounds`), the defaults (`commanders/simple/simple_commanders_flex_pca.f90`), map
reconstruction (`simple_rec3D_service`; `reconstruct3D m_estimator=flex` for the delivered maps),
the weight-set format (`project/simple_state_weight_set.f90`) and `cls_expansion`'s I/O (its strategy).

## Rules

Silent ones first: nothing reports these.

- Workers never write the project. A non-worker writes the sigma2 bootstrap and the even/odd
  repair at the start and its products in `publish`: the state weight set, the `flex_embedding`
  entry, hard labels into `ptcl3D` (0 = unassigned), each state's `vol`/`fsc`. The label write
  mutates the live field, so it comes after every stage that reads the selection.
  `refine3D_states` runs flex with `mkdir=no` on its own project and reads the labels and `vol`
  entries back (`simple_refine3D_states_controller.f90`, `run_flex_pca`): change both together. This
  in-place write is the design.
- `ml_reg` is a phase switch: `load_and_validate_sigma` turns it on (the fit whitens its planes by
  `params%l_ml_reg`), `infer_states` turns it off before any trial map.
- The embedding cache holds raw latents: write it before `apply_latent_deconvolution`, or a resume
  (`infile=`) deconvolves twice.
- Under `preimage_auto=yes` the commander leaves `npreimages` undefined: the application reads
  `cline%defined('npreimages')` to raise the ceiling to `AUTO_NSTATES`.
- Every selected particle is in an even/odd half: `validate_covariance_inputs` stops on an `eo = -1` row
  (`oris%get_half`) instead of the fits dropping it.
- `reconstruct_state_halves` keeps the trial maps raw through its `rec3D_request` (`nstates`, `one_part`,
  `no_nonuniform`, `no_trail_rec`, `all_states`, `no_automsk`): the service applies them to `params` for the
  reconstruction and restores the caller's values (and `vols`) after. A new setting the reconstruction must
  take goes into the request, never into a save-and-restore list of the caller. The command line the
  service gives volume assembly (`cline_rec`) still carries the same settings as keys.
- Per iteration (`fit_engine_iterate`): `iter_begin`, E-step pass, `iter_reduce`, `merge_stash`,
  `iter_finish`. The stash precedes `iter_finish`, which solves in place and frees the accumulators.
  `history%mix_*` survive iterations; only `kill` or the rank-change resize in `iter_begin` frees them.
- Threads: per-thread scratch (`fit%iter`, `history`'s `*th` arrays) is sized by `omp_get_max_threads()`
  at `fit_engine_iterate` entry, indexed by `omp_get_thread_num()+1`. Only `estep_particle` runs in
  the E-step region; bank build and `mstep_insert_batch` are called outside it. No `image%new` or
  `read_and_crop` (FFTW plans) in a region. `simple_flex_pca_basis` sets `omp_set_max_active_levels(1)`
  for its pass and restores the caller's level on the way out.
- The mod-4 split (`flex_pca_half_of`: rows {0,1} fit A, {2,3} fit B) is shared by master, workers
  and `check_flex_fits` in `commanders/simple/simple_commanders_refine3D.f90`, which repeats the
  100-row minimum per half (the paired run stops below it, or with `n_probe_iters < 2`).
- A new production `simple_flex_*` module goes into `LAYERS` of `scripts/check_flex_dag.py` (one unit
  per source, no upward edge or cycle, `simple_core_module_api` with `only:`) or the fast gate fails;
  `cls_*`, `traj_denoise` and testers are exempt. `flex_pca` keeps no run state at module scope: the
  application owns the session, plane store and PCG environment.

## Where to start

- `run/simple_flex_pca_application.f90` (phase order, the worker's one stage),
  `run/simple_flex_pca_project_gateway.f90` (every project read and write),
  `fit/simple_flex_probe_fit.f90` (bodies in `_estep`, `_update`, `_engine`, `_crossfsc`).
- `fit/simple_flex_pca_pcg.f90`: `flex_pcg_t` extends `volume/simple_pcg_lattice`, solves with
  `opt/simple_pcg_solver`. `rec_backend` (commander default pcg, `maxits_pcg=4`; `<=0` ships gridding)
  is the basis M-step only; the state maps follow `rec_states_backend` (gridding).
- Particle reads use `strategies/search/simple_matcher_ptcl_io.f90` and `simple_matcher_3Drec.f90`
  (`prep_imgs4rec`), a known deviation (layering policy, section 8).

## Tests

- `build/bin/simple_test_exec test=unit_heterogeneity` (2 s; `suite=` `flex_pca`, `flex_pcg_operator`,
  `flex_cls_expansion`, `flex_cls_expansion_kernel`); `python3 scripts/check_flex_dag.py . --verbose`.
- `test=lib_heterogeneity` (nightly, 35 s; `suite=` `flex_pca_two_state_phantom`,
  `flex_pca_deconvolution_20k`, `flex_pcg_operator_64`, `flex_pcg_solve_sweep`); with a
  `pcg_lattice` or `pcg_solver` change also `unit_reconstruction suite=pcg_lattice`.
- `test=flex_pca_blobs nthr=4` (7 min, 8.4 GiB): shared memory and two workers, label accuracy
  floor 0.60 (profile median 0.99); its `cls_expansion` route is a profile at chance level, not a gate.
  Flex initializing `refine3D_states`: `test=multistate_workflows suite=states`.
- The phantom of `run/simple_flex_pca_application_tester.f90` also feeds `lib_reconstruction`
  (`fractional_reconstruction`, `fraction_dose_export`) and the high-level state fixture.
- `trajectory_denoise`: `single_workflow_fcc` and `_wurtzite` run it, which is its coverage.

## See also

`src/main/flex/README.md` (module map), `doc/algorithms/heterogeneity_analysis/flex_pca.md`;
skills `simple-refine3d`, `simple-main-volume`, `simple-main-opt`.
