---
name: simple-main-strategies-search
description: Use when working in SIMPLE's src/main/strategies/search: the refine3D_exec and refine2D_exec matchers (a worker's or a shared-memory pass), particle batches (build_batch_particles2D/3D, read sessions, the particle and record buffers), the probabilistic batch pipeline the matchers run (per-batch fill and keyed draw, the part's readout record), reference preparation (reprojection model, 2D references, the band), a pass's sigma2 preparation, the per-particle 2D/3D strategies (greedy, shc, snhc, inpl, neigh, prob, cont, tseries) and their per-thread state, sampling and low-pass helpers, the classic polar shift map, 3D partial reconstruction, class-average registration utilities.
---

# `strategies/search`: the matchers and the per-particle searches

Owns one pass over a partition's particles: `refine3D_exec` (`simple_strategy3D_matcher.f90`) and
`refine2D_exec` (`simple_strategy2D_matcher.f90`) sample, prepare the references and particle batches, run
one `strategy3D`/`strategy2D` object per particle and write the part's outputs; the 3D pass then writes its
partial reconstructions, the 2D pass its class sums. Under `refine=prob*` it also fills each batch's table
rows and draws each particle's pose or class from its row before the loop. It does not draw a distributed run's
sample or merge (`parallelization/`), score a pair (`pftc/`), restore class averages (`class/`) or commit sigma2
(`sigma2/`). The tables themselves are `probabilistic/` (`simple-main-strategies-probabilistic`); the matchers
drive them per batch.

## Rules

Silent ones first.

- 3D pass: `ensure_field_halves` -> `adopt_reprojection_model_range` (Cartesian: `read_cart_refvols`) ->
  sample -> `read_reprojection_model` (`pftc%new`) -> `prep_sigmas_objfun` -> `alloc_ptcl_imgs` ->
  `memoize_refs` -> `prep_strategy3D` -> `begin_pipeline` (prob) -> `fill_pass_buffer` or
  `open_record_buffer` -> per batch `build_batch_particles3D`, `pipeline_batch` (prob), the loop ->
  `close_assignment` (prob) -> `write_sigma2`, orientations -> `pftc%kill` (stops if still allocated) ->
  partial reconstructions -> `release_pass_buffer` -> `qsys_declare_part_finished`. 2D: `cavger_new`,
  `cavger_read_all`, `cavger_init_online` -> `set_bp_range2D` -> `prep_pftc4align2D` (`pftc%new`,
  `prep_sigmas_objfun`, references, centring offsets, `memoize_refs`) -> `begin_pipeline` (prob) -> per
  batch `build_batch_particles2D`, `pipeline_batch` (prob), `prep_strategy2D_batch`, the loop,
  `cavger_update_sums` -> `close_assignment` (prob) -> sigma2 range, orientations, a worker's class-sum
  contribution.
- `pipeline_batch` fills the rows of the batch just prepared (`fill_tab_range`) and draws from them
  (`assign_batch`) before the batch's loop; the strategy objects read the drawn pose. `begin_pipeline` sets
  the iteration's model from the complete incoming field (`set_model_from_field`) before any batch:
  nothing a batch computes may come from the field it is changing.
- `prep_sigmas_objfun` after every `pftc%new`: the calculator points at `build%esig`'s spectra and `new`
  drops the pointer; without them `memoize_sigma2_weights` returns with no weights and a Euclidean
  score stops (`check_euclid_weights`).
- A worker never makes its 3D band: `adopt_reprojection_model_range` takes `kfromto` and `lp` from the
  reprojection model's header (Cartesian: `read_cart_refvols`) and stops on an incompatible one.
- In the particle loop a contained procedure sees the host's variables, not the loop's `firstprivate`
  copies: declare per-particle objects in the procedure (`o_sigma` in `choose_and_run_strategy`). Per-thread
  state is `s3D`/`s2D` columns `ithr`. The classic modes' draws in the loop (`ran3`, the `s3D%rts*` orders,
  `shcloc`, `irnd_uni`) use the thread's stream, so under the dynamic schedule a particle's draws depend on the
  thread. Every `prob*` draw (sample, in-plane, candidate order, assignment, 2D start class) is keyed on
  (seed, pass, iteration, particle, draw) (`rnd_draw_key`, `prob_rnd_scope`) and depends on no thread, batch or
  part; the project owner fixes the seed first (`adopt_run_seed`). A new `prob*` draw takes a key, not `ran3`.
- A worker reproduces the sample (`sample4update_reprod`) whenever the master drew it (probabilistic
  modes, `master_draws_sample*`, the polish) and draws only otherwise.
- Preparation records each particle's mask (`note_obs_mask`); class sums and 3D reconstruction normalise
  under the same one (`obs_on_box`). A cropped 3D record (`put_record` in `polarize_batch_particles3D`) is
  the reconstruction's observation bit for bit (`prep_rec_observation`): change both together.
  The 3D partial reconstruction, the reconstruction's state halves (`rec_state_half`, also PCG's) and the
  batch scoring take each used particle's half from `oris%get_half`, which stops on an `eo = -1`
  particle (in no half; serially, before the parallel regions); only an invalid state is skipped with a
  count.
- An empty partition still writes its sigma2 range, its orientations, zero PCG accumulators (3D), empty class
  sums (2D), a zero readout record (`prob*`, `write_empty_readout`) and `JOB_FINISHED`; without them the
  master's `merge_algndocs`, sigma2 commit or `merge_prob_readouts` stops.
- The matchers are the only callers of `build_batch_particles*`: no other program prepares the alignment's
  particles. Workers write `outfile` (rows `fromp:top`), sigma2 ranges, partials, contributions and, under
  `prob*`, their readout record `prob_readout_part<NN>.txt` (merged by the master in part order); never the
  project.

## Where to start

- The matchers; `simple_matcher_ptcl_batch.f90` (batches, `prep_sigmas_objfun`), `simple_matcher_ptcl_io.f90`
  (read sessions, particle and record buffers), `simple_matcher_refvol_utils.f90` (3D references),
  `simple_matcher_pftc_prep.f90` (2D references), `simple_matcher_smpl_and_lplims.f90` (sampling, band),
  `simple_strategy3D_srch.f90`/`simple_strategy2D_srch.f90` (the per-particle base), `simple_strategy3D_cont.f90`
  (Cartesian `refine=cont`), `simple_strategy2D_utils.f90` (class-average registration, clustering) and
  `simple_pspecs.f90` (the class averages' power spectra it compares).

## Tests

- `build/bin/simple_test_exec test=unit_pftc_align2D3D` (about 7 s): `suite=refine3d_in_plane_state`,
  `cavg_registration`, `2d_probability_table`, `3d_probability_table`; `unit_cart_align3D suite=pose_strategy`
  (about 5 s); `lib_cart_align3D` (1JYX recovery).
- The pipeline's draws: `test=polar_matching_bench suite=fills` (each route's batch pipeline at one thread
  against the full team, other batch sizes and three parts; fails on any differing entry).
- A change in what a pass computes: the workflow checks of `simple-main-strategies`.

## See also

`doc/policies/particle_read_path_policy.md`, `doc/policies/3D/refine3D_policy.md` sections 6-8,
`doc/policies/importance_sampling_fractional_update_policy.md`; skills `simple-main-pftc`, `simple-main-class`.
