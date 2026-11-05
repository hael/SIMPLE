---
name: simple-main-strategies-search
description: Use when working in SIMPLE's src/main/strategies/search: the refine3D_exec and refine2D_exec matchers (a worker's or a shared-memory pass), particle batches (build_batch_particles2D/3D, read sessions, the particle and record buffers), reference preparation (reprojection model, 2D references, the band), a pass's sigma2 preparation, the per-particle 2D/3D strategies (greedy, shc, snhc, inpl, neigh, prob, cont, tseries) and their per-thread state, sampling and low-pass helpers, the classic polar shift map, 3D partial reconstruction, class-average registration utilities.
---

# `strategies/search`: the matchers and the per-particle searches

Owns one pass over a partition's particles: `refine3D_exec` (`simple_strategy3D_matcher.f90`) and
`refine2D_exec` (`simple_strategy2D_matcher.f90`) sample, prepare the references and particle batches, run
one `strategy3D`/`strategy2D` object per particle and write the part's outputs; the 3D pass then writes its
partial reconstructions, the 2D pass its class sums. It does not draw a distributed run's sample or merge
(`parallelization/`), score a pair (`pftc/`), restore class averages (`class/`) or commit sigma2 (`sigma2/`).
The probability tables are `probabilistic/` (`simple-main-strategies-probabilistic`).

## Rules

Silent ones first.

- 3D pass: `ensure_field_halves` -> `read_reprojection_model` (`pftc%new`) -> `prep_sigmas_objfun` ->
  `alloc_ptcl_imgs` -> `memoize_refs` -> `prep_strategy3D` -> per batch `build_batch_particles3D`, the loop ->
  `write_sigma2`, orientations -> `pftc%kill` (stops if still allocated) -> partial reconstructions ->
  `release_pass_buffer` -> `qsys_declare_part_finished`. 2D: `cavger_new`, `cavger_read_all`,
  `cavger_init_online` -> `set_bp_range2D` -> `prep_pftc4align2D` (`pftc%new`, `prep_sigmas_objfun`,
  references, centring offsets, `memoize_refs`) -> per batch `build_batch_particles2D`, the loop,
  `cavger_update_sums` -> sigma2 range, orientations, a worker's class-sum contribution.
- `prep_sigmas_objfun` after every `pftc%new`: the calculator points at `build%esig`'s spectra and `new`
  drops the pointer; without them `memoize_sigma2_weights` returns with no weights.
- A worker never makes its 3D band: `adopt_reprojection_model_range` takes `kfromto` and `lp` from the
  reprojection model's header (Cartesian: `read_cart_refvols`) and stops on an incompatible one.
- In the particle loop a contained procedure sees the host's variables, not the loop's `firstprivate`
  copies: declare per-particle objects in the procedure (`o_sigma` in `choose_and_run_strategy`). Per-thread
  state is `s3D`/`s2D` columns `ithr`; draws in the loop (`ran3`, `ran_tabu`) use the thread's stream, so
  under the dynamic schedule a particle's draws depend on the thread.
- A worker reproduces the sample (`sample4update_reprod`) whenever the master drew it (probabilistic
  modes, `master_draws_sample*`, the polish) and draws only otherwise.
- Preparation records each particle's mask (`note_obs_mask`); class sums and 3D reconstruction normalise
  under the same one (`obs_on_box`). A cropped 3D record (`put_record` in `polarize_batch_particles3D`) is
  the reconstruction's observation bit for bit (`prep_rec_observation`): change both together.
  `normalized_eo` (`simple_matcher_3Drec.f90`) puts an `eo = -1` particle (in no half) in the even half:
  known defect, ruled to stop instead; rely on it nowhere.
- An empty partition still writes its sigma2 range, its orientations, zero PCG accumulators (3D), empty class
  sums (2D) and `JOB_FINISHED`; without them the master's `merge_algndocs` or sigma2 commit stops.
- The table workers in `commanders/simple/simple_commanders_prob.f90` repeat the matchers' preparation
  (references, sigma2, `build_batch_particles*`): change both. Workers write `outfile` (rows `fromp:top`),
  sigma2 ranges, partials, contributions; never the project.

## Where to start

- The matchers; `simple_matcher_ptcl_batch.f90` (batches, `prep_sigmas_objfun`), `simple_matcher_ptcl_io.f90`
  (read sessions, particle and record buffers), `simple_matcher_refvol_utils.f90` (3D references),
  `simple_matcher_pftc_prep.f90` (2D references), `simple_matcher_smpl_and_lplims.f90` (sampling, band),
  `simple_strategy3D_srch.f90`/`simple_strategy2D_srch.f90` (the per-particle base), `simple_strategy3D_cont.f90`
  (Cartesian `refine=cont`), `simple_strategy2D_utils.f90` (class-average registration, clustering).

## Tests

- `build/bin/simple_test_exec test=unit_pftc_align2D3D` (about 7 s): `suite=refine3d_in_plane_state`,
  `cavg_registration`; `unit_cart_align3D suite=pose_strategy` (about 5 s); `lib_cart_align3D` (1JYX recovery).
- A change in what a pass computes: the workflow checks of `simple-main-strategies`.

## See also

`doc/policies/particle_read_path_policy.md`, `doc/policies/3D/refine3D_policy.md` sections 6-8,
`doc/policies/importance_sampling_fractional_update_policy.md`; skills `simple-main-pftc`, `simple-main-class`.
