---
name: simple-main-strategies-probabilistic
description: Use when working in SIMPLE's src/main/strategies/search/probabilistic: the probability tables of refine=prob and prob_neigh (eul_prob_tab dense 3D, eul_prob_tab_neigh, eul_prob_tab2D dense and sparse) as the batch pipeline inside refine3D_exec and refine2D_exec (the iteration's model, fill_tab_range, the per-row draw assign_batch, the weighted insertion's extras), the keyed draws, particle temperatures and prob_beta, the parts' readout records, the strategy3D_prob and strategy2D_prob consumers, and the shared table utilities (shift prior, insertion cone, n_eff readout, the classic modes' samplers).
---

# `strategies/search/probabilistic`: the probability tables

Owns the probabilistic assignment as a batch pipeline inside the matcher pass (`refine3D_exec`, `refine2D_exec`),
the same in shared memory and in every distributed part. Per pass: the table (`new_pipeline`, `new_neigh_pipeline`),
the iteration's model (`set_model_from_field`), `open_assignment`; per batch: `fill_tab_range` on the prepared
batch, `assign_batch` (one keyed draw per particle from its row), then `strategy3D_prob`/`strategy2D_prob` polish the
drawn pose; `close_assignment` reports, or in a distributed part writes its readout record. Pair scores, in-plane and
disk draws are `pftc/`; the sample is the project owner's (`sample_ptcls4prob3D/2D`). There is no global table, part
file, assignment or insertion file, and no table program.

## Rules

Silent ones first.

- Every draw is keyed (`rnd_keyed(rnd_draw_key(scope, particle, family, ...))`, `simple_rnd.f90`) on
  `prob_rnd_scope(which_iter, l_2D)`: families `INPL` (pftc), `DISK` (shift map), `CAND` (`keyed_candidate_order`),
  `ASSIGN` (the row draw), `SAMPLE` (the owner's sample). A new draw in a fill or a selection takes a key, never
  `ran3`: an unkeyed draw depends on threads, batches and parts again. The fill sets the scope on the calculator
  (`pftc%set_rnd_scope`); a posterior or disk draw without it stops.
- The key holds the run's seed, `which_iter` and the particle's row in the project: a part with another
  `SIMPLE_SEED` or `which_iter` than the master's draws differently without error (an unset seed stops).
  `adopt_run_seed` records the seed and `RND_KEYED_VERSION` in `projinfo`; `continue=yes` stops on a version mismatch.
- A fill reads whole-field quantities only from the iteration's model, set once per pass before the first batch
  (`set_model_from_field`: the dense 3D budget `calc_num2sample` over `dist`; under `shift_prior=yes`
  `shift_prior_sigma2` over the sampled `shincarg`); otherwise only its particle's row. A whole-field read inside a
  fill makes a batch depend on the commits of earlier batches. The fill stops when the model is not set.
- Fill and draw alternate per batch: dense 3D `keep_batch_rows` and 2D `batch_rows_from_buffers` stop when the
  previous batch is not assigned; `assign_batch` stops on rows or buffered candidates of another batch.
- Complete runs still differ across thread counts (the preparations' arithmetic, summation order;
  `refine3D_policy.md` section 8): byte identity is judged at one thread.
- Fallbacks without an error: a dense row without a valid distance takes the first reference (class); a 3D
  particle with `istate < 1` keeps its pose; a neighbourhood particle without a candidate is re-seeded at its
  previous orientation (`seed_fallback_if_empty`), or keeps it (one warning with the count); a row without a
  temperature takes its argmin (warned). Dense 3D stops when every row is invalid or flat; sparse 2D on a particle
  without a candidate.
- Readout records: a distributed part writes `prob_readout_part<NN>.txt` at `close_assignment` (a zero record when
  it holds no sampled particle); the master's `merge_prob_readouts` reads them in part order, sums counts and sums
  (never means), reports and deletes them. The reader stops on a missing file and reads the keys in a fixed order:
  a new field changes `write_prob_readout`, `read_prob_readout` and `add_prob_readout` together.
- A part reproduces the owner's sample (`sample4update_reprod`, may be empty) and stops when the field was never
  sampled. The tables need `objfun=euclid`; `update_missing` with a probabilistic mode stops; so does a
  probabilistic part in the stream's 2D pool.
- `topk > 1`: the matcher sets the insertion cone from the states' FSC files (`insert_cone_from_fsc`; a state without
  one gets no cone; the master, not the parts, reports it); `open_assignment` stops without the direction space; the
  extras go per batch into `prob_insert_poses` (`setup_range`), resolved by `strategy3D_prob` and inserted by the
  partial reconstruction.
- The pass's context memory is checked against `ptcl_buffer_budget` with a warning, not a stop.

## Where to start

- `simple_strategy3D_matcher.f90`, `simple_strategy2D_matcher.f90` (`sample_particles_for_update`, `begin_pipeline`,
  `pipeline_batch`): the pipeline's order.
- `simple_eul_prob_tab.f90` (dense 3D: model, fill, `assign_rows`, extras, readout), `simple_eul_prob_tab_neigh.f90`
  (extends it: `assign_rows_neigh`, the re-seed), `simple_eul_prob_tab2D.f90` (`assign_rows_dense`/`_sparse`),
  `simple_eul_prob_tab_utils.f90` (`draw_posterior_row`, temperatures, `prob_beta_at`, shift prior, insertion cone
  and table, n_eff and readout records, `prob_rnd_scope`, `keyed_candidate_order`, the classic modes' samplers on
  `ran3`). Generator: `src/utils/math/simple_rnd.f90`; sample and seed: `search/simple_matcher_smpl_and_lplims.f90`.

## Tests

- `build/bin/simple_test_exec test=unit_pftc_align2D3D suite=3d_probability_table` and `suite=2d_probability_table`
  (selection rules on a pipeline batch, draws, extras, n_eff); `test=unit_numerics suite=random_draws` (the
  generator's pinned vectors).
- `test=polar_matching_bench suite=fills`: each route's pipeline at one thread against the team, in batches and parts,
  compared by particle (`polar_fill_identity.tsv`). Without `suite` it times the kernels.
- Workflows: `simulated_workflow` (solve2D and solve3D run the probabilistic modes; 1jxy also over two parts), judged
  as a profile. No unit tester covers the readout records or `merge_prob_readouts`; only the 1jxy two-part run
  exercises them.

## See also

`doc/policies/3D/refine3D_policy.md` section 8, `doc/policies/2D/solve2D_policy.md`,
`doc/policies/importance_sampling_fractional_update_policy.md`,
`doc/implementation_notes/completed/prob_align_batch_pipeline.md` (the pipeline and the key schema, sections 3-6) and
its `_report.md`; skills `simple-main-strategies-search`, `simple-solve3d-importance-sampling`, `simple-main-pftc`.
