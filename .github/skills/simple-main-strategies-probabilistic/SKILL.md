---
name: simple-main-strategies-probabilistic
description: Use when working in SIMPLE's src/main/strategies/search/probabilistic: the probability tables of refine=prob and prob_neigh (eul_prob_tab dense 3D, eul_prob_tab_neigh, eul_prob_tab2D dense and sparse), their partition part files and format words, the master's read and per-particle posterior draw (ref_assign, draw_posterior_row, particle temperatures, prob_beta), the assignment and weighted-insertion (topk) files, the strategy3D_prob and strategy2D_prob consumers, and the shared table utilities (shift prior, insertion cone, SHC samplers).
---

# `strategies/search/probabilistic`: the probability tables

Owns the two-phase probabilistic assignment. Table workers (`prob_tab`, `prob_tab_neigh`, `prob_tab2D`,
their frames in `commanders/simple/simple_commanders_prob.f90`) fill one part file each, `fill_tab_range`
per particle batch; the master (`prob_align*_exec`, `parallelization/`) reads every part into one table,
draws each particle's reference (`ref_assign`) and writes the assignment; the matcher's
`strategy3D_prob`/`strategy2D_prob` read it and polish the drawn pose. Pair scores and in-plane draws are
`pftc/`; the sample is the master's.

## Rules

Silent ones first.

- The fills read the incoming particle field: at each `fill_tab_range` dense 3D derives the refined-reference
  budget (`calc_num2sample` over the sampled particles' `dist`) and, under `shift_prior=yes`, every shifted fill
  the shift prior's width (`shift_prior_sigma2` over `shincarg`). Every fill of an iteration runs before any
  of its decisions reaches the field; otherwise the budget and the prior depend on the batch.
- Each fill calls `seed_rnd` and fills run under a dynamic schedule: in-plane, disk and candidate draws
  depend on the thread, so byte identity holds at one thread only. Every `ref_assign` draws in particle order
  on the master's stream, so shared-memory and distributed runs draw alike there.
- Fallbacks without an error: a row without a valid distance takes the first reference (class); a 3D
  particle whose assignment has `istate < 1` keeps its pose; neighbourhood particles without a candidate keep
  their orientation (one warning with the count). `read_assignment` keeps only the rows of the worker's `pinds`.
- Part files: `DIST_FBODY<part>.dat` (dense 3D, format word -4; 2D dense -7, 2D candidate stream -6) and
  `DIST_FBODY_neigh_<part>.dat` (no format word). A layout change changes writer, reader and format word
  together. The dense 3D and 2D reads stop on a missing part, a wrong format, size or reference count, or a
  particle in no part or in two; the neighbourhood read checks the reference and candidate counts only (its
  missing coverage check is left to the batch-pipeline refactoring, which retires these readers).
- Every part writes a file, empty when it holds no sampled particle. A table worker reproduces the master's
  sample and stops when the field was never sampled.
- The master allocates its table after the workers finish (peak memory). The 3D `prob_align*_exec` delete the
  previous part files and assignment first; distributed 2D does it in `cleanup_distributed_iteration_artifacts`.
- `topk > 1`: the master needs the direction space (`run_prob_align` builds it; `ref_assign` stops without
  it), sets the insertion cone from the states' FSC files (`insert_cone_from_fsc`) and writes
  `INSERT_FBODY.dat`; the 3D matcher stops when that file is missing.
- In shared memory the 2D table fills the particle buffer (cropped: the records) that the 2D matcher of the
  same iteration copies from and releases: a change to either side's reads changes the other.

## Where to start

- `simple_eul_prob_tab.f90` (dense 3D: fill, part I/O, `ref_assign`, assignment and insert-table I/O),
  `simple_eul_prob_tab_neigh.f90` (extends it), `simple_eul_prob_tab2D.f90`, `simple_eul_prob_tab_utils.f90`
  (`draw_posterior_row`, temperatures, shift prior, insertion cone and table, SHC samplers).

## Tests

- `build/bin/simple_test_exec test=unit_pftc_align2D3D suite=3d_probability_table` and
  `suite=2d_probability_table_io` (under a second each; the whole suite about 7 s).
- Draws or scores: `simulated_workflow` (solve2D and solve3D run the probabilistic modes), judged as a
  profile; `test=polar_matching_bench` for the table's kernels.

## See also

`doc/policies/3D/refine3D_policy.md` section 8, `doc/policies/2D/solve2D_policy.md`,
`doc/implementation_notes/completed/prob_align_batch_pipeline.md` (the refactoring that removes the global table,
the part files and the per-fill reseed); skills `simple-solve3d-importance-sampling`, `simple-main-pftc`.
