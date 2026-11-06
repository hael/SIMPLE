---
name: simple-solve2d
description: Use when working on SIMPLE's solve2D workflow and its refine2D stages, including the public refine2D program and its private worker, 2D search strategies, the probabilistic batch pipeline of refine2D_exec (eul_prob_tab2D, keyed class draws, readout records), class-sum accumulation in the 2D matcher, class-average restoration and convergence in the 2D strategy, distributed class-average assembly (cavgassemble_exec), carried class sums, cropped class sums, project updates, and streaming/rejection helpers.
---

# SIMPLE Solve2D

`solve2D` is de novo 2D class-average determination, the 2D equivalent of
`solve3D`: one staged schedule couples ab initio 2D classification with
initial 2D refinement of the classes, mirroring the 3D methodology. Its stages
are `refine2D` runs, the 2D counterpart of `refine3D`.

Use this skill for the orchestrating `solve2D` workflow and the `refine2D`
stages it drives. 2D has the 3D structure: the matcher searches and
accumulates the class sums, the strategy restores the class averages and
decides convergence.

## Read First

1. `doc/policies/2D/solve2D_policy.md` (with
   `doc/policies/particle_read_path_policy.md` for particle reads)
2. `src/main/ui/simple/simple_ui_refine2D.f90`
3. `src/main/exec/simple_exec_refine2D.f90`
4. `src/main/commanders/simple/simple_commanders_refine2D.f90`
5. `src/main/strategies/parallelization/simple_refine2D_strategy.f90`
6. `src/main/strategies/search/simple_strategy2D_matcher.f90`
7. `src/main/strategies/parallelization/simple_cavgassemble_strategy.f90`
8. `src/main/class/simple_classaverager.f90` (`_restore`, `_core`),
   `simple_cavg_sums.f90`, `simple_class_frcs.f90`

If fractional update or probabilistic sampling is involved, also read
`doc/policies/importance_sampling_fractional_update_policy.md`.

## Ownership Map

- `solve2D` (`simple_commanders_solve2D.f90`): orchestrating workflow, staged
  `refine2D` calls, terminal greedy pass, final class averages; stage rules in
  `src/main/solve/simple_solve2D_controller.f90`.
- `solve2D_chunks`: `commander_solve2D_chunks` in `simple_commanders_solve2D` runs one `solve2D` per
  stack-bound subset; its chunk type is `solve2D_chunk` (`src/main/solve/simple_solve2D_chunk.f90`).
- `refine2D`: public `simple_exec` program (`commander_refine2D`), shared
  memory or the distributed master with `nparts`; `simple_private_exec
  prg=refine2D_distr` is the same commander. `simple_private_exec prg=refine2D`
  is the per-partition worker (`commander_refine2D_distr_worker`).
- `simple_refine2D_strategy`: execution mode, iteration control, shared-memory
  restoration (`restore_cavgs_inmem`), the one convergence decision of both
  modes (`check_refine2D_convergence`), and the distributed master's own
  project: it draws the iteration's sample (`sample_ptcls4prob2D` under
  `refine=prob*`, else `sample_ptcls4update2D` when `master_draws_sample2D`
  holds), merges the parts' readout records (`merge_prob_readouts`) and runs
  `cavgassemble_exec` and `commit_sigma2_exec` on its own `(params, build)` at
  the thread count each needs and restores its own after.
- `simple_strategy2D_matcher` (`refine2D_exec`): draws the sample (shared
  memory) or reproduces the master's (a part), under `refine=prob*` fills and
  draws per batch (`pipeline_batch`), searches, writes orientations and sigma2 rows, accumulates the class sums.
  A worker writes its contribution (`cavger_write_contribution`); shared
  memory returns with the sums in the class averager and `l_frac`.
- `simple_cavgassemble_strategy` (`cavgassemble_exec`; `run_cavgassemble`
  behind `commander_cavgassemble`): reduces the contributions, restores,
  writes, makes the class document; updates the caller's project in memory.
- Probabilistic pipeline: `eul_prob_tab2D`
  (`strategies/search/probabilistic/simple_eul_prob_tab2D.f90`) per batch
  inside `refine2D_exec` (`new_pipeline`, `fill_tab_range`, `assign_batch`,
  `close_assignment`); a distributed part writes `prob_readout_part<NN>.txt`.
- `simple_strategy_sigma2`: shared sigma2 helpers
  (`prepare_canonical_sigma_update`, `commit_sigma2_exec`).
- `src/defs/simple_refine2D_fnames.f90`: class-average file names
  (`cavgs_iterNNN`, `_even`, `_odd`, `_ranked`); never build them by hand.
- `simple_matcher_ptcl_batch`: batch loading, preparation, `crop_particles2D`.
- `class/` and `simple_commanders_mkcavgs.f90`: class-average products and the
  `make_cavgs`/`cavgassemble` command-line wrappers.
- `project/` and `ori/`: particle/class metadata persistence.
- `stream/`: online and mini-stream 2D variants (the pool runs
  `refine2D_distr` with `stream2d=yes`).

## Working Rules

- Keep command/controller code orchestration-focused and matcher code focused
  on particle-domain alignment and class-sum accumulation; restoration and
  convergence belong to the strategy.
- Preserve the sample-once-and-reuse contract: the project owner
  (distributed master or shared-memory `refine2D_exec`) draws the subset once
  per iteration (keyed under `refine=prob*`) and distributed parts reproduce
  it (`sample4update_reprod`) rather than resampling.
- Preserve `sampled` as the current-round marker and `updatecnt` as cumulative
  update history when 2D workflows share the fractional-update machinery.
- Preserve the online single-read particle I/O contract: the matcher reads
  each batch once and sums the class averages from the images already held
  (full size, or Fourier-cropped). The pipeline fills the table on the
  matcher's prepared batch, so a probabilistic iteration reads and prepares
  each particle once; 2D has no particle buffer.
- Keep the cropped-particle decision single: `crop_particles2D` (not
  `stream2d`, `box_crop < box`, `ptcl2D`) decides how `build_batch_particles2D`
  prepares the batch, cropped or not, that the fill and the search then share.
  The particle disk cache was retired 2026-11-02.
- Under the polar shift map's key (the default) a particle on the disk is
  prepared under the box's mask (`ptcl_on_disk`, `prepimg4align`); keep the
  fill's rule and the preparation's reading the same field (`shsig`).
- `cavger_update_sums` buckets a batch by (class, half) and runs the buckets
  under a dynamic schedule, each in batch order, so the sums do not depend on
  the thread count. Keep that order in any change to the sums.
- The 2D table keeps no file: the fill's thread buffers become the batch's
  dense rows (`refine=prob`, `assign_rows_dense`) or candidate store
  (`refine=prob_neigh`, `assign_rows_sparse`), released after the batch's
  draws; each class draw is keyed on the particle (`RND_FAM_ASSIGN`). Change
  the table and `simple_eul_prob_tab2D_tester` together.
- `check_conv2D` prints no verdict: `check_refine2D_convergence` prints one
  `>>> CONVERGED: yes|no` (`convergence%print_verdict`) after `minits`,
  `maxits` and the full-assignment coverage, and runs after restoration,
  whose class document writes the particles' `res`.
- Only the assembly owner (distributed assembly, shared-memory restoration)
  reads, blends and publishes the carried sums `cavg_state.bin`; workers
  accumulate from zero. Without a usable carried set the iteration is a full
  update (`carryover_needs_full_update`).
- After distributed assembly the master writes only the `cls2D`, `cls3D`,
  `out` and particle segments; do not reread the particle segment there. It
  writes `ptcl2D` more than once per iteration (once the sample is drawn, in
  `merge_algndocs`, after assembly and at the iteration's end; `doc/policies/2D/solve2D_policy.md`,
  "Project segments"): do not add a write.
- Treat explicit offline or terminal class-average assembly commands as
  separate workflow stages. They may perform their own reads, but that does not
  relax the matcher worker's online single-read contract.
- For workflows that take a project file, set `mkdir=yes` before
  `params%new(cline)` unless the command is already a child command. Child
  workflow calls from inside the execution directory should use `mkdir=no` and
  the copied local project file.
- Treat previous 2D `eo` values as input project state. Do not regenerate them
  with `partition_eo` in a workflow whose purpose is to consume prior 2D class
  assignments (`ensure_field_halves` labels only rows without halves);
  validate `eo` with `get_eo()` equal to `0` or `1`.

## Tests

`simple_test_exec test=unit_pftc_align2D3D` (`2D probability table`, `3D probability table`),
`test=unit_numerics` (`random draws`) and `test=unit_reconstruction` (`class-average carry-over`); for behaviour, a
`solve2D` run in shared memory and distributed.

## Nearby Skills

- For stream particle sieving, read the Particle Sieving section of `.github/skills/simple-main-stream/SKILL.md`.
- For shared fractional-update contracts, read `.github/skills/simple-frac-update-trailing/SKILL.md`.

## Common Traps

- `refine2D` logic is spread across many `simple_strategy2D_*` modules.
- Streaming `refine2D` code is related but not identical to batch execution.
- Restoration takes the cut-off of the even/odd low-resolution averaging from
  the previous iteration's class FRCs (`params%frcs`): the shared-memory matcher
  reads them, and the distributed master passes `prev_frcs=.true.` to
  `cavgassemble_exec`. Keep the two modes equivalent; `make_cavgs` starts from
  none in both.
- Class averaging, selection, rejection, and project updates often live in
  different files than the search loop itself.
