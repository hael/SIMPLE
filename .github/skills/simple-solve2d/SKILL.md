---
name: simple-solve2d
description: Use when working on SIMPLE's solve2D workflow and its refine2D stages, including the public refine2D program and its private worker, 2D search strategies, prob_tab2D/prob_align2D part files, class-sum accumulation in the 2D matcher, class-average restoration and convergence in the 2D strategy, distributed class-average assembly (cavgassemble_exec), carried class sums, the shared-memory particle buffer, cropped class sums and the particle disk cache, project updates, and streaming/rejection helpers.
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

1. `doc/policies/2D/solve2D_policy.md` (with `particle_cache_policy.md` and
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
  project: it runs `prob_align2D_exec`, `cavgassemble_exec` and
  `commit_sigma2_exec` on its own `(params, build)` at the thread count each
  needs and restores its own after.
- `simple_strategy2D_matcher` (`refine2D_exec`): reproduces the sample,
  searches, writes orientations and sigma2 rows, accumulates the class sums.
  A worker writes its contribution (`cavger_write_contribution`); shared
  memory returns with the sums in the class averager and `l_frac`.
- `simple_cavgassemble_strategy` (`cavgassemble_exec`; `run_cavgassemble`
  behind `commander_cavgassemble`): reduces the contributions, restores,
  writes, makes the class document; updates the caller's project in memory.
- `prob_tab2D` (`exec_prob_tab2D` in `simple_commanders_prob.f90`) and
  `prob_align2D_exec` (`simple_prob_align_strategy.f90`); tables in
  `strategies/search/probabilistic/simple_eul_prob_tab2D.f90`.
- `simple_strategy_sigma2`: shared sigma2 helpers
  (`prepare_canonical_sigma_update`, `commit_sigma2_exec`).
- `src/defs/simple_refine2D_fnames.f90`: class-average file names
  (`cavgs_iterNNN`, `_even`, `_odd`, `_ranked`); never build them by hand.
- `simple_matcher_ptcl_batch`: batch loading, preparation, `crop_particles2D`.
- `simple_ptcl_cache`: downscaled particle disk cache (`cache=yes`) shared by
  `prob_tab2D` and `refine2D_exec`; owns its on-disk lifecycle.
- `class/` and `simple_commanders_mkcavgs.f90`: class-average products and the
  `make_cavgs`/`cavgassemble` command-line wrappers.
- `project/` and `ori/`: particle/class metadata persistence.
- `stream/`: online and mini-stream 2D variants (the pool runs
  `refine2D_distr` with `stream2d=yes`).

## Working Rules

- Keep command/controller code orchestration-focused and matcher code focused
  on particle-domain alignment and class-sum accumulation; restoration and
  convergence belong to the strategy.
- Preserve the probabilistic sample-once-and-reuse contract:
  `prob_align2D` chooses the sampled subset, and `prob_tab2D` plus the matcher
  reproduce that subset rather than resampling independently.
- Preserve `sampled` as the current-round marker and `updatecnt` as cumulative
  update history when 2D workflows share the fractional-update machinery.
- Preserve the online single-read particle I/O contract: the matcher reads
  each batch once and sums the class averages from the images already held
  (full size, or Fourier-cropped). In shared memory `prob_tab2D` fills the
  particle buffer, and the matcher copies its batches from it and releases it
  (`release_pass_buffer`), so a table iteration reads each particle once.
- Keep the cropped-particle decision single: `crop_particles2D` (uncached,
  not `stream2d`, `box_crop < box`, `ptcl2D`) decides for the table and the
  matcher, which then prepare particles from the in-memory cache entry (in
  shared memory the table's records feed the matcher). `cache=no` equals
  `cache=yes` byte for byte; keep it so.
- `cavger_update_sums` buckets a batch by (class, half) and runs the buckets
  under a dynamic schedule, each in batch order, so the sums do not depend on
  the thread count. Keep that order in any change to the sums.
- The 2D part files carry a format word: -4 for the candidate stream
  (`prob_snhc`, with per-particle counts), -5 for the dense rows (`prob`:
  distance and in-plane index per entry, then the refined entries).
  `prob_align2D` reads them in parallel and checks size, format and that each
  particle is in exactly one part. Change writer, readers and the tester
  together.
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
  writes `ptcl2D` more than once per iteration (in `merge_algndocs`, after
  assembly and at the iteration's end; `doc/policies/2D/solve2D_policy.md`,
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

## Particle Disk Cache Lifecycle

- `ptcl_cache_ensure` runs master-side only, before workers are scheduled, and
  must stay ahead of `gen_job_descr`/cline copies: fallback decisions propagate
  to worker command lines by flipping `cache=no` on the cline. Cache
  availability must be uniform across ranks; never make it best-effort per rank.
- The building (or adopting) process owns the cache files and removes them on
  exit — normal or hard exception — via the `cache_cleanup_glob` hook in
  `simple_defs` (armed by `ptcl_cache_own`, fired from `simple_exception`, module
  `simple_error`, and
  the tails of `simple_exec`/`single_exec`). Workers never take ownership.
- The cache basename hashes the execution directory so concurrent runs sharing
  a `cache_dir` cannot collide; all ranks recompute the name locally because
  qsys workers cd into the master's directory. The key file records the
  directory verbatim, geometry, and per-stack fingerprints; leftovers with the
  same name are validated and adopted or rebuilt, never trusted by name alone.
- The build refuses to claim more than 25% of the free space at `cache_dir`
  and falls back to uncached execution with a log message.
- Delete-on-exit means no cross-run cache reuse; within-run reuse across
  solve2D stages is the point (box_crop is stage-invariant), so do not add
  per-stage invalidation or per-stage names.

## Tests

`simple_test_exec test=unit_pftc_align2D3D` (`2D probability table I/O`) and
`test=unit_reconstruction` (`class-average carry-over`); for behaviour, a
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
