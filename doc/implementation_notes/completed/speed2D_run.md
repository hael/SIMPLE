# 2D speed and finishing the identical speed changes: the Mac run `speed2d_mac`

Implementation note, 2026-10-18. **Carried out 2026-10-18/19** (progress
below; report in `speed2D_run_report.md`). The
execution plan of the unattended night run `speed2d_mac` on the maintainer's
Mac. The design is in `speed_architecture_proposals.md` (section 1.1, the
finishing fixes, and section 1.2, the 2D workflow) and in
`identical_speed_run_review.md`, both under `doc/implementation_notes/planned/`
in `~/src/SIMPLE` (read them there; the proposals note there may be newer than
the copy in the worktree). This note holds the run's items, tests, phases and
progress.

## Rulings of 2026-10-18

- Tonight's chunk: what is left of the finishing fixes, then the 2D speed
  items. Every item leaves results byte-identical at one thread.
- Testing: one-thread identity on the simulated cases and on new 2D cases,
  and, for every new or changed parallel region, a thread stress test at the
  Mac's full thread count, repeated; one multi-threaded workflow run per
  phase to catch crashes and races.
- Not tonight: the particle cache stays off by default for `solve2D`, and the
  `solve2D` fractional-update default stays (ruled). The changes that alter
  results (the master drawing the stage 1-2 sample, partition balancing by
  active particles, the `maxits`/`minits`, `fillin` and `prob_snhc` renames)
  wait for the truth suite of the run `test_truth`. The 2D look-and-feel items
  and the compact 2D table record are the next night. "One project owner on
  the master" is a later run.

## Why thread stress testing

The previous run's parallel read of the dense table parts passed one-thread
identity and a one-against-four-thread tester, yet held a race the maintainer
fixed in commit 43de9b340: string temporaries built inside the OpenMP loop
and handed to `fopen` were not thread-safe, a coverage check needed `atomic
capture`, and `ref_assign` sized a work array by `nthr_glob` while its loop
ran on `omp_get_max_threads()` threads. One thread runs an OpenMP loop
serially, so identity can never show thread safety. Rules for every parallel
region this run adds or changes:

- No string or allocatable-character temporaries, and no calls into
  routines with module-level buffers, inside the region: build file names
  and messages before it.
- Work arrays sized by the team the region actually uses.
- Errors inside the region set a status per task, raised after the region.
- A tester check that runs the routine at one thread and at the Mac's full
  thread count (`omp_get_num_procs()`), at least 20 times, and compares bit
  for bit where the result is deterministic by construction; where it is
  not (dynamic schedules that move random draws), the check is that every
  run completes with finite results.

## The items

### Phase 1: finishing the identical speed changes

The maintainer's commit 43de9b340 already validates the dense part file's
header and counts, detects a particle in several parts with `atomic
capture`, builds the part file names before the parallel loop, and sizes the
assignment's work array by the OpenMP team. What is left:

1. **A short or unreadable dense part file stops the run.** The rows and the
   refined entries are still checked with `fileiochk`, which only prints at
   end of file; stop on any nonzero `iostat`, and compare the file size with
   the size its header implies, before reading rows.
2. **The active states in the dense part file.** Write the list of active
   states after the header and compare it on read, so that a row's reference
   index is checked as the old format allowed.
3. **Errors after the parallel read.** `read_dense_part` still stops from
   inside the parallel region; give it a status per part and raise the first
   failure, with its file name, after the loop.
4. **Volume assembly's handed-on project freed early.** Take the copy after
   the reconstructors are freed, make the strategy's `vol_proj` allocatable,
   and release it after its last use (or hand back only the fields assembly
   changes, whichever is simpler and identical).
5. **Dead code**: `gen_corr_grad_for_rot_8` (its only caller is a tester,
   which compares against a fresh load instead), the per-thread buffers
   `kcorrs` and `pft_r1_8` (never used), the kept builder's release in
   `distr_initialize` (cannot fire) and its comment; stale comments (the
   "candidate buffer" of the dense fill, the tester section header that
   names a run).
6. **The "load before evaluate" contract** of the shift search's `_aln`
   evaluators: the loaded reference's index and half kept per thread, a stop
   on mismatch.
7. **Table types**: the dense format word private to the dense table;
   `write_rows` refused with an error in the neighbourhood type; 64-bit
   offsets for the refined lists, or the size guard of
   `new_fixed_candidate_store`.
8. **Documents**: the planned note `next_release/calpha_nu_evidence_mask.md`
   names the deleted `nu_filter_vol`. Skill files are not edited; proposed
   text for `simple-main-nu-filt` and `simple-nonuniform-regularization` goes
   into the report.

### Phase 2: measuring the 2D workflow

All identical (benchmark output only).

9. `REFINE2D` written for every partition, not only the first.
10. Timers in `REFINE2D_STRATEGY` (shared memory and distributed): the
    probabilistic pre-step, scheduling and search, merging, class-average
    assembly, project rereads and writes, convergence, the star file and the
    JPEG, so its "accounted for" share means something.
11. A `PROB_ALIGN2D` report through the existing `write_prob_align_bench`.
12. Thread-seconds and load balance for the 2D search loop and the class
    sums, as the table fill reports them.
13. `scripts/plot_refine3d_bench.py`: the 2D critical path per iteration,
    like the 3D one.

Before Phase 3, run the 2D cases (below) on the Phase 2 build at several
thread counts (1, 4, 11; shared memory and two parts) and record where the
time goes. If the numbers contradict the order of Phase 3, reorder it and
say why.

### Phase 3: 2D speed

14. **Class sums balanced** (`simple_classaverager_restore.f90`, the static
    loop over classes): a dynamic schedule over (class, half) pairs, largest
    first; optionally the per-particle padding and FFT in a loop over
    particles before the accumulation per class. Each class and half still
    accumulates its particles in batch order, so the sums are identical at
    any thread count; the stress test compares them bit for bit at 1 and 11
    threads.
15. **One particle read per table iteration** in shared memory: fill the
    existing particle pass buffer (`fill_pass_buffer`,
    `simple_matcher_ptcl_io.f90`) in `prob_tab2D` and serve the search's
    batches from it; copy batches out in parallel.
16. **Search loops dynamic**: `schedule(dynamic)` in the 2D matcher and in
    the 3D matcher. Identical at one thread; at more threads only which
    random draws a particle meets changes, as accepted for the tables.
17. **The 2D dense assignment without rescans**: the `do while(
    any(ptcl_avail) )` of `simple_eul_prob_tab2D.f90` becomes a counter and
    only the heads that pointed at the assigned particle advance, as in the
    3D `ref_assign`; a tester keeps the present algorithm as the reference
    and compares assignments on synthetic tables, ties included.
18. **Populations counted once**: one `get_pops` pass instead of three
    `get_pop` calls per class in the reference preparation
    (`simple_matcher_pftc_prep.f90`).
19. **The last two heap trampolines**: array samplers for
    `gen_prob_euclid_val` and `gen_prob_power_euclid_val`
    (`simple_polarft_corr.f90`), with a sampler-equivalence tester like
    `test_likelihood_samplers`.
20. **The empty 2D partition.** A distributed 2D partition with no sampled
    particles is expected to crash (a batch size divides by zero; an empty
    table stops the run), where 3D handles it. Reproduce it first with a
    small case (for example many parts and a small sample); if it crashes,
    fix it as 3D does (a header-only table part, a guarded batch size).
    Partition balancing by active particles is not part of this item (it
    changes results).

## The identity test

As in the run `speed_mac` (`identical_speed_run.md`, section "The identity
test"), with its harness scripts in `~/speed_mac/bin` reusable after checking
them (`run_case.sh`, `compare.py`, `probtab_decode.py`, `item_src.sh`,
`snap_item.sh`, `bg`): one-thread runs with a fixed seed, each case in a
fixed directory, outputs moved aside after each run, every regular file
compared except exclusions declared in Phase 0 with reasons, project files
compared through `print_project_field` dumps (ptcl2D, ptcl3D, cls2D, cls3D,
out), table part files compared through the lossless decoder (accepted by the
maintainer). The fixture's checksums are checked after every case. The
cases run in the CTest environment with `OMP_NUM_THREADS=1`; `SIMPLE_SEED` is
set for `simple_exec` cases.

**Cases.** The fixture and every case are made again on the base build of
this run: master has moved (`solve3D` stage changes), so the previous run's
outputs are not a reference.

| case | what it covers |
| --- | --- |
| SW1, SW6 | the simulated workflows (1jxy with its two-part step, 6vxx), including their `solve2D` steps |
| MSI, CCI, PCI | the three-state `refine3D` (dense then neighbourhood table), the correlation objective with fractional update and trailing reconstruction, the PCG backend (as in the previous run, from the SW1 fixture) |
| E2 | `solve2D` with two parts at one thread on a larger simulated particle set: about 1,500 particles made by `simulate_particles` from the 1jxy truth map of the fixture (CTF, random poses and shifts), about 20 classes, all six stages (the `snhc_smpl`, `prob_snhc` and `prob` table paths, class-average assembly, the distributed master) |
| E2S | E2's command line in shared memory (the shared-memory iteration, the pass buffer) |
| R2 | `refine2D refine=prob` for two iterations from E2's final project, two parts (the dense 2D table, its reads and its assignment) |
| AD1, CR1 | `solve3D_addon` and `cont_refine3D_1jxy` at one thread, once at the end of the run |

E2, E2S and R2 must take at most about ten minutes each at one thread; size
the particle set to that. Check in Phase 0 that the 2D cases reach every
stage and table path they claim (from their logs).

**When to run what.** After every item: the item's Debug build and the
fast-gate suites, the identity cases the item can reach (all of SW1, E2, E2S
and R2 for 2D items; SW1, MSI, CCI, PCI for 3D and table items), and the
stress tests of any parallel region it touches. At the end of every phase:
SW6, the full `./compile_debug.sh` on an idle Mac, and one multi-threaded run
of E2S at 11 threads and of MSI with two parts of 5 threads (they must run to
the end with finite results; they are not byte-compared). AD1 and CR1 once,
at the end of Phase 3.

## Common rules

As in `identical_speed_run.md`, section "Common rules", with these changes:

- The base is master `43de9b340` or later, whatever `~/src/SIMPLE` holds at
  launch; the worktree is `~/src/SIMPLE-wt-speed2d`.
- Scratch is `~/speed2d_mac`; the previous run's `~/speed_mac` may be read,
  never written.
- No new test suites or registrations: new checks go into existing tester
  modules of registered suites (`unit_pftc_align2D3D` for the tables and
  kernels; for the class sums the closest registered tester of the class
  area, for example `simple_cavg_sums_tester.f90` if it is registered). The
  files the run `test_truth` edits stay untouched (the list in
  `identical_speed_run.md`).
- Fix as you go; an item that cannot be made identical is restored and
  recorded; per-item snapshots and diffs; nothing committed.

## Phases

**Phase 0: base and cases.** The Debug build with the fast gate and the base
Release build of the worktree as it is. The fixture (SW1 once, in a fixed
directory), then the E2 particle set, then every case twice on the base build
(AD1 and CR1 once each, as a reference for the end of the run), wall time and
peak memory recorded. The exclusion list from the pairs. Exit: every case
deterministic at one thread under the list; the 2D cases reach their paths.

**Phase 1** (items 1-8), **Phase 2** (items 9-13 and the 2D time budget),
**Phase 3** (items 14-20, AD1 and CR1 at its end). Exit of each: every kept
item identical on its cases and passing its stress tests, its diff written;
items not kept recorded with the cause; the phase-end runs done.

**Phase 4: documents and close.** `doc/policies/2D/solve2D_policy.md` (the
benchmark reports, the class sums, the pass buffer, the search schedule),
`doc/policies/3D/refine3D_policy.md` where Phase 1 changed what it
describes, this note's progress, the report beside it
(`speed2D_run_report.md`): what changed per item, the identity and stress
evidence, the 2D time budget from Phase 2 and the one-thread benchmark lines
before and after (as an indication), what was not kept and why, proposed
skill text, and a commit message. Then `full.diff` in the run directory.

## File table

| file | why | phase |
| --- | --- | --- |
| `doc/implementation_notes/planned/speed2D_run.md` | this plan, its progress, the report beside it | each |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90`, `_neigh.f90`, `_utils.f90`, `_tester.f90` | items 1-3, 7 | 1 |
| `src/main/strategies/parallelization/simple_refine3D_strategy.f90`, `src/main/commanders/simple/simple_commanders_rec_distr.f90` | items 4, 5 | 1 |
| `src/main/pftc/simple_polarft_calc.f90`, `_core.f90`, `_corr.f90`, `simple_pftc_shsrch_grad.f90`, `simple_pftc_inplane_tester.f90` | items 5, 6, 19 | 1 3 |
| `doc/implementation_notes/planned/next_release/calpha_nu_evidence_mask.md` | item 8 | 1 |
| `src/main/strategies/parallelization/simple_refine2D_strategy.f90` | items 10, 15, 20 | 2 3 |
| `src/main/strategies/search/simple_strategy2D_matcher.f90`, `simple_strategy3D_matcher.f90` | items 9, 12, 15, 16, 20 | 2 3 |
| `src/main/commanders/simple/simple_commanders_prob.f90` | items 11, 15, 20 | 2 3 |
| `src/main/class/simple_classaverager_restore.f90` and the class-area tester | items 12, 14 | 2 3 |
| `src/main/strategies/search/simple_matcher_ptcl_io.f90`, `simple_matcher_pftc_prep.f90` | items 15, 18 | 3 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab2D.f90` | items 17, 20 | 3 |
| `scripts/plot_refine3d_bench.py` | item 13 | 2 |
| `doc/policies/2D/solve2D_policy.md`, `doc/policies/3D/refine3D_policy.md` | documents | 4 |

## Progress

### Phase 0, 2026-10-18: base and cases

The base is master `43de9b340`; scratch is `~/speed2d_mac`, the run
directory `~/agent_runs/runs/speed2d_mac`.

- **Builds.** `./compile_debug.sh` passed its fast gate, 15 of 15 suites
  in 14.4 s. The base Release build is `~/speed2d_mac/builds/base`
  (gfortran 16 from Homebrew, Release, `BUILD_TESTS=ON`, OpenMP,
  `USE_ARCHOPT=ON`, installed into the build directory).
- **The 2D cases.** The particle set is 1,500 particles simulated from the
  fixture's 1jxy map (box 192, 1.3 A, SNR 0.1, CTF, random poses and
  shifts up to 3 pixels), imported into a project of its own. E2 is
  `solve2D ncls=20 nparts=2` on it (about 2 minutes at one thread), E2S the
  same in shared memory (about 3 minutes). Their logs show all six stages:
  `snhc_smpl` in stages 1-2, the sparse `prob_snhc` table in stages 3-5,
  the dense `prob` table in stage 6, class-average assembly by the master.
  With 1,500 particles the fractional update is off (the sample is all
  particles). R2 is `refine2D refine=prob maxits=2 nparts=2` from a fixed
  copy of E2's final project. `refine2D` is reached through
  `simple_private_exec prg=refine2D_distr`, which writes the project in
  place, so R2 runs on a copy in its case directory.
- **Determinism.** Every case ran twice on the base build with no
  differences under the declared exclusions (`exclusions.tsv` in the run
  directory): SW1 593 files, SW6 542, MSI 216, CCI 50, PCI 53, E2 140,
  E2S 137, R2 30. AD1 (44 check values, all pass) and CR1 (59, all pass)
  ran once as references for the end of Phase 3. SW1 and SW6 pass their
  truth checks.
- **MSI** on this master ends step 1 with 18, 28 and 41 particles per state
  and step 2 with 0, 20 and 67 (an empty state, whose previous volume is
  carried forward). It stays on the SW1 fixture.
- **A defect found and fixed (item 01).** The base crashed in the
  multi-threaded MSI run: a state of 9 particles got an all-zero FSC, so
  the post-processing's FSC=0.143 resolution was 0 A and the B-factor fit
  indexed outside its spectrum (SIGBUS). The post-processing now cuts at
  the first shell when no shell reaches FSC=0.143. It is identical on SW1,
  MSI, CCI and PCI. Why such a state gets an all-zero FSC is left for the
  maintainer.
- **Cost** (one thread, wall time, peak memory of the largest process):
  SW1 10 min, 3.7 GB; SW6 14 min, 3.8 GB; MSI 4.5 min, 4.3 GB; CCI 0.4 min,
  1.6 GB; PCI 2 min, 4.6 GB; E2 2 min, E2S 3 min, R2 0.3 min (all under
  0.3 GB in the master); AD1 19 min, 2.0 GB; CR1 33 min, 8.6 GB.

### Phase 1, 2026-10-18: finishing the identical speed changes

All eight items are kept. Each was built from its own snapshot and is
byte-identical to the base on SW1, MSI, CCI and PCI (the two polar-kernel
items also on E2, E2S and R2); SW6 is identical after the phase. Items 3,
1 and 2 were done in that order, so that the stops of item 1 go through
the per-part status of item 3.

- **Errors after the parallel read.** `read_dense_part` returns a status
  per part file and stops nowhere; it builds no strings inside the
  parallel loop. After the loop the first failing part stops the run, with
  its file name and the reason.
- **Short or unreadable part files.** Before reading past the header, the
  file size must be the size the header implies; every read checks its
  status.
- **The active states.** The dense part file holds the number and list of
  active states after the header (format word -3); the master compares
  them with its own, so a row's reference index is checked again.
- **Table types.** The format word is private to the dense table; the
  neighbourhood table refuses `write_rows`; the refined offsets stay
  default integers behind a size guard, as the candidate store does.
- **Dead code.** `gen_corr_grad_for_rot_8`, the per-thread buffers `kcorrs`
  and `pft_r1_8`, and the release of the kept reprojection builder in
  `distr_initialize` are deleted; stale comments fixed.
- **Load before evaluate.** `load_ref_aln_8` records the reference index
  and half per thread; the `_aln` evaluators take the reference index and
  stop on a mismatch.
- **Volume assembly's handed-on project.** Assembly copies its project out
  after freeing its reconstructors and volumes; the strategy holds it in
  an allocatable and frees it after the sigma2 commit, its last use.
- **Documents.** The planned C-alpha note uses `nu_filter_vols`. Proposed
  skill text for the two NU skills is in the report.
- **Tests.** The parallel read runs 24 times at every processor (11 on the
  Mac) against the read at one thread; the shift searches 20 times at
  every processor against one thread; the loaded reference is compared
  with a fresh load. A throwaway program (not in the repository) confirmed
  each new stop: a part cut by 4 bytes, a part with other active states,
  a particle in two parts.

### Phase 2, 2026-10-19: measuring the 2D workflow

All five items are kept. They change benchmark reports only; each was built
from its own snapshot and is byte-identical to the base on SW1, E2, E2S and
R2; SW6 is identical after the phase.

- **REFINE2D for every partition** (`REFINE2D_BENCH_PARTppp.txt`, as
  `REFINE3D`).
- **REFINE2D_STRATEGY** times the whole iteration up to the end of its
  finalisation: setup, pre-step, sigma2 preparation, matcher or scheduler,
  merge, class-average assembly, the master's project writes and rereads,
  sigma2 commit, convergence, star file and JPEG. Its "accounted for" share
  is now about 99.9%.
- **PROB_ALIGN2D** writes its record through `write_prob_align_bench`.
- **Thread-seconds and load balance** of the 2D search loop and of the
  class sums, in `REFINE2D` and in the offline `CAVGSUMS` record. The stress
  run for the two changed regions: E2S at 11 threads, 20 times, each to the
  end with finite class averages.
- **The summary script** builds the 2D critical path per iteration
  (`bench_iterations2D.csv`, `bench_profile2D.svg`, a report section by
  stage with the load balance).
- **The 2D time budget** (E2, 1,500 particles; wall time of the whole run):
  shared memory 144, 44 and 30 s at 1, 4 and 11 threads; two parts 95, 50
  and 51 s at 1, 4 and 5 threads per part. At one thread the matcher's
  class sums take 43 s, its search 38 s, particle preparation 19 s (and the
  table part's another 11 s), the table fill 25 s. At 11 threads the class
  sums keep their threads busy 63% of the time (47% in the worst iteration)
  and the search loop 85%; the table fill 99%. The distributed runs spend
  about 18 s in job launches and polling whatever the thread count. The
  order of Phase 3 stands: the class sums first, then the particle pass
  and the dynamic search schedule.

### Phase 3, 2026-10-19: 2D speed

All seven items are kept, and one defect was fixed on the way. Each item was
built from its own snapshot and is byte-identical to the base on SW1, E2,
E2S and R2 (and on MSI, CCI and PCI where it touches shared or 3D code).

- **Class sums balanced.** The batch is bucketed by class and half; the
  buckets run under a dynamic schedule, the largest first, and each sums its
  particles in batch order. The work arrays are sized by the thread team.
  Stress: `make_cavgs` on the E2 project at one thread and 20 times at 11
  threads, every class-average file byte-identical. Indication: the
  class-sum load balance at 11 threads went from 0.69 to 0.94. The optional
  per-particle FFT pass was not done (it needs a batch of padded images in
  memory).
- **One particle read per table iteration.** In shared memory the 2D table
  fills the particle buffer and the matcher's batches are copied from it;
  the copies run in parallel (`image_slab%get_many`, stress-tested in the
  stack I/O tester). The buffer holds raw images, so this saves the second
  read, not the second preparation, which costs more.
- **The particle buffer on macOS (defect).** The node's memory was read only
  from `/proc/meminfo`, so the buffer budget was zero and the buffer, 3D
  included, never ran on a Mac. It now reads `hw.memsize`. The cases then run
  the buffered paths and stay identical.
- **Search loops dynamic** in the 2D and 3D matchers; the 3D matcher's
  per-thread counters are sized by the team. Stress: E2S and a
  shared-memory `refine3D`, 20 times each at 11 threads, all to the end with
  finite results.
- **The 2D dense assignment without rescans.** A counter of the particles
  left, and only the classes whose head was the assigned particle advance.
  The tester keeps the old algorithm and compares assignments on synthetic
  tables with ties, and at one and at 11 threads.
- **Populations counted once** with `get_pops`, which takes an optional half.
- **The last two heap trampolines** are gone: array versions of the bounded
  and power samplers, with the old samplers kept in the tester as
  references.
- **The empty 2D partition.** Reproduced on the base (12 parts, 6 sampled
  particles: the table workers stopped). A part without sampled particles now
  writes an empty table and runs its matcher pass without batches, as 3D
  does. The same input runs to the end in dense and sparse mode.
- **For the maintainer:** other per-thread arrays still sized by the
  configured thread count rather than the team (`seed_bank` of the dense 3D
  fill, the 2D table worker's candidate buffers, `heap_vars`); `ref_assign`
  keeps one heap trampoline (`sample_likelihood_dist`).

### Phase 4, 2026-10-19: documents and close

- `doc/policies/2D/solve2D_policy.md`: the dense assignment's frontier, the
  empty partition, the search schedule and the particle buffer in shared
  memory, the balanced class sums, the benchmark reports.
- `doc/policies/3D/refine3D_policy.md`: the dense part file's state list,
  size check and errors after the parallel read; the load-before-evaluate
  check; the array samplers; the dynamic search schedule; the new and
  per-partition benchmark files.
- `doc/policies/particle_read_path_policy.md`: the 2D table fills the
  particle buffer in shared memory; the budget on macOS.
- The report: `speed2D_run_report.md`.
