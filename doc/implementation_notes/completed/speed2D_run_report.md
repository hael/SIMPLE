# 2D speed and finishing the identical speed changes: report of the run `speed2d_mac`

Report, 2026-10-19. The run carried out the plan in `speed2D_run.md` on the
maintainer's Mac, unattended, as one uncommitted diff. "Identical" means: the
outputs of the identity cases are byte for byte those of the build before
the run, at one thread with a fixed seed.

## What we learned

- **Every item of the plan is identical and kept.** Twenty changes in three
  areas (the probability tables' part files, the 2D benchmark reports, the
  2D speed items) leave every output of the one-thread identity cases
  unchanged. Two defects were found and fixed on the way, each identical on
  the cases.
- **The particle buffer never ran on a Mac.** The node's memory was read only
  from `/proc/meminfo`, so on macOS the buffer's budget was zero and the 3D
  particle buffer, and the new 2D use of it, silently stayed off. With the
  fix the identity cases run the buffered paths, and they stay identical.
- **A crash in the base on a sparse state.** In a multi-threaded three-state
  run a state of 9 particles got an all-zero FSC, and the post-processing
  indexed outside its spectrum (SIGBUS). The post-processing now cuts at the
  first shell when no shell reaches FSC=0.143.
- **Where 2D time goes** (E2, 1,500 particles, the Phase 2 build): at one
  thread the class sums take the largest part of the search step, then the
  search, then particle preparation; at 11 threads the class sums kept their
  threads busy only 63% of the time. Distributed runs spend about 18 s per
  run launching and polling jobs, whatever the thread count.
- **The 2D speed items work where the budget said** (an indication, one run
  each, alone on the Mac, shared memory at 11 threads): the class sums 7.9 to
  6.0 s (load balance 0.63 to 0.94), the search loop 6.0 to 5.4 s (0.85 to
  0.98), particle reading 0.8 to 0.3 s; the iterations' total 29.1 to
  26.6 s. The cases are small; the gains belong to the test machine to
  measure.
- **Thread stress testing found nothing**, and every parallel region the run
  added or changed has a repeated test at 11 threads.

## What changed

### The probability tables' part files (finishing the previous run)

- **Errors after the parallel read.** Reading a dense part file returns a
  status and stops nowhere inside the parallel loop, and builds no strings
  there. After the loop the first failing part stops the run, with its file
  name and the reason.
- **Short or unreadable part files stop the run.** The file size must be the
  size its header implies before anything past the header is read; every
  read checks its status.
- **The active states in the part file.** The dense part file holds the
  number and list of active states after the header (format word -3); the
  master compares them with its own, so a row's reference index is checked
  again.
- **Table types.** The format word is private to the dense table; the
  neighbourhood table refuses the dense row writer; the refined-entry
  offsets stay default integers behind a size guard.
- **Load before evaluate.** The shift search's evaluators take the reference
  index and stop unless the thread's loaded reference and half match.
- **Volume assembly's handed-on project** is copied after the reconstructors
  and volumes are freed, and the strategy frees its copy right after the
  sigma2 commit, its last use.
- **Dead code**: `gen_corr_grad_for_rot_8`, the per-thread buffers `kcorrs`
  and `pft_r1_8`, a release of the kept reprojection builder that could not
  run; stale comments. The planned C-alpha note uses `nu_filter_vols`.

### The 2D benchmark reports

- `REFINE2D` for every partition (it was partition 1 only).
- `REFINE2D_STRATEGY` times the whole iteration: setup, pre-step, matcher or
  scheduler, merge, class-average assembly, the master's project writes and
  rereads, sigma2, convergence, star file, JPEG. Its "accounted for" share is
  now about 99.9% (it covered only the sigma2 work).
- A `PROB_ALIGN2D` record for the 2D pre-step driver.
- Thread-seconds and load balance of the 2D search loop and of the class
  sums (also in the offline `CAVGSUMS`).
- `scripts/plot_refine3d_bench.py` builds the 2D critical path per iteration
  (`bench_iterations2D.csv`, `bench_profile2D.svg`, a report by stage).

### 2D speed

- **Balanced class sums.** The batch is bucketed by class and half; the
  buckets run under a dynamic schedule, the largest first, and each adds its
  particles in batch order, so the sums are the same at any thread count.
- **One particle read per table iteration.** In shared memory the 2D table
  fills the particle buffer and the matcher copies its batches from it, in
  parallel. The buffer holds raw images, so the second preparation remains.
- **Dynamic search loops** in the 2D and 3D matchers.
- **The 2D dense assignment without rescans.** Only the classes whose head was
  the assigned particle move on.
- **Class populations counted once** (`get_pops` with an optional half).
- **No heap trampolines left in the polar kernels**: array versions of the
  bounded and power samplers.
- **The empty 2D partition.** A distributed part with no sampled particle
  stopped the run (reproduced on the base: 12 parts, 6 sampled particles). It
  now writes an empty table and runs its matcher pass without batches, as 3D
  does.

### Fixes found on the way

- The post-processing cuts at the first shell when no shell reaches
  FSC=0.143 (it crashed).
- The node's physical memory on macOS (`hw.memsize`), so the particle buffer
  runs on a Mac.

## How identity was shown

- **Cases** (one thread, fixed seed, the CTest environment, made again on the
  base of this run): the simulated workflow of beta-galactosidase (SW1, with
  its two-part step) and of the spike (SW6); a three-state `refine3D` (MSI),
  correlation refinement with fractional update and trailing reconstruction
  (CCI), the PCG backend (PCI); and three new 2D cases on 1,500 particles
  simulated from the 1jxy map: `solve2D` with two parts (E2) and in shared
  memory (E2S), whose logs show all six stages and the three table paths,
  and `refine2D refine=prob` with two parts (R2). `solve3D_addon` (AD1) and
  `cont_refine3D_1jxy` (CR1) once at the start and once at the end, compared
  by their check values.
- **Determinism first**: every case twice on the base build, identical under
  the exclusions declared in Phase 0 (the previous run's list, plus the
  benchmark folder at R2's case root).
- **Every item built from its own snapshot**; the 3D and table items judged
  on SW1, MSI, CCI and PCI, the 2D items on SW1, E2, E2S and R2, items in
  shared code on all seven. Every comparison came out identical. SW6 was
  identical at the end of each phase; AD1 (44 check values) and CR1 (59) at
  the end of the run.
- **The part files** were compared through the previous run's lossless
  decoder (accepted by the maintainer), taught the new layout: every field
  both layouts hold, bit for bit; the state list is validated, not compared,
  since the base file lacks it.
- **Thread stress tests**, for every parallel region added or changed, at 11
  threads (`omp_get_num_procs()`), 20 times or more: the parallel part read
  (24 times, against one thread, bit for bit); the shift searches; the
  likelihood and sampler scores; the slab copy of the particle buffer; the
  dense 2D assignment; the populations; the class sums through `make_cavgs`
  (every class-average file byte-identical at 1 and 11 threads); the 2D and
  3D search loops, whose draws move with the thread, through E2S and a
  shared-memory `refine3D` (all runs to the end with finite results).
- **The new stops were checked** with a throwaway program outside the
  repository: a part cut by 4 bytes, a part with other active states, a
  particle in two parts each stop after the loop, naming the file.
- At each phase's end the clean Debug build passed the fast gate (15 of 15,
  14-15 s), and E2S at 11 threads and MSI with two parts of 5 threads ran to
  the end.

## The 2D time budget (Phase 2 build, alone on the Mac)

| run | threads | wall | 2D tables | search | outside the workers |
| --- | --- | --- | --- | --- | --- |
| E2S, shared memory | 1 | 144 s | 36.6 s | 102.5 s | 0.2 s |
| E2S | 4 | 44 s | 10.1 s | 31.5 s | 0.2 s |
| E2S | 11 | 30 s | 6.9 s | 21.0 s | 0.3 s |
| E2, two parts | 1 per part | 95 s | 18.9 s | 53.1 s | 18.1 s |
| E2, two parts | 5 per part | 51 s | 8.0 s | 20.9 s | 18.6 s |

Inside the shared-memory search step at one thread: class sums 43 s, search
38 s, particle preparation 19 s (the table's another 11 s), table fill 25 s;
reading about 1 s (page cache). Not in this run, for later: the distributed
launch and polling (the worker fold), and the matcher's particle buffers,
allocated every iteration with the team's size (2 s at 11 threads).

## Not done, and why

- **The per-particle FFT pass of the class sums** (optional in the plan): it
  needs a batch of padded images in memory; the balanced buckets already
  bring the load balance to 0.94.
- **The second particle preparation** in table iterations stays: the buffer
  holds raw images. Buffering prepared images would be a design change.
- **No negative unit tests of the new stops**: a stop ends the test program;
  the throwaway program checked them instead.

## Decisions for the maintainer

1. A three-state run gave a 9-particle state an all-zero FSC and a 2.6 A
   matching hand-off in volume assembly. The post-processing no longer
   crashes on it, but why such a state gets that FSC is not investigated.
2. Per-thread arrays still sized by the configured thread count rather than
   the team the loop runs on (the class of 43de9b340's fix): the dense 3D
   fill's `seed_bank`, the 2D table worker's candidate buffers, the polar
   kernels' `heap_vars`. Equal in production runs, where `nthr` sets both.
   The run sized by the team every array of a region it touched.
3. `ref_assign` (3D) keeps one heap trampoline (`sample_likelihood_dist`).
4. The particle buffer now runs on Macs (it was meant to); it holds up to
   half the node's memory per process.

## Proposed skill text (not applied)

- `simple-main-nu-filt`: `nu_filter_vols(vol_even_nu, vol_odd_nu,
  vol_src_even, vol_src_odd)`, with "the source pair must stay alive and
  unchanged from setup to the composition"; and "cleanup of module
  allocatables; the filter writes no scratch files".
- `simple-nonuniform-regularization` and its reference map: drop the cache
  file names and their clean-up check; "scores each candidate pair right after
  filtering it", "filters the source pair again per label".
- `simple-solve2d`: "In shared memory the 2D table fills the particle buffer
  and the matcher copies its batches from it; class sums are bucketed by
  (class, half) under a dynamic schedule, each in batch order."

## Files

Changed: `simple_eul_prob_tab.f90`, `_neigh.f90`, `_utils.f90`, `_tester.f90`,
`simple_eul_prob_tab2D.f90`, `_tester.f90`; `simple_polarft_calc.f90`,
`_core.f90`, `_corr.f90`, `simple_pftc_shsrch_grad.f90`,
`simple_pftc_inplane_tester.f90`; `simple_refine3D_strategy.f90`,
`simple_refine2D_strategy.f90`, `simple_strategy2D_matcher.f90`,
`simple_strategy3D_matcher.f90`, `simple_matcher_ptcl_io.f90`,
`simple_matcher_pftc_prep.f90`; `simple_commanders_rec_distr.f90`,
`simple_commanders_prob.f90`, `simple_commanders_volops.f90`;
`simple_classaverager.f90`, `simple_classaverager_restore.f90`;
`simple_oris.f90`, `simple_oris_getters.f90`, `simple_oris_tester.f90`;
`simple_image_slab.f90`, `simple_stack_io_tester.f90`, `simple_posix.c`;
`scripts/plot_refine3d_bench.py`; the policies `solve2D_policy.md`,
`refine3D_policy.md`, `particle_read_path_policy.md`; the planned note
`next_release/calpha_nu_evidence_mask.md`; this plan and report. None of the
files reserved for the run `test_truth` was touched. Files that run may also
edit: `simple_eul_prob_tab.f90` and `simple_refine3D_strategy.f90`.

## Evidence

The run directory `~/agent_runs/runs/speed2d_mac`: `journal.md` (every step,
decision and defect), `review/item_NN_*.diff` (one diff per item),
`full.diff`, `exclusions.tsv`, `phase0_determinism.txt`,
`phase0_resources.txt`, `notes/2D_time_budget.md`. Scratch `~/speed2d_mac`:
`out/<build>/<case>_<label>/compare.txt`, `logs/i*_summary*.txt`,
`logs/stress_*.txt`, `logs/repro_i24.txt`, the harness in `bin/`.

## Commit message

    Finish the identical speed changes; measure and speed up 2D

    Dense table part files: errors raised after the parallel read, size
    and active-state checks; the shift search checks its loaded reference;
    volume assembly frees its handed-on project early; dead code removed.
    2D benchmark reports for every partition and the whole iteration, with
    load balance, and a 2D critical path in plot_refine3d_bench.py.
    2D speed: balanced class sums, one particle read per table iteration,
    dynamic search loops, a rescan-free dense assignment, populations
    counted once, array samplers without trampolines, empty partitions
    handled. Byte-identical on the one-thread identity cases.

    Fixes: post-processing no longer crashes on a state without an FSC
    crossing; the particle buffer's memory budget works on macOS.
