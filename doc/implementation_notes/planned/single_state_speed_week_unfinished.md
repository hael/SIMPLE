# Single-state and multi-state speed week: the open work

**Status, 2026-10-20.** Most of the week's plan is done. What is left is
smaller than the plan: two profiling items, the last identical trims, one
larger restructuring (one worker process per iteration), and the measurements
on the large datasets that only the maintainer can run. Background, design,
the baseline and every finished item are in
[`../completed/single_state_speed_week_completed.md`](../completed/single_state_speed_week_completed.md).
The summary of days 1 and 2 is in
[`../completed/single_state_speed_week_report.md`](../completed/single_state_speed_week_report.md).

Terms used below. `solve3D` is de novo 3D map determination; `refine3D` is
3D refinement. The "table" is the probability table that scores every
particle against many reference projection directions (`prob`, and the
sparse neighbourhood version `prob_neigh`). A "part" is one partition of the
particles, run as its own process (`nparts`). PfCRT is the maintainer's
single-state benchmark (16,832 particles, box 300, baseline 4,567 s); msp1 is
the three-state benchmark (141,838 particles, baseline 10,415 s).
"Identical" means byte-identical outputs on the one-thread identity cases.

## What is done (for orientation)

Days 1 and 2 and the 2D table (run `prob_speed`); `prob_fast` retired; the
multi-state schedule and `ref_assign` without rescans (run `ms_speed`); the
identical parts of days 3 and 4 and the parallel table reads and compact table
record (run `speed_mac`). Details are in the completed note.

## What is open

1. **Score-loop profile at high resolution** (day 3, item 5). At shell 78 and
   box 160 the cost per reference is about twice what the shells and the box
   explain. The likely cause is the reference memo (about 4 GB at 5,000
   directions) streaming from memory for about 530 scattered references per
   particle. Profile one table part (perf counters, or two shell limits at the
   same number of directions) and decide whether to reorder the work, for
   example by blocking particles by neighbourhood so threads share references.
   Why open: `speed_mac` was limited to identical changes and did not profile;
   it only replaced the loop's per-call allocations by per-thread buffers.
2. **Single-precision weights in the score kernels** (day 3, item 6).
   `gen_euclid_crvec` and its bank version multiply single-precision data by a
   double-precision weight, which halves the vector width. This is a
   rounding-level change, so it is not identical: implement, measure, keep it
   only if the FSC check holds.
3. **Rest of day 4, item 3 (project reads).** Volume assembly and
   `commit_sigma2` still read the project file, because a copy of the
   master's in-memory project is not what they read. Decide whether to give
   them an in-memory project; the cost grows with particle count (a read of a
   751,729-particle project takes 1.3 s). Day 4's reprojection-model change
   covers the distributed master only; check whether the shared-memory
   strategy needs the same.
4. **The states' volume assembly in parallel** (multi-state item 5). Not an
   identical change, so `speed_mac` only added per-state timers and in-memory
   source halves. Several sums of a state's assembly depend on the thread
   count (radial FSC, the noise term of the ML regularisation, the smoothing
   weight of the filter's label update, the phase randomisation's per-thread
   random streams, FFTW plans). The nonuniform filter is module state with
   cache names that carry no state index. The multi-state peak-memory policy
   forbids holding several states' volumes in assembly. Parked; reopen only
   with a ruling that accepts non-identical output and a memory plan.
5. **One worker process per iteration** (day 5; the same item serves
   multi-state, estimated about 150 s on msp1). The table part and the search
   worker are two processes per partition and iteration, separated by the
   master's `ref_assign`; the second re-reads the project, reprojection model
   and particles and prepares them again. Design:
   - one worker per partition: fill the particle buffer, prepare references,
     fill and write its table part, signal "table written", wait for the
     master's assignment, then search and reconstruct with what it holds;
   - the master waits for the "table written" files, reads the tables, runs
     `ref_assign`, publishes the assignment by writing a temporary file and
     renaming it, then waits for the done files;
   - the worker waits in a sub-second `c_usleep` loop.
   Check first: `prob_align` and `refine3D` build their partitions through
   different `qenv%new` calls, and their first and last particles per part must
   be identical. Peak memory per worker rises (table, memoized references and
   particle buffer live together; the table part alone peaked at 0.6-0.7 GB on
   PfCRT), and the memory policy's budget must cover it. Implement behind the
   old path until the simulated tests (the distributed one included), a PfCRT
   run and the dataset with one stack per micrograph pass.
6. **Measurements only the maintainer can make.** One PfCRT run per block
   (days 3 and 4 together, then day 5, plus the one-stack-per-micrograph
   dataset), and one msp1 run after the merge, to compare state populations and
   per-state maps with the baseline. The estimates in the completed note
   (about 7,000 s for msp1 after `prob_speed`, 4,600-5,250 s with the schedule,
   4,000-4,700 s with the identical items) are not yet measured.
7. **A cap on refined references per particle across states** (138 of 7,500 in
   msp1's stage 5). Not ruled; it changes results and needs its own validation,
   as does a tighter `prob_athres` for the single-state `prob` table.

## Rulings that already bind this work

- Gains measured on the test machine are reported; compute-bound changes are
  timed there, PfCRT by the maintainer afterwards.
- The distributed simulated test is the `nparts=2` step in
  `simulated_workflow_1jxy` (no new CTest entry).
- `simulated_workflow_6vxx` uses root seed 20261003, `simulated_workflow_1jxy`
  20260923.
- No backwards compatibility: part files live for one iteration and need no
  format compatibility.
- The multi-state rulings (`prob_neigh` in `state` mode from 2,500 directions,
  the 2,500-direction cap, the stage-2 exit at 95% scanned) are implemented and
  stay.
- Parallel code is tested at 11-24 threads, repeatedly, not only at one thread.

## Acceptance checks

- One diff per item, compiled on its own; the fast gate passes.
- Identical items: byte-identical outputs of the two `simulated_workflow`
  suites (one thread, fixed seed) and of the one-thread identity cases.
- Changed results (items 2 and 4): every simulated test keeps passing its
  truth checks, with metrics recorded before and after; the PfCRT FSC 0.143
  stays within the run-to-run spread.
- Day 5: the distributed simulated step, `solve3D_addon`, a PfCRT run and the
  per-micrograph dataset pass; timing is read with
  `scripts/plot_refine3d_bench.py`.

## Not planned

Persistent workers across iterations (the existing backend starts one process
per task and does not help); routing `nparts=1` to the in-memory strategy
(the strategies differ in sampling and sigma paths); GPU offload of the table
fill (next release).
