# Single-state speed: the plan for the week of 2026-10-15

Implementation note, 2026-10-14. **Plan; nothing implemented.** Revised
2026-10-15 for the unattended run `prob_speed` on the Oracle Linux test
machine, which carries out days 1 and 2 and the same changes in the 2D table
(section "Execution").

## Goal and baseline

Make single-state `solve3D`/`refine3D` as fast as the CPU code can be for the
release; GPU work follows in the next release.

The baseline is the PfCRT run in `~/Processing/PfCRT/bench` (analysis in
`bench_summary/analysis.md`): 4,567 s, of which refinement 4,452 s, with 10
partitions of 8 threads, 16,832 particles, box 300 matched at 88-160. Every
change below is measured against it with `scripts/plot_refine3d_bench.py`.

The savings in this plan are estimates from that run's reports and from
reading the code, not measurements. Together they come to 1,350-2,050 s,
about 30-45% of the run.

| work | where the time is today (PfCRT) | estimated saving | result |
| --- | --- | --- | --- |
| job polling | 1.0-2.0 s per polled phase, two phases per iteration | 300-400 s | identical |
| dynamic thread scheduling in the table fill | 10-17% idle threads in the fill | 150-200 s | statistically the same |
| particle-frame shift in the `prob_neigh` score loop | score loop at high resolution, about 580 s | 400-500 s | changed as `prob_fast` already changes `prob` |
| exact trims in the kernels, the refinement and `ref_assign` | shift refinement about 460 s, `ref_assign` 0.6-1.6 s per iteration | 100-200 s | identical |
| nonuniform filter without the disk round trip | 330 s | 100-200 s | identical |
| reprojection model and project I/O on the master | 246 s plus several project reads per iteration | 50-150 s | identical |
| one worker process per iteration | a second process, project read, reference read, particle read and preparation per iteration | 250-400 s | identical |

"Statistically the same" means the same algorithm with a different
assignment of random draws to particles; multi-threaded runs were never
bit-reproducible (`simple_rnd.f90`: only a one-thread run is).

## Rulings of 2026-10-15

- **`prob_fast` is retired.** Its runs validated in 3D, and the same change
  will work in 2D. Its particle-frame paths become the only paths, in the
  `prob` and `prob_neigh` tables and in the 2D table; the key goes from the
  parameters, the parser and every user interface, and code reached only
  with `prob_fast=no` is deleted.
- **The 2D table gets the same changes** in the same run (section "The 2D
  table").
- **Gains measured on the test machine are reported.** These changes are
  compute-bound, so beta-galactosidase timings there are meaningful; PfCRT
  is measured by the maintainer afterwards.
- **The distributed simulated test** is a second `solve3D` step with
  `nparts=2` in `simulated_workflow_1jxy`, added as checks to that suite
  (no new CTest entry).

## Testing procedure

The work is tested in three blocks. Each block is checked on the simulated
high-level tests first, and only then on one PfCRT run:

| block | simulated tests | then |
| --- | --- | --- |
| days 1 and 2, the 2D table | after each item, with beta-galactosidase timings (the run `prob_speed`) | one PfCRT run for the block |
| days 3 and 4 | after each item | one PfCRT run for the block |
| day 5 | after the item | one PfCRT run, and the dataset with one stack per micrograph |

**The simulated tests** are the high-level CTest entries with a simulation
truth (`ctest -L highlevel`):

- `simulated_workflow_6vxx` and `simulated_workflow_1jxy`: `solve3D` on
  simulated movies, the final map checked against the truth (correlation and
  masked FSC). They run on one thread with `SIMPLE_SEED` set, so their
  outputs are reproducible: an identical change must leave their maps and
  orientations byte-identical.
- `solve3D_addon`: `solve3D` on simulated particles on eight threads, poses
  and map checked against the truth. This one exercises the dynamic
  scheduling.
- `cont_refine3D_1jxy`: the only one that runs distributed (`nparts`), but
  with `refine=cont`, not the probabilistic tables.

**The distributed gap.** `simulated_workflow` and `solve3D_addon` run
`solve3D` without `nparts`, so they take the shared-memory strategy, where
the table and the search run in the master process. The job polling (day 1),
the redundant segment write and the reprojection model (days 1 and 4) and the
worker fold (day 5) are only on the distributed path. `simulated_workflow_1jxy`
therefore gets a second `solve3D` step with `nparts=2` on the same project
(ruled 2026-10-15), whose final map must pass the same truth checks as the
shared-memory step.

**What passing means:**

- One diff per work item, reviewed and compiled on its own; the fast gate
  passes.
- **Identical changes:** byte-identical outputs of the two
  `simulated_workflow` suites before and after.
- **Changed results:** every simulated test keeps passing its truth checks,
  with the truth metrics recorded before and after. At the block's PfCRT run,
  the FSC 0.143 resolution stays within the run-to-run spread.
- **Timing:** the simulated tests are too small to time; the saving of a
  block is measured at its PfCRT run, read with the summary script, and
  recorded per block in the progress section of this note.

The particle I/O run, which edited the same matcher and prob files, finished
on 2026-10-15; its changes are merged before the run `prob_speed` starts.

## Day 1: measurement and the identical fixes on the master

1. **Benchmark layout and the missing report.** Compile the `BENCH/` change
   (one appended file per report family). Add timers and the `PROB_ALIGN`
   report to `exec_prob_align_neigh`, which writes none, so the `prob_neigh`
   pre-step can be split like the `prob` one.
2. **Job polling** (`src/utils/qsys/simple_qsys_ctrl.f90`, `schedule_jobs` and
   `schedule_array_jobs`).
   - The loop tests the done flags before refreshing them and then sleeps a
     whole second, so the master always waits 1 to 2 s after the last part
     has finished. Refresh before the test.
   - Poll every 100-200 ms with the existing `c_usleep` binding
     (`src/extlibs/unix/src/unix_unistd.F90`, already used with 100 ms in
     `simple_forked_process.f90`). Keep one second on batch queues, where the
     done files sit on a shared file system.
3. **Redundant project write.** The strategy writes the project segment
   before `prob_align` (`simple_refine3D_strategy.f90`, about line 1310)
   although nothing changed it since the read just above. Remove it.
4. **Dynamic scheduling in the table fill.** The single-state `prob` loop
   (`simple_eul_prob_tab.f90`, about line 281) and the `prob_neigh` loops
   (`simple_eul_prob_tab_neigh.f90`, about lines 188, 199, 402, 412) are
   static. The multi-state `prob_fast` loop is already dynamic, with the
   measured gain (load balance 0.90 to 0.98) in its comment. Each particle is
   processed whole by one thread, so only the random draw a particle gets
   changes.

Test: the simulated tests after each item. Expected on PfCRT: 450-600 s less,
measured at the block's PfCRT run after day 2.

## Day 2: the particle-frame shift in `prob_neigh`

At 5,000 directions the `prob_neigh` score loop is the largest single item.
For single-state runs it shifts the *reference* for every scored reference:
it rebuilds the shift matrix in double precision and runs one forward FFT per
shell, about 530 references per particle. The `prob` table with `prob_fast`
shifts the particle once (`memoize_shifted_ptcl`) and scores every reference
against the unshifted memoized references (`gen_likelihood_val_bank`). The
neighbourhood table never reads `prob_fast`.

- Port the same path to the neighbourhood fill: memoize the shifted particle
  after the seed, score with the bank kernel in the coarse and the fine
  passes, start the refinement from the converted shift
  (`shift_ptcl2ref`), and record the seed in `ref_assign_neigh`.
  Unconditionally: there is no key.
- Retire `prob_fast` (ruling of 2026-10-15): its paths in the `prob` table
  (the single-state particle-frame seed; the multi-state particle-frame grid
  and its quadratic step) become unconditional; the key leaves
  `simple_parameters.f90`, the parser, the `refine3D` and `solve3D`-family
  interfaces and the table part's benchmark lines; code reached only with
  `prob_fast=no` is deleted, and testers that compared the two paths test the
  remaining one.
- Estimated 3-5 times faster score loop: 400-500 s on PfCRT.
- Validation as for `prob_fast`: the same semantics (a particle-frame shift,
  the angular Nyquist bin included) already accepted for `prob`.

Test: the simulated tests, then the PfCRT run for days 1 and 2. Expected:
850-1,100 s less than the baseline.

## The 2D table: the same changes

The 2D table (`simple_eul_prob_tab2D.f90`) does per particle what the
neighbourhood table does: it estimates one shift at the previously assigned
class, scores every class with that shift through the shifted-reference
kernel (a shift matrix and one forward FFT per shell for every class), and
refines the shifts of a few classes around the best one. Both fill paths do
this, the plain table (`fill_tab_range`, used by `refine=prob` and the
`solve2D` seed pass) and `fill_tab_prob_snhc_range` (the `prob_snhc` stages
`solve2D` runs by default).

1. **Benchmark first.** The 2D table part's report has wall times but no
   thread-seconds split; give it the fill profile the 3D tables write (seed
   estimation, score loop, refinement, fill wall time, load balance).
2. **The particle-frame shift** in both fill paths: shift the particle once
   after the seed and score every class with the bank kernel; record the seed
   as the 3D table does.
3. **Dynamic scheduling** of the 2D fill loops (about lines 253, 305, 366 and
   376 today).

A 2D table scores every class for every particle (often 100-500 classes), so
the score loop should gain as much as in 3D; how much of `solve2D` that is,
the benchmark of item 1 tells.

## Day 3: exact trims in the table, the refinement and the assignment

1. **Shift refinement** (the `prob` stages; `src/main/pftc/simple_polarft_corr.f90`):
   every cost or gradient evaluation copies the reference to double precision,
   twice in the `gen_corr_grad_only_for_rot_8` path, and rebuilds the shift
   matrix with its sines and cosines. Copy the reference once per
   minimization; build the matrix from the shift without the second copy.
2. **The shift matrix range** (`simple_polarft_geom.f90`): the conversion to
   single precision covers the shells up to `interpklim` (about box_crop/2)
   instead of the searched `kfromto`. Convert only the searched shells.
3. **`ref_assign` for `prob`** (`simple_eul_prob_tab.f90`, about lines
   592-600): every assignment scans all particles (`any(ptcl_avail)`) and
   advances the head of every reference. Count the available particles and
   advance only the references whose head was the assigned particle, as the
   neighbourhood version already does.
4. **Neighbourhood table read** (`read_tabs_to_glob`): every part's file is
   opened and read twice, serially. Read it once, counting and filling in one
   pass.
5. **Score-loop profile at high resolution.** The cost per reference at shell
   78 and box 160 is about twice what the shells and the box explain. The
   likely cause is the reference memo (about 4 GB at 5,000 directions)
   streaming from memory for about 530 scattered references per particle.
   Profile one table part (perf counters, or two shell limits at the same
   number of directions) and decide whether to reorder the work, for example
   blocking particles by neighbourhood so threads share references.
6. **Single-precision weights in the score kernels** (`gen_euclid_crvec` and
   its bank version multiply single-precision data by a double-precision
   weight, which halves the vector width). Rounding-level change: implement,
   measure, keep only with the FSC check.

Test: the simulated tests after each item; items 1 to 4 must leave the
`simulated_workflow` outputs byte-identical.

## Day 4: the master's serial work

1. **Nonuniform filter** (`src/main/nu_filt/simple_nu_filter_state.f90`,
   `cache_filtered_vols`, and `simple_nu_filter_bank.f90`, `setup_nu_dmats`):
   every candidate low-pass (up to 8, cut at the FSC) writes the filtered even
   and odd volumes to disk as MRC; the objective reads them back, and the
   final composition by label (`simple_nu_filter_apply.f90`) reads them a
   third time.
   - Not by keeping them in memory: 2 x candidates x box^3 x 4 bytes, up to
     0.26 GB at box 160 but 1.7 GB at box 300 and 8.6 GB at box 512 in the
     final pass, about twice the filter's present peak.
   - Instead, without extra memory: compute the objective right after
     filtering each candidate (one pair alive at a time), and in the
     composition filter each candidate again from the Fourier transforms of
     the two input volumes (a multiply and two inverse FFTs per candidate)
     instead of reading it. Same numbers, no disk traffic.
   - First split the 6.6 s per iteration with timers; then try the candidate
     loop in parallel where memory allows.
2. **Reprojection model** (`materialize_reprojection_model`,
   `simple_refine3D_strategy.f90`, about lines 523-551): every iteration
   builds a new builder, reads the whole project and builds the toolboxes
   before projecting. Keep one reference builder in the strategy for the run;
   the projection itself is already parallel.
3. **Project reads.** `prob_align` and volume assembly each read the project
   again although the master holds it in memory (volume assembly reads the
   particle segment twice). Hand them the in-memory project. The cost grows
   with the particle count: small on PfCRT, large on a 750,000-particle
   project.

Test: the simulated tests, then the PfCRT run for days 3 and 4. Expected:
250-550 s less than the run after day 2.

## Day 5 and the next week: one worker process per iteration

The `prob_tab` part and the search worker are two processes per partition
and iteration, separated by the master's `ref_assign`. The second one reads
the project, the reprojection model and the particles again and prepares the
particles again. Both use the same sample, batch formula, shells, crop and
particle preparation.

- One worker per partition: fill the particle buffer, prepare the
  references, fill and write its table part, signal "table written", wait for
  the master's assignment, then search and reconstruct with the references and
  particles it already holds.
- The master: wait for the "table written" files, read the tables, run
  `ref_assign`, publish the assignment by writing a temporary file and
  renaming it, then wait for the done files.
- The worker waits with a sub-second `c_usleep` loop.
- Check before starting: `prob_align` and `refine3D` build their partitions
  through different `qenv%new` calls; their first and last particles per part
  must be identical.
- Peak memory per worker rises: the table, the memoized references and the
  particle buffer are live together (the table part alone peaked at 0.6-0.7
  GB on PfCRT); the memory policy's budget must cover it.
- Design and the partition check on day 5; implementation and validation the
  following week, behind the old path until the simulated tests (the
  distributed run among them), the PfCRT run and the per-micrograph dataset
  pass.

## Execution: the run `prob_speed`

Days 1 and 2 and the 2D table are carried out as an unattended run on the
Oracle Linux test machine (arun run `prob_speed`), in five phases. Nothing is
committed; the maintainer reviews the complete diff. Every phase records
itself in the progress section.

Common to every phase:

- **Data.** Beta-galactosidase (5,513 particles, box 256 at 1.275 Å per
  pixel, `mskdiam=180 pgrp=d2`), copied into a fresh directory under the
  run's scratch directory for each run; the originals are never written.
- **Execution settings.** `nparts=4 nthr=6` for distributed runs and
  `nthr=24` for shared-memory runs, identical in every phase. Timed runs use
  an out-of-tree Release build with tests built; nothing else heavy runs while
  a run is timed.
- **Timing (ruled 2026-10-15).** Gains are reported. Every timed comparison
  is the same case before and after, at least twice each way, read with
  `scripts/plot_refine3d_bench.py` (wall time, the critical path per
  iteration, the table parts' thread-seconds). A change is kept unless it is
  slower than the build without it beyond the run-to-run spread.
- **Validation.** Every `solve3D` run uses `nsample=2500` and starts from the
  one common `solve3D_cavgs` result with `cavg_ini_ext=yes`. Until Phase 2
  removes the key, every 3D case runs with `prob_fast=yes` on its command
  line, so that the baselines measure the path that becomes the only one.
- **Acceptance of a repeated case.**
  - `solve3D` (cases A and A-shm): the final FSC 0.143 resolution lies
    within the larger of the Phase 0 spread and 3 % of its Phase 0 value.
  - `solve2D` (cases E and E-shm): the last iteration's mean score is no
    lower than 0.97 times the lowest Phase 0 run, and the particle-weighted
    mean class resolution of the last iteration is no worse than 1.03 times
    the worst Phase 0 run.
  - **Simulated tests:** each by its own checks against its truth.
    `simulated_workflow_6vxx` and `simulated_workflow_1jxy` run on one thread
    with a fixed seed: after an identical change their maps and orientation
    documents must be byte-identical to the previous phase's (the
    distributed step included).
- **Gates.** The fast gate (run by `./compile_debug.sh`), the CTest unit
  suites of the touched areas, and `python3 scripts/check_test_registry.py .`
  whenever tests or their registration change.
- **Not in this run.** Days 3 to 5; the generated indexes under
  `doc/code_overview/`; skill files (proposed text goes into the report). The
  plan stays in `planned/`, since days 3 to 5 follow.

### Phase 0: baseline and the distributed simulated step

1. Build the tree as it is (`./compile_debug.sh`, with the fast gate) and the
   Release build for timed runs.
2. Add the second `solve3D` step with `nparts=2` to `simulated_workflow_1jxy`
   (ruling of 2026-10-15): from the same input project as the shared-memory
   step, in a directory of its own, with the same settings, its final map
   checked against the same truth. Note it in
   `test_environment_policy.md`; raise the CTest timeout of the entry if the
   suite needs it. This is the only source change of the phase.
3. Run the baseline cases on the Release build:
   - **S.** `solve3D_cavgs`, single state; its final project and maps are the
     common start.
   - **A.** `solve3D nsample=2500 cavg_ini_ext=yes prob_fast=yes` from S,
     distributed, twice.
   - **A-shm.** A's command line in shared memory, twice.
   - **E.** `solve2D` on the particles, distributed, twice.
   - **E-shm.** E's command line in shared memory, twice.
   For every case record the exact command line, the wall time, the final FSC
   0.5 and 0.143 resolutions (for E the score and class resolution of the
   acceptance) and the summary script's report, and keep its benchmark
   reports.
4. Run `simulated_workflow_1jxy` (with the new step), `simulated_workflow_6vxx`,
   `solve3D_addon` and `cont_refine3D_1jxy`; record their results and keep
   the maps and orientation documents of the two `simulated_workflow` suites
   for the byte comparisons.

Exit criteria: the build and the fast gate pass; the distributed step passes
its truth checks; every case is recorded with its command line and the
spreads of A, A-shm, E and E-shm; S is kept for the later phases.

### Phase 1: day 1

The four items of day 1, in this order: the benchmark report of the
`prob_neigh` pre-step, the job polling, the redundant project write, then the
dynamic scheduling of the 3D fill loops.

- After the first three (identical changes): the two `simulated_workflow`
  suites byte-identical to Phase 0.
- After the fourth: every simulated test passes its checks (the one-thread
  suites stay byte-identical, since one thread has nothing to reorder).
- Timed: A, A-shm and E twice each against Phase 0 (the polling also shortens
  the distributed 2D phases).

Exit criteria: the gates pass; the byte comparisons and the acceptance hold;
the timings and the critical-path change are recorded per item where the
spread allows, otherwise for the phase.

### Phase 2: day 2 and the retirement of `prob_fast`

The particle-frame shift in the neighbourhood fill, and the removal of the
`prob_fast` key with every path reached only with `prob_fast=no` (day 2).

- The simulated tests' outputs change: each passes its truth checks, and the
  truth metrics are recorded against Phase 1.
- The unit suites of the tables and the polar kernels (the table and in-plane
  testers) test the remaining paths.
- Timed: A and A-shm twice each against Phase 1, now without `prob_fast` on
  the command line.

Exit criteria: no `prob_fast` left in the sources, the interfaces or the
living documents; the gates pass; the acceptance holds; the timings and the
score-loop thread-seconds of the `prob_neigh` stages are recorded against
Phase 1.

### Phase 3: the 2D table

The three items of section "The 2D table": the benchmark profile first (and
a second timing baseline of E and E-shm with it if it changes the timers),
then the particle-frame shift in both fill paths, then the dynamic
scheduling. Extend the summary script where the 2D reports need it.

- The simulated tests pass their checks (their `solve2D` steps change).
- Timed: E and E-shm twice each against Phase 2.

Exit criteria: the gates pass; the acceptance of E and E-shm holds; the 2D
table's thread-seconds and the timings are recorded against Phase 2.

### Phase 4: documents and close

- The living documents: `refine3D_policy.md` (the probabilistic pre-step's
  benchmark reports, and the particle-frame shift as the only path in every
  table instead of the `prob_fast` paragraphs), `solve3D_chunks_policy.md`
  (no `prob_fast` to forward), `solve2D_policy.md` (the 2D table).
- Final repeats of A, A-shm, E and E-shm once each, and every simulated test.
- The progress section, then the report beside this plan, with a table of
  every case: wall time and resolution in Phase 0 and now.

Exit criteria: the gates pass; the report is written.

### File table

| file | why | phases |
| --- | --- | --- |
| `doc/implementation_notes/planned/single_state_speed_week.md` | this plan and its progress | each |
| `scripts/plot_refine3d_bench.py` | the summary script; fixes, and the 2D reports | each |
| `src/main/commanders/test/simple_commanders_test_highlevel_workflow.f90` | the distributed `solve3D` step of `simulated_workflow_1jxy` | 0 |
| `production/CMakeLists.txt` | the timeout and comment of `simulated_workflow_1jxy` | 0 |
| `doc/policies/test_environment_policy.md` | the distributed step | 0 4 |
| `src/utils/qsys/simple_qsys_ctrl.f90` | the job polling | 1 |
| `src/main/strategies/parallelization/simple_refine3D_strategy.f90` | the redundant project write | 1 |
| `src/main/commanders/simple/simple_commanders_prob.f90` | the `prob_neigh` pre-step report; `prob_fast` benchmark lines; the 2D table's profile | 1 2 3 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90` | dynamic scheduling; `prob_fast` retired | 1 2 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_neigh.f90` | dynamic scheduling; the particle-frame shift | 1 2 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_utils.f90` | shared fill and seed routines | 2 3 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab2D.f90` | the 2D table | 3 |
| `src/main/strategies/search/probabilistic/simple_strategy3D_prob.f90` | consumers of the recorded seed, if they need it | 2 |
| `src/main/strategies/search/probabilistic/simple_strategy2D_prob.f90` | consumers of the recorded seed, if they need it | 3 |
| `src/main/pftc/simple_pftc_shsrch_grad.f90` | `prob_fast` options made unconditional | 2 3 |
| `src/main/pftc/simple_pftc_inplane_tester.f90` | its `prob_fast` comparisons | 2 |
| `src/main/pftc/simple_polarft_corr.f90`, `src/main/pftc/simple_polarft_memo.f90` | bank kernels where a table needs a variant | 2 3 |
| `src/main/params/simple_parameters.f90`, `src/main/params/simple_parameters_parse.f90` | the `prob_fast` key removed | 2 |
| `src/main/ui/simple/simple_ui_refine3D.f90`, `src/main/ui/simple/simple_ui_solve3D.f90` | the `prob_fast` inputs removed | 2 |
| `doc/policies/3D/refine3D_policy.md` | the pre-step reports; the particle-frame shift as the only path | 1 2 4 |
| `doc/policies/3D/solve3D_chunks_policy.md` | no `prob_fast` to forward | 2 4 |
| `doc/policies/2D/solve2D_policy.md` | the 2D table | 3 4 |

## Not in this week

- Persistent workers that live across iterations: the existing persistent
  worker backend still starts one process per task and does not help.
- Routing `nparts=1` to the in-memory strategy: the two strategies differ in
  sampling and sigma paths; to be checked separately.
- Fewer refined references in `prob` (a tighter `prob_athres` or a cap):
  changes results and needs its own validation.
- GPU offload of the table fill: the next release.

## Progress

Nothing yet.
