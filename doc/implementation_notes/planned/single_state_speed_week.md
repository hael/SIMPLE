# Single-state speed: the plan for the week of 2026-10-15

Implementation note, 2026-10-14. Revised 2026-10-15 for the unattended run
`prob_speed` on the Oracle Linux test machine, which carried out days 1 and 2
and the same changes in the 2D table (section "Execution"; the record is in
"Progress" and the summary in `single_state_speed_week_report.md`). Days 3 to
5 are still planned. Multi-state runs added 2026-10-15 (section
"Multi-state").

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
- **After the Phase 3 stop:** the 2D particle-frame shift stays;
  `simulated_workflow_6vxx` runs with its own root seed, 20261003, and
  `simulated_workflow_1jxy` keeps 20260923 (floors unchanged); the crash of
  `solve3D`'s map symmetrization on very small data is fixed in this run
  (Progress, Phase 3).

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
   neighbourhood version already does. Carried out with the multi-state
   `ref_assign` by the run `ms_speed`.
4. **Neighbourhood table read** (`read_tabs_to_glob`): every part's file is
   opened and read twice, serially. Read it once, counting and filling in one
   pass. Carried out with the parallel reads of the `prob` table by the run
   `ms_speed`.
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
| `src/main/commanders/test/simple_commanders_test_highlevel_workflow.f90` | the distributed `solve3D` step of `simulated_workflow_1jxy`; a root seed per suite, 20261003 for 6vxx (ruling of 2026-10-15, after the Phase 3 stop) | 0 3 |
| `production/CMakeLists.txt` | the timeout and comment of `simulated_workflow_1jxy` | 0 |
| `src/main/simple_final_rec.f90` | the final reconstruction of a distributed one-thread `solve3D` left `nthr` off its workers' command line, which they require; found by the distributed step | 0 |
| `doc/policies/test_environment_policy.md` | the distributed step | 0 4 |
| `src/utils/qsys/simple_qsys_ctrl.f90` | the job polling | 1 |
| `src/main/strategies/parallelization/simple_refine3D_strategy.f90` | the redundant project write | 1 |
| `src/main/commanders/simple/simple_commanders_prob.f90` | the `prob_neigh` pre-step report; `prob_fast` benchmark lines; the 2D table's profile | 1 2 3 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90` | dynamic scheduling; `prob_fast` retired; `record_seed_shift` moves to the shared routines, for the 2D table | 1 2 3 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_neigh.f90` | dynamic scheduling; the particle-frame shift; imports `record_seed_shift` from its new place | 1 2 3 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_tester.f90` | its test of the seed recording with and without `prob_fast` tests the remaining one; imports `record_seed_shift` from its new place | 2 3 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_utils.f90` | shared fill and seed routines | 2 3 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab2D.f90` | the 2D table | 3 |
| `src/main/solve/simple_solve3D_utils.f90` | the map symmetrization's high-pass is never finer than its low-pass, so the symmetry search always has a Fourier shell (the crash on small data; ruling of 2026-10-15) | 3 |
| `src/main/volume/simple_volpft_corrcalc.f90` | refuses a Fourier band without a shell above the origin with an error instead of returning undefined correlations | 3 |
| `src/main/strategies/search/probabilistic/simple_strategy3D_prob.f90` | consumers of the recorded seed, if they need it | 2 |
| `src/main/strategies/search/probabilistic/simple_strategy2D_prob.f90` | consumers of the recorded seed, if they need it | 3 |
| `src/main/pftc/simple_pftc_shsrch_grad.f90` | `prob_fast` options made unconditional | 2 3 |
| `src/main/pftc/simple_pftc_inplane_tester.f90` | its `prob_fast` comparisons | 2 |
| `src/main/pftc/simple_polarft_corr.f90`, `src/main/pftc/simple_polarft_memo.f90` | bank kernels where a table needs a variant | 2 3 |
| `src/main/params/simple_parameters.f90`, `src/main/params/simple_parameters_parse.f90` | the `prob_fast` key removed | 2 |
| `src/main/ui/simple/simple_ui_refine3D.f90`, `src/main/ui/simple/simple_ui_solve3D.f90` | the `prob_fast` inputs removed | 2 |
| `doc/policies/3D/refine3D_policy.md` | the pre-step reports; the particle-frame shift as the only path | 1 2 4 |
| `doc/policies/3D/solve3D_policy.md` | the Fourier band of the map symmetrization's axis search (the crash fixed in Phase 3) | 4 |
| `doc/policies/3D/solve3D_chunks_policy.md` | no `prob_fast` to forward | 2 4 |
| `doc/policies/2D/solve2D_policy.md` | the 2D table | 3 4 |

## Multi-state

Added 2026-10-15. Multi-state runs matter even more than single-state ones
for the release (maintainer, 2026-10-15). This section measures them the same
way and adds their items to the blocks of the week.

### Baseline

The baseline is the three-state `solve3D` run in `~/Processing/msp1/bench`
(analysis in `bench_summary/analysis.md`; its per-iteration reports converted
to `BENCH/`, the originals kept): 10,415 s for 141,838 particles, 30,012 of
them searched per iteration, three independent states, box 288 matched at
88-160, 5 partitions of 16 threads. The default multi-state schedule: `shc`
in stages 1 and 2, `prob` in stages 3 to 5 (2,500 directions in stage 5).
The run used the multi-state table as it was before its particle-frame grid.

| where the time goes | seconds | share of the run |
| --- | --- | --- |
| probability tables (slowest part) | 7,075 | 68% |
| of which stage 5, 2,500 directions x 3 states | 4,632 | 44% |
| search workers (slowest); 789 s of it in the `shc` stages | 1,369 | 13% |
| the master around the tables: reading them, `ref_assign` | 789 | 8% |
| volume assembly | 374 | 4% |
| final reconstruction | 294 | 3% |

What differs from single state:

- The table is 68% of the run (43% on PfCRT). Every particle is scored
  against every direction of every state (7,500 references in stage 5) at
  zero shift, because independent states share no shift frame.
- In stage 5, two thirds of the fill's thread time is the shift refinement:
  138 refined references per particle, each with a 5 x 5 grid over every
  in-plane rotation and then L-BFGS-B.
- The table is dense, 28 bytes per particle and reference in the files. In
  stage 5 that is 6.3 GB written by the parts and read serially by the master
  every iteration, and 5.4 GB held by the master (peak 8.0 GB).
- `ref_assign` labels the states and then assigns within each state on one
  thread, rescanning every reference and particle per assigned particle: 8.6
  s per iteration in stage 5.
- `shc` stops at the first better reference. As the maps settle it scans
  most of the 3,000 references: 1 s per part in iteration 1, 26-28 s at the
  end of stage 2.

### What the run `prob_speed` brings

Retiring `prob_fast` makes the multi-state table always build the coarse grid
from 25 shifted copies of the particle, made once per particle and shared by
all its refined references. Where the grid fit holds, a quadratic step
replaces L-BFGS-B. The fill is scheduled dynamically, and the job polling is
faster.

On the two-state bgal `refine3D` test of 14 October this cut the refinement
3.4 times, from 141 to 41 thread-seconds, and the table step from 34 to 15 s.
Applied to msp1, it is estimated to save about 3,250 s in the tables and
150-200 s in polling, taking the run from 10,415 s to about 7,000 s. No
multi-state run has measured it yet. The first measurement is msp1 after the
merge (maintainer), plus case M below.

### Rulings of 2026-10-15 (multi-state)

- **`prob_neigh` in its `state` mode from 2,500 directions.** A multi-state
  run switches from `prob` to `prob_neigh` (`state`) at the first stage whose
  directions go from 1,000 to 2,500. In the default schedule that is stage 5,
  which so far ran `prob` over every direction of every state.
- **Directions capped at 2,500** in multi-state runs. Stages 7 and 8 of the
  eight-stage table, 5,000 directions today, stay at 2,500.
- **`shc` leaves stage 2 once it scans more than 95%** of the search space on
  average, whatever the state overlap. The measure is the mean over the
  particles of the percentage of references the search scanned, which it
  stores per particle as `frac` and the log prints as "% SEARCH SPACE
  SCANNED".
- **`refine3D_states` is capped at 2,500 directions too** (ruled later on
  2026-10-15). Its `prob_neigh` stage ran at 5,000 directions
  (`STAGE2_NSPACE`). The ruling caps the directions only; its 500 subspace
  directions (`STAGE2_NSPACE_SUB`) stay.
- Still open: a cap on the refined references per particle across states
  (138 of 7,500 in msp1's stage 5).

These rulings change `solve3D`'s multi-state schedule and `refine3D_states`.
`solve3D`'s particle route and its class-average route share the mode policy,
so both switch at stage 5. The stage-2 exit concerns the particle route, the
only one with `shc` stages. The run `ms_speed` carries them out
(`multistate_speed_run.md`, beside this plan).

### Testing multi-state

The simulated tests are single-state and do not reach these paths.

- **Identical changes** are checked in the unit suite of the 3D probability
  table (`simple_eul_prob_tab_tester.f90`). The present routine stays in the
  tester as the reference. The new one must give the same assignments on
  synthetic dense tables of several states with a fixed seed.
- **Timing on the test machine**, beta-galactosidase, run like the
  single-state cases (before and after, twice each way, read with the summary
  script). The cases are defined in `multistate_speed_run.md`:
  - a three-state `solve3D` from the common three-state `solve3D_cavgs`
    start (`solve3D_cavgs nstates=3` gives the three-state particle
    assignment that `cavg_ini_ext=yes` needs);
  - a three-state `solve3D` from the particles, which runs the `shc` stages
    as msp1 did;
  - a three-state `refine3D` at 1,000 directions, which times the dense
    table, its reads and `ref_assign`.
  - The highlevel tests with several states, `flex_pca_blobs` (which runs
    `refine3D_states_cavgs`) and `solve3D_pyramid_blobs`, check the changed
    paths against their planted truths.
- **Acceptance.** The schedule change alters results, and beta-galactosidase
  has no real states: its cases show that the schedule runs, every state
  keeps its particles and a map, and the resolutions stay close to the
  baseline. The check that counts is msp1: the maintainer compares the state
  populations and per-state maps with the baseline run.

### Items

| work | where the time is (msp1) | estimated saving | result | block |
| --- | --- | --- | --- | --- |
| the multi-state schedule (the rulings) | stage 5: 194 s of tables and 31 s around them per iteration after `prob_speed`; stage 2: 888 s | 1,750-2,400 s | changed | run `ms_speed`, then an msp1 run |
| `ref_assign` for several states | stages 3 and 4: 2.9-3.0 s per iteration at 3,000 references | about 70 s, growing with the particle count | identical | run `ms_speed` |
| reading the table parts in parallel | stages 3 and 4: 8.1 s per iteration | 130-180 s | identical | run `ms_speed` |
| a compact record for the dense table | stages 3 and 4: 2.5 GB per iteration written and read, 2.2 GB on the master | 80-120 s, and two thirds of the master's table memory | identical | days 3 and 4 |
| the states' volume assembly in parallel | 374 s | 100-200 s | identical | days 3 and 4 |
| one worker process per iteration | about 4 s per `prob` iteration | about 150 s | identical | day 5 |

With the schedule, stage 5 leaves the dense table. The dense-table items
therefore act on stages 3 and 4, and on the single-state `prob` stages. On top
of the run `prob_speed` (about 7,000 s), the schedule would take msp1 to
4,600-5,250 s and the identical items to roughly 4,000-4,700 s, against
10,415 s. These are estimates:

- **Stage 5:** the `prob_neigh` share assumes the neighbourhood table costs
  a third to a half of `prob` at the same directions. On PfCRT, single state,
  it cost a third (11.4 against 37.3 s per iteration).
- **Stage 2:** the exit assumes the scan passes 95% where the `shc` timings
  level off, around iterations 20-25 of the run. msp1's log has the exact
  percentages.

1. **The multi-state schedule** (rulings of 2026-10-15; results change).
   - `simple_solve3D_controller.f90`, `set_refine3D_mode_policy`, the
     multi-state branch:
     - `prob_neigh` with `prob_neigh_mode=state` from the first stage whose
       `nspace` exceeds 1,000, derived from the `NSPACE` table rather than
       from a second stage constant.
     - `NSTAGES_INDEPENDENT` keeps the default stop after stage 5 but gets a
       definition of its own, since `PROB_NEIGH_REFINE_STAGE - 1` no longer
       describes it.
     - With `nspace` 2,500, `nspace_sub` comes out at 126, as in
       single-state runs.
   - `set_refine3D_stage_controls`: when `nstates > 1`, `nspace` is at most
     2,500.
   - `simple_refine3D_states_controller.f90`: `STAGE2_NSPACE` 5,000 becomes
     2,500 (ruling); `STAGE2_NSPACE_SUB` stays 500.
   - The stage-2 exit: for a multi-state run's stage 2, the convergence rule
     (`simple_convergence.f90`, its multi-state branch) declares convergence
     once the mean scanned exceeds 95%, without the state-overlap term. The
     present rule (overlap and scan together) still ends the stage earlier if
     it is met. Stage 1 and every other user of the rule stay as they are.
   - How the controller tells the rule: through an internal key on stage 2's
     `refine3D` command line, set together with `fracsrch=95`. The key is
     declared like `cohort_sampling` (parameters and parser, marked
     internal, no interface entry).
   - The strategy's benchmark report records the mean percentage scanned, so
     the exit shows in the summary.
   - Documents: `solve3D_policy.md`, sections 3 and 6 (stage 5 in
     `prob_neigh`, the cap, the stage-2 exit), and the controller's comments.
   - Tested with MS and then msp1, as above.
2. **`ref_assign` for several states** (`simple_eul_prob_tab.f90`,
   `assign_greedy_state_labels`, `assign_refs_for_state`).
   - The state labels come from a global greedy frontier. A particle leaves
     the frontier at its smallest distance over all references, so its label
     is the state of its best reference. Compute that per particle, in
     parallel.
   - The frontier and the per-particle minimum agree except when two
     references of different states tie exactly for a particle's minimum, for
     example a flat particle whose distances are all equal. The tester checks
     both the general case and a tie, and the tie rule is written down (the
     lowest reference index).
   - Within a state, loop over that state's references only. Count the
     available particles instead of `any(ptcl_avail)`. Advance only the heads
     that pointed at the assigned particle, as the single-state item of day 3
     does. The random draws and their order are unchanged, so the
     assignments are identical.
3. **Reading the table parts** (`exec_prob_align`, `read_tab_to_glob`): read
   the part files in parallel; each fills only its own particles' columns.
   Skip the full initialisation of the master table, since every entry is
   overwritten and the read already stops on an incomplete part.
4. **A compact record for the dense `prob` table**, in the part files and on
   the master.
   - Per particle, one row of (distance, in-plane index) per reference, 8
     bytes. The reference and particle indices follow from the order.
   - The refined references go in a short list of (reference, shift) per
     particle.
   - An unrefined candidate has no shift of its own: zero in multi-state,
     the recorded seed in single state, applied as now.
   - In msp1's stages 3 and 4: 2.5 GB becomes about 0.7 GB per iteration,
     and the master's table about 0.7 GB instead of 2.2 GB. A multi-state run
     at 2,500 directions in `prob`, as before the rulings, would have gone
     from 6.3 to 1.9 GB.
   - The same change serves the single-state `prob` table. The
     neighbourhood and 2D tables are sparse and keep their records.
   - The part files live for one iteration, so the format needs no
     compatibility.
5. **The states' volume assembly** (`simple_commanders_rec_distr.f90`, the
   state loop). Time each state first; then restore the states in parallel
   where memory allows. *Not an identical change* (2026-10-17, see the status
   at the end): only the per-state timers were added.
6. **The worker fold** of day 5 applies unchanged; the table part and the
   search worker read and prepare the same particles.

Day 3's exact trims to the shift refinement (its items 1 and 2) also reach
the multi-state refinement wherever the quadratic step falls back to
L-BFGS-B.

## Not in this week

- Persistent workers that live across iterations: the existing persistent
  worker backend still starts one process per task and does not help.
- Routing `nparts=1` to the in-memory strategy: the two strategies differ in
  sampling and sigma paths; to be checked separately.
- Fewer refined references in `prob` (a tighter `prob_athres` or a cap):
  changes results and needs its own validation.
- GPU offload of the table fill: the next release.

## Progress

### Phase 0, 2026-10-15: baseline and the distributed simulated step

Done on the test machine `dell` (Intel Core Ultra 9 285, 24 cores, Oracle
Linux 8.10). Paths below are under the run directory
`/home/elmlundho/agent_runs/prob_speed`.

**Files changed.**

- `src/main/commanders/test/simple_commanders_test_highlevel_workflow.f90`:
  the `solve3D` step of `simulated_workflow` moved into an internal routine
  `run_solve3D_step`; for 1jxy a second call (step 11) runs `solve3D` with
  `nparts=2` from the same `solve2D` project (`mkdir=yes` works on a copy, so
  the input is the one the first step read), in its own directory
  (`dir_exec=solve3D_distr`), with the suite's settings (one thread, the same
  point group, mask and seed), and checks its final map against the truth
  with the same thresholds. The suite fails if either map fails.
- `src/main/simple_final_rec.f90` (added to the file table for this phase):
  the distributed step found a bug in the base code. The final
  reconstruction at the end of `solve3D` put `nthr` on its child command line
  only when `nthr > 1`; with `nthr=1` and `nparts=2` its two
  `reconstruct3D` worker processes stopped with "Missing key on command line:
  nthr" (`reconstruct3D` declares `nthr` required, and a worker process
  checks its command line against that). It is now set also whenever
  `nparts > 1`. Runs with more than one thread
  (every timed case) and shared-memory runs are unchanged: the one-thread
  suites gave byte-identical outputs before and after the fix.
- `production/CMakeLists.txt`: a comment on the 1jxy entry. Its timeout of
  7,200 s stays: the suite took 1,167-1,269 s on the Release build with
  three such runs side by side.
- `doc/policies/test_environment_policy.md`: the 1jxy row of the table of
  long-running CTest entries names the distributed step.

**Builds.** The machine's compilers come from the login shell
(`gcc-toolset-15`: gfortran 15.2.1; system FFTW; SIMPLE's own OpenBLAS;
`/usr/bin/python3` is 3.6, so Python 3.12 is put first on `PATH`). The
maintainer's checkout is built the same way (Release, `BUILD_TESTS=ON`,
OpenMP, `USE_ARCHOPT=ON`).

- Debug build with the fast gate, `./compile_debug.sh` (run by
  `scratch/keep/scripts/run_debug_build.sh`): 15 of 15 gate suites passed,
  after the fix as before it (`scratch/logs/debug_build_phase0_2.log`).
  `python3 scripts/check_test_registry.py .` passes.
- Release build for the timed runs, kept as `scratch/keep/build_phase0`:
  `cmake -S repo -B scratch/keep/build_phase0 -DCMAKE_BUILD_TYPE=Release
  -DBUILD_TESTS=ON -DPython3_EXECUTABLE=/usr/bin/python3.12`, then
  `make -j24 && make install`. Every run sets `SIMPLE_PATH` to the build and
  puts its `bin/` and `scripts/` first on `PATH`; `SIMPLE_QSYS=local`.

**Command lines of the cases.** Each case runs in a fresh directory holding
a copy of its input project, through `scratch/keep/scripts/run_case.sh <build>
<case dir> <input project> <arguments>` (wall time and peak memory of the
master process by `/usr/bin/time -v`); the whole series is
`scratch/keep/scripts/run_phase0_cases.sh` (variables `B` for the build and
`PH` for the phase directory of the kept reports). The input of S, E and
E-shm is a copy of `/data2/bgal_addon/control/bgal_addon.simple`
(beta-galactosidase, 5,513 particles, with a 50-class 2D solution). A and
A-shm start from S's final project, read from `scratch/cases/S`, which is
recreated from `scratch/keep/S` before each of them (the project points at
S's maps by absolute path).

| case | command line (run in the case directory) |
| --- | --- |
| S | `simple_exec prg=solve3D_cavgs pgrp=d2 mskdiam=180 nthr=24 projfile=bgal_addon.simple mkdir=yes` |
| A | `simple_exec prg=solve3D pgrp=d2 mskdiam=180 nsample=2500 cavg_ini_ext=yes prob_fast=yes nparts=4 nthr=6 projfile=bgal_addon.simple mkdir=yes` |
| A-shm | `simple_exec prg=solve3D pgrp=d2 mskdiam=180 nsample=2500 cavg_ini_ext=yes prob_fast=yes nthr=24 projfile=bgal_addon.simple mkdir=yes` |
| E | `simple_exec prg=solve2D mskdiam=180 ncls=50 nparts=4 nthr=6 projfile=bgal_addon.simple mkdir=yes` |
| E-shm | `simple_exec prg=solve2D mskdiam=180 ncls=50 nthr=24 projfile=bgal_addon.simple mkdir=yes` |

E repeats the settings of the 2D solution in the control project
(`mskdiam=180 ncls=50`, other keys at their defaults).

**Results on the base build.** The cases ran one after another with nothing
else on the machine, except S, which ran beside the four one-thread
simulated suites (S is the common start, not a timed comparison). FSC is the
Fourier shell correlation between the half maps, its resolution read at 0.5
and 0.143; for `solve2D` the mean score is the last "SCORE [0,1] AVG" of the
log and the class resolution is the population-weighted mean `res` of the
final class segment (`scratch/keep/e_metrics.sh`).

| case | run 1 | run 2 | resolution or score, run 1 / run 2 |
| --- | --- | --- | --- |
| S | 939.3 s | | FSC 0.5 7.591 Å, FSC 0.143 6.800 Å |
| A | 658.7 s | 663.2 s | FSC 0.5 4.533 / 4.533 Å, FSC 0.143 4.030 / 4.030 Å |
| A-shm | 430.0 s | 437.8 s | FSC 0.5 5.100 / 4.663 Å, FSC 0.143 4.239 / 4.185 Å |
| E | 148.6 s | 149.5 s | score 0.387 / 0.387, class resolution 9.511 / 9.894 Å |
| E-shm | 101.9 s | 95.8 s | score 0.387 / 0.387, class resolution 9.789 / 9.799 Å |

Spreads of the wall time (difference of the two runs): A 4.5 s (0.7 %),
A-shm 7.8 s (1.8 %), E 0.9 s (0.6 %), E-shm 6.1 s (6 %). The acceptance of
later phases follows from these: A's FSC 0.143 resolution within 0.121 Å
(3 % of 4.030 Å; the two runs agree exactly) of 4.030 Å; A-shm's within
0.126 Å (3 % of their mean, larger than the spread of 0.054 Å) of 4.212 Å;
E and E-shm a mean score of at least 0.375 (0.97 x 0.387) and a class
resolution no worse than 10.19 Å (E, 1.03 x 9.894) and 10.09 Å (E-shm,
1.03 x 9.799).

Critical path of the `refine3D` iterations from the summary script
(seconds, run 1 / run 2): A spends 603.9 / 609.1 s in the strategy, of which
the probability tables 241.1 / 248.3, the pre-step outside the tables (launch,
polling, table reads, `ref_assign`) 74.6 / 70.7, the search 78.2 / 77.6,
scheduling outside the search 71.9 / 72.6 and volume assembly 79.5 / 81.9:
the master's time around the slowest part is 152.3 / 149.2 s (25 %). A-shm
spends 392.0 / 398.8 s, of which the tables 174.0 / 186.9, the search
83.1 / 81.7, volume assembly 76.2 / 71.9 and the master's time around them
9.4 / 9.6 s. The table parts' thread-seconds, summed over all parts and
iterations (the summary script does not yet total them; summed from the
`THREAD-SECONDS` blocks of the reports):

| case | `prob` score loop | `prob` refinement | `prob` sum | `prob_neigh` score loop | `prob_neigh` refinement | `prob_neigh` sum | mean load balance `prob` / `prob_neigh` |
| --- | --- | --- | --- | --- | --- | --- | --- |
| A run 1 | 1894.8 | 337.8 | 2277.3 | 837.5 | 184.4 | 1088.5 | 0.76 / 0.76 |
| A run 2 | 1975.7 | 328.7 | 2350.0 | 867.6 | 188.1 | 1124.6 | 0.75 / 0.77 |
| A-shm run 1 | 1438.2 | 259.4 | 1733.4 | 823.0 | 182.7 | 1072.9 | 0.76 / 0.78 |
| A-shm run 2 | 1665.6 | 285.8 | 1989.8 | 717.0 | 165.1 | 941.4 | 0.72 / 0.78 |

The `prob_neigh` pre-step writes no `PROB_ALIGN` report yet (27 records,
all from the `prob` iterations), as the plan says for day 1.

Benchmark reports, summaries (`plot_refine3d_bench.py <run dir> <out dir>`),
the final projects and the case logs: `scratch/keep/bench/phase0/<case>_r<n>/`
and `scratch/keep/bench/phase0/S/`; S itself is kept in `scratch/keep/S`.

**Simulated tests on the Release build** (`scratch/keep/scripts/run_sim.sh
<build> <dir> 1 test=simulated_workflow suite=<suite>`, which sets the
environment of the CTest registration; `solve3D_addon` and
`cont_refine3D_1jxy` through `ctest`):

| test | result |
| --- | --- |
| `simulated_workflow_1jxy`, shared memory | pass: docking correlation 0.9659, whole-volume correlation 0.8796, FSC 0.143 against the truth 6.24 Å |
| `simulated_workflow_1jxy`, `nparts=2` | pass: docking correlation 0.9658, whole-volume correlation 0.8854, FSC 0.143 against the truth 5.67 Å |
| `simulated_workflow_6vxx` | pass: docking correlation 0.9528, whole-volume correlation 0.8276, FSC 0.143 against the truth 5.31 Å |
| `solve3D_addon` | pass in 567.5 s: union map against the truth 0.9389, FSC 0.143 4.83 Å; cohort pose excess 0.73 degrees |
| `cont_refine3D_1jxy` | pass in 397.1 s: median rotation errors 0.91 / 0.23 / 0.91 degrees in stages a / b / c, fractions in the 8 Å basin 0.998 / 0.997 / 0.998 |

Two runs of the same build of each `simulated_workflow` suite gave
byte-identical maps, FSC files and orientation documents, the distributed
step included. Two files differ between any two run directories by design:
`solve3D_manifest.txt` (run identifier, digests of absolute paths) and the
header word of `sigma2_state.bin` at bytes 121-128 (the layout digest, which
hashes the absolute paths of the particle stacks). The comparison script
`scratch/keep/compare_sim.sh <kept A> <kept B>` therefore compares the
sigma2 file outside that word. Kept for the comparisons of the later phases
(by `scratch/keep/keep_sim.sh`): `scratch/keep/sim/phase0/1jxy/7_solve3D`,
`scratch/keep/sim/phase0/1jxy/8_solve3D_distr` and
`scratch/keep/sim/phase0/6vxx/7_solve3D` (the final maps, half maps, FSC,
sigma2 state and a text dump of the particle orientations), and the
metrics of the two CTest entries in `scratch/keep/sim/phase0/ctest/`.

**Exit criteria.**

- The build and the fast gate pass: Debug build and gate 15/15 (twice), the
  Release build, the test registry check.
- The distributed step passes its truth checks (table above), in two runs.
- Every case is recorded with its command line, wall time, resolutions and
  summary, and the spreads of A, A-shm, E and E-shm are given above.
- S is kept in `scratch/keep/S`; the Release build in
  `scratch/keep/build_phase0`.

**Notes for later phases.** The summary script warns, for every run, that one
or two report keys occur twice and that the last iteration has no strategy
record (the closing sigma pass); the totals are not affected. In shared
memory the stage at which `prob` gives way to `prob_neigh` moved by one
iteration between the two A-shm runs, which shifts thread-seconds between
the two tables; compare their sum as well as each table.

### Phase 1, 2026-10-15: day 1 (report, polling, project write, dynamic scheduling)

The plan's line references were checked against the code first and held
(the job loops in `simple_qsys_ctrl.f90`, the write before the distributed
pre-step in the `refine3D` strategy, the static fill loops of the two 3D
tables).

**Files changed.**

- `src/main/commanders/simple/simple_commanders_prob.f90`: the driver of the
  neighbourhood pre-step (`prob_align_neigh`) now times its parts (setup,
  sampling and project write, scheduling of the table parts, table reads,
  `ref_assign`, assignment write) and writes them as a report of its own,
  `PROB_ALIGN_NEIGH_BENCH.txt`, in the format of the `prob` driver's
  `PROB_ALIGN_BENCH.txt` (the report writer takes the family name).
- `src/utils/qsys/simple_qsys_ctrl.f90`: the two blocking job loops
  (`schedule_jobs`, `schedule_array_jobs`) refresh the done flags before they
  test them, and pause 0.1 s between polls on the local queue (the
  `c_usleep` binding of the unix library); on batch queues and the
  persistent worker backend, whose done files may sit on a shared file
  system, they keep the one-second pause. Before, the master noticed the
  last part 1 to 2 s after it had finished.
- `src/main/strategies/parallelization/simple_refine3D_strategy.f90`: the
  distributed strategy no longer writes the particle segment before the
  pre-step; the master had just read the project file, so the write put back
  what the file held. The flag that guarded that read (`have_oris`) was set
  to true on every path of the strategy's initialization, so the read was
  unconditional in effect; the flag is deleted and the read made
  unconditional (dead code).
- `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90`: the
  table fill with a shift search always runs with dynamic OpenMP scheduling
  (it did only with the multi-state grid keys); the two copies of the loop
  are one.
- `src/main/strategies/search/probabilistic/simple_eul_prob_tab_neigh.f90`:
  the four fill loops (the stochastic and the subspace fill, each with and
  without shifts) run with dynamic scheduling.
- `scripts/plot_refine3d_bench.py`: the report totals the table parts'
  thread-seconds per table (seed estimation, score loop, refinement, their
  sum, the fill regions' wall time; load balance and cost per score averaged),
  and lists the `prob_neigh` driver's report next to the `prob` one.
- `doc/policies/3D/refine3D_policy.md`: the new report, and a paragraph on
  the dynamic scheduling of every shifted fill (the sentence that tied it to
  the grid keys is gone).

**Builds and gates.** Release builds: `scratch/keep/build_phase1` (all four
items, kept for Phase 2) and, for the split by item, a build of the first
three items alone (from a copy of the tree with the two table files at the
base revision; deleted after the runs). Debug build with the fast gate
(`./compile_debug.sh`): 15 of 15, before and after the last (comment-only)
edit (`scratch/logs/debug_build_phase1.log`,
`scratch/logs/debug_build_phase1_final.log`). The test registration did not
change.

**Byte comparisons.** `simulated_workflow_1jxy` (both `solve3D` steps, the
distributed one included) and `simulated_workflow_6vxx` on the build of the
first three items and on the build of all four: every kept file is
byte-identical to Phase 0 (`scratch/keep/compare_sim.sh`), with the same truth
metrics (1jxy: whole-volume correlation 0.8796 shared memory, 0.8854 with two
parts; 6vxx: 0.8276). The one-thread suites have nothing to reorder, so the
dynamic scheduling leaves them identical too, as the plan expected. Kept for
Phase 2: `scratch/keep/sim/phase1/`.

**The eight-thread simulated tests** on `build_phase1`
(`scratch/keep/sim/phase1/ctest/`): `solve3D_addon` passes in 505 s (union
map against the truth 0.9411, FSC 0.143 4.93 Å; Phase 0 0.9389 and 4.83 Å:
the dynamic order gives the particles other random draws);
`cont_refine3D_1jxy` passes in 382 s (it runs `refine=cont`, which uses no
table; only its last stage, the continuation of a polar `refine3D_auto`,
moved slightly: symmetric pose error 1.47 to 1.46 degrees).

**Timings.** Each case ran with the Phase 0 command lines
(`scratch/keep/scripts/run_phase1_cases.sh`), interleaved over three builds:
the build of the first three items, the build of all four, and the Phase 0
build for a third before-run in the same session; nothing else ran on the
machine. Reports and summaries: `scratch/keep/bench/phase1/<build>/<case>/`
(`phase0`, `phase1a` for the first three items, `phase1` for all four).

| case | Phase 0 build | first three items | all four items |
| --- | --- | --- | --- |
| A (`solve3D`, 4 parts x 6 threads) | 658.7, 663.2, 656.9 s | 549.4, 540.2 s | 521.3, 522.5 s |
| A-shm (`solve3D`, 24 threads) | 430.0, 437.8, 442.1 s | 446.6, 427.9 s | 416.0, 462.6 s |
| E (`solve2D`, 4 parts x 6 threads) | 148.6, 149.5, 148.4 s | 86.9, 87.9 s | 86.5, 85.9 s |

- **Job polling and the project write (first three items):** A 659.6 s to
  544.8 s on average, 115 s (17 %) less; E 148.8 s to 87.4 s, 61 s (41 %)
  less. In A the master's time around the slowest part (the pre-step outside
  the tables plus the scheduling outside the search plus the assembly
  outside volume assembly) fell from 149-155 s to 27-28 s; the strategy's wall
  time from 602-609 s to 489-498 s. The project write is too small to
  separate from the polling within the spread (a few seconds of A's run in
  total). Shared memory starts no jobs, and A-shm did not change beyond its
  spread.
- **Dynamic scheduling (the fourth item):** A 544.8 s to 521.9 s, 23 s (4 %)
  less, beyond the spread of either build (9 s and 1 s). The load balance of
  the table parts (slowest thread over the mean, averaged over the parts)
  rose from 0.75 to 0.96 in the `prob` table and from 0.77-0.80 to 0.94 in
  the `prob_neigh` table; the slowest table part on A's critical path from
  248-249 s to 222-225 s, and the summed wall time of the parts' fill
  regions from 517-519 s to 423-434 s (`prob`) and from 248-264 s to 214-218 s
  (`prob_neigh`). The thread-seconds summed over the parts rose a little
  (`prob` 2,277-2,350 to 2,450-2,516, `prob_neigh` 1,089-1,194 to
  1,205-1,229), which the shorter fill regions more than pay for. In A-shm
  the load balance rose the same way (0.72-0.77 to 0.96, 0.77-0.78 to 0.93) and
  the tables on its critical path fell from 174-188 s to 149 s in the run with
  35 iterations; its wall time is dominated by how many iterations a run
  takes to converge (34-35 before, 35 and 38 after: the second run's 38
  iterations explain its 462.6 s), so its total did not separate from the
  spread. E has no 3D table and did not change with the fourth item.
- **Together:** A 659.6 s to 521.9 s, 138 s (21 %) less; E 61 s (41 %) less.

**Acceptance.** A: FSC 0.143 at 4.030 Å in every run (Phase 0: 4.030 Å; FSC
0.5 4.533 Å throughout). A-shm: 4.185 and 4.030 Å (all four items) and 4.185
and 4.239 Å (first three). The third before-run on the Phase 0 build gave
4.030 Å, so the Phase 0 spread over its three runs is 0.209 Å (4.030-4.239 Å,
four Fourier shells at this box), larger than 3 % of their mean of 4.151 Å;
every Phase 1 run lies within 0.209 Å of that mean. E: last-iteration mean
score 0.387 in every run (limit 0.375); class resolution 9.39-10.13 Å (limit
10.19 Å).

**Exit criteria.**

- The gates pass: Debug build and fast gate 15/15.
- The byte comparisons hold after the three identical items and after the
  dynamic scheduling, for both `simulated_workflow` suites including the
  distributed step; every simulated test passes its checks.
- The acceptance of A, A-shm and E holds (above).
- The timings and the critical-path change are recorded per item for A and
  E, and for the phase in A-shm.

### Phase 2, 2026-10-15: the particle-frame shift in `prob_neigh`, `prob_fast` retired

**Inventory before the edits.** Every reference to `prob_fast` in the
tree, and what becomes of it:

| where | what it does | becomes |
| --- | --- | --- |
| `simple_parameters.f90`, `simple_parameters_parse.f90` | the key and its registration | removed: a command line passing it is refused as an unknown key |
| `simple_ui_refine3D.f90`; `simple_ui_solve3D.f90` (`add_prob_fast_input`, on `solve3D`, `solve3D_chunks`, `solve3D_addon`) | the developer input | removed, with the helper routine |
| `simple_eul_prob_tab.f90`, `fill_tab_range` | multi-state: the particle-frame coarse grid and its quadratic step (`grid_ptcl`, `grid_fit` options of the shift search); single state: the particle-frame seed (one shifted particle memo, scored with the bank kernel; refinement started from the seed converted at each reference's rotation) | unconditional; the reference-frame seed of the single-state table (the seed handed unconverted to the score kernel and to the refinement) was reached only with `prob_fast=no` and is deleted |
| `simple_eul_prob_tab.f90`, `ref_assign` and `record_seed_shift` | the seed of an assigned, unrefined candidate: as it is (particle frame) with the key, rotated by the candidate's in-plane angle (`materialize_seed_shift`) without it | particle frame only; the rotating branch and the routine's switch are deleted (in the multi-state table no seed is estimated, so both branches recorded nothing there) |
| `simple_eul_prob_tab_tester.f90` | tests `record_seed_shift` with and without the key | tests the remaining recording |
| `simple_commanders_prob.f90`, the table part's report | the key, the shifted particle memos and (multi-state) the grid fits, only with the key | the memos and the fits are always written (the key line is gone) |
| `simple_pftc_shsrch_grad.f90`, `simple_pftc_inplane_tester.f90`, `simple_eul_prob_tab_utils.f90` | comments naming the key on the grid options, their tests and the profile counters | comments name the multi-state table instead; the `grid_fit` and `grid_ptcl` options stay, off by default, because every other shift search (the 3D and 2D searches, the neighbourhood table) runs without them, so their tests keep comparing the two settings |
| `doc/policies/3D/refine3D_policy.md`, `doc/policies/3D/solve3D_chunks_policy.md` | the key's description; its forwarding to add-on runs | rewritten without the key (Phase 4 finishes the living documents) |
| `doc/refactoring_notes/completed/prob_search_speedup*.md` | the history of the key | kept: dated history |

The neighbourhood table never read the key: its fill shifted the reference
for every scored reference. It gets the particle-frame path as its only
path (below).

**Files changed.**

- `src/main/strategies/search/probabilistic/simple_eul_prob_tab_neigh.f90`:
  the neighbourhood fill takes the particle-frame path of the `prob` table,
  unconditionally. Wherever it seeds (`shc` mode, and the `geom` and `state`
  modes of single-state runs, with a shift search), the seed of
  `estimate_shift_seed` is converted once to the particle frame at the
  rotation its search ended at; the particle is shifted once by it
  (`memoize_shifted_ptcl` into a per-thread memo, unless the seed is below
  `SHERRSQ`) and every reference of the coarse pass (the subspace
  representatives) and of the fine pass (the neighbourhood, or the
  stochastic draws) is scored against that memo with the bank kernel
  (`gen_likelihood_val_bank`, which includes the angular Nyquist bin); a
  refined reference starts from the seed converted at its own rotation
  (`shift_ptcl2ref`) and clamped into the search box; `ref_assign_neigh`
  records the seed of an assigned, unrefined candidate as it is
  (`record_seed_shift`, as the `prob` table does) instead of rotating it by
  the candidate's in-plane angle. The table's report counts the shifted
  particle memos. The seed table and the table file keep their layout; the
  only consumer of the recorded shift, `simple_strategy3D_prob.f90`, reads it
  as the particle shift for both tables and needed no change.
- `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90`: the
  `prob_fast` paths unconditional (the multi-state grid options on, the
  single-state particle-frame seed), the reference-frame seed deleted;
  `record_seed_shift` lost its switch and the rotating branch.
- `src/main/strategies/search/probabilistic/simple_eul_prob_tab_tester.f90`
  (added to the file table): tests the remaining seed recording.
- `src/main/params/simple_parameters.f90`, `simple_parameters_parse.f90`,
  `src/main/ui/simple/simple_ui_refine3D.f90`, `simple_ui_solve3D.f90`: the
  key and its inputs removed. `simple_exec prg=refine3D ... prob_fast=yes` now
  stops with "unknown command-line argument: prob_fast".
- `src/main/commanders/simple/simple_commanders_prob.f90`: the table part's
  report writes the shifted particle memos always and the grid fits for every
  multi-state table.
- Comments only: `simple_pftc_shsrch_grad.f90`,
  `simple_pftc_inplane_tester.f90`, `simple_eul_prob_tab_utils.f90`.
- `doc/policies/3D/refine3D_policy.md`: the key's section rewritten as the
  description of the particle-frame shifts of every probabilistic table, the
  neighbourhood table included; `doc/policies/3D/solve3D_chunks_policy.md`:
  the key's row removed.

`materialize_seed_shift` (the rotation of a reference-frame seed) stays in
`simple_eul_prob_tab_utils.f90` for the 2D table, which Phase 3 changes.

**Builds and gates.** Release build `scratch/keep/build_phase2` (kept for
Phase 3); Debug build with the fast gate 15 of 15
(`scratch/logs/debug_build_phase2.log`), which runs the unit suites of the
tables and the polar kernels (`unit_pftc_align2D3D`: the 3D probability
table, the continuous in-plane search with the grid options, the polar
correlation). A review of the diff by a second agent found no defect.

**Simulated tests** (outputs change, as expected; truth metrics against
Phase 1; kept in `scratch/keep/sim/phase2/`):

| test | Phase 1 | Phase 2 |
| --- | --- | --- |
| `simulated_workflow_1jxy`, shared memory | correlation 0.8796, FSC 0.143 6.24 Å | 0.8801, 6.24 Å |
| `simulated_workflow_1jxy`, `nparts=2` | 0.8854, 5.67 Å | 0.8839, 4.54 Å |
| `simulated_workflow_6vxx` | 0.8276, 5.31 Å | 0.8287, 4.80 Å |
| `solve3D_addon` | union map 0.9411, FSC 4.93 Å, cohort pose excess 1.01 degrees | 0.9378, 5.03 Å, 2.60 degrees (floors 0.90, 6.36 Å, 5 degrees) |
| `cont_refine3D_1jxy` | pass | pass (no table; its last stage moved within its floors) |

**Timings** (`scratch/keep/scripts/run_phase2_cases.sh`, interleaved; the
before runs on `scratch/keep/build_phase1` with `prob_fast=yes`, the after
runs on `build_phase2` without the key; reports in
`scratch/keep/bench/phase2/<build>/<case>/`):

| case | Phase 1 build | Phase 2 build |
| --- | --- | --- |
| A | 521.3, 522.5 (in Phase 1), 521.0 s | 481.8, 480.3 s |
| A-shm | 416.0, 462.6 (in Phase 1), 408.2, 411.7 s | 380.7, 372.1, 371.1 s |

- A: 521.6 s to 481.1 s, 41 s (8 %) less. The `prob_neigh` score loop fell
  from 926-948 to 319-340 thread-seconds summed over the parts (2.8 to 2.9
  times fewer), its fill regions from 214-219 s to 107-112 s; the
  `prob_neigh` table on the critical path from 6.3-6.5 to 3.6-3.7 s per
  iteration at 5,000 directions and from 3.0 to 2.1 s at 2,500; the seed
  estimation (70-72 against 74-76 thread-seconds) and the refinement
  (190-192 against 205-208) changed little.
- A-shm: 410.0 s (the two runs of this session) to 374.6 s, 35 s (9 %) less,
  beyond both spreads (3.5 s and 9.6 s); the `prob_neigh` score loop from
  727-817 to 212-234 thread-seconds, its fill regions from 43-48 s to
  20-22 s, the tables on the critical path from 150 s to 122-123 s.

**Acceptance.** A: FSC 0.143 4.030 Å in every run. A-shm: 4.080, 4.239 and
4.132 Å, within 0.209 Å (the Phase 0 spread) of the Phase 0 mean of 4.151 Å;
the two before runs gave 4.239 Å.

**Exit criteria.**

- No `prob_fast` is left in the sources, the interfaces or the living
  documents (`rg prob_fast src production scripts nice doc/policies` finds
  nothing; it remains in this plan and in the dated history under
  `doc/refactoring_notes/completed/`).
- The gates pass (Debug and fast gate 15/15).
- The acceptance holds for A and A-shm; every simulated test passes its
  truth checks.
- The timings and the score-loop thread-seconds of the `prob_neigh` stages
  are recorded against Phase 1 (above).

### Phase 3, 2026-10-15: the 2D table

**Which 2D fill path each `solve2D` stage uses** (confirmed in
`simple_solve2D_controller.f90`, `set_refine2D_stage_search_policy` and
`set_cline_refine2D_seed_pass`, and in the logs of the Phase 0 runs of E):
`solve2D` runs six stages by default. Stages 1 and 2 (`refine=snhc_smpl`,
before `PROBREFINE_STAGE` = 3) use no probability table. Stages 3 to 5 run
`refine=prob_snhc`, whose table part fills with `fill_tab_prob_snhc_range`
(a random subset of the classes per particle, shrinking with the iteration).
The last stage, 6, runs `refine=prob`, which fills with the plain
`fill_tab_range` (every class for every particle). The seed pass, one dense
`refine=prob` iteration over all particles before stage 3, runs only on a
seeded run (`cls_init=prev`) and also uses `fill_tab_range`. In E (5,513
particles, 50 classes) that is 10 iterations without a table, 15 with
`fill_tab_prob_snhc_range` and 2 with `fill_tab_range`, in both the
distributed and the shared-memory run (in shared memory the pre-step runs
the same table part in the master's process).

**Files changed.**

- `src/main/strategies/search/probabilistic/simple_eul_prob_tab2D.f90`:
  - the fill profile: both fills record the thread-seconds of the seed
    estimation, the class score loop and the refinement, the wall time of the
    parallel region, the classes scored or refined, the shifted particle
    memos built and the shift searches' counters (`get_fill_profile`,
    `record_fill_profile2D`);
  - the particle-frame shift in both fills: the shift estimated at the
    particle's previous class is converted once to the particle frame at the
    rotation its search ended at, the particle is shifted once by it (unless
    it is below `SHERRSQ`), every class is scored against that memo with the
    bank kernel, a refined class starts from the seed converted at its own
    rotation and clamped into the box, and the assignment records an
    unrefined candidate's seed as it is;
  - dynamic scheduling of the four particle loops.
- `src/main/strategies/search/probabilistic/simple_eul_prob_tab_utils.f90`:
  `record_seed_shift` moved here from `simple_eul_prob_tab.f90` (now shared by
  the 3D, neighbourhood and 2D tables); `materialize_seed_shift`, the rotation
  of a reference-frame seed by the candidate's in-plane angle, had no caller
  left and is deleted. The imports of `simple_eul_prob_tab.f90`,
  `simple_eul_prob_tab_neigh.f90` and `simple_eul_prob_tab_tester.f90` follow
  (rows of the file table extended to this phase).
- `src/main/commanders/simple/simple_commanders_prob.f90`: the 2D table part
  hands its fill profile to the report; the report of a table that scores a
  different set per particle also writes the refined classes.
- `scripts/plot_refine3d_bench.py`: the share of the 2D table parts in the
  `refine2D` iterations that have one (slowest part, its fill, load balance).
- `doc/policies/2D/solve2D_policy.md`: the particle-frame shift, the dynamic
  scheduling and the fill profile of the 2D table.

A second agent reviewed the diff and found no defect; it pointed out that the
refinement call passed `irot` both as the search's in-out argument and inside
the expression of its start (not standard-conforming), now a temporary, and a
dead assignment, deleted.

**Builds and gates.** Release builds `scratch/keep/build_phase3a` (the fill
profile only, from the snapshot `scratch/src_phase3a`) and
`scratch/keep/build_phase3` (all three items). Debug build with the fast gate:
15 of 15 (`scratch/logs/debug_build_phase3.log`), including the unit suite of
the tables (`unit_pftc_align2D3D`).

**Timings** (`scratch/keep/scripts/run_phase3_cases.sh`, interleaved;
reports in `scratch/keep/bench/phase3/<build>/<case>/`):

| case | Phase 2 build | fill profile only | all three items |
| --- | --- | --- | --- |
| E (`solve2D`, 4 parts x 6 threads) | 86.1 s | 86.7, 85.1, 86.7 s | 83.7, 80.6, 80.6 s |
| E-shm (`solve2D`, 24 threads) | 94.9 s | 96.0, 93.8, 97.3 s | 87.8, 91.4, 87.7 s |

- The profile costs nothing measurable (the profile-only build against the
  Phase 2 build).
- The particle-frame shift and the dynamic scheduling: E 86.2 s to 81.6 s,
  4.6 s (5 %) less; E-shm 95.7 s to 89.0 s, 6.7 s (7 %) less, both beyond the
  spreads (1.6 and 3.6 s). The 2D score loop fell from 158 to 57-60
  thread-seconds summed over the parts in E (2.7 times fewer), and from 154
  to 53 in E-shm; seed estimation (about 65) and refinement (about 70) did not
  change; the load balance of the fill rose from 0.75-0.81 to 0.98; the
  slowest table part on the `refine2D` iterations with a table fell from
  24.8-27.1 s to 18.4-19.3 s (E) and from 27.8-28.1 s to 20.7 s (E-shm), from
  41-42 % to 34 % of those iterations.

**Acceptance of E and E-shm.** Last-iteration mean score 0.387 in every run
(limit 0.375); class resolution 9.81-9.91 Å (E, limit 10.19 Å) and
9.67-9.80 Å (E-shm, limit 10.09 Å).

**Simulated tests on `build_phase3`** (`scratch/keep/sim/phase3/`):

| test | Phase 2 | Phase 3 |
| --- | --- | --- |
| `simulated_workflow_1jxy`, shared memory | correlation 0.8801, FSC 0.143 6.24 Å | pass: 0.8863, 4.62 Å |
| `simulated_workflow_1jxy`, `nparts=2` | 0.8839, 4.54 Å | pass: 0.8837, 4.71 Å |
| `simulated_workflow_6vxx` | 0.8287, 4.80 Å | **fails**: correlation 0.5145 (minimum 0.80), FSC 0.143 against the truth 62.4 Å (maximum 40 Å) |
| `solve3D_addon` | pass, union map 0.9378, cohort pose excess 2.60 degrees | pass: 0.9425, 4.58 degrees (floor 5) |
| `cont_refine3D_1jxy` | pass | pass |

**Why the run stopped here (answered by the rulings below).** `simulated_workflow_6vxx` fails, run twice on
the same build with the same result (one-thread runs are deterministic).
Its `solve2D` step (80 particles, 4 classes) ends as well as before (final
score 0.407 against 0.414, mean class resolution 47 Å against 64 Å), but the
`solve3D` that starts from its class averages does not find the map (its
half-map FSC 0.143 reaches only 9.2 Å). To tell a regression from a fragile
test, the suite ran with three other root seeds on scratch-only copies of the
sources (an environment variable overriding the suite's seed;
`scratch/keep/bench/phase3/diag/`), on the Phase 2 behaviour (the profile-only
source) and on Phase 3:

| root seed | Phase 2 behaviour | Phase 3 |
| --- | --- | --- |
| 20260923 (the suite's) | pass | fail (correlation 0.51) |
| 20261001 | crash in `solve3D` (44 particles picked) | fail (correlation 0.55) |
| 20261002 | crash in `solve3D` | fail (correlation 0.53) |
| 20261003 | pass (0.8769, 4.71 Å) | pass (0.8763, 4.62 Å) |

With other seeds the suite picks only 44 particles, and the Phase 2 code
then stops in the map symmetrization after `solve3D`'s stage 3 with "i out of
range; get_ori" (an orientation index of 0; a defect of the present code on
very small data, outside this run's files). The evidence points to a suite
whose `solve3D` step is fragile at this particle count rather than to a
defect of the 2D change, but the suite's own seed now fails, and passing it
would need a change to the test (its seed, its size or its floors) that
neither the plan nor a ruling authorizes. The maintainer decides: accept the
2D particle-frame shift and change `simulated_workflow_6vxx` (for example more
particles, or another seed, with its floors kept), or keep the 2D table's
reference-frame seed.

**Rulings of 2026-10-15 after the stop.** (1) The 2D particle-frame shift
stays; `simulated_workflow_6vxx` gets its own root seed, 20261003 (the seed
that passed on both the Phase 2 behaviour and Phase 3), and
`simulated_workflow_1jxy` keeps 20260923; floors and every other setting stay.
(2) The crash in `solve3D`'s map symmetrization on small data is fixed in
this run.

**The seed per suite.** `simple_commanders_test_highlevel_workflow.f90` now
sets the root seed where it chooses each suite's files and symmetry
(`SEED_6VXX = 20261003`, `SEED_1JXY = 20260923`), and the log line
">>> Deterministic workflow seed" prints the seed used.

**The crash: cause and fix.** Reproduced with seeds 20261001 and 20261002 on
the Phase 2 behaviour (`scratch/diag_build_before`; 44 particles picked).
With diagnostic prints in a scratch copy: the map symmetrization after
`solve3D`'s third stage ran the symmetry-axis search with the high-pass
limit `hp` = 124.15 Å (the default, half the box: 0.5 x 191 x 1.3 Å) and the
low-pass limit 124.8 Å (the planned limit of the symmetry-search stage, which
the stage ladder takes from the class Fourier ring correlations; with 44
particles they never reach the thresholds, so the early stages are planned at
124.8 Å). At the stage's box (100 pixels of 2.496 Å) that leaves no Fourier
shell between the two limits (index range 2 to 1), so all 501 grid scores of
`volpft_srch4symaxis` were undefined (NaN), no grid point beat the initial
best, the best index stayed 0, and reading orientation 0 stopped the run
("i out of range; get_ori"). The fix is where the band is chosen and where it
is used:

- `src/main/solve/simple_solve3D_utils.f90`, `symmetrize`: the search's
  high-pass is `max(hp, lpsym)`, never finer than its low-pass, so the band
  holds at least the shell of the low-pass. It is set with the low-pass for
  each symmetrization (the setting at the command line's construction is
  gone). Whenever the low-pass is finer than `hp`, the normal case (12-20 Å
  against about 124 Å), nothing changes. A first version opened the
  high-pass to an octave below the low-pass; at this box that put the origin
  (index 0) into the band and made every score about 0.9996, so it was
  replaced.
- `src/main/volume/simple_volpft_corrcalc.f90`, `new`: a band without a shell
  above the origin is refused with an error that names both limits, instead
  of yielding undefined correlations.

With the fix both seeds run through (band: index 1 to 1, no undefined score;
`scratch/keep/sim/phase3/fix2_before_*.log`); their maps still miss the truth
check, as the Phase 3 runs with these seeds did, which is the data, not the
crash. This belongs in the final report.

**Final Phase 3 build** (`scratch/keep/build_phase3`, rebuilt;
`scratch/logs/final_phase3.log`):

- Debug build with the fast gate: 15 of 15
  (`scratch/logs/debug_build_phase3_final.log`).
- `simulated_workflow_1jxy` (seed 20260923): passes; both `solve3D` steps
  byte-identical to the Phase 3 run before the rulings (its `solve3D` runs
  with `pgrp_start` = `pgrp` and does no symmetry search, so the fix is off
  its path).
- `simulated_workflow_6vxx` (seed 20261003): passes, whole-volume
  correlation 0.8763 (minimum 0.80), FSC 0.143 against the truth 4.62 Å
  (maximum 40 Å); its symmetrization low-pass is 20.8 Å, finer than `hp`, so
  the fix does not change it.
- `solve3D_addon` passes (union map against the truth 0.9451, FSC 0.143
  4.93 Å, cohort pose excess 4.52 degrees against its floor of 5; it was 4.58
  on the build before the rulings and 2.60 in Phase 2, so it is worth
  watching) and `cont_refine3D_1jxy` passes
  (`scratch/keep/sim/phase3/ctest/`).
- The E and E-shm timings above stand: the rulings touch neither the 2D
  table nor `solve2D`.

**Exit criteria.**

- The gates pass (Debug and fast gate 15/15 on the final tree).
- The simulated tests pass their checks: `simulated_workflow_1jxy` (both
  steps) and `simulated_workflow_6vxx` on the final build, `solve3D_addon` and
  `cont_refine3D_1jxy` (above).
- The acceptance of E and E-shm holds.
- The 2D table's thread-seconds and the timings are recorded against
  Phase 2.
- The crash of the map symmetrization is reproduced, its cause found and
  fixed in the owning modules; the reproduction runs through both seeds.

### Phase 4, 2026-10-15: documents and close

**Files changed.** `doc/policies/test_environment_policy.md` (the root seed
of each `simulated_workflow` suite, from the ruling after the Phase 3 stop);
`doc/policies/3D/solve3D_policy.md` (added to the file table: the Fourier
band of the map symmetrization's axis search and the refusal of an empty
band); this plan (its header and this entry); the report
`single_state_speed_week_report.md` beside it. `refine3D_policy.md`
(the pre-step reports, the particle-frame shift as the only path of every
table), `solve3D_chunks_policy.md` (no `prob_fast` to forward) and
`solve2D_policy.md` (the 2D table) were brought up to date in Phases 1 to 3,
which changed their subjects, and needed nothing more.

**Final repeats** on the final Release build (`scratch/keep/build_phase3`,
whose sources have not changed since its last rebuild in Phase 3), once
each, with the Phase 0 command lines and the 3D cases without `prob_fast`
(`scratch/keep/scripts/run_phase4_cases.sh`; reports in
`scratch/keep/bench/phase4/<case>/`), on an otherwise idle machine:

| case | Phase 0 (mean) | final | FSC 0.5 / 0.143, or score and class resolution |
| --- | --- | --- | --- |
| A | 659.6 s | 478.7 s | 4.533 / 4.030 Å |
| A-shm | 436.6 s | 386.0 s | 4.663 / 4.080 Å |
| E | 149.0 s | 80.0 s | score 0.387, 10.03 Å |
| E-shm | 98.8 s | 95.1 s | score 0.387, 9.84 Å |

In A the strategy's wall time is 427.4 s (Phase 0: 603.9-609.1 s) and the
master's time around the slowest part 28.2 s (Phase 0: 149-155 s); the table
parts' load balance is 0.97 (`prob`) and 0.90 (`prob_neigh`), and the
`prob_neigh` score loop 299 thread-seconds (Phase 0: 838-868). In E the 2D
score loop takes 58 thread-seconds at a load balance of 0.98. E-shm's single
run lies above the three Phase 3 runs of the same build (87.7-91.4 s).

**Acceptance.** A: 4.030 Å, as Phase 0. A-shm: 4.080 Å, within the Phase 0
spread (0.209 Å) of its mean (4.151 Å). E: score 0.387 (limit 0.375), class
resolution 10.03 Å (limit 10.19 Å). E-shm: 0.387, 9.84 Å (limit 10.09 Å).

**Every simulated test** on the final build (`scratch/keep/sim/phase4/`):
`simulated_workflow_1jxy` passes (whole-volume correlation 0.8863 and, with
two parts, 0.8837, as in Phase 3), `simulated_workflow_6vxx` passes (0.8763),
`solve3D_addon` passes (union map 0.9460, cohort pose excess 3.68 degrees
against its floor of 5) and `cont_refine3D_1jxy` passes.

**Gates.** Debug build with the fast gate: 15 of 15
(`scratch/logs/debug_build_phase4.log`). No script lost its executable bit
(`git diff --summary` shows no mode change).

**Exit criteria.** The gates pass; the report is written
(`single_state_speed_week_report.md`).

### Status after the run `speed_mac`, 2026-10-17

The run carried out the identical parts of days 3 and 4 and of the section
"Multi-state" on the maintainer's Mac (`identical_speed_run.md` and its
report). Every item below was byte-identical to the base on the one-thread
identity cases.

- Day 3, item 1 (the shift refinement): done. One double-precision copy of
  the reference per minimisation; the gradient-only path and `shift_ref`,
  both unused, deleted. In addition, the all-rotation objective that
  L-BFGS-B's last accepted point repeats is reused.
- Day 3, item 2 (the shift-matrix range): done.
- Day 3, item 3 (`ref_assign`): done earlier by the run `ms_speed`.
- Day 3, item 4 (the neighbourhood read): done; the part writes its
  per-particle counts and the first pass reads only those.
- Day 3, items 5 and 6 (the score-loop profile, single-precision weights):
  not in the run. The score loop's per-call allocations and internal
  distance functions (heap trampolines on this compiler) were replaced by
  per-thread buffers.
- Day 4, item 1 (the nonuniform filter): done. Timers, the objective right
  after each candidate, and the composition filtering again from the input
  pair (timed not slower); no scratch files.
- Day 4, item 2 (the reprojection model): done for the distributed master;
  the builder is kept for the run and later iterations read only the
  segments the model uses.
- Day 4, item 3 (the project reads): `prob_align` gets a copy of the
  master's project (a copy of a 751,729-particle project takes 0.28 s, a read
  1.3 s); volume assembly's trailing reread and the one-record read of the
  matching low-pass are gone. Volume assembly and `commit_sigma2` still read
  the file: a copy of the master's project is not what they read.
- Multi-state, items 3 and 4 (the parallel reads, the compact record): done.
- Multi-state, item 5 (the states' assembly in parallel): not identical.
  Several sums of a state's assembly depend on the thread count (the radial
  FSC, the noise term of the ML regularisation, the smoothing weight of the
  filter's label update, the phase randomisation's per-thread random
  streams, FFTW plans made for the thread count), the filter is module state
  with cache names that carried no state index, and the policy on
  multi-state peak memory forbids holding several states' volumes in
  assembly. Per-state timers and the in-memory source halves were added
  instead.

