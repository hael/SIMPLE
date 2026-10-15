# Single-state speed: the plan for the week of 2026-10-15

Implementation note, 2026-10-14. **Plan; nothing implemented.**

## Goal and baseline

Make single-state `solve3D`/`refine3D` as fast as the CPU code can be for the
release; GPU work follows in the next release. `prob_fast=yes` is the default
from 2026-10-15.

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

## Testing procedure

The work is tested in three blocks. Each block is checked on the simulated
high-level tests first, and only then on one PfCRT run:

| block | simulated tests | then |
| --- | --- | --- |
| days 1 and 2 | after each item | one PfCRT run for the block |
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

**A gap: none of them runs the probabilistic tables distributed.**
`simulated_workflow` and `solve3D_addon` run `solve3D` without `nparts`, so
they take the shared-memory strategy, where the table and the search run in
the master process. The job polling (day 1), the redundant segment write and
the reprojection model (days 1 and 4) and the worker fold (day 5) are only on
the distributed path. Before day 1, give one simulated `solve3D` a distributed
run with `nparts=2`. This could be a second `solve3D` step in
`simulated_workflow_1jxy`, added as checks to that suite rather than as a new
CTest entry, or until that is decided a run by hand on the project the suite
leaves behind. Its final map must pass the same truth checks as the
shared-memory run.

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

The particle I/O run on the Dell edits the same matcher and prob files
(`simple_strategy3D_matcher.f90`, `simple_commanders_prob.f90`,
`simple_matcher_ptcl_io.f90`). Its remaining phases land before the items
that touch those files, or the items are rebased on them.

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
- Under `prob_fast` (the default), so `prob_fast=no` keeps today's path for
  comparison.
- Estimated 3-5 times faster score loop: 400-500 s on PfCRT.
- Validation as for `prob_fast`: the same semantics (a particle-frame shift,
  the angular Nyquist bin included) already accepted for `prob`.

Test: the simulated tests, then the PfCRT run for days 1 and 2. Expected:
850-1,100 s less than the baseline.

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
