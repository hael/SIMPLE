# Speed and architecture: before the release, and the next release

Proposal note, 2026-10-18. **Partly implemented: see "Status" below.** The maintainer
accepted the proposals on 2026-10-18, except two 2D defaults that stay as
they are (section 1.2 c). It collects the changes proposed after the run
`speed_mac`, split into what should be done before release 4 and what moves
to the next release. What was
done, and how it measures against the SIMPLE standards, is in
`identical_speed_run_review.md` (beside this note) and in the run's report,
`identical_speed_run_report.md`.

The frame, from the maintainer's rulings: release 4 is CPU-only and must be
as fast as possible; GPU work is the next release. Release 4 keeps no
backwards compatibility (except reading older project files), so changes to
user-visible names and meanings belong before it, not after. A high priority
(2026-10-18): **the 2D workflow as fast as possible, and consistent in look
and feel with the 3D workflow.**

Each proposal says what it is, why, whether results stay byte-identical
("identical"), and its size. Gains are estimates unless a measurement is
named. Items that change results need a ruling and a check with the truth
suite of the run `test_truth`.

## Status, 2026-10-20

Done in the run `structure2D_owner` (on the Oracle Linux test machine; plan
and report in `doc/refactoring_notes/completed/structure2D_and_project_owner*.md`):

- **The x86 check of the speed changes**: `aea5b1846`, `af13bc828` and
  `f078cf683` are byte-identical to their parents on the x86 build as on the
  Mac.
- **Section 1.5, the 2D structure aligned with 3D**: the 2D matcher searches
  and accumulates; the 2D strategy restores the class averages in shared memory
  and decides convergence in both modes; `refine2D` is a public `simple_exec`
  program. Byte-identical.
- **Section 1.3, one project owner on the master, and section 1.4**: the
  bodies of `prob_align`, `prob_align_neigh`, volume assembly, the sigma2
  commit, `prob_align2D` and class-average assembly are strategy routines on
  `(params, build)`; the commanders are thin wrappers; the typed entry points
  (`execute_from_project`, `execute_hand_on`) and the builder's `spproj_src`
  option are gone. The distributed `refine3D` master makes its reprojection
  model on its own builder and writes its particle field once before the
  workers and once per iteration; the distributed `refine2D` master writes
  only the segments class-average assembly changes. The moved steps no longer
  reseed the random generator (ruling): their check builds with the reseeds put
  back are byte-identical; `cfar` and `lp_est` reach the master one iteration
  earlier (ruling).

## Status, 2026-10-19

Done, all byte-identical at one thread unless said otherwise:

- **Section 1.1**: all eight items (run `speed2d_mac`, Phase 1;
  `doc/implementation_notes/completed/speed2D_run_report.md`).
- **Section 1.2 a, measuring**: done in `speed2d_mac` (Phase 2); the
  `PROB_ALIGN2D` report now also splits the master's reads and assignment
  (run `speed2d3d_mac`).
- **Section 1.2 b**: class sums balanced, one particle read per table
  iteration, dynamic search loops, the dense assignment without rescans,
  populations counted once, the last two heap trampolines (`speed2d_mac`);
  the compact dense 2D record, parallel part reads with a coverage check,
  per-particle counts in the sparse part files, and the sparse assignment's
  class lists sorted in parallel (`speed2d3d_mac`). Open: fewer project
  reads and writes on the master (with section 1.4).
- **Section 1.2 d, look and feel**: all items (`speed2d3d_mac`, Phase 2):
  3D-style log lines and one convergence verdict per iteration, the solve2D
  interface texts, class-average names through helpers, the shared sigma2
  helpers, the policy sections.
- **Section 1.2 f**: the empty 2D partition (`speed2d_mac`). Open:
  partition balancing by active particles (changes results).
- **Section 1.3**: the NU label optimisation and the master's dense 3D
  table stored by reference (`speed2d3d_mac`). The shift matrix fused into
  the alignment kernels was tried and is not identical (below). Open: one
  project owner on the master (done 2026-10-20, above), single-precision
  weights.
- **Also done** (`speed2d3d_mac`, Phase 4, changes uncached results by
  rounding, accepted): the 2D table prepares its particles from the same
  in-memory cropped image as the 2D matcher, and in shared memory the
  particle buffer keeps those cropped particles, so `cache=no` equals
  `cache=yes` byte for byte.

Measured by the maintainer, 2026-10-19: `solve2D` on msp1 (141,838
particles, 150 classes) is twice as fast overall on the merged code. PfCRT
ran well on the code before `f078cf683`.

Ruled 2026-10-19: the 2D code structure aligned with 3D moves from the next
release to before this one (section 1.5). The states' parallel assembly
(with the NU filter as an object, which existed for it) is left alone for
now (Part 2).

Two corrections to section 1.3, found by reading the code at `af13bc828`:

- The master's dense 3D table was already particle-major (each particle's
  references contiguous); the useful change was the opposite layout, each
  reference's distances contiguous, which is what was done.
- The correlation's shift-independent term cannot be computed once per
  minimisation without changing results: it is built from the shifted
  reference, and the shift matrix's modulus drifts by rounding along its
  recurrence. Only fusing the shift matrix into the kernels could be
  identical, and on this compiler it is not (the fused multiply-adds
  differ between the two contexts), so it was restored.

## Part 1. Before the release

### 1.1 Finish the current change

Small changes to the run's code, found by the review. None changes a result;
they are ready as Phase 7 of the run (`bash ~/agent_runs/runs/speed_mac/launch.sh`),
checked by the unit tests, and they produce `phase7.diff` to apply on top of
the run's diff.

1. Stop on a short or unreadable dense table part file (today the end of
   file only prints, and the master's table is no longer filled first).
2. Write and check the active state list in the dense part file.
3. Raise part-read errors after the parallel loop, not inside it.
4. Free volume assembly's handed-on project early (about 1 GB at 750,000
   particles at the master's peak).
5. Delete the dead code the review found (`gen_corr_grad_for_rot_8`, the
   buffers `kcorrs` and `pft_r1_8`, the kept builder's release in
   `distr_initialize`); fix stale comments.
6. Guard the "load the reference before evaluating" contract of the new
   shift-search evaluators.
7. Tidy the table types (format word private, `write_rows` refused in the
   neighbourhood type, 64-bit offsets).
8. Documents: the NU filter skills (needs approval) and the planned C-alpha
   note.

### 1.2 The 2D workflow: fast, and consistent with 3D

**Where 2D time goes.** Beta-galactosidase, 5,513 particles, 50 classes, on
the test machine: `solve2D` distributed (4 parts of 6 threads) went from 149
to 80 s with the speed-week changes, shared memory (24 threads) only from 99
to 95 s. The distributed gain was almost all job polling, which shared memory
never had. The table fill is now a small part of 2D: in the 17 iterations
with a table it is about 8 s of about 61 s. Most of the time is around it
(this split is inferred from the reports and the code, not measured line by
line):

- Particles are read and prepared twice per table iteration, once for the
  table and once for the search, serially in the shared-memory master
  (`prob_tab2D`, `simple_commanders_prob.f90` line 549; the matcher,
  `simple_strategy2D_matcher.f90` line 159). 3D reads them once per pass
  through its particle buffer.
- The class sums run in parallel over the 50 classes only, with a static
  schedule (`simple_classaverager_restore.f90`, lines 264-265), and each
  particle costs a forward FFT at the padded full box. At 24 threads the
  busiest thread carries at least 1.44 times the mean load.
- The search loop is statically scheduled (matcher, lines 163-165), while
  the cost per particle varies.
- References are prepared twice per table iteration; the table is written
  and read back in the same process; each of the six stages starts a fresh
  `refine2D` with new toolboxes and a new project read.
- The particle cache is off by default (`cache=no`), so particles are read
  and transformed at full box.

**a. Measure first** (identical, about 100 lines; before any 2D speed work):

- `REFINE2D` written for every partition, not only the first; timers for the
  probabilistic pre-step, scheduling, assembly, convergence and writes in
  `REFINE2D_STRATEGY` (today it times only the sigma2 work, so its
  "accounted for" share is meaningless); a `PROB_ALIGN2D` report like 3D's
  `PROB_ALIGN`; thread-seconds and load balance for the search loop and the
  class sums, as the table fill has; a critical-path summary for 2D in
  `scripts/plot_refine3d_bench.py`.
- Read the existing shared-memory reports of the 2D cases from the run
  `prob_speed` on the test machine; their particle-reading, class-sum and
  search lines settle the order of the items below before code is written.

**b. Identical speedups:**

| item | what | size |
| --- | --- | --- |
| class sums balanced | dynamic schedule over (class, half) pairs, largest first; optionally the per-particle padding and FFT in a loop over particles before the accumulation per class; each class still sums its particles in batch order, so results are identical at any thread count | 60-120 lines |
| one particle read per table iteration | in shared memory, fill the existing pass buffer in `prob_tab2D` and serve the search's batches from it; copy batches out in parallel | small to medium |
| search loop dynamic | `schedule(dynamic)` in the 2D matcher (and the 3D matcher, which is static too); identical at one thread, only the random draws move at more threads, as already accepted for the tables | 2 lines each |
| dense assignment without rescans | `do while( any(ptcl_avail) )` in the 2D table (`simple_eul_prob_tab2D.f90`, line 1000) becomes a counter, as in 3D; quadratic in the particle count today | small |
| populations counted once | one `get_pops` pass instead of `get_pop` three times per class in the reference preparation | small |
| last two heap trampolines | array samplers for the bounded and power samplers used by the first two `solve2D` stages (`snhc_smpl`) | about 120 lines |
| compact record for the dense 2D table | the 3D record (8 bytes per entry), parallel part reads, coverage check; per-particle counts in the sparse part files so the read stops scanning every index array twice | 300-400 lines |
| fewer project reads and writes on the master | `prob_align2D`, class-average assembly and the sigma commit work from the master's project; assembly writes only the segments it changes, not the whole project; done in the shape of section 1.4, not with new commander entry points | part of section 1.4 |

**c. Levers that change results** (ruled 2026-10-18):

- **The master draws the sample in stages 1 and 2**, as 3D does; today each
  worker draws within its own particles, so results depend on the number of
  parts. Accepted; checked with the truth suite.
- Not done (ruled): the particle cache stays off by default for `solve2D`,
  and its fractional-update default (`nsample` 200,000) stays as it is.

**d. Look and feel: what the user sees** (results unchanged):

- Log lines as in 3D: `CONVERGED: yes|no`, "% PARTICLES USED FOR UPDATING THE
  MODEL", the matching low-pass line, `STAGE n WITH LP x A`, and
  `SOLVE2D ENTERING/RETURNED FROM REFINE2D STAGE n`; `prob_align2D` prints
  one line per iteration instead of five; the stage table drops its dead
  `MAXCLSPOP` column. About 40 lines.
- Interface: correct the default texts of `center`, `lpstart` and `lpstop`;
  expose `nstages` and make `nsample` standard, as in `solve3D`; document or
  remove `refine=snhc_smpl` (the controller turns it into `prob_snhc` from
  stage 3); delete the unused `refine2D` program variable. About 30 lines.
- Output names through typed helpers, like `simple_refine3D_fnames.f90`,
  instead of the class-average file name built inline about 12 times. Small.
- Documents: give `solve2D_policy.md` the Defaults and Stage-table sections
  that `solve3D_policy.md` has, and correct its out-of-date sigma2 line.
- Code: the three canonical sigma2 helpers that the 2D strategy copies from
  3D go into one shared module (removes about 80 lines).

**e. Look and feel: names and meanings that change behaviour** (accepted
2026-10-18; before the release because release 4 breaks compatibility
anyway; the stream pool also runs `refine2D`, so stream is affected):

- `maxits` and `minits`: in 3D the length of a run, in 2D absolute
  iteration numbers; `solve2D` sets `minits=maxits` in stages 1 to 4.
- `fillin`: in 3D it samples never-updated particles every fifth iteration;
  in 2D it blocks convergence until every particle has a class.
- Mode names: `prob_snhc` (2D) against `prob_neigh_mode=snhc` (3D).

**f. Robustness:**

- A distributed 2D partition with no sampled particles crashes (a batch size
  divides by zero, and an empty table stops the run), where 3D handles it.
  Found by reading the code; reproduce, then fix. Then balance the 2D
  partitions by active particles, as 3D does. About 60-100 lines.

### 1.3 Shared architecture and 3D before the release

| item | what | gain | identical | size |
| --- | --- | --- | --- | --- |
| one project owner on the master (3D and 2D) | the master materialises its reprojection model on its own builder; the bodies of volume assembly, class-average assembly, `prob_align`, `prob_align2D` and the sigma commit become strategy or domain routines on `(params, build)`; the master reads only what workers write and writes the project once per iteration | about 3-6 s per iteration at 750,000 particles in 3D, similar in 2D; one to two fewer full projects at the peak; settles section 1.4 | yes | 500-900 lines in 3D, mostly moved; 2D smaller |
| faster label optimisation in the NU filter | done (2026-10-19): reuse the setup's smoothing at each level, skip all-zero lines and block the strided passes of the tent smoothing, packed voxel lists per colour in the Potts sweeps | the label optimisation was 8.7 of 14.9 s of the filter on the three-state identity case; likely 2-3 times on it | yes | about 150 lines |
| transposed master table, contiguous reads | done (2026-10-19): the master's dense 3D table stored by reference (it was particle-major already) so the assignment streams; consecutive particles read in one statement; refined entries placed directly | the assignment and reads on large multi-state tables | yes | about 150 lines |
| shift-independent work once per search | the correlation's shift-independent term once per minimisation; the shift matrix built inside the kernel loops. Tried 2026-10-19: neither is identical (see Status) | part of the refinement's memory traffic | no, rounding level | small to medium |
| single-precision weights in the score kernels | single-precision data times a single-precision weight | twice the vector width of the score loop | no (rounding level) | small; only if the truth suite is ready in time |

**Stretch, only if the items above land early: the worker fold.** One
worker per partition and iteration fills its table part, waits for the
master's assignment, then searches with the references and particles it
already holds (3D and 2D). Removes the second process start, project read,
reference build and particle preparation per iteration (about 9.5% on
PfCRT). Identical achievable; medium to large; it needs the job scheduler
split into schedule, poll and wait, and it is much easier after the project
owner. Otherwise it is the first item of the next release.

### 1.4 The commanders used as library routines

**What is there.** `prob_align`, volume assembly and the sigma commit gained
typed entry points beside `execute(cline)` (`execute_from_project`,
`execute_hand_on`), called by the distributed `refine3D` strategy with
in-memory projects.

**Why it should not stay.** A SIMPLE commander is a command object reached
through `execute(cline)`. Called from a strategy with typed arguments, three
commanders become library functions. The 2D path needs the same hand-offs
(section 1.2 b), so extending the pattern to class-average assembly and
`prob_align2D` would spread it. The repository already has the right shape:
the PCG assembly's body is `execute_rec3D_pcg_distributed_master(params,
build, ...)` in `strategies/parallelization`, called with the master's own
builder.

**Recommendation.** Keep the entry points through the current testing, then
remove them before the release as part of "one project owner on the master"
(section 1.3), done once for both 3D and 2D: the bodies become strategy
routines on the master's `(params, build)`, with the few settings that
differ today through the command line (threads, `combine_eo`,
`update_frac`, state) as explicit arguments, and `execute(cline)` becomes a
thin wrapper. There is then nothing to copy and nothing to hand back. The
quick alternative, if the pattern should not reach master even briefly, is
to revert the three hand-offs (about 4 s per iteration at 750,000 particles
given up).

### 1.5 Moved before the release (ruling of 2026-10-19)

**The 2D code structure aligned with 3D.** Class-average restoration and
the convergence decision move out of the 2D matcher (`refine2D_exec`) into
the 2D strategy, as 3D keeps reconstruction and convergence outside its
matcher; `refine2D` becomes a public program, as `refine3D` is. Results
unchanged; wide. It meets "one project owner on the master" in class-average
assembly, so the two are planned together: the 2D structure first, then the
project owner for 3D and 2D. The stream pool runs `refine2D` through its
commander and must behave as before (no edit of stream code expected; to be
checked).

## Part 2. The next release

| item | what | why later |
| --- | --- | --- |
| batched, GPU-ready kernels | pure kernels over explicit arrays (references, particle times CTF, shifts, rotations) scoring many references per call, first batched per particle on the CPU | the GPU release; FFT plans and precision change results |
| the worker fold (if not done) | as in section 1.3 | size and risk |
| persistent workers | worker processes that live for the run and exchange only model, sample and assignment; perhaps 5-15% per iteration on queued clusters | large; after the project owner and the worker fold |
| the NU filter as an object, the states assembled in parallel | filter state in an object, one per state, so that volume assembly can run states concurrently | left alone for now (ruling 2026-10-19); open when it returns: holding several states' volumes (the multi-state memory policy forbids it today) and the sums that depend on the thread count (radial FSC, ML noise term, label smoothing weight, per-thread random streams, threaded FFTW plans) |
| streaming assignment | keep only the best candidates per reference from each part; master table memory down by orders of magnitude on large multi-state runs | results change; needs a ruling; large |
| particle read-ahead in its own thread | read the next batch while the current one computes (today the read-ahead runs in the same thread) | medium; gain depends on the storage |
| table pre-step folded into the 2D search | in shared memory keep the table's prepared references and particle transforms for the search | medium to large; the one-read item of section 1.2 b takes most of the gain first |

## Part 3. Measurements that settle the order

- **PfCRT and msp1 on the merged code**, with the new split timers of volume
  assembly: the shares of tables, search, the NU filter and its label
  optimisation, the master's project handling and the assignment.
- **The 2D cases' existing shared-memory reports** from the run `prob_speed`
  on the test machine, then the 2D case at 1x24, 2x12 and 4x6 threads against
  shared memory, three runs each (its spread was 6%).
- **The master's peak memory** in a distributed `refine3D` at a large particle
  count.
- **Identity on the x86 build** of the test machine for the run's diff.

## Part 4. Suggested order

1. Section 1.1 now, while the first test jobs run on the run's diff.
2. The 2D measurements (section 1.2 a).
3. The identical 2D speedups and the look-and-feel items of 1.2 b, d and f,
   as one unattended run with one-thread identity tests (the 2D simulated
   workflow steps of the existing suites cover `solve2D`).
4. The 2D code structure aligned with 3D (section 1.5), then "one project
   owner on the master" for 3D and 2D (settles section 1.4).
5. The master-drawn sample of 1.2 c and the renamings of 1.2 e, checked
   with the truth suite.
6. The worker fold if time remains; everything in Part 2 after the release.

Items 1-3 and the NU label optimisation are done (Status).
