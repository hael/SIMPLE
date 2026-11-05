# Speed and architecture: before the release, and the next release

Proposal note, 2026-10-18. **Partly implemented: see "Status" below; the architecture review's
main item, the batch pipeline of probabilistic alignment, was implemented by the run
`prob_pipeline` (2026-11-04/05).** The maintainer
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

## Architecture review, 2026-11-04

Source review at `69a623cb8`. This section supersedes the remaining-work order
below where it conflicts with current code. The older measurements and release
estimates remain historical evidence, not predictions for the current search.
The architecture described here was proposed on 2026-11-04 and implemented by
the run `prob_pipeline` on 2026-11-04/05 (`prob_align_batch_pipeline.md`, its
Progress section holds the evidence and the measurements); the review's text
below is kept as the record of what was found.

### What remains, and what has changed

The polar shift-search refactor has already removed unused evidence work,
reused map geometry/plans/workspaces and avoided physical-width calculations
for candidates that cannot enter the polish shortlist. Its retained scope and
validation are in the existing
[polar shift-search record](../../refactoring_notes/completed/polar_shift_search_performance_refactoring.md).
Further gains should address the workflow's scientific-object lifetimes before
adding more numerical shortcuts.

The main remaining architectural cost was **probabilistic alignment implemented
as a table pre-step followed by a separate generic matcher pass** (resolved by
the batch pipeline, 2026-11-04; the description that follows is the state
reviewed):

- Shared-memory strategies still enter commanders that create another toolbox,
  publish/reload project fields, and serialize tables and assignments inside the
  same process. Distributed execution launches table and matcher workers
  separately. See `simple_prob_align_strategy.f90`,
  `simple_commanders_prob.f90` and the two refine strategies.
- Table workers prepare references/noise and particle batches, then destroy
  their scientific state. The matcher loads/prepares scientific state again and
  rebuilds particle batches. Existing canonical particle buffers save some reads
  and full-box preparation, not every polar transform.
- The 3D probabilistic matcher already skips discrete exhaustive scoring and
  skips reference memoization unless continuous in-plane polish needs it. Its
  selected-pose polish, sigma2 residual calculation and reconstruction are real
  consumers, not an accidental second exhaustive search. Optimize their shared
  preparation rather than deleting these operations.
- Dedicated particle-I/O prefetch is already implemented in
  `simple_stack_read_session`; the separate read-ahead thread in Part 2 is stale.
  The existing `_worker` backend retains a script dispatcher, not scientific
  objects inside the SIMPLE child processes.

Crucially, assignment is now **independent for each particle's candidate row**.
Dense 3D first chooses the state from that row's closest reference, then draws
within that state's block. Dense/sparse 2D and neighbourhood 3D likewise draw
from each particle's own eligible candidates. No reference/class population
balancing requires a global table (`ref_assign` in `simple_eul_prob_tab`,
`simple_eul_prob_tab2D` and `simple_eul_prob_tab_neigh`). The sampled `n_eff`
readout is diagnostic, not a feedback dependency between rows. Consequently,
bounded streaming of **complete rows** can preserve the posterior; the older
Part 2 proposal to retain only best candidates is a different, lossy proposal
and is not required for this architecture.

### Proposed execution contract

Make probabilistic alignment a complete strategy-owned pipeline on the caller's
`(params, build)`, with explicit begin/batch/end lifetimes. Reuse existing
preparation, table evaluation, assignment, polish and sigma2 domain routines;
do not create a second implementation of those algorithms. The commander stays
a command wrapper. Table/assignment files become transport or checkpoint
formats, rather than the mandatory interface between in-process operations.

The target batch path is:

```text
immutable iteration model + incoming particle metadata
  -> read/canonical observation -> prepare alignment batch
  -> evaluate complete candidate rows -> draw assignments + insertion extras
  -> polish selected poses when enabled -> sigma2 residuals when requested
  -> compact final decisions + canonical reconstruction observations
```

Ownership and lifetime:

| Lifetime | Owned state and boundary |
| --- | --- |
| Refine command | Partition identity, cache keys and small control metadata. No large alignment allocation is retained across the 3D reconstruction barrier merely because the process remains resident. |
| Alignment phase of one iteration | Immutable reference/half/state model, noise generation, band/grid, masks, sampled indices, incoming orientations/shifts and temperatures; compatible FFT plans and scratch capacity. One project owner draws the sample. Workers consume this snapshot and return changed records. |
| Bounded batch | Canonical observations, prepared polar particles/CTF terms, complete candidate rows and thread-local scratch. Keep them alive until selected-pose polish and sigma2 consumers finish, then recycle. |
| Across alignment/reconstruction barrier | Final poses, shift state/widths, frozen weighted-insertion poses/weights, noise contributions and the existing budgeted canonical observation buffer. No alignment PFTC/reference memo survives into 3D reconstruction. |

Represent the required consumers explicitly: candidate scoring, selected-pose
polish, sigma2 and reconstruction preparation. A path that only commits a
decision needs no polar particle preparation; a sigma2 or continuous-polish
consumer still does. Decide this from actual consumers, not a substring of the
mode name. A lightweight probabilistic finalizer can replace per-particle
allocation of a general search object once the shared pipeline works.

Preserve the distinctions between full-box input, canonical untapered cropped
reconstruction observations and masked/normalized alignment PFTs. Reusing one
as another requires a proven preparation contract, including CTF, band, half,
noise, centering and coordinate frame. Keep incoming metadata immutable until
each decision is complete; preserve the increment prior, shift-frame conversions,
invalid/empty-row fallbacks, temperature, candidate order and topk weights.
Continuous polish must not silently reweight insertion extras. The separate
Cartesian `pose_cont` pass is outside this initial fusion.

The memory target is references plus thread scratch plus bounded batch
observations/candidates plus compact particle results: no mandatory global
`Nparticles * Nreferences` table and no all-particle PFT cache. In 3D, release
alignment objects before allocating the existing one-state/one-half
reconstructor. Do not fuse all-state reconstruction into the alignment loop;
the [multi-state memory policy](../../policies/3D/separate_alignment_and_reconstruction_for_multistate_peak_mem_reduction.md)
is a hard boundary. Preserve 2D batch/class/half accumulation order initially.

### RNG is the dependency that must be resolved before batch fusion

Resolved by the keyed, stateless generator of `prob_align_batch_pipeline.md`
(section 4; phase 1 of the run, 2026-11-04): every draw of the probabilistic
routes is keyed on the run's seed, the pass, the iteration, the particle and
the draw. The reasoning below is kept as the record.

The mathematical row independence does not establish byte identity after
reordering. Candidate evaluation draws in-plane/shift samples; assignment later
draws in particle order on the master's intrinsic RNG stream. In shared memory
the main thread also evaluates candidates. Interleaving assignment with later
candidate evaluation changes that stream even at one thread. Longer-lived
workers also change `SIMPLE_SEED`'s process-local reseed counter.

First expose an explicit RNG contract for sample selection, candidate evaluation
and assignment, including restart and invalid rows that consume no assignment
draw. A later fusion can use independent logical streams keyed by iteration,
particle and phase, with any further candidate identity required by the draw
contract. This is a deliberate stochastic-trajectory change requiring truth
validation, not a byte-identical optimization. Merely fixing a seed or preserving
master assignment order does not settle it.

For an exact first phase, retain existing phase barriers, reseed semantics and
draw order. A resident worker can save reference setup across the barrier, but
its old batches have already gone: it cannot also eliminate every second polar
preparation without retaining/spilling those representations or changing RNG
ownership. Do not promise bounded memory, full batch reuse and the old RNG
schedule simultaneously without demonstrating how all three are met.

### Implementation order and review gates

1. **Measure the current boundary costs and expose typed entry points.** Split
   existing worker/matcher bodies into preparation, batch operations and finish
   routines. Use the caller's builder/project in shared memory; remove its
   command/toolbox and table/assignment file roundtrips. Initially retain table
   storage, barriers, numerical kernels and RNG behavior. Compare fixed-input
   rows, assignments, shift metadata, extras and noise output against the current
   path, not an older balanced-assignment implementation.
2. **Retain alignment state within an iteration.** Share compatible references,
   noise, plans and preparation capacity between evaluation and finalization.
   A distributed worker fold may retain these across assignment using the
   existing backend, with identical partition/sample coverage and bounded
   memory. It must not be presented as eliminating the second batch preparation
   automatically. Specify failure, cancellation and generation checks before
   introducing a waiting worker protocol.
3. **Separate RNG ownership, then fuse complete batches.** Extract one row
   selector shared by dense/sparse/neighbourhood adapters. Retain every eligible
   candidate, finish polish/noise while the prepared batch is live, and release
   its rows. If the master remains the selector, use a bounded in-flight window
   with backpressure and ordered decisions, not a gather of all rows. Cover
   missing/duplicate/stale records and empty partitions. Establish equivalence of
   the selection distribution and multi-seed scientific acceptance separately
   from the intentional trajectory change.
4. **Reuse the same pipeline in whole-command SPMD.** The existing
   [coarray SPMD proposal](coarray_spmd_execution.md) is the durable route to
   multi-iteration scientific residency and in-runtime publication/reduction.
   It needs the pipeline above; collectives alone do not fix object teardown or
   duplicated preparation. Keep transport/collectives in the execution adapter,
   not PFTC or scoring routines. Avoid a second independent persistent-worker
   framework. A per-iteration worker fold is an optional earlier adapter, not a
   prerequisite to implementing SPMD.

Acceptance combines unit assertions on complete row selection and metadata with
existing `unit_pftc_align2D3D`, fast gates, simulated workflows and multi-state
truth cases. Include sparse/dense/neighbourhood, topk, fractional sampling,
continuous polish on/off, sigma2 on/off, empty parts and multiple threads. The
exact first phases require observed fixed-seed identity; fused scheduling needs
the explicit RNG contract and repeated truth checks. Add counters for reference
loads, canonical and polar preparations, table bytes, optimizer/map work and
phase peak RSS. Target one preparation per batch in the fused path, unchanged
candidate/map work, and no shared-memory table transport. Unit/nightly runs can
establish these properties; representative PfCRT/msp1/difficult-sample wall time,
aggregate memory and NUMA measurements remain necessary to establish the gain.

Reference-neighbourhood tiling, SP weights, alternate contractions and GPU work
remain separate numerical/locality experiments. Search caps and refined-reference
budgets are scientific approximations and are not required to remove this
architectural overhead.

## Status, 2026-11-05

Done by the autonomous run `prob_pipeline` on the Dell, in the working tree
(`prob_align_batch_pipeline.md`, Progress: phases 0 to 4):

- Every random draw of the probabilistic routes (the sample, the in-plane draw
  per scored reference, the polar shift map's disk draw, the sparse fills'
  candidate order, the assignment draw, the diagnostics) comes from a keyed,
  stateless generator (a multiply-xor-shift hash of the run's seed, the pass,
  the iteration, the particle and the draw; `simple_rnd`), so a draw depends on
  none of the thread, the batch, the part or the process. The seed travels as
  `SIMPLE_SEED`; no new parameter.
- Probabilistic alignment is one batch pipeline inside `refine3D_exec` and
  `refine2D_exec`, in shared memory and in every distributed part alike: a
  batch is prepared once, its probability-table rows filled, each particle's
  pose or class drawn from its own row with the extras of the weighted
  insertion, then polished, its residual computed and committed. The master
  draws the sample, schedules its parts once, merges the parts' readout
  records in part order and assembles as before. The global table, the part
  files, the master's assignment step, the table programs and their commanders
  are retired.
- Measured (phase 4 of the note): the wall time per iteration and per phase of
  the trpv1 2D and 3D cases against the base build, the preparations per
  particle, and the peak resident memory of a worker, of the master and of the
  node; the identity of the drawn entries across thread counts, batch sizes and
  part counts holds (the fill harness), and a complete run across part counts
  agrees to the accumulators' tolerance frozen in phase 0.

## Status, 2026-11-03

Done on the Oracle Linux test machine, in the working tree (the maintainer's
request of 2026-11-03: sections 1.2 e, c and f, the worker fold, dropped, and
the single-precision weights):

- **Section 1.2 e, the 2D names and counts as 3D**: `refine2D` counts `minits`
  and `maxits` from `startit` (the `solve2D` controller plans in absolute
  iterations and emits run lengths; the stream's pool runs `maxits=1`); the 2D
  convergence guard is the new key `assign_all`, and `fillin` keeps its 3D
  meaning only, which 2D never had; 2D `refine=prob_snhc` is
  `refine=prob_neigh` with `prob_neigh_mode=snhc`, the only 2D mode.
  Byte-identical: distributed `solve2D` (two parts, one thread, fixed seed;
  full and fractional updates), the 1jxy and 6vxx workflows (path-derived
  bytes aside), `mini_stream`.
- **Section 1.2 c and f**: the distributed master draws the sample of the
  sampled stages without a table (`master_draws_sample2D`, as 3D's
  `master_draws_sample3D`), so the sample does not depend on `nparts` (one
  `refine2D` iteration at `update_frac=0.4` draws the same sample at two and
  three parts); the 2D partitions, the matcher's and the table's, balance the
  active particles. Changes results. Checked on the beta-galactosidase
  selection (2,756 of 5,513 particles active, 40 classes, `nsample=1000`, four
  seeds each): median class resolution 9.50 A before and 9.62 A after, the ten
  best classes 6.73 and 6.76 A, classes better than 10 A 21.8 and 20.8, within
  the seeds' spread.
- **The worker fold, 3D and 2D** (section 1.3, stretch): implemented, then
  dropped (ruling 2026-11-04). Where every part ran at once (the local queue, a
  computing unit per part), one job per part filled its table part, waited for
  the master's assignment and searched. It was byte-identical (one thread per
  part: distributed `solve2D` of 1jxy, distributed `solve3D` of 6vxx) and
  measured neutral on this machine: `refine2D` of trpv1 (35,645 particles, 100
  classes, box 256, four parts of six threads, two iterations) took 138.4 s
  folded and unfolded (three runs each, 124-146 s), `solve2D` of the
  beta-galactosidase selection 32.7 against 33.8 s. It saved a process start, a
  project read and (2D) a particle read per part and iteration, and it needed a
  second scheduling path, because a queue that cannot run every part at once
  would deadlock on a part waiting for the assignment. The gain counted above
  needs the table's prepared particles and references handed to the search;
  without that hand-over the fold is not worth its branching.
- **Single-precision weights in the score kernels** (section 1.3): the euclid
  rotation scan (`euclid_sums`) takes the single-precision weight
  (`sigma2_wk`); the `_8` routines of the shift search stay double. At the
  kernels (`polar_matching_bench suite=identity`, 32,000 rows per kernel): the
  distances within 1.6e-7 relative, one unit in the last place, and no row's
  best rotation moved. Per pair, pinned to the performance cores, 13.0 to
  11.8 us at one thread (9.7 %), 19.8 to 18.0 us at eight (9 %). The FSC
  check: three `refine3D` iterations of trpv1 (35,645 particles, c2, from its
  `solve3D` map, `refine=prob`, crop 160, four parts of six threads) give the
  same FSC=0.143 and FSC=0.5 resolutions at every iteration (4.04, 4.10,
  4.10 A; 6.10, 6.62, 6.77 A); the curves differ by at most 0.02 in a shell,
  as runs at six threads differ.

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
now (Part 2); rejected on 2026-11-04 as too memory-intensive.

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

The architecture review above is carried out by
[`prob_align_batch_pipeline.md`](../completed/prob_align_batch_pipeline.md) (planned
2026-11-04, ruled done): per-particle draws on keyed random streams, one batch
pipeline for shared memory and the distributed parts, no global table and no
part files. It absorbs the rows that stood here for the worker fold, the
streaming assignment, the read-ahead thread and the 2D table pre-step (the
first two superseded, the read-ahead already in the read session). What is left
for the next release:

| item | what | why later |
| --- | --- | --- |
| batched, GPU-ready kernels | pure kernels over explicit arrays (references, particle times CTF, shifts, rotations) scoring many references per call, first batched per particle on the CPU, inside the batch of the pipeline | the GPU release; FFT plans and precision change results |
| persistent workers, or the coarray SPMD route (`coarray_spmd_execution.md`) | worker processes or images that live for the run and exchange only model, sample and decisions; perhaps 5-15% per iteration on queued clusters | large; after the pipeline, which is the unit an image runs |

Rejected (ruling 2026-11-04): the NU filter as an object with the states
assembled in parallel. Holding several states' volumes and filter state at
once is too memory-intensive.

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
