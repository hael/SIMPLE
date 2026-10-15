# Particle I/O on per-micrograph stacks: review report

Report of the unattended run that carried out
[particle_io_few_particle_stacks.md](particle_io_few_particle_stacks.md),
written for the maintainer's review of the complete diff. It covers what
changed, the decisions taken without the maintainer, the deviations from the
plan, the evidence against the baseline, and what is left open. The plan's
progress section (section 11) holds the detail.

The run ran on 2026-10-13 and 2026-10-14 on the Oracle Linux test machine,
building with GCC 15.2 and CMake 3.26. Its data was beta-galactosidase: 5,513
particles in one stack, box 256 at 1.275 Å per pixel. Nothing is committed.
Each phase's diff is `/home/elmlundho/agent_runs/particle_io/review/phase_N.diff`,
and the cumulative diff after each phase is `full_after_N.diff`. The evidence
(scripts, logs and kept benchmark files) is under
`/home/elmlundho/agent_runs/particle_io/scratch/`.

## What changed

The layout on disk is unchanged and no particle is written to another file.
SIMPLE reads the same bytes in fewer, larger and earlier requests, labels its
half-sets per micrograph, and lets chained add-ons stop re-reading earlier
runs' particles.

- **Phase 0, the benchmark split** (`phase_0.diff`). Particle reading is timed
  apart from preparation in the 2D and 3D matchers, the reconstructions, the
  PCG worker (preconditioned conjugate-gradient reconstruction) and
  class-average restoration. The probability-table parts got a benchmark file
  of their own. Every reading pass reports the stacks it opened, its requests
  and the bytes requested and read; every benchmark file reports process and
  node memory. The baseline cases were then run on beta-galactosidase:
  - S, the `solve3D_cavgs` start;
  - A, `solve3D` from S;
  - B and C, `refine3D` distributed and in shared memory;
  - E, `solve2D`;
  - F, `emulate_solve3D_stream` with a base run and three add-ons;
  - P, A with the PCG backend.

  Four high-level tests were run with them.
- **Phase 1, the read session and coalesced spans** (`phase_1.diff`). A new
  pass-scoped read session (`src/fileio/simple_stack_read_session.f90`) opens
  each stack once per pass and keeps the stack at a batch boundary open. It
  reads the images a stack contributes in spans: a span bridges at most 2 MiB
  of unneeded bytes, and the scratch of a whole session is at most 64 MiB.
  Every particle reader moved to it, and the old batch reader was deleted.
  Read-only opens no longer check for existence twice.
- **Phase 2, the POSIX path and the prefetch window** (`phase_2.diff`). Stacks
  are opened concurrently on plain read-only descriptors, and the MRC header is
  decoded once per open. A producer thread reads up to 128 MiB, from up to 64
  stacks, ahead of the caller. The helpers live in `simple_posix.c`, with
  no new CMake target. Short reads, interrupted reads, unsupported advice,
  descriptor exhaustion and the fallback pixel modes are tested.
- **Phase 3, the particle buffer** (`phase_3.diff`). A 3D search fills a
  pass-sized slab of decoded images (`src/fileio/simple_image_slab.f90`). The
  same pass's reconstruction copies its images from the slab instead of reading
  them again. The slab stays within a budget: half the task memory a queue
  system grants, otherwise the node's memory over twice its processes. When it
  does not fit, the batches are read from the stacks as before.
- **Phase 4, project-owned half-sets** (`phase_4.diff`). One transaction in the
  master (`src/main/project/simple_eo_partition.f90`) labels the even and odd
  halves:
  - per micrograph group, largest group first, from its own random stream;
  - the same labels in `ptcl2D` and `ptcl3D`;
  - the provenance on the stack rows.

  Workers refuse to make a partition. Fields without micrograph rows, class
  averages and time series keep the per-index route. Every former call of the
  per-index routine is now a consumer of the transaction.
- **Phase 5, PCG out of the add-on chains** (`phase_5.diff`).
  - `solve3D_addon` refuses a frozen run that used PCG, naming its backend,
    before anything is written, and forces gridding.
  - `emulate_solve3D_stream` forces gridding.
  - The PCG inputs left both programs, and the PCG side of the frozen sets was
    deleted with its test.
  - Two rulings of 2026-10-14 were applied in this phase as well (below).
- **Phase 6, the add-on passes** (`phase_6.diff`).
  - The first add-on of a chain reads the base run's particles once at each of
    its two boxes, the stage box and the native box. It publishes a carried
    bundle: the union's raw sums at both boxes, the rows inside them, and a
    bundle manifest written last.
  - Later add-ons adopt the bundle and read no frozen particle.
  - The sigma2 noise model is a copy-on-write transaction. The frozen records
    are copied, the new particles start from the frozen run's half curves, and
    the means are republished over the union.
  - The final reconstruction is incremental: the add-on's own particles plus
    the carried sums.
  - A retired particle stays in the sums. A row inside the sums that is
    selected again stays retired.
- **Phase 7, documents and close** (this report).
  - A new short policy, `doc/policies/particle_read_path_policy.md`.
  - The reader policy in `project_stack_indexing_policy.md`.
  - The half-set text in `refine3D_policy.md`, `solve2D_policy.md` and
    `reconstruct3D_pcg_policy.md`.
  - The add-on and chunk policies, updated in Phases 5 and 6.
  - The final cases, and the plan moved to `completed/`.

## Decisions taken under delegation

Each decision is recorded, with its reason, in the plan's section 11 under the
phase named.

- *Case F uses 1,100 particles per add-on, not 1,200* (Phase 0). 1,200 gives
  two add-ons on this set, and the plan asks for three.
- *Case F starts from the control project* (Phase 0).
  `emulate_solve3D_stream` refuses the inputs that would start it from S.
- *The read session is a builder-free module in `src/fileio`* (Phase 1), so the
  unit tests drive it on synthetic stacks. Its 64 MiB scratch cap is for the
  whole session, shared by its read streams.
- *The `solve2D` acceptance is judged against every run of the Phase 0 code*
  (Phase 1). This gives a score of at least 0.377 and a class resolution of at
  most 12.11 Å. Read narrowly, against Phase 0's first two runs only, the rule
  is failed by the unchanged code itself.
- *A case with two Phase 0 values is accepted within the band between them,
  widened at each end* (Phase 3). For B this is 6.73 to 7.31 Å.
- *The prefetch buffer belongs to the C side, and the producer thread is
  unbound from the OpenMP master's core* (Phase 2).
- *Running out of descriptors stops the run* (Phase 2), naming the stack.
- *The buffer's budget under a queue is half the granted task memory*
  (Phase 3).
- *The half-set provenance lives on the stack rows* (Phase 4). A particle
  record cannot say "unlabelled": an unset half reads as even.
- *Labels from an earlier release are kept* (Phase 4) and recorded as
  per-index. Beta-galactosidase's control project already carries such a split.
- *The `make_cavgs` master makes the partition before its workers start*
  (Phase 4). This site was not in the plan's inventory; it became necessary
  once workers refuse.
- *A frozen run with no recorded backend counts as gridding* (Phase 5), the
  `solve3D` default. The backend and the PCG controls are no longer replayed
  from the manifest, and on the add-on's command line they are refused by
  name.
- *The carried route starts from the frozen run's final maps* (Phase 6). These
  are the restoration of the carried native sums.
- *A row reselected after its retirement stays retired while it is inside the
  sums* (Phase 6); searching it would count it twice.
- *The high-level `solve3D_addon` gate became a chain of two add-ons* (Phase 6),
  with every existing floor applied, unchanged, to the final union.
- *The read-count script counts images as bytes over the image size* (Phase 6).
  Since Phase 1 a request is a span, not one image.

## Deviations from the plan, and the rulings during the run

- **Rulings of 2026-10-14**, all recorded in the plan:
  - *The static review was deleted on purpose.* Its link left the plan's
    header.
  - *The add-on track serves the offline chains only.* The stream's 3D stage no
    longer runs add-ons. The stream ingestion policy left the file table, and
    the plan and both add-on policies describe the offline chains.
  - *No new environment variable.* The particle buffer's `SIMPLE_PTCL_BUFFER_MB`
    from Phase 3 was removed in Phase 5.
  - *The Phase 5 fence stop.* The ruling above had dropped the environment
    file's row from the file table, so the driver stopped after Phase 5. The
    row came back with phase 5, and the file stays at its base content.
- **Not done in this run, as the plan says:**
  - stacks kept open from the search into the reconstruction (they wait for
    the access-trace replay);
  - the production measurements on the ribosome project;
  - the generated indexes under `doc/code_overview/`, which still list a few
    deleted procedures and need regenerating;
  - skill files.
- **Tests the plan lists, met differently or not exercised:**
  - *Rollback and restart of a chain* are covered by construction: a run writes
    only its own directory and its bundle manifest last. No test run exercised
    them.
  - *`ptcl2D` and `ptcl3D` labels after prune, re-extraction and stream
    publication* are covered through the merge test: the stack rows travel with
    their particles.
- **The plan's `oris` setter files and `sp_project` files** were assigned to
  Phase 4 but needed no change: the transaction lives in its own module.

## Evidence against the baseline

No timing gain is claimed: the test set is one page-cached stack, which these
changes barely touch. Every read-path step was kept by the timing rule, since
none was slower than the build without it beyond the run-to-run spread. The
numbers are in section 11, per phase.

**Bit-exactness.**
- The stack-reader tests compare every image byte for byte with a single read,
  in pixel modes 0, 1, 2, 6 and 12.
- End-to-end single-threaded reconstructions are byte-identical before and
  after Phases 1 and 2.
- One iteration with and without the particle buffer is byte-identical.

**Final cases** (Phase 7, Release build; one run each against Phase 0):

| Case | Phase 0 | Phase 7 | Acceptance |
| --- | --- | --- | --- |
| A (`solve3D` from the class-average start, distributed) | 14:05, 14:09; FSC 0.143 4.35, 4.35 Å | 14:12; 4.41 Å | holds (limit 4.48 Å) |
| B (`refine3D`, five iterations, distributed) | 2:38, 2:38; 6.94, 7.10 Å | 2:48; 6.80 Å | holds (band 6.73 to 7.31 Å) |
| C (B in shared memory) | 1:46; 7.77 Å | 1:46; 7.77 Å | holds |
| E (`solve2D`) | 3:06, 3:08; score 0.389, class resolution 10.97, 10.94 Å | 3:09; 0.389, 11.05 Å | holds |
| F (`emulate_solve3D_stream`: a base run and three add-ons) | 39:22, 43:00; union FSC 0.143 4.03 / 3.89 / 3.80, 3.84 / 3.75, 3.80 Å | 30:21; 4.03 / 3.93 / 3.89 / 3.84 Å | holds (within 5 %) |

Peak memory is unchanged, except C in shared memory: 3.92 GB against 3.63 GB,
which is the particle buffer. Case F's add-ons read 4,000 / 0 / 0 frozen
particle images, against 12,000 / 18,600 / 25,200 in Phase 0. That is the first
add-on reading the base run's 2,000 particles once at each of its two boxes.
The add-ons took 261 / 310 / 252 s per 1,000 added particles, against
360, 322 / 578, 572 / 306, 520 s.

**High-level tests** (Phase 7, each judged by its own checks; `simulated_workflow_6vxx`
is recorded only):

| Test | Phase 0 | Phase 7 |
| --- | --- | --- |
| `solve3D_addon` | passed: union map against the truth, correlation 0.942, FSC 0.143 4.93 Å | passed (now a chain of two add-ons): 0.948, 4.74 Å |
| `simulated_workflow_1jxy` | passed: 0.898, 4.54 Å | passed: 0.879, 4.46 Å |
| `simulated_workflow_6vxx` | passed: 0.888, 4.71 Å | passed: 0.878, 4.71 Å |
| `cont_refine3D_1jxy` | passed: median rotation errors 0.91 / 0.23 / 0.91° | passed: 0.91 / 0.23 / 0.91° |

**Gates.**
- `./compile_debug.sh` builds with no compiler warning, and the fast gate
  passes 15 of 15.
- The Release build passes, and every unit suite passes.
- `scripts/check_test_registry.py` passes, and no file mode changed.

## What the maintainer measures next

On the ribosome project (751,729 particle rows in 19,520 stacks), section 4 of
the plan:
- the access-trace replay, which decides whether stacks should stay open from
  the search into the reconstruction;
- the stream-count sweep on spinning disks;
- the span gap budget and the scratch cap, named constants still at the plan's
  starting values (2 MiB and 64 MiB);
- the benchmark lines of section 5, cold and warm cache, read against the
  Phase 0 split of the same project.

## Open for the maintainer

- **Proposed skill text** (not applied):
  - `simple-main-strategies` and `simple-modern-fortran`: every particle
    reader goes through the pass's read session (`ptcl_read_session`); none
    opens a stack itself.
  - `simple-main-project`: the even/odd partition is the project's
    (`ensure_field_halves`), made by the master; workers never make one.
  - `simple-refine3d`: `solve3D_addon` carries its sums; chains are gridding
    only.
- **Regenerate** the generated indexes under `doc/code_overview/`.
- **A published add-on project names its sigma2 state by a bare name** that
  resolves only in its run directory. This is unchanged behaviour: consumers
  take the state from the run manifest's absolute path. It was noticed in
  Phase 6.
- **Disk:** each adopted add-on keeps its own complete carried bundle, plus
  the copy it adopted from its parent. Each set holds the even and odd raw sums
  and densities at the stage box and the native box, a few hundred MB per state
  at box 256. A chain's earlier run directories can be removed once the next
  add-on has adopted the bundle.
