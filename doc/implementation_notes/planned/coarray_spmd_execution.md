# SPMD execution of refine2D and refine3D on coarray images

Planned 2026-11-01, on master `eefc7e1a1`. **Whole-command SPMD is not implemented.** It records what the coarray backend
does today, the maintainer's goal for it, the proposed design, the rules that go with it, and the measurement run
to make on the Dell before any production code is written.

Reviewed 2026-11-01 against `20d07e53f`, including proposal commit `d0c7adecf`. The review corrections and
implementation gates below are recommendations, not additional maintainer rulings. No SPMD code or performance
results are supplied by this review.

The maintainer's statement of the goal (2026-11-01): coarray support is in SIMPLE, built by `compile_coarrays.sh`,
tested and working, with no serious performance benchmark yet. The script-based distributed execution we have
used for years should give way to a coarray-based distributed execution model that broadcasts more and relies less
on the file system, with proper coarray reductions of the partial class averages and the partial volumes. Cluster
support is not required in the first stage; the target is multi-core workstations, where OpenMP alone scales
poorly. Two rulings from that conversation are already folded in: the SPMD unit is what we launch as a
multi-iteration `refine2D` or `refine3D` command (not a `solve2D`/`solve3D` stage, which is too specific to those
applications), and no coarray is ever an allocatable component of a derived type.

## Execution architecture follow-up (2026-11-06)

The batch pipeline in `ea4c1c20a` supersedes this note's probability-table collection and assignment broadcasts.
Every image can fill, draw, refine and retain insertion poses locally. The shared iteration input is the frozen
model and sample; the outputs crossing images are changed particle/sigma2 records, raw additive accumulators and
small diagnostic counts. Batches do not synchronize with one another. The whole-command coarray backend remains
planned; the preparatory pipeline is implemented and validated by the maintainer.

Maintainer scope: clean up the architecture first; no coarray implementation or scoring-kernel redesign in this
refactoring. The current implementation work preserves the scientific settings and the alignment/reconstruction memory
boundary. Its contract, order and validation are kept here rather than in a second plan:

1. Give the 2D matcher explicit ownership of each allocated search strategy and release its target after use.
2. Restrict sparse 2D candidate indexing and storage to the current batch, preserving each particle's candidate
   order, global particle ID, keyed draw, seed shift and diagnostic sampling rank. Add uneven-batch regressions.
3. Use existing owner routines for shared-memory assembly and sigma2 commits; keep the authoritative project in
   the caller and separate project ownership from construction of replaceable execution resources. Preserve
   state-at-a-time restoration, crop/band changes and the repeated follow-up pass.
4. Select batch capacity before allocation from a common phase budget, accounting for fixed reference storage,
   per-particle work and the retained 3D observation buffer. Thread count bounds concurrency, not memory capacity.
5. Preserve the updated local-selection SPMD protocol in this note. Whole-command launch, collective payloads,
   placement and failure handling require the pinned coarray runtime and multi-image validation before release.
6. Profile the resulting scorer before introducing reference tiles/batched transforms. Keying the simplex start
   requires a deliberate random-stream version transition and independent truth checks.

Validation: build the current baseline before source edits; preserve its executables and frozen-input outputs.
Build incrementally, run the fast gate and the affected probability-table, reconstruction and numerical suites;
compare frozen-input fills across threads, batches and parts. Run simulated shared/distributed workflows and
Cartesian follow-up where ownership/lifecycle changes apply. Record actual results below, including unrun platform
and performance checks. No commit or push.

Progress: the baseline was rebuilt and its production executables saved. Implemented in the working tree:
allocatable 2D strategy ownership, batch-local sparse indexing, procedure-local matcher state, a resident
shared-memory refine3D project/general toolbox, owner assembly/sigma2 calls, and removal of redundant 2D master
particle writes. Intrinsic random-stream phase resets are explicit in the shared-memory strategies, preserving
the previous constructor-triggered reset sequence without reseeding inside owner routines. The existing keyed
generator/version is untouched. The follow-up has one lazily constructed parameters object; its iteration and
sampling settings change per pass. The alignment/reconstruction teardown boundary remains.
The strategy also restores the requested shift window before each alignment pass: parameter reconstruction
previously discarded convergence's widened window, and the resident lifecycle preserves that behavior explicitly.

The common memory planner caps preferred probabilistic batch capacity before allocation and deducts estimated
alignment storage from the 3D observation-buffer budget. These are dominant-array estimates, not a total-RSS
guarantee; fixed storage over budget still warns, with one particle slot and no retained buffer. The existing
fallback reconstructs from stack reads.

Verified on the maintainer's Mac with the Debug GNU Fortran 16 build (tests on, coarrays off):

- Incremental build and installation pass. The final fast gate passes 15/15 entries, including 1554 polar-alignment
  and 660 reconstruction assertions. Its 109.8 s under concurrent workflow load exceeds the warning-only 60 s budget.
- The 2D table regression passes 90 assertions, including uneven batches with nontrivial global IDs and seeded or
  refined shifts. The memory sub-suite passes 13 closed-form capacity and geometry checks.
- Baseline/current frozen-input fills agree in every non-timing field and probabilistic readout across 48
  configurations (threads 1/4, uneven batches, parts 1/3, states 1/2 and four fill routes).
  A final-source repeat at threads 1/11 also passes all 48 configurations with zero differences.
- A matched two-iteration serial 3D replay produces byte-identical maps, FSC, partial reconstructions and
  reprojection models (12 products). It exposed and verified the explicit shift-window reset above.
- A matched two-iteration serial 2D replay produces byte-identical class averages, FRC and class-state data
  (8 products) and all 79 particle records agree. The 3D replay's 1000 pose/sampling records agree; `cfar` now
  retains the current assembly diagnostic instead of the stale file-handoff value.
- The same 2D replay over two distributed parts agrees in all 79 particle records and all 8 class/FRC/state
  products, validating the reduced master publication path.
- A four-thread, two-iteration 2D replay completes under bounds checking; all 79 particle records are finite,
  have valid class IDs and show the two updates. Two further runs pass the same checks with actual 11-thread
  search teams, verified in the per-pass benchmark records.
- The full simulated workflow's `solve2D` portion matches the baseline in all 79 particle records and 84 binary
  artifacts. The canonical sigma2 file differs only in its path-dependent layout digest (header word 14);
  spectra and all other header/transaction bytes agree.
- Registry, FLEX DAG, source descriptors, navigation regeneration and document-link checks pass.
- `multistate_workflows suite=states nthr=4` exits 0 after all map, pose, population and negative-control checks.
  The planted populations (1500/900/600) are recovered without a state-label error; median rotation/shift errors
  are 1.25242 degrees / 0.0471645 pixels, within the existing 2-degree / 0.1-pixel gates.
- Complete `cont_refine3D_1jxy nthr=4` exits 0, including shared/distributed Cartesian refinement, two-iteration
  sampled follow-up and trailing populations, automatic follow-up, in-place continuation and NU bundle handoff.
- The full 1JXY shared-memory `solve3D` agrees with baseline in all 327 scientific readouts, all 79 final
  particle records and 41 binary artifacts. Canonical sigma2 differs only in its path-dependent layout digest;
  all spectra and transaction metadata agree. Truth-map checks pass with the same reported values.

- Complete `simulated_workflow suite=1jxy nthr=1` exits 0 on baseline and final source, including shared-memory
  and two-part 3D execution. All 28 workflow metrics agree and pass. All 204 binary/image/box products are
  accounted for: 201 are byte-identical; the three canonical sigma2 files differ only in layout digest, with
  every other byte equal. All 79 two-part particle records and all 250 two-part scientific readouts also agree.
  No product is missing, extra or numerically different.

The pre-rebase implementation completed this local validation; its source hashes matched the frozen
executables used for those comparisons.
The initially installed `build/bin` executables were stale; their early runs were stopped and excluded from
evidence. Comparisons use frozen copies of freshly linked `build/production` executables. Commands, hashes and
comparison reports are under `/tmp/simple-execution-architecture/`; final logs use `/tmp/simple-execution-bounds-*`.
No Linux/BOX scaling or production performance measurements have been made. No commit or push.

Rebase validation (2026-11-06): replaying `bd4d6513d` onto `1bb56ecd9` conflicted in the test suite imports,
test UI help and test policy row. The resolution keeps upstream's parent-module tester interfaces and polar
shift-map kernel suite together with the matcher-memory suite. The automatically merged refine3D strategy
retains both the execution-ownership refactor and upstream's staging-test interface. Navigation was regenerated.
The combined GNU Fortran 16 Release build/install and all 15 fast-gate entries pass (35.85 s), including 139 kernel
and 13 memory checks; registry, FLEX DAG, descriptors and diff checks pass. The earlier full workflow/identity
runs above were not repeated for this registration-only conflict resolution. Rebase continuation is left to
the maintainer under the repository's no-commit rule.

## Terms

An **image** is one process of a coarray program; `num_images()` of them run the same executable and are numbered
from 1. **SPMD** (single program, multiple data) means every image runs the same code path over its own share of
the data and the images meet at agreed points. A **collective** is an operation all images call together:
`co_sum` adds an array across images (into every image, or into one `result_image`), `co_broadcast` copies an
array from one image to all, `co_max` takes the maximum; `sync all` is a barrier. **OpenCoarrays** is the runtime
library gfortran's `-fcoarray=lib` builds link against; it implements images as processes of the **MPI** (Message
Passing Interface) library beneath it, and its launcher `cafrun` is a thin wrapper around `mpirun`. **NUMA**
(non-uniform memory access) describes a machine whose cores reach some of its memory faster than the rest; a
two-socket workstation has at least two NUMA nodes. The **controller** is the ordinary `simple_exec` process that
runs a command; a **part** is one worker's share of the particles (`params%part` of `params%nparts`, particles
`params%fromp` to `params%top`); the **master image** is image 1, the scientific output owner in the initial design
below. The **transport backend** is the coarray backend as it exists today.

## Why

Three reasons, in the maintainer's order. First, OpenMP alone scales poorly on large workstations: one process
with 64 threads can spread its working set over NUMA nodes and leaves other cores idle during serial sections.
SIMPLE's `NTHR_SHMEM_MAX = 20` limits the master thread helper; it is not a general cap on the shared-memory
matcher. Several processes with a smaller thread team each, placed within NUMA nodes, are worth measuring, and SIMPLE's
distributed mode already gives that shape, at the price of the file traffic below. Second, the distributed mode
moves everything through files, which on the production server means one spinning 12 TB disk. Third, the
distributed mode starts fresh worker processes for every step of every iteration: the project and prepared
references are re-read and local toolboxes rebuilt. The master already prepares the reprojection model once per
iteration; workers load it. Persistent images also make local buffer reuse possible, but existing teardown must
change explicitly before any such reuse is obtained.

The measured shares from the ribosome tally of 2026-10-13 (29,101 particles, 5 parts of 16 threads, 78 iterations,
6,572 s, `doc/implementation_notes/rejected/particle_caching_3D.md`): the probabilistic pre-step 46%, volume
assembly 11%, the master's reprojection model preparation 11%, scheduling residual about 6.7 s per iteration or
8%, the matcher's slowest part 11%. This is historical evidence, not a profile of the current code. Its pre-step
timer did not separate scoring, setup, table I/O and assignment. Collectives target the transport portions of both
assembly and the pre-step; model preparation and scoring compute remain and need separate current measurements.

## What the backend does today

The coarray backend is a transport for one partitioned step (`doc/policies/coarray_execution_policy.md`). The
controller runs the distributed strategy loop as before; for each step it stages the per-part command strings and
launches one child, `cafrun -np <nimages> simple_private_exec --coarray <from_part> <to_part> <strings>`, with
`nimages = min(computing units, parts)`. Each image runs the partition commanders of parts `this_image()`,
`this_image() + num_images()`, ... in-process, redirecting its output per part, then all images `sync all` and
exit. The data between steps still travels through files: alignment documents, partial
reconstructions, partial class sums. Per refine3D iteration that is one launch of the `refine3D`
workers (two with the follow-up pass), each of which reads the project, builds its local toolbox, loads the prepared
reprojection model and reads the sigma2 state again. Only `simple_private_exec` links the coarray copy of the library;
`this_image()`, `num_images()` and `sync all` appear only in `production/simple_private_exec_driver.f90` under
`#ifdef USE_COARRAYS`. The default launcher string is a bare `cafrun -np`, with `SIMPLE_COARRAY_SUBMIT_CMD` as the
override. The policy's own section 11 names in-runtime reductions as the intended follow-up of this foundation.

A worker-team reduction could already replace many part files with one aggregate handoff. The ordinary
controller would still need a file or another explicit interprocess boundary, because it is outside the image
team. The full in-memory iteration handoff requires the whole multi-iteration command to live inside that team.

## What the file system carries per iteration

| Payload | Shape | Size at typical settings | Written by | Read by |
| --- | --- | --- | --- | --- |
| partial reconstruction, gridding (`write_state_half_partial`) | per part, state and half: compressed native-grid Fourier numerator and real density, `(box_crop/2+1, box_crop, box_crop)` for an even box | box 256: 67.63 MB complex + 33.82 MB real = 101.45 MB per half; 202.90 MB per part/state; 8 parts, 2 states: 3.25 GB | each worker, once per iteration | the master, summed part by part (`reduce_partials` in `simple_volassemble_strategy.f90`) |
| partial reconstruction, PCG (`accumulate_worker_state_half`) | per part, state and half: raw full-range complex B and real D, with backend-specific layout/provenance | measure the PCG payload separately; do not substitute the gridding dimensions | each worker | the PCG master |
| partial class sums (refine2D) | per part: even and odd Fourier numerators, CTF-squared denominators, counts and offsets for `ncls` classes | 500 classes at box 128: 99.84 MB arrays plus about 0.01 MB metadata | each worker | `cavgassemble` on the master |
| alignment documents | per part: the updated orientation rows | a few MB | each worker | `merge_algndocs` on the master |
| sigma2 part files | per part: the updated particles' noise spectra | a few MB | each worker | `commit_sigma2` on the master |
| project file | whole | 100 MB to 1 GB for 10^6 particles | the master | every worker, every step |
| volumes, FSCs, class averages | per state or per iteration | 64 MB at 256 pixels | the master | every worker |

Sizes use decimal MB/GB and omit small file headers. The gridding writer calls `compress_exp` before writing:
`new_accumulator` and `alloc_rho` in `simple_reconstructor.f90` define the native shape, not `2*box_crop` cubed.
Each payload is written and read, so disk traffic is roughly twice its stored size. Partial reconstructions remain
a strong first reduction target; their dominance over tables and references depends on box, states and sample size.

## The design

### The SPMD unit and how it is launched

A multi-iteration `refine3D` or `refine2D` command run with `qsys=coarray` executes as one multi-image program for
its whole life: all iterations, the probabilistic pre-step, the search, the reconstruction or class-sum assembly,
the sigma2 commit and the convergence test. `nparts` images are launched, one part per image; the transport
backend's loop over several parts per image does not apply, since every image is a worker for the whole run.

The launch structure stays as it is at the outer level. The controller (`simple_exec`, normal build) runs the
command's commander; under `qsys=coarray` the commander submits one job for the whole command instead of one per
step, reusing qsys resource/partition staging with a distinct SPMD dispatch, and reads the project back when the child
returns. The `solve2D` and `solve3D` controllers therefore launch one `cafrun` per `refine2D`/`refine3D` command
they issue and run everything else (the start reconstruction, postprocessing, reseeding, state assignment) in their
own process as today. Moving a whole controller inside the images is possible later and nothing below prevents it;
the per-command unit already gets most of the gain. The stream pipeline and its 2D pool keep the script backend and
are not touched by this note.

The existing `--coarray` entry cannot be reused unchanged: private `prg=refine2D`/`refine3D` dispatches to a
single-iteration worker, and `submit_coarray_jobs` caps image count below `nparts` when resources are fewer. A
distinct SPMD entry must select the multi-iteration strategy without recursive submission and verify exactly one
image per part before entering any iteration. Reject an unavailable image/thread budget; do not silently run
collective commands inside the transport backend's strided partition loop.

Use an SPMD strategy lifecycle that reuses existing matcher/domain operations. Neither current in-memory strategy
is directly safe to run on every image: initialization, project writes/rereads, assembly and finalization contain
owner actions. Worker `part` also selects `l_distr_worker` behavior that writes and destroys partial results.
Keep explicit local partition bounds alongside the owner's global project view; master convergence, population
counts and restoration must use the global range, not image 1's worker range. No strategy, commander or domain
object calls `this_image()`; the driver and communication adapter own image identity.

The outer controller settles the execution directory/project once; the child receives `mkdir=no` and settled
paths. Image 1 creates/verifies shared halves, initial poses/states, references, sample history and sigma2 state;
other images adopt that state and construct local resources. Initially allow one startup project read per image,
then exchange changes in memory. Only image 1 owns progress/finished markers, shared cleanup and scientific output;
per-image diagnostic logs are uniquely named. Initialization and shutdown need the same ownership audit as iteration.

### One refine3D iteration on images

Reading `inmem_execute_iteration` and `distr_execute_iteration` in `simple_refine3D_strategy.f90` side by side,
the parts meet at these points, in this order within an iteration:

1. **Model, sample and local pipeline.** The master image prepares the iteration model once, then broadcasts its header
   and prepared arrays. Include frequency band, dimensions, state/direction layout, crop sampling, prepared Cartesian
   volumes when needed and optional independent-half diagnostic sections. Receivers allocate matching payloads and
   build process-local plans/memoization. The master draws the fractional-update sample over the whole project (today
   `master_draws_sample3D` on the master, reproduced by the workers with `sample4update_reprod`) and broadcasts the
   sampled indices, counters and sample size; workers adopt that exact sample without drawing or incrementing it
   again. Publish the frozen refined-reference budget, shift-prior width, active-state layout and temperature
   controls with that model. Every image runs the existing batch pipeline locally: fill, keyed row selection,
   extras, continuous refinement and residual. Rows and insertion poses remain local; no table collection or assignment
   broadcast is needed. Images can run different batch counts and capacities without synchronizing per batch.
2. **Search and merge.** Each image runs the local registration kernel with its own threads. At registration
   completion, gather changed numeric orientation fields by particle ID and verify coverage/overlap before applying
   them to the master's field. This must be an explicit hook/phase boundary: `refine3D_exec` currently performs
   reconstruction before it returns. Keep the pass buffer alive across the boundary, without a second stack read.
   The initial recommendation is to broadcast the committed field changes to all images; otherwise explicitly
   publish every global population/state-layout quantity that reference preparation and the next table derive.
   Locally current rows plus stale remote rows are not a valid global field.
3. **Reconstruction.** `calc_3Drec` (gridding) and `execute_rec3D_pcg_worker` (PCG) already build one half of one
   state at a time, write it and kill it. The reduction goes exactly where the write is: for each state and half, in
   the same order on every image, the half's numerator and density are summed into one image, which accumulates
   them into its even or odd reconstructor; after a state's two halves the owning image restores that state
   (`restore_state_from_parts` without `reduce_partials`), writes the volume and FSC and the trailing chain. With
   `result_image = 1` everything restores on the master image, as the master does today. Sum raw numerators and
   densities before trailing blend, density floors, folding, FSC/regularization or restoration; keep states and
   halves separate. Global merged rows and realized sample/population counts must already be available. PCG reduces
   its raw full-range B/D and contributor/mass metadata before its existing master finalization.
   Rotating the reduction root alone does not parallelize restoration: that owner blocks the next collective while
   it restores. Parallel state restoration is later work requiring an explicit staging/overlap schedule, a memory
   budget and a decision about which image may publish each state's outputs.
4. **Sigma2 and convergence.** Gather updated noise-spectrum rows by particle ID, validate transaction coverage,
   and commit once on the master; these rows are disjoint records rather than additive statistics. The convergence
   test runs on the master image's global field; its outcome, the matching low-pass limit and other per-iteration
   scalars the workers need are broadcast as small typed arrays.
5. **Writing.** Only the master image publishes scientific output: the changed project segments at the checkpoint,
   the volumes and half maps, the FSC and its plot, the heat maps, the iteration states and strategy benchmarks. The per-part
   benchmark records (`REFINE3D_BENCH_PARTppp.txt`, `REC3D_BENCH_PARTppp.txt`, `PROB_TAB_BENCH_PARTppp.txt`) keep
   their names, one file per image.

For `pose_cont=yes`, repeat the registration/merge phase: discrete registration without reconstruction, merge,
Cartesian follow-up of the same sample, final merge, then one reconstruction from refined poses. Retain one sigma2
transaction, with refined residuals replacing the discrete rows. This uses repeated communication phases, not a
promise that the follow-up requires no communication.

Persistent processes first save launches and intermediate file traffic. Buffer persistence is separate work:
the shared-memory strategy now retains its general toolbox/project, but the matcher kills its search toolbox,
PFTC, Cartesian references and sigma2 object and releases its pass buffer each pass. Keeping references alive
through reconstruction would undo a deliberate peak-memory reduction. Define reusable resource lifetimes and
invalidation by geometry, band, reference generation, particle preprocessing and sigma2 generation before claiming
caching gains. The 2D cropped-particle disk cache is retired. Its batch pipeline already shares prepared observations between
scoring and search; the 3D pass buffer remains subject to an explicit per-process budget.

### One refine2D iteration on images

The same shape with fewer parts: before sampling, the master checks carried-state usability and broadcasts the
full-update/fractional decision. Missing or incompatible carried sums require a full update. The master draws the
sample and frozen model; each image fills and selects its own batch rows or candidates. Each image searches its
part and returns changed rows before assembly. Sum even/odd Fourier numerators, CTF-squared denominators and integer
half populations to the master. Validate common geometry, centering offsets and `l_frac`; do not sum offsets or
flags. When `l_frac` and compatible carried sums permit blending, the master shifts previous sums once and applies
the existing represented/sample population blend; otherwise it uses the full current sums. It commits the carried
set, then restores with the previous FRC cutoff, estimates FRCs and makes the class document.
Broadcast merged/even/odd references plus class/half populations and centering metadata. Current reference
preparation reads global populations from the particle field, so canonical metadata or synchronized rows are
mandatory. Sigma2 commit and convergence use the global field. One restoring image followed by a reference broadcast is the initial design.

### Collectives, not coarray variables

Every exchange above is a `co_sum`, `co_broadcast` or `co_max` over an ordinary array, plus `sync all`. No coarray
variable is needed. The rules:

1. No coarray is ever an allocatable component of a derived type (ruled 2026-11-01). Coarray variables of any kind
   are avoided altogether in the first implementation; if a gather of variable-length pieces ever needs one, it is a
   single intrinsic-type buffer in the collectives module, never in a domain object.
2. Every image calls the same collectives in the same order. An image whose part has no particle in a state or
   half still calls the reduction for that state and half with a zero array. Existing empty-sample returns and
   empty-state skips occur before the write seam; they must reach a common publish phase instead of returning.
3. Reductions are made into one image (`result_image`) unless every image needs the result; broadcasts carry
   intrinsic-type arrays only, never characters or derived types.
4. Use bounded contiguous chunks and 64-bit shape/byte arithmetic before allocation or packing. Classic MPI counts
   limit individual messages; runtime/datatype byte limits can be lower than an element-count estimate. Choose and
   test a conservative chunk size on the pinned runtime; do not infer safety from the nominal box alone.
5. Failure on any image is `error stop` (what `THROW_HARD` already does through `simple_exception`), which
   requests termination of the team. Test nonroot failure and use a launcher deadline: mismatched collectives,
   killed processes and runtime failures are not proved hang-free by the language-level error path alone.

One module owns the collectives, with the intrinsic calls under `#ifdef USE_COARRAYS` and a single-image body
otherwise (`num_images() = 1`, reductions and broadcasts leave their arguments as they are), so every source goes
on compiling in both library copies as today. Its interface is in SIMPLE's terms, not MPI's:
`reduce_half_accumulator(numerator, density, result_image)`, `collect_particle_rows(...)`,
`broadcast_iteration_model(...)`, `broadcast_scalars(...)`, `barrier()`, `is_master_image()`, `nimages()`. An
alternative to the `#ifdef` is compiling the normal library with `-fcoarray=single`, which makes the module plain
Fortran in both builds; whether that costs compile time or binary size is measured before it is chosen
(`doc/policies/compile_time_policy.md`).

The seam in the matcher is the existing write: `write_state_half_partial` (gridding) and the raw write of
`accumulate_worker_state_half` (PCG) become a publish step that either writes the part file (script and transport
backends, unchanged) or reduces (SPMD), chosen by the strategy and passed down as it passes `l_write_partial_recs`
today; `calc_3Drec` does not learn about images. The same holds for the class sums, whose seam is
`cavger_write_contribution` (one `cavg_sums` contribution per part, read back by `cavger_assemble_sums_from_parts`),
and for the tables. The adapter is explicitly owned and passed through the workflow; numerical modules gain no
hidden communication globals. Domain objects expose bounded numerical payloads through their owning APIs, not
whole-object serialization or diagnostic getters that copy full volumes.

### Policy amendments this needs

`doc/policies/coarray_execution_policy.md` section 9 confines image identity to the private-exec driver. The
collectives module becomes the second allowed home, and the policy gains the rules above, the SPMD unit, and the
statement that the transport backend remains for the partitioned steps that are not multi-iteration commands
(`make_cavgs`, `reconstruct3D`, the preprocessing steps), which keep their part files.
`doc/policies/layering_and_project_ownership_policy.md` section 5 gets the stronger form of its ownership rule for
this mode: the master image is the one owner, and no image reads the project file after the start of the command.
The modes-equivalent rule (shared-memory and distributed runs behave the same; differences are fixed, not
documented) extends to three modes.

### Determinism, memory, threads, placement

With `result_image`, `co_sum` is a root-directed reduction; without it, the result is returned to every image.
OpenCoarrays' MPI implementation distinguishes `MPI_Reduce` and `MPI_Allreduce`. Addition order can differ from
the file backend and across image counts. Compare fixed-input accumulators with justified numerical tolerances;
when exact part order is required, collect parts and add in that order. Nonroot arguments become undefined after
a root-directed sum: reset or reload them before reuse, even if a particular runtime leaves their bytes intact.

Stochastic behavior needs a separate contract. `SIMPLE_SEED` uses a process-local seed-call counter, so persistent
images change streams compared with fresh workers even at one thread. The master alone draws the shared sample; probabilistic selections already use particle-keyed draws. The
remaining intrinsic simplex start needs its own stream contract before changing process lifetimes. Tolerances on sums do not establish equivalence of stochastic trajectories.

Budget memory per phase and in aggregate across images, including the root's simultaneous worker and owner roles.
Replicated global fields, padded tables, gather buffers, references, retained particles and runtime scratch are
additional costs; removing disk files does not remove their resident payloads. Reference memory depends on the
actual `pftsz`, frequency band, states, halves and memoization, not box alone. Prefer bounded table slabs, direct
contiguous accumulator access and one state/half at a time. Measure peaks at registration, reduction and restoration
before enabling persistence; do not assume the new mode has the same peak as a current worker.

The thread budget is images times threads at most the cores, with the master image using the same `nthr` as the
others; do not automatically carry over the script mode's master-side thread boost. First measure whether the
waiting images spin, sleep or yield and whether their assigned cores can usefully be borrowed. Runtime progress
policy and an eventual state-restoration schedule are measured tuning choices, not guaranteed speedups. Run
collectives on the image's main thread outside OpenMP regions; no concurrent collective calls from worker threads.

Placement decides whether a benchmark means anything. Binding defaults differ by MPI version and launcher.
Use a runtime-specific explicit mapping to disjoint core sets large enough for each image's OpenMP team, and record
both process and thread affinity. Identical `OMP_PLACES=cores` on unbound images does not by itself partition cores
between them. `--bind-to none` is only a control case unless each image also receives a distinct CPU set. Keep MPI
specific launcher syntax in the backend configuration and preserve `SIMPLE_COARRAY_SUBMIT_CMD` overrides.

### Restart and failure

A command run as one program loses the per-part restart of the script backend by design. What remains is the
iteration-level checkpoint the master image writes (the project segments, the volumes or references, the trailing
chains and canonical sigma2 state). Define how a complete iteration is recognized: in-place writes of several
artifacts are not an atomic checkpoint. Publish a completion record only after all required artifacts are closed
and consistent, and reject an incomplete generation on restart. A checkpoint must retain sampled/update counters,
geometry, state layout, carried-population provenance and the chosen RNG restart contract. The in-memory refine3D
strategy rejects `continue=yes` today; the SPMD strategy should take the distributed path's continuation semantics,
since the users of `qsys=coarray` are the users of the distributed mode. A failed image ends the whole command with
the per-image logs (`simple_private_exec_coarray_part_*.out`) as the diagnostics, as the policy already says.

### Clusters later

The communication contract can run across nodes without node-local shared memory assumptions. A later cluster
mode also needs consistent particle/source paths, executable/runtime deployment, cache placement and input
availability checks; changing launcher/binding alone is insufficient. Keep cross-image exchange collective-based,
without relying on memory-mapped files or `/dev/shm` shared between images.

### What this does not fix

Collectives do not reduce the arithmetic of scoring or model construction. NUMA placement can improve locality
but does not add hardware memory bandwidth. The pre-step needs a fresh compute/I/O/assignment breakdown; related
kernel work is in `doc/implementation_notes/completed/polar_matching_layout_bench.md` and
`doc/refactoring_notes/completed/prob_search_speedup.md`. The master still prepares the canonical model once;
parallel construction is optional later work. The initial wins are file round trips and launches. Buffer retention
and locality are separately measured opportunities, not gains implied by the payload table alone.

### Testing

Only `simple_private_exec` links the coarray library, so its thin platform entry must call ordinary tester bodies
using `simple_test_utils` assertions. Extend the existing capability-gated `coarrays` CTest entry; no new standalone
test program or executable. Adding subcases does not change the CTest budget; a new CTest entry would require the
budget, registry and policy count changes. The current smoke test establishes launch/synchronization only.

Validate fixed-payload transport first: root-directed and all-image sums, complex arrays, broadcast headers and
payloads, chunk tails, zero contributions, repeated buffer reset, exact-once coverage and missing/duplicate row
rejection. Then compare a fixed-reference/full-update iteration's accumulators and restored products with a
serial or file-backend oracle before adding probabilistic sampling. Cover fractional carry-over/centering,
globally empty states, sparse/dense tables, sigma2 row coverage, crop/band changes, continuation and refined poses.
Exercise 1, 2 and uneven 3-image particle partitions with multiple threads, plus deadline-bounded nonroot failure.

Use existing truth fixtures for scientific acceptance in all three modes and repeat over fixed seeds where streams
differ. `test=identity` is bytewise; truth-floor checks alone are not mode equality. Add explicit numerical
comparisons with independently justified tolerances and preserve exact metadata/coverage checks. A one-image SPMD
run is especially useful for isolating lifecycle errors before reduction order changes enter.

## Order of work

1. The measurement run below. It decides whether the gain on a workstation justifies the rest and sizes it.
2. Settle the launch/reentry, global-versus-local state and lifecycle contracts below. Add the communication
   adapter, single-image fallback, tests within the existing platform entry and runtime-specific placement.
3. A one-image SPMD refine3D lifecycle, then a fixed-reference/full-update gridding iteration on multiple images:
   reference broadcast, local registration, changed-row collection, raw half reductions and existing root restoration.
   Keep memory bounded and preserve the single-read particle path; amend policies when implementation lands.
4. Verify the local batch pipeline on the broadcast frozen model and sample, canonical sigma2 commit and convergence
   across iterations. Gate unsupported combinations explicitly until validated; do not silently change their route.
5. PCG raw B/D and metadata; Cartesian follow-up using the repeated registration/merge contract.
6. The SPMD refine2D lifecycle, local dense/sparse pipelines, raw class-sum reduction, owner carry-over and reference broadcast.
7. Measured resource retention, optional distributed model construction, and genuinely overlapping state restoration.
   Each needs its own invalidation/scheduling and peak-memory evidence within this same living note.

## Review findings and implementation gates

The per-command SPMD unit, one owner and raw-accumulator reductions are sound. The following issues must be
resolved before treating the proposal as implementation-ready; the body above incorporates the review corrections.

| Priority | Finding and required gate | Current source anchors |
| --- | --- | --- |
| P1 | A new launch/reentry and owner/local lifecycle are required. Current private commands are workers, while in-memory strategies perform global setup, writes and teardown. Prove a one-image SPMD command first. | `production/simple_private_exec_driver.f90:300,315`; `src/utils/qsys/simple_qsys_ctrl.f90:715`; `simple_refine3D_strategy.f90:1119,1223`; `simple_refine2D_strategy.f90:156,185` |
| P1 | Registration must expose committed rows before any root restoration. Empty samples/states must still reach the same publish sequence. The follow-up repeats registration/merge on the same sample. | `simple_strategy3D_matcher.f90:129,295,319`; `simple_matcher_3Drec.f90:53,85`; `simple_refine3D_strategy.f90:1794` |
| P1 | Global metadata and complete table schemas are part of transport. Preserve exact-once coverage and authoritative active-state/class layout; local-only row updates are insufficient. | `simple_eul_prob_tab.f90:160,1423`; `simple_eul_prob_tab2D.f90:1178,1334`; `simple_matcher_pftc_prep.f90:48`; `simple_refine2D_strategy.f90:566` |
| P1 | Reduce raw statistics, then apply population carry-over/trailing and nonlinear restoration once. Validate common offsets/geometry and PCG provenance separately. | `simple_classaverager_restore.f90:720,821`; `simple_volassemble_strategy.f90:119,140`; `simple_rec3D_pcg_strategy.f90:252,1735` |
| P2 | Persistence, numerical equivalence and speedups each need their own evidence. Current teardown is intentional, RNG streams change with process lifetime, and timing residuals do not isolate I/O. | `simple_strategy3D_matcher.f90:303,311,326`; `simple_rnd.f90:57`; `simple_refine3D_strategy.f90:1732`; `simple_ptcl_cache.f90:1` |

Basenames above refer to the existing refinement, matcher, probability, class and volume modules under `src/`;
these anchors were checked at the review revision, not inferred from the older ribosome note.

**Recommended communication contract.** Agree a small descriptor before payload allocation: command generation,
iteration/pass, model generation, sample identity, geometry/band, canonical state/class layout and payload counts.
Give every payload an owner, numerical kind, shape, lifetime and validation rule. Common metadata is verified or
broadcast; additive statistics are summed; disjoint particle records are collected by ID. No whole-`parameters`
copy or builder serialization. Keep the owner and local worker roles explicit on image 1 without creating a second
builder for the same project. Stable partition identity remains independent of transport image identity.

**Choices still to settle.** The review recommends synchronized changed particle rows for the initial version;
a smaller metadata-only replication scheme needs a complete consumer audit. Decide raw-payload chunk/gather budgets,
RNG stream/restart policy and initial supported option matrix (`prob_neigh`, sparse 2D, PCG, continuous follow-up, fractional
updates and continuation). Keep expensive buffer retention and parallel state restoration deferred until measured.
An unsupported combination must fail before collective execution, with the supported route stated explicitly.

**Runtime evidence.** Standard semantics are sufficient for this design, but runtime performance is a separate
gate. [GNU's broadcast documentation](https://gcc.gnu.org/onlinedocs/gfortran/CO_005fBROADCAST.html) requires equal
array shapes; [its sum documentation](https://gcc.gnu.org/onlinedocs/gfortran/CO_005fSUM.html) specifies root-only
results and undefined nonroot arguments. The inspected
[OpenCoarrays MPI source](https://raw.githubusercontent.com/sourceryinstitute/OpenCoarrays/main/src/runtime-libraries/mpi/mpi_caf.c)
distinguishes reduce/all-reduce and contains an elementwise array-broadcast path. Benchmark the exact production
buffers on the installed, pinned build before assuming a bulk broadcast. If cost scales with element count as
individual messages, resolve the runtime/adapter choice before relying on reference broadcasts for performance.
[Open MPI 5's launcher guidance](https://docs.open-mpi.org/en/v5.0.x/man-openmpi/man1/mpirun.1.html) calls out
thread-aware binding; defaults are not a portable placement contract.

Review validation: source/caller inspection, independent payload arithmetic and documentation checks only. No
build, numerical collective test, Dell run, Linux/BOX result or speedup is claimed. The measurement run below
remains proposed work.

## The measurement run on the Dell

The Dell (24 cores, 62 GB, Oracle Linux 8, gcc-toolset-15) is the performance reference (x86 is the target; the
Mac only has to build and pass). The run is an `arun` run; it changes no production code except the one test
entry it needs, and it ends in a report under `doc/implementation_notes/planned/`. Everything below runs
unattended once the data and the launcher are in place.

**Preconditions, recorded first.** The machine topology (`lscpu`, `numactl --hardware`: whether the 24 cores are
one or two NUMA nodes, which decides how much of a NUMA effect can show at all). Whether OpenCoarrays is installed
(`cafrun`, `caf --show`, the MPI beneath it and its version) and whether its version builds with gcc-toolset-15; if
it is not installed, a user-space build of MPICH or Open MPI plus OpenCoarrays in the run's scratch is the first
phase. The default binding of the installed `mpirun` (`mpirun -np 4 --report-bindings true` or MPICH's
equivalent), and the launcher string the run will use for every coarray case. The data: the simulated vcRNAP set
(`/data/simulated_vcrnapol`, 30,000 particles, box 180, 2.055 Å, two states with truth maps) is the main set,
large enough that per-iteration fixed costs do not dominate and small enough for many configurations in a night;
bgal (5,513 particles) is the quick check only.

**What is measured, and from where.** Every case runs with the benchmark output that is always on (`L_BENCH_GLOB`,
files under `BENCH/`): `REFINE3D_STRATEGY_BENCH.txt` for the iteration's phases on the master, `REFINE3D_BENCH_PARTppp.txt`
and `REC3D_BENCH_PARTppp.txt` for each part's search and reconstruction, `PROB_TAB_BENCH_PARTppp.txt` and
`PROB_ALIGN_BENCH.txt` for the pre-step, `VOLASSEMBLE_BENCH.txt` for the assembly; the 2D cases their counterparts.
From these, per iteration: wall time; each sequential worker phase's critical path (the maximum across its parts);
nonoverlapping master phases; directly timed reads/writes and launches/waits; and the remaining **unattributed
overhead**. Do not label wall time minus worker/master timers as transport: matcher totals include partial writes,
the pre-step includes setup/I/O/assignment, and model time can be nested inside initialization. State timer nesting
and phase boundaries before summing. The bytes of part files per iteration are taken from
the run directory before the iteration's cleanup (or computed from `box_crop`, `nparts`, `nstates`), and the peak
resident memory per process from the existing `report_rec_master_phase` line and process-tree sampling. Report
simultaneous aggregate RSS as well; a launcher's `/usr/bin/time -v` value is not the sum of all images' peaks.

**Scaling of the existing modes.** The same `refine3D` on the single-state half of the vcRNAP set (the 17,000
particles of state 1, started from the truth map low-passed to 15 Å, `refine=prob`, `objfun=euclid`, a fixed
sample size, 8 iterations, no continuous follow-up), on 24 threads in total except the optional single-NUMA-node control:

| Case | Mode | Images or processes x threads |
| --- | --- | --- |
| shared memory, default | in-process | 1 x 24; record actual matcher and assembly thread counts |
| shared memory, pinned | in-process | 1 x 24 with `OMP_PROC_BIND=close OMP_PLACES=cores` |
| shared memory, one node | in-process | 1 x (cores of one NUMA node) under `numactl --cpunodebind --membind`, if the box has two nodes |
| script backend | `qsys=local` | 2 x 12, 4 x 6, 8 x 3, with matching disjoint worker affinity |
| transport backend | `qsys=coarray` | 2 x 12, 4 x 6, 8 x 3, with the explicit binding |

For the launch comparison, match script and transport affinity, thread counts, starting checkpoint, seed, sample,
crop/band, storage and completed iteration count. Record physical cores versus SMT and runtime versions. Separate
warm and cold file-cache conditions and vary run order; different stochastic trajectories can change the work.
Their end-to-end difference then estimates the backend effect, not an isolated launch timer. Shared versus
distributed comparisons also change orchestration and work partitioning, so do not attribute the entire gap to
NUMA. Each configuration is repeated; report all times and repeat noisy pairs rather than selecting the faster run.

**The same in 2D.** `refine2D` on the same 17,000 particles, 100 classes, 10 iterations, in the shared-memory
default, the pinned shared-memory, and the 4 x 6 and 8 x 3 script and transport configurations. 2D is where OpenMP
scaling is weakest (short per-particle work, the polar references streamed from memory), so the multi-process gain
should show most clearly there.

**The collectives themselves.** A thin platform benchmark entry in `simple_private_exec`, with its checks in a
normal tester, run under `cafrun` with 1, 2, 4 and 8 images: `co_sum` into image 1 of native gridding arrays at
`box_crop` 180 (`91 x 180 x 180`, 23.59 MB complex plus 11.79 MB real per half), then stress cases at box 256
(67.63 MB plus 33.82 MB) and box 512 (539 MB plus 269 MB). Compare root-directed and all-image sums. Broadcast
actual polar-model, Cartesian-volume and 2D reference payloads, including complex types and the ranks/chunks used
by production; a 64 MB real array alone is not a sufficient proxy. Record contiguous versus packed access,
warm-up, allocation/packing separately, and the maximum per-image elapsed time over five repetitions. Restore
inputs before each reduction; compare against a closed-form or serial oracle with absolute and relative tolerances.
Record the pinned OpenCoarrays build: implementation of a standard collective need not give bulk-bandwidth
performance for every array rank/type. Record the CPU time of images waiting at the final barrier while image 1
sleeps for a few seconds, which shows whether
blocked images spin. Against these numbers stand the file equivalents from the refine3D cases: the time
`reduce_partials` takes to read and sum the parts, and the time the workers take to write them.

**What the report has to answer.** Whether several processes beat one on this machine, and by how much of it is
placement (the pinned row) versus other process-shape effects; measured I/O/launch time and unattributed overhead,
in 3D and 2D; what the same payload transfer and reduction cost as collectives; and
therefore the expected gain of the SPMD design on a 24-core box with a local disk, with the payload table scaling it
to larger boxes, more states and the production server's disk. The report leads with those four numbers and keeps
the per-case tables after them. If directly measured transport is under 10% in every 3D case and the pinned shared-memory
row closes most of the gap to the script rows, the design's workstation case rests on the persistent memory and the
2D results alone, and the report says so.

**Budget.** Pilot one case before committing to a night. The table has eight base 3D configurations, plus an optional
NUMA-local case. At two minutes per iteration, eight iterations and two repeats, that is about 4.3-4.8 hours for
3D alone; allow 2D, warm-up, repeated noisy cases and runtime setup separately. These are estimates, not observed runs.
