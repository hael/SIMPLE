# SPMD execution of refine2D and refine3D on coarray images

Planned 2026-11-01, on master `eefc7e1a1`. **Nothing here is implemented.** It records what the coarray backend
does today, the maintainer's goal for it, the proposed design, the rules that go with it, and the measurement run
to make on the Dell before any production code is written.

The maintainer's statement of the goal (2026-11-01): coarray support is in SIMPLE, built by `compile_coarrays.sh`,
tested and working, with no serious performance benchmark yet. The script-based distributed execution we have
used for years should give way to a coarray-based distributed execution model that broadcasts more and relies less
on the file system, with proper coarray reductions of the partial class averages and the partial volumes. Cluster
support is not required in the first stage; the target is multi-core workstations, where OpenMP alone scales
poorly. Two rulings from that conversation are already folded in: the SPMD unit is what we launch as a
multi-iteration `refine2D` or `refine3D` command (not a `solve2D`/`solve3D` stage, which is too specific to those
applications), and no coarray is ever an allocatable component of a derived type.

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
`params%fromp` to `params%top`); the **master image** is image 1, the only image that writes files in the design
below. The **transport backend** is the coarray backend as it exists today.

## Why

Three reasons, in the maintainer's order. First, OpenMP alone scales poorly on large workstations: one process
with 64 threads spreads one heap over every NUMA node, runs every serial section on one core while the others
idle, and SIMPLE caps the shared-memory master at `NTHR_SHMEM_MAX = 20` threads for that reason. Several
processes with a dozen threads each, one per NUMA node or a few per node, are the standard answer, and SIMPLE's
distributed mode already gives that shape, at the price of the file traffic below. Second, the distributed mode
moves everything through files, which on the production server means one spinning 12 TB disk. Third, the
distributed mode starts fresh processes for every step of every iteration, so nothing survives from one iteration
to the next: the project is re-read, the toolbox and the reprojection model rebuilt, the particle images re-read,
and the in-memory particle cache that the shared-memory mode can hold is unavailable.

The measured shares from the ribosome tally of 2026-10-13 (29,101 particles, 5 parts of 16 threads, 78 iterations,
6,572 s, `doc/implementation_notes/rejected/particle_caching_3D.md`): the probabilistic pre-step 46%, volume
assembly 11%, the reprojection model 11% (built again on every worker), scheduling overhead about 6.7 s per
iteration or 8%, the matcher's slowest part 11%. Coarray reductions address the assembly, the scheduling and the
redundant rebuilds; they do not touch the 46%, which is compute.

## What the backend does today

The coarray backend is a transport for one partitioned step (`doc/policies/coarray_execution_policy.md`). The
controller runs the distributed strategy loop as before; for each step it stages the per-part command strings and
launches one child, `cafrun -np <nimages> simple_private_exec --coarray <from_part> <to_part> <strings>`, with
`nimages = min(computing units, parts)`. Each image runs the partition commanders of parts `this_image()`,
`this_image() + num_images()`, ... in-process, redirecting its output per part, then all images `sync all` and
exit. The data between steps still travels through files: alignment documents, probability tables, partial
reconstructions, partial class sums. Per refine3D iteration that is two launches (`prob_tab`, then the `refine3D`
workers; three with the polish pass), each of which reads the project, builds the toolbox, builds the reprojection
model and reads the sigma2 state again. Only `simple_private_exec` links the coarray copy of the library;
`this_image()`, `num_images()` and `sync all` appear only in `production/simple_private_exec_driver.f90` under
`#ifdef USE_COARRAYS`. The default launcher string is a bare `cafrun -np`, with `SIMPLE_COARRAY_SUBMIT_CMD` as the
override. The policy's own section 11 names in-runtime reductions as the intended follow-up of this foundation.

Adding `co_sum` to this backend as it stands would gain little: the images that hold the partial sums exit at the
end of the step, and the process that needs the sum, the controller, is not an image. The reduction becomes natural
only when the unit of a coarray job is the whole multi-iteration command.

## What the file system carries per iteration

| Payload | Shape | Size at typical settings | Written by | Read by |
| --- | --- | --- | --- | --- |
| partial reconstruction, gridding (`write_state_half_partial`) | per part, state and half: the Fourier numerator (complex) and the density (real) of the padded volume, padded box `2 * box_crop` | `box_crop` 256: 257 x 512 x 512 complex = 540 MB plus 270 MB real; 1.6 GB per part and state; 8 parts, 2 states: 26 GB | each worker, once per iteration | the master, serially, summed part by part (`reduce_partials` in `simple_volassemble_strategy.f90`) |
| partial reconstruction, PCG (`accumulate_worker_state_half`) | per part, state and half: the raw accumulators of the solver | the same order of magnitude | each worker | the PCG master |
| partial class sums (refine2D) | per part: even and odd class sums and CTF-squared sums for `ncls` classes at `box_crop` | 500 classes at 128 pixels: about 130 MB per part | each worker | `cavgassemble` on the master |
| probability table (`prob_tab`, `prob_tab2D`) | per part: the sampled particles' scores against every reference | 10,000 particles x 2,500 references: 100 MB in total | each worker | `read_tabs_to_glob` on the master |
| assignment and insertion tables | per iteration | a few MB | the master | every worker |
| alignment documents | per part: the updated orientation rows | a few MB | each worker | `merge_algndocs` on the master |
| sigma2 part files | per part: the updated particles' noise spectra | a few MB | each worker | `commit_sigma2` on the master |
| project file | whole | 100 MB to 1 GB for 10^6 particles | the master | every worker, every step |
| volumes, FSCs, class averages | per state or per iteration | 64 MB at 256 pixels | the master | every worker |

The partial reconstructions dominate by two orders of magnitude, and they are also summed serially on one process.
That is the first reduction to move into memory.

## The design

### The SPMD unit and how it is launched

A multi-iteration `refine3D` or `refine2D` command run with `qsys=coarray` executes as one multi-image program for
its whole life: all iterations, the probabilistic pre-step, the search, the reconstruction or class-sum assembly,
the sigma2 commit and the convergence test. `nparts` images are launched, one part per image; the transport
backend's loop over several parts per image does not apply, since every image is a worker for the whole run.

The launch structure stays as it is at the outer level. The controller (`simple_exec`, normal build) runs the
command's commander; under `qsys=coarray` the commander submits one job for the whole command instead of one per
step, through the existing `qsys_ctrl` staging and `submit_coarray_jobs`, and reads the project back when the child
returns. The `solve2D` and `solve3D` controllers therefore launch one `cafrun` per `refine2D`/`refine3D` command
they issue and run everything else (the start reconstruction, postprocessing, reseeding, state assignment) in their
own process as today. Moving a whole controller inside the images is possible later and nothing below prevents it;
the per-command unit already gets most of the gain. The stream pipeline and its 2D pool keep the script backend and
are not touched by this note.

Inside the images, what runs is the in-memory strategy (`refine3D_inmem_strategy`, `refine2D_inmem_strategy`) over
the image's own part, plus collectives at the points where parts meet. The image's part comes to it as it does
today: the driver stages the command string with `part`, `nparts`, `fromp`, `top`, the commander builds
`parameters` from it, and the strategy sees `params%fromp` and `params%top`. No strategy, commander or domain
object calls `this_image()`.

### One refine3D iteration on images

Reading `inmem_execute_iteration` and `distr_execute_iteration` in `simple_refine3D_strategy.f90` side by side,
the parts meet at these points, in this order within an iteration:

1. **Sample and pre-step.** The master image draws the fractional-update sample over the whole project (today
   `master_draws_sample3D` on the master, reproduced by the workers with `sample4update_reprod`) and broadcasts the
   sampled indices and the sample size; the reproduction contract disappears. Under `refine=prob` each image fills
   its part of the probability table in memory; the table slices are summed into the master image as a
   zero-padded `co_sum` (every image holds a full-size table with zeros outside its part; 100 MB per image is
   cheap); the master image runs `ref_assign` as today and broadcasts the per-particle assignment (and the
   insertion table when `topk > 1`) as integer and real arrays. No table or assignment file.
2. **Search.** Each image runs `refine3D_exec` over its part with its own threads. The orientation rows it updated
   are packed as a real matrix, zero-padded to the full sample, and summed into the master image, which applies them
   to its field; the workers keep their own rows. Alternatively the master image broadcasts the merged rows back so
   every image holds an identical field; the first implementation should take whichever the sigma2 commit and the
   convergence report need with the least code, since both run on the master image.
3. **Reconstruction.** `calc_3Drec` (gridding) and `execute_rec3D_pcg_worker` (PCG) already build one half of one
   state at a time, write it and kill it. The reduction goes exactly where the write is: for each state and half, in
   the same order on every image, the half's numerator and density are summed into one image, which accumulates
   them into its even or odd reconstructor; after a state's two halves the owning image restores that state
   (`restore_state_from_parts` without `reduce_partials`), writes the volume and FSC and the trailing chain. With
   `result_image = 1` everything restores on the master image, as the master does today. With
   `result_image = 1 + mod(state - 1, num_images())` the states restore on different images in parallel, which the
   script mode cannot do; the master image then registers the volumes in the project after a `sync all`. Start
   with the first, measure, move to the second when the assembly share says so.
4. **Sigma2 and convergence.** The updated particles' noise spectra are summed into the master image the same way
   and committed there; the convergence test runs on the master image's field and its outcome, the matching
   low-pass limit and any other per-iteration scalar the workers need are broadcast as one small array.
5. **Writing.** Only the master image writes: the project segments once per iteration, the volumes and half maps,
   the FSC and its plot, the heat maps, the iteration states, the benchmark records of the strategy. The per-part
   benchmark records (`REFINE3D_BENCH_PARTppp.txt`, `REC3D_BENCH_PARTppp.txt`, `PROB_TAB_BENCH_PARTppp.txt`) keep
   their names, one file per image.

What the images keep between iterations is the point of the exercise: the particle field, the sigma2 state, the
`pftcc` allocation, the toolbox and, once the 3D in-memory strategy stops rebuilding its toolbox every iteration
(it does today, `kill_strategy3D_tbox` then `init_params_and_build_strategy3D_tbox`, because the crop and the
limits change between iterations), the particle records of a cropped pass. The 2D particle cache (`cache=yes`),
which `particle_caching_3D.md` found unavailable in the distributed mode because workers run from scripts, becomes
available with no further work. The reprojection model is still built on every image, as every worker builds it
today; splitting the references across images and gathering them is a later step (section "Order of work").

### One refine2D iteration on images

The same shape with fewer parts: the master image draws the sample and, under `refine=prob`, the 2D table is
summed into it and the assignment broadcast; each image searches its part; the partial class sums (even, odd and
CTF-squared, per class) are summed into the master image, which runs the restoration (`cavger_restore_cavgs`, the
FRCs, the class document) as `cavgassemble` does today and broadcasts the restored references; the sigma2 commit
and the convergence test run on the master image. The class sums are small enough that summing into every image
and restoring redundantly would also work, but restoration is threaded and not free, so one restoring image plus a
broadcast of `ncls` images is the default.

### Collectives, not coarray variables

Every exchange above is a `co_sum`, `co_broadcast` or `co_max` over an ordinary array, plus `sync all`. No coarray
variable is needed. The rules:

1. No coarray is ever an allocatable component of a derived type (ruled 2026-11-01). Coarray variables of any kind
   are avoided altogether in the first implementation; if a gather of variable-length pieces ever needs one, it is a
   single intrinsic-type buffer in the collectives module, never in a domain object.
2. Every image calls the same collectives in the same order. An image whose part has no particle in a state or
   half still calls the reduction for that state and half with a zero array; `mark_empty_state` and the paired
   partial contract already express this on the file side and carry over unchanged.
3. Reductions are made into one image (`result_image`) unless every image needs the result; broadcasts carry
   intrinsic-type arrays only, never characters or derived types.
4. Arrays larger than 2^31 elements are reduced in chunks (MPI counts are 32-bit integers). A padded box of 1024
   (513 x 1024 x 1024 = 538 million) fits; a padded box of 2048 does not.
5. Failure on any image is `error stop` (what `THROW_HARD` already does through `simple_exception`), which
   OpenCoarrays turns into termination of every image, so no image waits forever at a barrier for one that died.

One module owns the collectives, with the intrinsic calls under `#ifdef USE_COARRAYS` and a single-image body
otherwise (`num_images() = 1`, reductions and broadcasts leave their arguments as they are), so every source goes
on compiling in both library copies as today. Its interface is in SIMPLE's terms, not MPI's:
`reduce_half_accumulator(numerator, density, result_image)`, `reduce_table(tab, result_image)`,
`broadcast_assignment(...)`, `broadcast_scalars(...)`, `barrier()`, `is_master_image()`, `nimages()`. An
alternative to the `#ifdef` is compiling the normal library with `-fcoarray=single`, which makes the module plain
Fortran in both builds; whether that costs compile time or binary size is measured before it is chosen
(`doc/policies/compile_time_policy.md`).

The seam in the matcher is the existing write: `write_state_half_partial` (gridding) and the raw write of
`accumulate_worker_state_half` (PCG) become a publish step that either writes the part file (script and transport
backends, unchanged) or reduces (SPMD), chosen by the strategy and passed down as it passes `l_write_partial_recs`
today; `calc_3Drec` does not learn about images. The same holds for the class sums, whose seam is
`cavger_write_contribution` (one `cavg_sums` contribution per part, read back by `cavger_assemble_sums_from_parts`),
and for the tables.

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

`co_sum` runs as an MPI all-reduce whose summation order is neither the part order of today's file sum nor fixed
across image counts, so volumes and class sums will not be bitwise equal to the script mode's, and may differ in the
last bits between runs with different image counts. The existing mode comparisons are tolerance-based for the same
reason across OpenMP schedules, so this is accepted; if a test ever needs today's exact order, the master image can
pull the parts in fixed order instead, at the cost of serializing the sum.

Memory per image is what a worker holds today (its part's data, the full reprojection model, one half
accumulator at a time); the master image adds what the master holds today (the global table, the even and odd
reconstructors of the state being restored). Nothing doubles. Replicating the reprojection model on every image is
the price of the distributed shape, as it is for the script mode: 2,500 references at `box_crop` 256 are about
650 MB per image.

The thread budget is images times threads at most the cores, with the master image using the same `nthr` as the
others; the script mode's master-side thread boost for the assembly does not carry over, because the waiting images
hold their cores. Two consequences to measure rather than assume: MPI processes blocked in a collective spin by
default (Open MPI's `mpi_yield_when_idle` and the equivalent in MPICH change that), which would steal cycles from
the restoring image; and distributing the restoration by state (point 3 above) is the structural answer to the idle
images, not a thread boost.

Placement decides whether a benchmark means anything. Open MPI by default binds a process to one core when two or
fewer are launched and to a socket otherwise; the first cramps an image's OpenMP threads onto one core, the second
lets several images' threads contend for one socket. MPICH's launcher does not bind by default. The default launcher
string must therefore carry an explicit mapping (`--bind-to none` with `OMP_PLACES`/`OMP_PROC_BIND` set per image,
or a `--map-by <unit>:PE=<nthr>` form), and the measurement run records the binding it used.

### Restart and failure

A command run as one program loses the per-part restart of the script backend by design. What remains is the
iteration-level checkpoint the master image writes (the project segments, the volumes or references, the trailing
chains), so a crashed run restarts from the last iteration as a shared-memory run does. The in-memory refine3D
strategy rejects `continue=yes` today; the SPMD strategy should take the distributed path's continuation semantics,
since the users of `qsys=coarray` are the users of the distributed mode. A failed image ends the whole command with
the per-image logs (`simple_private_exec_coarray_part_*.out`) as the diagnostics, as the policy already says.

### Clusters later

Nothing above is node-local: collectives over MPI run across nodes, and the master image is the only writer, so a
multi-node run differs only in the launcher and the binding. The one thing to keep out of the design is anything
that assumes images share a machine, such as memory-mapped files or `/dev/shm` between images; collectives only.

### What this does not fix

The probabilistic pre-step is compute (46% of the ribosome run), and more processes add no memory bandwidth to a
bandwidth-bound kernel; its speedups are the matter of `doc/implementation_notes/planned/polar_matching_layout_bench.md`
and `doc/refactoring_notes/completed/prob_search_speedup.md`. The reprojection model is still built on every image until it is split. The wins are
the file round trips, the per-step launches and rebuilds, the persistent per-image memory and NUMA locality, and the
sizes in the payload table bound them.

### Testing

Only `simple_private_exec` links the coarray library, so a multi-image unit test needs an entry beside
`--coarray-smoke`, for instance `--coarray-test <name>`, dispatching to testers that call the collectives module;
the capability-gated `coarrays` CTest entry grows into a small suite (the CTest budget and the policy count change
with every added entry; `scripts/check_test_registry.py` runs after the tables change). Mode equivalence is tested
as it is today for the two existing modes: the same fixture through shared memory, the script backend and the SPMD
backend, compared with the existing tolerances, on the small simulated sets of the truth suite.

## Order of work

1. The measurement run below. It decides whether the gain on a workstation justifies the rest and sizes it.
2. The collectives module, its single-image fallback, the `--coarray-test` entry and the launcher binding default.
3. The SPMD refine3D strategy on the gridding backend: the sample broadcast, the table sum and assignment broadcast,
   the search over the image's part, the half-accumulator reduction into the master image, the sigma2 and
   convergence points; one job per command from the commander; the policy amendments.
4. PCG partials on the same seam; the polish pass (`pose_cont=yes`), whose second search pass runs on the images
   with no new collectives.
5. The SPMD refine2D strategy: table, class-sum reduction, reference broadcast.
6. Restoration distributed by state; the reprojection model split across images and gathered; the toolbox rebuild
   narrowed to what an iteration changes, so the particle records persist.

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
From these, per iteration: the wall time; the slowest part's compute; the master's compute (assembly, assignment,
commit); and the **transport share**, the wall time less the slowest part's compute less the master's compute,
which is what the files, the launches and the polling cost. The bytes of part files per iteration are taken from
the run directory before the iteration's cleanup (or computed from `box_crop`, `nparts`, `nstates`), and the peak
resident memory per process from the existing `report_rec_master_phase` line and `/usr/bin/time -v`.

**Scaling of the existing modes.** The same `refine3D` on the single-state half of the vcRNAP set (the 17,000
particles of state 1, started from the truth map low-passed to 15 Å, `refine=prob`, `objfun=euclid`, a fixed
sample size, 8 iterations, no polish), in these configurations, all on 24 threads in total:

| Case | Mode | Images or processes x threads |
| --- | --- | --- |
| shared memory, default | in-process | 1 x 24 (the `NTHR_SHMEM_MAX` cap applies where it applies; record what ran) |
| shared memory, pinned | in-process | 1 x 24 with `OMP_PROC_BIND=close OMP_PLACES=cores` |
| shared memory, one node | in-process | 1 x (cores of one NUMA node) under `numactl --cpunodebind --membind`, if the box has two nodes |
| script backend | `qsys=local` | 2 x 12, 4 x 6, 8 x 3 |
| transport backend | `qsys=coarray` | 2 x 12, 4 x 6, 8 x 3, with the explicit binding |

The script and transport rows with the same shape differ only in how the parts are launched, so their difference is
the launch and polling overhead; the script rows against the shared-memory rows give the NUMA and serial-fraction
gain that any multi-process scheme gets before any reduction is moved into memory; the pinned and one-node
shared-memory rows tell how much of the poor OpenMP scaling is placement alone. Each case is run twice; the report
gives both and flags a pair that differs by more than 10%.

**The same in 2D.** `refine2D` on the same 17,000 particles, 100 classes, 10 iterations, in the shared-memory
default, the pinned shared-memory, and the 4 x 6 and 8 x 3 script and transport configurations. 2D is where OpenMP
scaling is weakest (short per-particle work, the polar references streamed from memory), so the multi-process gain
should show most clearly there.

**The collectives themselves.** A `--coarray-bench` entry in `simple_private_exec` (SIMPLE-native, no driver
scripts), run under `cafrun` with 2, 4 and 8 images: `co_sum` into image 1 of a complex array and a real array of
the padded volume at `box_crop` 180 (181 x 360 x 360, 188 MB plus 94 MB, the size the refine3D cases reduce) and
at a padded box of 512 (540 MB plus 270 MB); the same summed into every image; `co_broadcast` of a 64 MB real
array; each timed over five repetitions, with the result checked against the serial sum of the same inputs to a
relative tolerance, which also confirms that the installed OpenCoarrays reduces complex arrays. The CPU time of the
images that wait at the final barrier while image 1 sleeps for a few seconds is recorded, which shows whether
blocked images spin. Against these numbers stand the file equivalents from the refine3D cases: the time
`reduce_partials` takes to read and sum the parts, and the time the workers take to write them.

**What the report has to answer.** Whether several processes beat one on this machine, and by how much of it is
placement (the pinned row) versus the serial fraction; what share of a distributed iteration is transport, in 3D and
in 2D, which bounds what any reduction scheme can recover here; what the same reduction costs as a collective; and
therefore the expected gain of the SPMD design on a 24-core box with a local disk, with the payload table scaling it
to larger boxes, more states and the production server's disk. The report leads with those four numbers and keeps
the per-case tables after them. If the transport share is under 10% in every 3D case and the pinned shared-memory
row closes most of the gap to the script rows, the design's workstation case rests on the persistent memory and the
2D results alone, and the report says so.

**Budget.** At about two minutes per refine3D iteration for 17,000 particles on this machine, the nine 3D cases
run twice take about five hours; the 2D cases about two; the collectives bench minutes; a user-space MPI and
OpenCoarrays build, if needed, under an hour. One night.
