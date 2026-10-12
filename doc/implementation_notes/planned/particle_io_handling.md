# Particle I/O handling for SIMPLE's iterative workflows

Status: draft plan for review, 2026-10-12. No code yet. Replaces the
single-stack cache plan of the same day, which is set aside (section 8).

## 1. The problem

Most of our projects hold thousands of per-micrograph particle stacks. Every
workflow built on the refine2D or refine3D machinery updates a small, freshly
drawn sample of particles in each iteration, and reads each drawn particle
from its stack again. On slow storage these workflows wait on the disk, not
on the CPUs.

The ribosome base run of `solve3D_chunks` on 2026-10-11 ran before its runs
were cut into whole stacks. It held 29097 particles at box 336 (about 450 KB
per image) on a spinning disk. It sampled 4998 particles per iteration over
84 iterations: about 420000 reads, roughly 14 per particle. In a late
iteration the probabilistic pre-step took 140 of the 204 seconds, and the
CPUs were often idle.

The stream adds a second, larger cost. Every `solve3D_addon` run reads all
frozen particles, the particles aligned by earlier runs, several times.

- **Frozen accumulations** (`calc_frozen_rec`). The add-on runs one stage,
  the last one the base run ran. It accumulates the frozen particles once at
  that stage's box and once at the native box: two passes, or one when the
  two boxes are equal.
- **The final reconstruction** (`calc_final_rec`, which runs
  `bootstrap_rec3D`). It re-estimates the noise model (sigma2) of the whole
  union, frozen and new particles together. That reads every particle four
  times: for the image power spectra, for the bootstrap map, for the
  residual sigma2 pass, and for the shipped map.
- **The size of a pass.** One pass over 500000 particles at box 336 reads
  about 226 GB of float32: 20 to 25 minutes on one spinning disk.
- **Over a session.** An add-on starts each time the new particles reach 10%
  of the frozen ones. With up to six passes over the union per add-on, a
  session growing from 100000 to 750000 particles reads its final dataset
  about fifty times.

These stream figures are arithmetic, not measurements.

## 2. Principles

- **The stacks stay as written.** No particle is copied or rewritten into
  another file. A rewrite costs a full read and a full write of the data,
  which no stream can afford. The existing cropped 2D cache (`cache=yes` in
  `solve2D`) is left as it is, and so is the format particles are extracted
  in.
- **Better handling of the stacks we have.** Fewer opens and fewer requests,
  requests issued ahead of need and served in disk order, and a page cache
  that keeps a run's particles once they have been read.
- **Fewer reads.** Reads that repeat work are removed, not made faster.
- **One reader for every workflow.** `discrete_read_imgbatch`
  (`simple_matcher_ptcl_io.f90`) serves the 2D and 3D matchers, the
  probability tables, reconstruction, class averages and FLEX. A change there
  reaches all of them.
- **Measured.** Each change is timed in SIMPLE's benchmark files, on the same
  runs, before and after.

## 3. How particles are read today

- **Per batch.** A batch is 50 images per thread in the matchers and 500 in
  reconstruction.
  - `discrete_read_imgbatch` groups the batch by stack and works through the
    stacks, up to `io_read_nstreams` at a time, each in its own thread: 2 on
    a spinning disk, 16 on an SSD, 8 elsewhere.
  - For each stack it opens the file, reads the header, reads each image
    with its own request, and closes the file.
- **Across batches.** The next batch opens the same stacks again. A stack
  with particles in several batches is opened, has its header read and is
  closed once per batch.
- **Per iteration.** The sample is read twice: by the probability-table pass
  and by the search pass. The second read usually comes from the page cache.
- **Per reconstruction** (`calc_3Drec` in `simple_matcher_3Drec.f90`).
  - Particles are read by state and half (even and odd), with one
    reconstructor in memory at a time.
  - Each particle is read once, but every stack is visited once per state and
    half: six times with three states.
  - Under soft state weights, a particle is read once for each state it has
    weight in.
- **What a read costs.**
  - From the page cache: tens of microseconds.
  - From a cold spinning disk: a seek of about 10 ms per stack visit, plus
    the file system's metadata on open.
  - On SSDs and network file systems: the latency of each open and each
    request.

So the levers are:

- fewer visits and requests per stack;
- many requests in flight, served in disk order while the CPUs compute;
- a working set that stays in memory;
- fewer passes.

## 4. Reader changes

These change only the read path: `simple_matcher_ptcl_io.f90`, the stack
reader (`simple_discrete_stack_io.f90`) and the image file layer. The
workflows do not change. The cropped 2D cache keeps its own reader
(`ptcl_cache_read_batch`). Its build reads the stacks through
`discrete_read_imgbatch`, so the build gains too.

### Read ahead

At the start of a pass the reader tells the kernel which byte ranges the
pass will read (`posix_fadvise` with `WILLNEED`). The kernel queues them all,
sorts them by position on disk, keeps the drive's queue full, and loads them
while the CPUs work. A few thousand blocking reads in program order become
one sorted sweep in the background.

- **Matchers and probability tables.** The iteration's sample is known
  before the first batch, from the sampler or its sample file. The whole
  sample is hinted at once: 5000 images at box 336 are 2.3 GB.
- **Reconstruction.** Its scans can exceed memory, so it hints a window
  ahead: the next batch or two while the current one is gridded.
- **Offsets.** An image's offset needs its stack's first data byte and bytes
  per pixel, both from the header. The reader remembers them per stack for
  the life of the process (next subsection). A stack not yet seen is hinted
  whole; per-micrograph stacks are small.
- **Bounded.** Never more than a fraction of the memory the system reports
  available, so the hints do not evict what the run already holds.
- **Platform.** `posix_fadvise` is Linux-only. A small C helper in
  `simple_posix.c` makes it a no-op elsewhere, as `simple_fs_is_rotational`
  does. On network file systems the effect depends on the file system.

### Stacks kept open

A process-wide pool keeps stacks open across batches and passes, keyed by
stack name, and remembers each stack's header fields. It replaces the open,
header read and close per stack per batch.

- The pool is bounded by the open-file limit. When it is full, the least
  recently used stack is closed.
- It is emptied at the end of each workflow step and on `kill`.
- Each stack has one unit, as Fortran requires. The reader already gives
  each stack to one thread at a time.
- In distributed runs each worker has its own pool.

### One request per stack

The reader reads from the first to the last needed image of a stack in one
request, then converts the images it needs and drops the rest. It splits the
read only at gaps of more than a few MB. On a spinning disk, passing over a
short gap costs about as much as reading it, while a second request costs a
rotation or a seek. float16 stacks convert with the existing
`simple_float16_to_float32`.

### Stream count

The cap of 2 on a spinning disk was set in August 2026 to protect long
sequential reads from head thrashing; it was not measured. Scattered reads
of a few images per stack are the opposite case, and arrays of spinning
disks need more requests in flight than 2.

`SIMPLE_IO_NSTREAMS` already overrides the cap, so a sweep (1, 2, 4, 8, 16)
on the ribosome settles it before any code changes. With read-ahead in place
the kernel's queue does the reordering, and the threads matter mostly for
decoding.

## 5. Fewer passes in add-on runs

This is the stream's largest particle I/O. Your rulings of 2026-10-12 set
the direction:

- **Frozen contributions are never recomputed.** They keep the noise
  weighting they were accumulated with, and that error is accepted.
- **The final reconstruction of an add-on is incremental.**

These change the `solve3D_addon` contract
(`doc/policies/3D/solve3D_addon_policy.md`, sections 3, 5, 7 and 8).

### Why the frozen sums are recomputed today

The frozen particles' poses and labels never change during a session (stream
3D ingestion policy, item 10). What changes is their noise model.

- The frozen term is accumulated with the frozen run's committed sigma2
  state, while the cohort's sigmas are seeded and estimated during the
  stage.
- Every add-on's final reconstruction fits a new sigma2 state to the union,
  and that state becomes the next add-on's frozen sigma2.
- So the weight of every frozen particle changes from one add-on to the
  next, and the frozen sums are recomputed from the particles.

### Frozen sums carried forward

- **What is carried.** The frozen sums at the stage's box and the native
  box, as they were accumulated. They are never rescaled or recomputed when
  sigmas change.
- **What is added.** The next add-on's frozen sums are this add-on's plus
  this run's particles at their final poses, with their best sigmas (the
  last pass of the run, below).
- **The boxes** are the same in every add-on of a chain, since each runs the
  frozen run's last stage.
- **The first add-on** starts from the base run's frozen sets and its
  committed sigma2 state, computed once as today.

### The noise model

- **Global.** Each particle's sigma2 spectrum is estimated and kept in the
  sigma2 state (`sigma2_part` in `simple_euclid_sigma2.f90`). Alignment and
  reconstruction use one global sigma, the average of those spectra. There
  are no groups.
- **The best estimate is the one used.** Estimating the noise is part of
  refinement. As soon as there is a better estimate of a particle's sigma,
  that is the sigma used from then on.
- **New particles need a sigma before they can be aligned.** They start
  from the union's committed global sigma, which is estimated from
  residuals. Each new particle's spectrum starts as that sigma, and the
  cohort's are then updated from their residuals during the stage.
  - Today the add-on drops the inherited registration and seeds the cohort
    from image power (`calc_pspec`, through `ensure_sigma2_for_iteration`).
    Image power includes the signal, so it is the worse estimate. That
    seeding pass goes.
  - The base run, with no union before it, still seeds from image power.
- **Frozen particles keep their committed sigmas,** and their contributions
  keep the weighting they were accumulated with.

### Incremental final reconstruction

At the end of the run, today's bootstrap sequence runs over the run's own
particles instead of the union. It has no image power seed, since the
cohort started from the union's committed global sigma.

1. **Bootstrap map:** the carried frozen sums at the native box, plus the
   run's particles at their final poses with their sigmas from the stage.
2. **Residual pass:** the run's particles against the bootstrap map at
   native sampling. Their new sigmas are committed.
3. **Shipped map and new carried sums:** the run's particles accumulated
   again with the new sigmas, at the native box and at the stage's box.
   - The shipped map is the carried frozen sums at the native box plus these.
   - The next add-on's frozen sums are the carried sums plus these.

Then:

- **The output** carries the carried sums and a sigma2 state that merges the
  frozen particles' committed sigmas with the run's new ones. The next
  add-on consumes both.
- **The full fit** is left to the session's final run (stream 3D ingestion
  policy).
- **What goes:** every pass over the union. That covers the two frozen
  accumulations and the image power spectra, bootstrap map, residual pass
  and full reconstruction of `bootstrap_rec3D`.

### What an add-on then reads

- **Its own particles:** the sampled stage iterations, and three passes at
  the end (bootstrap map, residual sigmas, shipped map and carried sums).
- **Nothing of the frozen particles.**

The frozen and final-reconstruction passes over the union, up to six per
add-on, go. A session then reads each particle a bounded number of times,
however many add-ons follow.

### Contract changes

- **The manifest** records the carried sums and the merged sigma2 state. A
  further add-on takes them as its frozen input.
- **Entry** (section 5 of the add-on policy). The cohort starts from the
  frozen output's committed global sigma. The add-on no longer drops the
  inherited sigma2 registration or seeds the cohort from image power.
- **The carried sums** are validated against the run context before a byte
  is read, as the frozen sets are today (`simple_frozen_accum.f90`).
- **The frozen term** keeps the weighting it was accumulated with. The
  add-on policy's section 7 records the accepted difference.
- **The final reconstruction** (section 8 of the add-on policy) runs the
  bootstrap sequence over the run's own particles on top of the carried
  sums. It merges the run's new sigmas with the frozen particles' committed
  sigmas, and no longer bootstraps the union's sigma2.

### Reconstruction scans per state and half

The reconstruction's scans per state and half (section 3) are a separate,
smaller point. Where every state's halves fit in memory, as at cropped boxes,
one scan can feed them all.

## 6. An option that needs your ruling: grouping the sample by stack

The sampler takes the particles with the lowest update count first and
shuffles the ties at random. It could break ties by stack instead, for
example by taking the stack-ordered list of tied particles from a random
starting point.

- **What stays the same.** Update counts stay level, so every particle is
  still updated equally often.
- **What changes.** Each iteration reads about 170 whole stacks instead of a
  few images from each of about 1000. Each iteration's subset is then
  correlated by micrograph, and the rule competes with class balancing.
- **What it would need.** It changes the statistics of the sample, so it
  would need a test against ground truth before adoption.

## 7. Order of work

1. **Measure today.**
   - Split particle reading out of particle preparation in the benchmark
     files (below).
   - Run the ribosome with whole-stack runs on its spinning disk, with
     `free -g` and `iostat -x` alongside.
   - Run msp1 on fast storage.
   - Run the `SIMPLE_IO_NSTREAMS` sweep.
2. **Change the reader** in three measured steps: stacks kept open, one
   request per stack, read-ahead.
3. **Add-on passes:** frozen sums carried forward, the run's sigmas merged
   into the union's sigma2 state, and the incremental final reconstruction.

Each step ends with your review of its numbers.

### Benchmarks

- **Where.** The benchmark files: `REFINE2D_BENCH_ITERnnn.txt`, the refine3D
  matcher's per-partition files, and `SOLVE3D_BENCH.txt`.
- **Timings.** A particle-read timing, separate from particle preparation.
- **Counts,** in the context sections, for the pass: stacks opened, read
  requests, bytes read, bytes hinted.
- **Reconstruction.** Its read timing, kept today under `DEBUG` in
  `calc_3Drec`, moves to the benchmark flag.

## 8. Set aside

- **A raw single-stack cache.**
  - Copying the active particles into one file makes later reads cheap, but
    costs a full read and a full write of the data.
  - For the whole ribosome project's 751729 particle rows that is 30 to 55
    minutes on one spinning disk.
  - A stream cannot pay that, and the particles would exist twice.
- **Memory mapping.** It pays for one large shared file. With per-micrograph
  stacks, mapping each stack costs what opening it costs.

## 9. Where the change would go

| File | Change |
| --- | --- |
| `src/main/strategies/search/simple_matcher_ptcl_io.f90` | `discrete_read_imgbatch` on the open-stack pool and one request per stack; a read-ahead routine for a pass's particles; read counters and timing for the benchmarks |
| `src/fileio/simple_discrete_stack_io.f90` | the process-wide pool of open stacks with remembered header fields |
| `src/fileio/simple_imgfile.f90`, `src/main/image/simple_image_io.f90` | read a span of images in one request and convert the needed ones |
| `src/fileio/simple_posix.c`, `src/fileio/simple_syslib.f90` | the `posix_fadvise` helper, a no-op off Linux; the memory-available query; the rotational cap per the sweep |
| `simple_strategy2D_matcher.f90`, `simple_strategy3D_matcher.f90`, `simple_commanders_prob.f90`, `simple_matcher_3Drec.f90` | read-ahead at the start of a pass, or a window ahead in reconstruction; read timing in the benchmark files |
| `simple_commanders_solve3D.f90`, `simple_solve3D_utils.f90`, `simple_frozen_accum.f90` | frozen sums carried forward and recorded in the manifest; the cohort started from the union's committed global sigma, with no image power seeding |
| `src/main/simple_final_rec.f90`, `exec_bootstrap_rec3D` (`simple_commanders_refine3D.f90`) | the add-on's bootstrap sequence over the run's own particles on top of the carried sums |
| `src/main/sigma2/` | the output's sigma2 state merged from the frozen particles' committed sigmas and the run's |
| `doc/policies/3D/solve3D_addon_policy.md` | sections 3, 5, 7 and 8 for the carried sums, the merged sigma2 state and the incremental final reconstruction |
| `doc/policies/` | a short policy for the particle read path |

## 10. Questions for you

- **Grouping the sample by stack** (section 6): worth a ground-truth test?
