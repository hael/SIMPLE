# Particle I/O on per-micrograph stacks: reading, half-sets and add-on passes

Implementation note, 2026-10-12. **Plan; nothing implemented.**

Note (2026-10-14): `emulate_solve3D_stream` has been removed. Case F and its
command line in the run table need restating with `solve3D_chunks`.

**Reviewed 2026-10-13** (a static review, deleted once its changes were in
this plan; ruling of 2026-10-14) and revised the same day to meet it. The
review found the diagnosis sound but the plan not ready: two scientific
contracts and several designs were missing. Hans ruled on the contracts the same day (section 8). The designs
are in sections 3.1 to 3.3.

It merges two notes of the same day:

- the few-particle-stacks proposal (`665ac2731`), written against master
  `32a495cef`;
- `particle_io_handling.md` from the agent checkout, now retired.

RELION references are to its master at `05cd1ad9`.

**Sampling strategies are decoupled from this note** (Hans, 2026-10-12). The
goal is a data path that lets any sampling scheme read cheaply. Workflows
should not have to be designed around disk latency. Nothing here changes
which particles are drawn.

**When.** The global-sigma refactoring is on master (`2307ea96f`). The whole
plan is carried out by the unattended run of section 10 on the Oracle Linux
test machine, on beta-galactosidase and without timing claims (Hans,
2026-10-14). The production measurements of section 4 follow on the
ribosome project.

**The consolidated-stack experiment is done** (2026-10-13, section 4). One
file instead of 19,520 shortened the refine3D stages of the ribosome base
run by 18.7%. No material layout gain was demonstrated on this dataset and
machine. Hans's conclusion: consolidating stacks, in the stream or in batch
extraction, is not worth pursuing; the work is SIMPLE's low-level I/O
performance.

The plan was prompted by the stream ribosome project: thousands of
per-micrograph stacks with a handful of particles each, on which `solve2D`
and `solve3D` seemed to spend their time in file metadata. The experiment
did not show that. Companion documents:

- [particle_cache_policy.md](../../policies/2D/particle_cache_policy.md);
- [solve3D_policy.md](../../policies/3D/solve3D_policy.md);
- [solve3D_addon_policy.md](../../policies/3D/solve3D_addon_policy.md);
- the global-sigma refactoring outline (`global_sigma2_only.md`).

Rulings so far (Hans):

- **No particle is rewritten into another file** (2026-10-12). No cache, raw
  or cropped, is added. The existing cropped 2D cache and the extraction's
  pixel format stay as they are.
- **No layout with fewer files** (2026-10-13). Consolidating stacks, in the
  stream or in batch extraction, is set aside on the experiment's result
  (section 3.4). The work is the low-level read path (section 3.1).
- **Half-sets per stack** for new partitions (section 3.2). The review moved
  their creation from a per-field routine to one project-owned partition.
- **Add-on passes** (section 3.3):
  - frozen contributions are never recomputed, and their old noise weighting
    is an accepted error;
  - the add-on's final reconstruction is incremental;
  - the noise model is global;
  - new particles start from the union's committed global sigma;
  - a frozen particle can be retired in the project but stays in the frozen
    volumetric components already generated (2026-10-13);
  - add-on chains (`solve3D_addon`, `solve3D_chunks`,
    `emulate_solve3D_stream`) support only the gridding backend: they look
    for initial solutions and speed, not PCG-level precision (2026-10-13);
  - the first add-on of a chain reads the base particles once and writes the
    first carried sums; `solve3D` does not change (2026-10-13);
  - the offline chains publish their last incremental add-on as it is, with
    the old weighting of the frozen contributions (2026-10-13).

## 0. Summary

- **Fewer files bought little** (measured, section 4). One file instead of
  19,520 shortened the refine3D stages by 18.7%. That file was still reopened
  once per batch, and the experiment did not separate metadata, bytes,
  conversion, copying and compute. It shows direction, not an upper bound.
  What the per-file work is today:
  - Each sampled particle typically sits in a stack that no other sampled
    particle touches.
  - The first touch of a stack in a part process costs two existence checks,
    two opens, two 1 KB header reads and two closes, all before one particle
    is read.
  - The opens run one after another in the master thread, deliberately so
    since a June race fix. Only the particle reads run in parallel.
- **A distributed `solve3D` iteration in probabilistic mode reads its sample
  three times**, in two freshly launched processes: the probability table,
  then the search and the reconstruction.
  - Reconstruction reads once per state and half.
  - The half label alternates by particle index, so a stack is usually opened
    in both halves.
- **The stream adds a larger cost.** Every `solve3D_addon` reads all frozen
  particles up to six times: two frozen accumulations, and four passes of
  `bootstrap_rec3D` in the final reconstruction.
- **The plan, none of which stores anything derived from the particles on
  disk:**
  - **Read path, bit-exact:**
    - a pass-scoped read session, so each stack is opened once per pass;
    - bounded coalesced reads within a stack run;
    - a POSIX read-only path that opens stacks concurrently and prefetches
      a bounded window ahead of the reads;
    - a pass-sized particle buffer in 3D, so the reconstruction reads from
      memory what the search already read, within an explicit memory budget;
    - optionally, stacks kept open from the search into the reconstruction.
  - **Half-sets per stack** (ruled): one project-owned partition, created
    once in the master. Reconstruction opens each sampled stack once per
    state instead of once per half, and overlapping neighbouring picks stay
    in one half.
  - **Add-on passes** (ruled): an add-on reads only its own particles; the
    first add-on of a chain reads the base particles once.
  - **Set aside:** a layout on disk with fewer files (section 3.4).
- **Counts.** Counted from the code for two stream-sized examples
  (section 2.5), before any gain from concurrent opens or prefetching, the
  opens per `solve3D` iteration fall:
  - 1.6–1.7× with one open per stack per pass;
  - 1.7–1.8× adding half-sets per stack.

  Opens per iteration stay of the order of the number of sampled particles.
- **`solve2D` on a stream pool** below `nsample`=200,000 reads every particle
  every iteration. One open per stack per pass halves its opens, and
  concurrent opens and prefetching overlap them.

## 1. Why this is not a particle cache

**The removed 3D cache.** `doc/policies/2D/particle_cache_policy.md` says
only that caching "never worked reliably for 3D" and was removed from those
paths (`c5dae62b7`). The policy and the history describe a store of
noise-normalized, Fourier-cropped particles at `box_crop`, which:

- tied a run to one crop (`666e95960` fixed `solve3D` at the final ladder
  crop to avoid rebuilding at every stage);
- changed the reconstruction numerics (the removed 3D text: "the only
  numerics delta versus uncached execution is the edge taper acting at
  `box_crop`");
- needed a validity key, an owner, cleanup, a space budget, rank uniformity
  and a `cache_dir` every node can see.

**A raw single-stack copy** was considered on 2026-10-12 and set aside.
Copying the active particles into one file makes later reads cheap, but it
costs a full read and a full write of the data. For the whole ribosome
project's 751,729 particle rows that is 30 to 55 minutes on one spinning
disk. A stream cannot pay that, and the particles would exist twice.
Memory-mapping was set aside with it: it pays for one large shared file, and
mapping each per-micrograph stack costs what opening it costs.

**What this plan does instead.** No file is written: the same bytes reach the
same code at the same box, with no key, owner, budget or uniformity rule.

- The read-path changes alter how the bytes are fetched.
- Half-sets per stack alter the half labels of new partitions.
- The add-on changes alter which passes run.

## 2. What an iteration does to the file system today

### 2.1 Layout

Extraction writes one stack per micrograph (`ptcls_from_<mic>`,
`simple_extract_strategy.f90:360`, registered by `add_stk` at :383).
Particles reach their pixels through `stkind`/`indstk` in
`os_ptcl2D`/`os_ptcl3D` and the `os_stk` row (`map_ptcl_ind2stk_ind`,
`simple_sp_project_ptcl.f90:8`).

### 2.2 Passes per iteration (distributed, probabilistic mode)

`refine3D_distr_strategy%execute_iteration` (`simple_refine3D_strategy.f90:1187`):

1. **The probability table. Read 1.**
   - `exec_prob_align` (`simple_commanders_prob.f90:185`) draws the sample
     over the whole project in the master (:213-215), writes it to the
     project, and launches `prob_tab` parts (:228-231).
   - Each part reads its sampled particles (`build_batch_particles3D`, :93).
2. **The search. Read 2.**
   - `gen_scripts_and_schedule_jobs` (:1289) launches the `refine3D` parts.
   - Each part reproduces its share (`sample4update_reprod`) and reads it in
     batches of `nthr*BATCHTHRSZ` (`simple_strategy3D_matcher.f90:370`).
3. **The reconstruction. Read 3,** once per state-and-half group in which a
   stack has sampled particles.
   - The same process reconstructs (`write_partial_recs`, :475), through
     `calc_3Drec` (`simple_matcher_3Drec.f90:67-79`) or the PCG worker
     (`simple_rec3D_pcg_strategy.f90:291`, :403).
   - Either way it loops over states and halves and reads each group
     separately.
   - Under fractional state weights (`group_pinds_by_weights`,
     `simple_matcher_3Drec.f90:189`) a particle is read for every state it
     has weight in.
4. **The polish pass**, with `pose_cont=yes` (:1298-1306).
   - The discrete pass skips reconstruction (`volrec=no`, :1288).
   - A further set of parts reads the sample again and reconstructs.

`solve2D` per iteration: the `prob_tab2D` parts make read 1. The `refine2D`
parts make read 2, whose search batch also feeds the class sums.

### 2.3 What one pass does per stack

`discrete_read_imgbatch` (`simple_matcher_ptcl_io.f90:189`) takes the sorted
path (`read_sorted_stack_runs`, :280), because every sampler returns sorted
indices. Per stack run:

- **First touch: the dimensions.** `stk_dims` (:38, called at :319)
  memoizes dimensions per process. On first touch it calls
  `find_ldim_nptcls` → `get_mrcfile_info` (`simple_imghead.f90:1269-1277`):
  an existence check, an open, a 1 KB header read and a close.
- **The open.** `dstack_io%open` (`simple_discrete_stack_io.f90:94`) calls
  `imgfile%open` (`simple_imgfile.f90:92`): an open, another existence check
  (:141), and the header read again (:146).
- **The particles.** One `read(pos=)` per particle (`rSlices(i,i)` via
  `read_single_mrc_image`), then a close.
- **Threading.**
  - The open loop (:331) and the close loop (:343) run in the master thread.
  - Only the particle reads run in the `omp parallel do` (:333), on at most
    `io_read_nstreams` threads: 16 on local non-rotational disks, 8 on
    network file systems and macOS, 2 on spinning disks, or
    `SIMPLE_IO_NSTREAMS` (`simple_syslib.f90:436-484`).
  - The serial open is deliberate. `c1f7470ed` ("fixed parallel i/o race",
    2026-06-08) moved the open before the parallel loop and put a
    `critical(dstack_io_open)` around lazy opens inside parallel regions.
- **Cold memos.** In distributed execution the part processes are new for
  every pass, so the memo is cold at the start of every `prob_tab` and
  `refine3D` process. That makes **two opens per touched stack in reads 1
  and 2, and one per state-and-half group in read 3.** With m = 5 and about
  one sampled particle per touched stack, the particle read is a small
  fraction of the file operations around it.
- **What an operation costs.** From the page cache, tens of microseconds per
  read. From a cold spinning disk, a seek of about 10 ms per stack visit plus
  the file system's metadata on open. On SSDs and network file systems, the
  latency of each open and each request.

### 2.4 Sampling and half-sets

- **`solve3D`.**
  - `nsample` is 10,000 by default (`NSAMPLE_SOLVE3D_DEFAULT`,
    `simple_solve3D_controller.f90:55`).
  - `balance` defaults to `cavg` when the project has a 2D solution and no
    `vol1`, to `class` with `vol1`, and to `none` without a 2D solution
    (`set_balance_default`, `simple_commanders_solve3D.f90:1833`). The
    stream's 3D stage now balances by default too.
  - On stream data the sampler is therefore the nested quota of
    `sample4update_class` (`simple_oris_sampling.f90:242`). `balance=none`
    uses `sample4update_cnt` (:183).
  - Both exhaust the lower update-count tiers and draw the cutoff tier
    uniformly over particles.
  - The result comes out in index order (`hpsort` in `_cnt`, packed in index
    order in `_class`), so reads already go stack by stack. What is not
    coupled to the stacks is which particles are drawn.
- **`solve2D`.** `nsample` is 200,000 (`NSAMPLE_DEFAULT_2D`,
  `simple_solve2D_controller.f90:26`). Below that the update fraction is 1
  and every particle is read every iteration.
- **Half-sets.** `partition_eo` (`simple_oris_setters.f90:791`) alternates
  `eo` by particle index, so any stack with two or more consecutive sampled
  particles has both halves.

### 2.5 Counts per iteration

S stacks of m particles, sample n, f = n/(Sm).

- Stacks touched per pass: T = S(1-(1-f)^m).
- Reconstruction opens: R = S[(1-(1-f)^e) + (1-(1-f)^o)], with e and o the
  even and odd counts in a stack.

Distributed, probabilistic mode, one state, no polish; opens are summed over
all parts. Each strategy is named in the table's first column.

| | A: 5,000 stacks × 6 = 30k ptcls | B: 40,000 stacks × 5 = 200k ptcls |
| --- | --- | --- |
| f (n = 10,000) | 0.33 | 0.05 |
| stacks touched per pass, today | 4,561 (91%), 2.2 ptcl each | 9,049 (23%), 1.1 ptcl each |
| opens per `solve3D` iteration, today (2T + 2T + R) | 25,281 | 45,800 |
| one open per stack (T + T + R) | 16,159 | 27,703 |
| ... + half-sets per stack (T + T + T) | 13,683 | 27,147 |
| particle read calls today (3n) | 30,000 | 30,000 |
| with at most one request per stack run (coalesced) | 16,159 | 27,703 |
| ... + half-sets per stack | 13,683 | 27,147 |
| `solve2D`, full sampling: opens today → one open per stack | 20,000 → 10,000 | 160,000 → 80,000 |

Some cases add to these counts:

- More states add reconstruction opens wherever a stack's sampled particles
  fall in several states.
- The polish pass adds one pass.
- Each open also carries an existence check and a header read.

These are counts from reading the code, not measurements. The file-operation
measurement of section 4 checks them.

### 2.6 The add-on passes

Every `solve3D_addon` run reads all frozen particles, the particles aligned
by earlier runs, several times.

- **The frozen accumulations** (`calc_frozen_rec`).
  - The add-on runs one stage, the last one its frozen run ran.
  - It accumulates the frozen particles at that stage's box and at the native
    box: two passes, or one when the two boxes are equal.
- **The final reconstruction** (`calc_final_rec` → `bootstrap_rec3D`).
  - It re-fits the noise model of the whole union.
  - That reads every particle four times: the image power spectra, the
    bootstrap map, the residual sigma2 pass and the shipped map.
- **The size of a pass.** One pass over 500,000 particles at box 336 reads
  about 226 GB of float32: 20 to 25 minutes on one spinning disk.
- **Over a session.** An add-on starts each time the new particles reach 10%
  of the frozen ones. A session growing from 100,000 to 750,000 particles
  therefore reads its final dataset about fifty times. This is arithmetic,
  not a measurement.

## 3. Plan

The read path (3.1) is the work: SIMPLE's low-level I/O performance (Hans,
2026-10-13). The measurements of section 4 set the order of its steps. The
add-on passes depend on the global-sigma refactoring in any case.

### 3.1 Read path (bit-exact)

These change only how bytes are fetched: `simple_matcher_ptcl_io.f90`, the
stack reader (`simple_discrete_stack_io.f90`) and the image file layer. The
cropped 2D cache keeps its own reader (`ptcl_cache_read_batch`). Its build
reads through `discrete_read_imgbatch`, so the build gains too.

**One open per stack per pass: a pass-scoped read session.**

- **The duplicate first-touch work goes.**
  - Take `ldim` and `nptcls` from the header that `dstack_io%open` reads
    anyway, and drop the separate `find_ldim_nptcls` on first touch.
  - The rule in the memo comment is kept: the physical header, not `os_stk`'s
    `box`/`nptcls_stk`, decides record offsets. The check moves into the
    open.
  - Drop the existence check after the open in `imgfile%open` (:141) for
    read-only opens. It is a second metadata operation on every open.
- **That alone gives one open per batch call, not per pass.**
  `discrete_read_imgbatch` builds and destroys its readers on every call, and
  the probability table, the search and the reconstruction call it batch
  after batch. A stack split across a batch boundary is opened again.
- **The session.** A read session owns the state of one pass:
  - the decoded headers;
  - reusable per-stream scratch, including the mode-12 conversion buffer,
    instead of the small reader and header arrays allocated per stack today;
  - exact counters (stacks opened, requests, bytes);
  - at least the descriptor of the stack at the current batch boundary, so
    the next batch continues it.

  Batches that end on stack boundaries are an alternative where they keep
  within the batch memory bound.
- **No change to threading.** Validation counts opens over a complete pass,
  not one reader call.

**Coalesced reads within a stack run.**

- Group a run's needed slices into one or more spans, each read with one
  positioned request into scratch, and copy the needed slices into
  `build%imgbatch`.
- A span bridges a gap only when that costs less than a separate request on
  the already-open descriptor. The decision weighs:
  - the bytes read beyond the needed slices;
  - a hard cap on scratch;
  - the conversion and copy cost of the extra slices, high for mode 12;
  - the measured cost of one more positioned read on an open descriptor.
    Open latency is not the right figure: the saved operation is a request,
    not an open.
- One span from the lowest to the highest index is not the rule. Slices 1
  and 5 of a five-image stack are two slices today, and one span would read
  and convert all five.
- On a spinning disk, passing over a short gap costs about as much as
  reading it, while a second request costs a rotation or a seek.
- This is second priority on stream data (one or two particles per stack per
  pass today), and first-order on large stacks and full sampling.

**A POSIX read-only path: concurrent opens and a bounded prefetch window.**

- **The path.** One small set of C helpers in `simple_posix.c`, bound in
  `simple_syslib.f90`, with no Fortran unit involved:
  - open read-only;
  - `pread`;
  - `posix_fadvise`;
  - close.
- **What the path must define**, in one place:
  - one authoritative MRC header decode;
  - 64-bit data offsets;
  - extended headers (`nsymbt`);
  - mode and signedness conversion;
  - short reads and `EINTR` retried;
  - end of file and truncation reported as errors;
  - the fallback to today's reader for a mode the path does not handle.
- **Concurrent opens.** Each stream opens, reads and closes its own stacks
  inside the parallel region, so up to `io_read_nstreams` opens are in flight
  instead of one.
  - The path touches no Fortran unit, so it does not depend on what the June
    race was.
  - `simple_stack_io_tester.f90`, which already reads in parallel
    (`test_dstack_parallel_read_*`), gains a case that opens in parallel,
    many times over, on more stacks than streams.
- **Prefetch in a bounded window.** The sample is known before the first
  batch, so reads can be issued ahead of use. `posix_fadvise` takes an open
  descriptor, not a path, so hinting a whole pass up front would mean
  opening every file twice or holding thousands of descriptors.
  - Instead a producer keeps a bounded window of stacks ahead of the
    consumer: it opens them and fills a double buffer with `pread`, and
    `WILLNEED` may be issued inside the window as a hint.
  - The advice is not binding, and the kernel is not relied on to order
    reads across files.
  - **Matchers and probability tables** run the window over the sorted
    sample. **Reconstruction** keeps the next batch or two ahead while the
    current one is gridded.
  - **Bounds.** The window has a fixed byte budget, so prefetched data never
    evict what the run already holds.
  - **Full sampling.** This also serves `solve2D` on stream pools, where
    sampling cannot help.
- **Conversion.** It is the existing code: `simple_float16_to_float32` for
  mode 12 and the current conversions for the other modes.
- **Platform.**
  - `posix_fadvise` is Linux-only and is a no-op elsewhere, as
    `simple_fs_is_rotational` already is.
  - On network file systems the effect of the hint depends on the file
    system.
  - The Windows build is not a concern.

**A pass-sized particle buffer: read once, use twice in 3D.**

- **Today** the 3D search reads its share in batches of 50 images per thread,
  each read, prepared and searched before the next (`refine3D_exec`,
  `simple_strategy3D_matcher.f90:170`).
  - The reconstruction then runs in the same process (`write_partial_recs`,
    :225, into `calc_3Drec` or the PCG worker). It reads the same particles
    from disk again, 500 at a time, once per state and half, and once per
    state with weight under fractional state weights.
  - 2D already reads once: the search batch also feeds the class sums.
- **The buffer.** Before the search loop the read layer fills one buffer
  with the images of the pass's share, in one sorted sweep through the read
  session.
  - **Representation.** One contiguous slab of decoded `real32` pixels plus
    an index map from particle to slab position. It is immutable once
    filled. Not an array of `image` objects: search preparation changes its
    input pixels, and every `image%new` allocates FFT storage and creates
    FFT plans.
  - **Use.** The search batches and the reconstruction copy their images
    from the slab into the existing bounded image batches. The results stay
    bit-identical. Read 3 disappears.
  - **Mode 12.** The slab holds the decoded float32 values, so conversion
    runs once. The memory figures below assume that.
- **Memory.** At box 336 in float32 a particle is about 450 KB.
  - A distributed part with 1,250 particles (10,000 sampled over eight parts)
    needs about 560 MB.
  - A shared-memory run with 10,000 particles needs about 4.5 GB.
- **The memory budget.** The buffer is allocated against an explicit budget,
  not against the memory the system reports available: several workers
  reading that figure at once would each claim the same memory.
  - The budget comes from the scheduler or the user where given, and
    otherwise from a conservative share per process or per node.
  - The allocation checks its status. A share that does not fit the budget,
    or an allocation that fails, falls back to today's batch reads, and the
    pool below then applies.
  - Benchmarks report the aggregate RSS of the node as well as each
    process's.
- **Scope.**
  - The buffer lives for one iteration in one process, and nothing is
    written to disk.
  - The probability table and the polish pass run in other processes and
    keep their own reads.
- **Overlap with compute** comes from the prefetch window while the buffer
  fills, not from the buffer itself. Filling a buffer before the loop on its own
  changes the order of the reads, not their number or cost.
- **The memory policy.** The second read exists because of
  `separate_alignment_and_reconstruction_for_multistate_peak_mem_reduction.md`
  (July 2026, `7a979b086` and `90d4ed283`). Inserting each batch into the
  reconstructors during the search kept every state's even and odd
  reconstructor allocated next to the all-state reprojection model. The
  policy separated the phases to bound the peak at the larger of the two, and
  accepted the reread as the price.
  - The buffer keeps the phases separate but raises both peaks. With A the
    alignment working set, R the reconstruction working set and B the
    buffer, the peak goes from max(A, R) to max(A + B, R + B). B is small
    next to A and R at the sizes above, but it is not free, hence the budget.
  - The policy's lifetime table does not list a buffer living across the
    assignment barrier. It gains that line, with the budget rule, when the
    buffer lands.
- **Measure** with the time-share measurement of section 4, which times
  reading in the search and in the reconstruction separately.

**Stacks kept open from the search into the reconstruction (optional).**

- In a `refine3D` part, reads 2 and 3 run in the same process. A pool of the
  part's open descriptors, kept from the search into the reconstruction,
  makes read 3's opens free.
- With the particle buffer this matters only for shares that exceed the
  buffer's budget.
- It is bounded by the open-file limit. When the pool is full, the least
  recently used stack is closed.
- A bounded LRU pool may reuse little: the search ends holding the tail of
  its order, while the reconstruction starts again near the beginning.
  Before any code, the real access trace of a pass is replayed against
  realistic pool sizes.
- It is worth adding only if that replay and the file-operation measurement
  show the metadata cost still matters after the steps above.

**The stream-count cap.**

- The cap of 2 on spinning disks was set in August 2026 to protect long
  sequential reads from head thrashing. It was not measured.
- Scattered reads of a few images per stack are the opposite case, and
  arrays of spinning disks need more requests in flight than 2.
- `SIMPLE_IO_NSTREAMS` overrides the cap already, so a sweep settles it
  before code changes (section 4).

### 3.2 Half-sets per stack (ruled 2026-10-12)

- **Two reasons:**
  - With per-index alternation, a stack with two or more sampled particles
    usually has both halves, so reconstruction opens it once per half. Per
    stack, it opens each touched stack once per state.
    - The gain grows with the number of sampled particles per stack.
    - In example A of section 2.5, reconstruction opens fall from 7,037 to
      4,561.
    - In example B, with about one sampled particle per stack, they fall
      only from 9,604 to 9,049.
  - Overlapping neighbouring picks share pixels. With per-index alternation
    they sit in opposite halves, which couples the half-maps and inflates
    the FSC. Keeping a micrograph's particles in one half removes that
    coupling. Overlap is also why RELION keeps all segments of a helical
    tube in one half-set (`divideParticlesInRandomHalves`,
    `exp_model.cpp:261-317`; its error text at :400: "half-sets are always
    per-filament").
- **The partition belongs to the project, not to a field.** Today
  `partition_eo` (`simple_oris_setters.f90:791`) labels one field wherever it
  is called. A randomized per-stack version of it cannot stay that way:
  - SIMPLE's particle sampling uses the same process-wide random state, so
    drawing a split would change the next sample, against this note's
    promise to change no draw.
  - Distributed call sites can each create a partition while only one
    result is kept.
  - `ptcl2D` and `ptcl3D` can be partitioned separately and disagree once
    the stack order is random.
  - In a partly labelled field, callers that test only
    `get_nevenodd() == 0` never label the unlabelled appended rows.
  - Stream merges and the final ingestion combine labelled and unlabelled
    rows.
- **One transaction, in the master, before any worker runs:**
  1. validate the grouping provenance and how complete the current
     partition is;
  2. derive one deterministic split from an isolated random stream;
  3. balance the active (`state > 0`) populations while every row of a
     group gets one label;
  4. write the same labels at once to `ptcl2D` and `ptcl3D`, which must be
     congruent;
  5. record enough provenance to tell a complete, a missing and an
     inconsistent split apart;
  6. workers consume the partition, and stop rather than invent one.

  The call sites that run `partition_eo` today (the refine2D and refine3D
  strategies, the 2D and 3D matchers, reconstruction, `calc_pspec`,
  `make_cavgs` and `solve2D`) become consumers. Existing labels never
  change; unlabelled rows of a partly labelled project are placed by the
  same rule. In the stream, chunk projects are partitioned as they are
  created, and the pool keeps their labels.
- **Groups are scientific, not storage.** `stkind` is a storage row. Native
  extraction makes it one micrograph, but imports can lay particles out
  otherwise, and several stack rows can share one physical file
  (`consolidate_particles`). The grouping follows the micrograph or
  filament identity where the project has it, not the physical path:
  - a stack row is a group when it maps to a micrograph row
    (`get_mic2stk_inds`), as native extraction writes;
  - a field in which an active stack row has no micrograph row (imported
    particle stacks, such as the beta-galactosidase test set) has no known
    scientific grouping and takes the per-index route, as today. The route
    is logged. This is a rule of provenance, not of size.
- **Placement.** Largest group first, to the half with fewer active
  particles, with deterministic random tie-breaking. The imbalance is
  bounded by the largest group. It is logged, and the transaction warns, or
  refuses, when there are too few independent groups.
- **No split by size.** A large micrograph or filament stays whole: its
  overlapping particles are exactly what the rule is for. A finer split
  follows only an explicit scientific identity, such as a filament. This
  replaces the earlier fallback to per-index labels for any stack above 1%
  of the field.
- **Routes without micrograph groups.** Class-average fields (no `stkind`)
  and time series (`ini_tseries`, in `make_cavgs`) take an explicit policy
  argument, or a separate per-index routine. `oris%partition_eo` alone has no
  workflow context to infer them.
- **FLEX.** Its repair of a degenerate split
  (`simple_flex_pca_project_gateway.f90:104`) uses the same project-level
  transaction, so it follows the same rule.
- **Not chosen.** Reading each batch once and routing it to both halves
  would avoid the label change. But it needs both half accumulators live,
  which the multistate peak-memory policy avoids.

### 3.3 Add-on passes (ruled 2026-10-12)

This is the stream's largest particle I/O. It changes the `solve3D_addon`
contract (`solve3D_addon_policy.md`, sections 3, 5, 7 and 8), and it depends
on the global-sigma refactoring (`global_sigma2_only.md`). The review found
four contracts the first rulings did not settle; all four are now ruled
(section 8).

**Why the frozen sums are recomputed today.** The frozen particles' poses and
labels never change during a session (stream 3D ingestion policy, item 10).
Their noise weighting does:

- The frozen term is accumulated with the frozen run's committed sigma2
  state.
- Every add-on's final reconstruction fits a new state to the union, and
  that state becomes the next add-on's frozen sigma2.

**Frozen sums carried forward.**

- The frozen sums at the stage's box and the native box are kept as they
  were accumulated. They are never rescaled or recomputed when sigmas
  change.
- The next add-on's frozen sums are this add-on's plus this run's particles
  at their final poses, with their best sigmas (the run's last pass, below).
- The boxes are the same in every add-on of a chain, since each runs the
  frozen run's last stage.
- **Retirement (ruled 2026-10-13).** Today a frozen particle rejected in the
  current project's `ptcl2D` retires: it leaves the frozen term and stays
  inactive in the union (`solve3D_addon_policy.md`, frozen membership;
  `simple_project_superset.f90`, `retire_from_frozen`). Once a particle is
  inside a carried sum it cannot be taken out, and that is accepted: a
  retired particle is inactive in the project but stays in the frozen
  volumetric components already generated. The add-on policy's
  frozen-membership text changes accordingly.
- **The first carried sums (ruled 2026-10-13).** The base run publishes no
  carry-ready sums, and its finished maps and fractional trailing
  accumulators are not the full-particle raw sums required. The first add-on
  of a chain therefore reads the base particles once, at the stage's box and
  the native box, as `addon_starting_state` (`simple_commanders_solve3D.f90`)
  does today, and writes the first carried sums. `solve3D` does not change.
- **Gridding only (ruled 2026-10-13).** Add-on chains look for initial
  solutions and speed, not PCG-level precision.
  - `solve3D_addon` replays the base run's backend today and offers PCG
    controls (`maxits_pcg`, `maxits_ml`, `pcg_solvent_check`). It drops
    them and refuses a frozen project whose base run used PCG, naming the
    base run's backend.
  - `emulate_solve3D_stream` drops its PCG inputs
    (`add_solve3D_base_inputs` with `l_include_pcg=.false.`) and forces
    gridding, as `solve3D_chunks` already does.
  - The PCG side of the frozen sets (the `B`/`D` payloads in
    `simple_frozen_accum.f90` and the frozen add in the PCG backend) is then
    dead and is deleted. A base `solve3D` keeps both backends.
  - This step does not wait for the global-sigma refactoring. It also
    removes the question of how a PCG chain would seed `bootstrap_rec3D`,
    which forces its bootstrap map through gridding.

**The noise model.**

- It is global: each particle's spectrum is kept in the sigma2 state, and
  alignment and reconstruction use one global sigma, their average.
- Estimating the noise is part of refinement. A better estimate of a
  particle's sigma is used as soon as it exists.
- New particles start from the union's committed global sigma: each new
  particle's spectrum starts as that sigma. The cohort's spectra are then
  updated from their residuals during the stage.
- The add-on no longer drops the inherited registration or seeds the cohort
  from image power (`calc_pspec` through `ensure_sigma2_for_iteration`).
  Only the base run, with no union before it, seeds from image power.
- Frozen particles keep their committed sigmas, and their contributions keep
  the weighting they were accumulated with.
- **The transaction is copy-on-write.** The inherited committed state belongs
  to the frozen output and must stay intact for rollback, and its rows do
  not describe the growing union. The add-on therefore:
  1. creates a candidate state of its own;
  2. copies the retained frozen records into it;
  3. initializes every cohort row from the right even or odd global curve;
  4. updates only the cohort during the stage;
  5. restores the union membership;
  6. recomputes both global means over the union;
  7. publishes atomically.

  The cohort comes from the membership mask, not from a suffix of appended
  rows: never-updated rows and reselected rows inside the frozen range can
  join it. Means committed while the frozen rows are still masked would be
  cohort-only means.

**Incremental final reconstruction.** Today's bootstrap sequence runs over
the run's own particles instead of the union, with no image-power seed:

1. **Bootstrap map:** the carried frozen sums at the native box, plus the
   run's particles at their final poses with their sigmas from the stage.
2. **Residual pass:** the run's particles against the bootstrap map at native
   sampling. Their new sigmas are committed.
3. **Shipped map and new carried sums:** the run's particles accumulated
   again with the new sigmas, at the native box and the stage's box.
   - The shipped map is the carried frozen sums at the native box plus
     these.
   - The next add-on's frozen sums are the carried sums plus these.

**Then:**

- **The output** carries the sums and a sigma2 state that merges the frozen
  particles' committed sigmas with the run's new ones.
- **The full fit.** Only the stream's p07 stage ends with a full `refine3D`
  over every selected particle (stream 3D ingestion policy). The offline
  chains do not: `emulate_solve3D_stream` ends at its last adopted add-on,
  and `solve3D_chunks` publishes its last add-on directly. They publish the
  incremental result as it is, and the historical weighting of the frozen
  contributions stays in offline results (ruled 2026-10-13: the fastest
  option). No terminal full reconstruction is added.
- **What an add-on reads:** its own particles only, in the sampled stage
  iterations and at the end. Nothing of the frozen particles, once valid
  carried sums exist. One reconstruction call owns one grid, so the passes
  at the end are:
  - three when the native box and the stage's box are equal: bootstrap map,
    residual pass, shipped map with the new carried sums;
  - four when they differ, since step 3 then reads once per box, unless a
    dual-grid accumulator or the cohort buffer serves both;
  - one more with `addon_diag=yes`, the default of `solve3D_chunks`, for the
    cohort-only reconstruction.
- **Over a session,** each particle is read a bounded number of times,
  however many add-ons follow.
- **The carried sums are one bundle, validated cheaply.** They must not go
  through the manifest's generic digest of non-sigma artifacts, which reads
  the whole file. As with today's frozen sets (`simple_frozen_accum.f90`):
  - components are written under temporary names;
  - a bundle manifest records the schema, dimensions, sampling, backend,
    state and half, byte sizes, the membership or lineage identity, the
    sigma2 generation and the parent and run lineage;
  - the bundle manifest is published last, and the solve manifest and the
    project only after the bundle is complete;
  - validation checks headers, sizes, context, lineage and generation, not
    full-file digests. The sigma2 identity keeps its structural and
    generation signature.

  Each adopted run owns its complete bundle, because stream cleanup can
  remove earlier run directories. A failed or rolled-back run leaves the
  previous bundle untouched.

### 3.4 A layout on disk with fewer files (set aside, 2026-10-13)

Extraction writing containers, large files that hold many micrographs'
stacks, was proposed so that any sampling would touch few files. In the
consolidated-stack experiment (section 4) one file shortened the refine3D
stages of the ribosome base run by 18.7% and gained nothing on passes that
read in order: no material layout gain on this dataset and machine. Hans set
it aside on that result, together with consolidating stacks in the stream or
in batch extraction.

`consolidate_particles` stays as a developer program for measurements, and
with it the paragraph "Stack Rows Sharing One File" of
`project_stack_indexing_policy.md`.

The prior art below is kept for reference.

**Prior art: cisTEM** (its master at `455de9a1`). cisTEM never has
per-micrograph stack files.

- **One stack per refinement package.**
  - The stack (`particle_stack_N.mrc`) is cut straight from the micrographs
    when the package is created: each micrograph is read once, and its
    particles are written consecutively (`MyNewRefinementPackageWizard.cpp`).
  - A package made from a subset of another one copies the selected
    particles into a new stack, unless all of them carry over.
- **Reading is a forward scan of one file.**
  - `MRCFile` (`src/core/mrc_file.cpp`) is one `std::fstream` that seeks
    only when it is not already at the wanted slice.
  - Reads are serialized under `omp critical` within a process.
  - Distributed jobs take contiguous particle ranges
    (`FirstLastParticleForJob`).
  - There is no mmap, `pread` or read-ahead.
- **Subsets are chosen in stack order.**
  - `refine2d` decides per particle before the read.
  - `refine3d` reads every particle of its range and decides afterwards, so
    a 10% update still reads the whole range.
  - The noise-spectrum estimate and `refine2d`'s starting averages sample
    random blocks of consecutive particles "to avoid seeking a lot in large
    files".
- **What it shows:**
  - With one file per package the per-file cost disappears.
  - Writing that file costs no extra pass, because it is written at the moment
    particles are cut from the micrographs. That is the moment SIMPLE's
    extraction would write containers.
  - Arbitrary sampling is still not free. Where seeks mattered, cisTEM shaped
    the sampling to the layout (blocks) or read everything.

**Prior art: cryoSPARC.** It is closed source. This comes from its public
guide and forum posts, which were not read in full.

- **Same indexing as SIMPLE.** Extraction writes one stack per micrograph.
  Each particle's metadata row holds a file path and an index in that file.
  A selection removes rows and rewrites no file.
- **Same cost.** Particles stay in many small files, and a subset still goes
  through the original stacks. No timings are published.
- **Two remedies.**
  - An SSD cache copies whole stack files to a local SSD at job start. Jobs of
    a project share it, and old data is evicted when space runs out.
  - The Restack Particles job writes the selected particles into fewer, larger
    files, with the number of particles per file as a parameter. It is a
    manual step. Users report it as slow, memory-hungry with large batches,
    and prone to open-file errors.
- **Reading.** Refinement loads as many images onto the GPU as fit and
  normally runs two threads, so reading and computing overlap. Particles
  come from the cache or from memory. No prefetch setting is documented.
- **What it shows.** Consolidation is the remedy users reach for. On the
  ribosome project it bought 18.7% of the stage time (section 4).

### 3.5 Not addressed

Every part process re-reads the project file and references at launch, and
the master writes the project segment each iteration. Persistent workers
would remove that, and the cold memo with it. That is out of scope here.

## 4. Measurements: the basis for the next decision

The benchmark split of section 5 and a regression baseline on
beta-galactosidase are Phase 0 of the unattended run (section 10). The
measurements below that need the production machine and the ribosome project
are the maintainer's, after the run:

- the open latency on the cluster file system and the positioned-read
  latency on the spinning disk, which set the production gap budget for
  coalesced spans;
- the access-trace replay, which decides whether stacks are kept open from
  the search into the reconstruction;
- the stream-count sweep, which decides whether the cap on spinning disks
  changes;
- the whole-stack ribosome run.

The beta-galactosidase set on the test machine is one stack, so none of
these can be measured there.

**The consolidated-stack experiment (done, 2026-10-13).**

- **Setup.** `consolidate_particles` copied the ribosome project's 436,455
  active particles (197 GB in 19,520 stacks of about 22) byte for byte into
  one file on the same disk. The copy is disk-bound, about 53 minutes on the
  production server: one spinning disk holding the stacks and the outputs,
  and 252 GiB of memory. `solve3D_chunks` then ran with identical settings
  and the same 15 runs of about 29,100 particles on the original project and
  on the consolidated one.
- **Result.**

  | | 19,520 stacks | One file | Change |
  | --- | --- | --- | --- |
  | base step, total | 8,638 s | 7,103 s | −17.8% |
  | refine3D stages | 8,011 s | 6,514 s | −18.7% |
  | final reconstruction | 413 s | 398 s | −3.6% |
  | add-ons 1–3, s per 1,000 particles added | 147, 141, 148 | 125, 129, 134 | −9% to −15% |

- **What it shows.** No material layout gain was demonstrated on this
  dataset and machine.
  - The gain is in the stage iterations, which read a random sample. The
    final reconstruction, which reads in order, gained nothing.
  - It is direction, not an upper bound on what layout or metadata work can
    buy. The one file was still reopened once per batch, and the experiment
    does not separate metadata, bytes, conversion, copying and compute.
    About half of the wall time of a matcher partition is outside its
    timers.
- **Caveats.** Ab initio is stochastic: the two base runs reached different
  three-state solutions, and compute differs with the solution, so the
  numbers are indicative. One machine, and the consolidated chain had
  finished 4 of its 14 steps.
- **A firm attribution**, if one is ever needed: capture one pass's selected
  particles and access trace and replay it against both layouts, or repeat
  fixed-sample runs to completion.
- **Consequence** (Hans): no layout with fewer files (section 3.4); the work
  is the read path (section 3.1).

The measurements:

- **File operations.** Run
  `strace -f -c -e trace=openat,newfstatat,statx,read,pread64,lseek,close` on
  one `prob_tab` part and one `refine3D` part of one iteration. Compare with
  2T, 2T and R for that part's range.
- **Time share.**
  - Particle reading gets its own timing in the benchmark files, separate
    from particle preparation (section 5). The read timer of `calc_3Drec`
    (`t_read`, behind `DEBUG`) moves to the benchmark flag.
  - Read the shares on the ribosome project for `solve3D` and `solve2D`.
- **Open latency** on the cluster file system. Time open, a 1 KB read and a
  close over a few thousand stacks not read recently, with one stream and
  with 8 streams. This says how much concurrent opens buy.
- **Positioned-read latency** on an already-open descriptor, cold and warm.
  This, not the open latency, sets the gap budget for coalesced spans.
- **The access trace of one pass**, replayed against realistic descriptor
  pool sizes before any pool is written (section 3.1).
- **Stream-count sweep.** `SIMPLE_IO_NSTREAMS` at 1, 2, 4, 8 and 16 on the
  ribosome's spinning disk, with `free -g` and `iostat -x` alongside.
- **The ribosome with whole-stack runs** (`solve3D_chunks`, `76c37002b`).
  Does the page cache hold a run's particles after their first read?

## 5. Benchmarks

- **Where:** `REFINE2D_BENCH_ITERnnn.txt`, the refine3D matcher's
  per-partition files, the refine3D strategy file, `SOLVE3D_BENCH.txt`, and
  a per-partition file of the probability-table parts (`prob_tab`,
  `prob_tab2D`), which make the first read of every iteration and have no
  benchmark file today; all in SIMPLE's format (CONTEXT, TIMINGS,
  % accounted for).
- **Timings:** particle reading, separate from particle preparation;
  conversion and copy time; buffer allocation time.
- **Counts, in the context sections, for each pass:**
  - stacks opened;
  - read requests;
  - bytes requested, and bytes actually read including over-read;
  - bytes prefetched, and how many of them were used;
  - the cache state, warm or cold.
- **Memory:** process RSS and the aggregate RSS of the node.
- **Rule:** a gain is claimed only from the same run measured before and
  after.

## 6. Tests

- **The particle buffer.**
  - With the buffer and without it, the search and the reconstruction
    produce bit-identical outputs.
  - The slab is unchanged after the search; fractional multi-state
    reconstruction reuses it once per state.
  - A share above the budget, or a failed allocation, falls back to batch
    reads, with the same outputs.
  - The reconstruction makes no read from disk when the buffer holds the
    share, checked by the open counter.
- **The read path.**
  - `build%imgbatch` is bit-identical between the old and new paths on a
    fixture with stacks of 1, 5 and 200 particles, on sorted, unsorted and
    duplicate index lists, in mode 2 and mode 12, with the output order kept.
  - One open per stack per pass, counted over a complete pass, including a
    stack split across consecutive batches.
  - Two logical stack rows sharing one physical file.
  - Several coalesced spans in one stack, with the over-read within its
    budget.
  - A dimension mismatch still stops.
  - The POSIX path: nonzero extended headers, 64-bit offsets, missing and
    truncated files, short reads, `EINTR`, end of file, advice the platform
    does not support, and the fallback modes.
  - The POSIX path opens in parallel, many times over, on more stacks than
    streams. Descriptor exhaustion and leaks are covered (pool eviction too,
    once a descriptor pool exists).
  - Prefetching changes no byte read.
- **Half-sets per stack.**
  - The halves differ by at most the largest group.
  - Every particle of a group carries one label, and no group is split for
    its size.
  - Class-average fields and time series take their explicit route.
  - `ptcl2D` and `ptcl3D` carry the same labels after merge, prune,
    reextraction and stream publication.
  - A partly labelled appended project gets its missing labels; existing
    labels are untouched.
  - Workers refuse a missing or inconsistent partition instead of creating
    one.
  - Creating the split leaves the particle sample unchanged (isolated random
    stream), and a fixed seed gives the same split.
- **Add-on passes.**
  - When no sigma changes, the carried sums equal sums recomputed from the
    frozen particles.
  - The shipped map from carried sums equals a full reconstruction of the
    union under the same weights, within tolerances justified by the changed
    accumulation order; drift is checked over the depth of a chain.
  - At least two chained add-ons, with rollback, restart and truncated or
    mixed bundles.
  - A frozen particle retired after at least one carried add-on stays in the
    carried sums and is inactive in the project; reselection after
    retirement.
  - Gridding, with equal and unequal native and stage boxes;
    `solve3D_addon` refuses `rec_backend=pcg`.
  - The first add-on of a chain reads the base particles once and writes
    the first bundle.
  - No frozen particle stack is read once a valid bundle exists.
  - Both offline chains publish their last incremental add-on as it is, with
    no terminal reconstruction.
  - The union's even and odd sigma means after the membership is restored.
  - A bundle of another run context is refused.
  - The merged sigma2 state holds the frozen particles' committed records
    unchanged.
- **Lightweight gates.** `scripts/check_test_registry.py` runs whenever a
  suite table or the test UI changes.

## 7. Order of work

The global-sigma refactoring, the prerequisite of the add-on track, is on
master (`2307ea96f`). The rest is carried out by the unattended run of
section 10, in this order (Hans, 2026-10-14: everything in one run, every
read-path step implemented and kept unless it is slower):

1. **The benchmark split** and a regression baseline on beta-galactosidase.
2. **The read path:** the pass-scoped read session and bounded coalesced
   reads; then the POSIX read-only path with concurrent opens and the
   bounded prefetch window; then the pass-sized particle buffer in 3D, with
   its budget and the policy line.
3. **Half-sets per stack:** the project-owned partition transaction and its
   consumers.
4. **PCG out of the add-on chains.**
5. **The add-on passes.**
6. **Documents and close.**

Each step ends with the regression cases of section 10 against the
baseline. Stacks kept open from the search into the reconstruction are not
in the run: they wait for the access-trace replay on the ribosome project
(section 4).

## 8. Decisions

None is open. The review's contracts were ruled by Hans on 2026-10-13:

- **Retirement.** A frozen particle can be retired in the project but stays
  in the frozen volumetric components already generated (section 3.3).
- **PCG.** Add-on chains support only the gridding backend (section 3.3).
- **The first carried sums.** The first add-on of a chain reads the base
  particles once, at the stage's box and the native box, as every add-on
  does today, and writes the first carried sums. `solve3D` does not change.
  The alternative, the base run writing the sums, would read the same
  particles and touch more code.
- **What the offline chains publish.** `solve3D_chunks` and
  `emulate_solve3D_stream` publish their last incremental add-on as it is,
  with the old noise weighting of the frozen contributions: the fastest
  option. In the stream, p07's full refinement still recomputes everything
  at the end of a session. Rejected: full final reconstructions in every
  offline add-on, and one terminal full reconstruction per offline chain.
- **The half-set fallback.** The question on its 1% threshold is gone: the
  fallback is removed (section 3.2).

## 9. Code map

| Area | Files |
| --- | --- |
| Read path | `src/main/strategies/search/simple_matcher_ptcl_io.f90` (`discrete_read_imgbatch`, `stk_dims`, read counters, the read session), `src/fileio/simple_discrete_stack_io.f90`, `src/fileio/simple_imgfile.f90`, `src/fileio/simple_imghead.f90`, `src/main/image/simple_image_io.f90` |
| Particle buffer | `simple_matcher_ptcl_io.f90` (the slab and its index map, filled in one sweep), `simple_strategy3D_matcher.f90` (filled before the search loop), `simple_matcher_3Drec.f90` and `simple_rec3D_pcg_strategy.f90` (copy from it) |
| POSIX path | `src/fileio/simple_posix.c`, `src/fileio/simple_syslib.f90` (`io_read_nstreams`, the new helpers) |
| Prefetch and benchmarks | `simple_strategy2D_matcher.f90`, `simple_strategy3D_matcher.f90`, `simple_commanders_prob.f90`, `simple_matcher_3Drec.f90`, `simple_refine2D_strategy.f90`, `simple_refine3D_strategy.f90` |
| Half-sets | the project-level partition transaction (with `sp_project`), `src/main/ori/simple_oris_setters.f90` (`partition_eo`, the per-index routine), and its call sites as consumers |
| Add-on passes | `simple_commanders_solve3D.f90` (`addon_starting_state`), `simple_solve3D_utils.f90`, `simple_frozen_accum.f90`, `src/main/project/simple_project_superset.f90` (retirement), `src/main/simple_final_rec.f90`, `exec_bootstrap_rec3D` (`simple_commanders_refine3D.f90`), `simple_commanders_solve3D_chunks.f90`, `src/main/sigma2/` |
| Policies | `doc/policies/3D/solve3D_addon_policy.md` (sections 3, 5, 7, 8; frozen membership under retirement; the PCG frozen add removed), `doc/policies/3D/separate_alignment_and_reconstruction_for_multistate_peak_mem_reduction.md` (the raw-image buffer across the assignment barrier), the documents that describe half-set assignment, a short policy for the particle read path |
| Tests | `src/fileio/simple_stack_io_tester.f90`, the `oris` setter testers, `src/main/volume/simple_frozen_accum_tester.f90` |

## 10. Execution

This plan is carried out as an unattended run on the Oracle Linux test
machine (arun run `particle_io`), in eight phases. Nothing is committed; the
maintainer reviews the complete diff when all phases are done. Every phase
records itself in section 11.

Common to every phase:

- **Data.** Beta-galactosidase (5,513 particles, box 256 at 1.275 Å per
  pixel, `mskdiam=180 pgrp=d2`, 50 selected class averages), copied into a
  fresh directory under the run's scratch directory for each run; the
  originals are never written. The set is one imported stack without
  micrograph rows, so it takes the per-index half-set route and exercises
  the read path only lightly. Few-particle stacks, micrograph groups and
  shared physical files are covered by the unit-test fixtures of section 6.
- **Execution settings.** `nparts=4 nthr=6` for distributed runs and
  `nthr=24` for shared-memory runs, identical in every phase. Timed runs
  use an out-of-tree Release build with tests built; nothing else heavy
  runs while a run is timed.
- **No timing claims.** No gain is claimed from this machine (Hans,
  2026-10-14): its test set is one stack, which these changes barely touch.
  The I/O gains are measured later by the maintainer on the ribosome project
  (section 4). Timings here guard against regressions only:
  - every read-path step is implemented, in the order of section 7;
  - a step is kept unless the case that exercises it is slower than the
    build without it beyond the run-to-run spread, measured at least twice
    each way when the difference is within the spread;
  - a dropped step is recorded in section 11 with its numbers.
- **Validation.** Every `solve3D` run uses `nsample=2500` and starts from the
  one common `solve3D_cavgs` result with `cavg_ini_ext=yes`. Beta-galactosidase
  is large enough for the Fourier shell correlation (FSC). The simulated
  workflow tests are judged by their ground-truth orientation and shift
  errors.
- **Acceptance of a repeated case.**
  - `solve3D` and `refine3D` (cases A, B, C): the final FSC 0.143
    resolution of every state lies within the larger of the Phase 0 spread
    and 3 % of its Phase 0 value.
  - `solve2D` (case E), as in the global-sigma run: the last iteration's
    mean score is no lower than 0.97 times the lowest Phase 0 run, and the
    particle-weighted mean class resolution of the last iteration is no
    worse than 1.03 times the worst Phase 0 run.
  - The add-on chain (case F): every step completes, and the FSC 0.143
    resolution of every step's union map lies within the larger of the
    Phase 0 spread and 5 % of its Phase 0 value. The wider margin is for
    the accepted old weighting of frozen contributions.
  - High-level tests: each by its own checks. `simulated_workflow_6vxx` is
    recorded, not judged (flaky on unchanged code).
- **Gates.** The fast gate (run by `./compile_debug.sh`), the CTest unit
  suites of the touched areas, and `python3 scripts/check_test_registry.py .`
  whenever tests or their registration change.
- **Not in this run.** The generated indexes under `doc/code_overview/`
  (the maintainer regenerates them); stacks kept open from the search into
  the reconstruction (they wait for the access-trace replay); the production
  measurements of section 4; skill files (proposed text goes into the
  report).

### 10.1 Phase 0: benchmark split and baseline

1. Find how SIMPLE is built on this machine and build the tree as it is
   (`./compile_debug.sh`, which also runs the fast gate), plus the Release
   build for timed runs.
2. Add the benchmark split of section 5, the only source change of this
   phase, in SIMPLE's benchmark conventions (`L_BENCH_GLOB`, `tic`/`toc`):
   - particle reading timed separately from particle preparation, in the
     refine2D and refine3D matcher files, in the reconstruction (the read
     timer of `calc_3Drec` moves from `DEBUG` to the benchmark flag), in the
     PCG worker and in class-average restoration;
   - a per-partition benchmark file for the probability-table parts;
   - per pass, in the context sections: stacks opened, read requests, bytes
     requested and bytes read (prefetch lines join in Phase 2);
   - the process's peak resident memory and the node's memory in use at the
     same point.
3. Run the baseline cases on the Release build with the split:
   - **S.** `solve3D_cavgs`, single state; its final project and maps are
     kept as the common start.
   - **A.** `solve3D nsample=2500 cavg_ini_ext=yes` from S, distributed,
     twice (the spread).
   - **B.** `refine3D` from A's final project with `objfun=euclid` and
     `update_frac=0.2`, five iterations, distributed, twice.
   - **C.** B's command line in shared memory, once.
   - **E.** `solve2D` on the particles, distributed, twice.
   - **F.** `emulate_solve3D_stream nptcls_base=2000 nptcls_addon=1100
     rollback=no addon_diag=yes` (a base run and three add-ons; 1,200 gave
     two, see section 11), with A's `solve3D` settings where the program
     takes them, gridding, twice.
   - **P.** A's command line with `rec_backend=pcg`, once; its final
     project is kept for the refusal check of Phase 5.
   For every case record the exact command line, the wall time, the final
   FSC 0.5 and 0.143 resolutions per state (for E the score and class
   resolution of the acceptance), the peak memory, and the benchmark
   timings with the read time and counts. For F also each step's wall
   time, seconds per 1,000 added particles, union FSC, and the frozen
   particle images each add-on reads.
4. Run `solve3D_addon`, `simulated_workflow_1jxy`, `simulated_workflow_6vxx`
   and `cont_refine3D_1jxy`, and record their results.

Exit criteria:

- the build and the fast gate pass;
- the read time and the counts appear in every benchmark file listed in
  step 2;
- every case is recorded with its command line in section 11, with the
  spreads of A, B, E and F;
- S and P are kept for the later phases.

### 10.2 Phase 1: the read session and coalesced spans

The first two read-path steps of section 3.1: the pass-scoped read session
(one header decode per open, read-only opens without the second existence
check, reusable scratch, counters, the boundary descriptor carried to the
next batch) and the bounded coalesced spans. Every reader of
`discrete_read_imgbatch` moves to the session. The gap budget starts at 2 MB
and the scratch cap at 64 MB, as named constants; the over-read of a span is
at most the bytes it needs. The production values follow the maintainer's
measurement.

Exit criteria:

- the Release and the debug builds, the fast gate, the unit suites of the
  touched areas and the registry check pass;
- the read-path tests of section 6 that cover the session and the spans
  pass, including one open per stack per pass across a batch boundary and
  two stack rows sharing one file;
- cases A and E repeated: the acceptance holds, and each step is kept or
  dropped by the timing rule.

### 10.3 Phase 2: the POSIX path and the prefetch window

The POSIX read-only path with concurrent opens and the bounded prefetch
window of section 3.1, with its edge cases, and the prefetch counts in the
benchmark files.

Exit criteria:

- the builds, the fast gate, the unit suites and the registry check pass;
- the POSIX-path tests of section 6 pass, including descriptor exhaustion,
  leaks and the fallback modes;
- cases A and E repeated: the acceptance holds, and the step is kept or
  dropped by the timing rule.

### 10.4 Phase 3: the particle buffer

The pass-sized particle buffer of section 3.1: the immutable slab with its
index map, copied into the existing image batches; the budget (the
computing environment's per-task memory where defined, otherwise physical
memory divided by twice the number of processes the node runs); the
fallback to batch reads; the line in the memory policy.

Exit criteria:

- the builds, the fast gate, the unit suites and the registry check pass;
- the particle-buffer tests of section 6 pass;
- cases A, B and C repeated: the acceptance holds, the peak memory of each
  is recorded against Phase 0 and stays within the budget plus the Phase 0
  peak, and the step is kept or dropped by the timing rule.

### 10.5 Phase 4: project-owned half-sets

The partition transaction of section 3.2 and its consumers; the per-index
route for fields without micrograph rows, class averages and time series.

Exit criteria:

- the builds, the fast gate, the unit suites and the registry check pass;
- the half-set tests of section 6 pass;
- on beta-galactosidase the log shows the per-index route, and cases A and
  E repeated meet the acceptance;
- `solve3D_addon`, `simulated_workflow_1jxy`, `simulated_workflow_6vxx` and
  `cont_refine3D_1jxy` run, each judged as in the common rules.

### 10.6 Phase 5: PCG out of the add-on chains

The "Gridding only" ruling of section 3.3: the PCG inputs leave
`solve3D_addon` and `emulate_solve3D_stream`, `solve3D_addon` refuses a PCG
base run, and the PCG side of the frozen sets is deleted with its tests.

Exit criteria:

- the builds, the fast gate, the unit suites and the registry check pass;
- `solve3D_addon` on a copy of a PCG base run stops with the refusal before
  any work; on a gridding base it runs and passes its own checks;
- outside history documents, the tree holds no frozen PCG payload.

### 10.7 Phase 6: the add-on passes

Section 3.3 as ruled: the first add-on reads the base particles once and
writes the first bundle; later add-ons read only their own particles; the
incremental final reconstruction; the copy-on-write sigma2 transaction;
retirement with carried sums; the bundle manifest; the offline chains
publish their last add-on as it is. The add-on policy changes with it.

Exit criteria:

- the builds, the fast gate, the unit suites and the registry check pass;
- the add-on tests of section 6 pass;
- `solve3D_addon` passes its own checks;
- case F repeated twice: the acceptance holds; add-ons after the first read
  no frozen particle image (by the counts); each step's wall time and
  seconds per 1,000 added particles are reported against Phase 0.

### 10.8 Phase 7: documents and close

1. The documents of section 9: the short read-path policy, the reader
   policy in the stack-indexing policy, the half-set documents, the add-on
   and chunk policies.
2. Cases A, B, C, E and F once each, and the four high-level tests, against
   Phase 0: the acceptance holds.
3. Move this plan to `doc/implementation_notes/completed/` with `mv`, and
   write the review report beside it. The report reads on its own: what
   changed and why, in plain words; the evidence against Phase 0 in tables,
   with no timing claim beyond regression; every decision taken during the
   run and where it is recorded; which read-path steps were kept or
   dropped; what the maintainer measures next on the ribosome project
   (section 4). SIMPLE terms are explained and acronyms spelled out on first
   use; no run-internal labels.

Exit criteria:

- the build, the fast gate and the registry check pass;
- the acceptance of step 2 holds;
- steps 1 and 3 are done.

### File table

| File | Change | Phases |
| --- | --- | --- |
| `src/main/strategies/search/simple_matcher_ptcl_io.f90` | Read timers and counts; the read session; coalesced spans; the POSIX path; the particle buffer. | 0 1 2 3 |
| `src/main/strategies/search/simple_matcher_ptcl_batch.f90`, `simple_strategy3D_matcher.f90`, `simple_strategy2D_matcher.f90`, `simple_matcher_3Drec.f90` | Benchmark split; sessions per pass; prefetch; the buffer filled before the search and read by the reconstruction; consumers of the partition. | 0 1 2 3 4 |
| `src/main/strategies/search/simple_ptcl_cache.f90` | Its build reads through the session. | 1 2 |
| `src/main/commanders/simple/simple_commanders_prob.f90` | The probability-table benchmark file; sessions; consumer of the partition. | 0 1 2 4 |
| `src/main/strategies/parallelization/simple_refine3D_strategy.f90`, `simple_refine2D_strategy.f90` | Benchmark split; the partition preflight; frozen PCG plumbing out; carried sums. | 0 4 5 6 |
| `src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90` | Read timers; sessions; the buffer; the frozen add removed. | 0 1 2 3 5 |
| `src/main/strategies/parallelization/simple_calc_pspec_strategy.f90`, `simple_make_cavgs_strategy.f90`, `simple_rec3D_strategy.f90` | Sessions; consumers of the partition. | 1 2 4 |
| `src/main/strategies/parallelization/simple_rec3D_service.f90` | Frozen PCG path out; carried sums. | 5 6 |
| `src/main/class/simple_classaverager.f90`, `simple_classaverager_restore.f90` | Read timers; sessions. | 0 1 2 |
| `src/main/flex/fit/simple_flex_pca_basis.f90`, `simple_flex_pca_plane_cache.f90`, `simple_flex_reconstructor_latent_ops.f90` | Readers move to the session. | 1 2 |
| `src/main/flex/fit/simple_flex_probe_fit_estep.f90`, `simple_flex_pca_embed.f90` | Their batch loops call the FLEX batch loader, so they own the pass's read session. | 1 2 |
| `src/main/volume/simple_reconstructor_openmpoffload.f90`, `src/main/commanders/single/single_commanders_experimental.f90` | Readers move to the session. | 1 2 |
| `src/fileio/simple_discrete_stack_io.f90`, `simple_imgfile.f90`, `simple_imghead.f90`, `simple_stack_io_tester.f90` | One header decode per open; read-only opens; span reads; the POSIX path; tests. | 1 2 |
| `src/fileio/simple_stack_io_tester.f90` | Phase 3: the tests of the image slab behind the particle buffer. | 3 |
| `src/fileio/simple_image_slab.f90` | New: the immutable slab of decoded images behind the particle buffer, with its budget check and failed-allocation fallback, kept free of the builder so the stack tester can drive it. | 3 |
| `src/fileio/simple_stack_read_session.f90` | New: the pass-scoped read session over many stacks (one open per stack per pass, the descriptor carried across batch boundaries, coalesced spans, reusable scratch, exact counts), kept free of the builder so the stack tester can drive it. | 1 2 |
| `src/main/image/simple_image_io.f90` | Span reads into scratch. | 1 2 |
| `src/fileio/simple_posix.c`, `src/fileio/simple_syslib.f90` | The node's memory in use for the benchmark files (Phase 0); the POSIX read helpers; the memory budget. | 0 2 3 |
| `src/fileio/simple_discrete_stack_io.f90`, `simple_imgfile.f90` | Phase 0: the bytes of one image of the open stack, so the benchmark files count the bytes read. | 0 |
| `src/main/commanders/simple/simple_commanders_rec_distr.f90` | Read timers in the gridding workers; frozen PCG add removed; carried sums. | 0 5 6 |
| `src/main/commanders/simple/simple_commanders_solve2D.f90`, `src/main/flex/run/simple_flex_pca_project_gateway.f90` | Consumers of the partition. | 4 |
| `src/main/ori/simple_oris_setters.f90`, `simple_oris.f90` | The per-index routine kept for its explicit routes. | 4 |
| `src/main/project/simple_eo_partition.f90` | New: the project-owned partition transaction. | 4 |
| `src/main/project/simple_sp_project.f90`, `simple_sp_project_ptcl.f90` | Entry of the partition transaction; micrograph mapping. | 4 |
| `src/main/commanders/test/simple_commanders_test_highlevel_frac.f90` | Follows the session and the partition. | 1 4 |
| `src/main/commanders/test/simple_commanders_test_class.f90`, `src/main/ui/simple_test/simple_test_ui_class.f90` | New unit testers in the suite tables; frozen PCG cases out. | 1 2 3 4 5 6 |
| `src/main/ui/simple/simple_ui_solve3D.f90` | PCG inputs out of the add-on and stream-emulation programs. | 5 |
| `src/main/commanders/simple/simple_commanders_solve3D.f90`, `simple_commanders_solve3D_chunks.f90` | Refusal of a PCG base; frozen PCG plumbing out; first carried sums; incremental final; publication. | 5 6 |
| `src/main/solve/simple_solve3D_utils.f90`, `simple_solve3D_manifest.f90`, `simple_solve3D_controller.f90` | Frozen PCG plumbing out; bundle manifest; carried sums. | 5 6 |
| `src/main/volume/simple_frozen_accum.f90` | PCG payloads deleted; the carried bundle. | 5 6 |
| `src/defs/simple_refine3D_fnames.f90` | Benchmark file names; frozen PCG names out; bundle names. | 0 5 6 |
| `src/main/params/simple_parameters.f90`, `simple_parameters_parse.f90` | Keys of the frozen PCG path out, if any. | 5 6 |
| `src/main/simple_final_rec.f90`, `src/main/commanders/simple/simple_commanders_refine3D.f90` | Incremental final reconstruction; the bootstrap over the cohort. | 6 |
| `src/main/project/simple_project_superset.f90` | Retirement with carried sums. | 6 |
| `src/main/sigma2/simple_sigma2_state.f90`, `simple_euclid_sigma2.f90`, `simple_sigma2_bootstrap.f90` | The copy-on-write sigma2 transaction of the growing union. | 6 |
| `src/fileio/simple_sigma2_state_file.f90`, `simple_sigma2_files.f90` | The candidate extended from the frozen state. | 6 |
| `src/main/commanders/test/simple_commanders_test_highlevel_addon.f90`, `src/main/ui/simple_test/simple_test_ui_highlevel.f90` | Chain, rollback and bundle cases. | 6 |
| `doc/policies/3D/separate_alignment_and_reconstruction_for_multistate_peak_mem_reduction.md` | The buffer across the assignment barrier. | 3 |
| `doc/policies/3D/solve3D_addon_policy.md`, `solve3D_chunks_policy.md`, `doc/policies/stream/stream_3D_ingestion_policy.md` | Gridding only; retirement; carried sums; incremental final; publication. | 5 6 7 |
| `doc/policies/3D/refine3D_policy.md`, `reconstruct3D_pcg_policy.md`, `doc/policies/2D/solve2D_policy.md` | Half-set assignment; frozen PCG add removed. | 5 7 |
| `doc/policies/project_stack_indexing_policy.md`, `doc/policies/particle_read_path_policy.md` | The reader policy; new: the short policy of the read path. | 7 |
| `doc/implementation_notes/completed/particle_io_few_particle_stacks.md`, `particle_io_few_particle_stacks_report.md` | The plan moved on completion, and the report. | 7 |

## 11. Progress

Each phase adds its entry here: date, files changed, evidence with log
paths, and how each exit item was met.

### Phase 0: benchmark split and baseline (2026-10-13 to 2026-10-14)

Paths below are relative to the run directory `/home/elmlundho/agent_runs/particle_io`
(`scratch/` holds the scripts, logs and kept projects).

**How SIMPLE is built on the test machine.** GCC and gfortran 15.2 from
`gcc-toolset-15` (enabled in `~/.bashrc` before its interactive check, so
background jobs get it), CMake 3.26, FFTW from the system, no environment
modules. The maintainer's own tree (`~/src/SIMPLE/build`) is a Release build
with tests, installed into its build directory; the login shell points
`SIMPLE_PATH` and `PATH` at it, so every run here sets both to the run's own
build. The debug build is `./compile_debug.sh` in the repository (it also runs
the fast gate, the build-time test gate of the `fast` CTest label). The timed
runs use an out-of-tree Release build with tests, `scratch/build_rel`, made by
`scratch/build_release.sh` (CMake with `-DCMAKE_BUILD_TYPE=Release
-DBUILD_TESTS=ON`, installed into the build directory, fast gate before the
install).

**The source change: the benchmark split.** SIMPLE writes benchmark files
(`*_BENCH*.txt`, a CONTEXT section and a TIMINGS section) whenever its
benchmark flag `L_BENCH_GLOB` is on, which it is in every build. This phase
adds the particle reads to them:

- The particle reader (`simple_matcher_ptcl_io`) counts, for the whole
  process, the header reads that find a stack's dimensions, the stacks it
  opens, its read requests, the bytes it was asked for and the bytes it read
  (equal until a later phase reads spans), and the wall time spent reading;
  the reads of the downscaled 2D particle cache are timed, not counted (the
  cache has its own reader). A caller takes a snapshot before and after a
  pass, so every count in a file belongs to one pass.
- `REFINE3D_BENCH_ITERnnn_PARTppp.txt` (the 3D search worker of one
  partition) reports the search pass and the reconstruction pass apart;
  "particle reading" is now separate from "particle preparation", and the
  reconstruction's reading has its own line.
- `REFINE2D_BENCH_ITERnnn.txt` (the 2D search, written by partition 1 only,
  as before) reports its pass the same way.
- New `REC3D_BENCH_ITERnnn_PARTppp.txt`: every reconstruction pass of a
  process appends one section, for gridding (`calc_3Drec`, whose debug-only
  read timer moved to the benchmark flag) and for the PCG worker (the
  reconstruction by preconditioned conjugate gradients): reading,
  preparation and accumulation timed apart, with the pass's counts. The
  reconstructions of the 3D search workers appear here too.
- New `PROB_TAB_BENCH_…`, `PROB_TAB_NEIGH_BENCH_…` and
  `PROB_TAB2D_BENCH_ITERnnn_PARTppp.txt`: the probability-table parts, which
  make the first read of every probabilistic iteration and had no file.
- New `CAVGSUMS_BENCH_PARTppp.txt`: the class sums of a `make_cavgs`
  partition (class-average restoration from the particles).
- Every new file, the two search files and the two strategy files of the
  masters (`REFINE3D_STRATEGY_BENCH`, `REFINE2D_STRATEGY_BENCH`) report the
  process's peak resident memory and, at the same moment, the memory in use on
  the node (MemTotal minus MemAvailable).

Files changed: `src/main/strategies/search/simple_matcher_ptcl_io.f90`,
`simple_matcher_ptcl_batch.f90`, `simple_strategy3D_matcher.f90`,
`simple_strategy2D_matcher.f90`, `simple_matcher_3Drec.f90`;
`src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90`,
`simple_refine3D_strategy.f90`, `simple_refine2D_strategy.f90`;
`src/main/commanders/simple/simple_commanders_prob.f90`;
`src/main/class/simple_classaverager_restore.f90`;
`src/defs/simple_refine3D_fnames.f90` (the two new file names);
`src/fileio/simple_posix.c` and `simple_syslib.f90` (the node's memory in
use); `src/fileio/simple_discrete_stack_io.f90` and `simple_imgfile.f90`
(the bytes of one image of an open stack); this plan (file table and this
entry). For this phase the file table added Phase 0 to the row of
`simple_posix.c` and `simple_syslib.f90`, and gained the row of
`simple_discrete_stack_io.f90` and `simple_imgfile.f90`.

**Data.** `/data2/bgal_addon/control/bgal_addon.simple` was confirmed to be
one imported stack (`/data/bgal/DAT/sumstack.mrc`, 5,513 particles in
32-bit floats, box 256) with no micrograph rows, every particle selected,
as section 10 assumes.

**Decisions taken.**

- **Case F runs three add-ons with `nptcls_addon=1100`, not 1200.** The
  emulation cuts `max(1, remainder / nptcls_addon)` add-ons, so 1,200
  particles per add-on on the 3,513 left after a base of 2,000 gives two
  add-ons (1,200 and 2,313), while the plan asks for a base run and three
  add-ons. 1,100 gives three (1,100, 1,100 and 1,313), so two add-ons follow
  the first, which is what Phase 6 must show reads no frozen particle.
- **Case F starts from the control project.** `emulate_solve3D_stream`
  refuses `cavg_ini_ext` and `vol1` (they are not inputs of the stream's base
  run), so it cannot start from S; it takes A's other `solve3D` settings.
- **Case C keeps `mkdir` at its default.** The earlier defect that nested
  shared-memory `refine3D` directories is fixed on this base; C is B's
  command line with `nthr=24` instead of `nparts=4 nthr=6`.
- **S ran on the first Release build of the split;** a review then changed
  only the benchmark output (zero guards on two "% accounted" lines, and
  comments), and every other case ran on the rebuilt tree.

**Command lines.** Every case runs through `scratch/run_case.sh <name>
<kind>` (in sequence via `scratch/queue.sh`), in a fresh
`scratch/runs/<name>`, on a copy named `bgal.simple`, with `SIMPLE_PATH` set
to the Release build and its `bin` and `scripts` first on `PATH`. A kept
project refers to its maps by absolute path inside the directory that wrote
it, so the script first recreates those directories from `scratch/keep/`.
Later phases repeat these lines exactly:

| Case | Starts from | Command line |
| --- | --- | --- |
| S | the control project | `simple_exec prg=solve3D_cavgs projfile=bgal.simple mskdiam=180 pgrp=d2 nparts=4 nthr=6` |
| A1, A2 | `keep/S/1_solve3D_cavgs/bgal.simple` | `simple_exec prg=solve3D projfile=bgal.simple mskdiam=180 pgrp=d2 nsample=2500 cavg_ini_ext=yes nparts=4 nthr=6` |
| B1, B2 | `keep/A/1_solve3D/bgal.simple` | `simple_exec prg=refine3D projfile=bgal.simple vol1=<scratch>/runs/A1/1_solve3D/recvol_state01.mrc mskdiam=180 pgrp=d2 objfun=euclid update_frac=0.2 maxits=5 nparts=4 nthr=6` |
| C | as B | `simple_exec prg=refine3D projfile=bgal.simple vol1=<same> mskdiam=180 pgrp=d2 objfun=euclid update_frac=0.2 maxits=5 nthr=24` |
| E1, E2 | the control project | `simple_exec prg=solve2D projfile=bgal.simple mskdiam=180 ncls=100 nparts=4 nthr=6` |
| F1, F2 | the control project | `simple_exec prg=emulate_solve3D_stream projfile=bgal.simple nptcls_base=2000 nptcls_addon=1100 rollback=no addon_diag=yes mskdiam=180 pgrp=d2 nsample=2500 rec_backend=gridding nparts=4 nthr=6` |
| P | `keep/S/1_solve3D_cavgs/bgal.simple` | A's command line with `rec_backend=pcg` |

The high-level tests run as CTest entries in the Release build,
`ctest -R '^<test>$' --output-on-failure -V` (`scratch/hl_tests.sh`).

**Results.** All runs exited with status 0, one at a time on the otherwise
idle 24-core machine (the load average at each start is the tail of the case
before it; nothing else ran). Wall time and peak memory are from
`/usr/bin/time -v` (peak memory: the largest resident set of any one process
of the run). FSC is the Fourier shell correlation of the half-maps at 0.5
and 0.143: `RESOLUTION_FINAL_STATE01` for S, A and P, and
`RESOLUTION_STATE01_ITER005` (the last iteration) for B and C, whose
iterations each reconstruct from a 20 % sample. For E: the last mean score
and the last particle-weighted mean class resolution of the log.

| Case | Wall time | Peak memory | FSC 0.5 / 0.143 (Å), or score and class resolution |
| --- | --- | --- | --- |
| S | 10:58 | 3.07 GB | 7.59 / 6.80 |
| A1 | 14:05 | 3.19 GB | 6.66 / 4.35 |
| A2 | 14:09 | 3.25 GB | 6.80 / 4.35 |
| B1 | 2:38 | 3.03 GB | 9.60 / 6.94 |
| B2 | 2:38 | 3.03 GB | 9.60 / 7.10 |
| C | 1:46 | 3.63 GB | 10.20 / 7.77 |
| E1 | 3:06 | 0.41 GB | score 0.389, 10.97 Å |
| E2 | 3:08 | 0.41 GB | score 0.389, 10.94 Å |
| F1 | 39:22 | 3.28 GB | see below |
| F2 | 43:00 | 3.30 GB | see below |
| P | 22:56 | 6.67 GB | 6.80 / 4.41 |

Spreads between the two runs: A 0.14 / 0.00 Å; B 0.00 / 0.16 Å; E score
0.000, class resolution 0.035 Å.

The add-on chain F, per step: wall time, seconds per 1,000 added particles,
the FSC 0.143 (and 0.5) of the union map, and the frozen particle images the
add-on reads. The frozen images are counted from the new files
(`scratch/addon_reads.py`): two frozen accumulations, at the stage box 162
and the native box 256, then three passes of the final reconstruction over
the whole union (bootstrap map, residual pass, shipped map), each reading
every frozen particle once; the fourth final pass, the image power spectra
(`calc_pspec`), writes no benchmark file and is counted from its log line.

| Step | Frozen / added / union | F1: s, s per 1,000, FSC 0.143 (0.5) | F2: s, s per 1,000, FSC 0.143 (0.5) | Frozen images read |
| --- | --- | --- | --- | --- |
| base | 0 / 2,000 / 2,000 | 928, –, 4.03 (6.40) | 913, –, 4.03 (6.40) | – |
| add-on 1 | 2,000 / 1,100 / 3,100 | 396, 360, 3.89 (4.41) | 354, 322, 3.89 (4.41) | 12,000 (6 × 2,000) |
| add-on 2 | 3,100 / 1,100 / 4,200 | 636, 578, 3.80 (4.35) | 629, 572, 3.84 (4.35) | 18,600 (6 × 3,100) |
| add-on 3 | 4,200 / 1,313 / 5,513 | 402, 306, 3.75 (4.29) | 683, 520, 3.80 (4.29) | 25,200 (6 × 4,200) |

F spreads (FSC 0.143): 0.00, 0.00, 0.04 and 0.05 Å. Add-on 3's wall time
differs because the two runs searched a different number of iterations.

**Particle reads in the benchmark files** (summed over iterations,
partitions and passes; read time is summed process time). On this one-stack
set every batch call opens the stack once, and the stack sits in the page
cache, so reading is about 1 % of the work; these figures are the baseline
the later phases' counts are compared with, not a measure of the I/O cost.

| Case | Pass | Read time (s) | Stacks opened | Read requests |
| --- | --- | --- | --- | --- |
| A1 | probability tables | 6.95 | 435 | 102,780 |
| A1 | search | 6.60 | 455 | 108,293 |
| A1 | search workers' reconstruction | 6.15 | 328 | 102,780 |
| A1 | all reconstructions (`REC3D`) | 7.35 | 368 | 116,313 |
| B1 | search / reconstruction | 0.35 / 0.35 | 29 / 40 | 5,515 / 5,515 |
| C | search / reconstruction | 0.27 / 0.27 | 5 / 20 | 5,515 / 5,515 |
| E1 | 2D probability tables (all parts) | 5.98 | 340 | 93,721 |
| E1 | 2D search (partition 1 only) | 2.31 | 135 | 37,233 |
| E1 | `make_cavgs` class sums | 0.32 | 9 | 5,513 |
| F1 | all reconstructions (`REC3D`) | 19.60 | 1,316 | 289,221 |
| P | all reconstructions (`REC3D`, PCG from stage 3) | 8.91 | 352 | 109,874 |

A2, B2 and E2 match their first runs' counts exactly.

**Evidence.** Case logs `scratch/logs/case_<case>.log` (command, start
load, `/usr/bin/time -v`), summaries `scratch/logs/summary_<case>.txt`, the
queue `scratch/logs/queue_phase0.log`, read tables
`scratch/logs/read_table_phase0_*.txt`, add-on reads
`scratch/logs/addon_reads_F*.txt`; every benchmark file of every case, and
F's `emulate_solve3D_stream_report.txt`, in
`scratch/keep/bench/phase0/<case>/`. Kept for later phases:
`scratch/keep/S/1_solve3D_cavgs` (the common start), `scratch/keep/A/1_solve3D`
(input of B and C, with `final_map_name`), `scratch/keep/P/1_solve3D` (the
PCG base for Phase 5's refusal check). Builds: `scratch/logs/build_debug_base.log`
(the base, fast gate 15 of 15), `scratch/logs/build_release_split2.log`.

**High-level tests** (Release build, run one at a time after the timed
cases; logs `scratch/logs/hl_phase0_<test>.log`):

| Test | Result | Key measures |
| --- | --- | --- |
| `solve3D_addon` | passed, 9:57 | union map against the truth: correlation 0.942 (floor 0.90), FSC 0.143 4.93 Å (floor 6.03 Å); new cohort's pose-error excess over the frozen particles' 2.34° (limit 5°); union against base map 0.9989 |
| `simulated_workflow_1jxy` | passed, 4:21 | whole-volume correlation with the truth 0.898, FSC 0.143 4.54 Å |
| `simulated_workflow_6vxx` | passed, 13:57 (recorded, not judged) | correlation 0.888, FSC 0.143 4.71 Å |
| `cont_refine3D_1jxy` | passed, 8:50 | median rotation error of stages a / b / c 0.91 / 0.23 / 0.91° (limit 1.5°) |

**Exit criteria.**

- *The build and the fast gate pass.* `./compile_debug.sh` on the final tree:
  no compiler warning (the six linker notes about executable stacks are the
  base's), fast gate 15 of 15 (`scratch/logs/build_debug_phase0_final.log`);
  the Release build likewise (`scratch/logs/build_release_split2.log`);
  `scripts/check_test_registry.py` passes; no file mode changed.
- *The read time and the counts appear in every benchmark file listed in
  step 2.* The 3D and 2D search files, the reconstruction file for gridding
  and for PCG (P), the class sums of `make_cavgs` (run inside `solve2D`, E)
  and the three probability-table files, all with the peak memory and the
  node's memory; see the files under `scratch/keep/bench/phase0/`.
- *Every case is recorded with its command line, with the spreads of A, B, E
  and F.* The tables above.
- *S and P are kept for the later phases.* `scratch/keep/S`, `scratch/keep/P`,
  and `scratch/keep/A` for B and C; a Release build of this phase's code is
  kept in `scratch/keep/build_phase0` for Phase 1's before-runs.

### Phase 1: the read session and coalesced spans (2026-10-14)

**The read path as found** (confirmed against the code before editing; the
plan's description holds, with these line numbers). `discrete_read_imgbatch`
(`simple_matcher_ptcl_io.f90`) took the sorted path whenever the batch's
particle indices and stack rows ascend, and otherwise grouped the batch by
stack name. Either way, per call: the dimensions of every touched stack came
from a separate header read on first touch per process (`stk_dims` memo,
`find_ldim_nptcls`: an existence check, an open, a 1 KB header read and a
close), then `dstack_io%open` opened the stack again and read its header
again (`imgfile%open` adds a second existence check), then one positioned
read per particle, then a close. The readers were built and destroyed in
every call, so a stack split across two batches was opened twice. Opens and
closes ran in the master thread; the reads ran in parallel, one thread per
open stack, at most `io_read_nstreams` at once (which also stats the first
stack's file on every call).

**Inventory of the particle readers.** Callers of `discrete_read_imgbatch`,
all of which move to the session:

| Caller | What the pass is |
| --- | --- |
| `simple_matcher_ptcl_batch.f90`, `build_batch_particles3D` | one batch of the 3D search (`simple_strategy3D_matcher.f90`) and of the 3D probability tables (`prob_tab`, `prob_tab_neigh` in `simple_commanders_prob.f90`) |
| `simple_matcher_ptcl_batch.f90`, `build_batch_particles2D` | one batch of the 2D search (`simple_strategy2D_matcher.f90`) and of `prob_tab2D`, when the downscaled cache is not in use |
| `simple_matcher_3Drec.f90`, `calc_3Drec` | gridding reconstruction, batches of every state and half |
| `simple_rec3D_pcg_strategy.f90`, `execute_rec3D_pcg_worker` | PCG reconstruction worker, batches of every state and half |
| `simple_classaverager_restore.f90`, `cavger_assemble_sums` and `transform_ptcls` | offline class sums (`make_cavgs`), and the particles of one class |
| `simple_calc_pspec_strategy.f90` | the particle power spectra that seed the noise model |
| `simple_ptcl_cache.f90` | building the downscaled 2D cache |
| `simple_flex_pca_basis.f90`, `simple_flex_pca_plane_cache.f90`, `simple_flex_reconstructor_latent_ops.f90` | the FLEX principal-component passes |
| `simple_reconstructor_openmpoffload.f90` | the GPU-offload reconstruction (not built here) |
| `single_commanders_experimental.f90` (two places) | nanoparticle particle processing and trajectory averages |
| `simple_commanders_test_highlevel_frac.f90` | the fractional-update high-level test's own reconstruction |

Other places that open particle stacks, which keep their readers:

- the downscaled 2D cache's reader (`ptcl_cache_read_batch`), one cache file
  read through `dstack_io`, as section 3.1 rules;
- `read_imgbatch` (`stack_io`, in the same module), used by `ctf_correct`
  (`simple_commanders_imgproc.f90`) and the nanoparticle commanders: whole
  ranges in stack order, not sampled passes;
- `dstack_io` read without a session by the project's stack splitting and
  sub-stack writers (`split_stk`, `write_substk` in
  `simple_sp_project_stk.f90`): export paths, one pass over the stacks in order;
- `stack_io` readers in image-processing commanders, extraction and export
  (`simple_procimgstk.f90`, `exec_write_classes`, the STAR exporters), which
  read whole stacks or write new ones.

**What changed.** A particle pass now reads through a *read session* that
lives for the whole pass (one search, one reconstruction, one probability
table, one class-sum pass) instead of building and destroying its readers in
every batch call.

- The session itself is builder-free and lives in a new file,
  `src/fileio/simple_stack_read_session.f90`, so the stack tester can drive
  it directly. For each batch it groups the requested images by stack file
  (stack rows sharing one file are one group), opens each group's stack once
  in the calling thread, reads the groups in parallel (one thread per open
  stack, at most `io_read_nstreams` at once, now probed once per session
  instead of once per call), and closes them. The stack of the batch's last
  group stays open for the next batch, so a stack split across a batch
  boundary is opened once. `simple_matcher_ptcl_io.f90` wraps it as
  `ptcl_read_session` (new, read, kill), which maps particle indices to stack
  and in-stack index through the project and feeds the benchmark counters of
  Phase 0. `discrete_read_imgbatch`, its per-process dimension memo and the
  separate header read on first touch are gone.
- **One header decode per open.** `dstack_io%open`
  (`simple_discrete_stack_io.f90`) takes the dimensions and image count from
  the header its open reads, and still stops on a box mismatch. Read-only opens
  in `imgfile%open` (`simple_imgfile.f90`) open with status "old" and skip the
  second existence check (a missing file fails the open itself; TIFF, which
  opens no unit, keeps its check). Per touched stack and pass, the file
  operations go from existence check, open, header read and close, plus
  existence check, open, existence check and header read, to one open with one
  header read.
- **Coalesced spans.** Within a group the needed images, sorted by their index
  in the stack, are cut into spans by `plan_spans`: a gap is bridged when it is
  at most `SPAN_GAP_BUDGET` (2 MB of file), the span stays within its share of
  `SPAN_SCRATCH_CAP` (64 MB of decoded scratch for the whole session, split
  over the read streams, the 16-bit buffer of a mode-12 stack counted), and its
  over-read stays at most the bytes it needs. Each span is one positioned read
  through the existing decoder (`rSlices`), into a per-stream scratch that is
  reused for the whole pass together with its mode-12 conversion buffer; the
  needed images are copied out, and a repeated index is read once. Pixel modes
  0, 1, 2, 6 and 12 read in spans; any other mode (Fourier modes 3 and 4)
  keeps one read per image.
- **Every caller moved.** All fifteen call sites of the inventory own a
  session spanning their batch loop. The FLEX batch loader takes the caller's
  session, which brought `simple_flex_probe_fit_estep.f90` and
  `simple_flex_pca_embed.f90` into the file table for this phase. The 2D cache
  reader keeps its own reader and only lost its now-redundant size hint.
- **Rule made explicit.** While a session lives it may hold a stack open, and
  a file connects to one Fortran unit, so nothing else in the process may open
  that stack in the meantime. The only place that did, the nanoparticle
  `ptclsproc_nano`, read every batch twice (once by particle range, then
  correctly by particle index); the redundant first read is removed.
- Removed with the old reader: the verbose per-batch log of `cls_expansion`
  and the "header probes" benchmark line (there are no separate header reads
  left to count).

Files changed in this phase: `src/fileio/simple_stack_read_session.f90`
(new), `simple_discrete_stack_io.f90`, `simple_imgfile.f90`,
`simple_stack_io_tester.f90`; `src/main/strategies/search/simple_matcher_ptcl_io.f90`,
`simple_matcher_ptcl_batch.f90`, `simple_strategy3D_matcher.f90`,
`simple_strategy2D_matcher.f90`, `simple_matcher_3Drec.f90`,
`simple_ptcl_cache.f90`; `src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90`,
`simple_calc_pspec_strategy.f90`; `src/main/commanders/simple/simple_commanders_prob.f90`;
`src/main/class/simple_classaverager.f90` (an unused import in an interface),
`simple_classaverager_restore.f90`; `src/main/flex/fit/simple_flex_pca_basis.f90`,
`simple_flex_pca_plane_cache.f90`, `simple_flex_reconstructor_latent_ops.f90`,
`simple_flex_probe_fit_estep.f90`, `simple_flex_pca_embed.f90`;
`src/main/volume/simple_reconstructor_openmpoffload.f90`;
`src/main/commanders/single/single_commanders_experimental.f90`;
`src/main/commanders/test/simple_commanders_test_highlevel_frac.f90`; this plan.

**Decisions taken.**

- *The session core is a new builder-free module in `src/fileio`*, not part of
  `simple_matcher_ptcl_io.f90`, so that the unit tests can drive it on
  synthetic stacks without a project; the particle-level wrapper stays where
  the plan puts the session.
- *The 64 MB scratch cap is for the whole session*, split over its read
  streams (4 MB per stream with 16 streams). A per-stream cap would allow
  1 GB of scratch with 16 streams. Because the share depends on the device's
  stream count, the span rule is a pure routine (`plan_spans`) that the tests
  check with explicit caps.
- *Images still pass through the decoder they always used*: a span is read and
  decoded by `rSlices` into scratch, then copied into the image with the same
  zero padding a single read leaves. The copy is the price of one request per
  span; Phase 2's own decoder writes into the images directly.
- *The `solve2D` acceptance is judged against all runs of the Phase 0 code.*
  The plan applies the rule "as in the global-sigma run", whose ruling judged
  against "the four `solve2D` runs on unchanged code (Phase 0's E1 and E2, and
  the two reruns of the base build)". Here the kept Phase 0 build was rerun six
  times for the timing rule anyway, and those runs show the measure's own
  spread: their class resolutions are 10.88, 11.20, 11.09, 11.76, 11.21 and
  11.04 Å. Read narrowly (Phase 0's E1 and E2 only, limit 1.03 × 10.97 =
  11.30 Å), the rule is failed by the unchanged code itself (11.76 Å). With
  all eight runs of the unchanged code the limit is 1.03 × 11.76 = 12.11 Å and
  the score floor 0.97 × 0.389 = 0.377. This reference set holds for the later
  phases. The maintainer may want to confirm the reading.

**Bit-exactness.** Besides the unit tests below, two single-threaded runs
gave byte-identical outputs on the Phase 0 and the Phase 1 builds
(`scratch/bitexact.sh`, logs `scratch/logs/bitexact_phase{0,1}_{rec,cavgs}.log`):
`reconstruct3D` of the kept project A (gridding, all 13 volumes and the FSC
file identical; one stack open instead of 12, 835 read requests instead of
5,513, 2.27 GB read for 1.45 GB needed) and `make_cavgs` of the control
project in two parts (the class averages, the per-half averages, the class
sums state and the FRCs identical; per part one open instead of three or four,
173 requests instead of about 2,757).

**Tests** (in the stack I/O suite of the `unit_core` fast suite, 217 checks,
all passing):

- the images of the session equal, bit for bit and with their storage padding,
  a single read of the same image, for stacks of 1, 5 and 200 images, in pixel
  modes 0, 1, 2, 6 and 12, for sorted, unsorted and repeated index lists, read
  in batches of four so that stacks are split across batches, in the order
  asked;
- one open per stack per pass: three stacks, one split across two batches,
  three opens; twenty stacks (more than the read streams of any device) in
  batches of four, twenty opens, every image exact;
- a file reached by two stack rows in one batch is one group, one open, one
  request for its two adjacent images;
- spans: a run with one gap is one request over five images; every other image
  is one request with an over-read (4 images) within the needed (5); distant
  groups are four requests with no over-read; repeated, unsorted indices are
  one request, each image read once;
- the span planner: at 576 KiB per image a gap of three images is bridged and
  a gap of four is not (the 2 MB budget); a span never exceeds its image cap;
  repeats stay in their span;
- a session whose box differs from the stack's stops the process (checked in a
  forked child).

**Timing rule** (Release builds, `nparts=4 nthr=6`, one case at a time on the
otherwise idle machine; the load average at each start is the tail of the case
before it; logs `scratch/logs/case_<case>.log`, queue
`scratch/logs/queue_phase1.log`, summaries `scratch/logs/p1_timing_summary.txt`
and `p1_E_spread.txt`). Three builds: the Phase 0 code (`scratch/keep/build_phase0`),
the session without coalescing (the Phase 1 source with the span planner's cap
set to one image, built for this measurement and then deleted), and the full
Phase 1 code.

| Case | Phase 0 | Session only | Session and spans |
| --- | --- | --- | --- |
| A (`solve3D`) wall time | 14:19, 14:21 | 14:17, 14:15 | 14:18, 14:13 |
| A FSC 0.5 / 0.143 (Å) | 6.94 / 4.35, 6.66 / 4.41 | 6.80 / 4.29, 6.80 / 4.41 | 6.66 / 4.41, 6.94 / 4.41 |
| E (`solve2D`) wall time | 3:07, 2:59, 2:56, 3:05, 3:09, 3:06 | 3:05, 3:12 | 3:09, 3:08, 3:14, 3:12, 3:13, 3:10 |
| E score, class resolution (Å) | 0.389 throughout; 10.88, 11.20, 11.09, 11.76, 11.21, 11.04 | 0.389; 11.50, 10.91 | 0.389, 0.389, 0.389, 0.389, 0.389, 0.388; 10.70, 11.88, 11.45, 10.97, 11.01, 11.26 |

- **A:** neither step is slower; the A acceptance holds (FSC 0.143 of 4.41 Å
  in both runs, within 3 % of Phase 0's 4.35 Å).
- **E:** the full Phase 1 build averages 3:10.8 against 3:03.8 for the Phase 0
  build. That is 7 s (3.8 %), within the Phase 0 build's own run-to-run spread
  of 12.7 s; the session-only build lies between. By the plan's rule both steps
  are kept. The acceptance holds on every Phase 1 run (scores at least 0.388,
  class resolutions at most 11.88 Å against the limit of 12.11 Å above).
- **Where the time goes on this data.** The benchmark files show where the E
  difference comes from. This set is one stack in the page cache, so requests
  and opens cost almost nothing, while every image now also passes through the
  scratch copy, and spans decode the bridged images too. Summed over a run's
  processes, the measured particle reading of E's probability tables rises
  from 6.0 to 9.4 s, and of A's probability tables from 6.8 to 13.7 s (A also
  reads 52 % more bytes than it needs, because its sample of 45 % of the
  particles leaves short gaps that spans bridge). Meanwhile the requests fall
  (A's tables 102,780 to 18,652; E's 93,721 to 6,120) and so do the opens (A's
  tables 435 to 164; E's 340 to 68). These counts are what the production
  measurements of section 4 will price.

**Exit criteria.**

- *The Release and the debug builds, the fast gate, the unit suites of the
  touched areas and the registry check pass:* `./compile_debug.sh` with no
  compiler warning (the six linker notes are the base's) and the fast gate 15
  of 15 (`scratch/logs/p1_build_debug_final.log`); the Release build and its
  fast gate (`scratch/logs/p1_build_release2.log`);
  `scripts/check_test_registry.py` passes; no file mode changed.
- *The read-path tests of section 6 that cover the session and the spans pass,
  including one open per stack per pass across a batch boundary and two stack
  rows sharing one file:* the tests above.
- *Cases A and E repeated, the acceptance holds, and each step is kept or
  dropped by the timing rule:* the table above; both steps kept.

Benchmark files of every timed run: `scratch/keep/bench/phase1/<case>/`. A
Release build of this phase's code is kept in `scratch/keep/build_phase1` for
Phase 2's before-runs.

### Phase 2: the POSIX read path and the prefetch window (2026-10-14)

**Ruling applied.** The 2026-10-14 ruling that the static review was deleted
on purpose: its link is gone from this plan's header, the sentence that the
plan was reviewed and revised stays.

**What changed.**

- **Plain read-only descriptors** (`src/fileio/simple_posix.c`, bound in
  `simple_syslib.f90`): open read-only (close-on-exec), a positioned read
  that continues short reads and retries interrupted ones (`EINTR`) and
  reports end of file by returning fewer bytes, a readahead hint
  (`posix_fadvise` `WILLNEED`; reported as unavailable, and ignored, where the
  platform has none), close, and every failure returned as an error number.
  No Fortran unit is involved, so the one-unit-per-file restriction of
  Phase 1 no longer applies to these reads.
- **The session reads through them** (`simple_stack_read_session.f90`, now
  the whole read path for MRC pixel modes 0, 1, 2, 6 and 12). Each read
  stream opens, reads and closes its own stacks inside the parallel region, so
  up to `io_read_nstreams` opens are in flight at once, and the stack at a
  batch boundary stays open for the next batch as in Phase 1. A stack's
  1024-byte header is read once per open and decoded by the same routine as
  the image file reader (`MrcImgHead%decode`, new in `simple_imghead.f90`):
  the dimensions, mode, signedness, image count and the data offset, 1024
  bytes plus the extended header, as 64-bit byte offsets. Spans are read as in
  Phase 1, into a per-stream byte buffer, and each image is decoded straight
  into its image with the conversions of the image file reader (float32
  copied, float16 through `simple_float16_to_float32` row by row, bytes and
  16-bit integers by their signedness), padding zeroed and the Fourier flag
  cleared. Missing, truncated and wrongly sized stacks, an index out of range
  and a failed read stop the run with the stack's name and the reason,
  reported after the parallel region. Any other pixel mode (the Fourier modes 3
  and 4) falls back to the image file reader.
- **The prefetch window.** After a batch the session is told the next one: by
  the pass's whole particle list when the caller reads in slices of it (the
  reconstructions, the power spectra, the FLEX passes, the cache build), or by
  the list a caller gives `expect` before its loop (the 3D and 2D searches and
  the three probability tables). It opens the next batch's stacks
  concurrently (at most `PREFETCH_MAX_STACKS`, 64, held open), plans their
  spans and hands the reads to a producer thread, which reads them into a
  window buffer it owns (`PREFETCH_WINDOW_BYTES`, 128 MiB) while the caller
  computes; spans beyond the window get a `WILLNEED` hint. The next read
  decodes every span the window covers from it and reads the others itself, so
  a wrong guess costs reads, never bytes. The window and the decoded images
  are the double buffer: the producer fills the window while the caller works
  on the images of the batch before.
- **Benchmark counts.** Every benchmark file with the particle read lines now
  also reports the bytes prefetched, the prefetched bytes images were decoded
  from, and the bytes advised.
- Removed as dead: the span reader, mode queries, stack-name getter and
  image-size getter that Phase 1 had added to `dstack_io`, the reusable
  mode-12 buffer Phase 1 had added to `rSlices`, and the bytes-per-pixel
  getter of `imgfile` (Phase 0).

Files changed in this phase: `src/fileio/simple_posix.c`, `simple_syslib.f90`,
`simple_imghead.f90`, `simple_imgfile.f90`, `simple_discrete_stack_io.f90`,
`simple_stack_read_session.f90`, `simple_stack_io_tester.f90`;
`src/main/strategies/search/simple_matcher_ptcl_io.f90`,
`simple_strategy3D_matcher.f90`, `simple_strategy2D_matcher.f90`;
`src/main/commanders/simple/simple_commanders_prob.f90`;
`src/main/strategies/parallelization/simple_calc_pspec_strategy.f90`; this
plan.

**Decisions taken.**

- *The producer reads into a buffer the C side owns.* Its lifetime then
  follows the producer, not a Fortran variable that could leave scope while a
  read is pending. The producer thread is created free of the binding of the
  OpenMP master thread, which would otherwise share one core with it.
- *Header decode is the image file reader's own*, applied to the bytes the
  descriptor read, so there is one authoritative decode. The header's own
  checks (box, image count) are the session's, as in Phase 1.
- *How the next batch is known.* The window works over the order the caller
  reads; callers that read slices of a whole list need nothing, the others
  declare their pass with `expect`. A batch off the expected order finds its
  place again in the list; the class-sum pass of `make_cavgs` and the
  nanoparticle commanders read without a window.
- *Out of descriptors, the run stops* with the error number (24, too many open
  files) and the stack's name. The session holds at most the window's stacks,
  the boundary stack and one stack per read stream.
- *Test hooks in the C layer* (at most a given number of bytes per read call,
  an `EINTR` every few calls, advice reported unavailable, and a lower limit
  on open files in a forked child) are how short reads, interrupted reads,
  missing advice and descriptor exhaustion are tested; only the tests set them.

**Bit-exactness.** The single-threaded end-to-end runs of Phase 1
(`scratch/bitexact.sh`: `reconstruct3D` of the kept project A and
`make_cavgs` of the control project in two parts) gave all 18 output files
byte-identical between the Phase 0 build and the final Phase 2 build
(`scratch/logs/bitexact_phase{0,2}_{rec,cavgs}.log`). In the reconstruction,
1.32 of the 1.42 GB read ahead were decoded from the window.

**Tests** (stack I/O suite of `unit_core`, 240 checks, all passing; every image
is compared bit for bit with a single read, storage padding and the Fourier flag
included, both images first filled with junk and flagged as transforms so that
a read must reset them):

- all Phase 1 cases, now on the descriptor path;
- prefetching: a pass read in batches, each followed by the prefetch of the
  next, is exact, and every prefetched byte is used; a wrong guess is read
  ahead, never used, and changes nothing;
- an extended header of 384 bytes; a stack whose last image lies beyond 4 GiB
  (a sparse file of 4.3 million images), read on both paths;
- a missing stack, a truncated one (end of file before the requested image), a
  stack of another box and running out of descriptors each stop the process
  with the reason (forked children; the reasons were checked in the log);
- reads cut to at most 100 bytes per call with an `EINTR` every third call,
  and advice reported unavailable, read exactly, with and without the window;
- mode 4 takes the image file reader and reads exactly as before;
- ten passes over 20 stacks (more than the read streams), opened in parallel
  each time: no descriptor outlives its session, and a session holds at most
  the window's stacks and the boundary stack.

**Timing rule** (Release builds, `nparts=4 nthr=6`, one case at a time on the
otherwise idle machine; logs `scratch/logs/case_<case>.log`, queue
`scratch/logs/queue_phase2.log`, summary `scratch/logs/p2_timing_summary.txt`;
before: the Phase 1 build kept in `scratch/keep/build_phase1`):

| Case | Phase 1 | Phase 2 |
| --- | --- | --- |
| A (`solve3D`) wall time | 14:02, 14:18, 14:16 | 14:11, 14:34, 14:34 |
| A FSC 0.5 / 0.143 (Å) | 6.80 / 4.35, 6.66 / 4.41, 6.66 / 4.41 | 6.66 / 4.35, 6.66 / 4.35, 6.66 / 4.41 |
| E (`solve2D`) wall time | 3:09, 3:10, 3:13 | 3:06, 3:07, 3:06 |
| E score, class resolution (Å) | 0.389; 10.95, 10.81, 11.05 | 0.389; 10.85, 11.06, 11.01 |

- **E** is not slower (means 3:06 against 3:11).
- **A** averages 14 s (1.6 %) longer, within the run-to-run spread (16 s for
  the Phase 1 build, 23 s for Phase 2), measured three times each way: by the
  plan's rule the step is kept. The benchmark files place the variation in
  the compute-bound probability-table fill, which varies by 75 s between runs
  of either build. The measured particle reading halves (A's probability
  tables 13.4 to 7.3 s and search 13.1 to 7.0 s summed over processes; E's
  tables 9.1 to 5.3 s), while particle preparation rises by about 2 to 3 s,
  probably the producer thread competing for memory bandwidth on this
  page-cached, one-stack set.
- **The window works.** In A, 57 % of the bytes the tables, the search and
  its reconstruction read came through the window (the search's 26.5 GB
  prefetched, all but 11 MB used); in E, 80 % (the 2D search's prefetch is
  used entirely).
- **Acceptance** holds: A's FSC 0.143 is 4.35 or 4.41 Å (within 3 % of
  Phase 0); E's scores are 0.389 and class resolutions at most 11.06 Å
  (limit 12.11 Å, Phase 1's reference set).

**Exit criteria.**

- *The builds, the fast gate, the unit suites and the registry check pass:*
  `./compile_debug.sh` with no compiler warning (the six linker notes are the
  base's) and the fast gate 15 of 15 (`scratch/logs/p2_build_debug_final.log`);
  the Release build and its fast gate (`scratch/logs/p2_build_release2.log`);
  `scripts/check_test_registry.py` passes; no file mode changed.
- *The POSIX-path tests of section 6 pass, including descriptor exhaustion,
  leaks and the fallback modes:* the tests above.
- *Cases A and E repeated, the acceptance holds, and the step is kept or
  dropped by the timing rule:* the table above; kept.

Benchmark files of every timed run: `scratch/keep/bench/phase2/<case>/`. A
Release build of this phase's code is kept in `scratch/keep/build_phase2`.

### Phase 3: the particle buffer (2026-10-14)

**What changed.** A 3D search now reads its share of the sample once, before
its first batch, into a particle buffer, and both the search batches and the
reconstruction that follows in the same process take their images from it: the
third read of an iteration (the search workers' reconstruction rereading the
particles) is gone.

- **The slab** (`src/fileio/simple_image_slab.f90`, new): one contiguous
  allocation of the share's decoded images, `box` by `box` 32-bit floats each,
  filled slot by slot, sealed, and from then on only read by copy into the
  existing bounded image batches. A copy equals the image a read delivers bit
  for bit, padding zeroed and flagged as real space; changing a copy never
  touches the slab. It is allocated only within its budget; a refused or failed
  allocation leaves it empty and never stops the run.
- **The pass buffer** (`simple_matcher_ptcl_io.f90`): `fill_pass_buffer` reads
  the share through the pass's read session (so with the prefetch window of
  Phase 2) into the slab and records each particle's slot; from then on, a batch
  whose particles are all in the buffer is copied from it without touching a
  file, and any other batch is read from the stacks as before. The search
  (`simple_strategy3D_matcher.f90`) fills it before its batch loop and releases
  it after its partial reconstructions, so it lives for one pass in one process.
  The gridding and PCG reconstructions need no change: they read through the
  same read session. Probability tables, power spectra, class sums and the
  final reconstruction keep their reads, as the plan says.
- **The budget** (`ptcl_buffer_budget`): half the
  task memory a queue system grants (`job_memory_per_task` of the project's
  computing environment, when the queue is not local); otherwise the node's
  physical memory (`MemTotal`, a new helper in `simple_posix.c`) divided by twice
  the number of processes the node runs, the parts and their master in a
  distributed run, one process in shared memory. The run first also read a
  user budget from an environment variable, `SIMPLE_PTCL_BUFFER_MB`; the
  maintainer ruled it out on 2026-10-14 (SIMPLE takes no new environment
  variables), and it was removed.
- **Benchmark lines.** The search's file reports the budget, the buffer's size,
  the outcome (filled, off, over budget or allocation failed), the fill's reads
  and time apart, and in every pass's counts the images copied from the buffer.
- **Policy.** The lifetime table of
  `doc/policies/3D/separate_alignment_and_reconstruction_for_multistate_peak_mem_reduction.md`
  gains the buffer, living across the assignment barrier, with its budget rule;
  the worker memory target becomes the larger working set plus the buffer.

Files changed in this phase: `src/fileio/simple_image_slab.f90` (new),
`simple_posix.c`, `simple_syslib.f90`, `simple_stack_io_tester.f90`;
`src/main/strategies/search/simple_matcher_ptcl_io.f90`,
`simple_strategy3D_matcher.f90`; the memory policy above; this plan (three rows
added to the file table for this phase).

**Decisions taken.**

- *The buffer is consulted by the read session*, so every reader inside the
  pass (the search batches, the gridding reconstruction and the PCG worker)
  takes it without new arguments, and a batch it does not wholly cover reads
  from the stacks.
- *Under a queue system the budget is half the granted task memory*, not all
  of it: `job_memory_per_task` is set (16 GB by default) for every project, and
  the rest of the task needs the other half, as the factor two of the node rule
  says. The local queue ignores that number, so here the node rule applies:
  6.7 GB for each of the four parts of a distributed run, 33.5 GB in shared
  memory.
- *The user's budget is an environment variable*, as the stream count's is, so
  no program's interface changes; it also turns the buffer off for the
  comparisons below.
- *Acceptance of a case with two Phase 0 values* is read as the band between
  them, widened at each end by the larger of the Phase 0 spread and 3 % of that
  end: for B, 6.94 to 7.10 Å becomes 6.73 to 7.31 Å.

**Equal outputs with and without the buffer** (`scratch/bufexact.sh`, logs
`scratch/logs/bufexact_*.log`; one `refine3D` iteration of the kept project A
in shared memory, one thread, a fixed seed, `objfun=euclid update_frac=0.2`).
The polar search is not reproducible from run to run even on unchanged code
and one thread: its transforms are planned with FFTW's `MEASURE`, which picks
algorithms by timing, and three unbuffered runs gave two different results.
The comparison therefore uses the Cartesian refinement (`refine=cont`), whose
transforms are planned deterministically: three unbuffered runs agree byte
for byte, and the run with the buffer, and the run with a 1 MiB budget (over
budget, batch reads), give the same 13 volume and Fourier shell correlation
files and the same orientations, byte for byte, on the final code too. With
the buffer, the search and the reconstruction open no stack and copy their
1,103 images from it.

**Tests** (stack I/O suite of `unit_core`, 248 checks, all passing): images
read, put in a slab and sealed come back bit for bit, padding zeroed and flagged
as real space, into junk-filled images, three times over (as a weighted
reconstruction takes a particle once per state), with every copy changed in
between; float32 pixels with a signalling NaN payload pass the read and the
slab bit for bit; a slab one byte over its budget is refused and stays empty;
a slab of 67 PB fails to allocate and reports it.

**Timing rule and memory** (Release builds, `nparts=4 nthr=6` for A and B,
`nthr=24` for C, one case at a time on the otherwise idle machine; logs
`scratch/logs/case_<case>.log`, queue `scratch/logs/queue_phase3.log`, summary
`scratch/logs/p3_timing_summary.txt`; before: the Phase 2 build kept in
`scratch/keep/build_phase2`; peak memory is the largest resident set of any
process of the run):

| Case | Phase 2: wall time | Phase 3: wall time | Phase 3: peak memory (Phase 0) | Phase 3: FSC 0.143 (Å) |
| --- | --- | --- | --- | --- |
| A (`solve3D`) | 14:21, 14:32 | 14:19, 14:23 | 3.21, 3.22 GB (3.19, 3.25 GB) | 4.29, 4.41 |
| B (`refine3D`, distributed) | 2:43, 2:41, 2:37 | 2:45, 2:43, 2:38 | 3.03 GB each (3.03 GB) | 6.94, 7.25, 7.10 |
| C (`refine3D`, shared memory) | 1:46, 1:48, 1:49 | 1:48, 1:48, 1:44 | 3.92 GB each (3.63 GB) | 7.77 each |

- The buffer is kept: no case is slower beyond the run-to-run spread (B and
  C averages within 1 s, A faster).
- Peak memory rises only where the buffer is a large part of a process: C, in
  shared memory, holds all 1,103 sampled particles (289 MB) on top of its
  Phase 0 peak (3.63 to 3.92 GB); a distributed part holds about 140 MB in A and
  80 MB in B, inside the master's larger peak. Every peak stays within the
  budget plus the Phase 0 peak.
- The reads it removes: in A, the search workers' reconstruction opened 164
  stacks for 44,715 requests in Phase 2 and opens none now; the search copies
  its batches from the buffer, and the fill reads the share once (with the
  window). The reconstruction reads that remain in A are the final
  reconstruction's, outside the search.
- Acceptance: A within 3 % of Phase 0's 4.35 Å; B within the band above (the
  Phase 2 build also gave 7.25 Å in its third run); C equal to Phase 0's 7.77 Å.

**Exit criteria.**

- *The builds, the fast gate, the unit suites and the registry check pass:*
  `./compile_debug.sh` with no compiler warning and the fast gate 15 of 15
  (`scratch/logs/p3_build_debug_final.log`); the Release build and its fast gate
  (`scratch/logs/p3_build_release2.log`); `scripts/check_test_registry.py`
  passes; no file mode changed.
- *The particle-buffer tests of section 6 pass:* bit-identical outputs with and
  without the buffer, the fallback over budget with the same outputs, no read
  from disk in the reconstruction (no stack opened), the slab unchanged and
  reused (above); the failed allocation (unit test).
- *Cases A, B and C repeated: the acceptance holds, the peak memory of each is
  recorded against Phase 0 and stays within the budget plus the Phase 0 peak,
  and the step is kept or dropped by the timing rule:* the table above; kept.

Benchmark files of every timed run: `scratch/keep/bench/phase3/<case>/`.
