# Particle I/O on per-micrograph stacks: reading, half-sets and add-on passes

Implementation note, 2026-10-12. **Plan; nothing implemented.**

**Reviewed 2026-10-13** ([static review](../completed/particle_io_few_particle_stacks_review_2026-10-13.md)),
and revised the same day to meet it. The review found the diagnosis sound but
the plan not ready: two scientific contracts and several designs were
missing. Hans ruled on the contracts the same day (section 8). The designs
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

**When.** The measurements of section 4 run once the global-sigma
refactoring has been validated (Hans, 2026-10-12). Their numbers set the
order of the read-path steps of section 3.1.

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
  filament identity where the project has it, not the physical path.
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
  solutions and speed, not PCG-level precision. `solve3D_addon` refuses
  `rec_backend=pcg`, as `solve3D_chunks` already forces gridding. The PCG
  side of the frozen sets (the `B`/`D` payloads in `simple_frozen_accum.f90`
  and the frozen add in the PCG backend) is then dead and is deleted. This
  step does not wait for the global-sigma refactoring. It also removes the
  question of how a PCG chain would seed `bootstrap_rec3D`, which forces its
  bootstrap map through gridding.

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

These run after the global-sigma refactoring has been validated, on the code
as it then stands. They need no change to the read path. The only code is
the benchmark split of section 5, which the time-share measurement needs.
Their numbers then decide:

- the order of the read-path steps, and which of them are worth doing;
- whether stacks need to be kept open from the search into the
  reconstruction;
- whether the stream-count cap on spinning disks changes.

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
  per-partition files, the refine3D strategy file and `SOLVE3D_BENCH.txt`,
  in SIMPLE's format (CONTEXT, TIMINGS, % accounted for).
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
    streams. Descriptor exhaustion, leaks and pool eviction are covered.
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

1. **The global-sigma refactoring is validated** first
   (`global_sigma2_only.md`).
2. **Measurements** (section 4): the benchmark split, then the file
   operations, the time share, the open and positioned-read latencies, the
   access-trace replay, the stream-count sweep and the whole-stack ribosome
   run.
3. **The read path.** The numbers set the order of its steps. The order shown
   is provisional:
   - the pass-scoped read session (one open per stack per pass);
   - bounded coalesced reads within a stack run;
   - the POSIX read-only path: concurrent opens, then the bounded prefetch
     window;
   - the pass-sized particle buffer in 3D, with its budget and the policy
     line;
   - stacks kept open from the search into the reconstruction, only if the
     trace replay supports it.
4. **Half-sets per stack** are ruled and do not wait for the numbers. The
   project-owned partition transaction is written first; the labels land
   with the first read-path step, and their effect on reconstruction opens is
   measured with it.
5. **Add-on passes,** as a separate track that builds on the global-sigma
   refactoring. Removing PCG from the add-on chains can land before, on its
   own.

Each step ends with a review of its numbers against the measurements of
step 2.

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
