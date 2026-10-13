# Particle I/O on per-micrograph stacks: reading, half-sets, add-on passes and data layout

Implementation note, 2026-10-12. **Plan; nothing implemented.**

It merges two notes of the same day:

- the few-particle-stacks proposal (`665ac2731`), written against master
  `32a495cef`;
- `particle_io_handling.md` from the agent checkout, now retired.

RELION references are to its master at `05cd1ad9`.

**Sampling strategies are decoupled from this note** (Hans, 2026-10-12). The
goal is a data path, and in time a layout on disk, that lets any sampling
scheme read cheaply. Workflows should not have to be designed around disk
latency. Nothing here changes which particles are drawn.

**When.** The measurements of section 4 run once the global-sigma
refactoring has been validated (Hans, 2026-10-12). Their numbers are the
basis for deciding which steps of section 3 come next, and whether the
container layout is worth its own design note.

**The consolidated-stack experiment runs first** (Hans, 2026-10-13). It
measures on the ribosome project, before any further change, what one file
instead of thousands is worth (section 4).

The plan was prompted by the stream ribosome project: thousands of
per-micrograph stacks with a handful of particles each, on which `solve2D`
and `solve3D` spend their time in file metadata. Companion documents:

- [particle_cache_policy.md](../../policies/2D/particle_cache_policy.md);
- [solve3D_policy.md](../../policies/3D/solve3D_policy.md);
- [solve3D_addon_policy.md](../../policies/3D/solve3D_addon_policy.md);
- the global-sigma refactoring outline (`global_sigma2_only.md`).

Rulings so far (Hans, 2026-10-12):

- **No particle is rewritten into another file.** No cache, raw or cropped,
  is added. The existing cropped 2D cache and the extraction's pixel format
  stay as they are.
- **Half-sets per stack** for new partitions (section 3.2).
- **Add-on passes** (section 3.3):
  - frozen contributions are never recomputed, and their old noise weighting
    is an accepted error;
  - the add-on's final reconstruction is incremental;
  - the noise model is global;
  - new particles start from the union's committed global sigma.

## 0. Summary

- **On this data the cost is per stack, not per byte.**
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
- **The plan, none of which stores anything derived from the particles:**
  - **Read path, bit-exact:**
    - one open per stack per pass;
    - one read call per stack run;
    - a POSIX read-only path that opens stacks concurrently and tells the
      kernel ahead of time what the pass will read;
    - a pass-sized particle buffer in 3D, so the reconstruction reads from
      memory what the search already read;
    - optionally, stacks kept open from the search into the reconstruction.
  - **Half-sets per stack** (ruled). Reconstruction opens each sampled stack
    once per state instead of once per half, and overlapping neighbouring
    picks stay in one half.
  - **Add-on passes** (ruled): an add-on reads only its own particles.
  - **A direction for discussion:** extraction writes containers, large
    files holding many micrographs' stacks, so that any sampling touches few
    files (section 3.4). New data only; nothing is rewritten.
- **Counts.** Counted from the code for two stream-sized examples
  (section 2.5), before any gain from concurrent opens or read-ahead, the
  opens per `solve3D` iteration fall:
  - 1.6–1.7× with one open per stack;
  - 1.7–1.8× adding half-sets per stack.

  Opens per iteration stay of the order of the number of sampled particles.
  Going further means a layout on disk with fewer files (section 3.4).
- **`solve2D` on a stream pool** below `nsample`=200,000 reads every particle
  every iteration. One open per stack halves its opens, and concurrent opens
  and read-ahead overlap them.

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
| with one read call per stack run | 16,159 | 27,703 |
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

The steps below are candidates. Which come next, and in what order, is
decided from the measurements of section 4. The add-on passes depend on the
global-sigma refactoring in any case.

### 3.1 Read path (bit-exact)

These change only how bytes are fetched: `simple_matcher_ptcl_io.f90`, the
stack reader (`simple_discrete_stack_io.f90`) and the image file layer. The
cropped 2D cache keeps its own reader (`ptcl_cache_read_batch`). Its build
reads through `discrete_read_imgbatch`, so the build gains too.

**One open per stack per pass.**

- Take `ldim` and `nptcls` from the header that `dstack_io%open` reads
  anyway, and drop the separate `find_ldim_nptcls` on first touch.
- The rule in the memo comment is kept: the physical header, not `os_stk`'s
  `box`/`nptcls_stk`, decides record offsets. The check moves into the open.
- Drop the existence check after the open in `imgfile%open` (:141) for
  read-only opens. It is a second metadata operation on every open.
- No change to threading. This halves the opens in reads 1 and 2.

**One read call per stack run.**

- In each run, read the span of slices from the lowest to the highest needed
  index with one `rSlices(first, last)` into a scratch buffer. Then copy the
  needed slices into `build%imgbatch`.
- Do so when the span costs less than separate requests: the slices read
  beyond the needed ones fit within a gap budget (latency × bandwidth).
  Start at 2 MB and set it from the open-latency measurement. Otherwise read
  per particle, as now.
- On a spinning disk, passing over a short gap costs about as much as
  reading it, while a second request costs a rotation or a seek.
- On stacks of 5 the span is the whole stack. The bytes are the same, and so
  are the mode-2 and mode-12 conversions.
- This is second priority on stream data (one or two particles per stack per
  pass today), and first-order on large stacks and full sampling.

**A POSIX read-only path: concurrent opens and read-ahead.**

- **The path.** One small set of C helpers in `simple_posix.c`, bound in
  `simple_syslib.f90`, with no Fortran unit involved:
  - open read-only;
  - `pread`;
  - `posix_fadvise`;
  - close.
- **Concurrent opens.** Each stream opens, reads and closes its own stacks
  inside the parallel region, so up to `io_read_nstreams` opens are in flight
  instead of one.
  - The path touches no Fortran unit, so it does not depend on what the June
    race was.
  - `simple_stack_io_tester.f90`, which already reads in parallel
    (`test_dstack_parallel_read_*`), gains a case that opens in parallel,
    many times over, on more stacks than streams.
- **Read-ahead.** At the start of a pass the reader hands the kernel the
  pass's files (`posix_fadvise` with `WILLNEED`). The kernel queues the
  reads, orders them by position on disk and runs them while the CPUs work.
  - **Matchers and probability tables.** The sample is known before the
    first batch. Per-micrograph stacks are small (5 particles at box 336 is
    about 2 MB), so whole files are hinted and no header is needed first. A
    large stack is hinted by range once its header is known.
  - **Reconstruction.** It hints a window ahead: the next batch or two while
    the current one is gridded.
  - **Bounds.** Never more than a fraction of the memory the system reports
    available, so the hints do not evict what the run already holds.
  - **Full sampling.** This also serves `solve2D` on stream pools, where
    sampling cannot help.
- **Conversion.** It is the existing code: `simple_float16_to_float32` for
  mode 12 and the current conversions for the other modes. A mode the path
  does not handle falls back to today's reader.
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
  with the raw images of the pass's share:
  - in one sorted sweep, each stack opened once per pass;
  - indexed by particle, and held exactly as read, so the results stay
    bit-identical;
  - with the read-ahead hint issued for the whole list before the sweep
    starts.

  The search batches and the reconstruction both take their images from the
  buffer. Read 3 disappears.
- **Memory.** At box 336 in float32 a particle is about 450 KB.
  - A distributed part with 1,250 particles (10,000 sampled over eight parts)
    needs about 560 MB.
  - A shared-memory run with 10,000 particles needs about 4.5 GB.
  - The buffer is capped by the memory the system reports available. A share
    that does not fit falls back to today's batch reads, and the pool below
    then applies.
- **Scope.**
  - The buffer lives for one iteration in one process, and nothing is
    written to disk.
  - The probability table and the polish pass run in other processes and
    keep their own reads.
- **Overlap with compute** comes from read-ahead while the buffer fills, not
  from the buffer itself. Filling a buffer before the loop on its own
  changes the order of the reads, not their number or cost.
- **Measure** with the time-share measurement of section 4, which times
  reading in the search and in the reconstruction separately.

**Stacks kept open from the search into the reconstruction (optional).**

- In a `refine3D` part, reads 2 and 3 run in the same process. A pool of the
  part's open descriptors, kept from the search into the reconstruction,
  makes read 3's opens free.
- With the particle buffer this matters only for shares that exceed the
  buffer's cap.
- It is bounded by the open-file limit. When the pool is full, the least
  recently used stack is closed.
- It is worth adding only if the file-operation and open-latency
  measurements show the metadata cost still matters after the steps above.

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
- **The rule.**
  - `partition_eo` (`simple_oris_setters.f90:791`) orders the stacks at
    random with the run's random stream (`seed_rnd`).
  - It gives each stack's particles to the half with fewer particles so far,
    so the halves are equal to within one stack.
- **When it applies.**
  - Only to new partitions: `partition_eo` runs only when
    `get_nevenodd() == 0`, at the call sites in the refine2D and refine3D
    strategies, the 2D and 3D matchers, reconstruction, `calc_pspec`,
    `make_cavgs` and `solve2D`. Existing labels never change.
  - In the stream, chunk projects are partitioned as they are created, and
    the pool keeps their labels.
- **The fallback.** Per index, as now, when per stack cannot balance the
  halves:
  - the field carries no `stkind` (class-average fields);
  - a time series (`ini_tseries`, in `make_cavgs`);
  - the largest stack holds more than 1% of the field's particles (imported
    projects with a few large stacks).

  The fallback is logged.
- **FLEX.** Its repair of a degenerate split
  (`simple_flex_pca_project_gateway.f90:104`) uses the same routine, so it
  follows the same rule.
- **Not chosen.** Reading each batch once and routing it to both halves
  would avoid the label change. But it needs both half accumulators live,
  which the multistate peak-memory policy avoids.

### 3.3 Add-on passes (ruled 2026-10-12)

This is the stream's largest particle I/O. It changes the `solve3D_addon`
contract (`solve3D_addon_policy.md`, sections 3, 5, 7 and 8), and it depends
on the global-sigma refactoring (`global_sigma2_only.md`).

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
- The first add-on starts from the base run's frozen sets, computed once as
  today.

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
- **The full fit** is left to the session's final run (stream 3D ingestion
  policy).
- **What an add-on reads:** its own particles only, in the sampled stage
  iterations and the three passes at the end. Nothing of the frozen
  particles.
- **Over a session,** each particle is read a bounded number of times,
  however many add-ons follow.
- **The carried sums** are validated against the run context before a byte
  is read, as the frozen sets are today (`simple_frozen_accum.f90`). The
  manifest records them with the merged sigma2 state.

### 3.4 A layout on disk that serves any sampling (direction, for discussion)

The read-path changes make each open and each request cheaper. But opens per
iteration stay of the order of the number of sampled particles, because the
sample decides which files are touched. A layout that lets any sampling
scheme read cheaply has to remove three costs:

- **The per-file cost** (open, metadata, header): fewer files.
- **The per-request latency** (a seek, or a network round trip): the
  particles in memory or on flash. On a spinning disk no layout makes
  scattered single-particle reads cheap; only memory does.
- **The serial wait:** requests issued ahead (read-ahead), since the sample
  is known before the pass.

**Containers.** One layout removes the first cost without a rewrite of
existing data, and keeps one stack row per micrograph.

- **Writing.**
  - Extraction appends each micrograph's particles to a large container
    file instead of writing one file per micrograph.
  - There is one container per extraction job, so each file has one writer.
  - The MRC header's image count is set when the container is closed.
- **Indexing.**
  - Each micrograph keeps its own `os_stk` row, so the 1:1 relation with
    `os_mic` holds. The row's path is the container, and its particles'
    `indstk` are their image positions in it.
  - The stack-indexing policy already separates the stack row from the
    physical file (`indstk` is the physical image index). Rows sharing one
    file need a paragraph in that policy, not a format change.
- **Reading.**
  - A part process opens each container once and keeps it open. Any particle
    is then one positioned read at a known offset.
  - Several streams read one container through separate descriptors, on the
    POSIX path. A Fortran file can be connected to only one unit.
  - Opens per iteration fall from thousands to the number of containers
    touched, whatever the sampling scheme.
- **On a spinning disk,** a container written in micrograph order is
  contiguous. Scattered reads inside it are short seeks within one region,
  not jumps between thousands of files interleaved with other outputs.
- **New data only.** Existing projects keep their files. Reextraction and
  motion-model refinement, which rewrite particles anyway, write containers.

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
    the sampling to the layout (blocks) or read everything. A good layout
    gets most of the way; read-ahead and memory residency do the rest.

**Costs and open points.**

- **Shared files.** Code that treats a stack file as one micrograph's must
  be audited. Removing or rejecting a micrograph deselects its rows and
  leaves its images in the container, and the space comes back only through
  compaction.
- **Ownership.** Extraction is the writer, so this changes extraction and
  the stream, which is Joe's area.
- **Order.** It needs the POSIX path first. Today's sorted read opens one
  reader per stack row in parallel, which would connect a shared container
  to several Fortran units.
- **When to decide.** From the measurements of section 4: how much of an
  iteration is opens, and what an open costs on the file systems in use.

**Set aside:**

- one stack per stream chunk, which breaks the 1:1 `os_stk`/`os_mic`
  relation that `simple_sp_project_mic.f90`, `simple_reextract_strategy.f90`,
  `simple_refine_motion_model_strategy.f90` and `prune_project` rely on;
- a repacked copy (section 1).

### 3.5 Not addressed

Every part process re-reads the project file and references at launch, and
the master writes the project segment each iteration. Persistent workers
would remove that, and the cold memo with it. That is out of scope here.

## 4. Measurements: the basis for the next decision

These run after the global-sigma refactoring has been validated, on the code
as it then stands. They need no change to the read path. The only code is
the benchmark split of section 5, which the time-share measurement needs.
Their numbers then decide:

- which of the read-path steps are worth doing, and in what order;
- whether stacks need to be kept open from the search into the
  reconstruction;
- whether the stream-count cap on spinning disks changes;
- whether the container layout gets its own design note.

**The consolidated-stack experiment** comes before the rest.

- **The program.** `consolidate_particles` (developer visibility) drops the
  inactive rows and copies the remaining images, byte for byte, into one
  stack in the execution directory. Every stack row then points at it, which
  is the container indexing of section 3.4
  (`project_stack_indexing_policy.md`, "Stack Rows Sharing One File").
- **The reader fix it needs.** The sorted read path now keeps consecutive
  stack rows that share a file in one run, read through one unit.
- **What it measures.** Same particles, same sampling, same pixels, one file
  instead of thousands: the most any layout change can buy. It measures the
  per-file cost and the locality on disk together.
- **How to run it.**
  - Start with the share of the ribosome base run (about 29,000
    particles, 13 GB), where the probabilistic pre-step took 140 of 204 s.
  - Run the same `solve3D` settings on the original project and on the
    consolidated one, each from a cold page cache.
  - Compare `SOLVE3D_BENCH.txt` and the refine3D strategy and matcher
    benchmark files, with `iostat -x` and `free -g` alongside.
  - The program logs the total to copy at the start, then every 30 seconds
    the gigabytes copied, the share done, the rate and the minutes left. At
    the end it logs the gigabytes copied and the time the copy took.

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
  with 8 streams. This sets the span-read gap budget and says how much
  concurrent opens buy.
- **Stream-count sweep.** `SIMPLE_IO_NSTREAMS` at 1, 2, 4, 8 and 16 on the
  ribosome's spinning disk, with `free -g` and `iostat -x` alongside.
- **The ribosome with whole-stack runs** (`solve3D_chunks`, `76c37002b`).
  Does the page cache hold a run's particles after their first read?

## 5. Benchmarks

- **Where:** `REFINE2D_BENCH_ITERnnn.txt`, the refine3D matcher's
  per-partition files, the refine3D strategy file and `SOLVE3D_BENCH.txt`,
  in SIMPLE's format (CONTEXT, TIMINGS, % accounted for).
- **Timings:** particle reading, separate from particle preparation.
- **Counts, in the context sections, for each pass:**
  - stacks opened;
  - read requests;
  - bytes read;
  - bytes hinted.
- **Rule:** a gain is claimed only from the same run measured before and
  after.

## 6. Tests

- **The particle buffer.**
  - With the buffer and without it, the search and the reconstruction
    produce bit-identical outputs.
  - A share above the cap falls back to batch reads, with the same outputs.
  - The reconstruction makes no read from disk when the buffer holds the
    share, checked by the open counter.
- **The read path.**
  - `build%imgbatch` is bit-identical between the old and new paths on a
    fixture with stacks of 1, 5 and 200 particles, on sorted and unsorted
    index lists, in mode 2 and mode 12.
  - One open per stack per call, checked by an open counter in `dstack_io`.
  - A dimension mismatch still stops.
  - The POSIX path opens in parallel, many times over, on more stacks than
    streams.
  - Read-ahead changes no byte read.
- **Half-sets per stack.**
  - The halves differ by at most the largest stack.
  - Every particle of a stack carries one label.
  - Each fallback applies: no `stkind`, a time series, a stack above 1% of
    the field.
  - Existing labels are untouched.
  - A fixed seed gives the same split.
- **Add-on passes.**
  - When no sigma changes, the carried sums equal sums recomputed from the
    frozen particles.
  - The shipped map from carried sums equals a full reconstruction of the
    union under the same weights.
  - A carried-sum file of another run context is refused.
  - The merged sigma2 state holds the frozen particles' committed records
    unchanged.
- **Lightweight gates.** `scripts/check_test_registry.py` runs whenever a
  suite table or the test UI changes.

## 7. Order of work

1. **The global-sigma refactoring is validated** first
   (`global_sigma2_only.md`).
2. **Measurements** (section 4): the benchmark split, then the file
   operations, the time share, the open latency, the stream-count sweep and
   the whole-stack ribosome run.
3. **Decision.** The numbers are reviewed, and the next steps are chosen from
   the candidates below. The order shown is provisional:
   - one open per stack per pass;
   - one read call per stack run;
   - the POSIX read-only path: concurrent opens, then read-ahead;
   - the pass-sized particle buffer in 3D;
   - stacks kept open from the search into the reconstruction;
   - the container layout, as its own design note.
4. **Half-sets per stack** are ruled and do not wait for the numbers. They
   land with the first step chosen, and their effect on reconstruction opens
   is measured with it.
5. **Add-on passes,** as a separate track that builds on the global-sigma
   refactoring.

Each step ends with a review of its numbers against the measurements of
step 2.

## 8. Open decisions

1. **Half-sets per stack:** is the 1% largest-stack threshold for the
   per-index fallback right?
2. **The container layout** (section 3.4): worth its own design note, with
   Joe, once the measurements are in?

## 9. Code map

| Area | Files |
| --- | --- |
| Read path | `src/main/strategies/search/simple_matcher_ptcl_io.f90` (`discrete_read_imgbatch`, `stk_dims`, read counters), `src/fileio/simple_discrete_stack_io.f90`, `src/fileio/simple_imgfile.f90`, `src/fileio/simple_imghead.f90`, `src/main/image/simple_image_io.f90` |
| Particle buffer | `simple_matcher_ptcl_io.f90` (the buffer, filled in one sweep), `simple_strategy3D_matcher.f90` (filled before the search loop), `simple_matcher_3Drec.f90` and `simple_rec3D_pcg_strategy.f90` (read from it) |
| POSIX path | `src/fileio/simple_posix.c`, `src/fileio/simple_syslib.f90` (`io_read_nstreams`, the memory-available query, the new helpers) |
| Read-ahead and benchmarks | `simple_strategy2D_matcher.f90`, `simple_strategy3D_matcher.f90`, `simple_commanders_prob.f90`, `simple_matcher_3Drec.f90`, `simple_refine2D_strategy.f90`, `simple_refine3D_strategy.f90` |
| Half-sets | `src/main/ori/simple_oris_setters.f90` (`partition_eo`) |
| Add-on passes | `simple_commanders_solve3D.f90`, `simple_solve3D_utils.f90`, `simple_frozen_accum.f90`, `src/main/simple_final_rec.f90`, `exec_bootstrap_rec3D` (`simple_commanders_refine3D.f90`), `src/main/sigma2/` |
| Policies | `doc/policies/3D/solve3D_addon_policy.md` (sections 3, 5, 7, 8), the documents that describe half-set assignment, a short policy for the particle read path |
| Tests | `src/fileio/simple_stack_io_tester.f90`, the `oris` setter testers, `src/main/volume/simple_frozen_accum_tester.f90` |
