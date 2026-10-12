# Particle I/O on per-micrograph stacks: reading, sampling, half-sets and add-on passes

Implementation note, 2026-10-12. **Plan; nothing implemented.**

It merges two notes of the same day:

- the few-particle-stacks proposal (`665ac2731`), written against master
  `32a495cef`;
- `particle_io_handling.md` from the agent checkout, now retired.

RELION references are to its master at `05cd1ad9`. Nothing here touches how
sampling units are built (`simple_view_partition_sampling.f90`).

The plan was prompted by the stream ribosome project: thousands of
per-micrograph stacks with a handful of particles each, on which `solve2D`
and `solve3D` spend their time in file metadata. Companion documents:

- [particle_cache_policy.md](../../policies/2D/particle_cache_policy.md);
- [solve3D_policy.md](../../policies/3D/solve3D_policy.md);
- [solve3D_addon_policy.md](../../policies/3D/solve3D_addon_policy.md);
- the global-sigma refactoring outline (`global_sigma2_only.md`).

Rulings so far (Hans, 2026-10-12):

- **No particle is rewritten into another file.** No cache, raw or cropped,
  is added. The existing cropped 2D cache and the extraction format stay as
  they are.
- **Half-sets per stack** for new partitions (section 3.3).
- **Add-on passes** (section 3.4):
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
    - optionally, stacks kept open from the search into the reconstruction.
  - **Half-sets per stack** (ruled). Reconstruction opens each sampled stack
    once per state instead of once per half, and overlapping neighbouring
    picks stay in one half.
  - **Stack-coherent sampling**, behind a parameter. A sample of n particles
    then comes from about n/m stacks (m particles per stack) instead of close
    to n. This changes which particles are drawn, not how many, nor how often
    over an epoch.
  - **Add-on passes** (ruled): an add-on reads only its own particles.
- **Counts.** Counted from the code for two stream-sized examples
  (section 2.5), under `balance=none`, before any gain from concurrent opens
  or read-ahead, the opens per `solve3D` iteration fall:
  - 1.6–1.7× with one open per stack;
  - 3.8–5.7× adding stack-coherent sampling;
  - 5.1–7.6× adding half-sets per stack.

  With the default `balance=cavg`, stack-coherent sampling gives less. A
  sampler dry run on the real project decides.
- **`solve2D` on a stream pool** below `nsample`=200,000 samples every
  particle every iteration, so sampling cannot help there. One open per stack
  halves its opens, concurrent opens and read-ahead overlap them, and
  anything beyond that needs fewer files (not proposed).

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
- Stack-coherent sampling alters which particles are drawn.
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
| stack-coherent sampling alone, `balance=none`, whole stacks (n/m per pass, both halves) | 10,000 | 12,000 |
| one open per stack + stack-coherent sampling | 6,667 | 8,000 |
| ... + half-sets per stack | 5,000 | 6,000 |
| particle read calls today (3n) | 30,000 | 30,000 |
| with one read call per stack run | 16,159 | 27,703 |
| ... + stack-coherent sampling (+ half-sets per stack) | 6,667 (5,000) | 8,000 (6,000) |
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

**Stacks kept open from the search into the reconstruction (optional).**

- In a `refine3D` part, reads 2 and 3 run in the same process. A pool of the
  part's open descriptors, kept from the search into the reconstruction,
  makes read 3's opens free.
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

### 3.2 Stack-coherent sampling (sampler only; behind a parameter)

This couples sampling to stacks. The quotas and tier rules stay. Only the
choice inside the cutoff tier changes, from particle-uniform to stack-major:

1. **Forced picks.** Per unit (one unit for `balance=none`), compute the
   quota and cutoff tier as now, and take every particle below the cutoff
   tier, as now.
2. **Stack order.** Order the stacks: those already touched by forced picks
   first, then the rest, each group shuffled with the run's seed. `stkind` is
   in the particle segment, so the sampler needs nothing from `os_stk`.
3. **The walk.** Walk the stacks. In each, take the cutoff-tier particles of
   units whose quota is not yet filled, at most c per stack per pass.
4. **Repeat** the walk until every quota is filled. With large stacks and a
   small c this takes more than one walk, which falls back towards today's
   spread.

**What stays the same.**

- The per-unit quotas.
- Lower tiers before higher ones. Under `balance=none`, every active particle
  is drawn once before any is drawn twice: the epoch property of
  `sample4update_cnt`.
- The sample is drawn once in the master (`exec_prob_align`,
  `master_draws_sample3D`), and the parts reproduce it from the `sampled`
  flags, so nothing changes in the parts.

**Where it applies.**

- `sample4update_cnt`, `sample4update_rnd` (the sticky `solve2D` stage), the
  non-greedy `sample4update_class` and the fill-in sampler.
- A cohort drawn at a stage's first iteration (`cohort_sampling`) stays
  stack-coherent for the whole stage.
- With `frac_best` the fill works inside the eligible top fraction.
- Greedy selection already fixes the particles and is left alone.

**Expected effect.**

- About n/m stacks per pass under `balance=none` with c ≥ m.
- Under `balance=class` or `cavg`, a stack's particles fall in different
  units, so fewer can be taken together.
- A mean-field estimate gives about 45–80% of today's stack count for
  example A, and 25–50% for B. It assumes the units are filled in turn, each
  preferring touched stacks, with a fraction φ of each unit in its cutoff
  tier (φ from 1 at an epoch start towards 0.5).
- The sampler dry run replaces the estimate with the count on the real
  project.

**Statistics.**

- Within one iteration the sample is clustered by micrograph: shared
  defocus, ice and drift.
- With m ≈ 5 the clusters are small, and n/m = 2,000 micrographs per
  iteration keeps the defocus spread.
- The per-stack cap c bounds the cluster on large stacks. c = 16 would mean
  at least 625 micrographs per 10,000 particles, and 4 MB per stack at box
  256 in float32.
- Coverage over an epoch is unchanged.
- Within an iteration, particles of dense micrographs are somewhat more
  likely to be drawn. The tiers even this out across iterations.

The sampler is driven by stack membership only, never by maps, poses or
projection directions. It stays off by default until the dry run, a
benchmark run and a ruling.

### 3.3 Half-sets per stack (ruled 2026-10-12)

- **Two reasons:**
  - With per-index alternation every sampled stack has both halves, so
    reconstruction opens it once per half: 2n/m of the 4n/m opens under one
    open per stack plus stack-coherent sampling. Per stack, reconstruction
    opens each stack once per state.
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

### 3.4 Add-on passes (ruled 2026-10-12)

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

### 3.5 Fewer files (not proposed now)

At full sampling (`solve2D` on pools below `nsample`), the only way to cut
opens further is fewer files. Three routes, none proposed:

- **One stack per stream chunk** at pool import. This breaks the 1:1
  `os_stk`/`os_mic` relation that `simple_sp_project_mic.f90`,
  `simple_reextract_strategy.f90`, `simple_refine_motion_model_strategy.f90`
  and `prune_project` rely on.
- **A repacked copy:** set aside (section 1).
- **A fractional `solve2D` update for stream pools.** This would let
  stack-coherent sampling apply, at the cost of a change to the 2D sampling
  policy.

Revisit once the read path is in and the time-share measurement shows what
part of a `solve2D` iteration is still opens.

### 3.6 Not addressed

Every part process re-reads the project file and references at launch, and
the master writes the project segment each iteration. Persistent workers
would remove that, and the cold memo with it. That is out of scope here.

## 4. Measurements before code

- **Sampler dry run.** Run the current and the stack-major samplers on the
  ribosome project file for a stage's worth of iterations under
  `balance=none`, `class` and `cavg`. Record the distinct `stkind` per
  iteration, a quota check per unit, and the coverage over an epoch. No image
  is read.
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

- **The read path.**
  - `build%imgbatch` is bit-identical between the old and new paths on a
    fixture with stacks of 1, 5 and 200 particles, on sorted and unsorted
    index lists, in mode 2 and mode 12.
  - One open per stack per call, checked by an open counter in `dstack_io`.
  - A dimension mismatch still stops.
  - The POSIX path opens in parallel, many times over, on more stacks than
    streams.
  - Read-ahead changes no byte read.
- **Stack-coherent sampling.**
  - The per-unit quotas are exact.
  - The lower tiers are exhausted before the cutoff tier.
  - Under `balance=none`, every active particle is drawn once per ⌈1/f⌉
    iterations.
  - The cap c is respected.
  - Under `balance=none` at an epoch start, with stacks of m and c ≥ m, the
    sample covers ⌈n/m⌉ stacks.
  - A fixed seed gives the same sample, and the parts reproduce the master's
    sample (`sample4update_reprod`) unchanged.
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

1. **Measurements:** the benchmark split first, then the dry run, the file
   operations, the time share, the open latency, the stream-count sweep and
   the whole-stack ribosome run.
2. **One open per stack per pass.**
3. **One read call per stack run.**
4. **Half-sets per stack.** It is ruled, small, and independent of the read
   path.
5. **The POSIX read-only path:** concurrent opens, then read-ahead.
6. **Stacks kept open from the search into the reconstruction,** if step 1's
   numbers still show the metadata cost.
7. **Stack-coherent sampling** behind a parameter, after the dry run and a
   ruling.
8. **Add-on passes,** as a separate track after the global-sigma refactoring
   lands.

Each step ends with a review of its numbers.

## 8. Open decisions

1. **Stack-coherent sampling:** is clustering by micrograph within an
   iteration acceptable, and what value of c?
2. **Fewer files:** should stream `solve2D` sample fractionally, so that
   stack-coherent sampling applies?
3. **Half-sets per stack:** is the 1% largest-stack threshold for the
   per-index fallback right?

## 9. Code map

| Area | Files |
| --- | --- |
| Read path | `src/main/strategies/search/simple_matcher_ptcl_io.f90` (`discrete_read_imgbatch`, `stk_dims`, read counters), `src/fileio/simple_discrete_stack_io.f90`, `src/fileio/simple_imgfile.f90`, `src/fileio/simple_imghead.f90`, `src/main/image/simple_image_io.f90` |
| POSIX path | `src/fileio/simple_posix.c`, `src/fileio/simple_syslib.f90` (`io_read_nstreams`, the memory-available query, the new helpers) |
| Read-ahead and benchmarks | `simple_strategy2D_matcher.f90`, `simple_strategy3D_matcher.f90`, `simple_commanders_prob.f90`, `simple_matcher_3Drec.f90`, `simple_refine2D_strategy.f90`, `simple_refine3D_strategy.f90` |
| Sampling | `src/main/ori/simple_oris_sampling.f90` (`sample4update_cnt`, `_rnd`, `_class`, fill-in), the sampling parameter in `simple_parameters*.f90` and its UI |
| Half-sets | `src/main/ori/simple_oris_setters.f90` (`partition_eo`) |
| Add-on passes | `simple_commanders_solve3D.f90`, `simple_solve3D_utils.f90`, `simple_frozen_accum.f90`, `src/main/simple_final_rec.f90`, `exec_bootstrap_rec3D` (`simple_commanders_refine3D.f90`), `src/main/sigma2/` |
| Policies | `doc/policies/3D/solve3D_addon_policy.md` (sections 3, 5, 7, 8), the documents that describe half-set assignment, a short policy for the particle read path |
| Tests | `src/fileio/simple_stack_io_tester.f90`, the `oris` sampling and setter testers, `src/main/volume/simple_frozen_accum_tester.f90` |
