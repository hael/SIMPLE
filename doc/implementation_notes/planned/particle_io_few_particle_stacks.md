# Particle reads on few-particle stacks: one open per stack, stack-coherent sampling (proposal, 2026-10-12)

Implementation note, 2026-10-12. **Proposal for review; nothing implemented.** Written against
master `32a495cef`; RELION references are to its master at `05cd1ad9`. Nothing here touches how
sampling units are built (`simple_view_partition_sampling.f90`). Prompted by the stream ribosome
project: thousands of per-micrograph stacks with a handful of particles each, on which `solve2D` and
`solve3D` spend their time in file metadata. Companion policies:
[particle_cache_policy.md](../../policies/2D/particle_cache_policy.md),
[solve3D_policy.md](../../policies/3D/solve3D_policy.md).

## 0. Summary

- On this data the cost is per stack, not per byte. Each sampled particle typically sits in a stack
  no other sampled particle touches, and the first touch of a stack in a part process costs two
  `inquire` stats, two opens, two 1 KB header reads and two closes before one particle is read. The
  opens run serially in the master thread (deliberately, since a June race fix); only the particle
  reads use the parallel streams.
- A distributed `solve3D` iteration in probabilistic mode reads its sample three times in two
  freshly launched processes (prob_tab, then search and reconstruction). Reconstruction reads once per
  (state, half), and `eo` alternates by particle index, so a stack is usually opened in both halves.
- Four changes, none of which stores anything derived from the particles:
  - **P1** one open per stack per pass (read path only, bit-exact), then concurrent opens through a
    thread-safe read-only path;
  - **P2** one read call per stack run (read path only, bit-exact);
  - **P3** stack-coherent sampling: the cutoff updatecnt tier is filled stack by stack, so a sample of
    n particles comes from about n/m stacks instead of close to n when the sample is sparse (changes
    which particles are drawn, not how many or how often over an epoch);
  - **P4** half-sets per stack instead of per particle index (a gold-standard decision, optional).
- Counted from the code for two stream-sized examples (section 2.5): opens per `solve3D` iteration
  fall 1.6-1.7x with P1a, 3.8-5.7x with P1a+P3, 5.1-7.6x with P1a+P3+P4 under `balance=none`, before
  any gain from overlapping opens (P1b). With the default `balance=cavg` P3 gives less; a sampler dry
  run on the real project decides (M1).
- `solve2D` on a stream pool below `nsample`=200,000 samples everything every iteration. Sampling
  cannot help there; P1a halves the opens, P1b overlaps them, and anything beyond that needs fewer
  files (P5, not proposed now).

## 1. Why this is not the particle cache

`doc/policies/2D/particle_cache_policy.md` says only that caching "never worked reliably for 3D"
and was removed from those paths (`c5dae62b7`). What the policy and the history describe is a store
of noise-normalized, Fourier-cropped particles at `box_crop`: it tied a run to one crop (`666e95960`
fixed `solve3D` at the final ladder crop to avoid rebuilding at every stage), changed the
reconstruction numerics (the removed 3D text: "the only numerics delta versus uncached execution is
the edge taper acting at `box_crop`"), and needs a validity key, an owner, cleanup, a space budget,
rank uniformity and a `cache_dir` every node can see.

Nothing below has any of that. No file is written; the same bytes reach the same code at the same
box; there is no key, owner, budget or uniformity rule. P1 and P2 change how the bytes are fetched,
P3 changes which particles are drawn, P4 changes the half-set labels of new projects.

## 2. What an iteration does to the file system today

### 2.1 Layout

Extraction writes one stack per micrograph (`ptcls_from_<mic>`,
`simple_extract_strategy.f90:360`, registered by `add_stk` at :383). Particles reach their pixels
through `stkind`/`indstk` in `os_ptcl2D`/`os_ptcl3D` and the `os_stk` row (`map_ptcl_ind2stk_ind`,
`simple_sp_project_ptcl.f90:8`).

### 2.2 Passes per iteration (distributed, probabilistic mode)

`refine3D_distr_strategy%execute_iteration` (`simple_refine3D_strategy.f90:1187`):

1. `exec_prob_align` (`simple_commanders_prob.f90:185`) draws the sample over the whole project in
   the master (:213-215), writes it to the project, and launches `prob_tab` parts
   (:228-231). Each part reads its sampled particles (`build_batch_particles3D`, :93). **Read 1.**
2. `gen_scripts_and_schedule_jobs` (:1289) launches the `refine3D` parts. Each reproduces its share
   (`sample4update_reprod`) and reads it in batches of `nthr*BATCHTHRSZ`
   (`simple_strategy3D_matcher.f90:370`). **Read 2.**
3. The same process reconstructs (`write_partial_recs`, :475): `calc_3Drec`
   (`simple_matcher_3Drec.f90:67-79`) or the PCG worker (`simple_rec3D_pcg_strategy.f90:291`, :403)
   loops over states and halves and reads each (state, eo) group separately. **Read 3, once per
   (state, half) group a stack has sampled particles in.**
4. With the polish pass (`pose_cont=yes`, :1298-1306) the discrete pass skips reconstruction
   (`volrec=no`, :1288) and a further set of parts reads the sample again and reconstructs.

`solve2D` per iteration: `prob_tab2D` parts (read 1) and `refine2D` parts, whose search batch also
feeds the class sums (read 2).

### 2.3 What one pass does per stack

`discrete_read_imgbatch` (`simple_matcher_ptcl_io.f90:189`) takes the sorted path
(`read_sorted_stack_runs`, :280) because every sampler returns sorted indices. Per stack run:

- `stk_dims` (:38, called at :319) memoizes dimensions per process. On first touch it calls
  `find_ldim_nptcls` → `get_mrcfile_info` (`simple_imghead.f90:1269-1277`): `inquire(exist)`, open,
  1 KB header read, close.
- `dstack_io%open` (`simple_discrete_stack_io.f90:94`) → `imgfile%open`
  (`simple_imgfile.f90:92`): open, `inquire(exist)` (:141), header read again (:146).
- One `read(pos=)` per particle (`rSlices(i,i)` via `read_single_mrc_image`), then close.
- The open loop (:331) and the close loop (:343) run in the master thread. Only the particle reads
  run in the `omp parallel do` (:333), on at most `io_read_nstreams` threads: 16 on local
  non-rotational disks, 8 on network file systems and macOS, 2 on spinning disks, or
  `SIMPLE_IO_NSTREAMS` (`simple_syslib.f90:436-484`).
- The serial open is deliberate: `c1f7470ed` ("fixed parallel i/o race", 2026-06-08) added the
  open before the parallel loop, a `critical(dstack_io_open)` around lazy opens inside parallel
  regions, and a test that opens serially and reads in parallel
  (`production/tests/simple_test_discrete_stack_io.f90`).

In distributed execution the part processes are new for every pass, so the memo is cold at the start
of every `prob_tab` and `refine3D` process: **2 opens per touched stack in reads 1 and 2, 1 per
(state, half) group in read 3.** With m = 5 and ~1 sampled particle per touched stack, the particle
read is a small fraction of the file operations around it.

### 2.4 Sampling

- `solve3D`: `nsample` = 10,000 by default (`NSAMPLE_SOLVE3D_DEFAULT`,
  `simple_solve3D_controller.f90:55`); `balance` defaults to `cavg` when the project has a 2D
  solution and no `vol1` (`class` with `vol1`, `none` without a 2D solution; `set_balance_default`,
  `simple_commanders_solve3D.f90:1833`), so on stream data the sampler is the nested quota of
  `sample4update_class` (`simple_oris_sampling.f90:242`); `balance=none` uses `sample4update_cnt`
  (:183). Both exhaust lower updatecnt tiers and draw the cutoff tier uniformly over particles. The
  result comes out in index order (`hpsort` in `_cnt`, packed in index order in `_class`), so reads
  already go stack by stack: what is not coupled to the stacks is which particles are drawn.
- `solve2D`: `nsample` = 200,000 (`NSAMPLE_DEFAULT_2D`, `simple_solve2D_controller.f90:26`); below
  that the update fraction is 1 and every particle is read every iteration.
- Half-sets: `partition_eo` (`simple_oris_setters.f90:791`) alternates `eo` by particle index, so
  any stack with two or more consecutive sampled particles has both halves.

### 2.5 Counts per iteration

S stacks of m particles, sample n, f = n/(Sm). Stacks touched per pass T = S(1-(1-f)^m).
Reconstruction opens R = S[(1-(1-f)^e) + (1-(1-f)^o)] with e, o the even and odd counts in a stack.
Distributed, probabilistic mode, one state, no polish; opens summed over all parts.

| | A: 5,000 stacks × 6 = 30k ptcls | B: 40,000 stacks × 5 = 200k ptcls |
| --- | --- | --- |
| f (n = 10,000) | 0.33 | 0.05 |
| stacks touched per pass, today | 4,561 (91%), 2.2 ptcl each | 9,049 (23%), 1.1 ptcl each |
| opens per `solve3D` iteration, today (2T + 2T + R) | 25,281 | 45,800 |
| P1a (T + T + R) | 16,159 | 27,703 |
| P3 alone, `balance=none`, whole stacks (n/m per pass, both halves) | 10,000 | 12,000 |
| P1a + P3 | 6,667 | 8,000 |
| P1a + P3 + P4 | 5,000 | 6,000 |
| particle read calls today (3n) | 30,000 | 30,000 |
| with P2 | 16,159 | 27,703 |
| with P2 + P3 (+ P4) | 6,667 (5,000) | 8,000 (6,000) |
| `solve2D`, full sampling: opens today → P1a | 20,000 → 10,000 | 160,000 → 80,000 |

More states add reconstruction opens where a stack's sampled particles fall in several states; with
fractional state weights (`group_pinds_by_weights`, `simple_matcher_3Drec.f90:189`) a particle is
read for every state it has weight in. The polish pass adds one pass. Each open also carries an
`inquire` stat and a header read. These are counts from reading the code, not measurements; M2
checks them.

## 3. Proposal

### P1. One open per stack per pass, then concurrent opens (read path, bit-exact)

- **P1a.** Take `ldim`/`nptcls` from the header that `dstack_io%open` reads anyway, and drop the
  separate `find_ldim_nptcls` on first touch. The rule in the memo comment (the physical header, not
  `os_stk` `box`/`nptcls_stk`, decides record offsets) is kept: the check moves into the open. Check
  whether the `inquire(exist)` after the open in `imgfile%open` (:141) can be skipped for read-only
  opens; it is a second metadata operation on every open. No change to threading.
- **P1b.** Overlap the opens. They are serial on purpose (`c1f7470ed`, section 2.3), so this needs a
  read-only path that is safe to open from several threads: POSIX `open`/`pread`/`close` through
  `simple_posix.c`, one descriptor per stream, no Fortran unit involved. Then each stream opens,
  reads and closes its own stacks inside the parallel region, and up to `io_read_nstreams` opens are
  in flight instead of one. Whether the June race was in the Fortran unit handling or elsewhere has
  to be established before P1b.
- Effect: P1a halves the opens in reads 1 and 2; P1b divides open latency by up to the stream count.
  Both apply to every program that reads particles, `solve2D` included.

### P2. One read call per stack run (read path, bit-exact)

In each run, read the span of slices from the lowest to the highest needed index with one
`rSlices(first, last)` into a scratch buffer and copy the needed slices into `build%imgbatch`, when
the span costs less than separate requests: span − needed slices ≤ a gap budget (latency ×
bandwidth; start at 2 MB, set from M4). Otherwise read per particle, as now. On stacks of 5 the span
is the whole stack. Same bytes, same mode-2 and mode-12 conversion. Second priority on stream data
(1-2 particles per stack per pass today); first-order on large stacks and full sampling.

### P3. Stack-coherent sampling (sampler only)

The coupling of sampling to stacks. Quotas and tier rules stay; only the choice inside the cutoff
tier changes from particle-uniform to stack-major:

1. Per unit (one unit for `balance=none`), compute the quota and cutoff tier as now and take every
   particle below the cutoff tier (forced, as now).
2. Order the stacks: those already touched by forced picks first, then the rest, each group shuffled
   with the run's seed. `stkind` is in the particle segment, so the sampler needs nothing from
   `os_stk`.
3. Walk the stacks; in each, take the cutoff-tier particles of units whose quota is not yet filled,
   at most c per stack per pass.
4. Repeat the walk until every quota is filled (with large stacks and a small c it takes more than
   one pass, which falls back towards today's spread).

Invariants kept: the per-unit quotas; lower tiers before higher ones, so under `balance=none` every
active particle is drawn once before any is drawn twice (the epoch property of
`sample4update_cnt`); the sample is drawn once in the master (`exec_prob_align`,
`master_draws_sample3D`) and the parts reproduce it from the `sampled` flags, so nothing changes in
the parts. Applies to `sample4update_cnt`, `sample4update_rnd` (the sticky `solve2D` stage), the
non-greedy `sample4update_class` and the fill-in sampler; a cohort drawn at a stage's first iteration
(`cohort_sampling`) stays stack-coherent for the whole stage. With `frac_best` the fill works inside
the eligible top fraction; greedy selection already fixes the particles and is left alone.

Expected stacks per pass: about n/m under `balance=none` with c ≥ m. Under `balance=class`/`cavg` a
stack's particles fall in different units, so fewer can be taken together. A mean-field estimate
(units filled in turn, each preferring touched stacks, a fraction φ of each unit in its cutoff tier,
φ from 1 at an epoch start towards 0.5) gives about 45-80% of today's stack count for example A and
25-50% for B. M1 replaces the estimate with the count on the real project.

Statistics: within one iteration the sample is clustered by micrograph (shared defocus, ice,
drift). With m ≈ 5 the clusters are small, and n/m = 2,000 micrographs per iteration keeps the defocus
spread. The per-stack cap c bounds the cluster on large stacks (c = 16 would mean at least 625
micrographs per 10,000 particles and 4 MB per stack at box 256, float32). Coverage over an epoch is
unchanged. Within an iteration, particles of dense micrographs are somewhat more likely to be drawn;
the tiers even this out across iterations. Proposed behind a parameter, off by default until M1 and
a benchmark run.

### P4. Half-sets per stack (optional; a decision)

Under P3 every sampled stack has both halves, so reconstruction opens it once per half: 2n/m of the
4n/m opens under P1a+P3. Assigning `eo` per stack (stacks in random order, each to the half with
fewer particles; halves equal to within one stack) lets reconstruction open each stack once per
state. A second reason to consider it: overlapping neighbouring picks share pixels, and with
per-index alternation they sit in opposite halves; keeping them in one half removes that FSC
coupling. Overlap is also why RELION keeps all segments of a helical tube in one half-set
(`divideParticlesInRandomHalves`, `exp_model.cpp:261-317`; its error text at :400: "half-sets are
always per-filament"). Only for new partitions (`partition_eo` runs only when
`get_nevenodd() == 0`); existing labels stay. Reading each batch once and routing it to both halves
would avoid the label change but needs both half accumulators live, which the multistate
peak-memory policy avoids; not proposed.

### P5. Fewer files (not proposed now)

The only way to cut opens at full sampling (`solve2D` on pools below `nsample`) is fewer files.
(a) Writing one stack per stream chunk at pool import breaks the 1:1 `os_stk`/`os_mic` relation that
`simple_sp_project_mic.f90`, `simple_reextract_strategy.f90`,
`simple_refine_motion_model_strategy.f90` and `prune_project` rely on. (b) A raw repacked copy at
job start avoids the cache's numerics problems but brings back its lifecycle, ownership and
multi-node problems. (c) A fractional `solve2D` update for stream pools would let P3 apply, at the
cost of a 2D sampling-policy change. Revisit once P1 is in and M3 shows what share of a `solve2D`
iteration is still opens.

### Not addressed

Every part process re-reads the project file and references at launch, and the master writes the
project segment each iteration; persistent workers would remove that and the cold memo together. Out
of scope here.

## 4. Measurements before code

- **M1, sampler dry run.** Run the current and the stack-major samplers on the ribosome project
  file for a stage's worth of iterations under `balance=none`, `class`, `cavg`: distinct `stkind`
  per iteration, quota check per unit, coverage over an epoch. No image is read.
- **M2, file operations.** `strace -f -c -e trace=openat,newfstatat,statx,read,pread64,lseek,close`
  on one `prob_tab` part and one `refine3D` part of one iteration; compare with 2T, 2T and R for
  that part's range.
- **M3, time share.** The refine3D bench (`L_BENCH_GLOB`, `rt_build_batch_ptcls` against
  `rt_align`) and the `calc_3Drec` read timer (`t_read`, behind its `DEBUG` flag), on the ribosome
  project, for `solve3D` and `solve2D`.
- **M4, latency per open** on the cluster file system: open + 1 KB read + close over a few thousand
  stacks not read recently, single stream and 8 streams. Sets the P2 gap budget and says how much P1b
  buys.

## 5. Tests

- P1/P2: `build%imgbatch` bit-identical between the old and new paths on a fixture with stacks of
  1, 5 and 200 particles, sorted and unsorted index lists, mode 2 and mode 12; one open per stack per
  call (an open counter in `dstack_io`); a dimension mismatch still stops. P1b extends
  `simple_test_discrete_stack_io` to open in parallel, many times over, on more stacks than streams.
- P3: per-unit quotas exact; lower tiers exhausted before the cutoff tier; under `balance=none`
  every active particle drawn once per ⌈1/f⌉ iterations; cap respected; under `balance=none` at an
  epoch start, with stacks of m and c ≥ m, the sample covers ⌈n/m⌉ stacks; fixed seed → same sample;
  parts reproduce the master's sample (`sample4update_reprod`) unchanged.
- P4: halves differ by at most the largest stack; every particle of a stack has one label; existing
  labels untouched.

## 6. Order and open decisions

Order: M2-M4 → P1a → M1 → P3 behind a parameter → P2 → P1b → P4 if accepted.

For Hans:
1. P3: acceptable that within an iteration the sample is clustered by micrograph, and the value of c?
2. P4: half-sets per stack for new projects?
3. P5(c): should stream `solve2D` sample fractionally so P3 applies?
4. P1b: what was the race fixed in `c1f7470ed`, and is a POSIX read-only path acceptable?
