# Threads in SIMPLE: safety and the thread budget

Proposal note, 2026-10-18. **Nothing here is decided or implemented.** It
records a read-only inventory of how SIMPLE decides how many OpenMP threads
each process uses, and how per-thread storage is sized, on master
`9531dee62`; then it proposes one rule for safety and one policy for the
thread budget, with the questions the maintainer has to answer. Line numbers
are approximate.

Why now. Commit `43de9b340` fixed a work array in the dense table's
assignment that was sized by the global `nthr_glob` while its parallel loop
ran on `omp_get_max_threads()` threads: if the team is larger than the array,
threads write past its end. The run `speed2d_mac` then found two more places
of the same kind. The maintainer's observation: "the policy that we seem to
have adopted on local-style execution is to use the thread budget of nparts
* nthr and cap it to some upper bound. I am not sure what is the best
approach on cluster environments. This may have to be looked over across the
entire library."

Terms. The **team** is the number of threads a parallel region runs with:
the OpenMP default, set by `omp_set_num_threads`, unless the region has a
`num_threads` clause. **Per-thread storage** is any array or set of objects
indexed by the thread number (`omp_get_thread_num()+1`). The **master** is the
process that schedules the workers of a distributed run; a **part** is one
worker's share of the particles; the **local queue** runs the parts as
background processes on the same machine.

## What we learned

1. **No production path overflows today, but 36 places are safe only by
   coincidence.** Of 49 per-thread storage sites in production code, 13 size
   their storage by the team that runs them; 36 size it by `nthr_glob` or
   `params%nthr`, which equal the team in normal runs only because building
   the parameters sets both. Three windows already break that equality, and
   each is one edit away from an out-of-bounds write. Two test paths do
   overflow or race today.
2. **There is no single thread budget.** At least four rules decide how many
   threads a master uses, with three different caps (20, 32 and the number
   of cores), and every in-process sub-command leaves its own thread count
   behind when it returns.
3. **On clusters the master runs on the login node by default, and only
   SLURM is recognised.** Job scripts never set `OMP_NUM_THREADS`; the "is
   this local?" test reads a parameter that defaults to `local` even when the
   project's computing environment says otherwise; and the PBS, LSF and SGE
   backends have bugs in the resources they request.

The safety rule can be applied everywhere without changing any result. The
budget needs the maintainer's rulings, mostly for clusters.

## Facts

### Where the team is set

- **Every `params%new`** (`simple_parameters_phases.f90`, `derive_parallel_settings`,
  about lines 527-540) sets the team to `nthr` (or to `omp_get_max_threads()`
  when `nthr` is 0 or absent; the interface default is 0, "give 0 if unsure")
  and sets `nthr_glob` to the same value. `ncunits`, the number of parts that
  run at once, defaults to `nparts`.
- **`set_master_num_threads`** (`simple_exec_helpers.f90`, about lines
  250-265) sets the team of a distributed master to `SLURM_CPUS_PER_TASK` if
  set (uncapped), otherwise to `min(NTHR_SHMEM_MAX=20, omp_get_max_threads())`;
  it does not change `nthr_glob`. The refine3D, reconstruct3D, make_cavgs and
  flex_pca masters call it *before* `params%new`, which resets the team, so it
  has no effect on their own work; the number only reaches the command lines
  of their sub-commands. The refine2D and cls_expansion masters call it
  *after*, where it gives `min(20, nthr)` (no boost) unless SLURM sets more.
- **`rec3D_master_nthr`** (`simple_rec3D_pcg_strategy.f90`, about lines
  224-232) gives the reconstruction master phase, on the local queue only,
  `max(nthr, min(PCG_MASTER_NTHR_CAP=32, nparts*nthr))` capped by
  `omp_get_num_procs()`; elsewhere `nthr_master`. The PCG master raises the
  team to it and puts it back to `params%nthr` afterwards (not to the value
  before).
- **The `solve3D` class sampler** (`simple_commanders_solve3D.f90`, about
  lines 1176-1182) uses `nparts*nthr` capped by the cores, without the 32 cap.
- **Sub-commands run inside a master** (volume assembly, class-average
  assembly, the sigma2 commit, post-processing, solve3D's initial 3D) each
  build parameters, which set the team and `nthr_glob` to their own `nthr`;
  nothing restores them afterwards. Only `simple_view_partition_sampling.f90`
  (about lines 132-160) saves and restores both.
- **Stream** has two tables of threads per stage that disagree
  (`simple_stream_master_resources.f90`, about lines 40-56, against
  `simple_defs_stream.f90`, about lines 36-59, used by the interface): for
  example reference picking 8 against 4, sieving 16 against 8. Its master
  runs a persistent-worker pool of 16 workers of 16 threads.
- **FFTW** is built with its own thread library (pthreads, not OpenMP), so a
  plan's thread count is fixed when it is made: `nthr_glob` for large images
  (`simple_image_core.f90`, about lines 36-80), inside a critical section.
- **Environment variables read for threads**: only `SLURM_CPUS_PER_TASK`
  (and `OMP_NUM_THREADS` implicitly, through `omp_get_max_threads()` when
  `nthr=0`). `PBS_NUM_PPN`, `NCPUS`, `NSLOTS`, `LSB_DJOB_NUMPROC` are not read.
- `flex_pca_basis` sets the nesting depth to one (`omp_set_max_active_levels(1)`,
  about line 198) and never restores it.

### Per-thread storage

49 sites in production code: 13 sized by the team that runs them, 36 by
`nthr_glob` or `params%nthr`.

| subsystem | sized by the running team | sized by `nthr_glob` or `params%nthr` |
| --- | --- | --- |
| polar kernels | none | `heap_vars` (`simple_polarft_core.f90`), the memoisation buffers (`simple_polarft_memo.f90`), the CTF scratch (`simple_polarft_ctf.f90`) |
| probabilistic tables | the 3D dense assignment (since `43de9b340`) | the 3D dense fill and its seed bank, the neighbourhood table's buffers and work arrays, the 2D table's fill and its assignment's work array (the twin of the fixed bug, `simple_eul_prob_tab2D.f90`) |
| 3D and 2D search | the particle batch's work arrays, the reference-volume slots | the 3D and 2D search objects (`simple_strategy3D_alloc.f90`, `simple_strategy2D_alloc.f90`), the 3D matcher's counters, the particle images of a batch, the module-level scratch images (`simple_image_geom.f90`), the 2D utilities, the correlation matrix |
| class averaging and reconstruction | the PCG reconstructor's per-thread sums | the class averager's padded images and stacks (`simple_classaverager_restore.f90`), the padded image heap of the 3D reconstruction, the volume interpolation's images, the symmetry search, the power-spectrum sums |
| image and FFT | rotational averaging and the stack reader (both pinned with `num_threads`), FFT plans (nothing indexed) | the non-local-means filters, vector quantisation |
| preprocessing, extraction, FLEX, PCA, clustering | kernel PCA, difference-map graphs, the Gaussian mixture, the FLEX embedding, basis and probe fit | CTF estimation, motion correction, micrograph processing, power spectra, class-average quality features, nanoparticle utilities, the particle extractor, the FLEX PCG pool |

The three windows in which the team and `nthr_glob` differ in production:

- **The PCG master's boost** raises the team up to 32 on the local queue while
  `nthr_glob` stays at `nthr`; its nested even/odd sections run their own
  teams. The per-thread storage it reaches is sized at use, so it is safe
  today.
- **`set_master_num_threads` after the parameters** (refine2D, cls_expansion)
  with `SLURM_CPUS_PER_TASK` above `nthr`. Only scheduling and file merges run
  in that window today.
- **Sub-commands that never restore**: after volume assembly or the sigma2
  commit the master's team equals the sub-command's `nthr_glob`, not the
  master's own `params%nthr`; for example the kept reprojection builder of
  the 3D master builds `heap_vars` sized by the outer `nthr` next to buffers
  sized by the sub-command's.

Two test paths fail by this reading (to be reproduced before they are
fixed): the 2D table tester calls the assignment with parameters never built
(`nthr_glob` = 1), so with more than one thread it writes past the work
array when run alone (`suite=2d_probability_table_io`); the class-average
state-matching tester reprojects before any parameters exist, so all threads
share one padded image and race.

### Clusters

| backend | where the master runs | what a part requests | gaps |
| --- | --- | --- | --- |
| local | the user's machine, beside the parts | nothing (background processes); `ncunits` parts at once, each `nthr` threads | no budget: only a warning when `ncunits*nthr` exceeds the cores, not counting the master |
| SLURM | the user's shell (login node) unless the user wraps it (`script=yes`, or the NICE template, which gives the stream master `-c 4`) | `--cpus-per-task` (`qsys_nthr` or the master's resolved `nthr`), `--mem`, `--time` | scripts do not export `OMP_NUM_THREADS` |
| PBS | as SLURM | `-l nodes=1:ppn=N`, memory, walltime | memory unit bug: 16,000 MB is written as `16kb` (`simple_qsys_pbs.f90`, about lines 76-80); walltime in `0-H:M:S`, not a PBS format; Torque syntax only |
| LSF | as SLURM | `-n N`, `span[hosts=1]`, memory, `-W` | the queue `cryoem_cpu` is hard-coded (`simple_qsys_lsf.f90`, about line 61); `-W` in the wrong format |
| SGE | as SLURM | no CPU, memory or time request (the CPU key is not mapped, so its binding line is never written) | every part runs `nthr` threads on one slot |
| persistent workers | as their base backend | one job per worker of `nthr` slots; tasks share the job's single memory request | the pool is 16 x 16 in stream; a task's slots are not checked against its own `nthr` |
| coarray | the user's process, blocked on the launcher | nothing from a scheduler | no rule for images times `nthr` |

Two further cluster risks:

- **Two sources of truth for the backend.** The queue environment takes the
  backend from the project's computing environment, but the "is this local?"
  tests of the reconstruction master and the class sampler read
  `params%qsys_name`, which defaults to `local` and is not taken from the
  project. A SLURM project run without `qsys_name=` on the command line
  therefore gets the local boost (up to 32 threads) on the login node.
- **`nthr=0` is resolved on the master's host.** Each part then asks the
  scheduler for the login node's core count.

## Proposal

### A. One safety rule (identical results; before the release)

1. **One place sets the team.** A small set of procedures (for example
   `set_team(n)`, returning the previous value, and `restore_team(previous)`)
   sets the OpenMP default and `nthr_glob` together. `params%new`,
   `set_master_num_threads`, the PCG boost and the sub-commands run inside a
   master go through it, and every sub-command restores the master's team when
   it returns, as `simple_view_partition_sampling.f90` already does.
2. **Per-thread storage is sized by the team that will run it**, read once,
   outside parallel regions, from `omp_get_max_threads()` (one helper, for
   example `team_size()`), never from `nthr_glob` or `params%nthr`. About 35
   local scratch sites change their size in place.
3. **Persistent per-thread stores record their width** (about 12 types:
   `heap_vars` and the memoisation buffers, the 3D and 2D search objects, the
   two candidate buffers, the class averager, the padded image heap, the
   particle extractor, the module-level scratch images, the FLEX pool, the
   symmetry search) and check `omp_get_max_threads() <= width` before their
   parallel users, or grow outside parallel regions; their `kill` loops over
   what was allocated.
4. **FFTW plans** take their thread count from the team at planning time, as
   they do now through `nthr_glob`, which rule 1 keeps equal to the team.
5. **The two failing tests** build parameters (or set the team) before they
   run.
6. **Tests**: every changed region is checked at the machine's full thread
   count, repeated, bit for bit against one thread where the result is
   deterministic; one multi-threaded run per workflow catches races; one
   tester runs a master's sub-command sequence and checks that the team and
   `nthr_glob` are restored.

Size: about 45 allocation edits, 12 guards, the setter at about 8 call sites,
and the tests. One unattended run.

### B. One thread budget (needs rulings)

One function decides the threads of every process and phase, from the
parameters and the effective backend, and every caller asks it instead of
applying its own rule. A starting point for the rulings:

| process or phase | local queue | cluster |
| --- | --- | --- |
| shared-memory run (no `nparts`) | `nthr`, capped by the cores | `nthr`, capped by the job's CPU set (`omp_get_num_procs()`) |
| worker (part) | `nthr` | `nthr`; the job requests `nthr` CPUs and the script exports `OMP_NUM_THREADS=nthr`; `nthr=0` is resolved on the worker from its own CPU set, not on the master's host |
| master while parts run | its own `nthr` (it only polls) | the same, within its allocation |
| master-only phases (assignment, volume and class-average assembly, NU filter, sigma2 commit, post-processing) | the whole local budget: `min(CAP, nparts*nthr, cores)`, one constant `CAP` for every phase | the master's own allocation (`omp_get_num_procs()`, which follows the scheduler's CPU set; the scheduler variables only as an upper bound) |
| master on a login node (no allocation) | not applicable | a small fixed cap (for example 4), with a log line recommending `script=yes` or a template that runs the master inside a job |

And around it:

- **"Is this local?"** is decided from the effective backend of the queue
  environment, not from `params%qsys_name`.
- **One cap constant** replaces `NTHR_SHMEM_MAX` (20) and
  `PCG_MASTER_NTHR_CAP` (32), and the class sampler's copy of the rule.
- **The local queue**: either a hard limit (`ncunits*nthr` plus the master's
  threads at most the cores) or the present warning, counting the master.
- **Stream**: one table of stage resources; on a cluster the forked stages
  are capped by the master's allocation.

### C. Scheduler backends (bugs; before the release if they stay)

- PBS: memory in MB or GB with the right unit; walltime as `HH:MM:SS`;
  `select=1:ncpus=N` for PBS Pro as well as Torque's `nodes=1:ppn=N`.
- LSF: the queue from `qsys_partition` instead of the hard-coded name; `-W`
  in minutes or `HH:MM`.
- SGE: a CPU request (`-pe <environment> N`) and the queue from
  `qsys_partition`.
- Or drop PBS, LSF and SGE from release 4 if nobody uses them (their bugs
  suggest so).

## Decisions for the maintainer

1. May the master run on a login node? If yes, should it be capped small
   (rule above) and the master-only phases run there anyway, or should those
   phases be submitted as one job of `nparts*nthr` CPUs (next release)? If
   no, should SIMPLE refuse to run a cluster master outside an allocation?
2. The master's budget on a cluster: the job's CPU set (`omp_get_num_procs()`)
   for every backend, or each scheduler's variables?
3. The local queue: hard limit or warning?
4. One cap constant: which value (20, 32, or the cores)?
5. Do PBS, LSF and SGE stay in release 4?
6. Stream on clusters: should the master's job cover the forked stages, or
   should the stages be capped by it? Should a persistent worker's memory
   request scale with its slots?

## Before the release, and later

- **Before the release**: section A (identical; one run with stress tests);
  in section B the effective-backend test for "is this local?", the job
  scripts exporting `OMP_NUM_THREADS`, `nthr=0` resolved on the worker, one
  budget function with one cap; section C or the removal of the unused
  backends.
- **Next release**: master-only phases submitted as jobs on clusters;
  persistent-worker memory scaled by slots; stream resources tied to the
  master's allocation.
