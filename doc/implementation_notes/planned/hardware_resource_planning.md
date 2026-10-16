# SIMPLE hardware resource planning

**Status:** DRAFT for user and system-administrator review, 2026-09-29

## 1. Purpose

This document defines a reproducible way to translate available hardware and a
SIMPLE workload into batch and streaming execution settings:

- `nthr`: OpenMP threads used by each worker or partition;
- `ncunits`: maximum workers or partition jobs allowed to run concurrently;
- `nparts`: total number of pieces into which the work is divided;
- `job_memory_per_task`: scheduler memory request for one concurrent worker;
- stage-specific persistent-worker counts and threads for streaming.

The immediate goal is guidance and a future planning tool. This draft does not
change SIMPLE defaults or automatically submit jobs.

Memory alone cannot determine `nparts`. CPU capacity limits concurrency, the
dataset determines useful partition granularity, and the scheduler determines
placement. GPU capacity, local storage, network throughput, and wall-time
limits may impose additional constraints. For streaming, the acquisition rate
also imposes a service-rate requirement: processing must keep up with incoming
movies without an unbounded backlog.

## 2. Existing SIMPLE behavior

SIMPLE already contains the following pieces of the resource model:

- `memreport=yes memreport_interval=<seconds>` writes per-process current and
  peak RSS samples to `memory_usage_<pid>.csv`.
- `scripts/memory_estimator.py` provides calibrated memory estimates for
  `motion_correct`, `solve2D`, and `solve3D`. The calibration data,
  targets, safety factors, and limitations are described in
  [memory estimator and benchmark guide](../../how2s/memory_estimator.md).
- `nparts` is the total number of distributed work partitions.
- `ncunits` is the maximum number of partition jobs dispatched concurrently.
- When `ncunits` is omitted, the parameter layer currently sets
  `ncunits = nparts`.
- `nthr` is the OpenMP thread count per partition job.
- The local queue backend already warns when `ncunits * nthr` exceeds the
  processor count reported by OpenMP.
- `job_memory_per_task` becomes the scheduler memory request for a distributed
  job. Its current default is 16000 MB; it is not generally derived from the
  calibrated Python estimator.

The older Fortran `simple_mem_estimator` is not the planning authority. It is
active for only a few program paths, several of its functions are placeholders
or disabled, and its coefficients are separate from the calibrated models.

## 3. Terms that must remain separate

### 3.1 Partitions and concurrent jobs

`nparts` describes work granularity. `ncunits` describes resource concurrency.
They are equal only when every part should run at once.

For example, `nparts=32 ncunits=8` creates 32 parts but runs at most eight at a
time. As each part completes, another is dispatched. This can improve load
balance without requiring resources for 32 simultaneous workers.

### 3.2 Per-process and whole-job memory

The memory target must be identified before applying any formula:

- **Single-worker peak:** one worker's allocation. This can determine
  `job_memory_per_task` and memory-limited concurrency.
- **Whole-commander peak:** the complete process already includes the relevant
  work. Do not multiply it by `nparts`.
- **Process-tree bound:** parent plus workers has already been aggregated. Do
  not multiply it again by `ncunits`.

The current calibrated targets are:

| Commander | Estimator target |
| --- | --- |
| `motion_correct` | single-worker peak RSS |
| `solve2D` | whole-commander peak RSS |
| `solve3D` | conservative parent plus largest `nparts` worker peaks |

### 3.3 Total work and active work

Partition calculations should use objects that actually perform work:
micrographs, movies, selected particles, active classes, or other workflow
units. Inactive or rejected records must not create empty parts when that can
be avoided.

### 3.4 Shared- and distributed-memory environments

The execution environment is a model input, not merely a deployment detail.
SIMPLE must distinguish at least these cases:

- **Shared-memory workstation:** one node whose CPU cores, RAM, GPUs, and I/O
  paths are shared by all active SIMPLE processes and OpenMP threads.
- **Distributed-memory cluster:** multiple scheduler-managed nodes; memory on
  one node cannot satisfy an allocation on another node, and work is divided
  among processes, jobs, MPI/coarray images, or persistent workers.
- **Hybrid cluster:** the usual cluster configuration, with shared-memory
  OpenMP execution inside each node and distributed workers across nodes.

On a workstation, `nthr` controls threads inside a worker, while several local
workers may still compete for the same memory bandwidth, RAM, filesystem, and
GPU. The planner must reserve resources for the desktop, operating system, and
other users when the workstation is not dedicated. Aggregate capacity is
bounded by the one machine:

```text
ncunits * nthr <= allocated workstation CPU threads
M_parent + ncunits * M_worker + M_cache <= usable workstation RAM
```

On a cluster, total memory across all nodes is not a single interchangeable
pool. Every worker must fit its scheduler placement, and the concurrent-worker
calculation must be performed per node before it is summed:

```text
U_node(i) = min(U_cpu_node(i), U_mem_node(i), U_gpu_node(i), U_admin_node(i))
ncunits   <= sum(U_node(i), i=1..N_nodes)
```

The parent or controller may consume memory only on a designated launch node,
so `M_parent` must be charged to that node rather than divided across the
allocation. Likewise, a GPU on one node cannot satisfy a worker placed on a
different node.

The environment also changes performance behavior. A workstation often uses
local storage but experiences interactive contention; a cluster may provide
more aggregate compute while sharing a parallel filesystem and network with
other jobs. Launch latency, metadata traffic, inter-node communication,
per-node merge work, and scheduler queue policy must be represented explicitly
in cluster measurements. Models trained only on a workstation must not be used
as cluster models without cluster validation, or vice versa.

## 4. Hardware inventory

Before choosing execution settings, record:

| Symbol | Hardware or policy value |
| --- | --- |
| `C_node` | CPU threads allocated to SIMPLE on one node |
| `R_node` | RAM allocated to SIMPLE on one node, in MiB |
| `G_node` | GPUs allocated on one node |
| `R_gpu` | usable memory per GPU, in MiB |
| `N_nodes` | nodes that may be used concurrently |
| `U_admin` | administrator or queue limit on concurrent jobs |
| `W_limit` | scheduler wall-time limit |
| `R_reserve` | RAM reserved for the OS, filesystem cache, launcher, and master |
| `D_local` | usable local scratch capacity and throughput |
| `E_kind` | execution environment: workstation, cluster, or hybrid cluster |
| `S_kind` | scheduler family or `local` |
| `S_version` | scheduler and site-adapter version |
| `S_memscope` | memory accounting scope: job, node, task, slot, or unresolved |
| `P_place` | scheduler/process placement policy across nodes, sockets, GPUs, and NUMA domains |

The inventory schema must be extensible. In addition to these common values,
it should retain processor architecture and instruction sets, sockets and NUMA
domains, measured memory bandwidth, accelerator type and supported precision,
host-to-device and device-to-device interconnects, storage class, and network
topology. Unknown fields must survive a read/write cycle so adding a new
hardware capability does not invalidate older site profiles.

Use scheduler allocations rather than the physical machine totals. On a login
node or shared workstation, use only resources explicitly available to the
run. Record whether `C_node` counts physical cores or logical hardware threads;
the selected convention must match the scheduler's CPU accounting.

A provisional memory reserve for dedicated nodes is:

```text
R_reserve = max(4096 MiB, 0.10 * R_node)
```

This is an initial administrative policy, not a calibrated SIMPLE constant.
Shared machines and filesystem-heavy stages may require a larger reserve.

### 4.1 Scheduler capability profile

Cluster planning must not assume that every site uses SLURM. Initial coverage
should include at least:

- SLURM;
- IBM Spectrum LSF and sites retaining legacy IBM queue interfaces such as
  LoadLeveler or local wrappers;
- PBS Pro, OpenPBS, and Torque-family deployments;
- Sun/Grid Engine and Univa Grid Engine deployments;
- direct local execution; and
- an administrator-configurable adapter for other schedulers.

This list defines coverage targets, not a claim that their command-line flags
are interchangeable. CPU, memory, task, node, GPU, array, and wall-time
requests have different meanings across schedulers and even across local site
configurations of the same scheduler. For example, memory may be enforced per
job, per node, per task, or per requested slot. A planner must never translate
`job_memory_per_task` into a native option until the adapter has established
the site's accounting convention.

Every scheduler adapter should expose the same normalized capability contract:

- scheduler identity, version, queues or partitions, and site configuration;
- submission, status, cancellation, and exit-status propagation;
- physical-core, hardware-thread, task, node, socket, and NUMA allocation;
- memory request and enforcement scope;
- GPU type, count, sharing, and device-binding semantics;
- wall-time, job-array, dependency, retry, requeue, and preemption behavior;
- maximum active and queued jobs, array throttling, and site fair-use limits;
- environment modules, containers, launchers, local scratch, and exported
  environment variables;
- stdout/stderr locations and the accounting records needed by telemetry.

Capabilities should be discovered where possible and confirmed by the system
administrator. Site overrides must be explicit, versioned, and retained in the
site profile. An unknown scheduler or unresolved memory/CPU convention must
disable automatic submission and return a rendered recommendation for review;
it must not silently fall back to SLURM semantics.

Scheduler support should be tracked through a coverage matrix with three
levels:

1. **Render coverage:** SIMPLE can generate and display the native request.
2. **Execution coverage:** submission, monitoring, cancellation, dependencies,
   arrays, and exit codes have passed adapter qualification.
3. **Resource coverage:** CPU, memory, GPU, placement, and concurrency requests
   have been verified against scheduler accounting on a real site.

Only resource-qualified adapters may apply automatically generated production
settings. The matrix must record scheduler and adapter versions, the tested
site, supported capabilities, known limitations, and the date of the last
qualification.

## 5. Workload inventory

Record the inputs that can materially change cost:

- commander or workflow stage;
- number of active work objects `N_work`;
- detector manufacturer and model, sensor generation, and readout mode;
- movie format, compression, bit depth or event representation, and file size;
- physical and super-resolution movie dimensions, frame count, and effective
  downscaled dimensions;
- movie acquisition rate, normal and peak inter-arrival times, burst size, and
  permitted processing lag;
- particle box size and sampled-particle count;
- number of classes, references, states, and iterations;
- sampling distance, mask diameter, symmetry, and reconstruction backend;
- `nthr`, because threads can add workspaces and FFT plans;
- CPU-only, GPU, persistent-worker, coarray, or scheduler backend;
- expected bytes read and written per work object.

Inputs outside a memory model's calibration range must be treated as
extrapolations and validated with telemetry before production use.

### 5.1 Acquisition and detector profile

The detector name alone is not a sufficient predictor. The same detector can
produce different costs in counting, super-resolution, dose-fractionated, or
event-based modes. Record the values that determine actual input volume:

- detector model and active sensor area;
- output format, including MRC, TIFF, EER, or another event representation;
- stored width and height, physical-pixel width and height, frames or events,
  and sampling distance;
- average and high-percentile compressed file size;
- sustained movie rate `lambda_avg` and peak rate `lambda_peak`;
- burst duration, microscope pauses, and the largest acceptable backlog or
  end-to-end latency.

For a movie with stored dimensions `Nx` by `Ny`, `Nf` frames, and effective
storage `Bpp` bytes per pixel, the uncompressed input size is approximately:

```text
B_movie = Nx * Ny * Nf * Bpp
R_input = lambda * B_movie
```

Event formats and compression require measured `B_movie`; their cost must not
be inferred from decoded dimensions alone. Both decoded compute cost and
on-disk I/O rate belong in the model.

## 6. Determine concurrency

Choose a candidate thread count `T = nthr`. Measure several values when
possible; the largest thread count is not necessarily the fastest or most
memory-efficient.

### 6.1 CPU limit

For independent workers on one node:

```text
U_cpu_node = floor(C_node / T)
```

Across equivalent nodes:

```text
U_cpu = N_nodes * U_cpu_node
```

For heterogeneous allocations, calculate `U_cpu_node`, `U_mem_node`, and any
accelerator limit independently for each node and sum the resulting per-node
worker capacities. Do not multiply the capacity of the largest node by
`N_nodes`.

For the local backend, require:

```text
ncunits * nthr <= C_node
```

unless deliberate oversubscription has been benchmarked and approved.

### 6.2 Memory limit for a single-worker model

Let `M_worker` be the estimator's recommended worker memory, not its raw fitted
value. Let `M_parent` be the parent/master peak retained while workers run. If
the parent is not measured, reserve it explicitly rather than assuming zero.

```text
U_mem_node = floor((R_node - R_reserve - M_parent) / M_worker)
```

`U_mem_node` must be at least one. If it is zero, reduce `nthr`, reduce input
dimensions where scientifically permitted, request a larger-memory node, or
use a stage-specific strategy. Do not solve the problem by allowing swapping.

For scheduler jobs that are placed one worker per allocation, request:

```text
job_memory_per_task >= M_worker
```

When several workers can share a node, the scheduler request and placement
policy must guarantee their aggregate memory plus the reserve fits that node.

### 6.3 Whole-commander or process-tree model

Compare the estimator's recommended total directly with the memory allocated
to the entire command. Do not multiply that value by `nparts` or `ncunits`.

The current `solve3D` model takes `partitions` and conservatively sums the
largest `nparts` worker peaks. It does not take `ncunits`; therefore it can
overestimate a throttled run where `ncunits < nparts`. That limitation should
be removed in a future calibration before automatic planning is enabled.

### 6.4 GPU and administrative limits

For stages where one worker owns one GPU:

```text
U_gpu = N_nodes * G_node
```

This must be replaced with a stage-specific value when multiple workers share
a GPU or one worker uses multiple GPUs. SIMPLE currently has no calibrated GPU
memory model, so GPU runs require measured device-memory evidence.

The final concurrency is:

```text
ncunits = max(1, min(N_work, U_cpu, U_mem, U_gpu, U_admin))
```

Omit `U_gpu` for CPU-only stages and omit any other limit that is genuinely
not applicable. Never omit an unknown constraint by silently treating it as
unlimited; report it as unresolved.

### 6.5 Streaming rate and persistent-worker pools

Streaming is a staged service, not one batch command. Preprocessing,
reference-picking/extraction, particle sieving, pool 2D, and later 3D work have
different units of work and may use different persistent or partition-worker
pools. The planner must recommend workers per stage and identify the stage
that limits sustained throughput. A global `workers` or `ncunits` value must
not conceal a slower downstream stage.

For stage `s`, measure its effective service rate with `W` workers:

```text
mu_s(W) = completed input movies per second at stage s
```

The smallest acceptable worker count is the smallest measured `W` satisfying:

```text
mu_s(W) >= safety_rate * lambda_peak
```

where `safety_rate` covers runtime variance and short acquisition bursts. This
condition is necessary but not sufficient. The worker pool must also satisfy
the per-node CPU, memory, GPU, local-scratch, and I/O constraints from section
6. Because workers contend for memory bandwidth and storage, do not assume
`mu_s(W) = W * mu_s(1)`; benchmark nearby counts until marginal throughput no
longer justifies another worker.

The planner should model the backlog explicitly:

```text
backlog(t + dt) = max(0, backlog(t) + arrivals(dt) - completions(dt))
drain_time      = backlog / max(mu_bottleneck - lambda_avg, epsilon)
```

Recommendations should include steady-state utilization, peak backlog,
predicted drain time after a burst, and reserved idle capacity. A configuration
that eventually finishes but falls continuously behind acquisition is a
failed stream configuration.

The current stream implementation already exposes different controls at stage
boundaries: preprocessing and reference-picking use partition concurrency,
particle sieving derives persistent workers from `nchunks`, and the master can
start a persistent-worker server. Calibration must preserve these distinct
semantics and test restart behavior; it must not replace them with one batch
formula.

## 7. Determine the number of parts

After concurrency is known, choose `nparts`. The initial requirements are:

```text
ncunits <= nparts <= N_work
```

More parts than concurrent slots can improve load balance and allow shorter
retries, but every part adds launch, filesystem, merge, and bookkeeping cost.

Two useful estimates are:

```text
P_time = ceil(total_estimated_work_seconds / target_part_seconds)
P_size = ceil(N_work / target_objects_per_part)
```

An initial general-purpose heuristic is:

```text
nparts = min(N_work, max(ncunits, P_time, P_size))
```

Until per-stage overhead is measured, cap routine batch runs near two to four
waves of work:

```text
nparts <= 4 * ncunits
```

This cap is a starting heuristic, not a scientific or architectural limit.
Some stages constrain or reinterpret `nparts`, and continuation files in some
workflows currently retain partition-shaped state. Check the owning workflow
before changing `nparts` across a restart.

For workflows with `nparts_chunk`, `nparts_pool`, or concurrently processed
chunks, calculate the total possible worker count. For example:

```text
total concurrent workers = concurrent chunks * nparts_chunk
```

and apply the CPU and memory limits to that total, not to either input alone.

## 8. Worked local example

Assume a dedicated workstation provides:

```text
C_node    = 32 CPU threads
R_node    = 65536 MiB
N_nodes   = 1
N_work    = 100 movies
nthr      = 4
M_worker  = 4608 MiB
M_parent  = 1024 MiB
R_reserve = 6554 MiB
U_admin   = 32
```

The `M_worker` value is the current recommendation for the documented
4096-by-4096, 16-frame `motion_correct` example; it must be recalculated for
the actual movie dimensions, frames, sampling, and threads.

```text
U_cpu = floor(32 / 4) = 8
U_mem = floor((65536 - 6554 - 1024) / 4608) = 12
ncunits = min(100, 8, 12, 32) = 8
```

A two-wave starting point is:

```text
nparts=16 ncunits=8 nthr=4 job_memory_per_task=4608
```

This is a capacity plan, not a performance optimum. A pilot should compare
throughput and peak memory with nearby thread/worker combinations.

## 9. Calibration corpus and model training

The planner must be trained from measurements covering multiple datasets,
workflow parameters, machines, and software builds. A single reference dataset
cannot distinguish a genuine resource relationship from a property of that
particular specimen or acquisition.

The calibration corpus `D_cal` is therefore an input to model fitting. Each
observation should record the following feature groups and measured outcomes:

| Feature group | Examples |
| --- | --- |
| Dataset | image dimensions, frames, particles, box size, classes, sampling, states, active work objects |
| Acquisition | detector and mode, format, compressed bytes per movie, sustained and peak movie rate, burst length and latency target |
| SIMPLE parameters | commander, `nthr`, `nparts`, `ncunits`, backend, iterations, masks, filters, cache and GPU settings |
| Hardware and scheduler | environment kind, CPUs, RAM, GPUs and GPU RAM, nodes, NUMA layout, local storage, queue limits and placement |
| Software environment | SIMPLE revision, model version, compiler, FFT and math libraries, MPI/coarray runtime, operating system |
| Outcomes | success or failure, peak worker and process-tree memory, elapsed time, throughput, CPU/GPU utilization, I/O, merge time, backlog and latency |

Training runs should deliberately vary both the data and the parameters. They
must include small, medium, and large cases; parameter combinations near
expected resource limits; and repeated runs where runtime noise is material.
Invalid combinations and resource failures are useful observations and should
be retained with their failure reason rather than silently discarded.

Models should normally be fitted per commander and execution backend. A single
universal model is acceptable only if validation demonstrates that it predicts
each supported commander as well as the dedicated models. Memory, elapsed time,
and throughput are separate targets; the setting with the lowest memory is not
necessarily the setting with the shortest elapsed time.

Streaming models should additionally be fitted per stage and worker-pool
topology. They must predict both service rate and resource use as persistent
worker count changes. The optimization target is the smallest safe pool that
meets the acquisition-rate and latency objectives, not the largest pool the
machine can launch.

Workstation, distributed cluster, and hybrid-cluster observations must be
identified explicitly. A model may share portable workload-size terms across
them, but environment-specific placement, communication, I/O, contention, and
memory coefficients require independent validation.

To avoid optimistic validation, hold out complete datasets and complete
machines. Randomly splitting rows from the same dataset and host between
training and validation would allow nearly identical runs to appear on both
sides. Every fitted model must publish:

- its training ranges and held-out error;
- the SIMPLE and software versions represented;
- its safety margin and uncertainty estimate;
- the feature values required to make a recommendation;
- an out-of-distribution rule that returns "no recommendation".

### 9.1 Installation-time site profile

At installation, the system administrator supplies or confirms the hardware
and scheduler inventory from section 4. The shipped calibration models combine
that inventory with representative workload envelopes from `D_cal` to produce
a site profile, for example:

- supported execution backends and launchers;
- conservative concurrency ceilings per node;
- baseline `nthr`, `ncunits`, and `job_memory_per_task` ranges per commander;
- GPU ownership and memory constraints;
- configurations that require a pilot before production use.

These are site defaults and feasibility limits, not fixed `nparts` values for
all future jobs. The installer does not know the dimensions or number of work
objects in a future dataset, so it cannot produce a reliable final partition
count.

### 9.2 Workload-time recommendation and local learning

When a user supplies an actual dataset and command line, the planner combines
the site profile with those workload features to recommend `nthr`, `ncunits`,
`nparts`, and `job_memory_per_task`. It should explain the limiting constraint
and identify which inputs are extrapolations.

Sites may optionally add telemetry from successful pilot and production runs
to a local calibration corpus. Retraining must be explicit, versioned, and
reproducible; a new local model replaces a shipped model only after held-out
validation shows that it is at least as safe within the site's supported
range. Raw scientific data need not be retained when the recorded features and
resource telemetry are sufficient for fitting.

### 9.3 Administrator-run reference benchmark

SIMPLE should publish a versioned reference benchmark that a system
administrator can run before selecting site defaults. The benchmark needs two
related but distinct artifacts:

1. An immutable input bundle containing openly licensed real data, compact
   desktop-sized subsets, manifests, checksums, provenance, and expected
   scientific outputs. A DOI-backed archive such as Zenodo is appropriate for
   distributing and identifying each released bundle.
2. An append-only results corpus containing the input fingerprint, SIMPLE and
   model versions, hardware and software profile, acquisition metadata,
   parameter vector, telemetry, correctness result, and failure reason for
   every run.

The input bundle is not itself the results database. Keeping them separate
allows the data release to remain immutable while benchmark observations grow
across machines and software versions. Results contributed to a shared corpus
must exclude raw user data, credentials, hostnames, and other site-sensitive
metadata unless the administrator explicitly approves them.

The normal qualification workflow should be:

1. Download a named benchmark release and verify its checksums.
2. Run its local execution adapter on a dedicated desktop or workstation,
   without requiring a scheduler or cluster account.
3. Sweep a bounded set of thread, partition, concurrent-worker, and persistent
   stream-worker settings within declared time, memory, and storage budgets.
4. Verify the expected scientific outputs while collecting resource and
   throughput telemetry.
5. Compare the observations with the shipped model and generate a versioned
   site correction with uncertainty and supported ranges.
6. Optionally validate that correction on cluster nodes using the real
   scheduler, network, and shared filesystem before enabling cluster-wide
   recommendations.

The shipped corpus supplies broad prior evidence; the local benchmark supplies
the correction for the administrator's actual machine. Extrapolation is
permitted only inside declared workload and hardware ranges. A local benchmark
that falls outside them must request another calibration point or return "no
recommendation" rather than extending a fitted curve without evidence.

## 10. Measurement and parameter-optimization framework

The project needs one framework that can run controlled experiments, collect
comparable measurements, identify missing regions of the calibration corpus,
and fit or validate recommendation models. It should extend the existing
memory-estimator scripts rather than create an unrelated performance system.

### 10.1 Framework components

The proposed framework has seven parts:

1. **Experiment manifest.** A versioned YAML or JSON file identifies the
   dataset, SIMPLE revision, command, parameter ranges, detector and acquisition
   profile, hardware requirements, local resource budget, repetitions, timeout,
   and expected outputs.
2. **Campaign generator.** It expands the manifest into explicit runs while
   respecting invalid combinations and a maximum CPU-hour or GPU-hour budget.
3. **Execution adapters.** Direct local execution is the portable default. The
   same experiment can optionally run through SLURM, PBS, LSF, SGE, persistent
   workers, or coarrays without changing its scientific inputs.
4. **Telemetry collector.** It records process-tree memory, CPU and GPU use,
   elapsed time, I/O, scheduler placement, stream arrival and completion rates,
   stage backlog, persistent-worker occupancy, exit status, and SIMPLE's own
   metrics in a common schema.
5. **Reference-data and results manager.** It downloads and verifies immutable
   benchmark bundles and writes scrubbed observations to an append-only local
   or shared results store.
6. **Feature and model pipeline.** It converts measurements into `D_cal`, fits
   commander/backend and stream-stage service models, evaluates held-out
   datasets and machines, and publishes a versioned model bundle.
7. **Advisor and report.** It combines a model bundle, site profile, and actual
   workload to explain recommended values and unresolved constraints.

Every generated run needs a stable experiment ID derived from the manifest,
dataset fingerprint, software revision, and parameter vector. Results must be
append-only so interrupted campaigns can resume without silently replacing an
earlier observation.

### 10.2 Where measurements run

Different evidence belongs at different frequencies and on different hosts:

| Layer | Where and when | Purpose |
| --- | --- | --- |
| Fast probes | build CI and optional local installation qualification; seconds | detect large regressions and characterize basic CPU, FFT, memory, I/O, launcher, and thread behavior |
| Bullet runs | a dedicated local desktop, workstation, or administrator-selected host; seconds to a few minutes | sample a small number of nearby parameter settings and estimate local scaling slopes |
| Representative pilots | the target workstation by default, or the target queue and storage path when qualifying a cluster; minutes | correct the shipped model for the actual dataset and site |
| Calibration campaigns | dedicated desktop/workstation by default; an optional controlled cluster node for cluster-specific behavior | fill the multi-dataset corpus and refit released models |
| Production telemetry | opt-in, sampled, and scrubbed of scientific data | detect drift and propose future calibration points |

CI runners are shared and noisy, so their absolute timings must not determine
production resource requests. They are useful for detecting discontinuities
and verifying that the measurement machinery still works. Published models
must rely on controlled hosts plus representative site pilots.

The portable benchmark must complete through the local adapter without a
batch queue, privileged administrator operations, or shared cluster storage.
Its manifests must cap elapsed time, peak storage, memory, and CPU/GPU use so a
system administrator can select a subset appropriate for an ordinary desktop.
Queue-based runs are additional evidence only when SIMPLE will be deployed on
a cluster. They measure placement, network, launcher, and shared-filesystem
effects; queue delay and unrelated cluster load must not define portable
kernel or CPU coefficients.

For clusters, bullet runs should include both one-node shared-memory probes and
multi-node distributed probes. This separates thread scaling within a node
from worker scaling across nodes and exposes shared-filesystem or network
bottlenecks. Workstation qualification normally omits multi-node probes but
must measure contention between concurrent local workers.

### 10.3 Low-level resource qualification

Before running a full reference workflow, SIMPLE should run a short,
scheduler-free qualification campaign that establishes a conservative operating
envelope for the local machine. These are performance and capacity probes, not
unit tests: scientific correctness is still a mandatory gate, but a shared or
slower machine must not fail merely because an absolute timing differs from a
reference host.

The campaign should isolate one scaling dimension at a time:

| Probe | Controlled sweep | Quantities measured |
| --- | --- | --- |
| Baseline | one worker and one thread | fixed startup, memory, I/O, and elapsed-time costs |
| Thread scaling | `nthr=1,2,4,...` up to physical-core capacity | throughput, CPU efficiency, memory growth, and the thread-scaling knee |
| Worker scaling | 1, 2, 4, ... concurrent workers at candidate `nthr` values | aggregate throughput, process-tree memory, I/O contention, and oversubscription |
| Partition scaling | increasing `nparts` at fixed concurrency | launch, merge, imbalance, and minimum useful task duration |
| Stream scaling | controlled movie replay while varying one stage pool at a time | per-stage service rate, utilization, backlog, latency, and drain time |
| Scheduler adapter | one tiny job, a bounded array, a dependency, and a placement probe | request semantics, task placement, accounting, exit codes, cancellation, and launch overhead |

Each point should be repeated enough to estimate runtime variability. Resource
decisions must use conservative bounds rather than the best observed run:

```text
Q_low(c) = lower confidence bound for throughput at configuration c
M_high(c) = upper confidence bound for peak process-tree memory at c
```

A configuration is inside the safe envelope only when all applicable
conditions hold:

```text
scientific correctness checks pass
M_high(c) <= (1 - reserve_mem) * usable_memory
threads(c) <= allocated_physical_cores - reserve_cores
GPU memory and ownership limits are satisfied
temporary storage remains below its reserved capacity
measured I/O demand remains below the sustainable local limit
```

For a stream configuration, also require:

```text
mu_low,s(c) >= rate_headroom * lambda_peak    for every required stage s
predicted peak backlog <= backlog_limit
predicted drain time <= drain_time_limit
```

The reserve values are policies, not universal constants. A dedicated compute
node may reserve less CPU than an interactive workstation, while an unknown or
noisy environment should use larger memory, throughput, and rate margins. The
qualification report must state the selected margins.

For a batch stage, let `Q_best` be the best conservative throughput among safe
configurations. Rather than selecting the largest configuration, form a
near-optimal set:

```text
C_near = {c in C_safe : Q_low(c) >= (1 - throughput_tolerance) * Q_best}
```

Choose from `C_near` the configuration using the fewest cores, workers, GPU
slots, and memory. This identifies the scaling knee: it preserves nearly all
measured throughput while avoiding resources whose marginal benefit is small.
For a stream stage, choose the smallest safe worker pool that meets the peak
arrival and latency contract. Nearby larger settings may be reported as burst
options, but should not become defaults without a measured benefit.

The campaign should proceed incrementally. Start from the one-worker baseline,
increase one dimension, stop before a predicted capacity boundary, and retain
the last two safe points around each scaling knee. Never deliberately drive the
machine into swapping, GPU out-of-memory, filesystem exhaustion, or sustained
stream backlog merely to discover a failure point. Boundary behavior should be
inferred from telemetry and approached with bounded steps.

The resulting low-level site profile should contain:

- recommended and maximum-safe `nthr` per representative kernel or commander;
- recommended and maximum-safe local concurrency;
- per-worker and process-tree memory bounds;
- the useful `nparts` range and minimum useful part duration;
- sustainable sequential and concurrent I/O rates;
- recommended persistent-worker counts and measured service-rate headroom per
  stream stage;
- confidence intervals, reserve policies, environment fingerprint, and the
  conditions that require requalification.

When qualifying a cluster, scheduler probes should use negligible scientific
work and safe resource requests. They should verify the resources actually
granted from scheduler accounting and the job environment, rather than testing
enforcement by intentionally exceeding memory or wall time. Queue wait time
must be reported separately from launch overhead and execution time. A failed
adapter probe blocks automatic submission for that scheduler but does not
invalidate the scheduler-free local hardware profile.

### 10.3.1 Fast probes and bullet-run scaling theory

Very fast tests can provide clues about scaling even when they are too small
to predict an entire workflow. The framework should treat these as short
probes, or "bullet runs": fire a small number of controlled measurements at
the parameter space, observe their direction and limiting resource, and use
that evidence to choose the next measurement.

Useful probes include:

- FFT throughput at representative 2D and 3D box sizes with 1, 2, 4, and 8
  threads;
- memory allocation, copy, transpose, and bandwidth at representative array
  sizes;
- MRC stack sequential and concurrent read/write throughput;
- process, MPI/coarray image, and scheduler-task startup latency;
- one small partition at several `nthr` values;
- two or more concurrent small partitions to expose memory and I/O contention;
- replaying detector-specific movies at controlled sustained and burst arrival
  rates while varying each stream stage's persistent-worker pool;
- measuring stage service rate, worker occupancy, peak backlog, and backlog
  drain time for representative movie dimensions and formats;
- GPU initialization, transfer, kernel throughput, and device-memory peaks
  where applicable.

Existing correctness tests may expose some of these operations and can emit
diagnostic measurements when doing so adds negligible work. Their pass/fail
criteria must remain correctness-based: an absolute timing from a shared CI
runner must not make a unit test fail. Stable performance probes should have
their own manifests so their inputs and repetitions are explicit.

A first-order scaling model for interpreting the probes is:

```text
T_run = T_fixed
      + W_cpu / (nthr * ncunits * efficiency)
      + nparts * T_launch
      + T_io(ncunits)
      + T_merge(nparts)

M_node = M_parent + ncunits * M_worker(nthr, workload) + M_cache
```

This is a hypothesis to fit and test, not an assumed law. From short runs, the
framework can estimate thread speedup and efficiency:

```text
speedup(T)    = elapsed(1) / elapsed(T)
efficiency(T) = speedup(T) / T
```

A sharp loss of efficiency indicates that more threads are unlikely to help;
rising elapsed time with concurrent workers indicates memory-bandwidth or I/O
contention; nearly constant time per work object under proportional resource
and dataset growth suggests useful weak scaling. Full pilots remain necessary
because startup, merging, cache effects, and contention may be absent from a
small probe.

### 10.3.2 Minimum viable SIMPLE qualification

The first executable qualification should use one deterministic synthetic
movie and the existing `motion_correct` memory harness. Motion correction is a
useful first probe because it exercises SIMPLE startup, MRC I/O, FFT work,
OpenMP scaling, and native process-memory telemetry without requiring a user
dataset. It does not characterize every commander, GPU behavior, distributed
execution, or all stream stages.

The minimum procedure is:

1. Record CPU topology, NUMA nodes, physical and logical cores, available RAM,
   swap activity, storage capacity, operating system, SIMPLE revision, and the
   executable paths.
2. Pin the probe to physical cores in one NUMA domain. Record the binding and
   avoid hardware threads for the initial sweep.
3. Generate a deterministic movie large enough that the measured command lasts
   at least about one second at the highest candidate thread count. Very short
   cases are smoke tests only.
4. Run `nthr=1,2,4,...` up to the physical cores available in that NUMA domain,
   with at least three repetitions. Require successful completion and valid
   memory telemetry for every point.
5. Select the smallest `nthr` within the agreed throughput tolerance of the
   best median result. Reject higher settings whose marginal throughput does
   not justify their additional cores or memory.
6. Run 1, 2, and then the proposed number of workers concurrently. Bind every
   worker to disjoint physical cores and verify per-worker elapsed time,
   aggregate throughput, process-tree memory, and swap activity.
7. Choose `ncunits` no larger than the measured safe concurrency. On a shared
   workstation, retain an explicit CPU and RAM reserve for the operating system
   and other users.
8. Set an initial `nparts` to one or two waves of measured concurrency when the
   workload contains enough work objects. Increase it only when a separate
   partition probe demonstrates better load balance or restart behavior.
9. Round the observed upper memory bound upward using the declared uncertainty
   and reserve policy. Do not reuse the fixture's memory request for larger
   movies; query the calibrated workload model and verify it with a pilot.
10. Repeat the selected configuration once as a final verification and save
    the manifest, raw CSV, logs, telemetry, environment fingerprint, and
    derived recommendation.

The first result is a conservative local baseline for `nthr`, `ncunits`, and
fixture memory. It is not yet a universal SIMPLE configuration. A second-level
qualification must cover representative detector dimensions, concurrent I/O,
and the stream-stage service rates required by the intended acquisition.

### 10.3.3 Initial run on the reference workstation

The minimum qualification was exercised on 2026-09-29 on a two-socket Intel
Xeon Gold 6242R workstation with 20 physical cores per socket, two hardware
threads per core, two NUMA nodes, and approximately 251 GiB of RAM. The probe
was pinned to physical CPUs in one socket and used a deterministic
2048-by-2048, eight-frame movie at 1.3 Angstrom per pixel. Each thread setting
was repeated three times.

| `nthr` | Median elapsed (s) | Worst elapsed (s) | Median speedup | Thread efficiency | Highest peak RSS (MiB) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 4.355 | 4.405 | 1.00 | 1.000 | 415.6 |
| 2 | 2.656 | 2.698 | 1.64 | 0.820 | 409.6 |
| 4 | 1.672 | 1.687 | 2.60 | 0.651 | 415.9 |
| 8 | 1.149 | 1.427 | 3.79 | 0.474 | 414.0 |
| 16 | 1.122 | 1.160 | 3.88 | 0.243 | 523.2 |

Sixteen threads improved median elapsed time by only 2.4 percent over eight
threads while consuming twice as many cores and increasing the highest
observed peak RSS by about 26 percent. Eight threads are therefore the
provisional scaling knee for this fixture.

Two concurrent eight-thread workers completed in 1.363 and 1.466 seconds. Four
workers pinned to disjoint groups of eight physical cores completed in 1.125
to 1.205 seconds, with peak RSS between 411.7 and 412.1 MiB per worker. The
four-worker result supports this provisional shared-workstation baseline:

```text
nthr=8
ncunits=4
nparts=8 when at least eight work objects are available
job_memory_per_task=1024 MiB for this 2048x2048x8 fixture only
```

This uses 32 of 40 physical cores and leaves eight cores as an interactive
reserve. The memory request is deliberately rounded well above the observed
peak, but it must not be used for production detector-sized movies. The host
had a system load near 16 during the measurements, and each concurrent-worker
configuration was run as only one wave. These results demonstrate the
procedure and establish provisional parameters; repeated concurrency waves,
I/O probes, and representative-movie pilots remain required before promoting a
site profile.

### 10.4 Filling calibration gaps

After each campaign, the framework should create a coverage report over
commander, backend, dataset scale, parameter range, and hardware class. The
next experiments should prioritize:

1. unsupported combinations required by users or administrators;
2. regions where model uncertainty is largest;
3. predicted memory or wall-time safety boundaries;
4. disagreement between the shipped model, local pilots, and production
   telemetry;
5. settings where nearby bullet runs show a change in scaling regime;
6. held-out datasets and machines needed to test generalization.

This adaptive design avoids an exhaustive Cartesian product while still
measuring dangerous boundaries. A gap is closed only after the new region has
both training observations and independent validation observations. The
campaign stops when its measurement budget is exhausted or all required
regions meet declared error and safety targets.

### 10.5 Adapting to new hardware

The framework must expect hardware classes that were not represented when a
model was released: new CPU architectures and vector widths, GPUs or other
accelerators, unified-memory systems, faster interconnects, computational
storage, and different node topologies. Recommendation logic must therefore be
capability-based rather than a list of recognized product names.

A hardware profile should describe observable capabilities and measured
behavior, including:

- available execution backends and supported numerical precision;
- physical cores, hardware threads, sockets, NUMA domains, and vector features;
- host and accelerator memory capacity, bandwidth, and allocation behavior;
- accelerator count, sharing rules, and interconnect topology;
- local and shared storage bandwidth, latency, and capacity;
- network and collective-operation characteristics;
- compiler, driver, firmware, math-library, MPI, and coarray runtime versions.

New hardware follows a staged qualification workflow:

1. **Discover.** Record the capabilities without assuming that a known model
   applies because a device name looks similar.
2. **Verify correctness.** Run the applicable fast and platform tests for the
   enabled backend. Unsupported operations disable only the affected backend.
3. **Run bullet probes.** Measure startup, FFT, memory, I/O, threading,
   accelerator transfer, and representative kernel scaling.
4. **Fit a provisional correction.** Start from architecture-independent
   formulas or the nearest validated model only when feature compatibility is
   explicit; fit hardware-specific coefficients from the probes.
5. **Run representative pilots.** Exercise supported commanders on held-out
   workloads and verify memory and wall-time safety margins.
6. **Promote the profile.** An administrator approves a versioned site model
   with declared ranges, uncertainty, and expiration triggers.

Until qualification reaches step 6, the planner may report measured
capabilities and conservative feasibility bounds, but it must label parameter
recommendations provisional or return "no recommendation". New hardware must
not be silently mapped to an old GPU, CPU, or node class.

The fitting design should share relationships that are genuinely portable,
such as array-size growth or the distinction between fixed and per-work-unit
costs, while keeping hardware-specific correction terms. This allows a new
machine to benefit from the existing corpus without pretending that its
throughput, concurrency, or memory behavior is already known. As local
measurements accumulate, adaptive campaigns should prioritize the regions
where the inherited model and observations disagree.

A hardware profile must be requalified after changes that can alter execution
behavior: major compiler or runtime upgrades, accelerator drivers or firmware,
FFT and math libraries, MPI/coarray implementations, scheduler placement,
memory configuration, storage, or network topology. Historical observations
remain in the corpus with their environment identity; they are not rewritten
as if they came from the new configuration.

## 11. Validation procedure

Before using a plan for a large production run:

1. Run the appropriate calibrated estimator and retain its JSON output.
2. Reject or review every calibration-range warning.
3. Run a representative pilot with `memreport=yes` and a short reporting
   interval.
4. Collect telemetry from the parent and every concurrent worker.
5. Confirm that observed per-process and aggregate peaks remain below the
   allocation with the declared reserve.
6. Confirm `ncunits * nthr` does not oversubscribe allocated CPUs.
7. Measure wall time per work object, I/O throughput, and merge overhead.
8. Recalculate `P_time`, `P_size`, `ncunits`, and `nparts` from the pilot.
9. For streaming, replay the expected sustained and peak acquisition rates and
   confirm that every stage remains stable, peak backlog is bounded, and the
   backlog drains within the declared latency target.
10. Repeat after changes to SIMPLE, compiler, FFT library, allocator, operating
   system, GPU stack, or important workflow settings.

The plan must fail closed: an unsupported commander or missing memory target
should produce "no recommendation" rather than a plausible-looking default.

## 12. Proposed planner output

A future resource-planning command should report both its recommendation and
the limiting factor:

```text
Commander: motion_correct
Hardware: 32 CPU threads, 65536 MiB usable RAM, 0 GPUs
Work: 100 movies
Threads per worker (nthr): 4
CPU-limited workers: 8
Memory-limited workers: 12
Administrative limit: 32
Recommended concurrent workers (ncunits): 8 [CPU limited]
Recommended total partitions (nparts): 16 [two waves]
Memory per task: 4608 MiB
Warnings: none
```

A stream recommendation should additionally report the acquisition contract
and each stage's pool rather than collapsing them into one worker count:

```text
Acquisition: detector=<model/mode>, movie=4096x4096x40, peak=0.50 movies/s
Sustainable input rate: 0.65 movies/s [30% headroom]
Limiting stage: preprocessing
Persistent workers: preprocessing=4, picking=2, sieving=2, pool2D=2
Predicted peak backlog: 8 movies
Predicted backlog drain time: 54 s
Warnings: detector mode is locally calibrated; burst duration is extrapolated
```

Machine-readable output should also include the model version, calibration
range, safety factor, reserve policy, formulas, every candidate limit, and all
warnings. This makes the recommendation auditable by a system administrator.

## 13. Gaps before automatic selection

Automatic `nparts` or `ncunits` selection should not become a default until:

- per-worker models exist for the main distributed stages;
- parent and worker overlap is measured rather than inferred;
- memory models include `ncunits` where concurrency changes the process-tree
  bound;
- stage-specific minimum useful part sizes and launch/merge overheads are
  calibrated;
- restart behavior when `nparts` changes is safe for the owning workflows;
- scheduler placement and memory semantics are validated for SLURM, PBS, LSF,
  legacy IBM queue interfaces, SGE, local, persistent-worker, and coarray
  execution;
- scheduler coverage is reported separately for request rendering, execution
  lifecycle, and verified resource semantics, with unknown adapters failing
  closed;
- GPU count and device-memory rules are available for GPU stages;
- the planner distinguishes physical cores, logical threads, sockets, and
  NUMA locality;
- filesystem capacity and I/O contention can constrain concurrency;
- recommendations are validated on more than one machine and SIMPLE build;
- a versioned, multi-dataset calibration corpus and reproducible fitting
  pipeline exist for every model used to generate automatic settings;
- installation-time site profiles and workload-time recommendations are
  validated independently;
- workstation, distributed-cluster, and hybrid-cluster behavior is represented
  explicitly in the corpus and validation reports;
- cluster capacity is computed and validated per node rather than from total
  aggregate RAM and CPU counts;
- fast probes, representative pilots, and calibration campaigns use one
  versioned experiment and telemetry schema;
- scheduler-free low-level probes establish conservative confidence bounds,
  safe capacity limits, and the throughput-scaling knee before selecting local
  defaults;
- model coverage reports identify unsupported regions and drive the next
  measurements;
- detector, movie-format, dimension, sustained-rate, and burst-rate families
  are represented in the stream calibration corpus;
- stream models distinguish stage-specific persistent-worker pools, service
  rates, backlog growth, and drain time;
- a versioned public benchmark bundle provides checksummed real inputs and
  expected scientific outputs independently of the append-only results corpus;
- the reference benchmark and bounded calibration sweeps work on a local
  desktop or workstation without a scheduler;
- cluster-specific network, placement, launcher, and shared-filesystem effects
  are qualified separately from the scheduler-free portable baseline;
- hardware profiles use extensible capability descriptions rather than fixed
  device-name tables;
- an unknown hardware class fails closed until correctness, bullet probes, and
  representative pilots establish a validated operating range.

Until those conditions are met, the planner should recommend explicit command
line values but leave the final decision with the user or system administrator.

## 14. Review questions

1. Should the planner optimize primarily for shortest elapsed time, minimum
   resource use, scheduler throughput, or a selectable objective?
2. What reserve policy is acceptable on dedicated nodes, shared workstations,
   and persistent-worker nodes?
3. Should `nparts` be allowed to change across restarts before partition-shaped
   continuation state is removed?
4. Which workflows require one part per socket, one part per node, or one part
   per GPU rather than the general worker model?
5. What minimum task duration should be targeted to avoid scheduler overhead?
6. Should the first implementation only report recommendations, or may it also
   populate `nparts`, `ncunits`, `nthr`, and `job_memory_per_task`?
7. Which representative public or synthetic datasets should define the shipped
   calibration corpus, and which dataset families must be held out?
8. May a site train local models automatically from telemetry, or must an
   administrator review and promote every fitted model?
9. Which fast probes are stable enough to run in ordinary CI, and which require
   a dedicated benchmark host?
10. What CPU-hour and GPU-hour budgets should limit adaptive calibration
    campaigns?
11. Which hardware changes invalidate only a performance correction, and which
    require complete backend requalification?
12. Who is authorized to promote a provisional hardware profile for general
    use at a site?
13. Which SIMPLE workflows support true multi-node execution, and which should
    remain constrained to one shared-memory node?
14. Which placement policies should the planner support first: one worker per
    node, one per socket, or several workers per node?
15. Which detector generations, readout modes, movie dimensions, and storage
    formats must the first stream calibration release cover?
16. What acquisition-rate safety factor, maximum backlog, and backlog-drain
    target define a successful stream configuration?
17. Which real benchmark datasets can be redistributed with clear provenance,
    licenses, checksums, and long-term versioning through Zenodo or an
    equivalent archive?
18. What elapsed-time, memory, CPU/GPU, download, and temporary-storage budgets
    keep the reference qualification practical on an ordinary desktop?
19. Who maintains the shared results schema, reviews contributed observations,
    and promotes a new version of the shipped model?
20. Which stream stages expose independently tunable persistent-worker pools,
    and which counts must remain fixed by the current implementation?
21. What throughput tolerance defines the near-optimal set from which the
    least resource-intensive configuration is selected?
22. How many repetitions and what confidence level are required before a
    low-level measurement can define a conservative site limit?
23. Which scheduler families and versions must be resource-qualified for the
    first release, and which may initially provide render-only coverage?
24. Which scheduler capabilities are mandatory for stream operation, including
    persistent workers, dependencies, requeue behavior, and array throttling?

## 15. Campaign implementation: first slice, 2026-10-15

The current checkout now implements a Python campaign entry point in
`scripts/hardware_resource_planner.py`, shared execution in
`scripts/benchmark/campaign.py`, and commander capability/argv adapters in
`scripts/benchmark/adapters.py`. The earlier conversational prototype and its
binary profile writer were absent from this checkout. This implementation
starts with reproducible measurement records; it does not yet export a binary
resource profile or change SIMPLE runtime defaults.

### Administrator interface

```bash
python3 scripts/hardware_resource_planner.py discover
python3 scripts/hardware_resource_planner.py campaign campaign.json --output output/calibration
# Explicit opt-in: launches preparation, SIMPLE, and scientific validators.
python3 scripts/hardware_resource_planner.py campaign campaign.json --output output/calibration --execute
```

The default campaign command only reads/hashes inputs and prints the proposed
measured-command argv and repetition schedule. It creates no output directory
and launches no process. Discovery reports host logical CPUs, process affinity,
and Linux host RAM; cgroup, scheduler, NUMA, GPU and non-Linux RAM qualification
remain outstanding. An administrator must supply the usable CPU/time/disk
budget explicitly; physical host totals are not allocation grants.

`scripts/hardware_campaign.example.json` is a template for all three adapters.
It requires real executable paths and an administrator-provided
`benchmark-fixture` directory beside the manifest. It is deliberately not
bundled with fake scientific validation. Relative fixture/executable paths are
resolved against the manifest directory. Start with one stratum and a small
input bundle before using the full template's 24 attempts (six candidates,
each with one warmup and three measured repetitions).

### Adapter contract

| Commander | Supported modes | Candidate keys |
| --- | --- | --- |
| `motion_correct` | worker; master rendering | `nthr`; master additionally `nparts`, `ncunits`, `job_memory_per_task` |
| `ctf_estimate` | worker; master rendering | Same as motion |
| `ctf_correct` | shared-memory | `nthr` only |

Scientific settings are explicit per case and remain identical across its
candidates. Patched and global motion, patched/global estimation, phaseflip,
heuristic Wiener and FSC-based Wiener are separate cases. The adapter rejects
resource/execution keys inside `science`, unsupported candidate keys,
oversubscribed CPUs, duplicate candidates, and invalid partition counts.
Worker commands set `part`, `fromp`, `top`, and `nparts=1`; master commands omit
worker selectors so the owning SIMPLE strategy takes the distributed path.
Every argv places `prg=` first, as required by SIMPLE's parser. The first real
campaign exposed this requirement; the adapter regression checks now assert it.

Master execution is intentionally unavailable in this first slice:
`simple_qsys_local` launches workers using `setsid`, so terminating the Python
runner's process group cannot guarantee cancellation of distributed workers.
The real queue job cancellation/accounting contract must be integrated before
enabling master campaigns. Rendering their argv is supported now.

### Fixture preparation and validation protocol

Each case declares a relative project filename, preparation argv arrays, a
validation argv array, scientific parameters, active work count, and candidate
resource dictionaries. The template invokes two fixture-supplied Python hooks:

- `prepare.py <commander> <workdir>` creates the declared project and all
  inputs it can modify/delete beneath the private run directory. Use SIMPLE's
  existing project/import commands; do not parse or patch project binaries.
- `validate.py <commander> <workdir>` checks scientific outputs against
  independently specified truth or a reviewed reference. On stdout it emits
  one JSON object with `passed` (boolean), nonempty `validator_version`, and
  a nonempty `metrics` object. Diagnostics belong on stderr. Wrong outputs
  must report false or exit nonzero.

Hooks are trusted administrator code, not sandboxed plugins. Their executable
paths are resolved and hashed, their argv are recorded, and helper scripts
should live in the checksummed fixture. Scripts referenced elsewhere are the
administrator's responsibility. `{workdir}` and `{fixture}` are literal argv
substitutions; no shell expansion is used. The latter always names the private
fixture copy, not the source bundle.

Every repetition copies the entire fixture, with no symlinks or hardlinks,
and verifies the copied bytes against the planning fingerprint. Hooks must
recreate embedded project paths into that private workspace: copying a
project binary does not rewrite its internal absolute paths. This matters
because estimation can delete its `forctf` input and correction updates
stack/project metadata. Input generation, import and validation are timed as
separate stages, outside the command's measured end-to-end wall time.

The engine enforces the validator protocol, not the scientific truth of an
arbitrary validator. Built-in known-motion, known-CTF and analytic correction
fixtures/validators remain to be implemented and qualified. The existing
motion memory harness remains unchanged; its identical-frame isotropic fixture
is useful for memory/startup measurements but cannot qualify patched motion.

### Lifecycle, budgets and evidence

Warm up each candidate once, then run at least three fresh-process repetitions
in reproducible shuffled blocks. Seed propagation and thread settings are
recorded. Linux children inherit an affinity restricted to the selected logical
CPUs within the process allocation. This is not yet physical-core/NUMA binding;
on other POSIX systems the CPU budget is a planning constraint only.

Each attempt has one deadline covering fixture copying, preparation, measured
execution and validation. The runner checks deadlines between stages and polls
processes every 100 ms, killing their process group on timeout/interruption
or a disk-limit breach. Filesystem copy/hash calls cannot be preempted by this
polling loop, and disk limits can be temporarily exceeded between polls. SIGINT
and SIGTERM trigger cleanup; SIGKILL cannot be handled. Hooks must not detach
children into independent sessions.

Worst-case reserved wall time is attempts times attempt timeout; reserved CPU
time multiplies that bound by the CPU budget. Expansion rejects campaigns
outside either budget. These conservative reservations include preparation,
warmups and failed/interrupted attempts, and persist across resume. Disk checks
count retained artifacts and fixture copies. This slice does not impose a hard
RAM cap; choose bounded fixtures and use site/container enforcement as needed.
It does not interpret a missing memory measurement as zero.

`campaign.json` saves the exact plan, fixture hashes, executable and runner/adapter source hashes,
relevant environment, platform, and CPU binding. Stable experiment IDs derive
from this identity and the candidate/repetition. A process lock prevents two
runners writing the same campaign. `observations.jsonl` appends started/finished
events with stages, argv, exit reasons, wall times, validation and telemetry;
each retry uses a new directory. Resume skips only passed attempts and rejects
changed inputs, executables or environment. A malformed/truncated journal fails
closed rather than discarding evidence. Provision extra reservation budget for
retries; failure does not replenish the budget.

Native `memory_usage_*.csv` files in the run directory yield per-process
high-water RSS records. Completeness is explicitly unverified, node aggregate
and swap are null, and every observation remains `qualification=unqualified`,
even when its hook passes. No production recommendation is exported from these
records yet. CPU/I/O accounting, complete process-tree memory, true allocation
discovery, built-in fixtures, noise-aware ranking, scheduler adapters, and
versioned profile export/Fortran consumption are the next implementation steps.

### Lightweight validation

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=scripts python3 -m unittest benchmark.test_campaign -v
```

These Python checks run tiny Python stand-ins, never SIMPLE or a compiler.
They cover default-off execution, explicit execution/resume, changed-input
rejection, scientific failure propagation, CPU/time budgets, adapter routing,
symlink rejection, retry accounting, timeout and descendant cleanup, and the
render-only master boundary. Actual SIMPLE benchmarks and Linux/macOS/site
qualification remain administrator-run validation.

### First authorized execution, 2026-10-15

An incremental Release build succeeded. The fast gate passed all 15 suites
in 8.43 s after rerunning outside the socket-restricted sandbox (the first
attempt failed the IPC socket test). The new runner then completed four
warmups and twelve measured motion runs, with three repetitions per setting:

| `nthr` | Median command wall time (s) | Maximum process RSS (MiB) |
| --- | --- | --- |
| 1 | 3.9564 | 410.16 |
| 2 | 2.3358 | 410.65 |
| 4 | 1.4355 | 419.80 |
| 8 | 1.0046 | 417.62 |

The fixture was 2048x2048, eight identical frames, 1.3 A/pixel, isotropic
alignment with patches disabled. Runs used eight distinct physical cores in
one socket, fixed seed and explicit OpenMP binding. Every run passed the
declared zero-motion checks: sampled integrated-image correlation >=0.99 and
maximum recovered shift <=0.1 pixels. Observed values were 0.9999999458 and
0.007524 pixels. These checks qualify this smoke fixture only.

The separate existing memory harness also completed twelve successful runs.
Its internal elapsed times are retained separately from the new runner's
command wall times. Eight threads were fastest among the measured settings,
but this does not establish a production optimum or worker concurrency.
CTF estimation/correction, patched-motion truth fixtures and partition sweeps
remain unmeasured in this session. Evidence and machine-readable summaries are
retained locally under `output/resource_benchmark_20261015/`.

## 16. Linux desktop motion profile and launcher

`scripts/benchmark_motion_desktop.py` extends the existing memory fixture helpers
with a bounded concurrent-worker probe. It generates 4096- and 8192-square
float32 movies at 1.3 A/pixel, eight frames each, with known global x shifts
`[0,0,0,0,0,1,2,3]` pixels. It runs `algorithm=patch mcpatch=yes`, records native
RSS, and checks integrated images against the generated truth. Image correlation
must be at least 0.98; profile export additionally requires recovery of global
shifts within 0.1 pixels. These are synthetic global-motion cases, not a
qualification of recovery from local deformation or real detector noise.

```bash
python3 scripts/benchmark_motion_desktop.py \
  --output output/motion_desktop_4k8k \
  --sizes 4096 8192 --candidates 4x1 8x1 2x8 4x8 8x4 \
  --repeats 3 --discard-images --execute
```

Without `--execute`, this command prints its plan. Candidate notation is
threads per worker times concurrent workers. One warmup and three repetitions
are run per candidate, with shuffled order for each repetition. All workers
read the same prepared movie but have separate output directories; inputs are
cache-warm and this is not a measurement of streaming independent movies from
storage. Throughput is completed movies divided by whole-wave wall time.
These worker waves do not measure the SIMPLE master's dispatch/merge overhead.

The probe uses one logical CPU per physical core and leaves at least 20% of
cores outside the benchmark pool on this desktop (32 used, eight reserved).
Workers share that pool with `OMP_PROC_BIND=FALSE`, avoiding pinning every
independent worker to the same first few cores. Each wave has a timeout and
process-group cleanup. Failed image checks stop the sweep. `--discard-images`
removes validated bulky outputs except the first warmup worker per candidate;
original movies/truth, STAR shifts, models, logs, per-wave results and telemetry
remain. Memory is observed, not hard-limited by the Python probe.

When three repetitions complete, export selects the highest median throughput
for batch use, preferring fewer CPU threads and then less aggregate worker
memory among candidates within 5% of the best. Low-latency use selects the
fastest measured single-worker setting. Memory requests are 1.5 times observed
maximum worker RSS, rounded upward to 256 MiB. Export verifies that installed
SIMPLE executables match the benchmark build's ELF code/read-only/initialized
data sections, allowing CMake's build/install RPATH changes; the resulting
profile records full installed executable hashes for later checking.

The primary exported configuration is `desktop_profile.json`, a UTF-8 JSON
file written next to `summary.json`. No executable launcher is generated.
The repository's optional launcher reads this same file:

```bash
python3 scripts/run_motion_desktop.py \
  --profile output/motion_desktop_4k8k/desktop_profile.json \
  --project /path/to/project.simple --movies 100 \
  --size 4096 --frames 8 --smpd 1.3
# Add --execute to run the printed command.
# Add --objective latency to select the measured single-worker option.
```

The JSON contract identifies `profile_type=simple_resource_profile`,
`schema_version=1`, and `commander=motion_correct`. `hardware`, `cpu_ids`, and
`binary_sha256` identify the measured machine, CPU pool and installed binaries.
`workload` gives frame count, pixel size in Angstrom/pixel and scientific mode;
`sizes` keys are square image side lengths in pixels. Each size contains `batch`
and `latency` settings: `nthr`, `ncunits`, `job_memory_per_task_mib`, and measured
throughput/time. Memory is MiB, time is seconds, throughput is movies/second.
`resource_policy` provides machine-readable reserve and partition rules;
`qualification` and `evidence` preserve measurement scope and provenance.
CTF profiles will identify their own commander after calibration; this motion
profile supplies no CTF defaults.

Native SIMPLE consumption is **not implemented yet**. The intended reader must
be opt-in (`auto_apply=false`), reject unsupported versions or mismatched
commander/workload/build/hardware, and preserve explicit job parameters.
It must check current CPU/memory allocation before applying recommendations,
derive `nparts` from the active movie count and clamped worker count, and never
rewrite an existing restart partition layout. The recorded partition policy
is explicitly unmeasured. Missing or incompatible profiles leave existing
SIMPLE defaults intact. JSON replaces the earlier proposed binary handoff;
no Python interpreter will be required by the future native reader.

Dimensions, frame count, sampling and active movie count are declarations by
the operator; the launcher does not parse project binaries to infer them.
Frame counts or sampling outside the measured cases are rejected. It checks
CPU model/kernel identity and executable hashes, clamps concurrency against
current affinity and host available RAM, reserves at least 4 GiB of RAM, and
sets `nparts=min(active_movies,2*ncunits)`. This two-wave partition count is an
unmeasured heuristic. Cgroup/scheduler allocation semantics remain outside this
desktop-only launcher. The fixed benchmark seed is retained in its printed
environment, together with the measured scientific settings.

Profiles remain provisional until representative real movies and the master's
complete launch/merge path are qualified. Increasing frames, adding gain/dose
work, changing sampling/format, or moving to cold or remote storage requires
another measurement. A small movie count may favor the latency profile even
when the batch-throughput profile is faster on a sustained queue.

### Measured result on this desktop, 2026-10-15

Forty waves completed and all 176 worker outputs passed image validation.
Recovered global shifts were within 0.014 pixels of the generated truth.
The Intel Xeon Gold 6242R desktop has 40 physical cores/80 logical CPUs;
the probe used physical CPUs 0-15 and 20-35, with eight cores reserved.

| Movie size | Best tested batch `nthr` / `ncunits` | Median batch movies/s | Worker request (MiB) | Best measured single-worker time at `nthr=8` |
| --- | --- | --- | --- | --- |
| 4096x4096 | 4 / 8 | 0.7900 | 2560 | 5.0685 s |
| 8192x8192 | 4 / 8 | 0.2102 | 9216 | 18.9387 s |

The compared candidates were 4x1, 8x1, 2x8, 4x8 and 8x4. This establishes
the best of those configurations under the eight-worker/32-core desktop
budget, not a global optimum. More workers and alternative partition counts
were not swept. All timings exclude preparation and validation and include
worker process startup and output writes; inputs are synthetic/cache-warm.

The reusable configuration is
`output/motion_desktop_4k8k/desktop_profile.json`, with all candidate measurements
in `summary.json`. A previously exported standalone launcher is retained only
as historical smoke-test evidence.
A real one-movie 4K smoke run through the exported launcher's local master
and worker completed successfully in 7.085 s, and its integrated-image
correlation was 0.99989524. This verifies launcher execution and project merge,
not the performance of the proposed two-wave partition policy. Fifteen Python
checks passed at that milestone. JSON-only export and consumer checks are
validated separately after the output-contract update.

## 17. Desktop CTF estimation calibration

`scripts/benchmark_ctf_desktop.py` runs isolated `ctf_estimate` workers using
4K/8K synthetic micrographs. `scripts/benchmark/generate_ctf_micrograph.py`
filters seeded Gaussian noise with an independently evaluated conventional CTF:
300 kV, 2.7 mm Cs, 0.1 amplitude contrast, 1.3 Angstrom/pixel, isotropic 2-micron
underfocus and zero phase shift. NumPy is required only by the fixture generator;
`--fixture-python` selects its interpreter. The benchmark runner uses standard
Python. The generator records its NumPy version, seed and physical truth.

```bash
python3 scripts/benchmark_ctf_desktop.py \
  --output output/ctf_desktop_4k8k \
  --fixture-python /usr/bin/python3.6 \
  --sizes 4096 8192 --candidates 4x1 8x1 2x8 4x8 8x4 \
  --repeats 3 --discard-images --execute
```

The interpreter path above is specific to this desktop's existing NumPy install.
Omitting `--execute` prints the plan without creating fixtures or running SIMPLE.
Each candidate has a warmup and three measured repetitions. The CPU pool and
process timeout/cleanup follow the motion benchmark. Scientific settings are
fixed: `ctfpatch=yes pspecsz=512 fit_phshift=no dfmin=0.5 dfmax=5 hp=30 lp=5`.
Each output must contain a CTF document with finite defocus, resolution and
correlation; both defocus axes must recover the external truth within 0.1 micron.
This validates global defocus, not a spatially varying patch-defocus surface.

Outputs include `plan.json`, input MRC/truth, per-worker commands/logs/CTF
models/telemetry, `results.jsonl`, `summary.json`, and `desktop_profile.json`.
The JSON uses the resource-profile contract with `commander=ctf_estimate`,
`auto_apply=false`, explicit workload coverage, installed executable hashes,
and `median_micrographs_per_second`. Selection and memory margins follow the
motion profile. Export requires complete successful observations, three measured
repetitions and a single-worker candidate. No launcher is generated.

Recommendations remain provisional: identical cached synthetic input per worker,
no real detector noise or ice variation, no astigmatic/phase-shift sweep, no
master dispatch/merge performance measurement. `nparts` remains an unmeasured
two-wave policy. The existing motion-only launcher rejects this CTF profile;
native SIMPLE consumption is still a future integration. `ctf_correct` is not
run by this benchmark.

### Measured CTF result on this desktop, 2026-10-15

All 40 waves completed successfully (five candidates, two sizes, one warmup
plus three measured repetitions). All 176 outputs passed the global defocus
truth gate; maximum error was 0.001025 micron. The same 32-physical-core pool
as motion correction was used, leaving eight physical cores outside the pool.

| Micrograph | Best tested batch `nthr` / `ncunits` | Median micrographs/s | Worker request (MiB) | Single-worker seconds at `nthr=8` |
| --- | --- | --- | --- | --- |
| 4096x4096 | 4 / 8 | 1.4912 | 1024 | 2.9455 |
| 8192x8192 | 4 / 8 | 0.5248 | 2560 | 8.1588 |

The profile at `output/ctf_desktop_4k8k/desktop_profile.json` was verified against
all observations and the current installed executable hashes. Twenty Python
checks passed; default-off planning created no output directory. No rebuild or
commit was performed. The pilot initially completed SIMPLE successfully but
exposed a validator filename mismatch (`*_ctf.txt`); that was corrected before
the full sweep. The original pilot evidence remains in `output/ctf_desktop_pilot`.

## 18. Particle-to-class-average desktop calibration

The user selected 2,000 synthetic particles, 128-pixel boxes and eight classes.
The current command for de novo class averages is `solve2D` (formerly
`abinitio2D`); its internal `refine2D` stages were formerly called `cluster2D`.
The benchmark uses these current names and introduces no legacy CLI aliases.

`scripts/benchmark_solve2d_desktop.py` adds opt-in shared-memory calibration.
The generic campaign adapter also supports `solve2D` in shared mode, with
explicit scientific settings and no particle-correction output-stack argument.
A distributed solve2D mode is deliberately not enabled before process cleanup
and master/worker accounting are qualified.

The independent NumPy fixture generator writes eight asymmetric analytic
ring-and-lobe shapes with distinct radii at random in-plane angles, with additive
Gaussian noise, a fixed seed,
no shifts and no CTF. Labels and noiseless templates are retained as external
truth. Each run imports the same immutable particle stack into a fresh project;
no candidate inherits another candidate's assignments, classes or cache.

The schedule is fixed across candidates: six stages, `extr_lim=20`,
`nits_per_stage=5`, `eo_stage=no`, `autoscale=no`, `refine=prob_snhc`,
`nsample=2000`, `cache=no`. These are the normal iteration limits, but the synthetic fixture is not a
production convergence qualification. `--extr-lim` and `--nits-per-stage` can
explicitly request a shorter diagnostic schedule. Execution includes final class-average
assembly; fixture preparation and output validation are outside the timer.
The runner applies a per-run timeout, a 20-GiB campaign disk ceiling and the
same 32-core affinity pool used by motion and CTF probes.

Validation requires the expected MRC class-stack geometry, finite nonzero
image data and valid class assignments for all input particles. Adjusted Rand
index and purity compare assignments with the external labels, independent of
class numbering. Profile export requires adjusted Rand index >=0.8 and all eight classes occupied
in every run, in addition to the structural checks. Failure retains measurements
and class averages for inspection but produces no resource recommendation.
Real-data quality remains unqualified; `auto_apply=false` remains mandatory.
The resulting class averages, projects, assignment listings, logs and telemetry
are retained for review.

```bash
python3 scripts/benchmark_solve2d_desktop.py \
  --output output/solve2d_desktop_2000_128_8 \
  --particles 2000 --box 128 --threads 4 8 16 --repeats 3 \
  --fixture-python /usr/bin/python3.6 --execute
```

Without `--execute` the command prints a plan only. The NumPy interpreter path
is specific to this desktop. After one warmup and three repetitions per thread
count, export writes `summary.json` and `desktop_profile.json`. Settings include
`nthr`, `nparts=1`, `ncunits=1` and a memory request with the existing 50% margin.
Here `nparts=1` explicitly means shared-memory execution, not a partition-count
optimum. Selection prefers fewer threads among cases within 5% of the best
median wall time. Results cover the exact particle count, box, class count and
schedule; extrapolation to larger boxes, more classes, sampled updates, CTF or
distributed execution requires another calibration. Native SIMPLE JSON loading
remains a separate future integration.

### Pilot findings, 2026-10-15

Three diagnostic pilots ran, each with one warmup and one measured repetition:

| Fixture / schedule | Threads | Measured seconds | Peak RSS (MiB) | Adjusted Rand index | Purity |
| --- | --- | --- | --- | --- | --- |
| Similar Gaussian-lobe shapes / shortened | 8 | 48.31 | 552.8 | 0.223 | 0.250 |
| Distinct ring-and-lobe shapes / shortened | 8 | 42.79 | 550.4 | 0.380 | 0.454 |
| Distinct ring-and-lobe shapes / normal limits | 16 | 59.04 | 1022.5 | 0.394 | 0.471 |

All runs assigned all 2,000 particles and produced finite class-average stacks,
but none recovered all eight intended classes to the quality threshold. The
normal-limit warmup occupied only four classes (ARI 0.357). These are diagnostic
results, not comparable performance candidates: fixtures, schedules and thread
counts changed between pilots. No full thread sweep or recommended resource
profile was exported. Further fixture/solver investigation or a validated real
particle project is needed before performance calibration is meaningful.

The current helper emits the distinct ring-and-lobe fixture. Earlier generator
hashes and original truth files remain with their respective pilots:
`output/solve2d_desktop_pilot`, `output/solve2d_desktop_ring_pilot`, and
`output/solve2d_desktop_full_pilot`. Each contains `summary.json`, class averages,
projects, assignments and validation metrics. Twenty-six Python checks pass,
including default-off behavior, routing, class-label permutation invariance,
quality-failure export blocking, and qualified-profile selection. No rebuild or
commit was performed; no benchmark processes remain running.

## 19. Consolidated calibration entry point

Contract: `hardware_resource_planner.py calibrate` owns default-off execution
and aggregates the supported motion, CTF-estimation and solve2D profiles into
one JSON collection. Workflow adapters own fixtures, scientific argv, truth
validation and profile qualification. One runner owns process groups, timeouts,
disk limits, affinity, repetitions, warmups, telemetry and evidence persistence.
The three previous desktop CLI names remain thin compatibility entry points;
they retain their per-workflow output-directory semantics. Existing scientific
settings and qualification thresholds must not change during this refactor.

Validation will cover dry-run side effects, compatibility arguments, shared
process cleanup/failure handling, adapter command construction, and aggregated
qualified/unqualified results using fake executables. No SIMPLE build or new
scientific calibration is required for this structural change. Picking,
solve3D, scheduler/cgroup qualification and resume remain outside this slice.

### Implemented interface and validation

```bash
python3 scripts/hardware_resource_planner.py calibrate \
  --workflows motion_correct ctf_estimate solve2D \
  --fixture-python /path/to/python-with-numpy \
  --workdir /scratch/simple-calibration \
  --output simple_resource_profile.json
# Review the plan, then add --execute to run it.
```

`--simple-exec` and `--private-exec` override installation paths. `--sizes`,
`--candidates` (threads x workers), `--threads` (solve2D), `--repeats`,
`--timeout`, `--disk-gib`, `--cpus`, and the solve2D workload/schedule options
control the campaign. Default image candidates are 4x1, 8x1, 2x8, 4x8, 8x4;
solve2D uses 4/8/16 threads. Automatically chosen candidates are filtered to
fit the reserved physical-core pool; explicitly oversized candidates fail.
The legacy motion wrapper keeps its original candidate defaults and discovers
executables from the installed environment.
Explicit subcommands and compatibility runners remain dry-run unless `--execute`
is supplied. The bare administrator entry point runs standard calibration (section 21). NumPy is only
required when executing workflows whose fixture generators need it.

`benchmark/calibration.py` owns the common lifecycle. The registry in
`benchmark/calibration_adapters.py` connects the motion, CTF and solve2D
scientific modules; `benchmark/resources.py` owns the physical-core pool and
`benchmark/profiles.py` owns identity and motion profile export. The optional
motion launcher imports profile helpers from this package. None of these
modules imports the former desktop runner scripts. Those scripts now contain
only compatibility CLI delegation.

The combined file has `profile_type=simple_resource_profile_collection`,
`schema_version=1`, `auto_apply=false`, hardware identity, and a `workflows`
object keyed by commander. Each entry contains `status`, `reason`, `evidence`
and `profile`:

- `provisional`: a qualified synthetic profile is embedded.
- `unqualified`: measurements exist but the repetition/quality/export gates
  were not met; `profile` is null.
- `failed`: fixture generation, execution, validation or export failed;
  `profile` is null and the reason is retained.

The file is atomically replaced after each completed workflow, retaining
previous results if a later workflow fails. Execution failures produce exit
status 2 after the remaining workflows are attempted. Unqualified measurements
are reported in JSON without being promoted to recommendations. Per-workflow
`plan.json`, `results.jsonl`, `summary.json` and eligible `desktop_profile.json`
files remain under the evidence directory. Existing per-workflow schema-1
profiles remain compatible with the optional motion launcher; that launcher
does not consume the collection itself. Native SIMPLE loading is not added.

New output/work directories are required; existing results are never resumed
or overwritten. Process-group cleanup and timeouts cover worker waves and
external preparation/validation commands. The disk ceiling is polled across
the entire combined evidence directory (not a filesystem quota). This remains
a local Linux calibration runner: scheduler/cgroup memory accounting and
multi-node scheduling are not implemented. Scientific settings and thresholds
are retained from the desktop adapters, including the solve2D ARI/occupancy
export gate. No previous measurement is relabelled as a new calibration.

Verified on 2026-10-15: all 36 Python checks pass, including fake-worker end-to-end
collection export, mixed workflow failures, unqualified output, timeout reaping,
disk-limit failure, adapter argv, fresh solve2D projects, and the three legacy
CLIs. A unified dry run planned all three workflows without creating its output
file or evidence directory. `git diff --check` passed. No SIMPLE compilation,
scientific benchmark rerun, or commit was performed for this refactor.

## 20. Picking, extraction and solve3D adapters (original implementation)

The default consolidated calibration now includes all six requested workflows:
`motion_correct`, `ctf_estimate`, `pick`, `extract`, `solve2D`, and `solve3D`.
`pick` is the actual SIMPLE commander name for particle picking. These additions
use the existing common runner; no separate desktop scripts were introduced.
The generic campaign adapter and example manifest also recognize the three
new commanders in shared mode.

```bash
python3 scripts/hardware_resource_planner.py calibrate \
  --fixture-python /usr/bin/python3.6 \
  --workflows motion_correct ctf_estimate pick extract solve2D solve3D \
  --workdir output/calibration_all \
  --output output/simple_resource_profile.json
# Add --execute after reviewing the plan.
```

New adapter scope and validation:

- `pick`: one synthetic micrograph per run, at each requested `--sizes` value,
  with a sparse grid of dark Gaussian particles, known centers, 1.3 A/pixel,
  noise sigma 0.03 and spacing twice `--box`. Uses `picker=new` in Gaussian
  diameter mode, `pick_roi=no`, `pcontrast=black`, fixed `lp=20`, `ndev=2`.
  Output `.box` coordinates are converted from origins to centers and matched
  one-to-one to truth within 0.15 box widths. Precision and recall must both
  reach 0.8; duplicate picks count as false positives. No template-based or
  segmentation-picker optimum is implied.
- `extract`: uses the same synthetic micrograph design and externally generated
  EMAN box coordinates, independently of picking success. Each candidate gets
  a fresh project. Mode is integrated-micrograph extraction, `pcontrast=black`,
  `ctf=no`, `backgr_subtr=no`, `wfloat16=no`, `outside=no`. Validation checks
  particle count, stack geometry and finite pixels, then compares every image
  against its independently cropped, contrast-inverted input. Affine planes
  are projected out before correlation to accommodate extraction's background
  ramp subtraction and normalization. Every crop must correlate >=0.99.
- `solve3D`: generates a known asymmetric Gaussian volume, simulates particles
  at SNR 0.1 with no CTF or translational error. Supplied-volume initialization
  uses `balance=none` on fresh particles, without an unnecessary 2D preparation run. Every candidate starts with a
  fresh project. This first case measures one stage with `vol1` set to the
  known volume, one state, C1 symmetry, `force_lp_range=yes`, `lpstart=15`,
  `lpstop=6`, `filt_mode=none`, and `automsk=no`. It is explicitly a known-start
  screening case, not full de novo reconstruction or a production schedule.
  Validation requires an iteration in the log, a reconstructed state volume,
  finite particle orientations, all requested particles in state 1, and all
  particle update counters positive. The reconstructed volume must correlate
  >=0.8 with smooth external truth in the same frame, resampled to the output
  grid if needed. No rotation/handedness search is hidden in the validator.

All three initially sweep `--threads` in shared memory, not concurrent jobs or
partitions. Export records `nparts=1`, `ncunits=1`, measured wall time and the
existing conservative memory margin. `--sizes` controls picking/extraction
micrograph sizes; `--box` controls their particle box and the solve3D volume/
particle box; `--particles` controls solve3D simulated particle count. The
existing `--candidates` worker-concurrency sweep still applies to motion/CTF.
Combined JSON entries retain the common provisional/unqualified/failed status.
Three measured repetitions and successful quality gates in every run are
required before any new resource recommendation is emitted.

The implementation lives in `benchmark/calibration_spatial.py`; NumPy-only
fixture generation and spatial validators run through `fixture_spatial.py`
using `--fixture-python`. Preparation and validation commands use the common
runner's process groups, timeout and disk budget. Scientific kernels and
Fortran sources are unchanged.

Verified 2026-10-15: 41 Python checks passed, covering registration, command
routing, isolated projects, coordinate conversion, duplicate-pick penalties,
nonfinite rejection, complete-observation export and quality gating. Small
NumPy checks accepted exact extraction crops and the known volume, and
rejected zero crops and an inverted volume. A six-workflow default dry run
created no output artifacts; `git diff --check` passed. No SIMPLE build or
native execution of the new adapters was performed in this change, so their
native end-to-end calibration remains outstanding. Existing motion/CTF results
must not be interpreted as measurements of these new workflows.

### Installed environment discovery

Calibration uses the administrator's loaded SIMPLE environment: executable lookup
uses `PATH`, with `$SIMPLE_PATH/bin` as a fallback. No build directory or prefix
is required, including for the compatibility desktop runners. Optional explicit
executable/prefix overrides remain for development. The inherited `SIMPLE_PATH`
is preserved unless an explicit prefix is supplied. Motion profiles pin the actual
benchmarked executable paths and hashes without assuming a `bin/` layout.
Missing required executables fail before fixture generation; previews remain
available without an installation. Shell aliases/functions are not executables:
load the installation environment so the commands are available on `PATH`.

## 21. Standard calibration and Slurm execution

An installed environment now needs only (from the scripts directory):

```bash
./hardware_resource_planner.py
# Or submit from the cluster login node with its SIMPLE environment loaded:
python3 scripts/hardware_resource_planner.py calibrate --slurm csbdevel --execute
```

Both default to motion correction, CTF estimation, picking, extraction, solve2D
and solve3D; evidence goes to `calibration-results/` and the combined JSON to
`simple_resource_profile.json` in the current directory. Existing destinations
are refused. The bare executable starts standard calibration immediately, as explicitly requested
by the administrator interface contract. Use the explicit `calibrate` subcommand
without `--execute` to preview without submitting or running anything.
All workflow, workload and destination options remain optional overrides. Standard
calibration discards repeated motion/CTF image outputs after validation to stay
within its disk budget, retaining representative warmup images and all metrics;
`--no-discard-images` retains all images.

Slurm submission requests one node, one task, 16 CPUs, 65536 MiB and one hour by
default. Override these with `--slurm-cpus`, `--slurm-memory` (MiB),
`--slurm-time` and optionally `--slurm-account`. The script and evidence paths
must be visible on the compute node. Python 3.9+ runs the planner (the entry point automatically selects an installed
`python3.9`–`python3.14` if the default Python is older); fixture Python
must have NumPy (`--fixture-python` can select a separate interpreter).
The loaded installation environment is exported to the batch job. `srun` binds
the task to allocated cores before discovery; candidates respect task affinity
and `SLURM_CPUS_PER_TASK`, without the desktop 20% core reserve. Slurm enforces
the requested memory and wall-time limits. JSON records allocation identity and
limits in both the collection and each exported profile. This measures within
one node, not distributed scaling, queue wait time or scheduler partition sizing.

Submission prints the job ID and writes `<output>.submission.json` with the
exact submission command/script; stdout/stderr go to `<output>.slurm-<jobid>.log`.
Calibration must be submitted from the cluster, not a desktop with unrelated
paths. Inside an existing single-node allocation, use normal `calibrate` under
`srun` with explicit `--cpus-per-task`; nested submission is refused.

The executable starts with a compact options outline, including the Slurm
submission/resource flags. Completed calibration prints a table of workflow,
size/box, batch/latency/shared mode, `nthr`, workers, `nparts`, requested MiB per
worker and qualification status. Failed/unqualified workflows have no recommended
values. Partition heuristics are explicitly marked unmeasured. The final message
links the created JSON to its intended future use by SIMPLE and explicitly states
that automatic consumption has not been implemented. Submission reports the job
ID, log destination (where the table will appear) and eventual JSON location.

### CSB validation, 2026-10-15

- Job 53893089 submitted to `csbdevel`, 8 CPUs/32768 MiB on `cn132` (Xeon Gold
  6342), synthetic 4K workload, three measured repetitions, shortened 2D schedule.
  Motion and CTF profiles qualified provisionally. Matching 4-thread/1-worker
  medians were 7.5756 s/movie and 4.8619 s/micrograph. Desktop Xeon Gold 6242R
  references were 7.3090 s and 5.1261 s. Installed binaries differ: these are
  environment comparisons, not isolated CPU-performance claims.
- Original picking/extraction calls exposed mandatory `nparts` CLI validation.
  Fixed both adapters with `nparts=1,part=1` while retaining explicit worker
  ranges (no nested scheduler submission). Retry job 53893144 on `cn123` completed
  exit 0 in 58 seconds: extraction qualified, picking correctly withheld its
  profile because the external-truth quality gate failed.
- Bare executable/Slurm preview verified on CSB with its default Python 3.6
  automatically selecting Python 3.12. Fifty local Python checks pass. No build
  was performed. The original job finished exit 2: solve2D failed its accuracy qualification
  and solve3D exposed an invalid final-stage index. The supplied-volume workflow
  begins at stage 4, so its single-stage adapter now uses `nstages=4` (a final
  index, not a stage count); a focused smoke retry follows.
- Local evidence snapshots and desktop comparison are in
  `output/csb_calibration_20261015/`; retry evidence retains its separate node
  identity and must not be presented as measurements from `cn132`.

- Focused job 53893236 exercised the startup outline and final table in a real
  Slurm log; solve3D reported a sigma2 empty-half error after the old 2D prep.
  Removed that preparation in favor of the documented fresh-particle,
  supplied-volume `balance=none` route. Job 53893274 tests that route with
  200 particles, 128-pixel boxes, four CPUs and one repeat (smoke only, not a
  qualified resource profile).

- Job 53893274 also failed in the installed SIMPLE sigma2 bootstrap: all 200
  particles were in the even half and none in the odd half, including without
  the 2D preparation. Thus the earlier preparation was not established as the
  cause. The fresh-particle supplied-volume adapter is retained, but solve3D
  calibration remains unqualified/failed until this initialization issue is
  resolved; no fabricated recommendation is exported. No Fortran changes or
  builds were made. The completed CSB test jobs are no longer running.

## 22. Reduced-cost standard calibration

The observed desktop sweep spent about 24 minutes on motion and seven minutes
on CTF alone. The no-argument standard now retains all six workflows but uses
4K inputs only, 500 particles with 128-pixel boxes, image candidates 4x1/8x1/4x2
and shared-workflow thread candidates 4/8, filtered to available cores. It uses
`extr_lim=4`, `nits_per_stage=1` and a 180-second per-process timeout. Warmup plus
three measured repetitions and all scientific quality gates are retained.
Motion and CTF each require 12 waves instead of 40. This is a smaller synthetic
reference workload, not a qualified prediction for 8K or large particle jobs.

`calibrate --full --execute` restores 4K/8K, 2000 particles, the five image
candidates, 4/8/16 shared threads, 20/5 solve2D scheduling and 600-second timeout.
Explicit workload options override either preset. Slurm forwards the preset and
chooses CPU candidates inside its allocation; compatibility desktop wrappers
retain their previous defaults. Evidence plans and exported JSON identify the
preset. The startup outline describes both choices. Existing running campaigns
and evidence are not modified or restarted by this change.

Validation: 51 Python checks pass, including standard/full selection, explicit
overrides, single-core fallback and Slurm forwarding. No new native benchmark
or build was run; reduced-preset total wall time remains unmeasured.

### Deterministic calibration inputs (initial implementation)

Calibration fixes `SIMPLE_SEED=20261015`, generates native simulation inputs
with one thread, and reuses each fixture across all candidates and repetitions.
The solve3D adapter uses the existing `simulate_particles` commander; its volume,
particle stack, truth orientations and explicit balanced half-set import are
included in the fixture identity. The other adapters retain their fixed analytic
or seeded Python fixtures until their native simulation/truth contracts are
integrated. `simulate_movie` exists but is not yet wired into calibration.

Plans and exported profiles now include a `reproducibility` object containing
the seed, generation thread count, reuse policy and SHA-256 hashes of declared
scientific inputs. Hashes are checked again after each candidate sweep; changed
inputs prevent profile export. Generator source and executable hashes remain
in the existing evidence. This describes reproducibility within the recorded
generator/build environment, not bitwise agreement across compilers or NumPy
versions. SIMPLE's RNG explicitly warns that concurrent OpenMP seed calls are
order-dependent; multi-threaded solver results therefore retain tolerance-based
scientific gates. Timings remain measurements with warmup and repetitions.

Validation: 53 Python checks pass, including input-mutation rejection, JSON
provenance, fixed seed environment, and the single-threaded native simulator
command with balanced import halves. No native build or simulation was executed
for this change. Scientific qualification of solve2D/solve3D is still outstanding.

## 23. Molecular calibration from the simulated workflow

The default calibration now follows the recipe in
`simple_commanders_test_highlevel_workflow.f90`. `--fixture-source analytic`
retains the earlier independent fixtures; compatibility desktop entry points
still default to those fixtures. `--molecule 6vxx` selects bundled spike
coordinates; `--molecule 1jyx` selects bundled beta-galactosidase (the test suite
spells its selector `1jxy`). No network download or test executable is needed.

Calibration is entirely external Python. No SIMPLE source changes, new native
commands, native simulation calls or rebuild are required. The earlier native
calibration changes, including the external-reference initialization patch and
its test registration, have been removed. `workflow_fixture.py` generates volume,
projections, movies and truth using NumPy. `data/*.xyz.gz` contains coordinate-only
copies of the existing example models; `data/molecules.json` records their atom
counts, checksums and source checksum. These assets ship with the Python tool;
runtime does not read the SIMPLE source tree or download molecular structures.
The density model uses equal-weight atom counts with 2 A Gaussian smoothing,
not SIMPLE's native atomic scattering factors. Generator metadata identifies
that distinction and the NumPy version.

The Python adapter generates a molecular volume, 100 spiral projections and
seeded shuffled projection subsets, then eight-frame movies at SNR 0.2. Seed
20260923 and one preparation thread are fixed; successive movies use seed+index.
CTF is isotropic at 2 um with 300 kV, Cs=2.7 mm and amplitude contrast 0.1.
The stationary specimen background is CTF filtered; independent frame noise
sets SNR 0.2. Known integer frame shifts and inverse-shifted noisy averages are
written as external truth, along with exact particle boxes. This is a reproducible
recipe, not a byte-identical replay of the native test simulator. The standard box remains 128; sampling
is `max(1.3,240/box)` A/pixel so the molecular support is not cropped to fit it.
Movies have distinct basenames and at most 10% particle footprint coverage;
the requested particle count is distributed over as many movies as needed.
Actual picking/extraction counts can differ and are recorded in validations
and qualified profiles, separately from the requested count.

Preparation runs the preceding commanders serially, only as far as required
by the selected workflows. It saves each starting project and keeps its data
at absolute paths in `_simulation/<image-size>`. Each timed candidate receives
a fresh project copy; measured outputs never become the next workflow's inputs.
Shared preparation is reused across workflows. All prerequisites, validation
and hashing are outside timing. `--fixture-timeout` bounds each preparation or
validation command separately (900 seconds by default); `--timeout` still
bounds each measured command. Total preparation time has not been measured.
This pipeline can take longer than the independent analytic fixtures.

This first molecular version sweeps `--threads` with one local job and one
partition where supported. It measures complete master commands for image
stages; `--candidates` worker-concurrency sweeps remain specific to analytic
fixtures. Local masters use the installed `SIMPLE_PATH/bin/simple_private_exec`,
whose identity must match discovery. Memory requests conservatively sum process
peak RSS, including child workers. A Linux subreaper supervisor cleans up local
workers even if queue scripts detach into their own sessions. Slurm submits the
same preparation and measurement path inside the granted single-node allocation.

Validation compares motion images with the simulator's optimal average, CTF
with actual defocus, picking with truth boxes, and extraction with crops of the
immutable input micrographs. solve3D starts de novo from the serial solve2D
project; it uses production `dock_volpair` in both hands, then independent NumPy
whole-volume correlation and soft-masked shell FSC against molecular truth
(correlation >=0.8, FSC=0.143 resolution <=40 A). The NumPy FSC discretization is
not claimed numerically identical to the suite's production FSC implementation.
Native output contracts and thresholds still need end-to-end validation.

solve2D checks finite class averages and complete, valid assignments. Continuous
molecular projection directions have no exact eight-class ground truth: the
old synthetic-label ARI gate does not apply. Its timing evidence is retained,
but its recommendation remains unqualified until an independent molecular
classification accuracy criterion is established. This limitation is explicit
in the validation and final JSON; no structural-only success is exported as
an accuracy-qualified recommendation.

Validation so far: 59 Python checks pass, including candidate project isolation,
single execution of serial prerequisites, deterministic subset bytes, rejection
of noise reconstructions, acceptance of identical truth maps, Slurm option
forwarding, and timeout cleanup of a detached worker. Both external coordinate
models, projection stacks, movies and truth outputs reproduce byte-for-byte
in the fixture checks; inverse-shifted frames reproduce the saved optimal image.
`git diff --check` passes and `git diff -- src production` is empty. No native
build, test executable or molecular calibration was run. Required next validation
is an installed SIMPLE smoke run followed by scientific validation and the
three-repetition candidate sweep. Existing measured JSON files were not changed.
When no workflow qualifies, the console labels the JSON a diagnostic report and
omits recommendation legends.

### solve3D timeout follow-up (2026-10-16)

The desktop and CSB comparison reached solve3D stage 3 before the shared
180-second candidate limit terminated it. Preserve the full de novo workflow
and scientific gates; give solve3D its own 1800-second candidate default.
`--solve3d-timeout` overrides that value. An explicit general `--timeout`
still applies to every workflow unless the solve3D override is supplied.
Preview and evidence plans record the effective limit, and Slurm forwards the
option through the existing generic CLI forwarding. Other workflow defaults
remain unchanged. This changes only the external Python calibration.

Validation: 60 Python checks ran successfully (one environment-dependent
check skipped), including preset/override precedence; `git diff --check`
passed. A focused 8-thread desktop retry uses the existing serial solve2D
input under `output/calibration_compare_20261015_external`, with one warm-up
and one measured candidate. Evidence: `output/solve3d_retry_20261016`.
Both full native solve3D runs completed with normal stops on the installed
desktop build cd2229e1: warm-up 1019.417 s (16.99 min), measured repeat
1003.308 s (16.72 min), 8 threads and 927 assigned particles. Peak recorded
memory was 1853.5/1843.3 MiB. The accuracy gate remains unqualified:
truth correlation 0.41465/0.40810 is below 0.8 (FSC143 reported 40 A).
Completion is verified; neither an accurate reconstruction nor optimal thread
count is claimed. The picked-input fixture needs further scientific work.
No rebuild or SIMPLE source change was performed. CSB was not rerun for this
follow-up; Slurm forwarding of the override passed a mocked submission check.
A full multi-candidate sweep takes substantially longer than the previous
timeout-limited campaign; request sufficient Slurm wall time explicitly.

### Extraction removed as a calibration target (2026-10-16)

The current default and selectable workflows are `motion_correct`,
`ctf_estimate`, `pick`, `solve2D`, and `solve3D`. Extraction is no longer a
standalone benchmark, recommendation-table entry, or supported campaign case.
The CLI outline, analytic adapter registry, molecular adapter registry and
campaign example reflect this removal. Earlier six-workflow results above
are historical evidence and are not rewritten.

The molecular preparation chain still runs extraction internally to create
particle stacks for solve2D and solve3D. It is preparation, not a timed
calibration candidate, and does not produce a separate resource recommendation.

### Package organization (2026-10-16)

The public executable remains `scripts/hardware_resource_planner.py`; its
Python-version bootstrap delegates to `resource_planner.cli`. Implementation
now lives under `scripts/resource_planner/`, with `workflows/`, `fixtures/`
(including molecular data), `tests/`, and `examples/`. Hardware discovery is
in `hardware.py`; process limits and cleanup are in `execution.py` and
`supervise.py`. The shared workflow mechanics live in `workflows/common.py`.
Source fingerprints include relative paths throughout the nested package.
The package README is the usage guide; this file remains the design record.

The three desktop compatibility wrappers and old profile-applying launcher
are retired from the source layout and preserved locally under ignored
`output/resource_calibration/retired_tools/`. Their launcher-specific checks
are archived there too; profile-export checks remain in the active suite.
Analytic fixtures remain accessible through the main CLI. Existing
`scripts/memory/` code is unchanged.

New runs default to `output/resource_calibration/standard/profile.json` and
`output/resource_calibration/standard/measurements/`, relative to the launch
directory. Existing directories are still refused; use explicit paths for
additional named runs. Logs remain alongside each attempt's measurements.
The new output subtree is Git-ignored. Historical campaign evidence is not
moved or rewritten because projects can contain absolute input paths.

Validation: all 55 remaining Python checks passed with
`SIMPLE_FIXTURE_PYTHON=/usr/bin/python3.6`, including deterministic fixture
replay, adapter preparation and process cleanup. Discovery and Slurm preview
also passed; no job was submitted. The five removed tests applied only to the
retired motion launcher. `git diff --check` passed. No native build or SIMPLE
runtime tests were performed for this reorganization.

### Direct particle preparation; no extraction (2026-10-16)

Supersedes the earlier internal-extraction design. No extraction commander,
adapter or validator remains in the implementation. Image calibration uses
the motion-to-CTF-to-picking preparation chain. Particle calibration has its
own cached fixture: Python molecular projections are sampled with the fixed
seed into a stack of the requested particle count and imported directly with
`ctf=no` and balanced even/odd assignments. solve2D starts from this project;
solve3D starts from a serial solve2D snapshot. Image preparation failures do
not affect particle preparation. Validation counts assignments against the
direct stack, rather than picking or extraction output. Fixture version is 3;
old extracted-particle timings do not qualify this new clean-projection case.
The dataset remains synthetic and its scientific quality gates are unchanged.
