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
  `motion_correct`, `abinitio2D`, and `abinitio3D`. The calibration data,
  targets, safety factors, and limitations are described in
  [memory.md](memory.md).
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
| `abinitio2D` | whole-commander peak RSS |
| `abinitio3D` | conservative parent plus largest `nparts` worker peaks |

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

The current `abinitio3D` model takes `partitions` and conservatively sums the
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
