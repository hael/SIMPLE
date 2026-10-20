# SIMPLE hardware resource planning: open work

Status (2026-10-20): an external Python tool measures a machine and writes recommended settings
to a JSON profile. SIMPLE itself does not read that profile, nothing sets `nparts`, `ncunits`,
`nthr` or `job_memory_per_task` automatically, and streaming is not planned at all. Design,
definitions, formulas and the full implementation record are in
[`../completed/hardware_resource_planning_completed.md`](../completed/hardware_resource_planning_completed.md).
Usage of the tool is in `scripts/resource_planner/README.md`; the calibrated memory estimator is
described in [`../../how2s/memory_estimator.md`](../../how2s/memory_estimator.md).

Terms. `nthr` is OpenMP threads per worker. `ncunits` is the number of workers allowed to run at
once. `nparts` is the number of pieces the work is divided into. `job_memory_per_task` is the
scheduler memory request for one worker. A "profile" is the JSON file the planner writes
(`output/resource_calibration/standard/profile.json` by default). "Qualified" means a measurement
passed its scientific correctness gate; an unqualified workflow gets no recommended values.

The learning so far. The Python side works and runs on a desktop and on the CSB Slurm cluster, but
the recommendations are provisional. They come from identical synthetic inputs, no real detector
noise, no measurement of the master's dispatch and merge cost, and an unmeasured rule for `nparts`
(two waves: `nparts = min(active_movies, 2 * ncunits)`). Two of the five workflows (solve2D and
solve3D) do not yet pass their accuracy gates, so they export no recommendation.

## 1. What exists (do not redo)

- `scripts/memory_estimator.py` with `memory_estimator_models.json`: calibrated memory estimates for
  `motion_correct`, `solve2D`, `solve3D`.
- `scripts/hardware_resource_planner.py` and `scripts/resource_planner/`: discovery, standard
  calibration (motion correction, CTF estimation, picking, solve2D, solve3D) on fixed-seed
  synthetic and molecular fixtures with direct particle import (no extraction), execution limits and
  cleanup, Slurm submission, schema-versioned JSON profiles with `auto_apply=false`, 60 Python
  checks. Dry run by default; `--execute` runs.
- The existing SIMPLE behaviour it measures: `memreport=yes`, `nparts`/`ncunits`/`nthr`, the local
  queue backend warning when `ncunits * nthr` exceeds the processor count. The older Fortran
  `simple_mem_estimator` is not the planning authority.

## 2. Open work

1. **Native consumption of profiles in SIMPLE (the main open item; not implemented).** The reader
   must be opt-in (`auto_apply=false` stays the default), read the JSON without a Python
   interpreter, reject unsupported versions and any mismatch of commander, workload, build or
   hardware, and preserve every explicit job parameter. It must check the current CPU and memory
   allocation before applying, derive `nparts` from the active movie count and the clamped worker
   count, and never rewrite an existing restart partition layout. A missing or incompatible profile
   leaves today's defaults untouched. Needs a ruling on question 6 below (report-only first, or also
   populate parameters). Owner modules: parameters, the `cmdline` normalisation in the commanders,
   the scheduler layer. The `ui -> exec -> commander` flow applies; do not add a parallel
   configuration path.
2. **Streaming planning (not started).** Stage-specific persistent-worker pools, acquisition
   contract (detector, movie size and frames, sustained and peak rate), service rates, backlog and
   drain time, and the report format of section 12 of the completed note (limiting stage, predicted
   peak backlog and drain time, stage pools). Needs detector, rate and burst families in a stream
   calibration corpus.
3. **Real-data qualification.** Profiles are provisional until real movies, real ice and detector
   noise, and the master's complete launch and merge path are measured. Measure `nparts` policy
   (currently unmeasured), per-stage minimum useful part sizes, launch and merge overhead, and parent
   plus worker memory overlap. Frame count, sampling or format outside the measured cases need new
   measurements.
4. **Accuracy gates that fail.**
   - solve2D: continuous molecular projection directions have no exact eight-class truth, so the old
     label-agreement gate does not apply. Needs an independent molecular classification criterion.
     Until then its timing is kept but no recommendation is exported.
   - solve3D: the full de novo workflow completes (about 1000 s on 8 threads, peak memory about
     1.85 GiB, normal stops) but truth correlation is 0.41 against a gate of 0.8 (FSC 0.143 at 40 Å).
     The picked-input fixture needs scientific work. Also on the CSB cluster solve3D calibration was
     failed by an initial sigma2 bootstrap that put all 200 particles in one half; cause not
     established (the installed SIMPLE, not the tool).
   - picking: its recommendation was withheld when the external-truth gate failed on CSB; re-check
     on the current fixtures.
5. **Remaining wiring.** `simulate_movie` exists but is not wired into calibration. Native
   end-to-end validation of the output contracts and thresholds (the Python FSC discretisation is
   not claimed identical to SIMPLE's). `ctf_correct` is not calibrated. CPU binding is within the
   process allocation only: no physical-core, socket or NUMA binding. Cgroup and scheduler allocation
   semantics are outside the desktop launcher.
6. **Gaps before any automatic selection** (the list of section 13 of the completed note, still all
   gating): per-worker models for the main distributed stages; memory models that include
   `ncunits` where concurrency changes the process-tree bound; safe restart when `nparts` changes;
   scheduler semantics validated for SLURM, PBS, LSF, legacy IBM interfaces, SGE, local,
   persistent-worker and coarray execution, with request rendering, lifecycle and verified resource
   semantics reported separately and unknown adapters failing closed; GPU count and device memory;
   physical cores versus logical threads, sockets and NUMA; filesystem capacity and I/O contention;
   validation on more than one machine and SIMPLE build; a versioned multi-dataset calibration
   corpus and reproducible fitting pipeline; independent validation of installation-time site
   profiles and workload-time recommendations; per-node (not aggregate) cluster capacity; one
   experiment and telemetry schema for probes, pilots and campaigns; scheduler-free low-level probes
   that set conservative bounds and the throughput knee; coverage reports; a versioned public
   benchmark bundle with checksummed real inputs; qualification of cluster network, placement and
   launcher effects apart from the portable baseline; extensible hardware descriptions; unknown
   hardware failing closed. Until these hold, the planner recommends explicit command-line values
   and leaves the decision to the user or administrator.

## 3. Acceptance checks

- Before using any plan on a large run (completed note, section 11): keep the estimator JSON; review
  every calibration-range warning; run a pilot with `memreport=yes`; confirm per-process and
  aggregate peaks stay under the allocation with the reserve; confirm `ncunits * nthr` does not
  oversubscribe; recompute `ncunits` and `nparts` from the pilot; for streams, replay sustained and
  peak rates and confirm bounded backlog that drains within the latency target; repeat after any
  change of SIMPLE, compiler, FFT library, allocator, OS or GPU stack.
- The planner fails closed: an unsupported commander or missing memory target gives "no
  recommendation", never a plausible default.
- Native consumption: with `auto_apply=false` a run is bit-identical to today; with it on, explicit
  parameters always win, and an unmatched profile changes nothing.

## 4. Rulings and open questions

Bound: JSON, not a binary profile, is the handoff; no Python interpreter at run time; the Python
tool stays outside the Fortran build; `auto_apply=false` is mandatory until real-data
qualification; extraction is not a calibration target; SIMPLE defaults are unchanged until a ruling.

Open review questions (all still open): the optimisation objective (shortest time, least resources,
throughput, or selectable); the reserve policy for dedicated, shared and persistent-worker nodes;
whether `nparts` may change across restarts; which workflows need one part per socket, node or GPU;
the minimum task duration; whether the first native version only reports or also populates
`nparts`, `ncunits`, `nthr`, `job_memory_per_task`; the datasets of the shipped corpus and which
families are held out; whether a site may train models automatically or an administrator must
promote each; which probes are stable enough for CI; CPU-hour and GPU-hour budgets of adaptive
campaigns; which hardware changes need full requalification; who may promote a provisional hardware
profile; which workflows support true multi-node execution; placement policies to support first;
detector generations and formats for the first stream release; the acquisition safety factor,
maximum backlog and drain target; redistributable benchmark datasets with licences and checksums;
budgets of the reference qualification on a desktop; who maintains the results schema and promotes
model versions; which stream stages have tunable pools; the throughput tolerance defining the
near-optimal set; repetitions and confidence for a conservative site limit; scheduler families and
versions to qualify first; scheduler capabilities mandatory for streaming.
