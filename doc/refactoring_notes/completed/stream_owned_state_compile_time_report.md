# Stream Owned-State Refactoring: Compile-Time Report

Date: 2026-10-04

Status: implemented in commit 93e855ada, measured on Linux, and validated by the
user-executed fast test gate. The rules derived from this and the follow-up work
are in `doc/policies/compile_time_policy.md`.

## Executive Summary

Changing three heavyweight inline components to owned allocatable state reduced
the measured clean Release build from **225.3 to 138.4 seconds**: **86.9 seconds
saved, a 38.6% reduction**, or approximately **1.63 times the build throughput**.
The compiler optimization level and test inclusion were unchanged across these
source-level experiments. All 14 fast tests passed after each retained step.

The effective pattern was reducing the inline derived-type graph that GCC must
initialize, copy, and finalize,
particularly where tests repeatedly declare an entire stage object. The first
ownership change shortened the dominant dependency chain; subsequent measurements
showed other branches becoming the limiting paths.

This is a measured single-machine result, not a promise of identical improvement
on every compiler or host. Each configuration has one recorded clean-build sample.

## Scope and Exclusions

This report covers only:

- Initial-analysis stage ownership of its particle sieve.
- Sieve-stage ownership of its particle sieve.
- Initial-analysis stage ownership of its parameters.
- The associated lifecycle, pre-start status, and fixture assertions.
- The queue-environment cleanup required before an owned sieve is deallocated.

It excludes changes to compilation/build scripts, compiler flags, warning
suppression, profiler implementation, and Cartesian optimizer/unit-test fixes.
Those changes are not credited with the reductions reported here. The baseline
was established after the preceding build-configuration work, and the same
configuration was used throughout the source-level experiments.

The retained implementation affects five Fortran sources:

| Source | Retained change |
|---|---|
| [src/main/stream/stages/simple_stream_stage_initial_analysis.f90](../../../src/main/stream/stages/simple_stream_stage_initial_analysis.f90#L143) | Owned allocatable sieve and parameters; allocation and teardown ordering |
| [src/main/stream/stages/simple_stream_stage_initial_analysis_tester.f90](../../../src/main/stream/stages/simple_stream_stage_initial_analysis_tester.f90#L1) | Assertions for lazy allocation, cleanup, partial setup, and reinitialization |
| [src/main/stream/stages/simple_stream_stage_sieve.f90](../../../src/main/stream/stages/simple_stream_stage_sieve.f90#L71) | Owned allocatable sieve; allocation-free pre-start status counters |
| [src/main/stream/stages/simple_stream_stage_sieve_tester.f90](../../../src/main/stream/stages/simple_stream_stage_sieve_tester.f90#L1) | Ownership and pre-start status assertions in existing tests |
| [src/main/sieve/simple_ptcl_sieve.f90](../../../src/main/sieve/simple_ptcl_sieve.f90#L261) | Release the owned queue environment during sieve cleanup |

No new test executable or suite was introduced. The numerical algorithms, stream
job settings, publication paths, and GUI metadata record layouts were not changed.
The stage component storage attributes intentionally changed; this requires
recompilation and should not be treated as a binary-compatible object-layout change.

## Measurement Methodology

### Fixed Configuration

- Checkout baseline: `912784aca`, with uncommitted experimental changes.
- Host: `fwl-c142338.ncifcrf.gov`, Linux x86-64.
- Hardware: two Intel Xeon Gold 6242R CPUs; 40 physical cores, 80 hardware threads.
- Toolchain: GNU Fortran/GCC 16.2.0.
- Configuration: Release, retaining `-O3`, with tests enabled.
- Build parallelism: 80 jobs, Make generator.
- Source-level baseline: a successful clean build of 836 compilation jobs.

The existing profiler performed a fresh configure and clean build for each
experiment and recorded monotonic start/end timestamps and elapsed duration for
each compiler invocation. Its wall-time measurement covers the timed Make phase
after configuration, including compilation, archive/link steps, and build targets;
it excludes configuration, the subsequent CTest run, and installation. Every
profile listed below recorded build exit status zero.

**Summed compile elapsed time is not CPU time.** It sums overlapping per-command
wall durations and may include scheduling or I/O waits. It is useful for comparing
work distribution, but must not be interpreted as process CPU accounting or added
directly to build wall time.

An earlier profile with missing timestamps was discarded. It is not evidence for
any result in this report.

### Experimental Sequence

1. Record the first valid clean Release baseline and identify the longest late
   compilation chain using per-file start and finish times.
2. Change one owned component at a time, preserving the workflow and ownership
   semantics and explicitly auditing teardown and callback lifetimes.
3. Extend existing lifecycle tests rather than creating parallel test scaffolding.
4. Run static source/metadata checks, then have the user perform the clean profile
   and runtime gate on the same machine.
5. Retain changes supported by measurements; reassess the limiting branch after
   each step rather than assuming per-file savings add to wall-time savings.

The source implementation was checked without launching builds. Static checks
included editor diagnostics, exact-scope comparisons, Fortran description headers,
test-registry consistency, whitespace, and byte-for-byte generated UI-default parity.
Compiler acceptance and runtime results came from the user-executed builds/tests.

### Running the Compile-Time Profiler

Run these examples from the repository root with the normal project build
dependencies and Python 3 available. Use the same compiler, build type, test
inclusion, and job count for before/after comparisons. Explicit `-j 80` matches
the measurements in this report; choose an appropriate fixed count on another host.

**Clean profiling removes and recreates the build directory. It does not install
the executables or run the test gate.** Saved profile directories are retained.

Record a clean Release baseline before making the next experimental change:

```bash
./scripts/profile_build.sh clean --label ownership-before -j 80
```

After the source change, record the follow-up and run the runtime gate separately:

```bash
./scripts/profile_build.sh clean --label ownership-after -j 80
ctest --test-dir build -L fast --output-on-failure
```

Each run prints its timestamped profile directory. The summary lists whole-build
time and the slowest compilation files; the saved records retain per-file start,
finish, and elapsed times. Compare the two exact directory paths printed by the
runs. For the baseline and final measurements already recorded in this report:

```bash
./scripts/profile_build.sh compare \
  build_profile/release-fixed_20261004_003210 \
  build_profile/initial-analysis-lazy-params_20261004_021155
```

For a new experiment, substitute its printed before/after directories. Negative
per-file B-minus-A deltas mean the second build was faster. The comparison writes
a full delta table into the second profile directory. Its current `compile CPU`
output label refers to summed compiler-command elapsed time, not CPU accounting;
interpret it as described above and use build wall time for the overall result.

To inspect the recorded final profile directly:

```bash
cat build_profile/initial-analysis-lazy-params_20261004_021155/summary.txt
head -n 25 build_profile/initial-analysis-lazy-params_20261004_021155/compile_sorted.tsv
```

A timestamp-only incremental rebuild can be measured in a build tree previously
configured by the profiling script:

```bash
./scripts/profile_build.sh touch \
  src/main/stream/stages/simple_stream_stage_initial_analysis.f90 \
  --label initial-analysis-incremental -j 80
```

The `touch` mode first brings that build tree up to date, then times the forced
rebuild of the specified source and its resulting build work. It is a body-rebuild
cost check, not a replacement for the clean ownership comparisons. For optional
Debug profiling use `clean --debug`; for a build without tests use `clean --no-tests`.
Those configurations are not directly comparable with the Release, tests-enabled
measurements reported here.

## Compile-Time Results

### Whole Build

| Configuration | Compile jobs | Build wall time | Change from preceding retained state | Summed compile elapsed |
|---|---:|---:|---:|---:|
| Original source-level baseline | 836 | 225.3 s | - | 1736.6 s |
| Initial-analysis sieve made allocatable | 836 | 149.3 s | -76.0 s / -33.7% | 1658.2 s |
| Sieve-stage sieve made allocatable | 836 | 148.5 s | -0.8 s / -0.5% | 1630.1 s |
| Initial-analysis parameters made allocatable | 836 | 138.4 s | -10.1 s / -6.8% | 1619.4 s |

The retained source-level changes therefore save **86.9 seconds / 38.6%** against
the valid 225.3-second baseline. The number of compile jobs is unchanged at 836.
The summed compile elapsed reduction is 117.2 seconds / 6.7%, considerably smaller
than the wall-time reduction: shortening the limiting chain improved scheduling
and completion time, not merely the total amount of compiler work.

The 0.8-second sieve-stage wall difference is too small to establish a meaningful
whole-build benefit from one sample. That change is supported by its substantial
local compilation reduction and passing tests, not by a claimed wall-time win.

### Files Directly Affected

All values below are per-file elapsed seconds across the retained ownership changes.

| File | Baseline | Initial-analysis lazy sieve | Both stages lazy sieve | Initial-analysis lazy parameters |
|---|---:|---:|---:|---:|
| Initial-analysis stage | 91.625 | 25.642 | 25.917 | 16.024 |
| Initial-analysis tester | 49.549 | 40.026 | 39.102 | 30.655 |
| Sieve stage | 26.130 | 26.131 | 10.409 | 10.294 |
| Sieve-stage tester | 26.217 | 26.297 | 17.437 | 17.429 |

The initial-analysis stage itself is **82.5% faster to compile** than the original
baseline; its tester is **38.1% faster**. Its sequential stage/tester work fell
from 141.174 to 46.679 seconds. These durations describe dependent compilation
work, not an independently additive wall-time saving.

### Generated Type Layout

Object-symbol inspection gave these default-initialization template sizes:

| Type | Inline baseline | After lazy sieve | After lazy parameters |
|---|---:|---:|---:|
| Initial-analysis stage | 114624 bytes / 111.9 KiB | 73408 bytes / 71.7 KiB | 35104 bytes / 34.3 KiB |
| Sieve stage | 89688 bytes / 87.6 KiB | 48472 bytes / 47.3 KiB | Not changed further |

The initial-analysis template shrank by **69.4%** overall. These are emitted
template sizes, not peak compiler memory or total executable sizes.

Generated copy/finalizer machine-code sizes did not necessarily decrease: they
grew slightly after the first allocation change. The evidence supports reducing
the inline type graph and its initialization/optimization burden, not a simplistic
claim that fewer emitted bytes always mean faster compilation.

## Implementation and Lifecycle Safeguards

### Initial-Analysis Sieve

The sieve is allocated immediately before its first construction, when the first
extraction finishes. Successful combination kills and deallocates it. Stage
teardown releases allocated active or inactive storage, and remains idempotent.
The activation flag, two warm-up cycles, final-ingestion transition, and class
selection/publication order remain unchanged.

### Sieve-Stage Sieve

The sieve is allocated on first start, after mask-diameter validation, and released
during stage teardown. A necessary edge case was initial status reporting:
previously it queried zero counters in the default inline sieve. It now reports
zero accepted/rejected counts while unallocated, without constructing the sieve.
Once allocated, the same getters supply the counters.

### Initial-Analysis Parameters

Parameters are allocated immediately before `parameters%new`. The constructor's
existing complete reset, parsing, and derivation remain authoritative. Teardown
releases parameters after the sieve, asynchronous jobs, queue environment, and GUI
resources. Cleanup after parameter-only setup and subsequent reinitialization are
covered by added assertions.

Constructor registry bindings are local to parsing. The queue constructor reads
parameters without retaining that object. The parameters' UI-program pointer is
a borrowed reference to UI-owned data; deallocating the parameters record does not
destroy that target. Owned strings are released by normal Fortran deallocation.

### Queue Callback Safety

The sieve owns a queue environment that can register a persistent-worker callback
target. Its destructor previously did not release that environment. The retained
cleanup calls `qenv%kill` before the enclosing sieve is deallocated, clearing owned
controllers and any matching callback registration. It does not introduce job
cancellation; already submitted work remains governed by the existing stream policy.

Allocatable ownership was used rather than borrowed pointers. This preserves
intrinsic deep-copy semantics for owned state and avoids aliases to temporary
parameters or constructor-local objects. No OpenMP scratch allocation or numerical
workspace representation was changed.

## Validation and Limitations

| Retained step | User-executed fast gate | Total test time | unit_stream |
|---|---:|---:|---:|
| Initial-analysis lazy sieve | 14/14 passed | 13.33 s | 3.77 s |
| Sieve-stage replication | 14/14 passed | 13.60 s | 4.14 s |
| Initial-analysis lazy parameters | 14/14 passed | 13.95 s | 4.27 s |

Important limits:

- One clean profile per configuration was recorded; cache state, system load,
  frequency scaling, and timing overhead were not statistically controlled.
- The large first ownership improvement is corroborated by per-file timings and
  layout changes, but exact percentages remain host/compiler-specific.
- Passing the fast gate validates its covered cases, not every production stream
  workload. Real persistent-worker reuse and streamed job workflows should receive
  their normal higher-level validation before production adoption.
- No scientific-runtime performance benchmark was performed for these ownership
  changes. Numerical algorithms and their optimization settings were unchanged.
- Module-directory size remained approximately 292 MiB. This was not a dependency
  fan-out or global module-footprint reduction.

## Natural Next Targets

### 1. STAR Exporter's Owned Parameters Snapshot

[src/main/star/simple_starproject_stream.f90](../../../src/main/star/simple_starproject_stream.f90#L17)
embeds a full parameters snapshot. Almost its entire roughly 37.5 KiB default
template is configuration, and the exporter is embedded in both refpick and
preprocessing stages. This is the strongest shared-owner candidate for the next
single-component experiment.

Preserve the owned snapshot rather than borrowing the caller's parameters.
Audit every assignment/read path, including optics-writing paths outside the
main initializer; define release/reuse behavior before changing storage. Measure
the exporter, both stage/tester families, and full build, then run the fast gate
and the relevant STAR-export workflow checks.

### 2. Refpick's Remaining Inline State

[src/main/stream/stages/simple_stream_stage_refpick.f90](../../../src/main/stream/stages/simple_stream_stage_refpick.f90#L85)
still embeds full parameters and two command lines. Each command-line template
is approximately 15.8 KiB, and its tester declares 13 local stage objects.

The latest measured refpick stage/tester durations are 17.562 and 36.829 seconds;
their 54.391-second chain now exceeds initial analysis's 46.679 seconds. Refpick's
tester finishes later. Its owned parameters or an independently initialized job
command line are natural next local targets after the shared exporter experiment.
Preserve pre-initialization queries and release state only after jobs/controllers
that may use it have been torn down.

### 3. Shared Sieve Parameters and Other Stream Parameters

The sieve's own private parameters snapshot is approximately 37.4 KiB and is used
by multiple callers. Indirect ownership could reduce its template globally, but
requires a wider lifetime and reuse audit. Its input can be a temporary local
parameters object: replacing the owned copy with a borrowed pointer would be unsafe.

Optics, preprocessing, pool2D, and solve3D stages also embed parameters and have
repeated fixtures. Treat each as a separate candidate using the new baseline,
rather than performing a blanket conversion. Their urgency depends on their
position in the measured completion graph, not template size alone.

### 4. Repeated Test-Fixture Construction

If ownership changes leave expensive test initialization, consider allocated,
independent fixtures constructed through an out-of-line helper. Keep isolation,
teardown assertions, deterministic inputs, and file cleanup. Do not replace local
fixtures with mutable shared test state merely to avoid generated constructors.

### Selection Rule

For every follow-up: identify a large inline owned component in a measured hotspot,
audit default-state access and all borrowed references/callbacks, change one
component, retain optimization/test inclusion, profile again, and run the runtime
gate. Keep local improvements distinct from critical-path improvements.

The high-level test commander remains a competing late branch at 29.511 seconds.
Reducing refpick alone may move the finish-time bottleneck there. No proportional
future speedup should be extrapolated from the initial-analysis result. Do not
change binary GUI record layouts or hot numerical storage without an independent
correctness and runtime-performance justification.

## Evidence Locations

All profile directories are under the repository's ignored `build_profile/`:

| Profile directory | Purpose |
|---|---|
| `release-fixed_20261004_003210` | Valid source-level baseline |
| `initial-analysis-lazy-sieve_20261004_014458` | First retained ownership change |
| `sieve-stage-lazy-sieve_20261004_015829` | Replication in the sieve stage |
| `initial-analysis-lazy-params_20261004_021155` | Final measured state |

Each directory contains profile metadata, the complete build log, compiler timing
records, sorted timings, the summary, and the successful build exit code. The fast
test results were supplied by the user after the profiled builds.
