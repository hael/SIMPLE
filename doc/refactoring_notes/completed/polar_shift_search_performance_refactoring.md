# Polar shift search: performance refactoring

Status: implementation and local validation completed; integrated in
`69a623cb8` on 2026-11-04. Archived in `completed/` on the same date.
Target-data performance/memory measurements remain follow-up work, with no
production speedup or Linux/BOX validation claimed. The staged/uncommitted
handoff statements below describe the recorded development snapshots.
Initial baseline source reviewed on `master` at
`e13422f3bca0205a808146084de7bf43ef0d6194`. This is the living record for its contract,
implementation sequence, review findings and validation results. The rebase onto
`ac9032343` preserves the newer increment-prior and numerical fixes; section 8
records its validation separately from the original baseline's evidence.

## 1. Objective and scope

Reduce solve2D and solve3D execution time while preserving the recovery gained
from dynamic per-particle re-centring and exhaustive shift/in-plane search.
First remove redundant computation and improve reuse. Changes to numerical
precision, work ordering or search discretization are separate experiments,
accepted only against their own quality and performance evidence.

There is no per-particle cap on the number of exhaustive searches in this plan.
A particle's search history alone does not establish that its next large-domain
search is unnecessary. Never-centred particles, movers and particles with broad
posteriors retain their existing route to the disk.

The shared map kernel serves dense and sparse 2D tables, dense and neighbourhood
3D tables, and classic polar strategies. Cover all these callers when changing
its interface. Keep the candidate map and its bounded alternating shift/in-plane
polish distinct from the final continuous in-plane polish of a committed pose.
The Cartesian `refine=cont` optimizer is outside this refactor.

Related records:

- [Speed-week closeout](../../implementation_notes/completed/single_state_speed_week_closeout.md):
  background for score-loop locality, precision and worker-lifetime experiments.
  Its timings and estimates are historical, not the baseline of this refactor.
- [Polar shift search](../../implementation_notes/completed/polar_shift_search.md)
  and [dynamic centring](../../implementation_notes/completed/polar_shift_search_dynamic_centring.md):
  scientific decisions this refactor preserves.
- [Probability-table memory closeout](probabilistic_alignment_table_peak_memory_validation_closeout.md):
  historical memory/equivalence evidence remains unrecorded. The task was
  archived by maintainer disposition, without inferring those results from
  this polar refactor.

## 2. Contract

For the first, result-preserving phases, retain:

- the master's sampled rows, half labels, reference candidates and their order;
- the disk/window rule, sentinel reset, `shincr`/`shsig` updates, bound inflation,
  and the preparation mask selected from the same particle state;
- the matching bands, radial/directional grids, Bessel cutoff, folded harmonics,
  cell areas, increment prior and its first-centring-only application;
- the tempered greedy cell and the physical-temperature posterior width, with
  existing normalization, frame conversions, cutoffs and tie-breaking;
- candidate ranking, polish bounds, acceptance/fallback rules, final continuous
  polish, assignment, reconstruction and sigma2 ownership;
- random draws and their order within a particle; no new reseeding;
- batching and bounded memory; no all-particle or all-reference full-ring cache.

Target byte identity on fixed-seed one-thread cases for those phases. Verify
kernel outputs separately at multiple threads. Dynamic scheduling currently
assigns thread-local random streams to particles, so a changed multi-thread
trajectory is not by itself evidence of a numerical error or equivalence.

Keep geometry and numerical caches in the existing PFTC objects, and candidate
shortlisting in the search/table objects. Commanders and strategies retain
orchestration. Preserve project ownership, particle-buffer reuse and the
sample-once-and-reproduce contract. Introduce no public setting or persistent
format for the initial work.

## 3. Source findings and owners

Paths in the following table are relative to `src/main/`.
These findings describe the initial baseline. Upstream subsequently replaced
the absolute prior with an increment prior, eliminating particle-dependent
prior tables; the rebase keeps that newer contract.

| Finding in the reviewed source | Owning procedures/files | Consequence to measure |
| --- | --- | --- |
| Production map callers request the greedy cell and discard `logz`; `shiftmap_draw` still sums tempered mass for evidence. Callers also calculate a prior normalizer used by that discarded output. | `pftc/simple_polarft_shiftmap.f90`: `shiftmap_prior`, `shiftmap_draw`; table fills and `strategies/search/simple_shiftmap_srch.f90` | Avoidable normalization and exponential-sum work; physical-width sums remain required. |
| Map setup is reached from each batch fill; setup destroys and rebuilds the maps, workspaces and FFT plans. | `setup_shiftmaps`, `new_shiftmap`; `simple_commanders_prob.f90` batch loops and table `fill_tab_range` routines | Allocation, table construction and planning may be repeated for unchanged geometry. |
| Every pair rebuilds Bessel-weighted shell tables in thread-local `wj`. | `gen_shiftmap`; `shiftmap_work` in `pftc/simple_polarft_calc.f90` | Weights depend on geometry, objective and current half-specific sigma2, not the reference or particle. |
| Every pair reconstructs the full-ring particle spectrum. | `gen_shiftmap`; particle memo/batch lifecycle | Repeated particle work across its reference loop. |
| Prior trigonometry is repeated per pair, and contributing width cells obtain shifts through trigonometric evaluation. | `prior_tables`, `shiftmap_draw`, `polar_shiftmap%get_shift` | Geometry and particle-only terms can be reused. |
| Screened references all receive widths before only the best are retained for polish. | Dense 2D/3D and disk-screen sparse/neighbourhood paths; `shiftmap_srch%srch` | Some width calculations never reach assignment. |
| Shell sums repeatedly traverse product arrays for different radii. | `gen_shiftmap` contraction into `gsum` | Potential cache/SIMD improvement; bandwidth limitation is a hypothesis. |
| Plain Euclidean scoring still multiplies single-precision complex terms by double-precision shell weights. | `pftc/simple_polarft_corr.f90`: `euclid_sums`, `euclid_shell_sum` | Possible vector-width/conversion cost; changing it changes rounding. |

The kernel already uses batched complex-to-real FFTs, single-precision map
arithmetic, contiguous angular slices and FFT-friendly direction counts. Retain
these improvements. The dynamic-centring note records dense-2D window costs of
1.2 ms/particle without the map versus 3.6 ms with it, and 5.1 ms for the disk.
Those older observations motivate measuring window work as well as disk work;
they do not establish the current bottleneck or an achievable speedup.

## 4. Implementation sequence

Each phase is a small, independently reviewable change. Measure before moving
to the next. Phase numbers indicate priority, not a promised delivery schedule.

### Phase 0 — establish the current baseline

Extend the existing benchmark record rather than introduce another reporting
system. Separate map setup, shell construction/contraction, FFT, cell selection,
physical-width calculation and candidate polish. Nested phase timers remain
breakdowns of the existing fill time, not extra terms in the wall-time total.
Aggregate per-thread counters; do not time individual cells.

Report disk and window particle/map counts separately, actual screened/polished
reference counts, grid dimensions and peak process memory. The existing
`screens` counter also counts window maps despite its printed disk-map label;
correct the reporting and update any consumers together. Record both wall time
and thread-seconds without adding them together.

Benchmark an early disk-heavy pass, a settled window-heavy pass, a high-resolution
3D neighbourhood pass and a mixed pass with newly sampled particles. Include
dense/sparse 2D and shared/distributed execution. Use the same input, sample,
references, sigma2, compiler, flags, seed, partitions and threads for each pair.
Repeat matched runs at least three times and record median and spread. Separate
setup/planning from steady-state throughput; report complete iteration time too.

Exit: commands, source revisions, machine, configuration, counts, memory and
logs are recorded here. Attribute time to measured phases before selecting
larger kernel or scheduling changes.

### Phase 1 — separate greedy selection from unused evidence

Provide a PFTC greedy-mode operation that does not request evidence or its prior
normalizer. Share loss/prior/cell-mass arithmetic with the existing evaluator;
do not create divergent implementations of the objective. Migrate every
production map caller, including classic strategies, and audit actual readers
before removing any argument or state.

Preserve the prior and cell measure in cell selection. Preserve the width's
physical temperature, existing intermediate rounding, exponential cutoff and
accumulation order initially. Remove the tempered evidence sum and redundant
prior-normalizer traversal from the production path. Retain the evidence/draw
operation for its existing numerical contract and tests.

Exit: identical greedy cells, losses and widths for representative grids,
including ties; unchanged one-thread workflow outputs; measured phase and
iteration times. Do not claim removal of the physical-width normalization.

### Phase 2 — cache map geometry and half-specific shell tables

Reuse unchanged map tables, workspaces and plans across batch fills. Geometry
validity includes box, angular size, shell band, radial grid, direction count
and harmonic cutoff. Workspace validity also includes the available thread
capacity; plan reuse retains dimensions, strides, alignment and plan mode.
Fixed-seed identity runs retain the existing estimated-plan policy.

Build immutable Bessel-weighted shell tables before the parallel region, once
per map/objective/noise generation and half. Use an existing reliable sigma2
refresh boundary or explicit invalidation; pointer identity does not prove that
sigma2 values are unchanged. Remove redundant per-thread `wj` after migration.

Keep `new`/`kill` symmetric and idempotent. Rebuilding one map must not invalidate
another or retain stale pair state. Size new thread storage for the actual
possible team, and never resize shared caches inside a scoring region.

Exit: unchanged results across repeated setup, changed bands/grids, sigma2
updates, team-size changes, unequal/short final batches, and repeated teardown.
Report allocation/planning counts and memory as well as runtime.

### Phase 3 — reuse particle and posterior geometry work

Prepare the full-ring particle spectrum once per particle/reference loop.
Invalidate it on every new batch and observation/CTF preparation, even when
batch size or particle slot is unchanged. Include half and band dependencies
where relevant. Keep reuse bounded to batch or per-thread storage; do not
duplicate every reference's full spectrum.

Cache direction vectors, radii squared, log cell areas and rotation geometry.
The upstream increment prior has no particle-dependent prior terms, so its
former per-particle tables are removed during the rebase. A separate zero-prior
specialization would need to preserve cell areas and posterior temperatures.

First preserve the existing cell loop and sum order. Treat arithmetic
reassociation or a different route to physical log masses as a numerical
experiment rather than silently including it in this phase.

Exit: geometry/half/frame tests and same-sized batch-remapping tests pass;
unchanged physical widths and centring decisions; fewer preparation calls
and a measured gain without unbounded cache growth.

### Phase 4 — calculate widths only for candidates that can survive

Where screening maps precede a smaller polish shortlist, keep a running list
of its best losses without changing the existing final ranking or tie rules.
Pass its current worst loss as a scalar width cutoff to the map evaluator.
After selecting the greedy cell, the evaluator skips physical moment reduction
only if the cell's ranking loss is strictly above that cutoff. Update the running
losses in the search layer afterward. Candidates later displaced may have incurred
width work; strictly excluded candidates cannot re-enter as the cutoff improves.

Preserve final candidate order and widths. Use one selection traversal per map;
do not repeat the mode solve to request a width, rebuild selected maps, retain
maps for all candidates, or change which references are mapped. Paths already
mapping only the final shortlist need no such change.

Exit: identical shortlist/assignment on ties and adversarial rankings; widths
match an independent direct-grid moment calculation; reduced width-call count
on representative fills. Keep only if the additional bookkeeping pays for itself.

### Phase 5 — improve shell contraction and reference locality

Use Phase 0 evidence to choose the next experiment. First tile angular bins and
radii while preserving each accumulator's shell order. Keep the first Fortran
index contiguous in hot loops and inspect vectorization/temporary-array reports.
Use plain, non-overlapping array interfaces where needed to expose SIMD,
respecting [Fortran argument aliasing rules](https://flang.llvm.org/docs/Aliasing.html);
`contiguous` is a proven caller contract, not a substitute for fixing strides.

Benchmark a tuned matrix-multiplication contraction separately if the shapes
warrant it. BLAS/reassociated sums and single-precision score weights are
numerical changes, even if the mathematical formula is unchanged. Keep the
existing normalization precision until its own evidence supports changing it.
Do not add per-file flags, new subset-library targets or blanket fast-math.

For locality experiments, compare bounded particle/reference tiles or grouping
particles with overlapping neighbourhoods. Preserve each particle's candidate
semantics. Changing work order changes thread-local random-stream assignment;
state that explicitly and use kernel comparisons plus multi-seed quality checks.

On Linux NUMA hosts, compare disjoint process placement and thread layouts,
first-touch owning-thread workspaces, and avoid nested FFT/BLAS oversubscription.
Retain dynamic particle scheduling unless measurements justify another choice.
Record bandwidth/cache counters where available; do not assume more threads help.

Exit: a repeatable improvement beyond timing spread at representative sizes,
acceptable peak memory, and the numerical/quality gates in section 5. Mac
measurements do not establish Linux scaling.

### Phase 6 — optional search-domain/discretization experiments

Only after the preceding measurements, consider an intermediate radius between
window and disk, selected from the existing uncertainty state with immediate
escalation on boundary hits. Retain the full disk for never-centred particles.
The proposed Laplace width is separate scientific work; it is not implemented
and cannot replace the coarse map's multi-modality information by assumption.

Fewer directions on inner rings, or evaluating the origin once, may reduce
redundancy. Preserve the combined central-cell mass and rotation semantics;
different ring grids require new quadrature/transform and recovery validation.
These are discretization changes, not identity-preserving refactors.

Exit: a stated mathematical contract and matched multi-seed recovery evidence
before adoption. A search-count cap or tighter origin-ranked candidate pool is
not an implicit fallback if these experiments fail.

## 5. Validation and acceptance

Follow the [test policy](../../policies/test_environment_policy.md),
[compile-time policy](../../policies/compile_time_policy.md) and
[ownership policy](../../policies/layering_and_project_ownership_policy.md).
Extend existing tester modules and suites; no standalone test executable or
new CTest entry is planned. Expected values come from direct sums, independent
grid evaluations or existing truth fixtures, not the new cache/selector itself.

| Level | Required evidence |
| --- | --- |
| Numerical units | Existing `unit_pftc_align2D3D`, especially `polar_shift_map`, continuous in-plane and probability-table suites. Cover both halves, anisotropic CTF, changed sigma2, folded harmonics, zero/nonzero prior, wraparound/ties, fixed and physical temperatures, disk/window and empty/small candidate sets. |
| Cache/lifecycle | Same-sized batches containing different particles, short last batch, changed mask/CTF preparation, map/band/noise changes, team changes and repeated construction/destruction; no stale pair or particle reuse. |
| Parallel | Compare serial and forced-team deterministic kernels, then repeat at full-machine thread counts (11-24 where available). The capped fast gate alone is insufficient. Preserve per-thread work ownership and avoid planner calls inside parallel regions. |
| Workflow identity | For phases claiming identity: byte-identical one-thread fixed-seed outputs and unchanged decoded candidates, assignments, widths/centring state and random-draw behavior. Diagnose any difference before reclassifying the phase as numerical. |
| Integration | Existing `simulated_workflow_6vxx` and `simulated_workflow_1jxy`; the latter includes `nparts=2`. Exercise fractional sampling, late first visits, movers and classic map-enabled strategies through existing fixtures. If current fixtures lack a needed case, add assertions/cases within the owning suites. |
| Scientific quality | For numerical/order/discretization changes: retain all simulation truth gates; compare pose/shift recovery, class populations/quality, state populations and maps, FSC and difficult-sample reconstruction success across matched seeds. Establish baseline spread and acceptance limits before inspecting the changed results. |
| Performance/memory | Matched phase, iteration and complete-run timings with counts and peak memory. Include beta-gal, PfCRT and a difficult previously failing dataset; add a multi-state case when shared-kernel changes affect it. Kernel throughput alone is not end-to-end acceptance. |

For each code phase, build incrementally in the existing `build/`, run
`scripts/run_fast_gate.sh build` and the affected area/high-level entries, then
install after the gate passes. Use Debug checks for lifecycle/bounds changes.
Run long workflows in the background and retain their logs. Record the source
revision plus a diff identifier for uncommitted changes, exact commands, platform
and results. Do not attribute historical or unobserved Linux runs to this work.

## 6. Deferred work and review decisions

- The speed-week goal of one worker per partition per iteration is now owned by
  the [batch-pipeline plan](../../implementation_notes/completed/prob_align_batch_pipeline.md).
  The earlier waiting-worker fold was dropped; no second worker framework is
  part of this polar refactor.
- The old distributed-master project-read item has substantially advanced:
  assembly and sigma2 commit already operate on the master's builder. Verify
  remaining shared-memory costs rather than reinstate that refactor wholesale.
- Parallel state assembly and GPU offload are outside this plan. Existing
  memory and scientific ownership constraints continue to apply.
- No total speedup is promised. Keep a change only when its relevant acceptance
  gates pass and its benefit exceeds measured variability without moving cost
  or excessive memory into another phase.

## 7. Original implementation and evidence, 2026-11-03

This section records the implementation and measurements before the rebase.
Its fingerprints, workflow identities, work reductions and timings belong to
the original `e13422f3` baseline. They are not validation of the later upstream
scientific changes or of the combined revision; use section 8 for that status.

The retained scope is implemented and is narrower than the first
cache experiment because exact one-thread identity is an acceptance condition.
No exhaustive map, reference candidate or random draw is removed.

| Phase | Current decision | Remaining acceptance |
| --- | --- | --- |
| 0: measurement | Counters distinguish disk/window/full maps, widths, selections, expansions and geometry builds/reuses. Timers surround construction, selection and posterior reduction; setup is nested in fill. The existing parser reads the final records. Construction's contraction/FFT split is deferred to external profiling. | Matched real-data timing, memory and scaling baseline. |
| 1: greedy selection | Implemented one shared evaluator with an evidence-free production mode operation. The evidence/draw API remains available. Both standard fixed-seed workflows reproduce their original outputs. | Representative real-data performance. |
| 2: geometry/weights | Matching geometry, workspaces and FFT plans reuse setup. Debug/Release lifecycle and actual 11-thread checks pass. Shared shell-weight scoring is deferred; retain the original per-thread weights and numerical generator. | Representative memory and setup measurements. |
| 3: reuse | Direction/cell-log geometry is reused. The initial particle-dependent prior tables passed identity checks but became obsolete under upstream's increment prior and are removed in the rebase. Full-ring particle-spectrum caching and zero-prior arithmetic specialization are deferred. | Representative performance. |
| 4: selective widths | A running loss cutoff skips only strictly worse candidates' physical moments. Ties keep widths; final ranking is unchanged; one selection pass per map. Final records show fewer widths with every map retained; the complete multi-state entry passes. | Representative performance. |
| 5: contraction/locality | Deferred until representative measurements justify another experiment. | Real-data performance, memory and numerical acceptance. |
| 6: discretization | Deferred scientific experiment. | Its own mathematical and multi-seed recovery evidence. |

### Configuration and reproducibility

- Baseline: tracked source at `e13422f3bca0205a808146084de7bf43ef0d6194`.
  The original source was archived into `/tmp/simple-polar-release-baseline`
  because Git worktree metadata writes were restricted. No Git refs were changed.
- Host: macOS ARM64; GNU Fortran 16.2; FFTW 3.3.11; OpenMP enabled. Original and
  refactored Release options match, including `-O3 -mtune=native -fPIC` and the
  same system numerics. No fast-math or per-file compiler option was added.
- The existing build initially used `CMAKE_BUILD_TYPE=debug`. It was switched
  to Release for matched comparisons, then restored to Debug, gated and installed.
  A separate full Release build checks the corrected test fixture without
  replacing that Debug build. Installs and child-process paths are
  frozen beneath each run root so rebuilding cannot replace running workers.
- Working results, frozen executables, harnesses and logs are under
  `/tmp/simple-polar-refactor-20261103/`. Long runs use existing test entrypoints
  with their ordinary seeds, thresholds, threads and timeouts. No test entry or
  CTest budget was added.

### Review findings and rejected experiments

The first implementation shared half-specific Bessel/shell weights and retained
the full-ring particle spectrum across references. Its unit tests passed, and
both standard simulated workflows passed their truth gates. Their maps, class
averages, FSCs and all other declared non-project artifacts were identical:
202 artifacts for 1jxy, including its two-part distributed solve3D, and 172 for
6vxx. Decoding all eleven project segments nevertheless found nine/fourteen
`ptcl2D` `shsig` differences. Every other scientific field agreed, including all
`ptcl3D` records. These small width differences were not accepted by loosening
an identity rule or numerical tolerance.

A direct production replay isolates that issue in about eight seconds:

```text
SIMPLE_SEED=20260923 simple_exec prg=solve2D projfile=<scratch>/replay.simple mkdir=yes mskdiam=180 ncls=4 nthr=1
```

Each replay copies the original nightly project, uses the same absolute input
images, and verifies the source project's SHA-256 is unchanged. Three original
replays agree exactly. The first refactored replay differs in seven widths.
Restoring the original interleaved spectrum loop reduces this to one; restoring
private per-pair shell-weight calculation removes that last difference. An exact
comparison verifies the cached weights themselves have identical values, but
copying them into the changed generator reintroduces differences. This supports
a code-generation/arithmetic-layout explanation; it is not a compiler-bug
diagnosis. The retained revision restores the original numerical generator and
its per-thread arrays, with construction profiling in a wrapper.

Other probes ruled out unused evidence calculation, the zero-prior shortcut,
scalar versus cached geometry readout, the combined versus separated posterior
loops, plan reuse, fine timers, and prior-table reuse as isolated causes. The
zero-prior shortcut was also removed to keep the original arithmetic. These
probes are diagnostic evidence, not complete workflow passes. The original
combined/independent posterior reduction orders remain preserved.

The first width shortcut repeated mode selection when a width was needed. It
was replaced with one selection pass and a scalar width cutoff before physical
moment reduction. No map is retained for every reference or rebuilt for a
selected candidate. A rejected width returns `sdev=-1`; a candidate tied at the
cutoff receives its ordinary physical width.

### Final source and local validation

Source diff SHA-256 (`src/main` against the baseline):
`f83c796117d4db047c48c0e8a8fad7949ac16064bfaf45eddfca0f416925abf8`.
Patch, executable digests, commands and run roots are recorded in
`/tmp/simple-polar-refactor-20261103/provenance.json` and the adjacent harnesses.
The original numerical generator body was compared directly with the baseline;
only its private procedure name changes. Construction timing wraps that body.
Original per-thread weights and interleaved spectrum expansion remain present.

Debug bounds checking caught a shape error in the newly added off-origin test
fixture: its phase ramp covered the scoring band, while stored reference and
particle arrays cover a wider band. The corrected ramp uses the full stored
bounds. This is test-only; the high-level workflows do not call that fixture,
and their production patch is unchanged. Its earlier Release assertion results
are superseded. Failure log: `/tmp/simple-polar-retained-debug-fast.log`.

| Final check | Observed result / evidence |
| --- | --- |
| Corrected Debug incremental build and install | Passed, no compiler warnings/errors; `/tmp/simple-polar-corrected-fixture-debug-build.log`, `/tmp/simple-polar-corrected-debug-install.log`. The working build is Debug. |
| Corrected Debug fast gate | 15/15 entries, including 1520/1520 polar-alignment assertions across ten suites; 255.0 s under concurrent compilation/validation, warning-only 60 s budget overrun; `/tmp/simple-polar-corrected-fixture-debug-fast.log`. |
| Corrected Debug full-machine map repeats | Three separate 171/171 passes with actual 11-thread teams; exact serial/team grids, cells, losses and widths; `/tmp/simple-polar-stress-retained-debug-{1,2,3}.log`. |
| Corrected Release full build, area, gate and install | Passed, no compiler warnings/errors; 1520/1520 assertions across ten suites with an actual 11-thread comparison, then 15/15 fast entries in 52.52 s. Logs: `/tmp/simple-polar-corrected-release-build.log`, `/tmp/simple-polar-corrected-fixture-release-{area,fast,install}.log`. Build/install root: `/tmp/simple-polar-retained-release-check`. |
| Standard `simulated_workflow suite=6vxx` | Exit 0; 173 actual artifacts/decoded projects identical to the original; whole-map correlation 0.906779 and FSC 0.143 crossing 4.30345 A, same as original. |
| Standard `simulated_workflow suite=1jxy` | Exit 0, including `nparts=2`; 203 actual artifacts/decoded projects identical to the original. Shared correlation/FSC 0.8990/4.23 A; distributed 0.8968/4.46 A, same as original. Wrong-truth controls pass. |
| Classic map-enabled 2D replay | `refine=snhc_smpl`: every decoded particle field and all 84 maps/stacks, FRC and class-state files identical. The default production replay also matches every particle field. |
| Standard `multistate_workflows nthr=8` | Exit 0; all four routes (class averages, particle search, state refinement and consensus), 129 metrics, pass. Consensus median rotation error 1.67140 degrees (limit 2.5), shift error 0.0822335 pixels (limit 0.15), ARI 0.718615 (floor 0.4), accuracy 0.902667 (floor 0.7). |
| Static/navigation checks | Test registry, FLEX DAG, source descriptors and `git diff --check` pass. Regenerated code overview and Fortran indexes. The five tracked index files also refresh stale pre-existing source maps, so their diff is broader than this change. |

Unit coverage includes both halves; nonzero/zero priors; independent direct-grid
moments for flat, multimodal and concentrated physical posteriors; cutoff ties
and skipped-width sentinel readout; changed sigma2, CTF and observations;
same-size remaps and short batches; changed map geometry and repeated teardown;
independent kth-loss cutoffs; and mixed width choices in a forced team. A
concentrated off-origin case additionally checks that a sub-resolution prior
preserves the represented zero-prior result. That case is not claimed as the
regression oracle for the rejected cache's production identity failure.

Workflow logs are `release-final-{1jxy,6vxx}.log` beneath the scratch root. Existing
identity reports are `runs/identity-final-{1jxy,6vxx}/identity_report.txt`. Their
raw identity exit is 1 only for the known unmatched wildcard described below;
no actual differing or missing artifact remains. Frozen workflow binaries in
`refactored-release-final/` precede only the test-fixture correction, with the
same production patch. Classic replay log:
`/tmp/simple-polar-replay-final-classic.log`. Multi-state log:
`final-release-multistate_workflows.log` beneath the scratch root.

### Work removed and performance limits

Complete final 3D table records give the following work counts. Dense coarse-map
counts and all reported shift-search, objective, gradient and alternation counts
match their original records. Original neighbourhood reports did not count maps;
window-only map counts match their original shift-search counts. No historical
disk-map count is invented for those records. New selections and particle
expansions equal actual maps throughout.

| Completed test path | Maps retained | Width reductions performed | Width reductions avoided | Geometry reuses |
| --- | ---: | ---: | ---: | ---: |
| 1jxy shared (53 fill records) | 59,856 | 38,576 | 35.6% | 26 |
| 1jxy two-part distributed (104 records) | 57,420 | 37,419 | 34.8% | 0 |
| 6vxx shared (85 records) | 159,213 | 57,060 | 64.2% | 90 |

These are reductions in physical moment calculations, not fewer exhaustive
searches or measured complete-run speedups. Shared workers demonstrate geometry
reuse. The distributed fixture records no within-worker geometry reuse; benefits
from different batching or worker lifetimes require separate measurement.
Details: `work-count-comparison.json` beneath the scratch root. The existing
`scripts/plot_refine3d_bench.py` parses all three final workflow paths, including
construction/selection/posterior timers and the new counts, without parser changes.

Three interleaved final Release comparisons were repeated after the multi-state
job exited, with one timing thread and a four-thread cap for the separate
forced-team assertions. All six map-test processes pass. Raw construction
ranges overlap, so no reliable construction speedup is attributed to this
revision. Original/refactored medians and observed ranges, in milliseconds:

| Map case | Original median (range) | Retained median (range) |
| --- | ---: | ---: |
| Window, 8 A | 0.061 (0.049-0.061) | 0.049 (0.048-0.050) |
| Disk, 8 A | 0.560 (0.456-0.567) | 0.456 (0.455-0.465) |
| Window, 2.7 A | 0.559 (0.466-0.560) | 0.466 (0.465-0.466) |

In the final evaluator's eight-reference screen, selection plus physical-moment
readout takes a median 0.088 ms (0.084-0.092) for all widths versus 0.047 ms
(0.045-0.048) with selective widths, across six half/repeat observations. Each
selective case retains all eight maps and measures two widths; each paired time
falls by 44-51%. This is a focused fixture result, not a complete-iteration
speedup. The first full pass also prepares its prior geometry. The benchmark
tester's whole-process RSS is not a comparable allocation baseline because the
new executable runs additional test cases. Results and commands:
`final-kernel-benchmarks.json` in the scratch root.

The completed high-level runs overlapped compilation and other validation;
their elapsed times are correctness-run observations, not matched performance
comparisons. The original pair-run harness also waited for its processes in
sequence, making its elapsed values possible upper bounds for the faster
member. The first final-source timing repeats also overlapped the multi-state
consensus route and show substantial spread; they are preserved in
`final-kernel-benchmarks-under-load.json`. The final repeats above replace them.

The full shared-process peak RSS observations were higher in the final runs:
1jxy original/final 2.468/3.226 GB, 6vxx 2.777/3.006 GB. The distributed worker
peaks were 2.245/2.244 GB. These are process-lifetime high-water marks from
unmatched runs, not map-allocation measurements. A fresh-worker replay of the
same prepared 1jxy references, noise and copied project passes all three pairs,
with byte-identical probability tables. Repeating those pairs after the
multi-state job exits gives original/final median elapsed times of 4.673/4.572 s,
with overlapping ranges 4.554-4.693/4.554-4.735 s. Peak RSS medians are
2.001289/2.001043 GB, with ranges 2.001224-2.001551/2.000945-2.001158 GB.
This window-heavy table has 175 maps and 87 polishes in both builds; the final
build measures 92 widths. It establishes no complete-run speedup. The earlier
under-load pairs are retained in `worker-memory-probe-under-load.json`.
This does not explain or dismiss the shared-process difference; matched
lifetime/memory profiling remains required. The added cache
storage is bounded by map/thread geometry; no particle-spectrum or shared-weight
cache remains. Commands/results: `worker-memory-probe.json` and
`full-workflow-memory-observations.json` in the scratch root.

Beta-gal, PfCRT, a difficult previously failing experimental dataset,
controlled lifetime/peak-memory comparisons and Linux/BOX scaling remain outstanding.
No precision, domain, work-order or scientific-policy change is accepted by
inference from these small fixtures.

### Existing test-environment limitations

- `test=identity target=simulated_workflow_{1jxy,6vxx}` exits 1 because its
  registered per-stage `*.simple` pattern matches neither original nor current
  run. The project now remains at its supplied path. The identity tool's
  ordinary run-prefix and `solve3D_run_id` normalization are retained; no output
  rule or exclusion was relaxed. All 376 actual artifacts/decoded projects match.
- `test=solve3D_addon nthr=8` fails after its base solve3D because the fixture
  reads `strict/1_solve3D/addon_gate.simple`, while solve3D keeps the project at
  its supplied path. The original-source Release executable fails at the same
  path. Logs: `original-release-solve3D_addon.log` and
  `final-release-solve3D_addon.log` in the scratch root. No addon pass is claimed.
- Original Debug baselines and superseded Release/multi-state attempts were
  intentionally interrupted for implementation changes. Their logs/products
  remain in scratch; they are not passing integration evidence. Short bisection
  workflows were deliberately terminated after decoding 2D records (exit -15).

### Remaining acceptance and handoff

At the pre-rebase handoff, the retained source, tests and this record were ready
for maintainer review. Debug/Release gates, complete simulated workflows,
multi-state truth gates, final paired probes, source fingerprint and
branch-context closeout were complete. All edits were uncommitted; no push was
made. The working build was restored to Debug and installed after its passing gate.

The retained scope can be reviewed using local numerical and integration
evidence. Deployment-scale speed, memory and recovery acceptance still needs
matched experimental datasets and the target Linux/BOX hardware. Deferred
shell-weight/spectrum caches, contraction/precision/order changes and search
discretization require separate measured justification; they are not silently
included in this implementation.

## 8. Rebase onto upstream `ac9032343`, 2026-11-04

The maintainer's `git pull --rebase` stopped while replaying `88dcbdc41`.
Resolved the eight conflict files by retaining upstream's increment-prior API,
harmonic-index wrapping, angular Nyquist corrections, insertion-loss helpers
and eligible-set normalization. The performance mode, strict width cutoff,
geometry/plan reuse and existing BENCH extensions remain.

Removed the obsolete per-particle prior arrays, readiness state, helper and
unused direction getter. All production and test mode/draw/prior callers use
the increment-prior API. The independent physical-moment oracle now applies
the prior to the cell's increment. Upstream's added tests remain present.
The numerical generator body matches `ac9032343` exactly after its private
procedure rename, including wrapped harmonic indices.

Verified combined-source checks:

- Incremental Debug build and installation pass without compiler warnings/errors.
- Complete alignment area: 1588/1588 assertions in 11 groups, including upstream's
  Nyquist, wrapped-order, insertion-weight and stage-omission coverage; actual
  forced map team 11.
- Debug fast gate: 15/15 in 106.2 s under concurrent Release compilation; the
  60 s budget overrun is warning-only.
- Three further Debug map runs: 178/178 each, actual teams 11, identical serial
  and parallel grids, modes, losses and physical widths.
- Registry, FLEX DAG, descriptors, conflict-marker scan, document links and
  `git diff HEAD --check` pass. Navigation indexes regenerated for merged source.
- Release incremental build and installation pass without compiler warnings/errors;
  complete alignment area 1588/1588, actual map team 11; fast gate 15/15 in 30.08 s.

Both standard simulated workflows exit 0, including final reconstruction and
wrong-truth controls. Reported map correlation / FSC 0.143: 6vxx 0.9018 / 4.54 A;
1jxy shared 0.898921 / 4.30345 A, two-part distributed 0.899477 / 4.30345 A.
The 1jxy suite asserts both distributed parts ran. The complete multi-state
entry also exits 0: all four routes and all 129 metrics pass. Median errors:

| Route | Rotation degrees (limit) | Shift pixels (limit) |
| --- | ---: | ---: |
| Class averages | 2.14521 (4.5) | 0.110326 (0.4) |
| Particle search | 1.63267 (2.5) | 0.0698007 (0.15) |
| State refinement | 1.22509 (2) | 0.0467114 (0.1) |
| Consensus | 1.65965 (2.5) | 0.0702912 (0.15) |

All integration checks use the new Release installation. State refinement
recovers the planted populations exactly. No pre-rebase workflow identity or
timing is claimed for this source, and no Linux/BOX or experimental-data run
was performed. The prior identity/addon harness limitations remain historical
findings; those entries were not repeated for this rebase.

Commands, logs and initial three-way conflict snapshots are under
`/tmp/simple-polar-rebase-20261104/`. Resolved source patch SHA-256 against
`ac9032343`: `4c03753aeecf0bf72f3688f2ee38c445cf872a43f7e50d052b3b3e729795a031`.

At this validation handoff, the assistant had performed no rebase continuation,
commit or push; those actions were reserved for the maintainer. All resolved
source and navigation files were staged, with zero unmerged index paths. This
record, untracked when the original commit was made, was included in the handoff.
The working build was Debug and installed, with no remaining local validation
runs. The maintainer subsequently completed integration as `69a623cb8`.
