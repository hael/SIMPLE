# Reference-picker acceleration implementation plan

This plan implements the contract in
[reference_picker_flcf_spec.md](reference_picker_flcf_spec.md). It is a
planning artifact only; no production code has been changed.

## Recommended delivery shape

Deliver in two tiers, in this order.

**Tier 1 — exact structural removals.** A group of small, separately provable,
output-identical changes in `simple_pickref`. They need no summed-area table,
no BLAS, no new module, and no opt-in control. The largest deletes a
per-coarse-position local-standard-deviation calculation whose result no
shipped code path consumes. Measure after this tier; it may change the whole
attribution.

**Tier 2 — one exact batched correlation engine, serving both correlation
loops.** Summed-area local statistics plus bounded dense matrix products, used
by the coarse `match_boximgs` loop and the fine `refine_upscaled` loop. Keep
the exhaustive `13 x 13` fine search and all references. Establish a
scalar-versus-batched shadow mode at each correlation boundary and keep the
scalar result authoritative.

After both tiers are measured in the stream, make a new attribution decision:
if correlation no longer rate-limits, investigate separately specified
reference pruning or search-space reductions only if warranted; if another
stage dominates, stop optimizing the picker and address that stage.

A Fourier replacement of the coarse search is out of scope and is recorded in
[Appendix A](#appendix-a-fourier-coarse-search-not-in-scope) only so the
reasoning is not lost. No phase below depends on it.

Do not combine kernel work with reference-generation, rotation, mask,
threshold, or distance-filter changes.

## Why the coarse loop is in scope from the start

`match_boximgs` (`src/main/pick/simple_pickref.f90:285-338`) and
`refine_upscaled` (`src/main/pick/simple_pickref.f90:668-713`) have identical
inner structure: `window_slim`, then `prenorm4real_corr`, then a serial loop
over references, then `maxval`. One engine serves both with no additional
mathematics and no Fourier conventions to prove.

The split between them is workload-dependent. Coarse work is approximately
`n_coarse_positions * R * N_coarse` and does not depend on particle count.
Fine work is approximately `npeaks * 169 * R * N_fine`. At typical production
geometry, fine exceeds coarse by roughly six times on a dense micrograph but is
comparable to coarse on a sparse one.

Consequently a refinement-only optimization is Amdahl-bounded: a five-times
refinement speedup gives roughly three times overall on dense micrographs and
under two times on sparse ones. Sparse micrographs set stream tail latency.
Pairing the coarse loop exclusively with a deferred FFT rewrite would leave
that bound in place for no benefit, which is the main reason a Fourier coarse
search is not on the critical path.

## Performance model

Let:

- `K` be retained coarse peaks;
- `C` be fine candidates per peak, at most 169 with current settings;
- `N` be pixels per reference box at the relevant resolution;
- `R` be references, normally up to 240 in the stream;
- `P` be coarse positions, approximately `(nx/offset) * (ny/offset)`;
- `G` be the fine grid cell count, approximately `nx * ny` at unit stride.

Current cost is approximately

```text
coarse:     O(P * (window copy + normalization + R*N_coarse))
          + O(P * N_coarse * (2*offset+1)^2)      <- dead avg_loc_sdev
fine:       O(K * C * (window copy + normalization + R*N_fine))
          + O(K * G)                              <- remap
          + O(K * G)                              <- distance_filter scan
```

The two `O(K * G)` terms and the dead local-standard-deviation term are
removed by tier 1 without touching any score.

The exact batched route does not change the leading dot-product count, but
changes its execution:

```text
one-time SAT construction per prepared resolution
+ O(positions) local-statistic queries
+ O(positions*N) bounded window gathering
+ dense matrix products totaling O(positions*R*N)
+ O(positions*R) normalized maximum reduction
```

This removes repeated mean subtraction and division over window pixels,
eliminates general image lifecycle overhead in the inner loop, and presents
the dominant work to optimized cache-aware vector kernels. The expected gain
increases with `R` and batch size. Small cases may use the scalar kernel
through a measured crossover.

## Target architecture

### Shared picker policy

Keep these steps in `simple_pickref`:

1. prepare micrograph resolutions and references;
2. perform coarse scoring and select candidate peaks;
3. invoke a selected fine-refinement score producer;
4. update the fine score grid;
5. apply shared threshold, distance, and output policy.

Separate the scalar correlation at each resolution into a named oracle without
changing its traversal order. Both producers consume the same position list
and fill the same result contract: best score and position plus the equivalent
grid updates.

### Batched real-space correlation engine

Add a small picking-domain type, provisionally in
`src/main/pick/simple_refine_corr_batch.f90`, instantiated once per prepared
resolution. Its responsibilities are:

- prepare double-precision summed-area tables for the prepared micrograph and
  its square, after the global amplitude precondition below;
- flatten prepared references into a consistent `N x R` bank and retain each
  residual reference sum;
- enumerate positions in legacy order, whether strided coarse positions or
  fine candidates;
- query local mean and variance in constant time;
- gather valid raw windows into a bounded `N x Cbatch` workspace;
- evaluate `transpose(S) * W` through the project's linear-algebra boundary;
- apply mean correction and scalar standard-deviation normalization;
- reduce scores in legacy reference/position order;
- expose counters and release owned buffers symmetrically.

Use existing `image` lifecycle and sanctioned access. If the general
`window_slim` API cannot gather valid 2D windows without allocation and
clearing, add a narrow image-layer primitive or iterator rather than exposing
raw representation to commanders.

### Global amplitude precondition

Before building the summed-area tables, subtract the global mean from the
prepared micrograph and rescale it to a bounded amplitude. Both operations are
exactly score-invariant, because Pearson correlation is invariant under a
global shift and a global positive scale.

This is required rather than optional. `mic_shrink` is produced by multiplying
`mic_raw` by `product(ldim_raw)` before the forward transform
(`src/main/pick/simple_pickref.f90:204`), clipping in Fourier space, and
inverse transforming without the matching division; only `mic_raw` is divided
back (`src/main/pick/simple_pickref.f90:220-221`). The subsequent band-pass
uses `hp = 0.` and does not remove DC
(`src/main/pick/simple_pickref.f90:229`). The batched numerator
`G - mu*sum(S)` is then a difference of large single-precision quantities, and
the SPEC's error gates would be at risk for reasons unrelated to the kernel.

Phase 0 records the observed mean and range of `mic_shrink` at both
resolutions so the rescaling constant is chosen from measurement.

### Linear-algebra boundary

The library is already present; only the Fortran interface is new work.

Verified on the current development host: BLAS and LAPACK are required
dependencies (`cmake/Dependencies.cmake:112-148`), BLAS resolves to Apple
Accelerate, and Accelerate is already in the link line of the built
`simple_exec`. Its BLAS exports the classic LP64 `sgemm_` symbol, which is the
ABI `cmake/Dependencies.cmake` documents. Adding a dense product therefore
adds no dependency, no CMake change, and no new link flag.

What does not exist is a boundary to call it through: the repository contains
no `gemm`/`sgemm` call, and every dense product is `matmul` (`src/main/pca/*`,
`src/main/simple_symanalyzer.f90`). Write a narrow module with an explicit
LP64 interface block rather than relying on an implicit interface, so an
ILP64 build fails to link instead of silently corrupting arguments.
Accelerate also exports `sgemm$NEWLAPACK$ILP64` variants, which are opt-in by
macro and must not be selected.

The alternative is `matmul` with a compiler-provided external-BLAS mapping.
That is not enabled in the current release flags
(`cmake/CompilerConfig.cmake:103`), so it would itself be a build change.
Benchmark `sgemm` against plain `matmul` before committing to the wrapper;
gfortran's `matmul` is respectable for these shapes and may make the wrapper
unnecessary.

### Coordinate remapping and grid representation

The fine score grid uses unit stride. Derive a direct mapping from a valid
refined origin to its `inds_offset` element and prove it with edge and
round-trip tests. This removes the current full-grid distance calculation,
`minloc`, array-wide `where`, and the `dists` automatic array
(`src/main/pick/simple_pickref.f90:677,723-733`).

Preserve the old loop-order behavior if multiple refined peaks target the same
grid element. Retain a checked fallback until all supported geometry satisfies
the mapping invariant.

Before cutting the scalar seam, decide whether the fine instance should keep a
full unit-stride grid at all. It currently materializes `positions`,
`inds_offset`, `box_scores`, `box_scores_mem`, and `loc_sdevs` over the whole
micrograph to hold results for the retained peaks only, which is the shared
cause of the remap cost, the `distance_filter` scan, the `dists` stack array,
and most of the picker's fine-resolution memory. The only consumers at fine
resolution are `refine_upscaled`, `distance_filter`, `get_positions`, and
`report_boxfile`. A sparse peak list subsumes several tier-1 items rather than
optimizing them, but it moves the seam, so it is a checkpoint decision and not
an implementation detail.

### Memory and batching

For box area `N`, reference count `R`, and candidate batch `Cb`, account
for at least:

```text
reference bank: N*R reals
candidate bank: N*Cb reals
result block:   R*Cb reals
two SATs:       approximately 2*micrograph pixels in double precision
```

For example, `N=10,000`, `R=240`, and `Cb=169` imply roughly 9.6 MB of
references and 6.8 MB of candidates in single precision before results, SATs,
and thread replication. Do not allocate a maximum batch per OpenMP thread
without a measured memory budget. Start with a configurable internal batch
size and reuse buffers.

Account separately for the existing fine-instance grid allocations and the
reduction from the `nboxes` overcount fix, so kernel memory and structural
memory are not conflated.

### Thread topology

Benchmark two controlled configurations:

1. OpenMP across positions with one BLAS thread per worker;
2. fewer/larger candidate batches with BLAS owning a bounded thread team.

Do not allow both layers to consume the full configured thread count. Set the
vendor thread-count control explicitly rather than assuming the linked
library's default.

That control is vendor-specific and the development and deployment hosts do
not share one. On the current macOS development host BLAS is Apple
Accelerate, which has no `OPENBLAS_NUM_THREADS`; its knob is
`VECLIB_MAXIMUM_THREADS`. On Apple Silicon, Accelerate also dispatches small
and medium dense products to AMX units shared across a core cluster, so
concurrent `sgemm` calls from many OpenMP threads contend for hardware an x86
OpenBLAS build does not have. The development host here has 24 logical and 16
performance cores, which makes that contention easy to hit and easy to
misread as a batching problem.

Consequence for this plan: tune batch size and topology on a host
representative of the deployment target, and label every figure with the host
and the resolved BLAS. Do not carry a batch size chosen on the Mac into the
cluster configuration.

Comparison mode should normally run the backends sequentially. Record the
effective OpenMP and BLAS thread counts, the resolved BLAS, and the host with
every timing.

### Controls and propagation

Provisional advanced keys:

- `refpick_backend=scalar|batched|compare`, covering both correlation loops
  because they share one engine.

It defaults to legacy behavior. Add it through typed parameters, registry,
validation, UI metadata, defaults, commander/strategy, distributed child
commands, and the stream master's explicit `pick_extract` command. Do not
read it ad hoc from `cmdline`.

Tier-1 changes take no control. Each is either the removal of a dead
computation or an indexing identity, and each is separately proved
output-identical.

Add a shadow cadence for live streams and a bounded diagnostic-output policy.

## Work phases

### Phase 0: attribution

Instrument before changing numerical code, in a **Release build**. The only
complete build tree on the development host today is Debug
(`-O0 -g -fcheck=do,mem`), which distorts tight scalar loops and array
intrinsics in opposite directions; no attribution number from it is usable.
Configure and build Release first, and record the build type with every
figure.

Measure, separately at coarse and fine resolution:

- positions entering each correlation loop;
- fine candidate count distribution and clipped edge neighborhoods;
- box area and reference count;
- time in `avg_loc_sdev`, window extraction, normalization/statistics,
  reference dots, maximum reduction, coordinate remapping, `distance_filter`,
  and `setup_iterators`;
- total coarse, refinement, picker, and stream-stage time, with the coarse and
  fine shares reported explicitly so the Amdahl bound is visible;
- observed mean and range of `mic_shrink` at both resolutions;
- thread counts, resolved BLAS, host, build type, and peak RSS, including the
  fine-instance grid allocations.

Include the normal 240-reference stream case, smaller reference banks, and
both a dense and a sparse micrograph. The dense/sparse pair is required, not
optional: it is what determines whether refinement-only work would have been
sufficient.

Gate: the profile confirms the dominant subcosts and records a reproducible
baseline. The reported refinement dominance, and the expectation that
`avg_loc_sdev` is the single largest term, are working hypotheses and not
substitutes for this measurement.

### Phase 1: exact structural removals

Each item is a separate change with its own output-identity proof.

1. Remove or lazily gate the `avg_loc_sdev` call at
   `src/main/pick/simple_pickref.f90:306` and the `loc_sdevs` /
   `loc_sdevs_mem` state it feeds. Its only consumers are `remove_outliers`
   (`447-471`), which `refpick` explicitly disables at line 113, and
   `get_loc_sdevs` (`550-555`), which has no caller in `src/` or
   `production/`. Record a decision on whether those two routines are retained
   as future capability; if retained, gate them so the shipped path does not
   pay for them.
2. Correct the `setup_iterators` first-pass `nboxes` overcount
   (`src/main/pick/simple_pickref.f90:246,248`), which roughly doubles the
   `positions` allocation and the `dists` loop bound.
3. Derive and test direct fine-coordinate remapping; remove the `dists`
   automatic array and the array-wide `where`.
4. Replace the `distance_filter` inner linear search
   (`src/main/pick/simple_pickref.f90:427-441`) with a logical selection array
   indexed through `inds_offset`.
5. Guard the `self%t = minval(scores_refined)` degenerate case
   (`src/main/pick/simple_pickref.f90:723`). If a refined peak leaves
   `box_score = -1.`, the threshold becomes `-1.`, `npeaks` becomes the whole
   grid, and `distance_filter` packs the entire grid into a greedy `O(n^2)`
   loop. The failure mode is a stalled stream stage. Fix it rather than
   preserving it in the parity contract, and update the scalar oracle to
   match.
6. Isolate the existing correlation loops as explicitly named scalar backends
   and define their result contract without altering loop or tie order.
7. Add a specialized valid-window gather that has no per-call allocation,
   whole-output clear, generic dimensionality branch, or impossible outside
   path.
8. Capture scalar score/position summaries on deterministic fixtures.

Gate: items 1-4 and 6-7 are byte-identical in output; item 5 changes behavior
only on an input that triggers the degenerate case, demonstrated on such an
input. The tier is measured and reported on its own, so its contribution is
not later attributed to the kernel.

### Phase 2: exact local statistics

1. Apply the global amplitude precondition to the prepared micrograph.
2. Build summed-area tables of the micrograph and squared values once per
   prepared resolution, using double precision.
3. Query every position's sum and sum of squares with a verified inclusive
   coordinate convention.
4. Reproduce population variance, not sample variance.
5. Define behavior for small negative roundoff, materially negative values,
   zero variance, and near-zero variance.
6. Compare each result against a direct double-precision two-pass oracle over
   random, constant, near-constant, gradient, and large-offset inputs.

Gate: local-statistic tests pass at every valid boundary and over several box
sizes before the tables are used for picking scores.

### Phase 3: exact batched correlation kernel

1. Flatten the normalized references and record their residual sums.
2. Gather positions in legacy order into bounded columns.
3. Calculate dense reference/window dots through the new linear-algebra
   boundary; benchmark `sgemm` against `matmul` before committing.
4. Apply residual-mean correction and divide by `N*sigma`.
5. Reduce each score in the same logical order as the scalar backend,
   explicitly matching strict-greater/tie behavior.
6. Add a scalar crossover for batches too small to amortize packing.
7. Benchmark batch size and the two permitted thread topologies, with vendor
   thread counts set explicitly.

Gate: all exact-kernel gates in the SPEC pass, including adversarial matrix
layout, ties, variance edges, one thread, and the normal worker thread count.
If the error gates fail, investigate the amplitude precondition before
relaxing a tolerance.

### Phase 4: shadow integration at both resolutions

1. Add `batched` as an opt-in producer for the coarse and fine loops.
2. Add `compare`: prepare common inputs, run scalar and batched correlation,
   keep scalar outputs authoritative, and compare per-position results before
   shared downstream policy.
3. Compare best scores, winning coordinates, invalid classifications,
   duplicate-grid updates, thresholds, and final one-to-one matched picks.
4. Write one bounded CSV or JSON-lines record per sampled micrograph. Dense
   score grids are explicit debug artifacts only.
5. Run correctness shadowing in offline replay, then a sparse live cadence.
6. Run separate scalar and batched replays for valid performance results, on
   both a dense and a sparse micrograph.

Gate: picker, throughput, and memory gates in the SPEC pass. Because the
batched route is exact, every pick mismatch is triaged individually; there is
no percentage budget at this tier. Keep scalar as the default until the
requester explicitly approves a change.

### Phase 5: reassess remaining work

Repeat the stage profile after tiers 1 and 2, on both a dense and a sparse
micrograph, and decide whether the picker is still the rate-limiting stream
stage. If it is not, stop here.

If correlation remains dominant, do not hide an algorithmic change inside the
batch implementation. Draft a separate SPEC for one or more opt-in
experiments:

- coarse top-K or score-margin reference retention;
- neighbor-aware reference candidate sets;
- hierarchical or hill-climbing fine positions;
- local quadratic peak fitting.

Gate each against the exhaustive all-reference/all-position oracle.

### Phase 6: stream rollout

1. Run deterministic fixtures and offline stream replay with every micrograph
   shadowed.
2. Categorize mismatches as tie order, layout/indexing, local-statistic edge,
   roundoff, or genuine behavioral change.
3. Run sparse live shadowing with bounded diagnostics.
4. Compare separate jobs for p50/p95 coarse, refinement, and picker latency,
   throughput, queue depth, peak RSS, and downstream particle/class-average
   quality, on the deployment host rather than the development host.
5. Keep legacy kernel defaults for at least the evidence-gathering period and
   record the approved transition separately. Tier-1 removals are not held
   behind this gate, since they are output-identical.

## Test design

### Structural-removal tests

For each tier-1 item, compare `.box` output, pick counts, scores, and the
final score grid before and after on the deterministic fixtures. For the
degenerate-threshold guard, construct an input in which every candidate of at
least one refined peak fails the near-zero-variance test, and assert that the
picker completes with a sane peak count rather than admitting the whole grid.

For the direct coordinate mapping, enumerate edge and round-trip cases and
assert the stride-one and origin invariants, including the fallback path.

### Local-statistic tests

Enumerate every valid window in small random images and compare SAT results
against a double-precision two-pass oracle. Include:

- rectangular micrographs and several even/odd box sizes;
- all four valid boundaries;
- constant and near-constant data;
- a large DC offset with small variance, both with and without the global
  amplitude precondition applied, to demonstrate why the precondition is
  required;
- gradients and isolated impulses.

### Batched score tests

Construct `S` and `W` so expected matrix orientation is unambiguous.
Compare the complete `R x C` batched score matrix against scalar Pearson
scores for:

- `R` and `C` equal to 1 and non-square combinations;
- non-multiple vector lengths and batch tails;
- positive, negative, and tied scores;
- normalized references with a deliberately retained residual sum;
- invalid/near-zero-variance candidates;
- multiple batch sizes and one/multiple OpenMP workers;
- both coarse and fine box areas.

The direct oracle should accumulate in double precision where practical.
Provisional comparison gates are in the SPEC.

### Correlation-contract tests

Exercise `match_boximgs` with controlled micrographs and `refine_upscaled`
with controlled coarse peaks:

- an interior peak with all 169 candidates;
- edge peaks and edge coarse positions with clipped neighborhoods;
- known winning references and positions;
- exact ties in position and reference order;
- two peaks mapping to the same fine-grid cell;
- black/white contrast and noisy distractors;
- one, 20, 120, and 240 references.

Compare every per-position winning score/position and the final score grid
before downstream filtering.

### Picker and stream fixtures

Use Ruben Meana-Paneda's simulated workflow in
`src/main/commanders/test/simple_commanders_test_highlevel.f90`. It already
contains simulated 6VXX/1JYX generation and a retained reference-picking
branch. Replace the source switch during implementation with a runtime test
argument and add coordinate truth or one-to-one matching. Its current
reference stack checks do not establish picking accuracy.

Include one dense and one sparse simulated micrograph, since the coarse/fine
cost split differs materially between them.

`production/tests/simple_test_gencorrs_fft.f90` is a style/timing reference,
not an oracle for Cartesian local Pearson correlation.
`production/tests/simple_test_phasecorr.f90` is not needed by this plan; it
would only matter under Appendix A.

### Performance matrix

| Dimension | Values |
| --- | --- |
| Micrograph | small fixture, typical production, largest supported |
| Particle density | sparse and dense, both required |
| Picking box | small, typical, large |
| References | 1, 20, 120, 240 |
| Refined peaks | low, typical, high |
| Fine candidates | clipped edge, partial, full 169 |
| Threads | 1, normal worker, topology study |
| Tier | baseline, tier 1 only, tier 1 + batched |
| Backend | scalar, batched, compare |
| Host | development (Accelerate) and deployment (cluster BLAS) |
| Build | Release only |
| Scope | kernel, coarse loop, refinement, picker, stream |

Report cold and warm runs separately and repeat/randomize backend order.
Report tier 1 on its own so its contribution is not attributed to the kernel.

## Comparison output schema

Suggested per-micrograph fields:

```text
micrograph_id, dims, box_coarse, box_fine, smpd, nrefs, nthr_omp, nthr_blas,
n_coarse_positions, n_refine_peaks, n_candidates_total, candidate_batch,
scalar_coarse_ms, batched_coarse_ms,
scalar_refine_ms, batched_refine_ms,
locsdev_ms, gather_ms, stats_ms, gemm_ms, reduction_ms, remap_ms,
distfilt_ms, iterator_setup_ms,
mic_mean, mic_absmax,
scalar_peak_rss, batched_peak_rss,
score_max_abs, score_rms,
n_same_winner, n_tied_winner, n_displaced,
displacement_p50, displacement_p95,
n_picks_scalar, n_picks_batched, n_matched,
n_unmatched_scalar, n_unmatched_batched,
result
```

The job summary aggregates counts/quantiles and identifies worst cases.
Comparison timing is diagnostic only.

## Likely failure modes

| Failure | Symptom | Primary test |
| --- | --- | --- |
| Dead statistic silently retained | tier 1 shows no gain | separate `avg_loc_sdev` timer |
| `loc_sdevs` consumer revived later | removed state needed again | recorded capability decision |
| SAT off-by-one | errors only at edges or every window shifted | enumerated boundary windows |
| SAT cancellation | negative/zero variance after large DC offset | offset plus small variance |
| Amplitude precondition skipped | error gates fail on production data but pass on fixtures | large-offset test with and without it |
| Sample vs population variance | box-size-dependent score scale | multiple box sizes |
| Reference residual mean ignored | small systematic score error | deliberately nonzero `sum(S)` |
| Matrix transposed/flattened wrongly | wrong references or shifted candidates | asymmetric `R x C` fixture |
| Batch tail skipped | final positions differ | non-divisible batch sizes |
| Tie reduction reordered | equal maxima choose another cell | exact tie fixture |
| Direct grid mapping off by one | all refined picks displaced | edge and round-trip mapping |
| Duplicate update changed | crowded peaks differ | same-cell fixture |
| Degenerate threshold blowup | picker stalls, `npeaks` equals grid size | all-invalid refined peak fixture |
| `dists` stack overflow | crash on hosts without raised `ulimit -s` | large-geometry run at default stack |
| BLAS oversubscription | high CPU but worse wall time | topology study |
| Vendor thread default assumed | timings differ between hosts | explicit thread-count control in harness |
| Accelerate AMX contention | GEMM slows as OpenMP threads rise, only on macOS | topology study on both hosts |
| Tuned on the wrong host | batch size regresses on the cluster | deployment-host performance run |
| Profiled in a Debug build | attribution ranks the wrong subcost | Release build asserted in the harness |
| Per-thread batch too large | RSS grows with thread count | memory matrix |
| Scalar crossover unstable | discontinuity near cutoff | sweep `R*C*N` |
| Coarse loop left scalar | speedup collapses on sparse micrographs | sparse/dense performance pair |
| Threshold amplification | tiny score drift changes picks | upper-tail and final matching |

## File impact map

Tier-1 production files:

- `src/main/pick/simple_pickref.f90`: `avg_loc_sdev` removal/gating,
  `nboxes` overcount, direct remapping, `distance_filter` selection array,
  degenerate-threshold guard, scalar seams;
- `src/main/image/simple_image*.f90`: only if `avg_loc_sdev` or its
  declaration in `simple_image.f90` is retired.

Tier-2 production files:

- `src/main/pick/simple_pickref.f90`: batched dispatch at both correlation
  loops, compare orchestration;
- new picking-domain batched-correlation module;
- `src/main/image/simple_image*.f90`: only a narrow valid-window/SAT access
  primitive if the picking module cannot use an existing sanctioned API;
- a new narrow linear-algebra wrapper with an explicit LP64 `sgemm` interface;
  none exists today, though the library is already linked;
- `src/main/pick/simple_picker_utils.f90`: pass the backend mode;
- `src/main/params/simple_parameters*.f90`, UI metadata, pick strategy,
  distributed child command, and stream master: typed opt-in propagation;
- `production/tests/` and high-level test commander/UI: kernel, correlation,
  picker, and simulated-workflow validation.

No correlation mathematics belongs in stream modules, UI metadata, parameter
parsing, or commanders.

## Review checkpoints

1. Approve the two-tier SPEC, provisional control name, and numerical parity
   gates.
2. Review Phase 0 attribution, including the dense/sparse split and the
   `avg_loc_sdev` share, and confirm the tier ordering against measurement.
3. Decide whether the fine instance keeps a full unit-stride grid or moves to a
   sparse peak list, before the scalar seam is cut.
4. Review tier-1 output-identity evidence and its measured contribution on its
   own.
5. Review local-statistic and batched-kernel parity before picker integration.
6. Review offline shadow accuracy and separate-run performance, on both dense
   and sparse micrographs, before any live stream trial.
7. Review the post-optimization profile before authorizing any
   behavior-changing pruning or search-space proposal.
8. Review stream throughput, memory, and downstream-quality evidence from the
   deployment host before considering a kernel default change.

## Validation constraints for implementation

Repository policy prohibits agents from compiling, linking, running CMake
builds, or executing test binaries unless the user explicitly requests it.
During implementation, use non-compiling source checks and language-server
diagnostics, then provide the exact build, unit, integration, stream replay,
Linux, and BOX commands for the user. Record only observed results as passed.

## Appendix A: Fourier coarse search, not in scope

Recorded so the reasoning is not lost. No phase, gate, control, or file in
this plan depends on anything below, and no work here is authorized.

Roseman's fast local correlation function would replace the strided
real-space coarse search with FFT-based score maps, sharing local sum and
sum-of-squares maps across references. The formulation and the conditions
under which it matches SIMPLE's complete-square Pearson contract are in
Appendix A of the SPEC.

It is off the critical path for two reasons. It addresses only the coarse
loop, which the tier-2 engine already covers at far lower risk and with an
exact scalar oracle. And it requires proving SIMPLE's transform scale,
template placement, conjugation, lag indexing, and valid-origin conventions
against direct computation before a single score can be trusted, which is a
larger verification burden than the entire batched kernel.

If a future profile shows coarse scoring dominating even after both tiers, the
starting points would be `production/tests/simple_test_phasecorr.f90`,
refactored from a printed demonstration into a self-checking Fourier
convention assertion, and a separate SPEC with its own opt-in control and a
percentage-agreement gate rather than the exact-parity gate used here.
