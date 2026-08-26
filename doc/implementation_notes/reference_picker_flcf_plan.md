# Reference-picker acceleration note

Status: optimized source branch implemented; 6VXX compare-mode numerical
validation observed; release benchmarking, optimized-only, and real-data
validation remain outstanding.

This is the single design and planning note for reference-picker acceleration.
It consolidates the former specification and plan. Git retains the earlier
documents if their full history is ever needed.

## Reviewed direction

Implement one new optimized reference-picker code path beside the untouched
legacy path. The legacy path remains the default and the validation oracle
while the optimized path is developed. Once the optimized implementation and
simulated tests are complete, compare both paths on real data and decide
whether the optimized path is ready to become the default.

## Implementation state

The source implementation now includes:

- one typed `refpick_backend=legacy|optimized|compare` key, defaulting to
  `legacy`, with validation and UI exposure for direct picking,
  reference-picker validation, simulated workflow validation, and streaming;
- propagation through picker construction, distributed job command lines, and
  the stream master/reference-picking child command;
- a narrow LP64 single-precision `transpose(A) * B` wrapper over the CPU SGEMM
  already linked by SIMPLE;
- `simple_pickref_corr_batch`, which prepares a globally centered/scaled
  micrograph copy, double-precision summed-area tables, flattened prepared
  references, bounded candidate/GEMM buffers, residual-reference-sum
  correction, and per-batch maximum reduction;
- optimized coarse and fine paths that do not allocate or calculate
  `loc_sdevs`, plus direct fine-grid remapping, invalid-neighborhood handling,
  corrected iterator allocation, and indexed distance-filter updates;
- `compare` diagnostics for coarse and fine score grids while retaining legacy
  output as authoritative;
- compare-mode wall timers around coarse/fine correlation and distance
  filtering, including per-stage and combined divergent-stage speedups; compare
  mode now also runs and checks the optimized distance filter before restoring
  legacy output;
- runtime `picker` and `refpick_backend` arguments in the existing simulated
  workflow, removing its compile-time picker switch;
- optimized-path candidate/reference/batch/workspace reporting and separate
  preparation, gather/statistics, SGEMM, and reduction timings.

The user-observed 6VXX compare run covered ten micrographs and 100 references.
All coarse and fine comparisons reported zero scores beyond the `2e-4`
tolerance. Worst coarse max/RMS errors were `1.788e-6`/`8.639e-8`; worst fine
max/RMS errors were `4.292e-6`/`5.222e-7`. The run completed the picking phase
normally. `legacy` remains the production default until optimized-only,
release-performance, second-system, and real-data validation are complete.

## Current optimized profile and next targets

The first 6VXX compare log was not a controlled release benchmark, so these
numbers are directional. Across ten micrographs, the instrumented optimized
correlation phases totaled about `7.03 s`: `3.40 s` SGEMM (48%), `2.60 s`
window gathering/statistics (37%), `0.78 s` score reduction (11%), and `0.25 s`
preparation (4%). Fine correlation accounted for about `6.06 s` (86%) of that
total; coarse correlation accounted for about `0.98 s`.

The release benchmark should determine optimization order:

1. Fine SGEMM is the largest measured component. Test bounded batch sizes and
   record BLAS/OpenMP thread settings before changing the algorithm.
2. Window gathering and summed-area statistics are the largest SIMPLE-owned
   kernel. If still material in release mode, tile or parallelize candidate
   gathering under a single explicit thread topology and avoid redundant
   candidate-buffer writes.
3. Score reduction is especially visible in the coarse pass. Compact invalid
   references once and use a vectorizable maximum reduction if profiling
   confirms the branch/reduction cost.
4. The new stage timers measure candidate construction, remapping, and both
   distance filters outside the correlation-engine phase counters. Optimize
   those only if their release timings remain significant.

Do not introduce search-space reductions until these exact-equivalence kernels
have been measured on the deployment BLAS.

The source review established these corrections and decisions:

1. SIMPLE already links CPU BLAS and LAPACK through
   `cmake/Dependencies.cmake` using the documented LP64 Fortran ABI. The picker
   only lacks a narrow CPU GEMM interface. Adding an explicit `sgemm` interface
   is manageable and does not require a new numerical dependency.
2. The optimized path must never calculate `avg_loc_sdev` at every coarse
   position. No current repository caller uses the public statistic and the
   internal outlier-removal call is disabled. If a concrete caller later needs
   it, calculate it once for the final selected boxes; otherwise the optimized
   path carries no local-standard-deviation state. The legacy path remains
   unchanged for reference and compatibility.
3. The `setup_iterators` first pass over-counts `nboxes` and over-allocates
   `positions`; it resets `nboxes` before fine refinement, so it does not
   double the later `dists` loop bound.
4. An all-invalid fine refinement can set the threshold to `-1`, admit the
   entire score grid, and pass zero entries from `inds_offset` to
   `positions(0,:)`. The optimized path must exclude invalid refinements.
5. The summed-area/BLAS formulation is mathematically equivalent to the scalar
   Pearson score but cannot be assumed bitwise identical because its operation
   order differs. Validation therefore compares numerical tolerances, final
   coordinates, and scores.
6. A sparse replacement for the full fine score grid is a separate broad
   redesign. It is not required for this implementation and is deferred.

## Goal and scope

Reduce reference-picking latency without reducing the current search space or
changing the scientific score. The production picker has two structurally
similar correlation loops:

- coarse search at 4 A/pixel with a three-pixel stride in `match_boximgs`;
- fine search at 2 A/pixel with unit stride in `refine_upscaled`.

With the current sampling constants, each retained coarse peak produces at
most a `13 x 13` fine neighborhood. Stream picking commonly expands ten
medoids into rotations and mirrors, producing as many as 240 references. The
optimized correlation engine must serve both resolutions.

The optimized path includes:

- no per-coarse-position `avg_loc_sdev` work;
- corrected iterator allocation, fine-coordinate mapping, filtering, and
  invalid-refinement handling;
- allocation-free gathering of known-valid 2D windows;
- summed-area local statistics;
- bounded BLAS matrix products across references and candidate windows;
- the same exhaustive candidates, references, score definition, and final
  selection policy as the legacy path;
- built-in comparison counters and timings suitable for simulated and real
  data validation.

The following are out of scope: reducing references or rotations, top-K
reference retention, shrinking the fine neighborhood, hill climbing,
threshold changes, a sparse fine-grid redesign, and a Fourier coarse search.

## Verified legacy behavior

Both `match_boximgs` and `refine_upscaled` currently perform this sequence for
every candidate position:

1. `window_slim` extracts a complete square, allocating small coordinate
   arrays and clearing the output image on every call;
2. `prenorm4real_corr` subtracts the window mean, computes population
   variance, and divides by its standard deviation in single precision;
3. a serial loop scores the normalized window against every normalized
   reference;
4. `maxval` retains the best score.

The complete square is part of the score contract. Zero-valued background in
an automasked reference still participates after complete-box mean
subtraction. For reference `S_r`, candidate window `T_c`, and `N` box pixels:

```text
score(r,c) = sum((S_r - mean(S_r)) * (T_c - mean(T_c)))
             / sqrt(sum((S_r - mean(S_r))^2)
                    * sum((T_c - mean(T_c))^2))
```

Other avoidable work in the legacy implementation is:

- `avg_loc_sdev` at every coarse position even though its internal consumer is
  disabled;
- an unused `loc_sdevs_mem` copy;
- oversized `positions` allocation from the first-pass `nboxes` over-count;
- a full `dists(self%nboxes)` calculation, `minloc`, and whole-grid `where`
  scan for every refined peak;
- a whole-grid `distance_filter` update with a linear peak search at every
  cell.

The legacy path remains unchanged. These observations define what the new
optimized path avoids rather than a sequence of edits to the oracle.

## Optimized correlation formulation

Let `S` be an `N x R` matrix of prepared reference columns and `W` an
`N x C` matrix of raw candidate-window columns. Obtain each candidate's sum
and sum of squares from double-precision summed-area tables:

```text
mu(c)       = sum(T_c) / N
var(c)      = sum(T_c*T_c) / N - mu(c)^2
sigma(c)    = sqrt(var(c))
G           = transpose(S) * W
score(r,c)  = (G(r,c) - mu(c) * sum(S(:,r))) / (N * sigma(c))
```

The residual-reference-sum correction is required because normalized
single-precision references need not sum to exactly zero. Candidate columns
must use the same pixel order as reference columns, and maximum reduction must
preserve legacy candidate/reference traversal and tie semantics.

Use the CPU BLAS already linked into SIMPLE through a narrow explicit LP64
`sgemm` interface. Keep BLAS details behind a small linear-algebra boundary;
the picker owns batching and score policy, not ABI declarations. No new CMake
dependency or link flag is expected. Verify the resolved BLAS on macOS and the
deployment host, and prevent nested OpenMP/BLAS oversubscription.

The prepared micrograph can carry a large DC value and amplitude after Fourier
resizing. Build the optimized tables and candidate bank from a temporary copy
with its global mean removed and a bounded positive scale. Pearson correlation
is invariant to those operations in exact arithmetic. The legacy oracle stays
on its original data, so comparisons still use numerical tolerances.

One bounded engine instance serves one prepared resolution and owns:

- double-precision summed-area tables for the micrograph and its square;
- the flattened normalized reference bank and residual reference sums;
- reusable candidate and result workspaces;
- best-score and best-position reductions;
- work and timing counters plus a correlation-workspace memory estimate;
- explicit prepare/evaluate/kill lifecycle.

## Simplified implementation workflow

### Step 1: preserve the oracle and add backend selection

Add one advanced key covering both resolutions:

```text
refpick_backend=legacy|optimized|compare
```

It defaults to `legacy`.

- `legacy` executes the current picker unchanged.
- `optimized` executes the new path.
- `compare` runs both paths on the same logical input, keeps legacy production
  output authoritative, and records bounded comparison diagnostics.

Add the key through typed parameters, parsing, validation, UI metadata, picker
utilities, and job serialization. The stream path copies its parent command
into `cline_pick_extract`; verify that the copied value survives queue-script
generation rather than adding a parallel propagation mechanism.

### Step 2: implement the complete optimized path

Implement the optimized coarse and fine loops as a coherent branch:

1. Prepare the reference bank and summed-area tables once per resolution.
2. Gather known-valid windows into bounded reusable batches.
3. Call the linked CPU BLAS through the explicit `sgemm` boundary.
4. Apply residual-mean correction and variance normalization.
5. Reduce scores in legacy order and preserve exact tie behavior.
6. Map fine coordinates directly through checked unit-stride indexing,
   preserving last-update-wins behavior for duplicate coordinates.
7. Replace the grid-by-peak filtering update with logical selection indexed by
   positive `inds_offset` values.
8. Exclude invalid fine results before thresholding or coordinate mapping.
9. Do not allocate or calculate `loc_sdevs`. If a demonstrated requirement
   appears, add a final-box-only calculation after selection.
10. Reuse bounded buffers and select a thread topology that does not give both
    OpenMP and BLAS unrestricted thread teams.

Instrumentation is part of this branch from the start. The correlation engine
records coarse/fine candidate counts, reference and batch counts, a bounded
workspace estimate, and preparation/gather-statistics/GEMM/reduction times.
In `compare` mode, outer wall timers now report legacy time, optimized time,
and speedup for coarse correlation, coarse distance filtering, fine
correlation/remapping, fine distance filtering, and their combined divergent
stages. Validation runs additionally record total picker time, thread settings,
and the resolved BLAS.

### Step 3: validate with existing simulated workflows

Use the existing simulated workflow in
`src/main/commanders/test/simple_commanders_test_highlevel.f90`. It generates
embedded 6VXX/1JXY-derived data and now accepts runtime `picker` and
`refpick_backend` arguments instead of a compile-time picker switch.

Exercise `legacy`, `optimized`, and `compare` for both simulated systems before
using real data. Compare score diagnostics, ordered coordinates, scores, pick
counts, and `.box` output. The existing workflow is the primary development
validation path; add a narrower deterministic fixture only if a discrepancy
cannot be isolated from its output.

### Step 4: validate the completed path on real data

Once the complete optimized path and simulated validation are in place, run
legacy and optimized jobs on representative real micrographs. Use `compare`
for correctness diagnostics and separate jobs for valid performance timing.

Report both dense and sparse data, typical and 240-reference banks, coarse and
fine timing, whole-picker latency, peak memory, pick matching, downstream
particle/class-average quality, and stream queue behavior. Measurements on
macOS Accelerate do not replace deployment-host OpenBLAS measurements.

Keep `legacy` as the default until the real-data evidence is reviewed and a
default change is explicitly approved.

## Validation gates

### Numerical and structural gates

- Summed-area sums, means, and population variances match a direct
  double-precision two-pass oracle at every boundary, for constant,
  near-constant, gradient, impulse, and large-offset inputs.
- Complete BLAS score blocks match scalar Pearson scores for asymmetric
  dimensions, batch tails, one and many references, deliberately nonzero
  reference residual sums, negative scores, invalid windows, and exact ties.
- Initial score tolerances are maximum absolute error `2e-4` and RMS error
  `5e-5`; tighten them from observed cross-platform results. Any relaxation
  needs a recorded numerical cause.
- Every non-tied case selects the same final coordinate. Equal candidate scores
  follow legacy traversal order.
- Direct mapping passes boundary and round-trip tests, preserves duplicate
  update order, and never indexes a zero `inds_offset` entry.
- Zero-valid and mixed-validity fine neighborhoods complete without selecting
  the full grid or indexing `positions(0,:)`.
- The optimized path performs zero per-coarse-position `avg_loc_sdev` calls.

### Simulated and real picker gates

- Cover black/white contrast, ROI and background subtraction on/off, edges,
  clipped and full fine neighborhoods, one and many references, duplicate
  mapped cells, and competing nearby coarse peaks.
- On simulated fixtures, compare ordered coordinates, scores, score grids,
  pick counts, and `.box` output. Every mismatch is either resolved or recorded
  as an explicitly demonstrated floating-point tie/rounding case.
- On real data, require one-to-one pick matching with every mismatch triaged;
  there is no unexplained percentage budget for this
  mathematically-equivalent implementation.
- The working performance target for the 240-reference workload is at least a
  fivefold correlation-kernel reduction, no material regression for small
  workloads, bounded memory, a non-growing stream queue, and picker p95 below
  70 percent of micrograph inter-arrival time.
- Label results by build type, host, resolved BLAS, OpenMP threads, and BLAS
  threads. Do not claim Linux or BOX validation without observed output.

## Expected file ownership

- `src/main/pick/simple_pickref.f90`: legacy dispatch and shared final picker
  policy; keep the legacy scoring path intact.
- `src/main/pick/`: optimized picker/correlation type with bounded state and
  explicit lifecycle.
- shared numerical utilities: narrow LP64 CPU `sgemm` interface using the BLAS
  already linked into SIMPLE.
- `src/main/image/`: only a narrow valid-window primitive if the optimized
  picker cannot use an existing sanctioned API.
- `src/main/params/`, `src/main/ui/`, `simple_picker_utils.f90`, picker
  strategy, and stream job creation: backend key definition and propagation.
- `production/tests/` and the high-level test commander/UI: numerical,
  simulated-picker, and comparison validation.

Correlation mathematics remains in the picking domain. Stream modules, UI
metadata, parameter parsing, and commanders only propagate policy and inputs.

## Deferred work

The full-grid-to-sparse-grid redesign, Fourier/Roseman coarse search, reference
pruning, and search-space reductions are not prerequisites and are not part of
this implementation. Reopen them only if the completed real-data profile
shows a remaining need.

## Source anchors

- Picker loops, iterator setup, filtering, and remapping:
  `src/main/pick/simple_pickref.f90`.
- Sampling constants and picker construction:
  `src/main/pick/simple_picker_utils.f90`.
- Local standard deviation and Pearson normalization:
  `src/main/image/simple_image_calc.f90`.
- Window extraction: `src/main/image/simple_image_geom.f90`.
- CPU BLAS/LAPACK discovery and ABI: `cmake/Dependencies.cmake`.
- Stream child-command creation:
  `src/main/stream/simple_stream_p04_refpick_extract_new.f90`.
- Simulated workflow:
  `src/main/commanders/test/simple_commanders_test_highlevel.f90`.

Repository policy leaves compilation and runtime testing to the user unless
explicitly requested. During implementation, agents should use source checks
and language-server diagnostics, provide the exact commands needed for user
validation, and record only observed results as passed.
