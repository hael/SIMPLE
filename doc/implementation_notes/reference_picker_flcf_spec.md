# Reference-picker acceleration specification

Status: planning assessment; no production behavior is authorized by this
document.

The implementation plan is in
[reference_picker_flcf_plan.md](reference_picker_flcf_plan.md).

## Decision summary

The work splits into two tiers, in this order:

1. **Exact structural removals.** A group of small, provably
   output-identical changes in `simple_pickref` that delete dead work and
   replace whole-grid bookkeeping with direct indexing. These need no new
   kernel, no summed-area table, no BLAS, and no opt-in control. The largest
   of them removes a per-coarse-position local-standard-deviation calculation
   whose result is never consumed by any caller.
2. **One exact batched correlation engine, applied to both correlation
   loops.** Summed-area local statistics plus bounded dense matrix products,
   used by the coarse `match_boximgs` loop and the fine `refine_upscaled`
   loop, which are structurally identical.

Tier 2 changes execution, not the score definition or the search space, which
makes the current scalar path an exact parity oracle for it.

A Fourier replacement of the coarse search along the lines of Roseman's fast
local correlation function is out of scope. It is recorded in
[Appendix A](#appendix-a-fourier-coarse-search-not-in-scope) so the reasoning
is not lost, but nothing in this SPEC depends on it and no phase is gated on
it.

All new kernel paths must be opt-in. The current scalar implementation remains
the default and authoritative result until the acceptance gates below pass.
Tier 1 changes are not opt-in, because each is either the removal of a dead
computation or an indexing identity; each must be separately proved
output-identical.

## Problem statement

Reference picking is reported to limit stream throughput. The production
picker has two resolutions:

- coarse search at 4 A/pixel and three-pixel stride;
- local refinement at 2 A/pixel and one-pixel stride.

The fine search radius is derived from the coarse offset and sampling ratio.
With the current constants it is six fine pixels in each direction, or at
most 169 candidate windows per retained peak
(`src/main/pick/simple_picker_utils.f90:15-17,47-54` and
`src/main/pick/simple_pickref.f90:668-713`).

Both correlation loops have the same shape. For every candidate position the
current path:

- calls the general `window_slim` routine, which allocates small indexing
  arrays and clears/copies an image;
- makes several full-box passes to calculate the window mean, subtract it,
  calculate variance, and divide by standard deviation;
- loops serially across every reference and performs a full-box dot product.

The stream normally expands ten medoids into twelve rotations plus mirrors,
giving up to 240 references
(`src/main/stream/simple_stream_p04_refpick_extract_new.f90:243-251` and
`src/main/strategies/parallelization/simple_pick_strategy.f90:352-382`).

Correlation is not the only cost, and on some workloads it is not the largest.
Three non-correlation costs must be accounted for before any kernel work.

### Dead local-standard-deviation calculation

`match_boximgs` calls `avg_loc_sdev(self%offset)` for every coarse position
(`src/main/pick/simple_pickref.f90:306`). That routine
(`src/main/image/simple_image_calc.f90:104-150`) is a naive `O(N*w^2)`
two-pass loop: for each box pixel it sums a `(2*offset+1)^2` window twice.

Its result, `self%loc_sdevs`, has exactly two consumers:

- `remove_outliers` (`src/main/pick/simple_pickref.f90:447-471`), which
  `refpick` explicitly disables (`src/main/pick/simple_pickref.f90:113`);
- `get_loc_sdevs` (`src/main/pick/simple_pickref.f90:550-555`), which has no
  caller anywhere in `src/` or `production/`.

The calculation is therefore dead in the shipped picking path, and its cost is
the same order as the entire coarse correlation. Phase 0 must time it
separately; the expected outcome is that removing or lazily gating it is the
single largest available win, at zero numerical risk.

### Whole-grid bookkeeping at fine resolution

The refinement `pickref` instance runs with `offset = 1`, so its score grid
spans the whole micrograph at unit stride while holding results for only the
retained peaks. Three consequences follow:

- `refine_upscaled` maps each refined position back to the grid by an OpenMP
  loop over `self%nboxes` computing `euclid`, then `minloc`, then an
  array-wide `where` (`src/main/pick/simple_pickref.f90:723-733`). This is
  `O(npeaks * grid)`.
- `distance_filter` updates `box_scores` with a `collapse(2)` loop over the
  whole grid containing an inner linear search over the retained peaks
  (`src/main/pick/simple_pickref.f90:427-441`). This is also
  `O(npeaks * grid)` and is not addressed anywhere else in this document's
  kernel work.
- `refine_upscaled` declares `real :: dists(self%nboxes)`
  (`src/main/pick/simple_pickref.f90:677`), an automatic array sized by the
  full grid. At production geometry this is tens of megabytes on the stack and
  is a stack-overflow exposure on any host that does not raise `ulimit -s`.

`setup_iterators` also over-counts `nboxes` in its first pass
(`src/main/pick/simple_pickref.f90:246,248` increment unconditionally and then
again under the mask), roughly doubling the `positions` allocation and the
`dists` loop bound.

### Attribution and Amdahl bound

Coarse correlation work is approximately
`n_coarse_positions * R * N_coarse` and is independent of particle count.
Fine correlation work is approximately `npeaks * 169 * R * N_fine`. With
typical production geometry, fine work exceeds coarse by roughly six times at
high particle counts but is comparable to coarse on sparse micrographs.

A refinement-only optimization is therefore bounded by the coarse share: a
five-times refinement speedup yields roughly three times overall on dense
micrographs and under two times on sparse ones. Sparse micrographs are exactly
the ones that set stream tail latency. Since `match_boximgs` and
`refine_upscaled` have identical structure, one engine must serve both; there
is no justification for optimizing only one and deferring the other to an
FFT-based rewrite.

### Latent threshold blowup

`refine_upscaled` sets `self%t = minval(scores_refined)`
(`src/main/pick/simple_pickref.f90:723`) immediately before setting the whole
grid to `-1.`. If any refined peak leaves `box_score = -1.` — reachable when
every one of its candidates fails the near-zero-variance test — then `t` is
`-1.`, and `npeaks = count(box_scores >= t)` returns the entire grid.
`distance_filter` then packs the whole grid, including `inds_offset == 0`
entries, into a greedy `O(n^2)` loop. The failure mode is a stalled stream
stage, not a wrong score. Such a peak also leaves `pos_refined = [0,0]`, which
writes the grid origin.

This must be fixed under its own change with the scalar oracle updated to
match. It must not be frozen into the batched backend's parity contract.

## Scientific score contract

For a square reference with `N` pixels and micrograph window `T_c`, SIMPLE
computes

```text
r(r,c) = sum((S_r - mean(S_r)) * (T_c - mean(T_c)))
         / sqrt(sum((S_r - mean(S_r))^2)
                * sum((T_c - mean(T_c))^2))
```

The prepared reference is normalized to zero mean and population standard
deviation one over the complete square by `prenorm4real_corr_3`. The current
window is normalized the same way, and `real_corr_prenorm_3` returns the dot
product divided by `N`
(`src/main/image/simple_image_calc.f90:1558-1572,1609-1614`).

The complete square is part of the contract. Automasked zero-valued template
background pixels still participate after complete-box mean subtraction. No
optimization in this SPEC may replace that footprint with the automask.

## Equivalent batched real-space formulation

Let `S` be an `N x R` matrix whose columns are the prepared references,
and let `W` be an `N x C` matrix of raw candidate windows. For candidate
`c`, calculate its raw sum, sum of squares, mean, and population standard
deviation:

```text
mu(c)    = sum(T_c) / N
var(c)   = sum(T_c*T_c) / N - mu(c)^2
sigma(c) = sqrt(var(c))
```

The sums can be read in constant time from double-precision summed-area tables
of the prepared micrograph and its square. Calculate all reference/window
numerators in a candidate batch with

```text
G = transpose(S) * W
numerator(r,c) = G(r,c) - mu(c) * sum(S(:,r))
score(r,c) = numerator(r,c) / (N * sigma(c))
```

The `mu * sum(S)` correction must be retained even though a normalized
reference sum should be close to zero. It preserves the current formula in the
presence of single-precision residual mean. A BLAS matrix multiplication, or
an equivalently verified dense kernel, replaces many small scalar dot-product
calls without changing the mathematical score.

Candidate windows must be gathered in the same pixel order as reference
flattening. Batches must be bounded; the implementation may maintain only the
current winning score/position and need not store all `R x C` scores.

The same formulation serves the coarse loop with `C` enumerating strided
coarse positions instead of fine candidates.

### Required amplitude precondition

The scalar path centers each window in single precision before taking the dot
product. The batched path takes the dot product against the raw window and
subtracts `mu * sum(S)` afterwards. These are identical in exact arithmetic
and materially different in single precision when the micrograph carries a
large DC level or a large amplitude scale.

`mic_shrink` carries both. It is produced by multiplying `mic_raw` by
`product(ldim_raw)` before the forward transform to avoid underflow
(`src/main/pick/simple_pickref.f90:204`), clipping in Fourier space, and
inverse transforming without the matching division — only `mic_raw` is divided
back (`src/main/pick/simple_pickref.f90:220-221`). The subsequent band-pass
uses `hp = 0.` and does not remove DC
(`src/main/pick/simple_pickref.f90:229`).

Because the Pearson score is invariant under a global shift and a global
positive scale, the implementation must, once per prepared micrograph and
before building the summed-area tables, subtract the global mean and rescale
to a bounded amplitude. This is exact, and without it the numerical gates
below are at risk for reasons unrelated to the kernel implementation.

Phase 0 must record the observed mean and range of `mic_shrink` at both
resolutions so the rescaling constant is chosen from measurement.

## Required behavior

### Exact structural removals

Each of the following is an independent change with its own proof of
output-identity. None takes an opt-in control.

1. Remove or lazily gate the `avg_loc_sdev` call and the `loc_sdevs` /
   `loc_sdevs_mem` state it feeds, subject to a decision on whether
   `remove_outliers` and `get_loc_sdevs` are retained as future capability.
   If retained, they must be gated so the shipped path does not pay for them.
2. Correct the `setup_iterators` first-pass `nboxes` overcount.
3. Replace the `refine_upscaled` distance/`minloc`/`where` remap with a direct
   coordinate-to-grid mapping, and remove the `dists` automatic array.
4. Replace the `distance_filter` inner linear search with a logical selection
   array indexed through `inds_offset`.
5. Guard the `self%t = minval(scores_refined)` degenerate case so an
   all-invalid refined peak cannot admit the entire grid.

A direct coordinate-to-grid mapping may replace the current distance/minimum
search only after the stride-one and origin invariants are asserted. If an
input violates those invariants, the code must use a correct fallback or fail
with a clear internal error.

Item 3 is a local patch for a representation problem: the fine instance
materializes a full unit-stride grid to hold results for the retained peaks
only. The alternative is a sparse peak list at fine resolution, whose only
consumers are `refine_upscaled`, `distance_filter`, `get_positions`, and
`report_boxfile`. That decision must be taken before the scalar seam is cut,
because it determines where the seam goes. If the sparse representation is
adopted, items 3 and 4 are subsumed by it.

### Correlation-backend behavior

The batched backend must preserve, at both resolutions:

- the preprocessing and exact candidate rectangle or coarse stride;
- evaluation of every current candidate position against every reference;
- complete-square normalization and Pearson scoring;
- current invalid/near-zero-variance handling;
- maximum-score selection and deterministic tie behavior in the legacy
  candidate/reference traversal order;
- coordinate origin and grid indexing;
- behavior when multiple refined peaks map to the same score-grid location;
- ROI, threshold estimation, peak selection, distance filtering, and all final
  outputs.

### Public selection

Use one advanced control covering both correlation loops, because they share
one engine:

- provisional `refpick_backend=scalar|batched|compare`.

It defaults to the legacy implementation. `compare` keeps the legacy result
authoritative and records bounded diagnostics.

The key must use the normal typed `parameters`, registry, validation, UI,
commander/strategy, distributed child-command, and stream propagation paths.
The stream master constructs the child picking command explicitly, so
propagation cannot be implicit.

### Ownership and lifecycle

`pickref` owns search policy and final reduction. Put batched numerical state
behind a narrow picking-domain type with explicit prepare/evaluate/kill
lifecycle. It should own:

- double-precision summed-area tables;
- the flattened normalized reference bank and per-reference residual sums;
- bounded candidate and result workspaces;
- timing and work counters.

One instance serves one prepared resolution, so the coarse and fine `pickref`
instances each own one.

A specialized valid 2D window gather may live in the image layer if it is a
general representation primitive. It must not expose image internals broadly.
BLAS layout and precision conversion should go through a narrow
linear-algebra boundary.

No such boundary exists today: the repository contains no `gemm`/`sgemm`
call, and every dense product is `matmul` (`src/main/pca/*`,
`src/main/simple_symanalyzer.f90`). The library itself is already there. BLAS
and LAPACK are required dependencies (`cmake/Dependencies.cmake:112-148`) and
are linked into every SIMPLE binary today; on the current development host
that resolution is Apple Accelerate, whose BLAS exports the classic LP64
`sgemm_` symbol matching the ABI that file documents. Writing the interface
module is therefore the only new work, and it must declare the LP64 signature
explicitly rather than relying on an implicit interface.

Do not add hidden global state. Reference caches must be keyed by dimensions,
sampling, preprocessing, and identity, with explicit invalidation and memory
bounds.

### Threading

The current outer OpenMP loops distribute positions or peaks. A threaded BLAS
call inside them can oversubscribe the machine. The implementation must select
and measure an explicit topology, initially favoring outer OpenMP with
single-threaded BLAS, and must set the vendor thread-count control explicitly
rather than assuming the linked library's default.

That control is vendor-specific, and the development and deployment hosts do
not share one. Apple Accelerate, which is what BLAS resolves to on the current
macOS development host, has no `OPENBLAS_NUM_THREADS`; its knob is
`VECLIB_MAXIMUM_THREADS`, and on Apple Silicon it dispatches small and medium
dense products to AMX units that are shared across a core cluster, so
concurrent `sgemm` calls from many OpenMP threads contend for hardware that an
x86 OpenBLAS build does not have. Batch size and thread topology must
therefore be tuned on a host representative of the deployment target, and any
figure obtained on the macOS development host must be labelled as such.

A larger-batch/threaded-BLAS alternative may be tested separately. Results
must be deterministic within the numerical gates.

## Compare-mode contract

The useful shadow comparison is at each correlation boundary:

1. prepare common micrograph, reference, and position inputs;
2. run scalar and batched correlation sequentially, or with partitioned
   resources only after thread safety is proved;
3. compare per-position best score, winning position, and final score-grid
   update;
4. run shared downstream policy with the scalar result authoritative;
5. emit bounded per-micrograph diagnostics.

Literal simultaneous execution in one process is not an initial requirement.
It risks OpenMP/BLAS oversubscription, races in picker module state, and biased
timings. On-the-fly means shadowing the same input before it leaves the stage.
Separate scalar and batched runs are the valid performance comparison.

Live-stream shadowing must have a configurable cadence. Diagnostics should
include geometry, reference/peak/candidate counts, stage times, score errors,
winning-position agreement, displacement, final pick matching, and peak
memory.

## Acceptance criteria

### Structural-removal gates

1. Each exact structural removal is separately shown to produce byte-identical
   `.box` output, pick counts, and scores on the deterministic fixtures, with
   the sole permitted exception of the degenerate-threshold guard, whose
   behavior change must be demonstrated only on an input that triggers the
   degenerate case.
2. The `avg_loc_sdev` removal is accompanied by a recorded decision on the
   future of `remove_outliers` and `get_loc_sdevs`.
3. Direct coordinate mapping is proved by edge and round-trip tests, with an
   asserted invariant and a correct fallback or clear internal error when the
   invariant does not hold.

### Exact-kernel gates

4. Summed-area mean and population variance match a double-precision direct
   oracle for every valid test window. Small negative variance from
   cancellation follows one documented clamp/error policy.
5. Every scalar reference/window score is compared with the batched score on
   deterministic random and adversarial inputs. Provisional single-precision
   gates are maximum absolute error at most `2e-4` and RMS error at most
   `5e-5`, measured with the global amplitude precondition applied. If those
   gates cannot be met, the first suspect is the precondition, not the
   tolerance.
6. Constant and near-constant references/windows agree on valid/invalid
   classification and score behavior.
7. Asymmetric data prove the candidate flattening and matrix orientation.
8. Multiple-reference maxima, exact ties, traversal order, and duplicate
   coordinate updates match legacy behavior.

These tolerances are provisional and should be tightened from observed
cross-platform results. Relaxation requires a recorded numerical explanation.

### Picker gates

9. `refpick_backend=scalar` is output-identical to the pre-refactor picker and
   remains the default.
10. On deterministic fixtures, batched and scalar correlation choose the same
    winning coordinate for every non-tied position at both resolutions;
    declared equal-score ties must follow legacy traversal semantics.
11. On representative production micrographs, the batched path reproduces the
    legacy pick set. Because the batched route is exact, every mismatch must be
    individually triaged and explained as a declared tie or a documented
    single-precision rounding event. There is no percentage budget for
    unexplained mismatches. Percentage-agreement gates apply only to the
    deferred behavior-changing optimizations, which are not authorized here.
12. Gates cover black/white contrast, one and many references, ROI on/off,
    background subtraction on/off, edges, density thresholds, and competing
    nearby coarse peaks.
13. `compare` never changes production `.box`, project, thumbnail, or stream
    output.

### Performance and operations gates

14. All attribution and performance figures come from a Release build.
    Debug builds (`-O0 -g -fcheck=do,mem`) distort tight scalar loops and
    array intrinsics in opposite directions and must not be used for any gate
    in this section.
15. Profiles separately report, at both resolutions: window gathering, local
    statistics, the dead-or-gated local-standard-deviation call, reference dot
    products, score reduction, coordinate remapping, distance filtering,
    iterator setup, and total time.
16. The structural-removal tier is measured and reported on its own, before
    any kernel work, so its contribution is not attributed to the kernel.
17. On the representative 240-reference stream workload, the batched
    correlation kernel should reduce correlation wall time by at least five
    times. The prior two-times target is too weak for replacing 240 scalar dot
    products with one well-shaped dense product, and would allow a poor packing
    implementation to pass. Release judgment additionally uses end-to-end
    picker time and queue stability: the target is a non-growing queue and
    picker p95 below 70 percent of micrograph inter-arrival time.
18. Reported speedups must be given for both a dense and a sparse micrograph,
    since the coarse/fine split and therefore the Amdahl bound differ
    materially between them.
19. Low-reference/small-candidate cases must not regress materially; a
    measured crossover may select the scalar kernel internally without
    changing results.
20. Candidate batches and reference storage remain within the configured
    worker memory budget. Peak RSS must be reported, including the reduction
    from the `nboxes` overcount fix and from any sparse-representation change.
21. No compilation, unit-test, integration-test, Linux, or BOX claim is made
    until its output has actually been observed.

## Deferred behavior-changing optimizations

These can reduce work further, but are not part of either exact tier:

- refine only the coarse winner or top-K references;
- retain references within a coarse-score margin;
- hill-climb, multiresolution, or fit a local surface instead of testing all
  169 fine positions;
- reduce rotations/mirrors;
- change thresholds or distance filtering.

A reference that loses at 4 A/pixel can win at 2 A/pixel, and a non-convex
local score surface can defeat hill climbing. Each proposal therefore needs a
separate opt-in, accuracy study, and SPEC revision.

## Evidence used for this assessment

- Coarse matching loop, including the dead local-standard-deviation call:
  `src/main/pick/simple_pickref.f90:285-338`, specifically line 306.
- `loc_sdevs` consumers: `src/main/pick/simple_pickref.f90:447-471` and
  `550-555`; `remove_outliers` disabled at
  `src/main/pick/simple_pickref.f90:113`.
- `avg_loc_sdev` implementation:
  `src/main/image/simple_image_calc.f90:104-150`.
- Fine refinement loop, remap, and degenerate threshold:
  `src/main/pick/simple_pickref.f90:668-740`, specifically 677 and 723-733.
- Whole-grid distance-filter update:
  `src/main/pick/simple_pickref.f90:427-441`.
- Iterator setup and `nboxes` overcount:
  `src/main/pick/simple_pickref.f90:237-283`, specifically 246 and 248.
- Micrograph amplitude handling:
  `src/main/pick/simple_pickref.f90:204,220-221,229`.
- Fine/coarse sampling constants:
  `src/main/pick/simple_picker_utils.f90:15-17,47-54`.
- Window extraction:
  `src/main/image/simple_image_geom.f90:88-122`.
- Complete-box normalization and correlation:
  `src/main/image/simple_image_calc.f90:1558-1572,1609-1614`.
- Existing summed-area implementation pattern:
  `src/main/image/simple_image_calc.f90:228-277`.
- Absence of a BLAS boundary, and the already-linked BLAS: repository-wide
  absence of `gemm` calls; `cmake/Dependencies.cmake:112-148`; Accelerate
  present in the link line of the built `simple_exec` on the development host,
  exporting LP64 `sgemm_`.
- Stream reference expansion:
  `src/main/stream/simple_stream_p04_refpick_extract_new.f90:240-264` and
  `src/main/strategies/parallelization/simple_pick_strategy.f90:323-382`.
- Ruben Meana-Paneda's simulated workflow and existing Fourier convention
  tests, described in the implementation plan.

## Appendix A: Fourier coarse search, not in scope

Recorded for completeness. No phase, gate, or control in this SPEC depends on
anything below, and no work here is authorized.

Roseman's fast local correlation function replaces the strided real-space
coarse search with FFT-based maps. For a complete-square mask `M`, it computes
shared local sum and sum-of-squares maps by convolution and a reference
numerator map by cross-correlation. With SIMPLE's current normalized
reference, the score would be

```text
local_mean(x) = corr(M,T)(x) / N
local_var(x)  = corr(M,T*T)(x) / N - local_mean(x)^2
score_r(x)    = corr(S_r,T)(x) / (N * sqrt(local_var(x)))
```

This is mathematically compatible with the complete-square Pearson contract in
this document, but only if transform scale, template placement, conjugation,
lag indexing, valid origins, and variance policy are proved by
direct-versus-FFT tests. SIMPLE's FFT conventions must not be assumed from a
textbook formula.

It is out of scope because it addresses only the coarse loop, and the coarse
loop is already covered by the tier-2 engine at far lower risk. Should a
future profile show coarse scoring dominating even after tier 2, this would be
the starting point, and it would need its own SPEC, its own opt-in control,
and a percentage-agreement gate rather than the exact-parity gate used here.

Reference: A. M. Roseman, "Particle finding in electron micrographs using a
fast local correlation algorithm," *Ultramicroscopy* 94 (2003), 225-236.
