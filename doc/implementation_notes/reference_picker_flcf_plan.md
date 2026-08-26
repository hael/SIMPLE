# Reference-picker acceleration note

Status: accepted and promoted to the sole production implementation on
2026-08-26. The former `legacy|optimized|compare` runtime selection has been
retired.

This is the single living design and implementation record for reference-picker
acceleration. Git history preserves the development-time scalar oracle,
comparison machinery, and earlier Roseman/FLCF assessment.

## Final contract

Reference-based picking has one implementation. There is no
`refpick_backend` CLI parameter and no runtime fallback to the retired scalar
path.

The production implementation must preserve:

- the complete-square Pearson correlation score;
- all coarse and fine candidate positions;
- all prepared references, including rotations and mirrors;
- existing reference preparation, contrast handling, filtering, ROI behavior,
  peak thresholds, distance filtering, and output coordinates;
- deterministic maximum and duplicate-coordinate update behavior;
- invalid-reference and zero/near-zero-window handling;
- existing `.box`, thumbnail, project, distributed, and stream outputs.

The acceleration changes how the exhaustive score is evaluated, not the
scientific search space.

## Production implementation

### Correlation engine

`simple_pickref_corr_batch` owns one bounded correlation workspace per
prepared picker resolution. It:

1. copies the real micrograph into owned storage;
2. removes its global mean and applies a bounded scale;
3. builds double-precision summed-area tables for the image and its square;
4. flattens the normalized reference bank and records residual reference sums;
5. gathers candidate windows into bounded column batches;
6. evaluates reference/window dot products through the existing LP64
   single-precision BLAS `SGEMM` boundary;
7. applies local mean and population-standard-deviation normalization;
8. reduces each candidate across valid references.

For reference matrix `S`, raw candidate matrix `W`, and `N` box pixels:

```text
mu(c)      = sum(T_c) / N
var(c)     = sum(T_c*T_c) / N - mu(c)^2
G          = transpose(S) * W
score(r,c) = (G(r,c) - mu(c) * sum(S(:,r))) / (N * sqrt(var(c)))
```

The residual-reference-sum term preserves the scalar Pearson definition when a
single-precision normalized reference does not sum to exactly zero.

### Coarse picking

`pickref%match_boximgs` always:

- prepares one batched correlation engine;
- scores the valid iterator positions;
- maps scores back through positive `inds_offset` entries;
- retains the unthresholded score copy used by peak/non-peak statistics.

The iterator allocation counts valid positions exactly. It no longer carries
the development-time over-allocation branch.

### Fine refinement

`pickref%refine_upscaled` always:

- enumerates the same clipped fine neighborhoods as the former scalar loop;
- scores every candidate against every valid reference in bounded batches;
- retains the first strict-greater maximum in candidate traversal order;
- excludes all-invalid neighborhoods;
- maps valid unit-stride coordinates directly into the fine score grid;
- preserves last-update-wins behavior when refined peaks share a grid cell.

### Distance filtering

Distance filtering retains the existing score-ordered greedy selection policy.
The final score-grid update uses an indexed logical keep set rather than
searching the retained peak list independently at every grid cell. Invalid
zero entries in `inds_offset` are never used as position indices.

### Removed development machinery

The production picker no longer contains:

- scalar coarse correlation or scalar fine-refinement implementations;
- legacy/optimized dispatch branches;
- compare-mode score grids, timers, or benchmark reports;
- backend state in `pickref`;
- per-position `avg_loc_sdev` calculation and its unused storage;
- the disabled local-standard-deviation outlier helper;
- backend propagation through parameters, UI, commanders, strategies,
  distributed commands, stream commands, or test commands.

The regular coarse/fine correlation engine counters remain useful operational
diagnostics and do not select behavior.

## Ownership

- `src/main/pick/simple_pickref.f90` owns picker policy, candidate
  enumeration, thresholding, refinement, and distance selection.
- `src/main/pick/simple_pickref_corr_batch.f90` owns bounded batched Pearson
  evaluation and its temporary numerical state.
- `src/utils/math/simple_linalg.f90` owns the narrow BLAS ABI wrapper.
- Picker utilities, commanders, stream stages, and UI code invoke the single
  behavior without duplicating numerical policy.

No hidden global backend state is introduced. The correlation engine retains
explicit `new`/score/`kill` lifecycle symmetry.

## Review findings

The implementation review established:

1. The dominant fine-refinement loop evaluates up to 169 positions per coarse
   peak and every expanded reference. Batching that exact work addresses the
   reported rate limiter without adopting Fourier local correlation.
2. Summed-area statistics eliminate repeated full-window mean and variance
   passes.
3. Dense matrix multiplication improves cache and vector utilization while
   retaining exhaustive scoring.
4. Direct fine-grid mapping removes the former per-peak full-grid
   `dists`/`minloc`/`where` remapping.
5. The accepted path does not need coarse-position local standard deviations:
   the internal consumer was disabled and no repository caller used the public
   getter.
6. Roseman's FLCF remains mathematically relevant but is unnecessary for the
   accepted production path and is not part of this implementation.

## Validation record

During development, the reported 6VXX comparison covered ten micrographs and
100 references. All coarse and fine comparisons were within the `2e-4`
acceptance tolerance. The recorded worst errors were:

| Stage | Maximum absolute error | RMS error |
| --- | ---: | ---: |
| Coarse | `1.788e-6` | `8.639e-8` |
| Fine | `4.292e-6` | `5.222e-7` |

The user subsequently confirmed that the picker work is validated and
authorized removal of the backend selector and scalar oracle. That
confirmation is the promotion decision; this cleanup task does not claim new
compiler, runtime, Linux, or BOX results that were not directly observed here.

## Validation criteria for future changes

Future changes to the sole picker must continue to check:

- summed-area local statistics against a double-precision direct oracle;
- batched score blocks against complete-square scalar Pearson scores;
- asymmetric layouts, batch tails, invalid references, constant and
  near-constant windows, negative scores, and exact ties;
- coarse and fine coordinate equality on deterministic fixtures;
- edge-clipped refinement neighborhoods and duplicate mapped cells;
- black/white contrast, ROI/background combinations, one and many references;
- representative real micrographs, stream queue behavior, peak memory, and
  resolved BLAS/OpenMP thread settings.

Any future search-space reduction, reference pruning, threshold change, or
Fourier/Roseman backend is a new scientific behavior and requires its own
explicitly approved contract.

## Source anchors

- Picker policy and mapping: `src/main/pick/simple_pickref.f90`.
- Batched Pearson engine:
  `src/main/pick/simple_pickref_corr_batch.f90`.
- BLAS boundary: `src/utils/math/simple_linalg.f90`.
- Picker construction: `src/main/pick/simple_picker_utils.f90`.
- Stream invocation:
  `src/main/stream/simple_stream_p04_refpick_extract_new.f90`.
- Simulated workflow:
  `src/main/commanders/test/simple_commanders_test_highlevel.f90`.
