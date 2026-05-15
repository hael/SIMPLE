# Microchunk Quality Refactor Notes

**Date:** May 15, 2026
**Context:** Integration note for using the current `chunk_default_v2` class-average quality model inside `cluster2D_microchunked` while preserving the streaming/microchunk orchestration.

## Short Summary

`cluster2D_microchunked` owns the streaming orchestration: it partitions the stream into bounded projects, runs pass-1 and pass-2 ab initio 2D tiers, builds a reference chunk, matches later chunks against that reference, tracks completion with sentinel files, and combines match chunks into the final project. The class-average quality decision should come from the shared `simple_cavg_quality` backend.

`simple_cavg_quality` provides the quality decision layer:

- extract a compact feature vector for every class average
- robustly normalize features within each run
- keep only conservative hard rejections as pre-clustering vetoes
- use k-medoids only inside normalized scalar-feature space to identify good and bad class-average groups
- derive an automatic score threshold from the cluster means
- preserve raw and normalized feature tables for validation

The main refactor opportunity is to keep the microchunk lifecycle as-is and route only the class-average rejection policy through the quality backend. This should be controlled by a single internal logical flag in the microchunk implementation, not by a new command-line option and not by a second rejection-type abstraction.

The quality backend is the current chunk model:

- instantiate `cavg_quality_model`
- initialize it with `CAVG_QUALITY_MODEL_CHUNK_DEFAULT`, currently `chunk_default_v2`
- call `evaluate_cavg_quality`
- propagate `quality%states`, `quality%scores`, and cluster annotations into the standard project fields for class-average selection

The selected model owns the decision boundary, feature policy, hard-reject behavior, and any later retrained chunk default. The microchunk code should only decide whether to use the quality backend.

## Current Quality Policy

### 1. Quality Evidence Is Scalar And Vectorized

The supported direction is scalar quality metrics as features in a normalized vectorized learning analysis. Pairwise histogram/spectrum matrices and relational clustering-style quality evidence are not part of the active policy. Clustering is used as an automatic thresholding helper within scalar feature space, not as a pairwise image-similarity model.

Each feature is allowed to be weak alone. The model learns how stable scalar evidence combines across population support, resolution, mask geometry, local variance, correlation/FRC proxy, center-edge signal, particle presence, and internal detail.

### 2. Robust Feature Normalization Is Required

The useful unit of comparison is not the raw scalar metric. It is the metric relative to the distribution of class averages in the current project or chunk. The successful implementation uses median/MAD normalization over non-hard-rejected classes, then clips normalized values to limit outlier leverage.

This makes feature values comparable across datasets with different particle counts, different absolute class-average contrast, and different failure modes.

### 3. Conservative Hard Rejection Stays Conservative

The quality module only hard-rejects pathological cases:

- zero population
- stored class resolution worse than 40 A
- no valid connected component after foreground segmentation
- foreground-component centroids outside the mask radius
- more than 10 pixels of the largest foreground component outside the mask disc
- essentially zero local variance in both foreground and background

Aggressive hard rejections are difficult to recover from because they bypass the multivariate decision. Most real examples should enter feature-space scoring, even if they look weak by one scalar metric.

Resolution is not otherwise removed from the model. The hard gate only handles catastrophic `res > 40 A` failures; ordinary resolution differences remain active scalar evidence through `neg_log_res`.

### 4. Feature Families Define Model Eligibility

The model infrastructure excludes histogram/spectrum pairwise matrices. The feature inventory has explicit families, and diagnostic-only entries are reserved for hard-reject explanation or active feature hypotheses. Local variance, internal detail, resolution, population support, and presence remain eligible model evidence.

The feature table keeps cheap scalar diagnostics such as `mask_inside`, `single_component`, connected-component shape, central presence, and internal-detail diagnostics from a broad 20-6 A band-pass.

### 5. The Chunk Profile Is Local-Signal Heavy

The chunk default is a scalar feature-space model with learned emphasis on population/resolution support, local foreground/background variance, correlation/FRC proxy, center-edge signal, particle presence, and connected-component diameter.

The practical interpretation is:

- local foreground variance is the strongest single contributor
- class population still matters, but should not dominate
- resolution helps, but is not a sufficient proxy for quality
- centering helps catch pathological classes without becoming a hard rule
- background/local context is useful when combined with foreground signal

### 6. The Boundary Margin Remains Internal

The signed boundary margin is useful for calibration, but it is an internal model parameter in `simple_cavg_quality`, not a command-line knob. Tuning should happen through the learn/promote workflow, with validation summaries checked against trusted manual selections.

## Microchunk Orchestration Responsibilities

`simple_stream_cluster2D_microchunked` and `simple_microchunked2D` have several strong design choices that we should preserve:

- chunk generation is isolated from chunk processing
- each microchunk is a self-contained project
- pass-1, pass-2, refchunk, and match tiers have explicit lifecycle flags
- sentinel files make restart behavior understandable
- project merging is centralized
- reference selection is tracked and exposed for downstream stream GUI use
- match chunks reuse refchunk class averages rather than reclustering globally

These are orchestration strengths. The quality-vector work should not be used as an excuse to rewrite the scheduler.

## Main Refactor Recommendation

Refactor class-average rejection in `simple_microchunked2D` so that the decision source is `simple_cavg_quality`, while the microchunk lifecycle remains unchanged.

The backend choice should be a module-local logical constant or private logical field, for example:

```text
USE_CAVG_QUALITY_BACKEND = .true.
```

When this flag is `.true.`, `microchunked2D%reject_cavgs` should use the current chunk model. This switch should not be exposed in the UI.

Concretely:

1. Add a small internal helper, for example `apply_cavg_quality_to_chunk`, that:
   - reads the chunk project and class-average stack
   - initializes a `cavg_quality_model` with `CAVG_QUALITY_MODEL_CHUNK_DEFAULT`
   - calls `evaluate_cavg_quality`
   - maps `quality%states` to `cls2D`, `cls3D`, and particle state/class fields
   - records `quality_score`, `quality_cluster`, model name, and model context where the project schema already supports it
   - writes selected/rejected class-average stacks
   - optionally writes a compact feature table for analysis/debug builds
   - updates `chunk%nptcls_selected`

2. Route `microchunked2D%reject_cavgs` through the helper when `USE_CAVG_QUALITY_BACKEND` is `.true.`.

3. Keep the existing sentinel behavior:
   - `REJECTION_FINISHED` on success
   - `REJECTION_FAILED` and `COMPLETE` on missing/mismatched class averages
   - refchunk still captures `self%refs` and `self%box`
   - match chunks still update accepted/rejected particle counters

4. Keep any non-quality-backend branch private and temporary, for debugging and side-by-side validation only.

This should make the change low-risk: the scheduler sees the same class states and particle-state side effects as before, but those states come from the validated feature-space selector.

## Suggested Implementation Stages

### Stage 1: Share Project Update Code

`exec_model_cavgs_rejection` and `microchunked2D%reject_cavgs` both need to propagate class-average decisions into project state. Extract the common state-mapping work into a small helper module or module procedure.

Candidate responsibilities:

- set `cls2D` and `cls3D` state/accept/quality fields
- call `map2ptcls_state`
- zero `class` and `class_match` for rejected particles
- write selected/rejected stacks

This reduces risk before changing the actual selection policy.

### Stage 2: Add A Logical Quality Backend Switch

Introduce a local logical choice inside `simple_microchunked2D`, not a user-facing CLI parameter:

```text
USE_CAVG_QUALITY_BACKEND = .true.
```

This keeps the backend choice explicit during validation without adding another operational knob. The true branch should use `CAVG_QUALITY_MODEL_CHUNK_DEFAULT`; it should not add `rejection_type`, and it should not carry a separate set of scalar thresholds in the microchunk layer.

### Stage 3: Use The Current Chunk Model As The Integrated Profile

The current `chunk_default_v2` profile is the default first integration candidate. Microchunk tiers differ:

- pass-1 chunks are smaller and noisier
- pass-2 chunks have already passed one quality gate
- refchunk quality is high-impact because it seeds all match chunks
- match chunks need stable, comparable decisions against a fixed reference

If validation shows the same profile is too strict for pass-1 or too permissive for ref/match, introduce stage-aware internal models in `simple_cavg_quality`, not scattered thresholds in `simple_microchunked2D`.

For example:

```text
CAVG_QUALITY_MODEL_CHUNK_PASS1
CAVG_QUALITY_MODEL_CHUNK_PASS2
CAVG_QUALITY_MODEL_CHUNK_REF
CAVG_QUALITY_MODEL_CHUNK_MATCH
```

Each model should specify weights, boundary margin, minimum score separation, feature policy, and hard-reject behavior through the normal model machinery. The extraction and normalization code should remain shared.

### Stage 4: Replace Duplicate Metric Code

`simple_cluster2D_rejector` and `simple_cavg_quality` currently duplicate several concepts:

- Otsu foreground segmentation
- connected-component filtering
- foreground/outside-mask geometry
- local variance
- population and resolution scoring

Once the quality-vector backend is validated, `simple_cluster2D_rejector` can be reduced to either:

- a compatibility wrapper around `simple_cavg_quality`, or
- a private scalar backend available only for targeted debugging/tests

This would avoid two diverging definitions of "bad class average."

### Stage 5: Standardize Diagnostics

The useful diagnostic outputs from `model_cavgs_rejection` are:

- `cavgs_quality_features.txt`
- `cavgs_quality_reference_summary.txt`
- `cavgs_quality_reference_threshold_scan.txt`
- selected/rejected class-average stacks

For microchunk streaming, always writing full diagnostic tables for every chunk may be noisy. A good compromise:

- always annotate `cls2D` with `quality`, `quality_cluster`, and `accept`
- always write selected/rejected stacks, preserving current behavior
- write feature tables only under an internal debug/analysis setting or for failed/refchunk cases
- consider writing one aggregate per-tier summary after completion

This keeps normal stream runs light while preserving the ability to debug difficult partitions.

## Integration Sketch

Current microchunk rejection:

```text
read chunk project
read class averages
cluster2D_rejector%new
reject_pop
reject_res
reject_mask
reject_local_variance
map rejected flags to project state
write selected/rejected stacks
touch REJECTION_FINISHED
```

Proposed with the logical backend flag enabled:

```text
read chunk project
read class averages
model%init_preset(CAVG_QUALITY_MODEL_CHUNK_DEFAULT)
evaluate_cavg_quality
map quality%states to project state
annotate quality fields
write selected/rejected stacks
touch REJECTION_FINISHED
```

The outer microchunk flow does not need to change: only the decision source for class-average rejection changes.

## Open Validation Questions

1. Does the current difficult-stream profile preserve behavior on easier batch-mode datasets?
2. Are pass-1 chunks too small for reliable two-cluster k-medoids decisions, or does the existing `nfit < 4` and low-separation fallback cover this?
3. Should refchunk and match chunks use the same profile, or should refchunk be slightly less permissive because it seeds all downstream matching?
4. Do we need to retain scalar population floors in pass-1 as hard pre-filters, or does the feature-vector population weight handle this safely?
5. How much diagnostic output is acceptable in long streaming runs?

## Recommended Next Step

Integrate `simple_cavg_quality` into `microchunked2D%reject_cavgs` behind the internal `USE_CAVG_QUALITY_BACKEND` logical flag and validate on difficult stream partitions. The enabled path should initialize `CAVG_QUALITY_MODEL_CHUNK_DEFAULT`, so future retraining and promotion of `chunk_default_v2` automatically changes the stream behavior without touching microchunk code. If the state propagation and sentinel behavior remain unchanged, this should be a local refactor rather than a scheduler rewrite.

The safest initial target is refchunk and match rejection, because those tiers most closely resemble the real-world validation sets. Pass-1/pass-2 can follow once we know the quality-vector behavior is stable on smaller chunks.
