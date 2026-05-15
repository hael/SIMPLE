# Microchunk Quality Refactor Notes

**Date:** May 13, 2026  
**Context:** Follow-up to the `cluster_cavgs_quality` experiments on difficult real-world stream partitions, with an eye toward improving `cluster2D_microchunked` without disturbing Joe's streaming/microchunk orchestration.

## Short Summary

Joe's `cluster2D_microchunked` code already solves the hard orchestration problem: it partitions the stream into bounded projects, runs pass-1 and pass-2 ab initio 2D tiers, builds a reference chunk, matches later chunks against that reference, tracks completion with sentinel files, and combines match chunks into the final project. The weaker piece is not the chunk scheduler. It is the class-average quality gate inside `microchunked2D%reject_cavgs`, which still uses a cascade of scalar threshold tests.

The new `simple_cavg_quality` implementation gives us a better quality decision layer:

- extract a compact feature vector for every class average
- robustly normalize features within each run
- keep only conservative hard rejections as pre-clustering vetoes
- use k-medoids in feature space to identify good and bad class-average groups
- derive an automatic score threshold from the cluster means
- preserve raw and normalized feature tables for validation

The main refactor opportunity is therefore to keep the microchunk lifecycle as-is and replace or wrap only the class-average rejection policy.

## What We Learned From The Experiments

### 1. Scalar Gates Are Useful But Too Brittle Alone

The original microchunk rejector applies independent decisions:

- population floor
- FSC resolution ceiling
- mask/connected-component rejection
- local-variance rejection

This is simple and fast, but each test is forced to be either a hard veto or irrelevant. That made tuning awkward because several metrics are weak alone but useful in combination. Population and resolution were especially tempting to over-weight because they often separate easy cases, but they do not generalize cleanly across difficult real-world stream partitions.

### 2. Robust Feature Normalization Was Crucial

The useful unit of comparison is not the raw scalar metric. It is the metric relative to the distribution of class averages in the current project or chunk. The successful implementation uses median/MAD normalization over non-hard-rejected classes, then clips normalized values to limit outlier leverage.

This made feature values comparable across datasets with different particle counts, different absolute class-average contrast, and different failure modes.

### 3. Conservative Hard Rejection Should Stay Conservative

The new quality module only hard-rejects pathological cases:

- zero population
- stored class resolution worse than 40 A
- no valid connected component after foreground segmentation
- foreground-component centroids outside the mask radius
- more than 10 pixels of the largest foreground component outside the mask disc
- essentially zero local variance in both foreground and background

This was important. Aggressive hard rejections are difficult to recover from because they bypass the multivariate decision. Most real examples should enter feature-space scoring, even if they look weak by one scalar metric.

### 4. Some Candidate Features Did Not Survive Representative Chunk Training

The current feature bank no longer includes `spectrum_dynrange`, `neg_ice_score`, or `log_fg_bg_locvar_ratio`, and the histogram/spectrum pairwise matrices have been removed from the model infrastructure. Those signals either failed to separate manual good/bad classes consistently, were redundant with retained scalar features, or inverted on difficult stream partitions.

The remaining feature table keeps cheap scalar diagnostics such as `mask_inside`, `single_component`, histogram entropy, connected-component shape, central presence, full-image contrast, and masked histogram variance. The current chunk default was promoted from the representative batch7chunk learning cycle with the `base12_pruned_plus_histvar` policy, which keeps masked histogram variance but zeros `mask_inside`, `single_component`, `cc_area_frac`, `presence`, and `log_contrast` in the model.

### 5. The Current Best General Profile Is Still Local-Signal Heavy

The chunk default has evolved into a scalar feature-space model with learned emphasis on population/resolution support, local foreground/background variance, correlation/FRC proxy, center-edge signal, histogram entropy, and connected-component diameter. The obsolete spectrum dynamic-range and ice-ring diagnostics are no longer part of the feature vector.

The practical interpretation is:

- local foreground variance is the strongest single contributor
- class population still matters, but should not dominate
- resolution helps, but is not a sufficient proxy for quality
- centering helps catch pathological classes without becoming a hard rule
- background/local context is useful when combined with foreground signal

On the representative stream-chunk learning set in batch7chunk, this became the promoted `chunk_default_v2` profile.

### 6. The Boundary Margin Should Remain Internal

The signed boundary margin is useful for calibration, but exposing it as `quality_recall_margin` on the command line was confusing. It is now an internal constant in `simple_cavg_quality`. Future tuning should happen in code, with validation summaries checked against trusted manual selections, not through casual runtime tweaking.

## What Joe's Microchunk Code Already Does Well

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

Concretely:

1. Add a small internal helper, for example `apply_cavg_quality_to_chunk`, that:
   - reads the chunk project and class-average stack
   - calls `evaluate_cavg_quality`
   - maps `quality%states` to `cls2D`, `cls3D`, and particle state/class fields
   - writes selected/rejected class-average stacks
   - optionally writes a compact feature table for analysis/debug builds
   - updates `chunk%nptcls_selected`

2. Replace the current scalar cascade inside `microchunked2D%reject_cavgs` with that helper.

3. Keep the existing sentinel behavior:
   - `REJECTION_FINISHED` on success
   - `REJECTION_FAILED` and `COMPLETE` on missing/mismatched class averages
   - refchunk still captures `self%refs` and `self%box`
   - match chunks still update accepted/rejected particle counters

4. Keep `cluster2D_rejector` temporarily as a fallback or comparison backend until the batch-mode validation confirms easy-dataset behavior.

This should make the change low-risk: the scheduler sees the same class states and particle-state side effects as before, but those states come from the validated feature-space selector.

## Suggested Implementation Stages

### Stage 1: Share Project Update Code

Right now `exec_cluster_cavgs_quality` and `microchunked2D%reject_cavgs` each know how to propagate class-average decisions into project state. Extract the common state-mapping work into a small helper module or module procedure.

Candidate responsibilities:

- set `cls2D` and `cls3D` state/accept/quality fields
- call `map2ptcls_state`
- zero `class` and `class_match` for rejected particles
- write selected/rejected stacks

This reduces risk before changing the actual selection policy.

### Stage 2: Add A Quality Backend To Microchunk Rejection

Introduce a local policy choice inside `simple_microchunked2D`, not a user-facing CLI parameter. For example:

```text
MICROCHUNK_REJECTION_BACKEND = 'quality_vector'
```

or, better, use a logical/module constant during transition:

```text
USE_QUALITY_VECTOR_REJECTION = .true.
```

This lets us keep the old scalar backend nearby for a short validation period without exposing another knob.

### Stage 3: Use Stage-Aware Profiles Only If Evidence Requires It

The current quality profile was tuned on difficult stream partitions and should be the default first integration candidate. However, microchunk tiers differ:

- pass-1 chunks are smaller and noisier
- pass-2 chunks have already passed one quality gate
- refchunk quality is high-impact because it seeds all match chunks
- match chunks need stable, comparable decisions against a fixed reference

If validation shows the same profile is too strict for pass-1 or too permissive for ref/match, introduce stage-aware internal profiles in `simple_cavg_quality`, not scattered thresholds in `simple_microchunked2D`.

For example:

```text
CAVG_QUALITY_PROFILE_PASS1
CAVG_QUALITY_PROFILE_PASS2
CAVG_QUALITY_PROFILE_REF
CAVG_QUALITY_PROFILE_MATCH
```

Each profile should specify weights, boundary margin, and minimum score separation. The extraction and normalization code should remain shared.

### Stage 4: Replace Duplicate Metric Code

`simple_cluster2D_rejector` and `simple_cavg_quality` currently duplicate several concepts:

- Otsu foreground segmentation
- connected-component filtering
- foreground/outside-mask geometry
- local variance
- population and resolution scoring

Once the quality-vector backend is validated, `simple_cluster2D_rejector` can be reduced to either:

- a compatibility wrapper around `simple_cavg_quality`, or
- a legacy scalar backend retained only for old workflows/tests

This would avoid two diverging definitions of "bad class average."

### Stage 5: Standardize Diagnostics

The most useful diagnostic outputs from `cluster_cavgs_quality` were:

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

Proposed:

```text
read chunk project
read class averages
evaluate_cavg_quality
map quality%states to project state
annotate quality fields
write selected/rejected stacks
touch REJECTION_FINISHED
```

The outer microchunk flow does not need to change.

## Open Validation Questions

1. Does the current difficult-stream profile preserve behavior on easier batch-mode datasets?
2. Are pass-1 chunks too small for reliable two-cluster k-medoids decisions, or does the existing `nfit < 4` and low-separation fallback cover this?
3. Should refchunk and match chunks use the same profile, or should refchunk be slightly less permissive because it seeds all downstream matching?
4. Do we need to retain scalar population floors in pass-1 as hard pre-filters, or does the feature-vector population weight handle this safely?
5. How much diagnostic output is acceptable in long streaming runs?

## Recommended Next Step

After the easy batch-mode validation finishes, integrate `simple_cavg_quality` into `microchunked2D%reject_cavgs` behind an internal backend constant and run the same difficult stream partitions again. If the state propagation and sentinel behavior remain unchanged, this should be a local refactor rather than a scheduler rewrite.

The safest initial target is refchunk and match rejection, because those tiers most closely resemble the real-world validation sets. Pass-1/pass-2 can follow once we know the quality-vector behavior is stable on smaller chunks.
