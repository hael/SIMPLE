# Distance-Transform Shape Metrics for Class-Average Rejection

## Status

This is a design and implementation-plan note for adding rotationally invariant
shape evidence to SIMPLE class-average rejection. It does not specify an
immediate code change.

The natural integration point is the existing `model_cavgs_rejection` backend in
`src/main/cavg_quality`, not the streaming `simple_cluster2D_rejector` path. The
current model backend already owns class-average feature extraction, robust
normalization, learned feature weights, thresholding, analysis output, learning,
and promotion.

## Current SIMPLE Context

`model_cavgs_rejection` currently extracts a 12-feature bank:

- class population and resolution metadata;
- Otsu foreground centering and connected-component geometry;
- foreground/background local variance;
- optional stored class score evidence;
- center/edge signal and presence evidence;
- medium-frequency interior/boundary texture evidence.

The existing foreground path already performs the most important preprocessing
for distance-transform shape features:

- low-pass filter the class average;
- zero the edge average;
- threshold with Otsu;
- identify connected components;
- discard image-spanning components;
- keep the main foreground component;
- measure centroid displacement and foreground area relative to the expected
  mask disc.

That means the distance-transform proposal should reuse the current segmentation
path rather than introduce a parallel thresholding policy.

Relevant implementation anchors:

- `src/main/cavg_quality/simple_cavg_quality_feats.f90`: feature definitions,
  current foreground segmentation, hard rejects, and normalization.
- `src/main/cavg_quality/simple_cavg_quality_types.f90`: fixed feature count and
  shared model/result types.
- `src/main/cavg_quality/simple_cavg_quality_model.f90`: built-in models, model
  file I/O, feature weights, scoring, clustering, and thresholding.
- `src/main/cavg_quality/simple_cavg_quality_analysis.f90`: feature tables,
  analysis report, feature AUC, and threshold scans.
- `src/main/cavg_quality/simple_cavg_quality_learn.f90`: feature-policy search
  and learned-model generation.
- `src/main/image/simple_image_bin.f90`: connected components, binary
  morphology, border masks, and hole/background helpers.
- `src/main/commanders/simple/simple_commanders_reproject.f90` and
  `src/main/volume/simple_volinterp.f90`: existing volume reprojection path for
  building reference projection stacks.

## Review of the Idea

A binary foreground mask from a high-SNR 2D class
average is often stable enough to support shape descriptors, and the interior
distance transform is a useful way to summarize thickness, support, and gross
geometry without in-plane alignment.

The distance-transform histogram is the best first descriptor. If the binary
mask is rotated without changing its foreground support, the count of pixels at
each distance from the boundary is unchanged except for pixel-grid and
interpolation effects. This makes the histogram naturally invariant to in-plane
rotation and translation. Normalizing distances by the mask radius, expected
molecular diameter, or pixel size also makes it robust to box-size and sampling
differences.

The proposed complementary metrics are also appropriate:

- PCA eigenvalue ratio captures elongation without depending on orientation.
- Compactness or circularity captures ragged boundaries and aggregates.
- Hole count or hole area captures simple 2D topology.
- Distance-transform moments and quantiles summarize thickness distribution.

There are several important caveats:

- A distance-transform histogram is not a unique fingerprint. Different shapes
  can share similar thickness distributions, so it should be treated as evidence
  rather than as a standalone identity test.
- The histogram alone does not encode topology. Holes, component count, and
  foreground/background connectivity should be measured separately.
- A min/max envelope over reference projections is usually too permissive and
  too sensitive to outlier projections. Robust percentiles, nearest-reference
  distances, or a low-dimensional reference distribution are safer.
- The proposed reference-generation workflow is not reference-free. It is a
  reference-calibrated geometric filter. SIMPLE can still support a
  reference-free shape feature family, but projection-library rejection is a
  separate mode.
- The current learned model assumes feature weights are non-negative and that
  larger normalized values mean better class averages. Raw descriptors like
  aspect ratio or maximum thickness are not universally "higher is better".
  They must either be converted into quality-oriented scores or kept out of
  default models until reference calibration or manual learning validates them.
- Segmentation is the largest practical failure mode. The feature should inherit
  the existing low-pass/Otsu/connected-component path and expose diagnostics
  before it becomes a hard rejection criterion.

## Recommended Direction

Implement this in two stages.

Stage 1 should add reference-free shape descriptors as diagnostics and an
optional `shape` evidence family in `model_cavgs_rejection`. The default built-in
models should keep zero weight on new shape features until manual
`quality_mode=analyze` runs show that the new family improves false-positive or
false-negative behavior.

Stage 2 should add a reference-calibrated shape score. Users can generate
reference projections with the existing `reproject` program, then pass the
projection stack to `model_cavgs_rejection`. The rejection backend can compute
the same distance-transform signatures for reference projections and class
averages, then add one or more "closeness to reference envelope" scalar features.

Do not wire this into the live microchunk stream rejector first. The stream path
has stricter lifecycle and particle-state propagation requirements. The safer
path is to validate shape evidence through `quality_mode=analyze`, learn a model,
promote a built-in only after validation, and then consider stream integration.

## Proposed Shape Signature

### Mask Preparation

For each class average or reference projection:

1. Copy the image into a binary-image work object.
2. Zero the edge average.
3. Low-pass filter with the same foreground segmentation cutoff currently used
   by `model_cavgs_rejection`.
4. Apply Otsu thresholding.
5. Identify connected components.
6. Remove invalid full-image or image-spanning components.
7. Keep the largest valid foreground component.
8. Optionally fill small internal holes for thickness descriptors, while also
   retaining hole metrics from the unfilled mask.

This should be a shared helper used by both the existing foreground geometry
features and the new shape features, so future segmentation changes remain
consistent.

### Distance Transform

Compute the 2D Euclidean distance transform of the foreground mask: each
foreground pixel receives its distance to the nearest background pixel.
Background pixels are zero.

Prefer an exact separable squared Euclidean distance transform over a naive
all-pairs search or a chamfer approximation. Class-average boxes are modest, but
the feature may run over many class averages and many reference projections.

Distances should be represented in pixels internally, then normalized for
features by one of:

- expected mask radius in pixels;
- `mskdiam / smpd`;
- foreground equivalent-disc radius;
- explicit reference-stack scale when using calibrated mode.

### Histogram Features

Build a normalized distance-transform histogram over foreground pixels only.
Use a fixed number of bins and fixed normalized distance bounds so class averages
and reference projections are comparable.

Useful scalar summaries:

- normalized maximum distance, the approximate inradius of the thickest support;
- normalized mean and median distance;
- upper quantiles of distance, such as 75th and 90th percentiles;
- entropy of the distance histogram;
- coefficient of variation of nonzero distances;
- foreground area fraction relative to the expected mask disc.

The full histogram should be available to the reference-calibrated score, but
the learned quality model should consume a compact set of scalar scores unless a
deliberate high-dimensional feature-bank version is introduced.

### Shape and Topology Features

Add compact scalar descriptors alongside the distance-transform histogram:

- PCA aspect ratio from foreground pixel coordinates.
- Compactness, for example `4*pi*area/perimeter^2`.
- Boundary density or perimeter normalized by foreground area.
- Hole count and hole-area fraction, derived from connected components in the
  background.
- Number of foreground components before largest-component reduction.

Care is needed with directionality. Compactness can be "higher is better" for
many junk-rejection cases, but aspect ratio and thickness are often
particle-dependent. For reference-free learned models, prefer either diagnostic
output or quality-oriented transforms such as "typicality relative to the
dataset" rather than raw values that may be good for one target and bad for
another.

### Radial Distance-Transform Profile

The centroid-based radial profile should be optional. It is useful for roughly
centered, symmetric, spherical, or toroidal targets, but less generally robust
for flexible, elongated, or strongly view-dependent particles. It should not be
part of the default shape family unless validation shows it helps.

## Reference-Calibrated Mode

Reference-calibrated rejection should compare each class average against a
library of projected reference signatures.

Reference-library generation can reuse the existing `reproject` command rather
than projecting volumes inside `model_cavgs_rejection` at first. A user would
produce a stack of model projections over a sufficiently dense angular grid,
using the appropriate symmetry, mask diameter, sampling distance, contrast, and
box size. The rejection backend would then read that stack and compute the same
shape signatures.

Recommended comparison strategy:

- Compute a normalized distance-transform histogram for every reference
  projection.
- Compute PCA aspect ratio, compactness, area fraction, and topology metrics for
  every reference projection.
- For each class average, compute the nearest or low-percentile distance to the
  reference library rather than applying independent min/max bounds to every
  feature.
- Use robust reference statistics: median, MAD, percentile envelopes, or
  nearest-reference distances.
- Convert distances into quality-oriented scalar features, for example
  `shape_ref_score = -robust_z(distance_to_reference_library)`, so larger values
  mean closer to expected model geometry.

Good histogram distances for this use case:

- 1D Wasserstein distance, implemented as the sum of absolute cumulative
  histogram differences.
- Chi-squared distance with a small epsilon for empty bins.
- Hellinger or Jensen-Shannon distance, already conceptually aligned with
  normalized histograms.

For shape vectors, avoid raw min/max gates as the primary decision. A
Mahalanobis distance in a low-dimensional descriptor space, a robust nearest
neighbor distance, or a percentile envelope over the final scalar distance is
preferable.

## Implementation Plan

### Phase 0: Documentation and Validation Targets

Define the exact descriptors, feature names, output tables, and validation
criteria before code changes.

Acceptance criteria for this phase:

- clear distinction between reference-free shape descriptors and
  reference-calibrated shape scores;
- agreement that shape evidence starts as analysis/learn evidence, not as a hard
  reject;
- decision on whether initial reference support accepts a projection stack only
  or also accepts a volume directly;
- expected output files and analysis columns named.

### Phase 1: Shape Helper Design

Introduce a small shape-analysis helper owned by `src/main/cavg_quality`.

Planned responsibilities:

- build or receive the cleaned foreground mask;
- compute foreground area, centroid, component count, and hole metrics;
- compute the exact 2D Euclidean distance transform;
- build a normalized distance-transform histogram;
- compute scalar shape descriptors;
- expose reference-library comparison helpers.

Keeping this helper under `src/main/cavg_quality` is lower risk than adding a
general image API immediately. If other subsystems later need distance
transforms, the helper can be promoted into `src/main/image`.

### Phase 2: Feature-Bank Integration

Add a `shape` evidence family to the quality feature bank.

Required updates:

- increase the feature count in `simple_cavg_quality_types.f90`;
- add feature indices and definitions in `simple_cavg_quality_feats.f90`;
- compute shape features from the cleaned foreground mask;
- update built-in weight arrays in `simple_cavg_quality_model.f90`;
- preserve existing built-in behavior by assigning zero shape weights to current
  defaults;
- bump model-file and analysis-table versions;
- support reading older model files by padding missing new weights with zero, or
  explicitly require retraining for shape-aware models;
- update `simple_cavg_quality_analysis.f90` so raw and normalized shape columns
  appear in analysis and feature tables;
- update `simple_cavg_quality_learn.f90` so feature-policy search can include
  shape-aware policies.

Feature-policy options should stay small. A reasonable first set is:

- current policies unchanged;
- `microchunk_plus_shape`;
- `microchunk_plus_score_shape`;
- `microchunk_plus_score_signal_shape`.

The default model should not silently start rejecting more classes because shape
features exist. Shape-aware rejection should require a learned model, a promoted
model, or an explicit reference-calibrated score.

### Phase 3: Reference Library Support

Add optional reference-stack support to `model_cavgs_rejection`.

Recommended first interface:

- user supplies a reference projection stack;
- backend computes reference shape signatures at apply/analyze time;
- backend emits a reference summary file;
- backend adds a scalar reference-closeness feature to the shape family.

Possible future interface:

- user supplies `vol1`, `nspace`, `pgrp`, and `mskdiam`;
- backend invokes or reuses the existing reprojection machinery internally.

The projection-stack interface is preferable for the first implementation
because it keeps volume projection, contrast handling, masking, and orientation
sampling outside the rejection backend.

### Phase 4: Diagnostics and Outputs

Analysis output should make shape behavior inspectable.

Recommended outputs:

- raw and normalized shape columns in the main analysis table;
- per-feature AUC, robust separation, current weight, and suggested weight, as
  already done for current features;
- optional `cavgs_shape_reference_summary.txt` with reference-library
  descriptor distributions, histogram-distance percentiles, and outlier counts;
- optional `cavgs_shape_features.txt` if the shape signature contains more detail
  than the compact model feature columns can show;
- selected/rejected stacks unchanged.

The report should explicitly distinguish:

- hard rejects caused by existing validity gates;
- soft rejects from the learned model;
- classes whose shape score is outside the reference envelope;
- classes where segmentation failed or produced degenerate shape descriptors.

### Phase 5: Validation

Validation should happen before any shape-aware built-in model is promoted.

Synthetic tests:

- small binary masks with known distance-transform values;
- discs, rectangles, rings, and two-component masks;
- invariance under 90-degree rotations and approximate stability under
  interpolated in-plane rotations;
- compactness and aspect-ratio sanity checks.

Reference-projection tests:

- use SIMPLE's existing `reproject` path to generate projections of a known
  volume;
- verify that projected references fall inside their own reference distribution;
- verify that obvious synthetic junk masks fall outside;
- test contrast and mask-radius sensitivity.

Manual-data validation:

- run `quality_mode=analyze` on existing manually selected projects;
- compare current defaults against shape-aware learned candidates;
- inspect false positives and false negatives with selected/rejected stacks;
- check whether shape evidence removes junk aggregates without rejecting
  legitimate rare projections;
- compare chunk and pool contexts separately.

Promotion criteria:

- shape-aware model improves balanced accuracy or false-positive rejection on
  multiple datasets;
- manually good hard rejects do not increase;
- gains survive leave-one-dataset-out diagnostics;
- behavior remains stable when no reference stack is supplied.

### Phase 6: Optional Stream Integration

Only after batch/analyze validation should stream integration be considered.

The stream path currently uses `simple_cluster2D_rejector` directly. A future
stream integration should decide whether to:

- replace the scalar rejector with `model_cavgs_rejection` decisions;
- call only the shape helper as an extra scalar rule;
- keep the existing stream rejector unchanged and reserve shape evidence for
  batch/pool cleanup.

This decision should account for stream sentinels, chunk lifecycle, selected and
rejected stack writing, particle-state propagation, and retry behavior.

## Open Questions

- Should the first implementation expose only a reference projection stack, or
  should it also accept a volume and generate projections internally?
- What distance normalization is most stable across projects: mask radius,
  equivalent-disc radius, or `mskdiam / smpd`?
- Should aspect ratio be a diagnostic-only descriptor unless a reference model is
  present?
- How should multi-state references be represented: one combined projection
  library, state-specific libraries, or best score over all states?
- How much hole filling should be applied before the distance transform, and
  should hole metrics be computed before or after cleanup?
- Should shape features be unavailable when the foreground mask area is too small
  but not already hard rejected?
- Are there target classes where compactness penalizes true projections, such as
  elongated complexes, filaments, or membrane-associated particles?

## Summary

Distance-transform shape evidence is worth pursuing for SIMPLE, but it should be
introduced as a measured extension of the existing class-average quality model.
The safest path is to reuse the current segmentation machinery, add
rotation-invariant shape diagnostics, validate them through `quality_mode=analyze`
and `quality_mode=learn`, then add reference-calibrated scoring using projection
stacks generated by the existing reprojection workflow.

The most important design rule is to avoid turning raw geometric descriptors
into hard rejection gates too early. Shape should first become inspectable,
learnable evidence; only validated reference-calibrated distances should become
strong rejection signals.
