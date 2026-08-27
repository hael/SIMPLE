# Particle Picking

## Problem

Convert a corrected micrograph into particle-center coordinates while keeping
candidate generation, scoring, spatial suppression, and extraction as explicit
stages. SIMPLE provides segmentation-based picking for reference-free opening
and reference-based picking once representative particle images exist.

## Segmentation-based opening

The segmentation picker builds a filtered micrograph representation at the
requested particle scale, separates particle-like foreground from background,
labels connected regions, and reduces accepted regions to center coordinates.
Diameter policy and segmentation mechanics are separate: the former selects
the physical scale, while the latter owns the image-to-components decision.

This route is used early in preprocessing and streaming because it does not
require class-average references. Its coordinates seed the initial 2D analysis
that later produces stronger picking references.

## Reference bank

Reference picking prepares every supplied class average under the selected
filter, mask, contrast, rotation, and mirror policy. Let a prepared reference
`S_r` and a candidate micrograph window `T_c` contain `N` pixels. Invalid
references and windows with zero or numerically negligible variance are not
eligible.

## Complete-square Pearson score

For every candidate coordinate and every valid reference, SIMPLE evaluates the
same exhaustive Pearson correlation:

```text
mu_c        = sum(T_c) / N
var_c       = sum(T_c*T_c) / N - mu_c^2
score(r,c)  = [S_r dot T_c - mu_c*sum(S_r)] / [N*sqrt(var_c)].
```

References are normalized, but their residual single-precision sum is retained
in the numerator; assuming it is exactly zero changes scores at the rounding
level.

The production evaluator batches candidate windows into columns, obtains
window sums and squared sums from double-precision summed-area tables, and
computes all reference/window dot products with SGEMM. This is an evaluation
acceleration only: it does not prune references, coordinates, rotations, or
mirrors.

## Coarse and fine search

The coarse pass scores every valid iterator position and retains an
unthresholded score grid for peak and background statistics. Peaks passing the
selection policy enter fine refinement.

For each coarse peak, fine refinement enumerates the same clipped local
neighborhood at unit stride and scores every coordinate against every valid
reference. A strict-greater comparison retains the first maximum in traversal
order. If two refined peaks map to the same fine-grid cell, the established
last-update rule is preserved.

## Spatial competition

Candidate peaks are ordered by score. Greedy distance filtering accepts a peak
only when it is sufficiently far from already accepted higher-scoring peaks.
The final coordinate set therefore represents score maxima under a physical
exclusion radius, not an independent threshold at every pixel.

ROI/background policy and peak thresholds are applied without changing the
definition of the Pearson score. Accepted coordinates feed the particle
extractor, `.box` output, project metadata, thumbnails, distributed batches,
and stream stages.

## Implementation

- Segmentation: `src/main/pick/simple_pickseg.f90` and
  `src/main/pick/simple_picksegdiam.f90`.
- Reference policy and coordinate selection:
  `src/main/pick/simple_pickref.f90`.
- Bounded Pearson engine: `src/main/pick/simple_pickref_corr_batch.f90`.
- Workflow construction: `src/main/pick/simple_picker_utils.f90` and
  `src/main/simple_particle_extractor.f90`.
- Accepted acceleration contract:
  `doc/implementation_notes/reference_picker_flcf_plan.md`.

