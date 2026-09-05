# Particle Picking

## Problem

Locate the centers of particles in a corrected micrograph. Two estimators are
used at different points in a project: a reference-free segmentation picker,
which needs only an approximate particle diameter and is used to bootstrap,
and a reference-based picker, which scores every position against a bank of
class averages once such averages exist. Both reduce to the same three steps:
score candidate positions, suppress non-maxima under a physical exclusion
radius, and threshold.

## Segmentation picking

The micrograph is low-passed to the requested resolution, denoised by a
total-variation filter, and binarized. The default threshold is Otsu's, the
threshold that maximizes the between-class variance of the intensity
histogram, applied in its tight variant; when a window size is supplied, a
Sauvola local-contrast threshold is used instead so that uneven ice does not
bias the segmentation. Two erosion passes remove thin bridges. Connected
components are labeled and filtered by size: either within
`mean +/- ndev * sd` of the component-size distribution (default 2.5 standard
deviations), or, when a molecular diameter is given, within a factor of four
of the expected area `pi (d/2)^2`. The center of mass of each surviving
component is the pick.

**Diameter estimation.** The streaming variant of this picker first
estimates the particle diameter across many micrographs. Each micrograph is
shrunk to 4 A per pixel, background-subtracted at 400 A, masked for carbon
and ice, denoised, and binarized at the intensity quantile that leaves 17
percent of the field as foreground. Component diameters are binned; bin
means are z-scored against the global median with a MAD scale and accepted
below three sigma. The consensus diameter then sets the box size and the
reference bank's scale.

## Reference picking

**Reference bank.** Each class average is automasked, damped in its
negative values, centered, and low-passed at `min(max(30, 0.15 d), 15)` A
for maximum diameter `d`. The bank is expanded over `nrots` in-plane rotations
(about 100 references in total by default, 12 rotations in streaming) and,
optionally, their mirrors. Rotations and mirrors are enumerated rather than
searched because a rotation of the template is a different template; nothing
in the score is rotation-invariant.

**Score.** For a prepared, normalized reference `S_r` with `N` pixels and a
micrograph window `T_c` at position `c`, the score is the Pearson
correlation, written so that the window statistics factor out:

```text
mu_c       = sum(T_c) / N,
var_c      = sum(T_c^2) / N - mu_c^2,
score(r,c) = [ S_r . T_c - mu_c sum(S_r) ] / [ N sqrt(var_c) ],
score(c)   = max_r score(r, c).
```

The window mean and variance come from summed-area tables in double
precision, and the reference-window dot products for a batch of positions
are one matrix product. This is an exact evaluation of the same correlation
at every position; nothing is pruned.

**Two-resolution search.** A coarse pass scores positions on a 3-pixel grid
at 4 A per pixel; a fine pass rescores the neighborhood of each coarse peak
at unit stride and 2 A per pixel. Ties are broken by first occurrence in
traversal order.

**Peak selection.** Peaks are selected in four stages:

1. keep the top-scoring positions up to a count cap equal to the number of
   grid cells at twice the stride;
2. greedy non-maximum suppression: walk positions in descending score and
   discard any within `d/3` of an already accepted peak. The accepted set is
   the set of score maxima under a physical exclusion radius, not an
   independent threshold at every pixel;
3. Otsu threshold on the surviving scores;
4. quantize the scores into five levels by sorted means and shift the Otsu
   level down by one, two, or three levels according to the requested
   particle density (low, optimal, high), then optionally cap the count.

**Background diagnostic.** Positions farther than `d/2` from every accepted
peak form a background set. The standardized mean difference and a
Kolmogorov-Smirnov test between peak and background scores are reported, and
a warning is issued when the two distributions are not separable
(`smd < 0.2` and `p > 0.5`), which indicates references that do not match
the data or a micrograph without particles.

## Rationale

- Segmentation needs no prior model and is robust to the absence of
  references, but its centers are biased by contrast and shape. It is used
  where nothing better exists and its picks are refined by 2D classification
  into references for the second picker.
- Pearson correlation normalizes each window by its own variance, so thick
  ice and carbon edges do not produce high scores merely by having high
  contrast. The mean subtraction and normalization are what distinguish it
  from a plain matched filter.
- Non-maximum suppression at one third of the diameter is what makes the
  output a set of particle centers rather than a correlation map.

## Implementation

- Segmentation: `src/main/pick/simple_pickseg.f90`, `src/main/pick/simple_picksegdiam.f90`;
  thresholds in `src/utils/simple_segmentation.f90`.
- Reference bank and coordinate selection: `src/main/pick/simple_pickref.f90`,
  `src/main/strategies/parallelization/simple_pick_strategy.f90`.
- Batched Pearson evaluation: `src/main/pick/simple_pickref_corr_batch.f90`.
- Workflow and extraction: `src/main/pick/simple_picker_utils.f90`,
  `src/main/simple_particle_extractor.f90`.
- Design note: `doc/implementation_notes/reference_picker_flcf_plan.md`.
