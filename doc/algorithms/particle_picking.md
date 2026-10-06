# Particle Picking

## Problem

Locate the centers of particles in a motion-corrected micrograph: the
weighted sum of the gain-corrected, aligned movie frames (warped by the local
deformation model where it was accepted), dose-weighted when the total dose
is given. No CTF correction is applied before picking.

Two estimators are used at different points in a project: a reference-free
segmentation-based picker for bootstrapping, which needs only a broad range
of plausible diameters, and a reference-based picker, used once 2D references
exist. The segmentation-based picker finds particles as connected components
of a denoised, binarized micrograph. The reference-based picker scores every
candidate position against a bank of references, suppresses non-maxima
within a physical exclusion distance, and thresholds. Its references are
either class averages or reprojections of a 3D map; the map route is the one
the streaming pipeline uses, and it gives much better picks.

## Segmentation picking

The segmentation-based picker bootstraps the streaming pipeline, where it is
run on the first 1000 micrographs as they are preprocessed. It treats
particles as compact bright blobs in a heavily denoised micrograph and
decides from the data which blob sizes are particles.

**Preparation.** The slowly varying background is removed by a 400 A
high-pass, the contrast is set so that particles are bright, and the
micrograph is Fourier-cropped to 4 A per pixel. Amorphous carbon is then
located and excluded. The gradient magnitude of the micrograph, low-passed to
15 A, is histogrammed in overlapping 256 A patches on a 128 A grid, and the
patch histograms are grouped by k-means into five clusters. The clusters are
ordered by their spread about the mode, and the ordered list is split where
the pooled histograms of the two groups differ most. When that
total-variation distance exceeds 0.2, the group with the wider spread is
masked as carbon. A micrograph that is more than 98 percent carbon is not
picked.

**Cascade filter.** The micrograph is then denoised by a fixed cascade in
which each stage removes a different kind of noise. The mean of the edge
pixels is subtracted and negative values are damped tenfold, so that dark
features cannot compete with the bright particles. The micrograph is
low-passed to 15 A and smoothed by a B-spline smoother (smoothing weight 7);
this intermediate image is the denoised view shown to the user. Non-local
means then averages each pixel with pixels whose neighborhoods look alike,
which suppresses noise without blurring particle edges. Iterated conditional
modes (weight 100), a Markov-random-field denoiser, pulls the image toward
piecewise-constant regions, and a 3 x 3 median filter removes the isolated
pixels that remain.

**Binarization and components.** The filtered micrograph is thresholded at
the intensity that leaves 17 percent of its pixels as foreground. It is
eroded twice (4 A per pass) to cut thin bridges between neighbors, and its
holes are filled: every background region except the largest becomes
foreground. The carbon mask is applied and the connected components are
labeled. Each component's diameter is measured, with 8 A added back for the
erosion; components between 20 and 500 A are kept, and the centroid of each,
the mean position of its pixels, is a pick.

**Diameter bins and box size.** A first batch of micrographs (100 in
streaming) decides which picks are particles. Their diameters are sorted into
five bins with edges at 20, 100, 200, 300, 400, and 500 A. A bin is accepted
when its mean lies within three robust standard deviations of the batch
median, the robust standard deviation being the median absolute deviation
scaled to match a Gaussian. The box is the largest diameter in the accepted
bins times a factor that falls linearly from 1.5 for the smallest particles
to 1.0 for particles 400 pixels across, rounded to the nearest FFT-friendly
size. Picks outside the accepted bins are dropped, and the rest are written
with that one box. The bins and the box are decided once, by the first
batch, so every later batch is picked alike.

The stand-alone picker (`prg=pick picker=segdiam`) runs the same
preparation, cascade, and binarization, but keeps components between 12 A
and an upper bound `moldiam_max` instead of selecting bins.

## Reference picking

**Where the references come from.** There are two sources: a set of selected
class averages, or a 3D map reprojected along directions spread evenly over
the sphere. In streaming, the map is built automatically from the bootstrap
data. The segmentation picks are classified in 2D by [solve2D](solve2d.md),
the good class averages are selected, and an ab initio map with three states
is determined from those class averages alone. A state is eligible when it is
a single compact object, holds at least a tenth of the population, and is
resolved to within 1.5 times the best state's FSC resolution. Among the
eligible states, the one covering the most distinct views is kept and
reprojected along 50 directions. Outside streaming, any map can be
reprojected and its projections given as the references.

**Reference bank.** When there are more references than requested (`ncls`,
10 in streaming), they are first reduced to that many medoids by
class-average clustering. Each reference is then automasked, damped in its
negative values, centered, and low-passed at `0.15 d` clamped to 15 to 30 A,
for `d` the largest automasked reference diameter. The bank is expanded over
`nrots` in-plane rotations and, optionally, their mirrors; by default the
number of rotations is chosen to give about 100 references before mirroring,
and streaming uses 12 rotations with mirrors. Rotations and mirrors are
enumerated rather than searched because a rotation of the template is a
different template; nothing in the score is rotation-invariant. Each
prepared reference is normalized to zero mean and unit variance. The
micrograph itself is low-passed at a fixed 15 A.

**Region of interest.** By default, positions on crystalline ice and on
amorphous carbon are excluded before scoring. Ice is found in overlapping
128-pixel boxes resampled to about 1.7 A per pixel: a box is masked when the
power around its strongest Fourier component near 3.7 A, the crystalline-ice
reflection, exceeds five times the power typical of the 15 to 6 A band. The
test needs a pixel size of about 1.7 A or finer and is skipped otherwise.
Carbon is found by the same patch-histogram test as in segmentation picking.
A micrograph with less than 1 percent of its area left is not picked.

**Score.** For a prepared reference `S_r` with `N` pixels and a micrograph
window `T_c` at position `c`, the score is the Pearson correlation. Because
the reference has zero mean and unit variance, the score can be written so
that the window statistics factor out:

```text
mu_c       = sum(T_c) / N,
var_c      = sum(T_c^2) / N - mu_c^2,
score(r,c) = [ S_r . T_c - mu_c sum(S_r) ] / [ N sqrt(var_c) ],
score(c)   = max_r score(r, c).
```

**Two-resolution search.** A coarse pass scores every third position at 4 A
per pixel; a fine pass rescores the 13 x 13 neighborhood (169 positions) of
each coarse peak at unit stride and 2 A per pixel. Ties are broken by first
occurrence in traversal order.

**Batched evaluation with BLAS.** The formula separates the window
statistics from the dot products, and the two are computed differently. The
window mean and variance come from summed-area tables of the micrograph and
of its square, built once in double precision, so each costs four table
lookups whatever the box size. The dot products `S_r . T_c`, which are
nearly all of the work, are computed as matrix products, 256 positions at a
time. Each window is unrolled into a column of a matrix `W` (`N` rows, one
column per position), the whole bank, rotations and mirrors included, forms
the columns of a matrix `R`, and `G = R^T W` is computed by one
single-precision general matrix multiply (SGEMM) from the BLAS library SIMPLE
is linked against. Each entry `G(r, c)` is one dot product; the score formula
is applied entry by entry, and the maximum over references is taken in each
column.

Before the tables are built, the micrograph is mean-subtracted and scaled to
a largest absolute value of one, which keeps the single-precision products
well conditioned. The `mu_c sum(S_r)` term keeps the score an exact Pearson
correlation when a single-precision reference does not sum to exactly zero.
Every candidate position is scored against every reference and nothing is
pruned; the correlation is computed directly in real space, not through
Fourier transforms. On test micrographs the batched scores agreed with a
direct pixel-by-pixel evaluation to within `5e-6`.

**Peak selection.** On the coarse grid, peaks are selected in four stages:

1. keep the top-scoring positions, up to as many as there are grid cells at
   twice the stride;
2. greedy non-maximum suppression: walk the positions in descending score
   and discard any within the exclusion distance of an already accepted
   peak. The exclusion distance is one third of the reference box width
   unless a distance in Angstrom is given (`thres`);
3. Otsu threshold on the surviving scores;
4. quantize the scores into five levels by sorted means, snap the Otsu
   threshold to the nearest level boundary, and lower it by zero, one, or two
   levels for the requested particle density (low, optimal, high); then
   optionally cap the count.

The fine pass then moves each surviving peak to the best-scoring position in
its neighborhood, repeats the non-maximum suppression at the finer sampling,
and applies the optional count cap again.

**Background diagnostic.** On the coarse grid, positions farther than half
the reference box width from every accepted peak form a background set. The
standardized mean difference and a Kolmogorov-Smirnov test between peak and
background scores are reported, and a warning is issued when the two
distributions are not separable (`smd < 0.2` and `p > 0.5`), which indicates
references that do not match the data or a micrograph without particles.

## Rationale

- Segmentation needs no prior model and works without references, but its
  centers are less precise. A pick is the unweighted centroid of a
  thresholded blob. Where contrast varies, the threshold trims the blob
  unevenly, and the centroid of a silhouette is the particle's center only
  for compact, roughly symmetric views: for elongated or asymmetric
  particles the offset changes with the view. Touching particles merge into
  one blob, which is either dropped by the size window or picked at a point
  between them. Segmentation is therefore used where nothing better exists,
  and its picks are turned by 2D classification, and in streaming by a 3D
  map built from the class averages, into references for the second picker.
- Map reprojections make better references than the class averages they are
  built from. They cover the views evenly, including views too rare in the
  bootstrap data to form a clean class. They are less noisy than any single
  class average, because the map pools the signal of every class average
  that went into it, although they keep the map's own noise and errors.
  And every one is a projection of the same density, so as long as the
  chosen state is the particle, the bank includes no junk classes or classes
  of a contaminant.
- Pearson correlation normalizes each window by its own variance, so thick
  ice and carbon edges do not produce high scores merely by having high
  contrast. The mean subtraction and normalization are what distinguish it
  from a plain matched filter.
- The correlation is evaluated in real space with BLAS because the search
  needs only some of the shifts: every third position in the coarse pass and
  small neighborhoods of the coarse peaks in the fine pass. A Fourier
  cross-correlation would compute every shift for every reference. Gathering
  exactly the needed windows and letting one matrix multiply form all their
  dot products with the bank uses the processor's cache and vector units far
  better than a loop of single dot products, and the search stays exhaustive
  over its candidates.
- Non-maximum suppression within a physical exclusion distance is what makes
  the output a set of particle centers rather than a correlation map.

## Implementation

- Segmentation: `src/main/pick/simple_picksegdiam.f90`; carbon and ice
  detection, the cascade filter, and binarization in
  `src/main/image_processing/simple_micproc.f90`.
- Diameter bins and box size: `src/main/pick/simple_segdiam_bin_picker.f90`;
  the streaming preprocessing pass in `src/main/stream/simple_mini_stream_utils.f90`.
- Streaming references from a 3D map:
  `src/main/stream/stages/simple_stream_stage_initial_analysis.f90`; policy in
  `doc/policies/stream/reference_generation_policy.md`.
- Reference bank and coordinate selection: `src/main/pick/simple_pickref.f90`,
  `src/main/strategies/parallelization/simple_pick_strategy.f90`.
- Batched Pearson evaluation: `src/main/pick/simple_pickref_corr_batch.f90`;
  the SGEMM wrapper `gemm_tn` in `src/utils/math/simple_linalg.f90`.
- Workflow and extraction: `src/main/pick/simple_picker_utils.f90`,
  `src/main/preprocess/simple_particle_extractor.f90`.
- Design note: `doc/implementation_notes/completed/reference_picker_flcf_plan.md`.
