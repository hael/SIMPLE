# Motion Model Particle Re-extraction

## Problem

After particle alignment, the extraction coordinates can be improved using
the translations already estimated for each particle. Returning to the raw
movie also makes it possible to produce separate particle datasets for short
intervals of the exposure, while correcting the motion within each interval.

`refine_motion_model` performs this re-extraction using a saved binary motion
model, then reconstructs each interval from the saved particle poses. Despite
its name, it does not compare movie particles with reference images, estimate
new trajectories, or write a newly refined motion model. Its polynomial refit
uses the existing patch measurements, not new image-alignment measurements.

This is distinct from [reference-free motion correction](motion_correction.md),
which estimates the initial global and local motion. The
[workflow policy](../policies/refine_motion_model_policy.md) owns the current
execution, scientific, project and output contracts; this chapter explains the
algorithm.

## Model

Let `g(t)` be the stored global drift and `d_p(t)` the stored local displacement
at patch center `p`. All displacements are in pixels of the scaled movie.
The reference frame `r` is the saved model's fixed frame, not the start of
each output block. Global drift is re-referenced as

```text
G(t) = g(t) - g(r).
```

When the saved patch model was accepted and its measurements are available,
each component of local motion is fitted to `d_p(t) - d_p(r)` by unweighted
least squares, solved by SVD. The 18 basis functions are

```text
{1, x, x^2, y, y^2, xy} multiplied by {tau, tau^2, tau^3},
tau = t - r,
x = (X - 1)/(Nx - 1) - 1/2,
y = (Y - 1)/(Ny - 1) - 1/2.
```

`X,Y` are one-based pixel coordinates and `Nx,Ny` are the scaled movie
dimensions. Denote the fitted vector field by `L(X,Y,t)`. The time basis
makes `L(X,Y,r) = 0`. If patch correction was rejected or is unavailable,
`L = 0` and only global drift is applied. Refitting updates the in-memory
coefficients and per-axis fit RMSDs; it does not alter the measured offsets,
change the saved reference-frame metadata, or perform a new acceptance test.

The saved coefficients are not reused because they follow an export
convention (fitted relative to the first frame) rather than the reference
frame `r`. Refitting the stored patch measurements at `r` gives the same
unweighted fit that produced the micrograph. More generally, the persisted
model does not depend on a particular reference frame: the stage-drift and
local offsets can be re-referenced to any convention and the polynomial
refitted from them.

For a particle centered at `c`, the applied motion comes from
`G(t) + L(c,t)`. Local deformation is evaluated once at the particle center:
each particle frame receives a translation, not a pixelwise deformation
across the particle box.

## Algorithm

### Recenter the extraction boxes

Only particles active in both the 2D and 3D records are retained. Select the
existing particle translation `s` from the requested orientation segment
(3D by default), and round it to an integer vector `k`:

```text
k = nint(s)
new_corner = old_corner - k
new_shift_2D = old_shift_2D - k
new_shift_3D = old_shift_3D - k.
```

Both residual translations receive the same correction. If the box size
changes, add `(old_box - new_box)/2` to the new corner to preserve its center.
These coordinates are zero-based in the scaled-movie coordinate system.
By convention, particle shifts are also in scaled-movie pixels, so a
translation and a box corner can be combined without conversion.

### Prepare the movie

Read the saved model, its movie, and its gain reference when present. Decode
EER with the stored grouping and upsampling settings. Apply gain correction
and cure the stored detector-outlier coordinates, also including zero-valued
EER gain defects. Defect replacement uses neighboring non-outlier pixels,
with random replacement for insufficient neighborhoods. Fourier-crop the
frames to the stored scaled dimensions when necessary.

The model contains dose metadata and frame-quality weights, but this
re-extraction path currently does not apply dose weighting or use the
frame-quality weights in the block sums.

### Extract motion-corrected particle frames

For every retained particle and every selected frame, subtract the global
and local displacement from the padded extraction corner. Round to an
integer extraction position, window the frame, apply the remaining
subpixel translation as a Fourier phase shift, and clip to the output box.
The corrected particle frames are retained for reuse across output blocks.

### Sum frame blocks

The `f` in `fromf`, `tof` and `stepf` stands for fraction: indices count the
movie fractions stored in the model, so they have the same meaning for every
movie, EER included, where a fraction is a group of raw EER frames.
With `tof=0`, the master resolves a common endpoint from the first eligible
movie model. This selects all fractions only when the selected movies have the
same saved fraction count; see the policy's current limitations.

For inclusive bounds `a = fromf`, `b = tof` and block width `s = stepf`, make
`ceil((b-a+1)/s)` blocks. For block index `j`, starting at one,

```text
end_j   = min(a - 1 + j*s, b)
start_j = max(a, end_j - s + 1).
```

Thus `a=1, b=17, s=5` gives `1-5`, `6-10`, `11-15`, and `13-17`.
Blocks are not cumulative. A short final remainder is extended backward to
retain the full width, so the last block overlaps its predecessor. When the
entire selected interval is shorter than the width, it is used once.

For particle `i`, sum the corrected frames with equal weights:

```text
S_i,j = sum over t=start_j..end_j of corrected_particle_i,t.
```

Each block starts from zero; there is no accumulation from the preceding
block. Afterwards, invert contrast for black input particles, subtract a
background ramp, and normalize background noise using a circular mask of
radius `0.75 * box/2`. Background-ramp subtraction is unconditional in this
path.

### Outputs

The default block width is three saved movie fractions (`stepf=3`). Explicit
`stepf` values override this default.

Each block `start-end` gets a directory `frames_<start>-<end>/` holding one
particle stack per movie, `ptcls_from_<movie>_<start>_<end>`, and a project
named after the input project. That project points to the block's stacks,
carries the recentered box coordinates and 2D/3D shifts, and keeps the other
segments of the input project. Movies with no retained particles are left out.

### Reconstruction after extraction

The distributed master finishes writing every interval project before it
reconstructs the intervals, sequentially, inside `distr_execute`. Extraction
workers only write their extraction products. The existing `finalize_run`
remains completion bookkeeping.

Every `frames_<start>-<end>` project is copied into
`rec_frames_<start>-<end>` with endpoints padded to at least two digits, for
example `rec_frames_01-03`. The copy owns its project metadata and registered
maps and FSC curves. Existing state-map/FSC registrations are removed from the
copy before publishing interval products; the source interval project stays
unchanged. Reconstruction directories must not already exist.

Reconstruction uses gridding, `objfun=cc`, `ml_reg=no`, and
`cc_emit_sigma=no`. There is no sigma2 bootstrap, loading, whitening, or
emission. Maps use the saved 3D orientations, residual translations, CTF
metadata, the interval stacks' native geometry, and the supplied `pgrp` and
`mskdiam`. Additional undocumented reconstruction-crop arguments are not a
supported cross-dispatch mode. Optional filtering, automatic envelope masking,
and postprocessing are disabled; ordinary gridding correction and spherical
soft support still apply. Each interval is reconstructed by a nested,
complete run of the `reconstruct3D` program on its own project copy, with its
own freshly parsed parameters; the master then registers the delivered merged
map and FSC in that copy.

The master establishes missing even/odd labels once before extraction jobs
start. Interval reconstruction verifies that selected particle rows and
half-set labels agree across all projects. Selection follows `sample4rec`:
active particles with `updatecnt>0` when any active particle was updated,
otherwise all active particles. The master's working directory and thread
settings are restored after each run; its own parameters are never edited.
`nparts=1` reconstructs in the master process; larger values run the
program's own distributed gridding master.

## Rationale

Integer recentering incorporates the existing particle alignment without
discarding its subpixel residual. Using one temporal reference for every
block keeps all extracted datasets in a common motion coordinate system.
The polynomial transfers patch measurements smoothly to each particle center
without re-estimating motion from weak individual-frame particle images.

Fixed-width blocks normally contain equal numbers of frames. Anchoring the
last block at the final frame avoids a short, lower-exposure remainder, at
the cost of shared data with the preceding block. Those overlapping outputs
must not be treated as statistically independent observations or combined
as though they were disjoint portions of the exposure.

## Implementation

- CLI registration: `src/main/ui/simple/simple_ui_preproc.f90`,
  `new_refine_motion_model`.
- Entry point: `src/main/commanders/simple/simple_commanders_motion.f90`,
  `exec_refine_motion_model`.
- Recentering, frame-block endpoints, worker execution, and project assembly:
  `src/main/strategies/parallelization/simple_refine_motion_model_strategy.f90`:
  `distr_execute` coordinates dispatch, `assemble_motion_fractions` merges and
  writes interval projects, then `reconstruct_motion_fractions` delivers maps.
- Model loading, detector correction, per-frame extraction, summation, and
  normalization: `src/main/preprocess/simple_particle_extractor.f90`.
  `generate_cumulative_particles_stack` retains its historical name but sums
  only the supplied inclusive interval.
- Binary persistence and SVD refitting:
  `src/main/motion/simple_motion_model.f90`.
- Shared polynomial coordinate mapping, basis, and evaluation:
  `src/main/motion/simple_motion_correct_utils.f90`.
- Existing numerical coverage:
  `src/main/motion/simple_motion_model_tester.f90` checks binary round trips,
  accepted/rejected patch metadata, known polynomial coefficients at reference
  frames 1 and 4, and zero/nonzero fit RMSDs. It does not exercise the complete
  extraction and distributed assembly workflow.
