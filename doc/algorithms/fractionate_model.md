# Fractionate Model: Particle Re-extraction

## Problem

After particle alignment, the extraction coordinates can be improved using
the translations already estimated for each particle. Returning to the raw
movie also makes it possible to produce separate particle datasets for short
intervals of the exposure, while correcting the motion within each interval.

`fractionate_model` performs this re-extraction using a saved binary motion
model, then reconstructs each interval from saved 3D poses and exports raw
FSC diagnostics by default. With `oritype=ptcl2D`, it instead generates class
averages from saved 2D alignments and classes. Its polynomial refit uses the
existing patch measurements.
The application does not compare movie particles with reference images,
estimate new trajectories, or write a newly refined motion model.

This is distinct from [reference-free motion correction](motion_correction.md),
which estimates the initial global and local motion. The
[workflow policy](../policies/fractionate_model_policy.md) owns the current
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

Both particle fields are read independently and keep their own metadata,
with congruent row counts. In 3D mode retain particles active in both fields;
in 2D mode use the `ptcl2D` selection. Select the existing particle translation
`s` from the requested orientation segment
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
Once every interval's maps (`ptcl3D`) or class averages (`ptcl2D`) are made,
the distributed master removes these folders with their stacks and projects
(`remove_fraction_stacks`); a failed run keeps them. The copies in
`rec_frames_*` or `cavgs_frames_*` then point at stacks that no longer exist.

### Reconstruction after extraction

This branch applies to the default `oritype=ptcl3D`.

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
supported cross-dispatch mode. Nonuniform filtering and automatic envelope
masking are disabled. Classical postprocessing is enabled with `bfac=0` and
`mirr=no`: `_pproc.mrc` receives FSC amplitude weighting and a Butterworth
low-pass at FSC=0.143, with no B-factor fitting or sharpening. The native
`_lp.mrc` companion is retained and no mirrored map is written. Ordinary
gridding correction applies to every map; spherical support applies to standard outputs after writing the unmasked
`_unfil` half maps. Each interval is reconstructed by a nested,
complete run of the `reconstruct3D` program on its own project copy, with its
own freshly parsed parameters; the master then registers the delivered merged
map and FSC in that copy.

Before moving to the next interval, the master removes only the reconstruction's
exact `recvol_stateNN_part...` and `rho_recvol_stateNN_part...` even/odd files,
after validating every populated state's map/FSC and writing output registration.
Shared filename helpers identify the targets; canonical directory and file
checks precede deletion. Failed reconstructions retain their partials. Fraction
stacks and projects in `frames_a-b/` are removed after the last interval (see
Outputs). Final maps, half maps, FSCs and the copied reconstruction project
remain.

The master establishes missing even/odd labels once before extraction jobs
start. Interval reconstruction verifies that selected particle rows and
half-set labels agree across all projects. Selection follows `sample4rec`:
active particles with `updatecnt>0` when any active particle was updated,
otherwise all active particles. The master's working directory and thread
settings are restored after each run; its own parameters are never edited.
`nparts=1` reconstructs in the master process; larger values run the
program's own distributed gridding master.

## Dose-series export

The `_even_unfil.mrc` and `_odd_unfil.mrc` maps are exported after density and
interpolation correction, before spherical support. They are unmasked,
unregularized real-space float32 volumes. Standard half maps and the merged
volume retain their spherical mask, so the existing FSC surface describes the
same standard half-map pair.

After all intervals complete, usable state-1 dose metadata produce a
`dose_series.txt` manifest at the run root. Each row names the unmasked even/odd
maps using paths relative to that file and records the original inclusive
fraction range plus its accumulated physical dose bounds. The header includes
native geometry, saved voltage, dose per fraction, contributor/half counts and
a fingerprint of the fixed state-1 partition. The producer tag is
`# program fractionate_model`; manifest readers require that program name.

For fraction dose `d` from the saved motion model and interval `[a,b]`, the
bounds are `(a-1)*d` and `b*d` e/Angstrom^2. This dose already includes EER
fraction grouping; it is not multiplied by that grouping again. Starting after
fraction 1 retains the preceding dose. Only micrographs contributing to the
actual state-1 reconstruction sample enter the compatibility checks.

The shared manifest class checks metadata and real MRC headers. Publication
uses a temporary file and final rename after the full series succeeds. Missing
or incompatible dose metadata prevents manifest publication, with a reason in
the log, while leaving successful reconstructions and FSC outputs available.
The consumer `dose_filter3D` applies its own statistical mask, prunes the
predecessor of an overlapping final interval and requires at least two disjoint
windows. Fitting is the separate application [`dose_filter3D`](dose_filter3D.md).

## Raw fraction FSC diagnostics

Before the reconstruction loop, the strategy allocates an FSC matrix indexed
by native Fourier shell, actual fraction interval and state. Each delivered
binary FSC is copied directly into its interval/state column after checking
its shell count and finite samples. Unpopulated-state curves remain `NaN`.
For native box size B, sampling s in Å/pixel and shell j=1,...,B/2,

```text
resolution_j = B*s/j
reciprocal_frequency_j = j/(B*s).
```

After all intervals finish, each state gets `fsc_fractions_stateNN.tsv` with
columns `shell`, `resolution_A`, `reciprocal_frequency_A^-1` and one
`FSC_fractions_<first>-<last>` column per actual inclusive interval. The decimal
precision preserves the raw binary single-precision FSC values on round trip.

Fortran-controlled gnuplot renders two panels: all raw interval FSC curves and
a 3D surface of the same raw matrix. The surface uses reciprocal frequency,
fraction-interval midpoint index `(first+last)/2`, and FSC as its x/y/z
coordinates. Exact fraction-interval labels preserve the correspondence,
including an overlapping final interval. The y coordinate needs no physical
dose or voltage metadata.

Outputs share the `fsc_fractions_stateNN` prefix: `_diagnostics.png`,
`.gp`, and `_surface.tsv`, alongside the matrix TSV. No PDF is generated.
Negative values, zero and endpoints are displayed without clipping; missing
samples stay missing, with no interpolation across them. Sparse and empty
states are displayed without inventing a surface. Gnuplot runs headlessly;
plotting failures warn and retain the numerical export and reconstructions.
There is no Python rendering dependency.

For `oritype=ptcl3D`, a second diagnostic compares the standard merged
`recvol_stateNN.mrc` maps between fraction intervals, within each state.
Each unordered pair is evaluated once on the native shell grid; no alignment,
additional masking or filtering is applied. These are merged-map similarity
curves, distinct from the even/odd FSCs above. They do not use `_pproc`, `_lp`
or half maps.

Every `rec_frames_aa-bb/` folder receives
`fsc_cross_fractions_stateNN.tsv`, `_diagnostics.png`, `.gp` and `_surface.tsv`
(the latter three share the `fsc_cross_fractions_stateNN` prefix). The TSV
contains shell, resolution and reciprocal frequency, followed by one labelled
column per other target interval in its original order; self-comparisons are
excluded. The existing plotting routine draws the curves on the left and the
frequency/target-interval/FSC surface on the right, using the same style as
`fsc_fractions_stateNN_diagnostics.png`. Its title identifies the local
reference interval and state. Missing maps remain `NaN`, sparse comparisons
are shown as points and empty comparisons are explicitly labelled. Plotting
failures retain the numerical files. The 2D branch does not run this step.

The 3D branch delivers the reconstruction maps, unsharpened FSC-weighted
postprocessed maps, raw half-map FSCs and cross-fraction merged-map FSCs. No SNR
transformation, fit selection, dose-model fitting or model-derived outputs are
performed. The raw diagnostics impose no dose-metadata or fitting-coverage
requirements. Standalone Python analysis scripts remain independent.

## Accumulated interval maps

The 3D execution also writes one set of totals per populated state in
`rec_frames_total/`. It sums every saved unmasked/unfiltered even interval map
and every odd interval map separately, without averaging or dose weighting.
Overlapping final intervals remain in the sum exactly as saved. These are
restored-map sums, not a joint reconstruction from unique frames.

The outputs are `recvol_stateNN_even_unfil.mrc`,
`recvol_stateNN_odd_unfil.mrc`, and their literal sum
`recvol_stateNN_unfil.mrc`. All three remain unfiltered and unmasked.
`image_msk%automask3D` generates `automask3D_stateNN.mrc` from a copy of that
merged sum using non-tight Otsu, largest-component selection and fixed
`amsklp=15` Å, `binwidth=1`, `edge=1`. At this edge width the saved mask is
binary on the native grid. Caller mask settings are preserved.

The same binary mask supplies `distance_to_surface.mrc` for state 1 and
`distance_to_surface_stateNN.mrc` for other populated states. These volumes
have the mask's box and sampling. Each voxel holds the Euclidean distance in
Å to the nearest opposite-valued voxel centre, positive where the mask is 1
and negative where it is 0. The calculation subtracts the existing distance
transform of the complementary mask from that of the mask; it does not wrap
at box edges or apply a half-voxel surface offset. Empty/full masks are rejected.

Copies of the even and odd totals receive the ordinary soft spherical mask
from `mskdiam` before FSC calculation. `fsc_stateNN.bin` and
`fsc_stateNN.tsv` retain that spherical-masked curve, without density masking
or phase randomization. The merged sum receives the usual FSC amplitude
weighting and Butterworth low-pass at FSC=0.143; its
`recvol_stateNN_pproc.mrc` output retains both the spherical support and the
density envelope. Sharpening and mirrored outputs are disabled. The 2D branch
does not generate totals.

## Exporting voxel features

The standalone [exporter](../../scripts/export_fraction_voxel_features.py)
writes a voxel-by-feature TSV for one state:

```sh
python3 scripts/export_fraction_voxel_features.py /path/to/master --mskdiam 100 -o voxel_features.tsv
```

It requires NumPy, SciPy and mrcfile. Use `--state NN` for another state,
`--threads N` for FFT threads, `--chunk-rows N` to size TSV row batches, and
`--tmpdir /scratch` to choose temporary storage. The default output is
`MASTER/voxel_features_stateNN.tsv`. A companion `<output>.json` records
geometry, units, source files, interval ordering and filtering settings.

Rows contain every voxel centre inside the hard sphere of the supplied
`mskdiam` in Å, centred on SIMPLE's origin voxel `(N//2,N//2,N//2)` in zero-based
coordinates. This selection includes both density-envelope and solvent voxels.
Rows run in Z/Y/X order, with X fastest. The first five columns are
`x_vox`, `y_vox`, `z_vox`, `distance_to_surface_vox`, and `in_envelope` (0/1).
Distance values are read from the total's Å-valued distance volume and divided
by the isotropic sampling to obtain voxel units.

Four intensity groups follow: all unfiltered even fractions, all unfiltered
odd fractions, all FSC-filtered even fractions, then all FSC-filtered odd
fractions. Each group is ordered by increasing interval end fraction, hence
increasing cumulative exposure; overlapping intervals are retained. Column
names identify the group and exact interval. There are `5 + 4*N` columns for
`N` intervals. Absolute dose metadata are not required for this ordering.

The unfiltered columns use the saved `_even_unfil.mrc` and `_odd_unfil.mrc`
files. The native application filters only the merged map, so the script
generates filtered halves independently using each interval's
`fsc_stateNN.bin`. Each raw half receives the native soft sphere before FSC
amplitude weighting and order-eight Butterworth low-pass at FSC=0.143, with
no sharpening or density-envelope masking. A half is never replaced by the
merged `_pproc` map. Temporary filtered maps are deleted after export; no
input is modified. Memory mapping and row batches avoid keeping the full
feature matrix or every volume in RAM. Geometry mismatches, missing inputs,
nonbinary envelopes and invalid values fail without publishing a partial TSV.

## Per-voxel even/odd analysis

```sh
python3 scripts/analyze_fraction_voxel_features.py voxel_features.tsv -o voxel_stats
```

The analyser reads the TSV and its `.tsv.json` geometry sidecar. Unfiltered
analysis and output are temporarily suspended by commenting out the unfiltered
entry in `GROUPS`. Unfiltered columns are neither required nor parsed, and old
unfiltered output files are left untouched. The exporter itself is unchanged.

Active outputs retain the filtered Pearson, Spearman and `even + odd` population
variance (`ddof=0`) maps, histogram figure, three-panel distance-scatter figure,
and progress JSON. New maps are `voxel_stats_filtered_pearson_pvalue.mrc` and
`voxel_stats_filtered_spearman_pvalue.mrc`. All maps retain full original XYZ
geometry. Correlations/p-values for constant series are `NaN`; background voxels
are `-1` in both correlation and both p-value maps; background variance stays
zero. These background values are written from initialization onward. Variance units are squared input intensity.

Pearson p-values use the two-sided exact beta null for independent bivariate
normal observations, with p=1 for two nonconstant pairs. Spearman p-values use
the two-sided absolute-tail pairing permutation test with average tied ranks.
All N! pairings are enumerated for N<=8, including tied-value multiplicities;
above eight, use 9,999 fixed-seed Monte Carlo permutations and the conservative
(count+1)/(B+1) estimate. `--spearman-resamples B` changes the latter count.
Null distributions are cached by sorted rank tie patterns, and results do not
depend on chunk size. No asymptotic Spearman p-value is substituted.

These p-values are nominal and uncorrected for multiple voxel tests. Interval
observations need independence/exchangeability under their respective nulls.
Shared frames in overlapping intervals, temporal dependence, and shared
FSC-derived filtering can violate those assumptions. The JSON records methods,
permutation count/seed, assumptions and interval overlap. Exact enumeration
removes the small-sample approximation; it does not remove dependence bias.
See [Pearson](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.pearsonr.html)
and [permutation tests](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.permutation_test.html).

Filtered histogram panels retain separately normalized inside/outside-envelope
Pearson and Spearman distributions. The filtered scatter retains Pearson,
Spearman and filtered sum variance vs signed distance in voxel units (positive
inside, negative outside). Its shared fixed-seed uniform voxel sample is capped
at 200,000 (`--scatter-max-points`); maps and histograms use every table voxel.

Default analysis chunks remain 65,536 rows (`--chunk-rows`). All five active
MRCs, including p-values, write directly to their final paths and flush with
header statistics after every chunk. Both figures refresh after the first/last
chunks and periodically (`--plot-every`, default ten seconds). Progress JSON
updates per chunk; failed/interrupted outputs remain partial, and only
`complete` certifies completion. Source files are preserved. Other options
remain `--metadata`, `--bins`, and `--overwrite`; dependencies remain NumPy,
SciPy, mrcfile and Matplotlib.

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
  `new_fractionate_model`.
- Entry point: `src/main/commanders/simple/simple_commanders_motion.f90`,
  `exec_fractionate_model`.
- Recentering, frame-block endpoints, worker execution, and project assembly:
  `src/main/strategies/parallelization/simple_fractionate_model_strategy.f90`:
  `distr_execute` coordinates dispatch, `assemble_motion_fractions` merges and
  writes interval projects, then `reconstruct_motion_fractions` delivers maps
  or `generate_fractions_cavgs` invokes `commander_make_cavgs_distr`.
- Model loading, detector correction, per-frame extraction, summation, and
  normalization: `src/main/preprocess/simple_particle_extractor.f90`.
  `generate_cumulative_particles_stack` retains its historical name but sums
  only the supplied inclusive interval.
- Binary persistence and SVD refitting:
  `src/main/motion/simple_motion_model.f90`.
- Shared polynomial coordinate mapping, basis, and evaluation:
  `src/main/motion/simple_motion_correct_utils.f90`.
- Raw FSC matrix export and native gnuplot diagnostics:
  `src/main/strategies/parallelization/simple_fractionate_model_strategy.f90`.
  The reconstruction coordinator keeps the collected matrix after interval
  partial-volume cleanup and writes each state's tables and plots.
- Existing numerical coverage:
  `src/main/motion/simple_motion_model_tester.f90` checks binary round trips,
  accepted/rejected patch metadata, known polynomial coefficients at reference
  frames 1 and 4, and zero/nonzero fit RMSDs. It does not exercise the complete
  extraction and distributed assembly workflow.

## Class averages from 2D particle metadata

`oritype=ptcl2D` selects particles, coordinates, saved in-plane alignment,
classes, particle CTF values and halves from `ptcl2D`. Both particle segments
are read independently and retain their own records; their row counts must
agree. The selected 2D shift supplies integer recentering, which is subtracted
from both fields' own residual shifts. The existing project half-partition
checks remain in force.

After extracting the same equal-weight, non-dose-weighted intervals, the
master copies each assembled project into `cavgs_frames_aa-bb/` and invokes
`commander_make_cavgs_distr`, using the solve2D command-copy pattern. It passes
`mskdiam` and the original compute settings, and preserves supplied averaging
options such as `ncls`, `remap_cls`, `objfun`, `ml_reg` and `pgrp`. Omitted
options retain the commander's normal defaults. Class remapping or expansion,
when requested, is performed by that existing workflow.

The interval caller always sets `refs=fraction_cavgs_aa-bb.mrc`, overriding
only the parent's output filename. The interval directory therefore contains
that merged stack, its `_even` and `_odd` stacks, class metadata and class
FRCs. The fraction stacks are removed after the last interval. This branch does not produce volume
FSC diagnostics or a volume dose-series manifest; the default `ptcl3D`
reconstruction path keeps its existing behavior.
