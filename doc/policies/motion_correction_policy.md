# Reference-free Motion Correction Policy

Status — 2026-10-07: In development. This is a working policy, not a declaration
that the implementation and all integration checks are complete.

## Purpose and scope

SIMPLE's reference-free motion correction estimates stage drift and, optionally,
local beam-induced motion from the frames of each movie. Alignment references
come from that movie, with the frame being aligned excluded from its own
reference. No reconstructed volume, particle orientations, or external particle
reference may enter this estimation path.

This policy governs `motion_correct` and the motion-correction stage of
`preprocess`, including shared-memory, distributed-worker, and streaming use of
that stage. It defines the scientific and data contracts that changes must
preserve. Current implementation limitations are identified separately below.

`refine_motion_model` is a separate, reference-based workflow. Particle
extraction from saved motion models consumes these outputs; it is not
reference-free motion estimation. These consumers must respect the producer's
metadata conventions without changing this workflow's objective.

The [algorithm description](../algorithms/motion_correction.md) explains the
alignment method. This policy owns its execution, representation, output, and
validation requirements. Gain-orientation analysis has its own
[policy](motion_gain_analysis_policy.md).

## Ownership and execution

| Owner | Responsibility |
|---|---|
| UI, parameters, commanders, and execution strategies | Register options, resolve defaults, prepare the gain, select movies, schedule work, and merge project records |
| `simple_motion_correct_iter` | Orchestrate one movie, choose patch acceptance or fallback, write products, and update its project row |
| `simple_motion_correct` | Prepare frames, coordinate global and local correction, generate weighted products, and manage movie-local state |
| `simple_motion_align_hybrid` | Estimate translations using discrete and continuous leave-one-out alignment |
| `simple_motion_patched` | Define patches, estimate residual local trajectories, fit the deformation field, and generate warped micrographs |
| `simple_motion_correct_utils` | Share gain handling, EER fraction calculation, spatial normalization, polynomial basis, and evaluation |
| `simple_motion_model` | Preserve motion metadata and serialize binary and STAR representations |

Scheduling must not change gain orientation, frame selection, motion conventions,
or output normalization. Strategies own orchestration; numerical algorithms stay
in the motion and image modules. Reuse the shared polynomial helpers rather than
introducing another basis or coordinate conversion in a consumer.

The motion-correction module holds movie state in module variables. Concurrent
movies must not share this state within one process; OpenMP parallelism inside a
movie is distinct from independent distributed worker processes. Complete the
per-movie cleanup before initializing the next movie, including fallback paths.

## Movie preparation

The required order is decoding or reading, gain correction, defect curation,
Fourier scaling, global alignment, and optional local alignment. Gain correction
and defect coordinates refer to the decoded movie grid before Fourier scaling.

- Use movie-row sampling and microscope metadata. Updating the output sampling
  must not overwrite the meaning of the original movie sampling.
- Alignment requires at least two frames, both before and after `max_dose`
  truncation. `motion_correct_init` stops with a hard error otherwise; invalid
  preparation must not proceed to spectrum generation or publish a movie record.
- Distinguish user-requested output scaling from the aligner's temporary Fourier
  crop for shift search. Return estimated shifts in the scaled movie's pixel
  units, not in the temporary search grid's units.
- Preserve both original and scaled dimensions and sampling distances. Dimension
  rounding and EER upsampling must remain consistent with the exported metadata.

### EER fractions and dose accounting

EER processing requires either an explicit `eer_fraction` or a supplied
`total_dose` from which grouping can be determined using `fraction_dose_target`.
Only complete groups are decoded; leftover raw EER frames are excluded and the
effective total dose is adjusted accordingly. The current supported upsampling
values are 1 and 2; the latter doubles each image dimension and halves the pixel
size.

After decoding, a frame index identifies a fraction, not a raw EER frame.
`eer_fraction` and `eer_upsampling` must travel with the movie model. Keep these
quantities distinct:

- `total_nframes`: available movie frames or decoded fractions before dose
  truncation;
- `nframes`: frames retained for alignment and integration;
- `dose_per_frame`: effective dose per retained frame or fraction;
- `target_dose_per_frame`: the grouping target, not a measured fraction dose;
- model `total_dose`: the supplied acquisition dose;
- model `accumulated_dose`: `nframes * dose_per_frame`.

Dose weighting is enabled by supplying `total_dose`. Supplying an EER grouping
alone does not enable it. `max_dose` acts on whole frames: the current selection
retains the first frame whose cumulative dose exceeds the limit, if there is
one. It is not a strict fractional-frame dose cap. Changes to this endpoint rule
must update dose metadata and endpoint tests together.

## Gain orientation and detector defects

### Apply the gain exactly once

The gain filename and `flipgain` form one operational pair. For a requested
flip, `flip_gain` materializes an oriented gain file and updates both the
caller's gain path and the command-line `gainref`. It does not reset either
flipping mode; the caller owns that reset.

The shared-memory and distributed initializers of `motion_correct` and
`preprocess` call this helper. When `gainref` and `flipgain` are defined, they
then reset both `params%flipgain` and the command-line `flipgain` to `no`, before
correction or worker dispatch. Direct shared-memory runs therefore prepare
their own gain; scheduled workers receive the prepared path and `no`, so their
helper call performs no gain I/O or flipping. Streaming prepares the gain before
constructing its worker command and explicitly sets that command's `flipgain`
to `no`. Workers must never repeat a flip already applied by their scheduler.

When a consumer instead loads a model containing a pending flip, it must pass
that saved mode into `correct_gain`. The helper loads or prepares the gain,
applies the pending flip in memory, then multiplies frames. The returned gain
must retain the same orientation for defect detection. Do not flip movie frames
to compensate for gain orientation, overwrite the input gain in this path, or
infer orientation from a filename suffix.

An ordinary image gain must match the movie dimensions. EER gain preparation
belongs to the EER decoder; `.gain` files are not ordinary non-EER gain images.
Automatic gain selection must be resolved before workers perform correction;
`auto` is not a saved binary-model flipping mode.

### Curate defects on the movie grid

Statistical outliers are detected on the gain-corrected frame sum using the
current six-standard-deviation threshold. Apply the resulting spatial mask to
every frame. Replacement estimates must exclude all masked pixels and be
computed before substituting that frame's values, so traversal order does not
feed repaired values into later estimates.

The current curation uses an 11 by 11 neighborhood: a median when sufficient
valid neighbors remain, and bounded random replacement for sparse or heavily
defective neighborhoods. Tests of stochastic replacement must use a fixed seed
and one thread.

Zero-valued pixels in the prepared EER gain are detector defects even when they
are not statistical outliers. Their orientation must match the gain actually
multiplied into the frames. Persisted outlier coordinates currently describe
the statistical detections, not the complete EER gain-defect mask; consumers
must recover the latter from the gain. The reference-free curation limitation
described below must not be mistaken for complete detector-defect coverage.

## Alignment and local model acceptance

### Global drift

Keep the discrete correlation search followed by continuous refinement. Both
must exclude the current frame's weighted contribution from its reference.
Alignment filtering and B-factor weighting belong to the search objective;
they must not become unintended extra filtering of the output micrographs.

Frame weights are initialized uniformly and may be updated from alignment
correlations. With positive correlation support, the default softmax criterion
produces nonnegative weights normalized to sum to one. Integration must use
the selected weights consistently, not introduce an implicit renormalization
in one output path. Patch alignment receives the global frame weights.

The working alignment reference is the central frame,
`nint(0.5 * nframes)`. `mcconvention=first` or `relion` re-references the applied
global shifts to frame one. A measured displacement is corrected by applying
its negative; do not reverse this sign in metadata merely because the image
operation uses a corrective shift.

### Residual local motion

Local estimation operates on globally corrected frames. It estimates residual
motion, not another copy of stage drift. Enable it only when the selected
algorithm permits patches and the effective grid has more than one patch.
`algorithm=iso` disables local correction.

Automatic patch selection uses a 200 Angstrom physical target and a minimum
patch size of 200 scaled pixels, with FFT-friendly patch dimensions. Explicit
patch counts and small or single-axis grids need geometry and fit-validity
checks; allocating patch arrays alone is not evidence of a usable model.

Fit X and Y independently by SVD using the same 18-term basis. For each spatial
monomial in the order `1, x, x^2, y, y^2, x*y`, store coefficients for
`t, t^2, t^3` consecutively. Spatial coordinates use
`pix2polycoords(p, N) = (p - 1) / (N - 1) - 0.5`, with the corresponding X or Y
dimension; time is the frame index minus the fit's reference frame. Preserve
the absence of a constant temporal term and the zero displacement at that
reference.

With `mcpatch_thres=yes`, both per-axis fit RMSDs must be strictly below the
acceptance threshold: 4 scaled pixels for `simple`, or 5 for `first`/`relion`.
On failure, retry once with each grid count rounded after halving and bounded
below by one, provided at least two patches remain. If the retry fails, generate
both products using global correction only and report the fallback.

`mcpatch_thres=no` explicitly bypasses the RMSD acceptance threshold, with a
warning for an unsatisfactory fit. Save whether local correction was actually
accepted; the presence of patch measurements or coefficients is not an acceptance
flag.

## Coordinate and reference conventions

Do not treat binary model arrays as STAR arrays. The current producer uses these
representations:

| Quantity | Binary motion model | STAR export |
|---|---|---|
| Stage drift | Scaled movie pixels; preserves the working convention supplied by global alignment | Original movie pixels, relative to frame one |
| Measured local offsets | Scaled movie pixels, relative to frame one through `get_poly4star` | The active writer exports the fitted local model rather than a measured-offset table |
| Persisted local coefficients | Scaled-pixel displacements, fitted with `t = iframe - 1` | Original-pixel displacements, with the same first-frame temporal origin |
| Patch geometry | One-based coordinates on the scaled movie grid | Not an interchangeable serialized patch-geometry record |
| Statistical outliers | One-based coordinates on the decoded, unscaled movie grid | Zero-based coordinates on that same grid |

Here `binning = smpd / smpd_movie`; converting a scaled-pixel displacement to
original movie pixels multiplies it by `binning`. EER's decoded, possibly
upsampled grid is the original movie grid for this conversion.

### Handling `motion_model` metadata

- Keep `drift_offsets_x/y` and `local_offsets_x/y` in pixels of the scaled
  working grid (`ldim`, `smpd`). `model_coeffs_x/y` must produce displacements
  in those same pixels when evaluated with the defined polynomial basis.
  Do not store STAR-scaled values in these components.
- Setters and binary I/O preserve the supplied values; they do not infer or
  convert units. Producers must supply working-grid values. Consumers applying
  motion to that grid must use them without another `binning`, image-scale or
  EER-upsampling factor. Use saved model geometry, not current project sampling.
- Convert units only at an explicit boundary. For the same reference frame,
  `offset_decoded = binning * offset_scaled`. With unchanged normalized spatial
  coordinates and time basis, every polynomial coefficient uses the same
  factor: `coeff_decoded = binning * coeff_scaled`, not powers of `binning`
  based on polynomial degree. For example, at `binning=2`, a displacement of
  0.5 working pixels is 1 decoded movie pixel.
- STAR export must convert temporary copies and must not rescale the stored
  offsets or coefficients in place. Writing STAR before the binary model must
  leave the binary values unchanged. Patch coordinates, defect coordinates,
  frame weights and dose metadata have their own conventions; do not apply
  displacement conversion indiscriminately to all model arrays.
- Keep scale and temporal reference separate. `refit_polynomial(ref_frame)`
  changes the coefficient reference, not the displacement units; its results
  remain in working-grid pixels. Evaluate them with `t = iframe - ref_frame`.
  It does not rewrite the stored measured offsets or `fixed_frame`. Do not
  serialize a consumer-refitted polynomial as a first-frame polynomial without
  restoring the required temporal convention.

`motion_model%new` accepts movie-header dimensions and physical sampling before
EER upsampling. It derives the stored decoded geometry using the same validated
`eer_scale_movie_convention(smpd, ldim, eer_upsampling, smpd_out, ldim_out)` as
`motion_correct_init`, then calculates `binning`. The subroutine supplies adjusted
sampling and dimensions through its output arguments; non-EER callers
keep their input values. Binary reading must not repeat this conversion.
The extractor retains original sampling as `smpd_physical`: saved movie sampling
for non-EER and EER mode 1, twice that sampling for EER mode 2.
This value initializes ordinary frames or the EER decoder and remains
unchanged during Fourier cropping. Do not pass decoded EER sampling to `eer%new`
or replace saved acquisition sampling with a current project's sampling. This
preserves the existing binary and STAR grid conventions.

The saved `fixed_frame` is initialized from the central working reference; it
does not redefine the first-frame temporal origin of exported coefficients.
A consumer requiring a different reference must explicitly re-reference both
axes and refit or transform the polynomial consistently. In particular, do not
evaluate first-frame coefficients at `iframe - fixed_frame` without that change.

STAR stage drift must be zero at frame one for both axes. Local polynomial
displacements must also vanish there, irrespective of the working correction
convention. STAR rows beyond `nframes` up to `total_nframes` use the current
unprocessed-frame marker `-9999` for shifts and zero frame weight; these are not
valid measured displacements.

## Dose weighting and intensity normalization

Keep two scientifically different output products:

| Project field | Integration contract |
|---|---|
| `forctf` | Motion-corrected, non-dose-weighted average with scalar weight `1 / nframes` for every retained frame |
| `intg` | Motion-corrected integration using the selected alignment frame weights, with frequency-dependent dose weights when enabled |

Generate `forctf` before applying dose weights to the reusable frame buffers.
The global-only and accepted-patch paths must use the same scalar-weight and
dose-weight conventions. Falling back to global correction must not change
the normalization or apply dose weighting twice.

`image%apply_dose_weighing` uses frequency-dependent geometric weights with
squared weights summing to one over the selected frames. At zero frequency,
each dose weight is `1 / sqrt(nframes)` for a full-movie selection. These weights
do not replace the scalar alignment weights: apply both exactly once when
forming `intg`. Do not add a frame-count multiplier to match the appearance or
mean intensity of another output.

Consequently, equal means or amplitudes between `intg` and `forctf` are not a
general correctness criterion. As an integration-stage control, for identical
constant frames with uniform alignment weights and dose weighting enabled, the
current full-movie contract
gives `forctf = constant` and `intg = constant / sqrt(nframes)`. Without dose
weighting, those products agree in this control case. Use this closed-form
control to detect omitted or duplicated scalar weighting.

## Products and persistence

For the ordinary movie workflow, publish the integrated micrograph (`intg`),
CTF micrograph (`forctf`), thumbnail (`thumb`), STAR document (`mc_starfile`), and
binary motion model (`mcmodel`, `.mmodel`). Publish patch-fit diagnostics
(`gofx`, `gofy` and the trajectory plot) when patch correction was attempted.
Update sampling and dimensions to describe the output images, while preserving
the source movie and model metadata needed to interpret them. Nano time-series
output exceptions are outside this policy's ordinary movie contract.

Write products before publishing their resolved paths in the project row.
The `bid` field is the RMS local offset over frames and the patch grid actually
fitted (the reduced grid after a retry), and zero unless patches were accepted.
Distributed assembly must preserve movie identity and the same product/model
associations as a shared-memory run.

The binary model stores the movie and gain paths, pending gain flip, original
and scaled geometry, frame counts, dose and EER metadata, global offsets, scalar
frame weights, patch geometry and measurements, fitted coefficients, acceptance
flag, reference metadata, and statistical outlier coordinates. Consumers must
honor `patch_accepted`; rejected patch arrays may still be present.

The current binary header is `FILE_VERSION=1`, `MODEL_VERSION=0`. Version 1
places a length-prefixed gain-flipping string immediately after the gain
filename. Readers reject other versions; version 0 files must be regenerated,
not interpreted by guessing a missing mode. The binary layout uses native
Fortran numeric representations and is not a promised cross-platform archival
format. Layout changes require versioning and an independent binary fixture.

Per-movie motion-correction `.star` files (`mc_starfile`) are provided for use
in other software. Within SIMPLE, the binary `motion_model` (`.mmodel`,
`mcmodel`) is the authoritative metadata source; these STAR files must not be
used as internal processing inputs or fallback metadata. STAR export remains
supported for interoperability, but does not preserve all binary model state.
Its gain path must identify a usable, correctly oriented gain because
the current STAR writer does not serialize `flipgain`. Its local-model version
indicates whether local correction was accepted. The legacy `.poly` output
is deprecated in favor of the binary `motion_model` object.

For `extractfrommov=yes`, extraction and re-extraction read only the binary
`mcmodel` metadata. A missing entry, blank path, or absent model file falls
back to integrated-micrograph extraction for that micrograph. Neither STAR
nor legacy `.poly` is a fallback metadata source. Invalid binary models and
missing movies referenced by an existing model remain errors.

## Current implementation limitations

These limitations qualify the requirements above; they are not alternative
scientific conventions to propagate into new code.

- In `cure_outliers` in `simple_motion_correct`, adding zero-valued EER gain pixels
  to the repair mask is inside the branch for at least one statistical outlier.
  A gain-only defect can therefore remain uncurated, and the early return for
  a nearly constant frame sum also bypasses this work. Complete EER defect
  handling must be tested independently of statistical detection.
- For EER, STAR dose rate currently comes from `target_dose_per_frame`, whereas
  the binary also stores the effective `dose_per_frame`. Grouping and truncation
  can make these quantities differ. Do not assume that STAR and binary dose
  metadata are numerically interchangeable.
- `corrs2weights` returns all-zero weights when every correlation is nonpositive.
  Global alignment warns on negative mean correlation but does not reject the
  movie at that point. Successful execution alone therefore does not establish
  a valid, nonzero integrated micrograph.

## Validation and change control

Follow the [test environment policy](test_environment_policy.md). Numerical
tests need independent expected values, fixed seeds where applicable, explicit
tolerances, and assertions that affect exit status. A successful program exit
or visually plausible micrograph alone is insufficient.

The `motion gain`, `motion model` and `particle extractor` sub-suites of
`unit_motion` provide gain-helper checks, binary fixtures, polynomial-refit
checks and movie-extraction initialization/fallback checks. The
image sub-suite of `unit_image` covers image flips. These are starting points,
not evidence of complete end-to-end coverage. Changes affecting this workflow
must cover the relevant contracts:

- Known global translations, corrective-shift sign, and leave-one-out reference
  construction, including the zero-motion control.
- Known local coefficients on a non-square grid and a non-first working
  reference; first-frame STAR export and scaled/original pixel conversion.
- Non-unit `binning` for both offset axes and polynomial coefficients;
  unchanged binary values before/after STAR export and no consumer-side
  double scaling.
- Accepted patches, rejected patches, the reduced-grid retry, forced threshold
  bypass, and global-only fallback.
- All gain modes, no-gain processing, prepared-gain reuse without double
  flipping, and direct versus worker execution.
- Statistical-only, gain-only, overlapping, and boundary defects; no valid
  neighborhood and nearly constant frames.
- EER grouping with leftover raw frames, upsampling, dose truncation endpoints,
  and separation of target, nominal, and accumulated dose.
- Closed-form `intg` and `forctf` intensity controls, with and without dose
  weighting, through both global and local integration paths.
- Binary round trips with and without patch/outlier arrays, rejected patches,
  differing runtime parameters, and explicit version rejection.
- Short or invalid movies, cleanup between movies, and project/product identity
  after distributed assembly.

Changes to signs, units, temporal origins, gain semantics, normalization,
acceptance rules, or binary layout require coordinated producer/consumer changes
and tests. Update this policy and the algorithm note when their respective
contracts change. New optional scientific behavior requires explicit opt-in;
do not silently change existing defaults or weaken acceptance tests to make a
new implementation pass.

## Implementation references

- [Movie iterator](../../src/main/motion/simple_motion_correct_iter.f90)
- [Frame preparation and integration](../../src/main/motion/simple_motion_correct.f90)
- [Hybrid alignment](../../src/main/motion/simple_motion_align_hybrid.f90)
- [Patched alignment and deformation](../../src/main/motion/simple_motion_patched.f90)
- [Shared motion utilities](../../src/main/motion/simple_motion_correct_utils.f90)
- [Motion model persistence](../../src/main/motion/simple_motion_model.f90)
- [Motion correction execution](../../src/main/strategies/parallelization/simple_motion_correct_strategy.f90)
- [Preprocessing execution](../../src/main/strategies/parallelization/simple_preprocess_strategy.f90)
- [Image dose weighting](../../src/main/image/simple_image_ops.f90)
- [Correlation to frame weights](../../src/utils/math/simple_stat.f90)
