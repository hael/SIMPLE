# Motion Model Fraction Re-extraction and Reconstruction Policy

## Purpose and scope

This policy governs the current `refine_motion_model` application: re-extract
particle fractions using saved movie motion, assemble a project for each
exposure interval, then reconstruct those intervals from saved particle poses.
It defines the execution, scientific, project, and output contracts that changes
must preserve. Current limitations and the extent of validation are stated
separately below.

The application does not perform reference-based particle-frame alignment,
estimate new trajectories, or publish a newly refined motion model. Its local
polynomial refit uses saved patch measurements. The reconstructed interval maps
do not feed back into motion estimation or particle alignment in this workflow.

The [algorithm description](../algorithms/refine_motion_model.md) explains the
coordinate operations and summation. The
[motion-correction policy](motion_correction_policy.md) owns the producer's
binary-model, EER, gain, and coordinate conventions. The
[layering and project-ownership policy](layering_and_project_ownership_policy.md),
[stack-indexing policy](project_stack_indexing_policy.md), and
[particle-read policy](particle_read_path_policy.md) apply to the data handoffs.

## Execution and ownership

The commander normalizes the command line inline — the missing-key defaults
and the fixed reconstruction settings below — then retains one strategy
lifecycle:

```text
initialize -> execute -> finalize_run -> cleanup
```

For the public master run, `execute` is `distr_execute`. It performs these tasks
in order, with a barrier between extraction and assembly:

1. Create interval directories, dispatch extraction workers, and wait for them.
2. Call private `assemble_motion_fractions` to merge part projects, write every
   complete interval project, and remove consumed per-part interval stack tables.
3. Call `reconstruct_motion_fractions` on the returned project paths and ranges.
   Complete each interval reconstruction before beginning the next.

Assembly and reconstruction are execution work. Keep them out of
`finalize_run`, which retains completion bookkeeping, and out of `cleanup`,
which releases resources and temporary queue files.

| Owner | Responsibility |
|---|---|
| UI and parameters | Declare inputs, own baseline defaults, derive and validate parsed settings |
| Motion commander | Normalize the command line, select the strategy and drive its existing lifecycle |
| Distributed master | Establish shared half sets, schedule work, own assembled projects, sequence interval reconstructions and register delivered products |
| Extraction worker | Read its assigned movies/particles, extract stacks, and write part projects |
| `assemble_motion_fractions` | Own assembly-local workspaces and return completed interval paths/ranges |
| `ptcl_extractor` and motion-model modules | Decode/correct movies, interpret saved motion, extract frames, and normalize sums |
| `reconstruct3D` program (strategy, service and volume assembly) | Insert saved poses, restore maps/half maps and calculate FSCs in each nested interval run |

Workers read the execution project and write their part files; they do not
assemble global projects or start interval reconstructions. Worker completion
is declared in the existing worker finalization path. A public `nparts=1` run
still uses the master strategy and performs the complete workflow. Direct
partition invocations are worker executions, not equivalent public runs.

The master is the sole writer of each assembled project. Each reconstruction
owns a separate project copy and runs as a nested, complete `reconstruct3D`
program on it, through its strategy routine (`run_rec3D` in
`simple_rec3D_strategy`) with its own freshly parsed parameters and builder. The master's own `parameters` object is never edited
or copied for a reconstruction. Heavy stored components and imports follow the
[compile-time policy](compile_time_policy.md).

## Inputs and interval selection

Input is a SIMPLE project with movie/stack associations, binary `mcmodel`
records, particle box coordinates, and usable 2D/3D particle records. Processed
movies, saved models, and required gain files must be accessible. Binary model
reading validates its supported format; STAR exports and legacy polynomial
sidecars are not interchangeable extraction inputs.

The public UI exposes `pgrp`, `mskdiam`, `nparts`, and `nthr` for the combined
workflow, along with extraction controls. Optional `box` changes the extraction
box; it must satisfy normal parameter checks and the extractor's even-box
requirement. This application always writes particle stacks as MRC mode 12
(float16), and does not expose `wfloat16` in its UI. Its commander forces
`wfloat16=yes`, so a supplied shared CLI key cannot override this storage rule.
Other applications retain the optional `wfloat16=yes|no` setting, defaulting
to `no`, including `extract` and `reextract`.
Storage precision does not change the integration rule. The default
recentering segment is `ptcl3D`; the existing
`oritype=ptcl2D` extraction path can instead supply the integer recentering
translation. Reconstruction always reads saved 3D poses.

`fromf`, `tof`, and `stepf` index the decoded fractions saved in the model. For
EER, one such fraction contains the saved group of raw EER frames. Do not
reinterpret these indices as individual raw EER frames or regroup the movie
using unrelated runtime dose settings.

- `fromf` defaults to 1 and must be positive.
- `stepf` defaults to 3 and must be positive. Explicit widths remain supported.
- `tof=0` requests automatic resolution of the final fraction. The current
  master resolves it once from the first eligible movie model and uses that
  common endpoint for the run; it is not resolved separately for each movie.
- After endpoint resolution, `tof` must not precede `fromf`. Extraction checks
  the selected range against each model it processes.

Keep the parameter baseline authoritative. The UI publishes generated defaults
under the [UI policy](ui_layer_policy.md); it must not own a separate interval
width. The reconstruction mode assignments below are fixed workflow settings,
not alternative UI defaults.

Output blocks are inclusive, normally full-width, and not cumulative. Anchor a
short final remainder at `tof`, extending it backward to the requested width
when the selected interval is long enough. Thus `fromf=1,tof=8,stepf=3` gives
`1-3`, `4-6`, `6-8`. A selected interval shorter than `stepf` is used once.
Overlapping outputs share observations and must not be treated as independent
exposure groups or summed as disjoint data.

## Extraction and intensity contract

Keep only particles active in both the 2D and 3D records. Recenter the extraction
corner by the rounded translation from the selected segment and apply the same
integer correction to both residual-shift records. Preserve the physical center
when changing the box. Reconstruction must use the resulting saved residual
translations, not reapply the integer shift or fit new particle poses.

The binary model owns movie geometry, acquisition sampling, EER decoding,
gain orientation, saved detector defects, and motion units. Respect the
producer/consumer conventions in the motion-correction policy:

- Apply gain correction and repair the recorded detector defects, including
  the supported EER gain-defect handling, on the appropriate decoded grid.
- Fourier-scale frames to the model's saved working geometry when needed.
- Re-reference global drift to the model's saved fixed frame. Refit accepted
  local patch measurements consistently at that reference; rejected or absent
  local correction uses global drift alone.
- Apply the displacement at the particle center as an extraction translation
  and a residual Fourier phase shift. Do not add another binning or EER scale
  factor to saved working-grid displacements.

Cache the corrected particle frames for the selected range. Zero each output
particle before summing only its requested block, with equal weights. This path
explicitly disables dose weighting. Saved frame-quality weights are loaded as
model metadata but are not applied to these block sums.

After each sum, apply the existing contrast convention, background-ramp
subtraction, and background-noise normalization. These are processed particle
images, not preserved detector-count sums. No dose-dependent or frame-quality
weight may be introduced silently under the same output contract.

## Assembly and project handoff

Assemble every interval from the completed worker parts using the existing
project adders and whole-row transfers. Preserve particle-to-stack indexing,
recentered coordinates, residual shifts, CTF metadata, state labels, and the
other retained particle fields. Movies with no retained particles are excluded.
Keep count consistency checks and consumed-part cleanup with assembly.

During assembly, every interval project must point to its own physical fraction stacks and name
its actual project file. Return the exact completed paths and bounds to the
reconstruction coordinator; do not rediscover intervals from a directory scan.
No reconstruction starts while an interval project is still being assembled.

The master calls the project-owned half-set helper before worker dispatch.
Preserve the shared partition through extraction and assembly. Do not generate
independent partitions for individual intervals or relabel particles from map
agreement. Reconstruction checks that selected particle rows and their half
labels agree across interval projects.

The execution project is the master's mutable state. With the default
`mkdir=yes`, normal execution-context setup copies the input project into the
run directory. With `mkdir=no`, shared half-set bookkeeping can be written into
the supplied execution project; do not promise that mode leaves it untouched.

## Reconstruction contract

Reconstruction is automatic for a public run with assembled particles. Use the
existing `reconstruct3D` gridding route with hard state assignments and these
fixed settings:

```text
objfun=cc ml_reg=no cc_emit_sigma=no rec_backend=gridding
filt_mode=none automsk=no postprocess=no unfil_unmasked=yes
m_estimator=no update_frac=1 trail_rec=no
```

There is no noise-model bootstrap, sigma2 loading/whitening, residual sigma2
emission, ML prior, or trailing update in this step. CC alone is not a universal
promise of no sigma2 in other applications; this workflow also fixes
`cc_emit_sigma=no` and does not call their sigma-estimation stages. An inherited
sigma2 registration may remain in the copied project's metadata, but it must
not be consumed or regenerated by this reconstruction.

Use existing orientations, residual shifts, state assignments, and CTF metadata.
The positive internal translation setting enables saved shifts; it introduces
no translation search. Select rows through `sample4rec`: active rows with
`updatecnt>0` when any active row has been updated, otherwise all active rows.
An extracted row therefore need not be a reconstruction contributor. Verify the
same selected rows and half labels for every interval.

Native interval-stack box/sampling are authoritative and must agree with the
prepared parameters. The supported application inputs reconstruct that native
interval geometry; changing the extraction `box` is a separate operation.
Additional hidden reconstruction-crop arguments are not a supported
cross-dispatch contract. Do not substitute current project sampling for
acquisition/model geometry in extraction.

Normal particle preprocessing, gridding density/interpolation correction and
FSC calculation remain active. Write the `_unfil` half-map exports before
spherical support; standard half maps and the merged reconstruction retain
that support and the supplied symmetry/mask size.

Run each interval as a nested `reconstruct3D` program on its own project copy:
a sparse command line carrying the fixed settings, symmetry, mask diameter and
thread count, without movie partition indices or unrelated refinement
controls. The nested run parses its own parameters and builds its own
toolboxes. One partition reconstructs in the master process (`nparts` is not
forwarded); multiple partitions run the program's own distributed gridding
master. Restore the master's working directory and thread context after each
run; retain its logging destination.

## Products and publication

For each interval `a-b`, write:

| Directory | Products |
|---|---|
| `frames_a-b/` | Assembled input-named `.simple` project and `ptcls_from_<movie>_<a>_<b>.mrcs` stacks; the directory and its contents are retained after reconstruction |
| `rec_frames_aa-bb/` | A separate input-named project copy, populated-state reconstruction maps, half maps, FSCs, and normal service artifacts |

Reconstruction endpoints are padded to at least two digits, for example
`rec_frames_01-03`. The current implementation rejects an existing
reconstruction directory instead of reusing its files.

Within each reconstruction iteration, clean up only after all populated states
have delivered valid maps/FSCs, their FSCs have been copied into the in-memory
matrix, and output registration has been written. Failed, incomplete or skipped
reconstructions retain their fraction folders and partial volumes.

Retain every `frames_a-b/` directory, its assembled project and its particle
stacks after both successful and unsuccessful reconstruction. These remain
available for subsequent analysis without re-extraction.

In the completed reconstruction directory, remove only the even/odd gridding
partials named by `refine3D_partial_rec_fname` and `refine3D_partial_rho_fname`
for the run's exact states and parts. Validate the canonical run root and exact
interval reconstruction directory before deletion. Require existing targets to be non-directory
file entries without redirection, and preflight them before deletion. Do not use
broad `recvol*` globs: final maps, half maps, FSCs and the copied project remain.
Verify cleanup and log the removed partial-file count.

The copied reconstruction project's stack/particle rows continue to reference
the retained fraction stacks. Raw FSC export and plotting use the collected
in-memory matrix.

Before the reconstruction loop, allocate an FSC array with dimensions
`(native Fourier shells, intervals, states)` using the assembled project's
state count and the native box/sampling. Read each delivered binary FSC directly
into its interval/state column, checking the shell count and finite values.
Retain signed correlations and endpoints unchanged. Missing curves, including
unpopulated states, remain `NaN`; they are not zero-correlation observations.

After all intervals finish, write `fsc_fractions_stateNN.tsv` in the master
execution directory for each state. Columns are `shell`, `resolution_A`,
`reciprocal_frequency_A^-1`, followed by `FSC_fractions_<first>-<last>` for
each actual inclusive interval. Shells run from 1 to `box/2`, without DC;
resolution is `box*smpd/shell` and reciprocal frequency is its inverse.
The TSV retains enough digits to round-trip the binary single-precision FSCs.
The raw diagnostics below use this same matrix.

Use the reconstruction service's canonical state-map, half-map and FSC names.
Update copied-project `projinfo` to its own project and directory. Remove
inherited state-associated map/FSC entries through the existing project helper,
then explicitly register delivered interval maps and FSCs in `out` through its
adders. Do not leave original full-exposure products advertised as delivered
interval maps, including for states without reconstructed particles. Unrelated
retained output metadata is not a declaration of a newly generated product.

Only publish reconstruction registrations for products actually delivered.
Write modified segments through their normal ownership paths. A failed
reconstruction must not be reported as a successful complete application run;
partial directories are not resumable checkpoints under the current contract.

## Dose-series half maps and manifest

Every populated state's `recvol_stateNN_even_unfil.mrc` and
`recvol_stateNN_odd_unfil.mrc` contains the density-corrected, deapodized,
unregularized half volume before spherical masking, written as real-space
float32 MRC (mode 2). The internal `unfil_unmasked=yes` setting is fixed on the
fraction reconstruction's child command; other reconstruction callers retain
the default masked `_unfil` behavior. Support provenance identifies the new
exports as `none`. Standard half maps, the merged map and their FSC retain
spherical support and the supplied `mskdiam`.

After successful reconstruction and raw FSC export, the master writes
`dose_series.txt` for state 1 when its physical dose metadata are usable.
The shared `dose_series_manifest` class owns its versioned reader/writer and
validates actual MRC headers before the consumer reads voxel arrays. The
manifest contains the box, sampling, saved voltage, dose per fraction,
selected-particle and half counts, a partition fingerprint, and one row per
actual interval with fraction bounds, accumulated physical dose bounds and
relative even/odd paths. Paths are relative to the manifest directory so the
complete run directory can be moved. Version 1 separates columns with spaces
or tabs; path tokens cannot contain whitespace.

Use only the actual `sample4rec` state-1 contributors to identify motion-model
metadata. Dose per fraction and voltage come from their saved binary motion
models; do not use a current CLI default or an EER target dose. The stored dose
already includes EER grouping. For an original inclusive interval `[a,b]` and
saved dose per decoded fraction `d`, write `dose_lo=(a-1)*d` and `dose_hi=b*d`
in double precision, in physical e/Angstrom^2. Retain the dose before `fromf`
when the selected interval starts after fraction 1. Voltage scaling belongs
to the consumer.

All contributing micrographs must have compatible positive finite doses,
supported voltage and sufficient saved fractions. Check the state-1 selected
rows and half labels across intervals and fingerprint their ordered identities.
A volume-only consumer can record that invariant but cannot independently
verify the underlying particle partition. Unused micrographs do not determine
the state-1 dose metadata.

Keep the final overlapping interval and its predecessor in the manifest. The
consumer validates the supported overlap pattern, drops the predecessor and
requires at least two disjoint surviving windows. Publish the complete
manifest by temporary-file write followed by rename only after validating its
metadata and all referenced unmasked halves. A failed interval must not leave
a partial manifest presented as complete.

Missing or incompatible physical metadata, unavailable state-1 halves, a
changed state-1 sample or insufficient disjoint intervals suppresses this new
manifest with an explicit reason. It must not discard valid maps, stacks or
raw FSC outputs. Existing fatal reconstruction and global partition checks
remain fatal. Report file-publication errors explicitly. No dose fitting is
performed by `refine_motion_model`; `dose_filter3D` consumes the manifest in a
separate application.

## Raw FSC diagnostics

After exporting the matrix, generate one two-panel figure for each state in
`fsc_fractions_stateNN_diagnostics.png`. Do not generate PDF figures.
The existing reconstruction strategy
owns this reporting step and consumes its in-memory FSC matrix directly.
Keep the command file `fsc_fractions_stateNN.gp` and numerical surface table
`fsc_fractions_stateNN_surface.tsv` beside the raw matrix TSV so the figures
can be reproduced. The [implementation record](../implementation_notes/completed/fraction_fsc_diagnostics.md)
records the surviving export and plotting work.

| Panel | Content |
|---|---|
| Raw FSC curves | One curve per actual inclusive fraction interval; reciprocal frequency in Å⁻¹ on x and FSC on y, with exact interval labels |
| Raw FSC surface | The complete raw matrix as a 3D surface; reciprocal frequency in Å⁻¹ on x, fraction-interval midpoint index on y, and FSC on z |

For an interval `[a,b]`, its midpoint index is `(a+b)/2`; preserve the original
fraction numbers. Label or explicitly map these coordinates to the exact
inclusive intervals, including overlapping final intervals. These are saved
fraction indices, not physical dose coordinates. The plotting routine uses
these indices and does not require a physical dose axis.

Preserve all raw finite correlations, including negative values, zero and
endpoints, without clipping or display caps. Missing curves remain `NaN` in
the tables and missing in both panels. Do not interpolate across missing
values or bridge missing intervals with an invented surface. A single interval,
single shell, or empty state must remain an honest sparse/empty display rather
than a fabricated grid. Keep the native frequency grid and actual interval
order; do not resample the observations for display.

Generate the raw FSC figure through Fortran-controlled
gnuplot, using the headless `pngcairo` terminal and the repository's existing native plotting
patterns. There is no Python or Matplotlib runtime dependency and no
interactive display or pause. Report unavailable gnuplot, renderer failures,
or missing figure outputs as warnings. Retain the raw matrix TSV and any
available command/data artifacts, and report which plots were not produced.
A plotting failure does not discard successful reconstructions or the raw
FSC export.

The application stops after unfiltered reconstruction and raw FSC export. It
does not preprocess SNR, select fitting points, fit decay or critical-exposure
models, or generate fitted parameters, predictions, residuals or new dose
weights. No fitting-coverage or dose-compatibility gate may reject the raw
output. Existing standalone Python analysis scripts remain independent.

## Current limitations

- Automatic `tof` assumes a common usable interval across the selected movies.
  A longer movie can be truncated to the first model's endpoint; a shorter one
  can fail during extraction. There is no up-front all-model range audit.
- The worker associates active movies and stacks in their existing row order.
  That association must already be consistent; arbitrary reordered or regrouped
  particle stacks are not mapped back to movies by an explicit identity join.
- Interval provenance is carried by directory/stack names and the execution's
  path/range list. The application does not yet add a dedicated integration
  range/rule record to every stack for downstream dose-history reconstruction.
- Detector-defect repair can draw replacement samples. Independent unseeded
  runs, or different process partitions, need not produce byte-identical images.
  This does not justify changing weighting or normalization conventions.
- Workers hold the decoded movie and cache all requested particle frames before
  forming the output blocks. The application does not stream one block at a
  time through a bounded particle-frame cache.
- Generic CLI keys can still reach the parameter layer even when absent from
  this program's UI. In particular, a hidden `box_crop` affects the master's
  prepared parameters but is not forwarded to the nested interval
  reconstructions. Such overrides are not validated as a supported mode. The
  fixed reconstruction settings above govern ordinary application use.

## Validation and change control

Follow the [test-environment policy](test_environment_policy.md). Expected
values must come from independent calculations or synthetic truth, with fixed
seeds and declared tolerances where needed. A normal stop or plausible image
alone is not evidence of numerical correctness.

Changes must cover the affected contracts:

- Baseline/default publication, explicit interval-width overrides, inclusive
  endpoints, non-first starts, short selections and the final overlap.
- Master/worker ordering, all-part assembly, stack/particle indexing, completed
  project handoffs and consistent selections/half labels across intervals.
- Preserved poses, residual shifts, CTF/native geometry, and extraction-box
  changes when that path is affected.
- No sigma2 dependency or emission, including a project with an unavailable
  inherited sigma2 registration; no allocation of a whitening spectrum.
- Copied-project ownership, removal of stale state products, correct map/FSC
  registration, missing/empty inputs, and caller-context restoration.
- Numerical map recovery against synthetic truth through the nested CC
  gridding route.
- Raw matrix export: binary-FSC round trips, native shell/frequency coordinates,
  exact interval labels, signed values/endpoints, and missing-state NaNs.
- Plot verification: both raw panels, uncapped FSC values, midpoint-index and
  interval mapping, missing-data gaps, sparse/empty states, and nonfatal
  renderer failures that retain the numerical export.
- Complete one- and multiple-partition runs. For structural changes, compare
  before/after binary products and normalized project contents with the same
  seed, partitioning and runtime settings.

Remaining automated coverage includes the motion-model/extractor tests in
`unit_motion`, the UI checks in `unit_ui`, and `fractional_reconstruction` in
`lib_reconstruction`. The dedicated interval and motion-fraction reconstruction
testers have been removed at the maintainer's request. These remaining tests
do not provide direct coverage of the complete interval reconstruction
workflow or all real detector, EER, movie-length, crop, or queue-system
combinations.

Changing the integration rule, temporal grouping, motion conventions, pose
handling, reconstruction objective, sigma2 behavior or artifact ownership
requires coordinated source/UI/documentation updates and appropriate numerical
and workflow checks. Keep the algorithm explanation aligned with this policy;
record implementation history and dated results in implementation notes.

## Implementation references

- [Application strategy](../../src/main/strategies/parallelization/simple_refine_motion_model_strategy.f90): range helpers, worker extraction, master dispatch, `assemble_motion_fractions`, and `reconstruct_motion_fractions`.
- [Motion commander](../../src/main/commanders/simple/simple_commanders_motion.f90) and [UI](../../src/main/ui/simple/simple_ui_preproc.f90): command-line defaults, lifecycle and public input contract.
- [Particle extractor](../../src/main/preprocess/simple_particle_extractor.f90) and [binary model](../../src/main/motion/simple_motion_model.f90): saved-motion extraction and equal-weight normalized sums.
- [Reconstruct3D strategy](../../src/main/strategies/parallelization/simple_rec3D_strategy.f90) (`run_rec3D`), [reconstruction service](../../src/main/strategies/parallelization/simple_rec3D_service.f90) and [volume assembly](../../src/main/strategies/parallelization/simple_volassemble_strategy.f90): the nested interval reconstruction runs.
- [Project half-set ownership](../../src/main/project/simple_eo_partition.f90) and [reconstruction row selection](../../src/main/ori/simple_oris_sampling.f90): shared labels and `sample4rec`.
