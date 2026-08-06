# Automasking Policy

## Scope

This document describes 3D automasking policy for `refine3D`, `abinitio3D`, `reconstruct3D`, and the shared `volassemble` step.

The current architecture treats automasking as:

- a workflow policy selected by `automsk`
- a per-state artifact named by convention
- a post-reconstruction operation owned by `volassemble`
- an opt-in input to phase-randomized FSC solvent correction
- an opt-in input to matching-reference solvent flattening

## Public policy

The user-facing control is:

- `automsk=no`: no automask generation; FSC uses the spherical mask
- `automsk=yes`: generate a per-state envelope mask from the even/odd volumes
- `automsk=tight`: same as `yes`, but request tighter Otsu-style thresholding

`envfsc=no` is the default and keeps the automask out of FSC estimation.
`envfsc=yes` requests phase-randomized solvent correction when a compatible
automask is available in a nonuniform refinement.

`envref=no` is the default and retains the spherical matching-reference mask.
`envref=yes` requests RELION-style solvent flattening of each matching
reference with its compatible state automask. It changes only the references
used to generate projections; particle images, FSC estimation, NU filtering,
and matching-bandwidth selection are unchanged.

NU filtering always uses the spherical support derived from `mskdiam`.
Automasks do not define or restrict the NU objective domain.

A future NU-evidence envelope is a separate correlation-derived artifact. It
must not overwrite `automask3D_stateNN.mrc`, feed FSC correction, or replace
spherical NU support. Its intended workflow consumer is matching-reference
masking after independent validation and temporal recovery safeguards.

`mskfile` is no longer part of the CLI policy. Passing `mskfile` is a hard error.

## Ownership

### `volassemble`

`volassemble` is the sole producer of workflow automasks.

Responsibilities:

- decide whether a state mask should exist
- generate the mask from the current even/odd pair
- persist the result as `automask3D_stateNN.mrc`

### `reconstructor_eo`

The reconstructor is a consumer, not a producer.

For a nonuniform refinement, the reconstructor first calculates the provisional
gold-standard FSC with the broad spherical mask. After `volassemble` has
generated or loaded the current automask, it optionally replaces that curve
with the phase-randomized solvent-corrected FSC when `envfsc=yes`.

Legacy non-NU reconstruction paths may still consume a compatible state mask
directly when `envfsc=yes`; they fall back to the spherical mask when it is
missing.

The reconstructor no longer regenerates missing masks on demand.

### Matcher reference preparation

With `envref=yes`, reference preparation loads the compatible
`automask3D_stateNN.mrc`, subtracts the median density in its soft transition,
and multiplies the selected reference volume by the mask before reprojection.
This is applied after choosing regular versus NU-derived and independent
versus merged references. A missing or incompatible mask falls back to the
spherical reference mask for that iteration.

## State-specific artifacts

Mask files are named:

- `automask3D_state01.mrc`
- `automask3D_state02.mrc`
- ...

This avoids the old single-file collision problem in multi-state workflows and keeps mask ownership aligned with the state model.

## Current workflow

1. `refine3D` or `abinitio3D` produces partial reconstructions.
2. `volassemble` builds even, odd, and merged state volumes.
3. If `automsk != 'no'`, `volassemble` may generate `automask3D_stateNN.mrc`.
4. If a NU `filt_mode` is active, the NU filter constructs spherical support from `mskdiam`.
5. With `envfsc=no`, the FSC remains the broad-sphere curve and the automask has
   no FSC role.
6. With `envfsc=yes`, `volassemble` calculates unmasked/broad-sphere, masked,
   and randomized-phase masked curves, then writes the corrected curve as the
   state FSC used for resolution and bandwidth decisions.
7. With `envref=yes`, the next iteration solvent-flattens its matching
   references with the compatible state mask before projection.

## Current implementation notes

- Multi-state automasking is supported.
- Masks are internal workflow artifacts; they are not recorded as explicit CLI outputs.
- `tight` should be preserved end-to-end as a policy value, not collapsed to a boolean.
- The implementation regenerates masks when they are missing or incompatible, at `startit`, and every `AMSK_FREQ` iterations.
- NU filtering does not consume these masks.

## Compatibility rules

A state mask is compatible only if:

- its box matches the current cropped box
- its sampling matches the current cropped sampling

Dimension-only checks are not sufficient once stage-dependent rescaling or recropping is in play.

## Recommended direction

Near-term improvements:

- Preserve `automsk=tight` all the way through staged `abinitio3D` control logic.
- Use the same compatibility rule in both `volassemble` and `reconstructor_eo`.
- Move mask-production decisions into a small helper or policy object so `volassemble` owns execution, but not every branch itself.

Longer-term architectural target:

- `volassemble` should remain the execution site for expensive volume-domain work.
- Mask policy should be represented explicitly and threaded into `volassemble`, rather than being reconstructed ad hoc from scattered local conditions.
