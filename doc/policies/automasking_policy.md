# Automasking Policy

## Scope

This document describes 3D automasking policy for `refine3D`, `abinitio3D`, `reconstruct3D`, and the shared `volassemble` step.

The current architecture treats automasking as:

- a workflow policy selected by `automsk`
- a per-state artifact named by convention
- a post-reconstruction operation owned by `volassemble`
- an optional FSC mask input consumed by the reconstructor

## Public policy

The user-facing control is:

- `automsk=no`: no automask generation; FSC falls back to the spherical mask
- `automsk=yes`: generate a per-state envelope mask from the even/odd volumes
- `automsk=tight`: same as `yes`, but request tighter Otsu-style thresholding

`mskfile` is no longer part of the CLI policy. Passing `mskfile` is a hard error.

## Ownership

### `volassemble`

`volassemble` is the sole producer of workflow automasks.

Responsibilities:

- decide whether a state mask should exist
- generate the mask from the current even/odd pair
- persist the result as `automask3D_stateNN.mrc`
- reuse the same state mask as input to nonuniform filtering when available

### `reconstructor_eo`

The reconstructor is a consumer, not a producer.

Responsibilities:

- try to load `automask3D_stateNN.mrc`
- verify that the mask matches the current box and sampling
- apply it to FSC volumes only when `envfsc=yes`
- fall back to the spherical mask when no compatible state mask exists

The reconstructor no longer regenerates missing masks on demand.

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
4. If `filt_mode=nonuniform`, `volassemble` prefers that state mask as the nonuniform-filter support mask.
5. During FSC calculation, the reconstructor attempts to load the same state-specific mask.
6. If no compatible mask exists, FSC uses the spherical fallback mask.

## Current implementation notes

- Multi-state automasking is supported.
- Masks are internal workflow artifacts; they are not recorded as explicit CLI outputs.
- `tight` should be preserved end-to-end as a policy value, not collapsed to a boolean.
- The implementation currently behaves as "generate if missing or incompatible". If periodic refresh every `AMSK_FREQ` iterations is desired, that must be enforced explicitly in `volassemble`.

## Compatibility rules

A state mask is compatible only if:

- its box matches the current cropped box
- its sampling matches the current cropped sampling

Dimension-only checks are not sufficient once stage-dependent rescaling or recropping is in play.

## Recommended direction

Near-term improvements:

- Preserve `automsk=tight` all the way through staged `abinitio3D` control logic.
- Regenerate masks on the intended cadence instead of only when files are missing.
- Use the same compatibility rule in both `volassemble` and `reconstructor_eo`.
- Move mask-production decisions into a small helper or policy object so `volassemble` owns execution, but not every branch itself.

Longer-term architectural target:

- `volassemble` should remain the execution site for expensive volume-domain work.
- Mask policy should be represented explicitly and threaded into `volassemble`, rather than being reconstructed ad hoc from scattered local conditions.
