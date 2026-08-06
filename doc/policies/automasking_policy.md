# Automasking Policy

## Scope

This document describes 3D envelope-mask policy for `refine3D`, staged
`abinitio3D`, and the shared `volassemble` step. Standalone density masking
commands remain available, but density/Otsu masks are no longer generated in
the refinement loop.

The current architecture treats refinement automasking as:

- a workflow policy selected by `automsk`
- a NU-evidence-derived per-state artifact named by convention
- a post-reconstruction operation owned by `volassemble`
- an input to matching-reference solvent flattening whenever `automsk != no`

## Public policy

The user-facing control is:

- `automsk=no`: no refinement-loop envelope generation
- `automsk=yes`: generate a per-state envelope from NU cross-half evidence
- `automsk=tight`: retained as a compatibility spelling of `yes`; the NU
  envelope has no density/Otsu tightness mode

In 3D refinement, `automsk=yes|tight` requires
`filt_mode=nonuniform|nonuniform_lpset`. `filt_mode=none|uniform|fsc` requires
`automsk=no`. Standalone density-mask utilities are outside this refinement
invariant.

`envfsc=no` is the default. **`envfsc=yes` is an unfinished branch.** The density
automask that used to feed it is no longer generated in the loop, and the
NU-evidence envelope must never reach the FSC. It currently exits with a hard
error. Re-enabling it requires a simple, fast density-based mask generated on
the fly from the current map.

Matching references are solvent-flattened with the current, lagged NU-evidence
envelope whenever `automsk != no`; there is no separate `envref` control (it has
been removed from the CLI and the library). A compatible legacy density mask is a
transition fallback; otherwise the first iteration uses the sphere. This changes
only the references used to generate projections. Particle images, FSC
estimation, NU filtering, and matching-bandwidth selection are unchanged.

NU filtering always uses the spherical support derived from `mskdiam`.
Envelope masks do not define or restrict the NU objective domain. In
particular, the correlation-derived NU envelope must never feed FSC correction
or replace spherical NU support.

`mskfile` is no longer part of the CLI policy. Passing `mskfile` is a hard error.

## Ownership

### `volassemble`

`volassemble` is the sole producer of refinement-loop envelope masks. It:

- decides whether a state NU envelope should exist
- derives it after static NU optimization and accepted `nu_refine` extensions
- regenerates it freely each cycle, allowing the envelope to shrink or grow with
  resolution; there is no monotonic recovery guard and no `pct_signal` failsafe
- persists the result as `nu_envmask3D_stateNN.mrc`

The envelope is generated before `nu_filter_vols` releases the mask-packed NU
unary storage. Iteration N therefore writes the mask consumed while iteration
N+1 materializes its reprojection model.

### FSC consumers

The NU-evidence envelope is never used for FSC correction because it is selected
from cross-half agreement. A nonuniform refinement first calculates the
provisional gold-standard FSC with the broad spherical mask. **`envfsc=yes` is an
unfinished branch:** parameter validation exits with a hard error before volume
processing begins. Supporting it requires generating a simple, fast density
mask on the fly and passing that transient mask directly to phase-randomized FSC
correction.

### Matcher reference preparation

Whenever `automsk != no`, reference preparation first loads the compatible
`nu_envmask3D_stateNN.mrc`, subtracts the median density in its soft transition,
and multiplies the selected reference volume by the mask before reprojection.
This is applied after choosing regular versus NU-derived and independent versus
merged references. The fallback chain is:

1. NU-evidence envelope
2. compatible legacy density envelope
3. spherical reference mask

## State-specific artifacts

Current refinement envelopes are named:

- `nu_envmask3D_state01.mrc`
- `nu_envmask3D_state02.mrc`
- ...

`automask3D_stateNN.mrc` is a legacy density-envelope name. It is no longer
produced in the refinement loop, but may be consumed as the second
reference-mask fallback. The FSC path does not consume it.

## Current workflow

1. `refine3D` or staged `abinitio3D` produces partial reconstructions.
2. `volassemble` builds even, odd, and merged state volumes.
3. The NU filter constructs spherical support from `mskdiam`, optimizes the
   static bank, and accepts any supported `nu_refine` extensions.
4. If `automsk != 'no'`, `volassemble` derives
   `nu_envmask3D_stateNN.mrc` before NU unary storage is released.
5. The NU envelope is used only by the next matching-reference preprocessing
   pass and never by FSC correction or NU support.
6. `envfsc=yes` exits with a hard error until a fast density mask can be
   generated on the fly for phase-randomized correction.

## Regeneration and recovery

Multi-state NU-evidence envelope generation is supported. The implementation
regenerates envelopes when they are missing or incompatible, at `startit`, and
every `AMSK_FREQ` iterations.

Regeneration overwrites the per-state envelope each cycle. The envelope is free
to **shrink** as resolution improves as well as to grow: there is no monotonic
grow-only guard and no `pct_signal > 50%` regeneration failsafe. Spherical NU
support keeps every omitted region observable to the evidence calculation, so a
domain that recovers reproducible signal re-enters the envelope on the next
regeneration. Incompatible box or sampling changes bootstrap a new envelope.

## Compatibility rules

A state mask is compatible only if:

- all three dimensions match the current cropped box
- its sampling matches the current cropped sampling

Dimension-only checks are not sufficient once stage-dependent rescaling or
recropping is in play.
