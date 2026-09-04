# Automasking Policy

## Scope

This document describes 3D envelope-mask policy for `refine3D`, staged
`abinitio3D`, and the shared `volassemble` step. Standalone density masking
commands remain available. Volume assembly can also generate a density/Otsu
mask for FSC correction independently of refinement-reference automasking.

The current architecture has two independent controls and artifacts:

- `automsk` requests a NU-evidence-derived per-state envelope that defines
  the filter-field background: outside it the NU filter takes the coarsest
  bank candidate, so matching references carry the excluded density heavily
  low-pass filtered (never removed), and the same field derives the Q_NU
  precisions on the pcg backend
- `envfsc` requests an on-the-fly density-derived per-state envelope for
  phase-randomized FSC correction and cFAR. Non-PCG derived final maps may
  reuse it; PCG maps are never masked after the solve.

Volume assembly owns both artifacts on the gridding path; the PCG strategy
owns their corresponding diagnostics, filter-field constraint, and solve
support. Neither mask may replace the broad spherical support used by the NU
objective.

## Public policy

The user-facing control is:

- `automsk=no`: no refinement-loop NU-evidence envelope generation
- `automsk=yes`: generate a per-state envelope from NU cross-half evidence
- `automsk=tight`: valid for standalone density-mask commands, but rejected in
  NU refinement because the NU-evidence envelope has no Otsu tight variant

In 3D refinement, `automsk=yes` requires
`filt_mode=nonuniform|nonuniform_lpset`; `automsk=tight` and non-NU filtering
with `automsk != no` are rejected. Standalone density-mask utilities are
outside this refinement invariant.

`envfsc=no` is the general default; `refine3D_auto` defaults it to `yes` unless
the user supplies a value. With `envfsc=yes`, volume assembly averages the
current even/odd half maps, low-pass filters the merged map at `envmsklp`,
applies non-tight Otsu segmentation, retains the largest connected component,
grows it by `binwidth`, and applies a cosine edge of width `edge`. `envmsklp`
defaults to `ENVMSKLP_DEFAULT` (20 A) and must be positive when `envfsc=yes`.
For `abinitio3D`, `binwidth` defaults to `ENVMSKWIDTH_DEFAULT` (7 voxels) and
is preserved in staged and final reconstruction commands.
The resulting density envelope is used for masked and randomized-masked FSC,
cFAR, and compatible final-map postprocessing. This path is independent of
`automsk` and NU filtering; `amsklp` remains the NU-evidence/standalone
automasking scale and does not control this density FSC envelope.

Matching references are never multiplied with an envelope before reprojection
(2026-09-02): hard-removing density that is present in the particle images
(e.g. a detergent micelle) destroys pose discrimination under the euclid
objective. With `automsk=yes` the down-weighting happens through the NU
filter field instead -- the background, defined as the complement of the NU
evidence envelope derived in the same evidence pass, takes the coarsest bank
candidate (cisTEM-style heavy background low-pass). The matcher applies the
spherical soft reference mask only; there is no separate `envref` control.
Particle images, FSC estimation, and matching-bandwidth selection are
unchanged by `automsk`.

NU filtering always uses the spherical support derived from `mskdiam`.
Envelope masks do not define or restrict the NU objective domain. In
particular, the correlation-derived NU envelope must never feed FSC correction
or replace spherical NU support.

`mskfile` is no longer part of the CLI policy. Passing `mskfile` is a hard error.

## Ownership

### Volume assembly

On the gridding backend, volume assembly owns both mask-production paths:

- half-map restoration generates `automask3D_stateNN.mrc` on every active
  `envfsc=yes` calculation, directly from the current half maps
- NU postprocessing decides whether `nu_envmask3D_stateNN.mrc` should exist
  and regenerates it according to the cadence below. With `automsk=yes`, it is
  deliberately derived from the static candidate bank before optimization and
  adaptive extension so it can constrain that same pass. A diagnostic envelope
  generated without arming the background may include accepted extensions.

The NU envelope is generated before `nu_filter_vols` releases the mask-packed
NU unary storage, and constrains the local filtering field in that same
assembly pass.

On the PCG backend, the reconstruction strategy builds the NU-evidence
envelope while the replay unaries are live and installs it as a fixed
coarsest-label boundary condition before constructing `Q_NU`. It independently
builds the conservative density mask used as solve support.

### FSC consumers

The NU-evidence envelope is never used for FSC correction because it is selected
from cross-half agreement. With `envfsc=no`, the reported radial FSC and cFAR
use the broad spherical FSC mask. With `envfsc=yes`, the density envelope is
passed to phase-randomized FSC correction and the same envelope is applied to
the cFAR copies.

The randomization onset is the first shell where the genuinely unmasked FSC is
below 0.8. The two half maps are independently phase-randomized beyond that
shell and then masked. The corrected curve starts two shells after the onset;
if no usable onset exists, the code retains the unmasked curve. The diagnostic
files are `fscu_stateNN.bin`, `fsct_stateNN.bin`, and `fscn_stateNN.bin`.

### Matcher reference preparation

Matcher reference preparation never reads or multiplies either envelope.
References receive only the broad spherical soft mask. With `automsk=yes`, the
NU-evidence envelope has already influenced the reference through the
coarsest-bank background assignment: through synthesized NU-filtered maps on
the gridding backend, or through the in-solve `Q_NU` precision on PCG.

## State-specific artifacts

The two state-specific artifacts are:

- `automask3D_stateNN.mrc`: current density/Otsu envelope produced by
  `envfsc=yes`; used by FSC/cFAR in memory and compatible non-PCG final
  postprocessing from disk
- `nu_envmask3D_stateNN.mrc`: NU-evidence envelope produced by `automsk=yes`;
  records the envelope that constrains the NU filter/prior background and its
  regeneration cadence; it is never applied directly to a reference

Both files are state-local. They have different statistical provenance and are
not interchangeable in the FSC or NU-objective paths.

## Current workflow

1. `refine3D` or staged `abinitio3D` produces partial reconstructions.
2. `volassemble` restores even and odd state volumes and calculates radial FSC
   and cFAR. `envfsc=yes` generates the density envelope during this step;
   `envfsc=no` uses the broad spherical FSC mask.
3. `volassemble` restores the merged state volume.
4. The NU filter constructs spherical support from `mskdiam` and evaluates the
   static candidate-bank unaries.
5. If `automsk=yes`, `volassemble` derives
   `nu_envmask3D_stateNN.mrc` from those live unaries and fixes the envelope
   background to the coarsest candidate; `automsk=no` leaves the spherical
   field unconstrained.
6. The NU filter optimizes the static field and accepts any supported
   `nu_refine` extensions inside that fixed background.
7. The NU envelope affects matching references only through the local filter
   field or `Q_NU`; it is never multiplied into a reference and never enters
   FSC correction or NU objective support.
8. Non-PCG final postprocessing may reuse a compatible
   `automask3D_stateNN.mrc` when `envfsc=yes`; PCG postprocessing applies no
   mask after the solve.

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
