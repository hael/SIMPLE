# Automasking Policy

## Scope

This document describes 3D envelope-mask policy for `refine3D`, staged
`abinitio3D`, and the shared `volassemble` step. Standalone density masking
commands remain available. Volume assembly can also generate a density/Otsu
mask for FSC correction independently of refinement-reference automasking.

The current architecture has two controls and two artifacts, coupled in one
direction (policy 2026-09-09):

- `automsk` requests a NU-evidence-derived per-state envelope that defines
  the filter-field background: outside it the NU filter takes the coarsest
  bank candidate, so matching references carry the excluded density heavily
  low-pass filtered (never removed), on both backends. It is the only
  envelope that excludes detergent.
- `envfsc` requests an on-the-fly density-derived per-state envelope (the
  conservative density envelope). On gridding it is applied post hoc to the
  FSC pair with phase-randomized solvent correction and to the cFAR copies;
  on PCG it is the solve support of both the base and the regularized solve,
  so the FSC pair is envelope-constrained inside the estimator. Non-PCG
  derived final maps may reuse it; PCG maps are never masked after the solve.
- `automsk=yes` implies `envfsc=yes` on both backends. The density envelope
  is derived in parameter validation, never requested separately, because a
  density-constrained estimate is the envfsc contract by construction.

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
the user supplies a value, and `automsk=yes` promotes it to `yes` everywhere
(logged when it overrides an explicit or defaulted `no`). With `envfsc=yes`,
the density envelope is built from the current even/odd average: low-pass at
`envmsklp`, non-tight Otsu segmentation, largest connected component, spherical
dilation by `binwidth` layers, and an outward cosine skirt of width `edge`.
`envmsklp` defaults to `ENVMSKLP_DEFAULT` (20 A) and must be positive when
`envfsc=yes`. The dilation has a shared physical minimum, `ENVMSKWIDTH_A_MIN`
(7.5 A, the former abinitio3D default of 7 layers at 1.075 A/pixel): whenever
the envelope is in use, `binwidth = max(binwidth, ceiling(7.5 A / smpd_crop))`
at the sampling the envelope is built at, so the same physical envelope comes
out at every crop level and in every program; an explicit larger `binwidth`
wins. The density envelope is used for masked and randomized-masked FSC, cFAR,
the PCG solve support, the Euclidean null shell of the NU evidence envelope,
and compatible final-map postprocessing. `amsklp` remains the
NU-evidence/standalone automasking scale and does not control this density
envelope.

Matching references are never multiplied with an envelope before reprojection
(2026-09-02): hard-removing density that is present in the particle images
(e.g. a detergent micelle) destroys pose discrimination under the euclid
objective. With `automsk=yes` the down-weighting happens through the NU
filter field instead -- the background, defined as the complement of the NU
evidence envelope derived in the same evidence pass, takes the coarsest bank
candidate (cisTEM-style heavy background low-pass). The matcher applies the
spherical soft reference mask only; there is no separate `envref` control.
Particle images and matching-bandwidth selection are unchanged by `automsk`;
FSC estimation follows the implied `envfsc=yes` (see below).

NU filtering always uses the spherical support derived from `mskdiam`.
Envelope masks do not define or restrict the NU objective domain. In
particular, the correlation-derived NU envelope must never feed FSC correction
or replace spherical NU support.

The NU evidence envelope's null model has two regimes (policy 2026-09-09),
keyed on how the base pair was solved:

- spherical base pair (gridding; PCG bootstrap without a lag-one reference):
  the robust median + `nu_msk_sig` MAD of the margin over the observed
  support, where the generous sphere makes solvent the majority; validity is
  the solvent majority
- envelope-constrained base pair (PCG under `automsk=yes`, or an explicit
  `pcg_mskfile`): the estimator has removed the far solvent, so the null is
  designated by Euclidean geometry instead of estimated from a mixture -- the
  median/MAD are taken on the density envelope's dilation ring (dilated minus
  core) at full weight of the base support; labels are free on the observed
  density envelope and fixed solvent outside it, nesting the evidence envelope
  inside the density envelope; validity is shell sufficiency

If the null is invalid or the envelope is empty, the density envelope itself
is armed as the filter-field background (logged as `EVIDENCE FALLBACK`); the
provenance string records which envelope ran.

`mskfile` is no longer part of the CLI policy. Passing `mskfile` is a hard error.

## Ownership

### Volume assembly

On the gridding backend, volume assembly owns both mask-production paths:

- half-map restoration generates `automask3D_stateNN.mrc` on every active
  `envfsc=yes` calculation, directly from the current half maps
- NU postprocessing regenerates `nu_envmask3D_stateNN.mrc` on every
  `automsk=yes` competition. It is deliberately derived from the static
  candidate bank before optimization and adaptive extension so it can
  constrain that same pass. A diagnostic envelope generated without arming the
  background may include accepted extensions.

The NU envelope is generated before `nu_filter_vols` releases the mask-packed
NU unary storage, and constrains the local filtering field in that same
assembly pass.

On the PCG backend, the PCG master runs the same assembly-owned NU competition
(`simple_nu_state_filter`), so the NU-evidence envelope is produced and
consumed exactly as on gridding, and hands over the support that constrained
the base pair so the evidence null takes the Euclidean-shell regime. It
independently builds the conservative density mask used as solve support, but
only under `automsk=yes`; with `automsk=no` no density mask is built and every
PCG solve runs on the spherical support (policy 2026-09-06). Under
`automsk=yes` the density envelope constrains both the base and the
regularized solve once a prior reconstruction exists (policy 2026-09-09).

### FSC consumers

The NU-evidence envelope is never used for FSC correction because it is selected
from cross-half agreement. With `envfsc=no`, the reported radial FSC and cFAR
are computed on the shipped half-maps, which carry the soft spherical support
at `msk_crop` from the reconstruction itself (no second mask, 2026-09-09). With
`envfsc=yes` on gridding, the density envelope is passed to phase-randomized
FSC correction and the same envelope is applied to the cFAR copies. With
`envfsc=yes` on a support-constrained pair (PCG under `automsk=yes`, or
`pcg_mskfile`), no post-hoc mask and no phase-randomized correction are
applied: the envelope is already inside the estimate. This is a deliberate
choice, not a claim that a constrained estimate is bias-free -- a common
window on both halves can still contribute correlated power -- and it is
therefore REPORTED rather than hidden: every diagnostic evaluation logs
`>>> FSC MODE` and writes the same line into the resolution text, naming one
of three modes (spherical support, envelope post hoc with solvent correction,
estimator-constrained without correction).

The randomization onset is the first shell where the genuinely unmasked FSC is
below 0.8. The two half maps are independently phase-randomized beyond that
shell and then masked. The corrected curve starts two shells after the onset;
if no usable onset exists, the code retains the unmasked curve. The diagnostic
files are `fscu_stateNN.bin`, `fsct_stateNN.bin`, and `fscn_stateNN.bin`.

### Matcher reference preparation

Matcher reference preparation never reads or multiplies either envelope.
References receive only the broad spherical soft mask. With `automsk=yes`, the
NU-evidence envelope has already influenced the reference through the
coarsest-bank background assignment in the synthesized NU-filtered maps, on
both backends.

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
   `envfsc=no` evaluates the shipped, support-masked halves without a second
   mask.
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
   field; it is never multiplied into a reference and never enters
   FSC correction or NU objective support.
8. Non-PCG final postprocessing may reuse a compatible
   `automask3D_stateNN.mrc` when `envfsc=yes`; PCG postprocessing applies no
   mask after the solve.

## Regeneration and recovery

Multi-state NU-evidence envelope generation is supported. The envelope is
regenerated from the live evidence on every NU competition under
`automsk=yes`; there is no cadence, and the artifact on disk has no
in-workflow reader (the former `AMSK_FREQ` planner was removed, 2026-09-09).

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
