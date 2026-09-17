# Automasking Policy

## Scope and public modes

This policy covers 3D refinement, staged `abinitio3D`, reconstruction workers,
and `volassemble`. The public refinement modes are:

- `automsk=no`: retain spherical reconstruction/reference support
- `automsk=yes`: use the conservative density/Otsu envelope
- `automsk=nu`: prefer the lagged NU-evidence envelope and use the density
  envelope whenever the NU artifact is absent, incompatible, empty, or has an
  invalid evidence null

`automsk=nu` requires `filt_mode=nonuniform|nonuniform_lpset` in public and
internal commands. Final native reconstruction preserves the refinement's NU
filter mode rather than pairing `automsk=nu` with `filt_mode=none`.
`automsk=tight` remains available to standalone density-mask utilities but is
rejected in 3D refinement. Active refinement automasking implies `envfsc=yes`.

The density envelope is made from the even/odd average low-pass filtered at
`envmsklp`, followed by non-tight Otsu segmentation, largest-component
selection, at least `binwidth` dilation, and an `edge` cosine skirt. Its
dilation has the shared physical minimum `ENVMSKWIDTH_A_MIN` (7.5 A).

The NU-evidence envelope is derived from cross-half NU evidence at `amsklp`
and thresholded with `nu_msk_sig`. It remains state-local and may shrink or
grow on each regeneration.

## Lag and fallback

Early consumers in iteration N use `nu_envmask3D_stateNN.mrc` written by
iteration N-1. This intentional one-iteration lag applies to the gridding
FSC mask and to the matcher's reference fallback (PCG solve support is
always the density envelope, see below). The current NU filtering pass
derives a new evidence envelope and applies it immediately to the current
`_nu_filt` references before publishing it for the following iteration. If
generation does not produce a usable mask (empty evidence field, no
component kept, or an invalid null), the density envelope is used for the
current products and the stale NU artifact is removed, so a missing,
empty or incompatible artifact can never reactivate an old mask.

When no usable lagged NU mask exists, `automsk=nu` uses the conservative
density envelope. A PCG base solve with neither a NU artifact nor a prior
reference bootstraps on the sphere; its current base pair then supplies the
density fallback for the regularized replay. Invalid current evidence removes
the stale NU artifact and uses density for current NU products.

`automsk=nu` is an opt-in for particles whose NU-evidence envelope is
tighter than the density envelope and honest -- soluble particles without
ordered solvent or detergent. It is not a default and must not become one
without this record: the evidence envelope excludes density that the
particle images contain whenever that density is not reproducible between
halves at the evidence candidates -- a detergent micelle is the standard
case -- and a reference that no longer explains the micelle loses pose
discrimination under the euclid objective (PfCRT collapse, 2026-09-02;
`pcg_priors_history.md`). That is why the evidence envelope was demoted to a
diagnostic on 2026-09-13 and the conservative density envelope, which
retains every density present at `envmsklp`, became the reference mask of
`automsk=yes`.

## Reconstruction and reference use

On PCG, the conservative density envelope is installed as support for both
the base and regularized solves under `automsk=yes` and `automsk=nu` alike.
The mask therefore constrains the estimator and is not applied again to the
diagnostic half maps. The NU-evidence envelope is never a solve support
(review 2026-09-17): the evidence null of a constrained base pair is
designated on the density envelope's dilation ring at full support weight,
which lies outside an evidence envelope, so a NU-supported solve empties its
own null, invalidates the next envelope and falls back to a density envelope
derived from a map that is zero outside the NU support -- an oscillating
support with no way back for density the evidence excluded. An explicit
`pcg_mskfile` remains a development override and is reported as `explicit`
support.

On gridding, reconstruction retains spherical support. The selected envelope
is applied after reconstruction where required.

Matching references always receive the broad spherical soft mask first. The
automatic envelope is then applied after filtering:

- NU assembly applies the current selected envelope to `_nu_filt` even/odd
  products before they are used for reprojection
- matcher-owned non-NU filtering applies the density envelope for
  `automsk=yes`
- if NU products have not yet been assembled, `automsk=nu` reference fallback
  prefers the lagged NU artifact and then the density artifact

The NU objective domain itself always remains the spherical `mskdiam` support;
neither automatic envelope replaces that domain.

## FSC and phase randomization

Phase-randomized solvent correction is performed only when
`rec_backend=gridding` and a post-hoc envelope is selected. `automsk=yes` uses
the density envelope. `automsk=nu` uses the compatible lagged NU envelope, or
the density fallback.

PCG never performs phase randomization. It reports the FSC of the half maps
estimated on their installed support. A shared mask can still bias that FSC;
the method is therefore reported explicitly rather than presented as a
solvent-corrected curve.

Every resolution report and corresponding log line includes an `FSC MODE`
record with:

- reconstruction backend
- support kind (`sphere`, `density`, `explicit`, or `mixed`; never `nu`)
- post-hoc mask kind (`none`, `density`, or `nu`)
- whether phase randomization was applied

For gridding phase randomization, the onset is the first shell where the
genuinely unmasked FSC is below 0.8. Both halves are independently randomized
beyond that shell, and correction begins two shells later. If no usable onset
exists, the unmasked curve is retained.

## NU evidence validity

A spherical base pair estimates the evidence null robustly over the observed
support and requires solvent to remain the majority. A support-constrained PCG
pair uses the density envelope's dilation ring, at full support weight, as
its Euclidean null shell -- which is why the solve support must be the
density envelope in every mode. An empty or insufficient null makes the NU
mask invalid and selects the density fallback.

## Artifacts and compatibility

- `automask3D_stateNN.mrc`: conservative density envelope and fallback
- `nu_envmask3D_stateNN.mrc`: NU-evidence envelope
- `*_pcg_support.txt`: lagged solve-support provenance

A state mask is compatible only when all three dimensions match the current
cropped box and its sampling matches `smpd_crop`. Dimension-only checks are not
sufficient across stage-dependent recropping or rescaling.

Final reconstruction metadata prefers the NU artifact in `automsk=nu` and
falls back to the density artifact. Compatible non-PCG final postprocessing may
reuse that selected mask; PCG final maps are not masked again after the solve.
Both final-reconstruction routes (direct, and the sigma-bootstrap route taken
when the registration box differs from the native box) keep the refinement's
NU `filt_mode` under `automsk=nu`, since the mode is valid only there; under
`automsk=yes` the shipped map is classical (`filt_mode=none`).

## Validation criteria (`automsk=nu`, 2026-09-17)

- UI and parameter parsing expose `no|yes|nu` for the 3D refinement
  programs and their internal reconstruction/assembly commands only;
  `automsk=nu` with a non-NU `filt_mode` is rejected.
- Gridding: NU -> density fallback for the post-hoc FSC mask, with
  phase-randomized correction.
- PCG: density support for both solves in `yes` and `nu`; sphere only for
  the bootstrap without a density source; no phase randomization.
- Resolution text and logs state the method actually used (`FSC MODE`).
- Current NU products carry the current valid NU-evidence envelope,
  otherwise the density fallback; the choice is logged per state.
- Shared and distributed PCG paths preserve identical policy.
- The mode is opt-in; the PfCRT collapse record above is the reason.
