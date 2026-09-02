# Refine3D Auto Policy

This document records the current policy for `refine3D_auto`. It describes the
automated wrapper around base `refine3D`; the base iteration contracts remain
in [refine3D_policy.md](refine3D_policy.md).

## 1. Scope

`refine3D_auto` is a single-state automated refinement workflow. It chooses
conservative defaults, prepares a starting reference when needed, runs base
`refine3D`, and then reconstructs a final all-particle map.

It is not a separate matcher implementation. Once startup material is ready,
the refinement iterations are delegated to `commander_refine3D`.

## 2. Defaults

`refine3D_auto` sets hard workflow defaults:

- `balance=no`
- `greedy_sampling=no`
- `trail_rec=yes`
- `refine=prob_neigh`
- `ml_reg=yes`
- `overlap=0.99`
- `nstates=1`
- `objfun=euclid`
- `lplim_crit=0.143`
- `incrreslim=no`

NU volume filtering is independent of `incrreslim`. On gridding the FSC does
not cap the NU candidate bank or shell extension. On PCG, where the field is
an empirical precision rather than a post-hoc filter, adaptive evidence-bank
discovery is bounded to two Fourier shells beyond the evidence pair's measured
FSC=0.143 crossing. Bootstrap NU filtering never uses the generic parsed
startup `lp` as a volume-filter ceiling because it is not evidence about the
resolution of supplied half maps.

It also supplies overridable defaults when the user has not provided them:

- `mkdir=yes`
- `center=no`
- `sigma_est=global`
- `combine_eo=no`
- `prob_inpl=yes`
- `nsample=25000`
- `autoscale=yes`
- `filt_mode=nonuniform`
- `nu_refine=yes` (gridding extends the filter bank; PCG extends the evidence
  bank used to construct `Q_NU`; the PCG extension is FSC-frontier bounded and
  spacing-aware)
- `automsk=yes`
- `envfsc=yes`
- `keepvol=no`

The default `envfsc=yes` is guarded, so an explicit user value remains
authoritative. It generates a density/Otsu envelope from the current merged
half maps, low-pass filtered at `envmsklp`, and uses it for phase-randomized
FSC correction and cFAR. On PCG it is also the conservative support for both
the base and regularized solves. `envmsklp` defaults to
`ENVMSKLP_DEFAULT` (20 A), while `amsklp` remains the separate NU-evidence
smoothing scale. The NU-evidence envelope never enters FSC correction or PCG
solve support.

With `automsk=yes`, the NU filter-field background is the complement of the
NU evidence envelope (`nu_envmask3D_stateNN.mrc`), derived from the static
candidate bank at the start of the same evidence pass. It is fixed before
adaptive candidates are challenged, so accepted probes refine precision inside
a causal boundary rather than redefining their own support. Background voxels take the
coarsest bank candidate, so matching references carry the envelope-excluded
density (detergent micelle, disordered belt) heavily low-pass filtered
rather than removed, and the same field derives the Q_NU precisions.
Matching references always take the spherical soft mask only -- they are
never multiplied with an envelope before projection, because a reference
must not hard-remove density that is present in the particle images. The
PCG solve support remains the conservative density envelope, never the
evidence envelope. There is no separate `envref` control.

If `filt_mode` is overridden to a non-NU mode, `automsk` must also be set to
`no`; other combinations are rejected.

## 3. Starting Reference

Explicit `vol1` takes precedence. If `vol1` is absent, `refine3D_auto` may use
the project `os_out` state-1 `vol` entry when the file exists and its native
box and sampling match the current run.

If no compatible starting volume is available, `refine3D_auto` runs a
`reconstruct3D` startup pass and uses `vol_state01.mrc` as the initial
reference.

For an explicit external `vol1`, `ref_pose_init=cc` selects the shared
external-reference transition. Before its fixed-reference CC pass, the service
runs the normal native-grid particle-image sigma2 bootstrap. In the initialized
cohort, the CC pass then replaces active matching shells with
reference-conditioned residuals and consolidates the resulting spectra for the
first Euclidean iteration.
Particles outside the capped cohort retain image-bootstrap records until they
are updated by refinement. `ref_pose_init=none` trusts the supplied reference
and does not invoke this wrapper-owned transition.

When NU filtering is active and an existing initializer is used, the workflow
requires a compatible same-stem raw native even/odd pair. It accepts
`_unfil` half maps when present and otherwise uses the same-stem even/odd
half maps. If the raw pair is missing or incompatible, the workflow falls back
to startup reconstruction instead of trusting stale derived NU products.

When the raw pair is compatible, `refine3D_auto` generates fresh same-stem
`_nu_filt` bootstrap references before the first matcher pass. With
`nu_refine=yes`, that bootstrap may run the sequential shell challenger from
the finest populated base-bank label.

Under `rec_backend=pcg`, no `_nu_filt` bootstrap reference is produced. The
startup reconstruction supplies the first `Q_NU`-regularized primary pair;
when an existing initializer is used, the first matching pass necessarily
precedes PCG evidence construction and uses that initializer with only the
spherical matcher support.

## 4. Autoscaling and Sampling

With `autoscale=yes`, native boxes larger than the minimum box are downscaled
toward the default target sampling of 1.3 A. Smaller boxes or `autoscale=no`
run at native sampling.

The translation search limit is derived from the active cropped sampling and
clamped to a practical range.

The workflow samples up to `nsample` active particles per iteration. If active
particles do not exceed `nsample`, update fraction is disabled and each
iteration is a full update. Otherwise `update_frac` is set to
`nsample / active_particles`.

Automatic iteration planning targets roughly four updates per active particle,
caps the run length, and enforces a minimum of ten iterations unless the user
explicitly supplied `maxits`.

## 5. Refinement and Final Reconstruction

After startup, `refine3D_auto` runs base `refine3D` with:

- `prg=refine3D`
- trailing reconstruction weighted from realized sampled-update bookkeeping
- the planned `maxits`
- the selected starting `vol1`

After refinement, it runs a final `reconstruct3D` pass from all particle
images. Final reconstruction sets `postprocess=yes`, sets `nu_refine=no`, and
turns `filt_mode` back to `none` when the refinement used NU filtering.

Final-map postprocessing is classical global FSC/B-factor postprocessing. NU
filtering is a refinement-reference feature, not a separate final-map
postprocess path.

## 6. Output Policy

The final reconstruction writes ordinary reconstruct3D products and then
`write_final_rec_outputs` records the final map products using the requested
resolution target.

`refine3D_auto` remains single-state. Multi-state automated refinement belongs
to [refine3D_states_policy.md](refine3D_states_policy.md), the in-development
[classify3D_refs_policy.md](classify3D_refs_policy.md), base `refine3D`, or the ab
initio workflows.
