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

NU volume filtering is independent of `incrreslim`: the FSC does not cap the
NU candidate bank or shell extension. Bootstrap NU filtering likewise never
uses the generic parsed startup `lp` as a volume-filter ceiling because it is
not evidence about the resolution of supplied half maps.

It also supplies overridable defaults when the user has not provided them:

- `mkdir=yes`
- `center=no`
- `sigma_est=global`
- `combine_eo=no`
- `prob_inpl=yes`
- `nsample=25000`
- `autoscale=yes`
- `filt_mode=nonuniform`
- `nu_refine=yes`
- `automsk=yes`
- `envfsc=no`
- `keepvol=no`

The default `envfsc=no` uses the broad spherical FSC mask. With `envfsc=yes`,
volume assembly generates a density/Otsu envelope from the current merged half
maps, low-pass filtered at `envmsklp`, and uses it for phase-randomized FSC
correction and cFAR. `envmsklp` defaults to `ENVMSKLP_DEFAULT` (20 A), while
`amsklp` remains the separate NU-evidence smoothing scale. The NU-evidence
envelope never enters FSC correction. `envfsc` and `envmsklp` remain advanced
filter controls and are independent of the default `automsk=yes` NU-reference
policy.

With `automsk=yes`, matching references consume the lagged
`nu_envmask3D_stateNN.mrc`. The first iteration uses the spherical
matching-reference mask; once volume assembly has written a compatible NU
envelope, later references are solvent-flattened before projection without
changing the NU filter, FSC policy, particles, or matching low-pass. There is no
separate `envref` control.

If `filt_mode` is overridden to a non-NU mode, `automsk` must also be set to
`no`; other combinations are rejected.

## 3. Starting Reference

Explicit `vol1` takes precedence. If `vol1` is absent, `refine3D_auto` may use
the project `os_out` state-1 `vol` entry when the file exists and its native
box and sampling match the current run.

If no compatible starting volume is available, `refine3D_auto` runs a
`reconstruct3D` startup pass and uses `vol_state01.mrc` as the initial
reference.

When NU filtering is active and an existing initializer is used, the workflow
requires a compatible same-stem raw native even/odd pair. It accepts
`_unfil` half maps when present and otherwise uses the same-stem even/odd
half maps. If the raw pair is missing or incompatible, the workflow falls back
to startup reconstruction instead of trusting stale derived NU products.

When the raw pair is compatible, `refine3D_auto` generates fresh same-stem
`_nu_filt` bootstrap references before the first matcher pass. With
`nu_refine=yes`, that bootstrap may run the sequential shell challenger from
the finest populated base-bank label.

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
to [refine3D_multi_policy.md](refine3D_multi_policy.md), base `refine3D`, or
the ab initio workflows.
