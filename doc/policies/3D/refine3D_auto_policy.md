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

NU volume filtering is independent of `incrreslim`. Since 2026-09-16 there
is one NU competition for every workflow (the `nu_refine` shell walk is
retired): the candidate bank is a generated ladder of hard Butterworth
rungs -- the coarse rungs 20, 15, 12, 10, 8 and 6 A plus fine rungs every
`NU_LADDER_FINE_STEP` (2) Fourier shells, widened only when the
`NU_BANK_MAX_MEMBERS` (16) budget requires it -- bounded one shell coarser
than the finest member, which is the ML-regularized (closed-form Wiener)
pair whenever `ml_reg=yes`. Without a regularized pair the ladder is
bounded at `fsc/1.5` of the base pair (else Nyquist). The regularized pair
competes for voxels like any other member, so the FSC never gates the
filter directly; it enters only through the regularizer `P_tau`, the
bound, resolution reporting, convergence and the `envfsc` correction. The
matching low-pass handoff is the content extent of the finest selected
label with at least 1% of the signal voxels at that label or finer: a hard
rung hands off its cutoff, the regularized member hands off the pair's
FSC=0.143 resolution, with no headroom (`nonuniform_filtering_policy.md`
sections 8, 10, 12). Bootstrap NU filtering never uses the generic parsed
startup `lp` as a volume-filter ceiling because it is not evidence about
the resolution of supplied half maps.

It also supplies overridable defaults when the user has not provided them:

- `mkdir=yes`
- `center=no`
- `sigma_est=global`
- `combine_eo=no`
- `prob_inpl=yes`
- `nsample=25000`
- `autoscale=yes`
- `filt_mode=nonuniform`
- `automsk=yes`
- `envfsc=yes`
- `keepvol=no`
- `regpass=yes`, `regpass_fsc=0.8` (section 5)

The default `envfsc=yes` is guarded, so an explicit user value remains
authoritative -- except that `automsk=yes` implies `envfsc=yes` (policy
2026-09-09), so with the default `automsk=yes` envfsc cannot be switched off.
It generates a density/Otsu envelope from the current merged half maps,
low-pass filtered at `envmsklp` and dilated by at least `ENVMSKWIDTH_A_MIN`
(7.5 A), and uses it for phase-randomized FSC correction and cFAR on
gridding. On PCG it is the conservative support for both the base and
regularized solves, and the FSC is then reported on the constrained pair
without post-hoc correction (the `>>> FSC MODE` line says which). `envmsklp`
defaults to `ENVMSKLP_DEFAULT` (20 A), while `amsklp` remains the separate
NU-evidence smoothing scale. The NU-evidence envelope never enters FSC
correction or PCG solve support; on PCG its null is designated on the density
envelope's dilation ring rather than estimated from the constrained pair.

With `automsk=yes`, the conservative density envelope is the support of the
signal model (policy 2026-09-13). Outside it there is no signal to filter --
unreconstructed under the PCG support projection, solvent on gridding -- so
the filter field takes the coarsest bank candidate there, fixed before
adaptive candidates are challenged, and the `_nu_filt` matching references
are multiplied by the envelope after filtering. The NU evidence envelope
(`nu_envmask3D_stateNN.mrc`) is still derived from the static candidate bank
at the start of the same evidence pass, as a diagnostic of the evidence
field; it is never armed and never multiplied into a reference, because it
cuts out detergent density that the particle images contain. The matcher
itself applies the spherical soft mask only. The PCG solve support remains
the conservative density envelope, never the evidence envelope. There is no
separate `envref` control.

`filt_mode` may be overridden to a non-NU mode with `automsk=yes` kept
(2026-09-14): the density envelope then reaches the references through the
matcher instead of the `_nu_filt` products, and the matching low-pass comes
from the FSC at `lplim_crit` instead of the NU handoff. Note that with
`ml_reg=yes` and a non-NU mode the references are the shipped ML-regularized
pair, unfiltered by the matcher.

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
`_nu_filt` bootstrap references before the first matcher pass, from the
same generated ladder as every later iteration (the regularized member
requires `ml_reg=yes`).

Under `rec_backend=pcg` the startup reconstruction runs the same NU
competition inside the PCG master and produces the same `_nu_filt` bootstrap
references and matching-lp handoff as gridding (policy 2026-09-06). The
initial volume (external `vol1` or the compatible project volume) is passed
to it as `vol1` (2026-09-11), so its PCG base pair is envelope-constrained
from the start instead of bootstrapping on the sphere.

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
caps the run length, and enforces a minimum of three iterations unless the user
explicitly supplied `maxits` (or a larger `minits`). The minimum was ten until
2026-09-11; on PfCRT the forced iterations degraded an already converged map
(cFAR 0.78 to 0.62, FSC=0.5 4.14 to 4.31 A over iterations 1-10). Once the
overlap criterion (0.99) is met after the third iteration the run stops.

## 5. Registration Pass, Refinement and Final Reconstruction

After startup, with `regpass=yes` (default), `refine3D_auto` runs one
global registration pass (2026-09-16): a single `refine3D` iteration with
`refine=prob`, `nspace=5000`, no subspace, every active particle regardless
of the sampling policy, matched against the masked startup references and
band-limited at the resolution where the startup pair's FSC falls below
`regpass_fsc` (default 0.8). The rationale: previous poses are a fixed point
of the previous objective at the full band; the solvent-mask constraint on
the references is a new objective, and a neighbourhood search at the full
band cannot leave the old basins (aldolase run 16: 0.998 orientation overlap
in iteration 1 and no motion after). Previous poses therefore only produce
the startup reference; the pass gives every particle a global search under
the constraint at a band where the reference is trustworthy. The band is
imposed with `lpstop`, never `lp`: an explicit `lp` sets `l_lpset`, which in
every non-NU `filt_mode` matches both halves against the merged reference
and would silently break gold standard. A user `lpstop` caps the pass as
well (the coarser of the two applies) and is restored for the main run. The
pass is skipped, logged, when no startup FSC exists or it never reaches
`regpass_fsc`. It logs `REGISTRATION PASS REASSIGNED`: the fraction of
directions that moved by more than the orientational basin width at the
pass band (`res / (mskdiam/2)`), by more than twice it, and the fraction of
shifts that moved by more than one pixel -- the number that says whether
re-basining happened. The main run then continues from the pass output as
iteration 2. `regpass_fsc` is the knob for testing coarser bands.

The main run is base `refine3D` with:

- `prg=refine3D`
- `refine=prob_neigh`, `nspace=20000`, `nspace_sub=500` unless given
- trailing reconstruction weighted from realized sampled-update bookkeeping
- the planned `maxits`, counted from iteration 2 after a registration pass
- the selected starting `vol1` (the pass output when the pass ran)

After refinement, it runs a final `reconstruct3D` pass from all particle
images. Final reconstruction sets `postprocess=yes` and
turns `filt_mode` back to `none` when the refinement used NU filtering.
`automsk` is inherited (2026-09-09): on PCG the shipped map is estimated on
the same density-envelope support as every refinement iteration, with the
same implied `envfsc=yes` and the same reported `>>> FSC MODE`, rather than
falling back to the sphere for the map that matters most. The PCG support
is derived from a reference volume (`build_pcg_state_support`), so the
final `reconstruct3D` receives `bootstrap_rec3D`'s own gridding bootstrap
map as `vol<state>` (2026-09-11); without it the base pair bootstrapped on
the sphere and the resolution doc reported the gridding wording, "density
envelope applied post hoc ... phase-randomized correction", for a PCG map
(bgal). The base pair and the reported FSC are now envelope-constrained.

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
