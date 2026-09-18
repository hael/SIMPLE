# Nonuniform Filtering Policy

This document records the current nonuniform-filtering contract. It focuses on
active code paths and should stay aligned with what the implementation does
today.

## 1. Scope

Nonuniform filtering is a volume-domain regularization feature. On both
reconstruction backends it selects a local low-pass limit inside spherical
`mskdiam` support and writes NU-filtered derived references
(`simple_nu_state_filter`, shared by the gridding volassemble and the PCG
master since 2026-09-06). Its candidate bank is not truncated by the FSC. The
finest selected label separately governs the matching bandwidth handed to
later iterations.

It is not a separate final-map postprocessing workflow. `postprocess` and the
automatic `reconstruct3D` postprocess step use the ordinary global
FSC/B-factor path, even when the reconstruction that produced their input
volumes used a NU `filt_mode`.

## 2. User Controls

The public filter selector is `filt_mode`.

Supported values:

- `none`
- `uniform`
- `fsc`
- `nonuniform`
- `nonuniform_lpset`

`filt_mode=nonuniform` enables NU-filtered volume products. After the first
iteration, the finest cutoff selected by the NU filter supplies the matching
bandwidth while the matcher retains independent half-map topology.

`filt_mode=nonuniform_lpset` enables the same NU filter and promotes the
selected NU bandwidth into an LP-set matching run. LP-set matching uses merged
registration-reference topology.

On `rec_backend=pcg` (2026-09-06) both NU modes run the same post-hoc NU
competition as gridding inside the PCG master: the unregularized `_unfil` pair
seeds the candidate bank, the `P_tau`-regularized pair is the finest member
of the competition (`ml_reg=yes`), and the `_nu_filt`/`_nu_locres`
products and the matching-lp handoff are written exactly as on gridding. The
former in-solve `Q_NU` replay precision and its controllers were removed.

The competition input is the unregularized (base) pair on both backends.
Two alternatives were tried and retired on 2026-09-09 (`nu_input=gridding|ml`,
records 2026-09-08b-g and 2026-09-09 in
`doc/implementation_notes/pcg_priors_history.md`): the gridding half of the PCG
accumulation became moot once the like-for-like selection removed the
footprint artifact, and the ML-regularized pair cannot seed the competition
because `P_tau` is a global per-shell shrinkage driven by the global FSC: it
removes exactly the local content beyond the global FSC that the competition
must find to advance the band, and the under-converged replay compounds the
shrinkage across iterations (PfCRT: labels drifted coarser every iteration,
the sampler settled on the smooth reference, early stopping fired, three runs
pinned at 8.4-9.1 A).

There is one NU competition (2026-09-16): the former `nu_refine` control
and the high-resolution shell walk it enabled are retired (section 10). The
bank is the coarse ladder plus fine rungs generated from the box (section
8), the same in every workflow and in `postprocess_nu`.

`mskdiam` controls the spherical NU support mask. `automsk` separately controls
NU-evidence envelope generation, but is valid only while NU filtering is active.
`ml_reg` provides the regularized even/odd pair, the finest member of the
competition.

Refinement automasking is independent of the filter mode (2026-09-14):
`automsk=yes` multiplies the matching references by the conservative density
envelope in every `filt_mode` -- at assembly on the `_nu_filt` products in NU
modes, in the matcher after its own filter otherwise
(`automasking_policy.md`). `automsk=tight` is rejected in 3D refinement.
`automsk=nu` requires NU filtering, applies a valid current evidence envelope
to `_nu_filt` products, and uses its lagged artifact for early consumers;
density is the fallback.

`envfsc=no` is the general default and the broad-sphere FSC remains the
reported curve. `refine3D_auto` overrides that default to `yes` unless the user
sets it explicitly.
With `envfsc=yes`, volume assembly generates a density envelope on the fly from
the merged current half maps, smoothing it at `envmsklp`. Gridding performs
phase-randomized solvent correction with the selected density or NU mask; PCG
never phase-randomizes. `envmsklp` defaults to 20 A through
`ENVMSKLP_DEFAULT`; it is separate from `amsklp`, which continues to set the NU
evidence smoothing scale. FSC
correction changes reported resolution metadata; it does not truncate the NU
filter bank or directly set the NU matching bandwidth. `envfsc` is independent
of `automsk` and can operate without NU filtering.

With `automsk=yes` the conservative density envelope is the
support of the signal model (policy 2026-09-13): outside it there is no signal
to filter -- unreconstructed under the PCG support projection, solvent on
gridding -- so the filter field takes the coarsest bank candidate there, and
the `_nu_filt` references are multiplied by the envelope after filtering, on
both backends. With `automsk=nu`, the current valid NU evidence envelope is
armed instead and density is used when it is invalid or unavailable. With
`automsk=no`, the entire spherical `mskdiam` support remains unconstrained.
The former `envref`
parameter has been removed.

## 3. Ownership

`simple_vol_pproc_policy.f90` owns per-state automask regeneration/reuse
decisions. NU support is fixed by the `simple_nu_filter` API.

`commander_volassemble` in `simple_commanders_rec_distr.f90` owns execution:
restoring half maps, planning postprocessing, running automask generation,
running NU filtering, writing derived products, and recording NU matching
bandwidth metadata.

`simple_nu_filter` owns the filter algorithm and its module-level working
state: candidate-bank setup, objective generation, ordered-label smoothing,
optional high-resolution extension, output synthesis, diagnostics, and cleanup.

`simple_matcher_refvol_utils.f90` owns matcher reference loading. It decides
whether to use NU-filtered even/odd references, a merged NU reference, or a
regular fallback reference.

`simple_matcher_smpl_and_lplims.f90` and
`simple_refine3D_strategy.f90` own the selected-LP handoff into matching.

The standalone `nu_filt3D` program also uses `simple_nu_filter`, but it is an
explicit filtering command rather than a workflow policy layer.

## 4. Volume Assembly Contract

Workflow NU filtering runs in Cartesian `volassemble` after the state
half-maps and merged map have been restored and after provisional
FSC/resolution metadata for the state has been calculated.

For each state, `volassemble` then:

1. calculates the radial FSC and cFAR, generating a density envelope and doing
   phase-randomized correction when `envfsc=yes`
2. restores and writes the merged/base state volumes
3. plans the NU-evidence envelope action
4. configures spherical NU support from `mskdiam`
5. configures the full static NU candidate bank
6. with active automasking, derives the NU-evidence envelope from the live
   static unaries and fixes the background outside the selected envelope to
   the coarsest candidate; `nu` selects valid evidence with density fallback
7. optimizes the local filter map
8. with active automasking, multiplies the NU-filtered even and odd references
   by the selected envelope, then writes NU-filtered even, odd, merged, and
   local-resolution products
9. records the finest locally selected NU low-pass limit for later handoff

Low-resolution even/odd insertion is a registration-reference preparation
trick. It must not feed `volassemble` FSC calculation, automasking, NU
filtering, ordinary half-map handoff, or on-disk half-map products.

When trailing reconstruction is active, the trailing blend is applied to the
restored half maps used by automasking and to the NU base/auxiliary inputs
before NU filtering.

## 5. Inputs

The NU filter consumes:

- the current unfiltered even volume
- the current unfiltered odd volume
- a spherical support diameter in Angstrom
- optionally, an auxiliary even/odd pair (the regularized pair), the finest
  member of the bank
- optionally, the auxiliary pair's effective resolution in Angstrom
- optionally, the base pair's FSC=0.143 resolution, which bounds the hard
  rungs when no auxiliary pair is supplied

When `ml_reg=yes`, `volassemble` uses the `_unfil` even/odd pair as the base NU
input and passes the ML-regularized even/odd pair as the auxiliary member.
The auxiliary effective resolution comes from the state FSC(0.143) resolution,
`res0143s(state)` (clamped by a set `lp` in `lpset` topology).

On `rec_backend=pcg` the base input is the unregularized solve pair.

On both backends the input halves are deapodized and carry the soft spherical
support at `msk_crop` exactly once (2026-09-09): the PCG solve support, and on
gridding the identical `mask3D_soft` applied after deapodization in
`restore_gridding_pair`. The NU machinery builds its own logical sphere from
`mskdiam` and does not mask the inputs again.

## 6. Spherical Support Contract

All NU entry paths use a spherical support mask derived from `mskdiam`.
`setup_nu_dmats` constructs that mask internally; callers cannot supply an
arbitrary logical envelope. This prevents density- or correlation-conditioned
masks from changing the normalized Huber objective domain and keeps a broad
solvent population available for future NU-evidence null estimation.

Spherical geometry alone does not guarantee a valid solvent-majority null.
`mskdiam` must be generous enough to include substantial solvent around the
particle. If a NU-evidence segmentation reports more than half of support as
signal, workflow integration must reject that envelope or use a statistically
different null estimator; a warning alone is not sufficient for automatic
reference masking.

The cost is memory: persistent objective storage scales as
`n_support_voxels * n_candidates`. Any future proposal to reduce support with a
dilated envelope must include a replacement null estimator, temporal recovery
guards, and a measured memory justification. It must not silently weaken this
API invariant.

Envelope generation and compatibility remain separate from NU support.
Standalone `nu_filt3D` therefore exposes `mskdiam`, not `automsk`, for NU
support.

The Huber unary is WHITENED by a radially-resolved raw E/O noise profile
(`image::nu_objective_noise_profile`: shell-wise Gaussian-scaled MAD of the
raw even-odd difference over real-space radius, gap-filled and smoothed, with
per-voxel linear interpolation between shell centres). Reconstruction noise is
not spatially stationary — deapodization amplifies the periphery and solve
supports taper it — and the earlier single global scale put peripheral
residuals in the wrong Huber regime, compressing their cost-improvement
margins and biasing both the filter competition and the evidence envelope
toward the centre. (That flaw was historically masked by the gridding
under-deapodization bug, whose radial fade approximately cancelled the true
sigma(r) rise; fixing deapodization exposed it as over-tight envelopes.
Measured on the neutral fixture: sigma(r) edge/centre 1.29; whitening raised
envelope recall of true density from 0.48 to 0.61 at unchanged component
count.) `>>> NU WHITENING PROFILE` reports shells, min/max and edge/centre
ratio at every setup.

When standalone NU-evidence envelope generation is enabled, its public shape
controls are limited to `nu_msk_sig` (robust evidence threshold) and `amsklp`
(physical evidence scale, in Angstrom). Production fixes the evidence form to
the radially-whitened Huber-cost margin, density weight to zero, MRF
beta to 1, and minimum component fraction to 0.1. It also fixes binary growth
at 1 A and the cosine edge at 6 A; `nu_filt3D` converts those physical lengths
to the nearest voxel counts at the input-map sampling, with a one-voxel
minimum. These values reproduce the prior 1-pixel and 6-pixel defaults at the
1 A/pixel reference sampling without creating additional public tuning knobs.
Their meanings remain part of the contract: beta controls boundary smoothness;
a positive density weight would retain strong but poorly ordered density; the
component fraction removes components smaller than that fraction of the
largest; and scale-free evidence would protect weak ordered density from being
outvoted by a high-contrast core.

Mask ownership is strict:

- spherical `mskdiam` support defines the NU objective domain
- `nu_envmask3D_stateNN.mrc` is generated from NU evidence under active
  automasking; `automsk=nu` uses it as the active envelope and lagged PCG/FSC
  mask when valid
- density-derived `automask3D_stateNN.mrc` is generated during volume assembly
  when `envfsc=yes`; the same generated mask feeds phase-randomized FSC and
  cFAR, is applied to the `_nu_filt` references under `automsk=yes`
  (2026-09-13), and may be reused by non-PCG final postprocessing
- PCG uses density support for `automsk=yes`; `automsk=nu` prefers the lagged
  NU artifact and falls back through density to a spherical bootstrap
- neither automatic envelope replaces the spherical NU objective support

## 7. Outputs

Workflow filtering writes derived products beside the primary reconstruction
outputs:

- `vol_state_even_nu_filt.mrc`
- `vol_state_odd_nu_filt.mrc`
- `vol_state_nu_filt.mrc`
- `vol_state_nu_locres.mrc`

Actual names append `NUFILT_SUFFIX`, currently `_nu_filt`, to the even, odd,
and merged state volume names. The local-resolution map appends
`NULOCRES_SUFFIX`, currently `_nu_locres`, to the merged state volume name.

The `_nu_locres` map stores the resolutions in Angstroms. Voxels outside the NU
support mask, and values above Nyquist, are written as zero.

Base even/odd and merged volumes remain the primary reconstruction outputs.
NU-filtered products are derived references and diagnostics.

## 8. Filter Algorithm

The current filter performs these steps:

1. build a retained low-pass bank from the base even/odd pair
2. optionally replace the finest discrete bank member with an auxiliary pair
3. cache low-pass-filtered bank volumes as local scratch files
4. compute mask-packed unary objective costs for retained candidates
5. smooth each candidate objective over a mask-normalized local support
6. select the label per in-mask voxel coarse to fine: a finer candidate
   replaces the incumbent only if it wins at its own smoothing scale with
   both smoothed alike (2026-09-08; see section 9)
7. apply ordered-label Potts smoothing to the candidate map, which never
   promotes a voxel beyond the label it entered with
8. synthesize filtered even/odd outputs from the selected labels
9. write the merged `_nu_filt` output as the even/odd average
10. write the same-grid `_nu_locres` map

The bank is the static ladder `[20, 15, 12, 10, 8, 6, 5, 4]` A (the
machinery of commit ed36eb4c's abinitio3D, the only NU mechanism since
2026-09-18; the generated dense ladder of 2026-09-16 and the shell walk of
section 10 are both gone). Given the pair's FSC=0.143 resolution the bank
is capped (2026-09-08): only candidates coarser than `fsc/1.5`
(`NU_BANK_FSC_HEADROOM`, about two ladder labels finer than the FSC) are
retained, never fewer than two. The unary prices a finer candidate by the
noise it admits from the other half, which holds for a gridding pair but
not for a spectrally regularized one: on PfCRT a PCG base pair populated
the finest label at a 6 A FSC and pinned the matching band there. The July
gridding runs populated one to three labels beyond the FSC (4.7/1.8/0.6% of
the sphere at 6.3 A), which the cap admits. With `ml_reg=yes` the
ML-regularized pair is appended as one more member of the bank BESIDE the
finest retained rung (2026-09-18; until then it took that rung's place):
the two compete voxel by voxel under the same unary, smoothing and prior,
and the auxiliary shares the finest rung's Potts coordinate, so replacing
a finest-rung voxel by the regularized pair costs the prior nothing -- the
unary alone decides, and the regularized pair is included exactly where it
wins. Its effective resolution is the pair's FSC=0.143 (clamped by a set
`lp`); its label hands off that resolution. The master logs one `>>> NU
BANK CAP:` line and one `>>> NU AUXILIARY MEMBER:` line per state. Absent
an FSC (the standalone `nu_filt3D` program) the bank is uncapped.

An opt-in replay-evidence API can compact this full unary bank before it is
released. Callers must tag the setup source as `base_unfil`; the API fingerprints
and rechecks the exact half pair and rejects the ML auxiliary-replacement path.
It adds a zero cross-half-prediction null to a separate ordered-label model.
Because raw zero prediction has a systematic Huber-loss offset relative to a
smoothed predictor even for independent noise, and selecting the best of several
signal candidates adds a multiple-comparison advantage, the null score
subtracts the robust median-plus-three-MAD offset of
`C_zero-min(C_signal bank)` over the generous spherical support. This calibrates
the actual competing bank while retaining sensitivity to genuinely coarse
shared signal whenever its candidate wins, rather than treating the 20-A label
as solvent. The API then freezes selected
cutoff, normalized label entropy, and nested support confidence through
20/12/8/5 A plus the spherical-support geometry in
`nu_evidence_state`. This evidence analysis
does not alter the NU filtering label map or outputs.
`expand_nu_evidence_band_weights` expands the frozen state into per-band
lack-of-evidence weight fields (`1 - a_b` inside the spherical evidence
support, 1 outside it), recreating the packed lexicographic order from the
frozen geometry alone so it works after `cleanup_nu_filter`. Before any
replay use, `assert_nu_evidence_replay_ready` enforces the readiness
contract: a state whose explicit null wins less than
`NU_EVIDENCE_MIN_NULL_FRAC` or more than `NU_EVIDENCE_MAX_NULL_FRAC` of the
OBSERVED part of the generous spherical support marks a failed null
calibration (starved and saturated null respectively) and hard-errors --
validity alone does not qualify evidence to parameterize a precision. The
observed part excludes exact zero/zero voxels that a density-constrained PCG
solve leaves inside the sphere (`nu_observed_mask`, set by `setup_nu_dmats`
with the same test as the whitening profile); every calibration statistic
(null-bias center, spatial beta, temperature, null/uncertain/band-support
fractions) is confined to it, unobserved voxels are frozen at the explicit
null with zero band support, and the summary reports `observed_fraction`.
The spherical NU support itself is unchanged. The compact evidence state is a
diagnostic and envelope input only; the in-solve `Q_NU` consumer was removed
on 2026-09-06 (`doc/implementation_notes/pcg_priors_history.md`).
With `automsk` enabled the NU-evidence envelope is regenerated from the static
candidate bank while the raw per-voxel evidence margins are live. It remains a
diagnostic under `automsk=yes`; under `automsk=nu` it is the current
coarsest-bank boundary and reference envelope when valid, with density
fallback. Accepted adaptive candidates do not redefine it in the same pass.
`write_nu_evidence_envmask` remains the single producer, called from the
shared `simple_nu_state_filter` on both backends. The full post-hoc NU
filtering path described in this document is production behavior on both
backends.

Auxiliary replacement is conservative. If supplied, the auxiliary pair replaces
the finest discrete label only when its effective resolution is finer than that
label. It is not appended as an extra sidecar candidate.

Persistent unary costs are mask-packed. Full-volume objective arrays are
temporary work buffers; values outside the NU mask must not influence in-mask
objective smoothing or label selection.

Like-for-like selection (2026-09-08). Each candidate's smoothed unary uses
its own radius (1.5 x LP, capped at 30 A), so an argmin over `dmats_mask`
compares differently smoothed fields: two candidates with near-identical raw
unaries do not tie, the smaller radius wins at local minima of the unary
field, the larger at maxima, an intermediate one almost never. An honest
gridding pair never exposes this (adjacent fine candidates differ by the
admitted noise band); a regularized pair does, and the populated fine label
then follows the radius table (PfCRT record 2026-09-08d in
`doc/implementation_notes/pcg_priors_history.md`). The selection is therefore
sequential, coarse to fine: at each level the incumbent and the candidate
are both smoothed at the candidate's radius from the raw unaries kept in
`raw_dmats_mask`, and the candidate wins only with a strictly lower cost.
Identical unaries tie exactly and the coarser label keeps. `dmats_mask`
(own-radius smoothing) remains the input of the Potts prior, the beta
estimate and the evidence envelope, and the Potts sweeps
may only move a voxel to a coarser label than the one it entered with. A
uniform 10% tie margin was tried first and rejected: the natural cost
differences are below 1% between the coarse labels and about 12% at the
fine end, so it collapsed honest pairs to the coarsest label.

## 9. Objective and Label Smoothing

Candidate objective maps are smoothed before voxelwise selection with a
normalized tent kernel over the NU mask. The support radius is candidate-scale:

```text
radius_A = 0.5 * AWF * LP(A)
AWF = 3.0
maximum radius = 30 A
```

The NU filter always applies ordered-label Potts smoothing after the initial
voxelwise selector. This is part of the algorithm, not a workflow switch.

The ordered-label prior:

- uses the 26-neighbor 3D voxel neighborhood
- updates with an 8-color schedule
- evaluates penalties on the ladder label index (integer coordinates of
  the static bank; the log-resolution reference-ladder coordinate of
  2026-09-16 went with the generated ladder); the auxiliary member shares
  the finest retained rung's coordinate, so a boundary between the two is
  free (2026-09-18)
- tolerates jumps of up to one ladder step
- penalizes larger jumps with a linear-quadratic hinge
- normalizes neighbor penalties by the number of in-mask neighbors
- preserves the current label on ties within a small tolerance

Degenerate implementation exits, such as a single label or numerical-zero
beta, may skip smoothing. Users do not select a no-smoothing mode.

## 10. High-Resolution Extension (removed)

The sequential shell walk (`nu_refine=yes`: one Fourier-shell challenge at a
time beyond the finest rung, frontier bookkeeping, majority-z acceptance,
thinned retention, walked-label cleanup, walk depth persisted for restarts)
was retired on 2026-09-16 together with the `nu_refine` control, and the
generated dense ladder that replaced it was withdrawn on 2026-09-18 after
the PfCRT regressions (`latest3`/`latest4`: abinitio3D climb stalls,
refine3D_auto 4.03/4.50 A against 3.93/4.14 A on 2026-09-11 from the same
particles). The static ladder plus the auxiliary pair of section 8 -- the
machinery of commit ed36eb4c's abinitio3D, which produced the best PfCRT
maps -- is the only NU mechanism, in abinitio3D, refine3D_auto and
postprocess_nu alike. Records: `doc/implementation_notes/pcg_decision_log.md`
(2026-09-16 to 2026-09-18).

## 11. Matching References

On the gridding backend, matcher reference loading first looks for NU products.

Plain `nonuniform` prefers independent `_nu_filt` even/odd references. If they
do not exist yet, it falls back to regular even/odd references, then to the
merged state volume if half-map references are unavailable.

`nonuniform_lpset` with active LP-set matching uses the merged registration
reference and prefers the merged `_nu_filt` product when it exists.

State count alone must not force merged-reference matching. The selected NU LP
does not choose reference topology; LP-set mode does.

The ordinary low-pass filter is not applied on top of either NU reference
path.
Reference preparation treats NU filtering, like ML regularization, as filtering
already done during assembly.

Reference preparation applies ordinary spherical soft support first. Assembly
then applies density under `automsk=yes` or valid NU evidence under
`automsk=nu`, with density fallback, before projection.

## 12. FSC Correction and Matching-Bandwidth Handoff

The active FSC is selected before NU filtering for resolution estimation and
reporting:

- `envfsc=no`: use the provisional broad-sphere FSC unchanged
- `envfsc=yes` with gridding: select the density/Otsu envelope, or lagged NU
  evidence with density fallback for `automsk=nu`; find the genuinely
  unmasked FSC 0.8 crossing, independently randomize both half-map phases beyond
  that shell, apply the density envelope, and use
  `(FSC_masked - FSC_randomized_masked) /
  (1 - FSC_randomized_masked)` starting two shells after the crossing
- if no usable crossing exists, retain the genuinely unmasked curve
- PCG: report the FSC on installed solve support without phase randomization

The unmasked, masked, and randomized-masked diagnostics are written as
`fscu_stateNN.bin`, `fsct_stateNN.bin`, and `fscn_stateNN.bin`. The corrected
curve replaces `fsc_stateNN.bin` and its text resolution report.

FSC estimation and NU filtering have separate bandwidth roles. The FSC
enters the bank as its cap (`fsc/1.5`, section 8) and as the auxiliary
member's resolution (its `P_tau` shrinkage and its label resolution). The
NU filter chooses the candidate applied at each volume voxel from its
retained bank. After that volume operation the handoff
(`record_nu_alignment_lowpass_limit`) is the finest selected label whose
cumulative population, that label or finer, reaches
`NU_ALIGN_LP_MIN_SIGNAL_PCT` (1%) of the SIGNAL voxels of the NU mask, i.e.
the mask minus the solvent/background clamp (2026-09-13); the auxiliary
member hands off its own resolution. No gate relative to the whole mask
applies (`min_assigned_pct=0`): the 5% whole-mask support gate of
2026-08-30 capped the PfCRT matching band at 5-6 A against a 4.1 A map
because the coarsest background clamp is a large share of the mask. The
raw finest label (2026-09-02 to 2026-09-13) let 54 voxels of 412k set the
band at 3.37 A against a 3.62 A map. The handoff is logged per state as
`>>> NU MATCHING LOW-PASS HANDOFF: ...` with the raw finest label beside
it. The `fsc/1.5` cap of the bank is what lets the abinitio3D
merged-reference climb (`nonuniform_lpset`) match ahead of its FSC; the
same cap applies to `refine3D_auto` and `postprocess_nu` since 2026-09-18
(the uncapped shell walk of `nu_refine=yes` is gone).

`incrreslim` retains its classical matcher meaning: on an FSC-driven matching
path it permits ten shells beyond the selected FSC criterion. The NU-selected
matching path takes precedence and does not repurpose `incrreslim` as a volume
filter control.

In multi-state runs, the populated state with the finest valid NU-selected
limit determines the single project-level matching bandwidth, matching the
classical global-bandwidth policy.

Staged `abinitio3D` passes an `lpstop` ceiling only in its non-NU stages:
the per-stage `lpstages` limit, or the FSC=0.5 stage-boundary promotion of it
(bounded by the ladder's hard fine bound `LPSTOP_BOUNDS(1)`, 4.5 A). In the
NU stages (`NU_FILTER_STAGE` onwards, `nonuniform_lpset`)
the controller passes no `lpstop` unless the user set one (2026-09-08, the
July policy restored): the handoff is already bounded by the bank (section 8:
the regularized member's FSC=0.143 extent), and a ceiling on top of it only pins the map -- with
the 4.5 A bound in place the PfCRT handoff asked for 4.14/3.98 A, matching
was clamped to 4.5 A and the FSC sat at exactly 4.50 A for 30 iterations. An
NU-selected project limit in those stages may therefore promote matching
beyond the stage plan and beyond 4.5 A. A user `lpstop` remains an explicit
ceiling in every stage (the coarser of it and any stage ceiling applies).
None of this changes the evidence-driven update policy of `refine3D_auto`.

That project `lp` is consumed as follows:

- every nonuniform mode may use it in later non-fresh iterations
- `nonuniform_lpset` also promotes it to command-line `lp`
- fresh stage starts do not consume it unless the run is continuing
- explicit user `lp` remains a hard override
- `lpstop` still caps promoted matching bandwidth

In `nonuniform_lpset`, promotion also activates LP-set topology. Plain
`nonuniform` updates bandwidth while preserving gold-standard half-map
matching.

## 13. Workflow Defaults

`refine3D_auto` defaults to:

- `filt_mode=nonuniform`
- `automsk=yes`
- `ml_reg=yes`
- `envfsc=yes`

With `automsk=yes`, density fixes background voxels of the local filter field.
With `automsk=nu`, valid NU evidence replaces it and also masks references;
density is the fallback.

The `envfsc` default is overridable. When enabled, it uses the independent
density-envelope path described above; it does not enable `automsk` or change NU
support.

`refine3D` exposes `filt_mode`, `automsk`, and `ml_reg` through the ordinary
UI/CLI definitions.

Staged `abinitio3D` defaults to `filt_mode=nonuniform` at the public interface,
but the controller only enables NU filtering from `NU_FILTER_STAGE`; the
bank is the static ladder of section 8. Because abinitio3D is not currently a
gold-standard workflow, staged `nonuniform` is promoted to
`nonuniform_lpset` before the disabled `GOLD_STD_STAGE`. The controller forces
`envfsc=no` before `ENVFSC_STAGE` and forwards the requested value at that stage;
scheduled stage `lp` remains on the refine3D command line.
The default `multivol_mode=independent` policy stops at stage 5, before this
NU-filtering and envfsc stage boundary, unless the user explicitly requests
later stages. The separate final original-sampling reconstruction still
inherits the parent `envfsc` request.

The abinitio3D cavgs route disables NU filtering and automasking.
