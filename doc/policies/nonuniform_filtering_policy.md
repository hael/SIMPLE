# Nonuniform Filtering Policy

This document records the current nonuniform-filtering contract. It focuses on
active code paths and should stay aligned with what the implementation does
today.

## 1. Scope

Nonuniform filtering is a volume-domain regularization feature. On the
gridding backend it selects a local low-pass limit inside spherical `mskdiam`
support and writes NU-filtered derived references. On the PCG backend the same
Potts-regularized local-resolution field supplies the in-solve `Q_NU`
precision and no derived NU-filtered volume is produced. Its candidate bank is
not truncated by the FSC. The finest evidenced cutoff separately governs the
matching bandwidth handed to later iterations.

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

On `rec_backend=pcg` (2026-09-02) both NU modes require the active `Q_NU`
replay, which is the SOLE NU regularization mechanism on that backend: no
post-hoc NU filter runs, no `_nu_filt`/`_nu_locres` products are written, and
the matcher consumes the primary Q_NU-regularized maps directly (per-half for
`nonuniform`, merged for `nonuniform_lpset`) with no additional filtering.
Only the evidence-derived scalar matching-lp handoff survives, for both
modes.

`nu_refine=yes` enables iterative high-resolution NU shell extension. This is
on by default in `refine3D_auto`, off by default elsewhere, and explicitly set
to `no` by staged `abinitio3D`. The static `nu_refine=no` PCG evidence and
four-band operator are a compatibility boundary: adaptive safeguards must not
alter their candidate geometry, posterior masses, or matching handoff.

`mskdiam` controls the spherical NU support mask. `automsk` separately controls
NU-evidence envelope generation, but is valid only while NU filtering is active.
`ml_reg` can
provide an auxiliary even/odd pair for static NU filtering, but not for
`nu_refine=yes` shell-extension runs.

Refinement automasking is available only in NU modes:
`automsk=yes` requires `filt_mode=nonuniform|nonuniform_lpset`;
`automsk=tight` is rejected because the NU-evidence envelope has no Otsu tight
variant. All other filtering modes require `automsk=no`.

`envfsc=no` is the general default and the broad-sphere FSC remains the
reported curve. `refine3D_auto` overrides that default to `yes` unless the user
sets it explicitly.
With `envfsc=yes`, volume assembly generates a density envelope on the fly from
the merged current half maps, smoothing it at `envmsklp`, and performs
phase-randomized solvent correction. `envmsklp` defaults to 20 A through
`ENVMSKLP_DEFAULT`; it is separate from `amsklp`, which continues to set the NU
evidence smoothing scale. The NU-evidence envelope never reaches the FSC. FSC
correction changes reported resolution metadata; it does not truncate the NU
filter bank or directly set the NU matching bandwidth. `envfsc` is independent
of `automsk` and can operate without NU filtering.

Matching references are never multiplied with an envelope before reprojection
(2026-09-02): a reference must not hard-remove density that is present in the
particle images. With `automsk=yes`, the filter-field background -- the
complement of the NU evidence envelope, derived in the same evidence pass that
builds the filter bank -- takes the coarsest bank candidate, so the
envelope-excluded density (detergent micelle, disordered belt) enters the
matching references heavily low-pass filtered rather than removed, and the
same field derives the Q_NU precisions on the pcg backend. With `automsk=no`,
the entire spherical `mskdiam` support remains unconstrained. The conservative
density envelope is not an NU filter-field input; it remains independently
owned by FSC correction and PCG solve support. This masking policy is
independent of `nu_refine`, which controls only candidate-bank extension. The
matcher applies the spherical soft mask only. The former `envref` parameter
has been removed.

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
6. with `automsk=yes`, derives the NU-evidence envelope from the live static
   unaries and fixes its background to the coarsest candidate; with
   `automsk=no`, leaves the spherical field unconstrained
7. optimizes the local filter map
8. optionally runs `nu_refine` high-resolution shell extension inside the
   fixed background
9. writes NU-filtered even, odd, merged, and local-resolution products
10. records the finest locally selected NU low-pass limit for later handoff

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
- optionally, an auxiliary even/odd replacement pair
- optionally, the auxiliary pair's effective resolution in Angstrom

When `ml_reg=yes`, `volassemble` uses the `_unfil` even/odd pair as the base NU
input. In static NU mode (`nu_refine=no`), it may also pass the
ML-regularized even/odd pair as an auxiliary replacement source. The auxiliary
effective resolution comes from the state FSC(0.143) resolution,
`res0143s(state)`.

When `nu_refine=yes`, the ML-regularized auxiliary replacement is not supplied;
the high-resolution shell challenger owns the resolution-extension experiment.

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
- `nu_envmask3D_stateNN.mrc` is generated from NU evidence when `automsk=yes`;
  it constrains the filter field's background and is never multiplied into a
  matching reference or installed as PCG solve support
- density-derived `automask3D_stateNN.mrc` is generated during volume assembly
  when `envfsc=yes`; the same generated mask feeds phase-randomized FSC and
  cFAR and may be reused by non-PCG final postprocessing
- PCG builds a conservative density support independently of `automsk`:
  `envfsc=yes` installs it in the base and replay solves, while `envfsc=no`
  leaves the base solve spherical and installs it in the replay only. Before
  any reconstruction exists, the base necessarily bootstraps on the sphere
  and supplies the replay's density support for either `envfsc` value
- the correlation-derived NU-evidence envelope never feeds FSC correction or
  NU support

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
6. select the best candidate per in-mask voxel
7. apply ordered-label Potts smoothing to the candidate map
8. synthesize filtered even/odd outputs from the selected labels
9. write the merged `_nu_filt` output as the even/odd average
10. write the same-grid `_nu_locres` map

The nominal static bank is `[20, 15, 12, 10, 8, 6, 5, 4]` Angstrom before any
high-resolution extension. The FSC does not remove candidates from this bank.

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
The spherical NU support itself is unchanged. The PCG replay consumes the expanded weights as
the `Q_NU` precision when `pcg_nu_lambda_rel > 0`
(`doc/implementation_notes/pcg_priors.md` Stage 6). Since 2026-08-28 this is
the DEFAULT whenever `rec_backend=pcg` runs with NU filtering and the euclid
ML replay; in that mode the post-hoc NU filter and its
`_nu_filt`/`_nu_locres` products are NOT generated -- the in-solve precision
already performs the local regularization, and the LP-set matching handoff
derives from the frozen replay evidence (finest evidenced local cutoff).
With `automsk` enabled the NU-evidence envelope is still produced. It is
regenerated from the static candidate bank while the raw per-voxel evidence
margins are live, then fixed as a coarsest-bank boundary condition during
ordered-label Potts regularization and adaptive shell challenges. Accepted
adaptive candidates refine the precision inside that preliminary boundary;
they do not redefine their own support in the same pass.
That clamped field supplies `Q_NU`; the artifact is not a second filtering or
reference-masking mechanism. `write_nu_evidence_envmask` remains the single
producer shared by gridding volume assembly and PCG evidence construction.
The full post-hoc NU filtering path described in this document is production
behavior only for the gridding backend. PCG+NU rejects a disabled `Q_NU`
rather than falling back to post-hoc filtering.

Auxiliary replacement is conservative. If supplied, the auxiliary pair replaces
the finest discrete label only when its effective resolution is finer than that
label. It is not appended as an extra sidecar candidate.

Persistent unary costs are mask-packed. Full-volume objective arrays are
temporary work buffers; values outside the NU mask must not influence in-mask
objective smoothing or label selection.

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
- evaluates penalties on retained-bank coordinates, not raw label numbers
- tolerates adjacent retained-bank coordinate steps
- penalizes larger jumps with a linear-quadratic hinge
- normalizes neighbor penalties by the number of in-mask neighbors
- preserves the current label on ties within a small tolerance

Degenerate implementation exits, such as a single label or numerical-zero
beta, may skip smoothing. Users do not select a no-smoothing mode.

## 10. High-Resolution Extension

`nu_refine=yes` enables a sequential high-resolution shell ratchet after the
static bank has been optimized.

The extension:

- starts from voxels assigned to the finest populated retained label
- challenges the next unrepresented Fourier shell
- evaluates the challenger only on that frontier mask
- accepts a challenger only when enough tested frontier voxels prefer it
- requires at least 5% challenger wins and a minimum absolute seed support
- may accept multiple contiguous shell steps in one iteration
- stops at the first unattempted, unsupported, or rejected challenger
- never proposes a shell beyond the Fourier grid

For PCG only, adaptive discovery also stops two Fourier shells beyond the
current evidence-pair FSC=0.143 crossing. This bounds repeated use of the same
half pair while leaving the established gridding shell walk unchanged.

The challenge test itself is unary-only. After one or more challengers are
accepted, the final expanded label field is cleaned with the same ordered-label
Potts prior used by the static bank.

Accepted shell steps are challenged at full Fourier sampling. The retained
extension bank is thinned for memory: every second extension shell is kept,
plus the current terminal shell. The active mask-packed objective bank is
capped; if compaction cannot free room, extension stops rather than growing
without bound.

Each accepted challenger also updates the raw best-candidate cost before
candidate-scale objective smoothing. A diagnostic envelope generated after the
walk can therefore describe the accepted bank. With `automsk=yes`, however,
the filtering/precision boundary is deliberately generated from the static
bank before the walk, because the boundary must constrain rather than be
selected by the adaptive candidates.

On the gridding path the accepted high-resolution depth for each state is persisted in
`nu_highres_depth_stateNN.txt` so the next iteration can seed the same
extension depth.

For PCG, accepted candidates are converted to an evidence posterior using
Fourier-shell coordinates anchored to the finest static-bank interval. The
static coordinates remain exactly 1 through 8; only extension candidates use
the shell-distance continuation. Voronoi cell widths supply the candidate
integration measure, normalized to the static bank's total signal mass.
Adding, removing, or thinning densely spaced shell probes therefore does not
create fine-band support merely by changing the label count. The static
eight-candidate bank retains integer coordinates and unit masses exactly.

The PCG adaptive matching handoff requires the established 5% assignment
support and then guarantees two Fourier shells beyond the FSC crossing. This
keeps the matching/evidence loop able to advance without letting a single
high-resolution voxel set the global bandwidth. `nu_refine=no` retains the
historical raw-finest handoff.

## 11. Matching References

On the gridding backend, matcher reference loading first looks for NU products.

Plain `nonuniform` prefers independent `_nu_filt` even/odd references. If they
do not exist yet, it falls back to regular even/odd references, then to the
merged state volume if half-map references are unavailable.

`nonuniform_lpset` with active LP-set matching uses the merged registration
reference and prefers the merged `_nu_filt` product when it exists.

State count alone must not force merged-reference matching. The selected NU LP
does not choose reference topology; LP-set mode does.

On the PCG backend, both NU modes use the primary `Q_NU`-regularized maps:
independent halves for `nonuniform`, and the merged map for
`nonuniform_lpset`. PCG never looks for `_nu_filt` products and applies no
additional NU filtering.

The ordinary low-pass filter is not applied on top of either NU reference
path.
Reference preparation treats NU filtering, like ML regularization, as filtering
already done during assembly.

Reference preparation never multiplies either the NU-evidence or density
envelope into a reference. The matcher applies only the ordinary spherical
soft support before projection. Envelope information reaches a PCG matching
reference solely through the Potts-regularized `Q_NU` solve.

## 12. FSC Correction and Matching-Bandwidth Handoff

The active FSC is selected before NU filtering for resolution estimation and
reporting:

- `envfsc=no`: use the provisional broad-sphere FSC unchanged
- `envfsc=yes`: generate a non-tight density/Otsu envelope on the fly from the
  merged half maps after low-pass filtering at `envmsklp`; find the genuinely
  unmasked FSC 0.8 crossing, independently randomize both half-map phases beyond
  that shell, apply the density envelope, and use
  `(FSC_masked - FSC_randomized_masked) /
  (1 - FSC_randomized_masked)` starting two shells after the crossing
- if no usable crossing exists, retain the genuinely unmasked curve

The unmasked, masked, and randomized-masked diagnostics are written as
`fscu_stateNN.bin`, `fsct_stateNN.bin`, and `fscn_stateNN.bin`. The corrected
curve replaces `fsc_stateNN.bin` and its text resolution report.

FSC estimation and NU filtering have separate bandwidth roles. The FSC never
caps the local filter bank or its high-resolution extension. The NU filter
chooses the cutoff applied at each volume voxel from its full bank. After that
volume operation, the finest cutoff selected anywhere inside the NU support
mask becomes the project-level matching `lp` for the next iteration.

`incrreslim` retains its classical matcher meaning: on an FSC-driven matching
path it permits ten shells beyond the selected FSC criterion. The NU-selected
matching path takes precedence and does not repurpose `incrreslim` as a volume
filter control.

In multi-state runs, the populated state with the finest valid NU-selected
limit determines the single project-level matching bandwidth, matching the
classical global-bandwidth policy.

Staged `abinitio3D` additionally passes its current planned `lpstages` value as
`lpstop`. Consequently, an NU-selected project limit may keep matching coarser
than the plan but cannot promote matching beyond the current ab-initio stage
boundary. This workflow constraint does not change the evidence-driven update
policy used by `refine3D_auto`. When the user explicitly supplies a coarser
`lpstop`, the staged workflow uses the coarser of that value and the planned
stage boundary.

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
- `nu_refine=yes`
- `automsk=yes`
- `ml_reg=yes`
- `envfsc=yes`

With `automsk=yes`, the NU-evidence envelope fixes background voxels of the
local filter field to the coarsest bank member. It is not a reference mask.

The `envfsc` default is overridable. When enabled, it uses the independent
density-envelope path described above; it does not enable `automsk` or change NU
support.

`refine3D` exposes `filt_mode`, `nu_refine`, `automsk`, and `ml_reg` through
the ordinary UI/CLI definitions. It does not force NU shell extension unless
requested.

Staged `abinitio3D` defaults to `filt_mode=nonuniform` at the public interface,
but the controller only enables NU filtering from `NU_FILTER_STAGE`. During NU
stages it emits `nu_refine=no`. Because abinitio3D is not currently a
gold-standard workflow, staged `nonuniform` is promoted to
`nonuniform_lpset` before the disabled `GOLD_STD_STAGE`. The controller forces
`envfsc=no` before `ENVFSC_STAGE` and forwards the requested value at that stage;
scheduled stage `lp` remains on the refine3D command line.
The default `multivol_mode=independent` policy stops at stage 5, before this
NU-filtering and envfsc stage boundary, unless the user explicitly requests
later stages. The separate final original-sampling reconstruction still
inherits the parent `envfsc` request.

The abinitio3D cavgs route disables NU filtering and automasking.
