# Nonuniform Filtering Policy

This document records the current nonuniform-filtering contract. It focuses on
active code paths and should stay aligned with what the implementation does
today.

## 1. Scope

Nonuniform filtering is a volume-domain reference-generation feature. It
selects a local low-pass limit inside a support mask, writes NU-filtered
derived volumes, and may hand a selected matching bandwidth to later
iterations.

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

`filt_mode=nonuniform` enables NU-filtered volume products and keeps ordinary
gold-standard half-map matching unless another policy explicitly changes the
matching bandwidth.

`filt_mode=nonuniform_lpset` enables the same NU filter and promotes the
selected NU bandwidth into an LP-set matching run. LP-set matching uses merged
registration-reference topology.

`nu_refine=yes` enables iterative high-resolution NU shell extension. This is
on by default in `refine3D_auto`, off by default elsewhere, and explicitly set
to `no` by staged `abinitio3D`.

`automsk` and `mskdiam` control the NU support mask. `ml_reg` can provide an
auxiliary even/odd pair for static NU filtering, but not for `nu_refine=yes`
shell-extension runs.

## 3. Ownership

`simple_vol_pproc_policy.f90` owns per-state postprocessing decisions:
automask regeneration/reuse and the source of the NU support mask.

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
half-maps and merged map have been restored and after FSC/resolution metadata
for the state has been calculated.

For each state, `volassemble` then:

1. plans the automask and NU mask source
2. regenerates a state automask if needed
3. builds the NU support mask
4. configures `simple_nu_filter`
5. optimizes the local filter map
6. optionally runs `nu_refine` high-resolution shell extension
7. writes NU-filtered even, odd, merged, and local-resolution products
8. records the selected NU matching low-pass limit for later handoff

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
- a logical support mask
- optionally, an auxiliary even/odd replacement pair
- optionally, the auxiliary pair's effective resolution in Angstrom

When `ml_reg=yes`, `volassemble` uses the `_unfil` even/odd pair as the base NU
input. In static NU mode (`nu_refine=no`), it may also pass the
ML-regularized even/odd pair as an auxiliary replacement source. The auxiliary
effective resolution comes from the state FSC(0.143) resolution,
`res0143s(state)`.

When `nu_refine=yes`, the ML-regularized auxiliary replacement is not supplied;
the high-resolution shell challenger owns the resolution-extension experiment.

## 6. Mask Selection

For workflow `volassemble`, NU mask precedence is:

1. a freshly regenerated state automask
2. an existing compatible state automask when `automsk != 'no'`
3. a spherical support mask derived from `mskdiam`

A reusable state automask must match the current cropped box and cropped
sampling. Missing masks, incompatible masks, the stage start iteration, and
`AMSK_FREQ` iterations trigger regeneration when automasking is enabled.

Standalone `nu_filt3D` follows the same practical mask choice: generate an
automask when `automsk != 'no'`, otherwise use the spherical `mskdiam` mask.

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

The `_nu_locres` map stores spatial frequency in inverse Angstrom. Voxels
outside the NU support mask, and values above Nyquist, are written as zero.

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

The static bank is `[20, 15, 12, 10, 8, 6, 5, 4]` Angstrom before any
high-resolution extension.

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

The challenge test itself is unary-only. After one or more challengers are
accepted, the final expanded label field is cleaned with the same ordered-label
Potts prior used by the static bank.

Accepted shell steps are challenged at full Fourier sampling. The retained
extension bank is thinned for memory: every second extension shell is kept,
plus the current terminal shell. The active mask-packed objective bank is
capped; if compaction cannot free room, extension stops rather than growing
without bound.

The accepted high-resolution depth for each state is persisted in
`nu_highres_depth_stateNN.txt` so the next iteration can seed the same
extension depth.

## 11. Matching References

When NU mode is active, matcher reference loading first looks for NU products.

Plain `nonuniform` prefers independent `_nu_filt` even/odd references. If they
do not exist yet, it falls back to regular even/odd references, then to the
merged state volume if half-map references are unavailable.

`nonuniform_lpset` with active LP-set matching uses the merged registration
reference and prefers the merged `_nu_filt` product when it exists.

State count alone must not force merged-reference matching. The selected NU LP
does not choose reference topology; LP-set mode does.

The ordinary low-pass filter is not applied on top of the NU reference path.
Reference preparation treats NU filtering, like ML regularization, as filtering
already done during assembly.

## 12. Matching Low-Pass Handoff

After writing NU products, `volassemble` records the finest selected NU
low-pass limit for each state. In multi-state runs, the populated state with
the finest selected NU limit determines the single project-level matching
bandwidth, matching the classical global-bandwidth policy.

That project `lp` is consumed narrowly:

- `nu_refine=yes` may use it in later non-fresh iterations
- `nonuniform_lpset` may promote it to command-line `lp`
- fresh stage starts do not consume it unless the run is continuing
- static plain `nonuniform` with `nu_refine=no` does not use it to override
  ordinary matching bandwidth
- explicit user `lp` remains a hard override
- `lpstop` still caps promoted matching bandwidth

In `nonuniform_lpset`, promotion also activates LP-set topology. In plain
`nonuniform`, a `nu_refine` handoff may update bandwidth while preserving
gold-standard half-map matching.

## 13. Workflow Defaults

`refine3D_auto` defaults to:

- `filt_mode=nonuniform`
- `nu_refine=yes`
- `automsk=yes`
- `ml_reg=yes`

`refine3D` exposes `filt_mode`, `nu_refine`, `automsk`, and `ml_reg` through
the ordinary UI/CLI definitions. It does not force NU shell extension unless
requested.

Staged `abinitio3D` defaults to `filt_mode=nonuniform` at the public interface,
but the controller only enables NU filtering from `NU_FILTER_STAGE`. During NU
stages it emits `nu_refine=no`. Because abinitio3D is not currently a
gold-standard workflow, staged `nonuniform` is promoted to
`nonuniform_lpset` before the disabled `GOLD_STD_STAGE`; `envfsc` remains
`no`, and scheduled stage `lp` remains on the refine3D command line.
The default `multivol_mode=independent` policy stops at stage 5, before this
NU-filtering stage boundary, unless the user explicitly requests later stages.

The abinitio3D cavgs route disables NU filtering and automasking.
