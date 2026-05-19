# Nonuniform Filtering Policy

## Scope

This document defines how nonuniform filtering is selected, where it runs, what mask it uses, and how it should evolve inside the current `refine3D` and `abinitio3D` architecture.

## Public policy

The user-facing control is `filt_mode`.

Supported values:

- `none`
- `uniform`
- `fsc`
- `nonuniform`

Only `filt_mode=nonuniform` activates the nonuniform volume filter path.

## Execution point

Nonuniform filtering is executed in Cartesian `volassemble`, after:

- even/odd assembly
- sampling-density correction and FSC estimation
- merged-volume restoration and write
- gridding correction
- optional trailing half-map blending
- per-state postprocessing planning and any required automask regeneration

This is intentional. The expensive shared-memory volume work happens in one place for both shared-memory and distributed `refine3D` paths, with `volassemble` owning the derived filtered references.

The filter must not consume even/odd maps after low-resolution insertion. Low-resolution even/odd blending is a registration-reference trick applied only after reference read/mask/filter and immediately before generating the reprojection model for 3D registration. `volassemble` must not write low-resolution-blended half-maps to disk, and those blended maps must not feed FSC, automasking, nonuniform filtering, or ordinary half-map handoffs.

If trailing reconstruction is active, the trailing blend is applied to the clean restored half-maps before automask generation and before `_nu_filt` products are generated. Because low-resolution-blended half-maps are never persisted, the ordinary previous even/odd half-maps remain valid trailing inputs.

## Inputs

The nonuniform filter consumes:

- the current unfiltered even volume
- the current unfiltered odd volume
- optionally, auxiliary pre-filtered even/odd candidate pairs supplied by `volassemble`
- an auxiliary effective resolution, in Angstrom, for each auxiliary candidate pair
- a logical support mask

Mask precedence is:

1. freshly regenerated state-specific automask from `volassemble`
2. existing compatible state-specific automask, when `automsk != 'no'`
3. spherical support mask derived from `mskdiam`

When `ml_reg=yes`, `volassemble` uses the `_unfil` even/odd pair as the base nonuniform input. In the static NU bank path (`nu_refine=no`), it may also supply the ML-regularized even/odd pair as an auxiliary candidate source. In the refinement-ratchet path (`nu_refine=yes`, as in `refine3D_auto`), ML-regularized auxiliary candidates are not supplied; the high-resolution shell challenger owns the refinement experiment. Any supplied pairs are clean of low-resolution insertion and receive any trailing blend before filtering. The auxiliary coordinate is derived from the state FSC(0.143) resolution, `res0143s(state)`.

The auxiliary source remains a distinct candidate source. It does not replace a
base low-pass bank member. When an auxiliary candidate wins, source provenance is
tracked separately from the effective local low-pass label: `srcmap` identifies
the auxiliary source used for output synthesis, while `filtmap` stores the
nearest base-bank label implied by the auxiliary effective-resolution coordinate
for diagnostics and map visualization.

## Outputs

For each state, the filter produces:

- `vol_state_even_nu_filt.mrc`
- `vol_state_odd_nu_filt.mrc`
- `vol_state_nu_filt.mrc`

Actual filenames are built by appending `NUFILT_SUFFIX`, currently `_nu_filt`, to the base even, odd, and merged state volume names. These are derived products. The base even/odd and merged volumes remain the primary reconstruction outputs.

## Current implementation strategy

The filter currently:

1. builds a bank of low-pass filtered even/odd volumes from the unfiltered pair
2. optionally appends auxiliary pre-filtered even/odd pairs to that candidate bank only for static-bank NU filtering
3. maps all candidates onto a shared filter-bank coordinate axis
4. caches the low-pass-filtered base bank on disk
5. computes voxelwise objective maps across all candidates
6. chooses the best candidate per voxel
7. optionally refines the candidate map with an ordered-label spatial smoothing prior
8. logs selected-low-pass statistics and neighbor-continuity diagnostics
9. synthesizes filtered even/odd outputs from the selected candidate map
10. writes the merged `_nu_filt` volume as the average of the filtered even/odd outputs

This design is scientifically reasonable and easy to debug, but it pays a large I/O and memory-traffic cost.

In nonuniform mode, matcher reference loading tries `_nu_filt` even/odd references first, falls back to the regular even/odd references before filtered products exist, and avoids applying the ordinary low-pass filter on top of the nonuniform reference path.

When `filt_mode=nonuniform`, `nu_refine=yes`, and the user has not set an
explicit `lp`, the 3D matching/reprojection low-pass limit follows the finest
active NU filter-bank limit from the previous `volassemble` pass instead of the
global FSC resolution. After writing the matching `_nu_filt` even/odd products,
`volassemble` updates the ordinary project `lp` field to the finest active NU
bank limit, using the same project metadata convention as `set_bp_range3D`.
The matcher reads that project `lp` value on the next iteration; a fresh first
iteration or missing project `lp` falls back to the ordinary FSC/project-`lp`
policy. Explicit `lp` remains a hard user override, and `lpstop` still caps the
selected matching bandwidth.

## Ordered-label smoothing regularization

The nonuniform filter always applies ordered-label smoothing after the unary
voxelwise selector described above and before writing filtered volumes.

The smoothing stage is intended to reduce abrupt local jumps in the selected filter-bank label. It initializes from the ordinary voxelwise argmin, then runs a small number of ICM-style passes over the label map using the fully connected 26-neighbor 3D voxel neighborhood. Updates use an 8-color parity schedule so voxels updated within the same pass are not neighbors under the full 3x3x3 neighborhood. The neighborhood penalty is evaluated on a candidate-coordinate axis rather than on raw label indices:

- base low-pass candidates use coordinates `1..n_base`
- auxiliary candidates are projected onto that same axis from their required effective resolutions
- one-step differences are tolerated by the current penalty setting
- jumps larger than one filter-bank step are penalized with a convex hinge to discourage larger jumps more strongly
- neighbor penalties are normalized by the number of in-mask neighbors, so boundary and thin-mask voxels do not receive systematically weaker or stronger regularization
- ties preserve the current label within a small tolerance instead of drifting to the lowest candidate index

Auxiliary candidate resolutions are mandatory whenever auxiliary candidate volumes are supplied. In the current `volassemble` ML-regularized static-bank path, the auxiliary candidate is the ML-regularized even/odd pair and its coordinate is derived from the state FSC(0.143) resolution, `res0143s(state)`. When `nu_refine=yes`, no ML-regularized auxiliary candidate is supplied to the initial NU competition.

Diagnostics log the estimated smoothing beta, candidate and auxiliary counts, jump-penalty settings, candidate coordinates, changed voxels per iteration, and mean site energy.
Auxiliary-candidate diagnostics also log the unary aux-vs-best-base margin and
the auxiliary source assignment counts before and after ordered-label smoothing,
so workflow logs can distinguish an auxiliary candidate that loses the unary
competition from one that is removed by the spatial prior.

### Ordered-label smoothing

Ordered-label smoothing is always active in the NU filter. Standalone
`nu_filt3D`, iterative `volassemble`, coupled high-resolution refinement, and
the `postprocess_nu` Fourier-shell workflow all use the same smoothing policy.

The Cartesian `volassemble` path and standalone `nu_filt3D` therefore call the
optimizer without any smoothing switch:

```fortran
call optimize_nu_cutoff_finds()
```

### High-resolution bank extension

The optional high-resolution extension path can add finer low-pass candidates after the initial candidate map has been selected. It first identifies voxels currently assigned to the finest base-bank label, evaluates the next high-resolution candidate only within that local extension mask, and then updates only those eligible voxels.

When requested, the extension step applies the same ordered-label prior to this constrained two-label decision. The old finest label and proposed new high-resolution label are the only labels that can change inside the extension mask, but neighborhood costs are evaluated against the surrounding frozen label field. This lets the extension avoid creating large discontinuities at the boundary of the finest-resolution region while preserving the local nature of the refinement.

Iterative workflows gate this behavior through `nu_refine`. The default is `nu_refine=no`, so staged `abinitio3D` nonuniform filtering does not expand the high-resolution bank unless a caller deliberately opts in. `refine3D_auto` defaults `nu_refine=yes` while still allowing an explicit override.

When `nu_refine=yes`, `volassemble` uses an evidence-driven conservative
ratchet: an iteration can evaluate and accept one or more speculative
high-resolution Fourier-shell candidates, but only sequentially. Each challenge
candidate is the next unrepresented shell above the current finest base-bank
label, not a coarse hard-coded Angstrom ladder. A challenge is attempted
whenever at least one base-bank voxel sits on the current finest label. The
candidate is applied to the current map and promoted into the next iteration's
starting bank only when the ordered-label extension refinement moves at least
5% of the tested frontier voxels inside the NU refinement mask to that
speculative shell. The denominator is the tested extension frontier for the
current challenger, not the total map volume. If a candidate is accepted, the
same iteration may challenge the next Fourier shell; the walk stops at the
first unattempted or rejected challenger.

The uncoupled postprocessing workflow is owned by `postprocess_nu`, not by
`volassemble` or the iterative refinement path. It estimates the NU filter map
from raw even/odd half maps without ML-regularized auxiliary candidates. It
starts from the static base bank, then challenges one unrepresented Fourier
shell at a time. A postprocess challenger is accepted only when the ordered-label
extension refinement assigns at least one tested frontier voxel to that shell;
only accepted shells become part of the active bank before the next shell is
challenged. The postprocess path may continue this sequential challenger loop
through Nyquist, but it must not pre-generate or pre-evaluate the full sampling
ladder.

The candidate objective bank should be stored only for voxels inside the NU
mask. Full-volume objective arrays are temporary work buffers for objective
generation and mask-normalized tent smoothing; persistent unary costs used for
candidate selection and ordered-label smoothing are mask-packed. Objective
values outside the NU mask must not contribute to smoothed in-mask unary costs.

Terminal original-sampling `abinitio3D` reconstruction follows this
postprocessing policy. If the top-level ab initio run used
`filt_mode=nonuniform` and final postprocessing is enabled, the final
`reconstruct3D` command preserves `filt_mode=nonuniform`, which dispatches to
`postprocess_nu` rather than the classical postprocessing-only path. The final
ab initio output copy stage must preserve the NU products under the
`rec_final_stateNN` names, including `_pproc_nu` and `_lp_nu`.

`postprocess_nu` first writes the standard classically postprocessed comparison
outputs using the ordinary FSC-derived filtering path: `_pproc` and `_lp`, plus
the mirrored `_pproc_mirr` map when mirroring is enabled. It then determines or
accepts the B factor in the same way and uses the NU filter map as a local
postprocessing transfer-function selector for the merged reconstruction. Bins
at or better than the global FSC resolution use the global B factor exactly.
Worse local-resolution bins transition through an asymmetric sigmoid toward a
positive damping plateau, with the default inflection near 8 A. For negative
sharpening B factors, `B_floor = alpha * (-B_global)`. The sigmoid fraction is
normalized to be zero at the global FSC resolution, so the B-factor field is
continuous at the boundary between classically sharpened and damped regions.
All bins are then filtered with the same 4-pixel Hann antialiasing window at
the finest active NU bin, rather than with the global FSC transfer, so locally
promoted high-resolution bins are not cut back to the global FSC limit.
`_lp_nu` is written as the corresponding unsharpened Hann-antialiased map. The NU products
are written separately as `_pproc_nu` and `_lp_nu`, plus `_pproc_nu_mirr` when
mirroring is enabled.

## Strengths of the current design

- Shared-memory volume work is centralized in one execution step.
- The mask used for automasking can also guide nonuniform filtering.
- The implementation is explicit and testable.
- Cached filtered volumes reduce repeated FFT/filter work when the same cutoff bank is reused.
- The ordered-label smoothing stage is orthogonal to candidate-bank construction and preserves the output synthesis path.

## Current limitations

- Disk-backed cache files add avoidable filesystem traffic inside an already expensive volume step.
- The cutoff search still performs repeated full-volume passes.
- Ordered-label smoothing regularization is controlled by the NU filter implementation rather than by workflow-level flags.
- Policy and execution are mixed together in `volassemble`, which makes the commander harder to reason about.
- The nonuniform filter depends on module-level cached state, which limits composability and future concurrency.

## Recommended direction

Medium-term improvements:

- Option to dynamically change the filter bank?
- Return both the selected cutoff map and summary statistics as explicit outputs.
- Promote ordered-label smoothing activation and strength to explicit workflow parameters if validation supports it.

Architectural target:

- `volassemble` should orchestrate nonuniform filtering.
- A dedicated volume postprocessing service should own the filter-bank setup, objective evaluation, and output synthesis.
- When `ml_reg=yes`, `volassemble` should feed the `_unfil` pair as the base nonuniform input. It may add the ML-regularized pair as an auxiliary candidate source only in static-bank NU filtering (`nu_refine=no`).
