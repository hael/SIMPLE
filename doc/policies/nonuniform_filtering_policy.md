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

`filt_mode=nonuniform` activates the nonuniform volume filter path.

In staged `abinitio3D`, `filt_mode=nonuniform` stays on the static
discrete-bank policy with `nu_refine=no`, lets the ML-regularized auxiliary
pair replace the finest discrete NU label only when its effective resolution is
finer than that label, and promotes the finest selected NU bandwidth back into
the matching `lp`. The older abinitio3D automatic low-pass modes,
`filt_mode=uniform` and `filt_mode=fsc`, are not supported. From
`GOLD_STD_STAGE=5`, single-state abinitio3D enables envelope-masked
gold-standard FSC reporting and removes the scheduled stage `lp`; the matching
bandwidth still comes from the NU-selected project `lp`. Multi-state
abinitio3D keeps gold-standard matching and `envfsc` off, while still allowing
automasking and the NU-selected matching `lp`. The high-resolution
`nu_refine=yes` ratchet is reserved for `refine3D_auto`.

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
- optionally, an auxiliary pre-filtered even/odd replacement pair supplied by `volassemble`
- an auxiliary effective resolution, in Angstrom, for the replacement test
- a logical support mask

Mask precedence is:

1. freshly regenerated state-specific automask from `volassemble`
2. existing compatible state-specific automask, when `automsk != 'no'`
3. spherical support mask derived from `mskdiam`

When `ml_reg=yes`, `volassemble` uses the `_unfil` even/odd pair as the base nonuniform input. In the static NU bank path (`nu_refine=no`), it may also supply the ML-regularized even/odd pair as an auxiliary replacement source. The auxiliary pair is ignored unless its effective resolution is finer than the finest discrete bank member. If it is finer, it replaces that finest label; it is not appended as an extra label beside the discrete bank. In the refinement-ratchet path (`nu_refine=yes`, as in `refine3D_auto`), ML-regularized auxiliary replacements are not supplied; the high-resolution shell challenger owns the refinement experiment. Any supplied pairs are clean of low-resolution insertion and receive any trailing blend before filtering. The auxiliary effective resolution is derived from the state FSC(0.143) resolution, `res0143s(state)`.

## Outputs

For each state, the filter produces:

- `vol_state_even_nu_filt.mrc`
- `vol_state_odd_nu_filt.mrc`
- `vol_state_nu_filt.mrc`

Actual filenames are built by appending `NUFILT_SUFFIX`, currently `_nu_filt`, to the base even, odd, and merged state volume names. These are derived products. The base even/odd and merged volumes remain the primary reconstruction outputs.

## Current implementation strategy

The filter currently:

1. builds a bank of low-pass filtered even/odd volumes from the unfiltered pair
2. optionally replaces the finest discrete label with an auxiliary pre-filtered even/odd pair only when the auxiliary resolution is finer
3. maps all retained labels onto a shared filter-bank coordinate axis
4. caches the low-pass-filtered base bank on disk
5. computes voxelwise objective maps across all candidates
6. applies candidate-scale, mask-normalized AWF-like objective smoothing
7. chooses the best candidate per voxel
8. optionally refines the candidate map with an ordered-label spatial smoothing prior
9. logs selected-low-pass statistics and neighbor-continuity diagnostics
10. synthesizes filtered even/odd outputs from the selected candidate map
11. writes the merged `_nu_filt` volume as the average of the filtered even/odd outputs

This design is scientifically reasonable and easy to debug, but it pays a large I/O and memory-traffic cost.

In nonuniform mode, single-state matcher reference loading tries `_nu_filt`
even/odd references first, falls back to the regular even/odd references before
filtered products exist, and avoids applying the ordinary low-pass filter on
top of the nonuniform reference path. Multi-state matcher reference loading
uses the merged state reference instead of independent even/odd half-map
references; when NU products exist, it uses the merged `_nu_filt` state volume.

When `filt_mode=nonuniform` and the user has not set an explicit `lp`, the 3D
matching/reprojection low-pass limit follows the finest selected NU
filter-bank limit from the previous `volassemble` pass instead of the global
FSC resolution. In `nu_refine=yes` shell-extension runs, candidate bins that
fail the 5% tested-frontier acceptance threshold do not advance the global
matching bandwidth. After writing the matching `_nu_filt` even/odd products,
`volassemble` updates the ordinary project `lp` field to that selected NU
limit.
The matcher reads that project `lp` value on the next iteration; a fresh first
iteration or missing project `lp` falls back to the ordinary FSC/project-`lp`
policy. Explicit `lp` remains a hard user override, and `lpstop` still caps the
selected matching bandwidth. The NU-derived matching LP is independent of the
reference topology: single-state gold-standard runs consume even/odd
NU-filtered references independently, while multi-state runs consume merged
state references and keep `envfsc` disabled.

## Candidate-scale objective smoothing

Before voxelwise selection, each candidate objective map is locally averaged
with a normalized tent kernel over the NU mask. The tent radius is tied to the
candidate's effective low-pass limit rather than a fixed physical width:
`radius_A = 0.5 * AWF * LP(A)`, currently with `AWF = 3.0` and an upper cap of
30 A. An auxiliary replacement label uses the supplied effective resolution for
objective smoothing and low-pass reporting while keeping the finest retained
label coordinate for Potts regularization.

This is an AWF-like evidence aggregation step: lower-resolution candidates are
judged from broader local evidence, while high-resolution challengers remain
more spatially local. The masked normalization is radius-specific and cached,
so each objective map only needs the current temporary full-volume work array
plus the persistent mask-packed objective bank.

## Ordered-label smoothing regularization

The nonuniform filter always applies ordered-label smoothing after the unary
voxelwise selector described above and before writing filtered volumes.

The smoothing stage is intended to reduce abrupt local jumps in the selected filter-bank label. It initializes from the ordinary voxelwise argmin, then runs a small number of ICM-style passes over the label map using the fully connected 26-neighbor 3D voxel neighborhood. Updates use an 8-color parity schedule so voxels updated within the same pass are not neighbors under the full 3x3x3 neighborhood. The neighborhood penalty is evaluated on a candidate-coordinate axis rather than on raw label indices:

- retained low-pass candidates use coordinates `1..n_base`
- an auxiliary replacement, when active, uses the finest retained coordinate
- one-step retained-bank coordinate differences are tolerated by the Potts
  hinge, allowing smooth local transitions between adjacent retained filters
- jumps larger than one retained-bank coordinate step receive a
  linear-quadratic penalty to discourage abrupt resolution jumps more strongly
- neighbor penalties are normalized by the number of in-mask neighbors, so boundary and thin-mask voxels do not receive systematically weaker or stronger regularization
- ties preserve the current label within a small tolerance instead of drifting to the lowest candidate index

Auxiliary replacement resolutions are mandatory whenever auxiliary replacement
volumes are supplied. In the current `volassemble` ML-regularized static-bank
path, the auxiliary pair is the ML-regularized even/odd pair and its effective
resolution is derived from the state FSC(0.143) resolution, `res0143s(state)`.
When `nu_refine=yes`, no ML-regularized auxiliary pair is supplied to the
initial NU competition.

Diagnostics log the estimated smoothing beta, retained label counts,
jump-penalty settings, candidate coordinates, changed voxels per iteration, and
mean site energy. Auxiliary replacement diagnostics also log whether the
replacement was accepted or ignored, the replacement effective resolution, and
the unary margin versus the best coarser retained label.
Neighbor-continuity diagnostics report unique 26-neighbor links and classify
retained-bank coordinate differences as identical links, Potts-tolerated
one-step transitions, or penalized larger jumps. One-step transitions remain
visible in the log as context, but the continuity health assessment is driven by
the larger jumps that the ordered-label prior actually penalizes.

### Ordered-label smoothing

Ordered-label smoothing is always active in the NU filter. Standalone
`nu_filt3D`, iterative `volassemble`, and coupled high-resolution refinement
all use the same smoothing policy.

The Cartesian `volassemble` path and standalone `nu_filt3D` therefore call the
optimizer without any smoothing switch:

```fortran
call optimize_nu_cutoff_finds()
```

### High-resolution bank extension

The optional high-resolution extension path can add finer low-pass candidates
after the initial candidate map has been selected. It first identifies voxels
currently assigned to the finest populated base-bank label, prunes any empty
finer base-bank labels, evaluates the next high-resolution candidate only
within that local extension mask, and then updates only those eligible voxels.
This challenge is unary-only: the full-bank Potts prior is not applied during
extension because the extension experiment is already constrained to a
one-shell step on the current finest populated frontier. After the sequential
extension walk stops, the accepted label field is cleaned with the ordinary
ordered-label Potts prior over the final accepted bank. This cleanup uses the
mask-packed unary costs for all active labels, including accepted shell
challengers, and therefore can suppress spatially isolated high-resolution
islands without changing the shell-by-shell acceptance test itself.

Iterative workflows gate this behavior through `nu_refine`. The default is
`nu_refine=no`. Staged `abinitio3D` uses the static discrete-bank policy with
`nu_refine=no` for all stages, including automasked stages. `refine3D_auto`
defaults `nu_refine=yes` and still allows an explicit override.

When `nu_refine=yes`, `volassemble` uses an evidence-driven shell ratchet: an
iteration can evaluate and accept one or more speculative high-resolution
Fourier-shell candidates, but only sequentially. Each challenge candidate is the
next unrepresented shell above the current finest populated base-bank label, not
a coarse hard-coded Angstrom ladder. A challenge is attempted whenever at least
one base-bank voxel sits on the current finest label. In iterative refinement,
the candidate is applied to the current map and promoted into the next
iteration's starting bank only when the unary objective moves at least 5% of the
tested frontier voxels to that speculative shell and the challenger has enough
absolute support to avoid one-voxel shell walking. If a candidate is accepted,
the same iteration may challenge the next Fourier shell; the walk stops at the
first unattempted challenger or the first challenger below this conservative
frontier rule. Each challenge logs the old/new Fourier shell, tested frontier
size, unary wins, and accepted-shell depth for the next iteration. The shell
extension helper still supports `accept_pct=0` for manual diagnostics.

The shell challenge itself stays at the full Fourier sampling rate, but the
retained high-resolution extension bank is thinned for memory. Accepted odd
shell steps are kept as temporary frontiers so the next contiguous shell can be
tested; once the next retained step is accepted, the lower-resolution temporary
label is dropped and its voxels fall back to the nearest retained coarser base
label. Persisted high-resolution depth is rebuilt with the same policy: keep
every second extension shell plus the current terminal shell. The retained
filter-bank coordinates are compacted after thinning, so adjacent retained
labels remain one coordinate step apart even when they are two Fourier shell samples
apart. This bounds the mask-packed unary bank for fine-sampling data sets
without skipping any shell challenge.

The refinement implementation keeps a hard cap on the number of mask-packed
distance-matrix candidates retained at once. When the cap is reached, selected
base labels are preserved and unselected high-resolution labels are compacted
out of the active objective bank before another shell can be accepted. If the
cap remains full after compaction, extension stops rather than allocating an
unbounded unary bank. Once labels are finalized and no further Potts cleanup is
needed, the retained unary bank and associated mask-index work arrays are
released before synthesizing the output volumes.

Large-volume memory ownership is deliberately biased toward mask-packed or
short-lived state. Label/source maps use a compact integer kind, the finest
frontier objective cache is stored only as a mask-packed vector, and the caller
mask can be released after `setup_nu_dmats` copies it into the filter state.
The normalized smoothing support is allocated lazily for the active radius. It
is released after static candidate-bank smoothing, but is kept across adjacent
high-resolution shell challenges so repeated equal-radius extension steps can
reuse it; `nu_filter_vols` releases it before output synthesis. Accepted shell
insertion and stride thinning are combined into one mask-packed bank rebuild so
the common extension path avoids appending the unary bank and immediately
copying it again for thinning.

Static discrete NU filtering writes a single matching low-pass limit back to
the project. This limit is the finest assigned NU label in the state mask; when
an auxiliary ML-regularized replacement is active and selected, its actual
auxiliary resolution participates in that minimum through the replaced finest
label. Staged `abinitio3D` uses this promoted NU limit with
`filt_mode=nonuniform`, so static NU reference generation influences the next
matching bandwidth without enabling the high-resolution shell ratchet.

Map postprocessing is classical-only. `postprocess` and the automatic
`reconstruct3D` postprocess step use the ordinary global FSC/B-factor path even
when earlier reconstruction stages used a NU `filt_mode`. Nonuniform filtering
remains a volume-domain reference-generation feature that writes `_nu_filt`
derived products during assembly; it does not define a separate final-map
postprocess workflow.

The candidate objective bank should be stored only for voxels inside the NU
mask. Full-volume objective arrays are temporary work buffers for objective
generation and candidate-scale, mask-normalized tent smoothing; persistent
unary costs used for candidate selection and ordered-label smoothing are
mask-packed. Objective values outside the NU mask must not contribute to
smoothed in-mask unary costs.

Terminal original-sampling `abinitio3D` reconstruction uses the same classical
postprocessing policy. If staged ab initio refinement used a NU `filt_mode`,
that policy controls the staged `_nu_filt` reference generation, but the
terminal postprocess map remains classical.

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
- When `ml_reg=yes`, `volassemble` should feed the `_unfil` pair as the base nonuniform input. It may supply the ML-regularized pair as an auxiliary replacement only in static-bank NU filtering (`nu_refine=no`), and only use it when it extends beyond the finest discrete label.
