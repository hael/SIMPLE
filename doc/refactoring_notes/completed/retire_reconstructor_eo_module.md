# Retire the `simple_reconstructor_eo` Module

## Status

Proposed on 2026-08-28; code-reviewed against `master` the same day (call-site
inventory, both PCG FSC/report duplicates, and the gridding module verified;
review findings folded into the sections below). This is the single living
design, review, implementation, and validation record for retiring
`src/main/volume/simple_reconstructor_eo.f90` after removal of the
`reconstructor_eo` composite type.

Implemented 2026-08-28 (all five stages; see the Completion Record). User-side
compilation and runtime validation are pending.

## Objective

Remove the obsolete `simple_reconstructor_eo` module without changing
reconstruction results, workflow behavior, artifact names, or public CLI
behavior.

The previous gridding-backend refactor removed the `reconstructor_eo` type,
its composite lifetime, and `builder%eorecvol`. The file that formerly owned
that type now contains unrelated residual procedures:

- workflow-specific even/odd accumulator I/O;
- gridding half-pair restoration and ML replay sequencing;
- half-map masking, FSC, cFAR, resolution, and report generation.

The first two groups have one workflow owner and do not require a public
library API. The diagnostic group is duplicated in the PCG shared and
distributed paths and should become a backend-neutral volume-domain service.

This is an ownership refactor, not a new feature. It requires no opt-in key and
must leave default and optional reconstruction behavior unchanged.

## Current Call-Site Inventory

As of `master` on 2026-08-28, the only import of
`simple_reconstructor_eo` is inside `restore_state_from_parts` in
`src/main/commanders/simple/simple_commanders_rec_distr.f90`.

| Procedure | Static call sites | Runtime behavior |
| --- | ---: | --- |
| `read_gridding_pair_accumulators` | 2 | Called once for every partial inside the `nparts` reduction loop and, when applicable, once for the persistent trailing chain. |
| `write_gridding_pair_accumulators` | 1 | Called by the contained `write_trail_chain_set`; its three callers represent mutually exclusive seed, blend, and replacement paths. |
| `restore_gridding_pair` | 3 | The calls are mutually exclusive bootstrap/conical, bootstrap/non-conical, and ordinary paths; exactly one executes for each restored state. |

No other source file constructs a `reconstructor_eo`, refers to
`builder%eorecvol`, or calls these procedures. Their static visibility is
therefore broader than their actual ownership.

## Required Behavioral Contract

### Single-accumulator mechanics

`simple_reconstructor` remains the owner of one gridding accumulator and its
numerical mechanics:

- numerator and sampling-density storage;
- raw numerator/rho reads and writes;
- zero padding of a compatible earlier accumulator;
- accumulation and weighting;
- sampling-density correction and ML-prior attachment;
- base and final restoration;
- deapodization and lifecycle management.

The refactor must not move even/odd workflow policy into the `reconstructor`
type or make a single accumulator aware of refinement artifact names.

### Gridding-pair accumulator I/O

The contained pair I/O must preserve all existing behavior:

- filenames remain `<fbody>_even.mrc`, `rho_<fbody>_even.mrc`,
  `<fbody>_odd.mrc`, and `rho_<fbody>_odd.mrc`;
- the four numerator/rho components retain the existing bounded parallel read;
- matching grids are read directly;
- a smaller previous grid is zero-padded only under the existing fractional-
  update contract;
- a larger previous grid is a hard error;
- an incomplete optional pair resets both destination accumulators;
- an incomplete required trailing pair is a hard error;
- trailing-chain manifests, generation checks, and weight policy remain owned
  by `restore_state_from_parts`.

### Gridding half-map restoration

The gridding restoration sequence must remain numerically and operationally
unchanged:

1. Restore the undeapodized base pair used by the legacy gridding FSC oracle.
2. Write deapodized `_unfil` maps when ML regularization is active.
3. Calculate or accept the bootstrap FSC/cFAR and optional conical result.
4. Attach the selected prior to both accumulator densities.
5. Restore and write the final deapodized half maps.
6. Preserve the non-ML path and its merged-accumulator ordering.
7. Preserve the FSC array, FSC=0.5 and FSC=0.143 resolution clamps, cFAR, and
   all current output filenames.

In particular, gridding FSC must continue to use the legacy undeapodized
half-map representation. Consolidating diagnostics must not silently switch it
to the final deapodized representation used by downstream consumers.

### Common half-map diagnostics

The shared diagnostic operation must preserve:

- ordinary spherical-mask FSC;
- optional density-envelope masking and phase-randomized FSC correction;
- cFAR calculation with the existing cone construction;
- FSC=0.5 and FSC=0.143 resolution determination and Nyquist clamping;
- non-destructive treatment of caller-owned half maps;
- existing gridding and PCG mask-radius choices;
- current mask, FSC-array, text-report, and log content.

The common evaluator must operate on explicitly prepared real-space half maps.
It must not infer whether an input came from gridding or PCG, and it must not
perform an implicit backend-dependent FFT, inverse FFT, or deapodization.

The gridding caller will retain a small adapter that reproduces its legacy
representation conversion before calling the common evaluator. PCG will pass
its restored real-space base half maps directly.

Two clarifications from the 2026-08-28 code review:

- `reconstructor%restore_base` already returns a real-space, clipped,
  undeapodized volume, and `phase_rand_fsc` performs its own guarded
  FFT/IFFT round trips on internal copies. The gridding adapter is therefore
  a representation *guarantee* (real-space, undeapodized, legacy clip), not
  a new conversion; today's `ifft` calls on copies inside
  `calc_masked_cfar`/`calculate_gridding_pair_fsc` are state-guarded no-ops
  on this input.
- "The evaluator does not write files" holds at the evaluator level only.
  `phase_rand_fsc` in `simple_fsc` writes the `fscu_state*`, `fsct_state*`,
  and `fscn_state*` binary arrays internally, identically for both backends
  today. These lower-level artifacts are preserved as-is; the no-write
  contract governs masks, FSC arrays, half maps, and text reports selected
  by the workflow.

## Target Ownership

### `simple_reconstructor`

Keep the existing single-accumulator type and `gridding_half_restore` result
mechanics in `src/main/volume/simple_reconstructor.f90`. Do not add pair
lifetime, refine3D naming, mask, FSC-reporting, or trailing-chain policy to
this module.

### `restore_state_from_parts`

Move the following procedures into `restore_state_from_parts` as sibling
contained procedures:

- `read_gridding_pair_accumulators`;
- `write_gridding_pair_accumulators`;
- `restore_gridding_pair`.

These procedures coordinate type-bound domain operations and workflow
artifacts. They do not provide reusable numerical kernels, and all their call
sites already share the host procedure's state. Their placement makes the
single workflow owner explicit without creating another one-caller module.

The contained procedures should continue to receive explicit arguments where
that makes mutation and lifecycle visible. Host association must not be used
to hide which reconstructor or image objects are modified.

### `simple_halfmap_diagnostics`

Add `src/main/volume/simple_halfmap_diagnostics.f90` with a narrow public API:

```fortran
type :: halfmap_diagnostics_result
    real, allocatable :: fsc(:)
    real              :: res_fsc05   = 0.
    real              :: res_fsc0143 = 0.
    real              :: cfar        = 0.
  contains
    procedure :: kill
end type halfmap_diagnostics_result

subroutine evaluate_halfmap_pair(params, state, even, odd, average, &
    spherical_mask_radius, result, envmask)

subroutine write_halfmap_diagnostics(result, box, smpd, fname)
```

The exact spelling may follow local conventions during implementation, but the
contract is fixed:

- `even`, `odd`, and `average` are prepared real-space images;
- `spherical_mask_radius` is explicit so gridding and PCG retain their current
  choices;
- `result` owns the FSC and scalar diagnostics;
- `envmask` is an optional returned artifact;
- the evaluator does not write files or select refine3D filenames;
- the writer receives an explicit filename and contains no backend-specific
  policy.

The module may use `simple_fsc` primitives, `image`, `image_msk`, and
`parameters`. It must not depend on `simple_refine3D_fnames`, a commander, or a
reconstruction strategy. `simple_fsc` remains the lower-level general FSC
utility and should not acquire volume-workflow dependencies merely to avoid a
new owning module.

### Callers

The gridding commander and PCG strategy remain responsible for:

- preparing backend-correct diagnostic half maps;
- selecting the existing spherical mask radius;
- resolving refine3D artifact filenames;
- writing the returned envelope mask and FSC array;
- execution-specific logging;
- deciding when diagnostics are required.

## Implementation Plan

### Stage 1: Freeze and expose the common contract

1. Add `simple_halfmap_diagnostics` and
   `halfmap_diagnostics_result`.
2. Move the backend-neutral mask, phase-randomized FSC, cFAR, resolution, and
   explicit-filename report code into it.
3. Keep evaluator inputs non-destructive by operating on owned work images.
4. Keep envelope-mask generation separate from artifact writing by returning
   the optional mask to the caller.

### Stage 2: Remove PCG duplication

1. Replace `calculate_state_fsc` and `calculate_distributed_fsc` in
   `simple_rec3D_pcg_strategy.f90` with calls to the common evaluator.
2. Replace `write_fsc_summary` and `write_distributed_fsc_summary` with the
   common explicit-filename writer.
3. Retain shared/distributed log prefixes and orchestration at the strategy
   call sites.
4. Confirm both paths pass identical scientific policy and differ only in
   execution mechanics.

Migrating PCG first establishes the common API using maps already represented
in its required real-space input domain.

### Stage 3: Internalize gridding-pair ownership

1. Move pair read, pair write, and pair restoration into
   `restore_state_from_parts`.
2. Preserve the existing argument lists initially; simplify them only after
   behavior parity is established.
3. Keep the bounded parallel pair read and all missing/required artifact
   semantics unchanged.
4. Keep trailing-chain manifest and weighting decisions in their existing
   contained workflow procedures.

### Stage 4: Migrate gridding diagnostics

1. Add a contained gridding adapter that produces the exact real-space
   representation corresponding to the legacy undeapodized FSC calculation.
   The adapter applies to the ordinary restored base pair only. The trailing
   bootstrap FSC oracle already operates on the previous iteration's final
   deapodized half maps read via `read_and_crop`; that call site passes its
   maps to the common evaluator directly, with no adapter, preserving the
   existing representation asymmetry between the two oracles.
2. Replace `calculate_gridding_pair_fsc` with the common evaluator through
   that adapter.
3. Replace `gridding_pair_diagnostics` and its writer with the common result
   and writer.
4. Keep FSC-array, automask, and resolution-report filenames explicit in the
   commander.
5. Preserve bootstrap FSC/cFAR/conical injection and ML replay ordering.

### Stage 5: Delete the obsolete module

1. Remove the final `use simple_reconstructor_eo` import.
2. Delete `src/main/volume/simple_reconstructor_eo.f90`. Also update the two
   comment-only references that would otherwise trip the validation grep:
   `src/main/flex/simple_flex_pca_rec3D.f90` (radius-convention comment) and
   `src/main/commanders/test/simple_commanders_test_highlevel.f90`
   (scope comment).
3. Update current architecture and algorithm documentation to point to
   `simple_reconstructor`, `simple_halfmap_diagnostics`, and
   `restore_state_from_parts`.
4. Update generated Fortran indexes after source changes.
5. Propose updates to `.github/skills` references that name the deleted module;
   apply those instruction/skill edits only with explicit approval.

## Review Findings and Risks

### Contained-procedure size

`read_gridding_pair_accumulators` and `restore_gridding_pair` are not trivial,
but another public one-caller module would preserve the same false abstraction.
Their bodies principally coordinate file artifacts and calls to existing
domain methods. Keeping numerical kernels type-bound on `reconstructor` and
moving reusable diagnostics into their own volume module preserves the
commander/domain boundary.

If either contained procedure gains a second independent workflow consumer,
that is the trigger to introduce a narrowly named pair-I/O or pair-restoration
module. Anticipated reuse alone is not sufficient.

### Representation drift

The highest scientific risk is accidentally evaluating gridding FSC on
deapodized final maps instead of its current undeapodized base representation.
The common evaluator's real-space input contract and the explicit gridding
adapter are required to make this transition reviewable.

### Mask-policy drift

The current gridding and PCG paths do not express their spherical mask radius
identically. The common module must accept the radius as an argument rather
than selecting one global value. Density-envelope generation must continue to
use the same `parameters` values and average-map construction as each current
path.

### Known residual backend divergences (2026-08-28 code review)

The refactor is behavior-preserving, so the following backend deltas survive
it. They are recorded here so the shared module makes them visible rather
than hiding them, and so a later coherency decision has a checklist:

- **Spherical mask radius (`envfsc=no`)** — the largest scientific delta.
  Gridding masks with the broad rim radius
  `box_crop/2 - COSMSKHALFWIDTH - 1` for both FSC and cFAR
  (`calculate_gridding_pair_fsc`); PCG masks with the user radius
  `params%msk_crop` for both (`calculate_state_fsc`,
  `calculate_distributed_fsc`). FSC curves and cFAR scores are therefore not
  directly comparable across backends in this mode. Unifying the radius is a
  deliberate, behavior-changing follow-up decision, out of scope here; the
  explicit `spherical_mask_radius` argument is what makes that follow-up a
  one-line policy change per caller.
- **Resolution-report filename policy.** The gridding commander resolves the
  text-report name through `resolve_fsc_txt_fname` (honoring `outfile` and
  `which_iter`); PCG always writes the plain
  `refine3D_resolution_txt_fbody(state)`. The explicit-filename writer
  preserves both; whether PCG should honor `which_iter` is a separate
  question.
- **Envelope-mask source map.** Gridding builds the density envelope from the
  undeapodized base average; PCG builds it from its base-solve average. This
  is inherent to the backends' representations and intentionally preserved.
- **Report-loop bound.** The gridding writer iterates `size(res)` and indexes
  `fsc(k)` directly; the PCG writers bound the loop with
  `min(size(res), size(fsc))`. The arrays are equal-length in practice
  (`get_resarr` yields `filtsz` entries); the common writer adopts the
  defensive `min` form.
- **Duplicated cone constants.** All three cFAR sites construct
  `fsc_area_score_result` with the literal tuple `(256, 20., 0.143, 1)`.
  The shared module owns these as named constants.

### Additional near-duplicate not migrated

`shipped_pair_res` in `simple_rec3D_pcg_strategy.f90` is a third mini
evaluator (soft spherical `msk_crop` mask, FSC, `get_resolution`, Nyquist
clamps) used only as the NU-replay over-regularization diagnostic. It is
deliberately out of scope: it computes no cFAR, writes no artifacts, and its
"never a resolution claim" semantics differ from the half-map diagnostic. If
the common evaluator later grows a cFAR-free spherical mode, it becomes a
natural consumer; do not force it into the initial API.

### Side-effect drift

The current routines mix calculation with mask, FSC, half-map, and report
writes. Moving calculation without auditing every write could remove an
artifact or write it at a different point in the lifecycle. Artifact ownership
and filename selection must be explicit at each migrated call site.

### Lifecycle and parallel-I/O regressions

The refactor must preserve `new`/`kill` symmetry for work images, masks,
conical FSC results, and allocatable FSC arrays. Internalizing pair reads must
not serialize the existing four-component bounded I/O or leave open image and
rho handles on an error-free path.

## Validation Matrix

The refactor is complete only after the following paths have been checked for
numerical and artifact parity:

| Dimension | Required cases |
| --- | --- |
| Backend | gridding, PCG |
| PCG execution | shared, distributed |
| Gridding regularization | ML on, ML off |
| FSC mask | `envfsc=no`, supported `envfsc=yes` path |
| Directional regularization | conical FSC off, conical FSC on |
| Trailing reconstruction | disabled, bootstrap without a chain, established chain |
| Partial accumulator compatibility | matching grid, smaller fractional-update grid, missing optional pair, incomplete required pair, larger incompatible grid |

For every applicable case, compare:

- FSC array values;
- FSC=0.5 and FSC=0.143 resolutions;
- cFAR;
- density-envelope mask content when generated;
- `_unfil`, final half-map, and merged-map content;
- accumulator and rho artifacts;
- resolution report text and artifact names;
- trailing-chain generation and manifest behavior.

## Validation Criteria

Agent-side lightweight checks:

- `rg` finds no remaining source reference to `simple_reconstructor_eo`,
  `reconstructor_eo`, or `builder%eorecvol`;
- the two PCG execution paths call the same half-map evaluator and writer;
- the gridding path has one explicit representation adapter;
- all current artifact writes have an identified caller after the move;
- the Fortran source-index parser completes against `src`;
- generated code-overview indexes are updated;
- `git diff --check` passes;
- editor or language-server diagnostics report no new syntax errors when the
  tooling is available.

User-side compilation and runtime checks:

- build `simple_exec`, `simple_test_phase_rand_fsc`,
  `simple_test_rec3D_backend`, and the continuous PCG half-set tests;
- run `simple_test_phase_rand_fsc`;
- run `simple_test_rec3D_backend` for gridding/PCG comparison;
- run the continuous PCG half-set validation paths;
- exercise the gridding cases in the validation matrix, including trailing and
  fractional-update accumulator restoration;
- compare pre-refactor and post-refactor artifacts within their established
  numerical tolerances;
- run representative Linux/BOX validation before claiming platform parity.

Compilation and runtime tests remain with the user under repository policy.
No Linux/BOX result may be claimed without observed output.

## Completion Record

All five stages implemented on 2026-08-28. Compilation and runtime validation
remain with the user; no build or test binary was run.

### Stage 1 (done)

`src/main/volume/simple_halfmap_diagnostics.f90` added with
`halfmap_diagnostics_result`, `evaluate_halfmap_pair`, and
`write_halfmap_diagnostics`. Deviations from the sketched API, all
contract-preserving:

- `state` is a parameter (needed by `phase_rand_fsc` artifact naming and the
  cone score log line);
- an optional `cones` argument returns the `fsc_area_score_result` needed by
  gridding's conical regularization (`add_conical_invtausq2rho`);
- the cone construction constants `(256, 20., 0.143, 1)`, previously
  duplicated at three sites, are named module parameters;
- the writer uses the defensive `min(size(res), size(fsc))` loop bound;
- the envelope mask is copied to the optional `envmask` before
  `phase_rand_fsc` runs, matching the legacy write-before-correction content
  (the mask is not modified afterwards, so content is identical either way).

Imports are exactly the sanctioned set: core API, `parameters`, `image`,
`image_msk`, `simple_fsc`. No refine3D-filename, commander, or strategy
dependency.

### Stage 2 (done)

Both PCG duplicate pairs deleted. A single module-level
`calculate_pcg_state_diagnostics(params, state, context, even, odd, avg,
diagnostics)` owns the PCG mask policy (`params%msk_crop`), the automask
artifact write, and the context-labelled resolution log lines
(`'RECONSTRUCT3D'` / `'DISTRIBUTED'`); shared and distributed call sites now
provably pass identical scientific policy. Both call sites write the
resolution report through the common explicit-filename writer with
`refine3D_resolution_txt_fbody(state)` resolved at the call site. The
automask artifact is now written after the FSC evaluation instead of before;
its content is unchanged (the mask is never modified by the evaluation). The
`simple_fsc` and `simple_image_msk` imports of the strategy module were
retired with the duplicates.

### Stage 3 (done)

`read_gridding_pair_accumulators`, `write_gridding_pair_accumulators`, and
`restore_gridding_pair` moved verbatim (argument lists preserved, explicit
arguments retained) into `restore_state_from_parts` as sibling contained
procedures. The bounded four-component parallel read, zero-padding contract,
missing/required semantics, and trailing-chain ownership are unchanged.

### Stage 4 (done)

The contained `calc_gridding_pair_diagnostics` adapter builds the merged
average, selects the legacy broad radius
(`box_crop/2 - COSMSKHALFWIDTH - 1`), writes the automask on the envfsc path,
and calls the common evaluator. Both restoration oracles route through it:
the ordinary path with the real-space undeapodized base pair, the trailing
bootstrap with the previous final half maps (no adapter semantics needed, as
reviewed). The redundant pre-`new` of the bootstrap conical result at the
call site was dropped (the evaluator constructs the cone result itself,
exactly as `calc_masked_cfar` already did). `restore_gridding_pair` retains
its trailing FSC/resolution block, `refine3D_fsc_fname` write, and Nyquist
clamps; the text report is written by the commander through the common
writer via `resolve_fsc_txt_fname`.

### Stage 5 (done)

`src/main/volume/simple_reconstructor_eo.f90` deleted; the final import
removed; the two comment-only references
(`simple_flex_pca_rec3D.f90`, `simple_commanders_test_highlevel.f90`)
updated; `doc/algorithms/reconstruction.md` now points to
`simple_reconstructor`, `simple_halfmap_diagnostics`, and
`restore_state_from_parts`. Generated indexes regenerated
(`generate_codeoverview.pl`, `gen_fortran_indexes.pl --root src`); both
parsers completed and the new module is indexed. Skill references in
`.github/skills/` that name the deleted module are listed for approval below
and were NOT edited.

### Lightweight checks performed (agent side)

- `rg` finds no source reference to `simple_reconstructor_eo`,
  `reconstructor_eo`, `eorecvol`, or any retired procedure name under `src`,
  `production`, `scripts`;
- both PCG execution paths call `calculate_pcg_state_diagnostics` →
  `evaluate_halfmap_pair` and the common writer;
- the gridding path has exactly one representation adapter;
- every artifact write from the old routines has an identified caller
  (automask: adapter and PCG helper; FSC array: `restore_gridding_pair`
  gridding / call-site `arr2file` PCG; text report: commander via
  `resolve_fsc_txt_fname` / PCG call sites; fscu/fsct/fscn: unchanged inside
  `phase_rand_fsc`);
- the Fortran source-index parser completed against `src`;
- `git diff --check` passes;
- procedure begin/end balance checked per edited file.

### Pending user-side validation

Build `simple_exec`, `simple_test_phase_rand_fsc`, `simple_test_rec3D_backend`
and the continuous PCG half-set tests; run the validation matrix above and
compare artifacts within established tolerances. No Linux/BOX result claimed.

### Skill edits (approved and applied 2026-08-28)

- `.github/skills/simple-main-volume/SKILL.md`: replace
  `simple_reconstructor_eo.f90` in Read First with
  `simple_halfmap_diagnostics.f90`.
- `.github/skills/simple-refine3d/SKILL.md`,
  `references/refine3d-code-map.md`, `references/volume-postprocessing.md`:
  repoint the module reference to `simple_halfmap_diagnostics.f90` and the
  pair I/O to `restore_state_from_parts`.
- `.github/skills/simple-abinitio3d-importance-sampling/SKILL.md`: same
  repointing.
- `.github/skills/simple-frac-update-trailing/SKILL.md` and
  `references/frac-update-contract.md`: `reconstructor_eo%read_eos_parallel_io`
  is stale twice over; the current owner is
  `read_gridding_pair_accumulators` contained in `restore_state_from_parts`.
