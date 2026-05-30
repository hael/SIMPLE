# Refactoring: Merge Projects With Distinct CTF Models

## Implementation Status

The first implementation phase is in place:

- `src/main/commanders/simple/simple_commanders_project_core.f90` routes
  `merge_projects` through the dedicated project table parameter `projtab` and
  requires `projfile_merged` as the explicit output project.
- `src/main/ui/simple/simple_ui_project.f90` exposes the file-table-only
  interface and no longer treats `merge_projects` as an in-place two-project
  project command.
- `src/fileio/simple_projfile_utils.f90` has a reusable
  `merge_selected_project_files` helper, following the same file-array merge
  shape as `merge_chunk_projfiles`.
- Input `.simple` files are inspected through existing `sp_project`
  segment-info and segment-read methods, avoiding the full project read path
  that rewrites `projinfo`.
- If a listed project file contains only metadata, the helper fails with the
  offending filename and, when possible, a hint pointing at a nearby
  data-bearing numbered-stage project file.

The helper is now a generic project-field merger. It accepts any data-bearing
project shape whose populated data segments match across all inputs:

- movie/micrograph-only projects through `os_mic`
- stack-only projects through `os_stk`
- particle projects through `os_ptcl2D` and/or `os_ptcl3D`
- class/output/optics-bearing projects when those segments are present in all
  inputs

Top-level project files that contain only metadata segments such as `projinfo`
and `compenv` are rejected because there is no data segment to merge.

The helper currently:

- validates all populated data segments as all-or-none across inputs
- validates stack box and sampling distance when stack rows carry those fields
- validates micrograph sampling distance when mic-only rows carry `smpd`
- validates row-level CTF-model fields on `os_mic` and `os_stk` when CTF is
  enabled
- validates particle `stkind` when particles are merged with stacks
- preserves complete rows with `transfer_ori`
- precomputes output row offsets and parallelizes independent row
  transfer/remap loops for mic, stack, class, output, optics, `ptcl2D`, and
  `ptcl3D` segments
- parallelizes large particle validation scans and row-level `ogid` discovery
  while reporting deterministic first failing rows
- preserves exact `os_ptcl2D%state` values when `ptcl2D` exists
- remaps stack `fromp` / `top` only when particle rows are present, remaps
  particle `stkind`, particle/class `class`, and row-level `ogid`
- remaps `ogid` independently of whether `os_optics` exists
- copies/remaps `os_optics` only when the optics segment is part of the matching
  input field set
- rebuilds `os_ptcl3D` from merged `os_ptcl2D` only for the historical
  ptcl2D-only input shape, then removes 2D-clustering fields from `ptcl3D`;
  this synthetic `ptcl3D` row construction is also parallelized

No `os_optics` backfill path was added. `sp_project%get_ctfparams` already
reads `smpd`, `ctf`, `kv`, `cs`, `fraca`, and `phaseplate` from `os_stk` for
particle analysis and particle-specific defocus from particle rows, so the
remaining resolver-related work is an audit of direct CTF construction outside
that project API plus dedicated merge tests.

## Goal

Support N-project merging for SIMPLE projects that may come from different
microscopes or CTF models but have compatible analysis data. The immediate use
case is merging independently processed projects that have already been scaled
to the same sampling distance, but the merge utility should not be hard-coded
to a 2D-selected particle project shape.

The general contract is:

1. Every input has the same populated project data segments.
2. Row counts may differ.
3. Segment rows are preserved verbatim except for required merge-local
   rewrites.
4. Cross-project local namespaces and foreign keys are remapped.
5. CTF-model values come from authoritative row-level fields, not from
   `os_optics`.

For CTF-aware rows, the merger must preserve:

- `smpd`
- `kv`
- `cs`
- `fraca`
- `ctf`
- `phaseplate`

`fraca` is deliberately included. It is not a microscope hardware constant, but
it changes the CTF and therefore must be merged and resolved with `kv`, `cs`,
and `smpd`.

## Non-Goals

Do not change heterogeneous movie, micrograph, or particle import in this
refactor. The current scenario assumes data have already been imported and
processed in separate internally consistent projects when needed.

Do not make `os_optics` required. Some datasets do not have optics-group
metadata, depending on collection and import path.

Do not use `os_optics` as a source of truth or a repair source for missing
authoritative row-level CTF fields.

Do not rescale or re-extract during merge. For stack/particle projects, the
current target requires input stacks to have already been scaled/extracted to
the same sampling distance and compatible box size. Mismatches should fail
validation with an actionable message.

Do not attempt to scientifically reconcile per-project class averages,
reconstructions, or output artifacts. If `cls2D`, `cls3D`, or `out` segments are
part of the matching input field set, the generic merger may concatenate and
remap their project rows, but it does not create a newly optimized combined
analysis product.

## Target Policy

The authoritative CTF-model representation is the row-level bundle used by each
processing stage:

- `os_mic`: source of truth for movie/micrograph CTF fitting history.
- `os_stk`: source of truth for particle-stack-level CTF constants used by
  CTF-aware 2D/3D analysis.
- `os_ptcl2D` / `os_ptcl3D`: source of truth for particle-level CTF values such
  as defocus, astigmatism angle, phase shift, `stkind`, and selection state.

`os_optics` rows are optional metadata. Row-level `ogid` values, when present,
are optics-group assignments and must be remapped during merge even if the input
project has no `os_optics` segment. The merge must still work from row-level
CTF-model fields when no durable optics rows exist.

Because this generic merger requires matching populated data segments, mixed
presence of `os_optics` is treated as a project-shape mismatch in this phase. If
none of the input projects have `os_optics`, the merge is still valid and any
row-level `ogid` values are remapped. If future workflows need mixed optional
optics metadata, that should be an explicit extension rather than an implicit
backfill from optics rows.

The merge must not use `os_optics` to repair missing CTF-model fields. For this
workflow the authoritative `mic`, `stk`, and particle rows are expected to
already contain the values needed for CTF fitting and downstream analysis. A
missing authoritative field is a validation error unless CTF is explicitly
disabled for that row.

All active stacks entering a merged 2D/3D project should have one effective
analysis sampling distance. Different sampling distances should be handled by
an explicit scaling or re-extraction step before this merge.

## Original Code Findings

Before this refactor, the `merge_projects` commander in
`src/main/commanders/simple/simple_commanders_project_core.f90` was not the
right abstraction:

- It supported exactly two projects by hard-coding `nprojs = 2`.
- It read `projfile` and `projfile_target`, appended the target into the first
  project with `sp_project%append_project`, and wrote back to `projfile`.
- It detected different box or sampling distance and then ran `reextract`.
  Heterogeneous CTF-model merging should instead require compatible data for
  the current scenario.
- It delegated row merging to `sp_project%append_project`, which mixed merge
  policy with project mutation.
- It remapped `ogid` only inside the `os_optics` branch, missing projects that
  carry row-level optics-group assignments without an `os_optics` table.

The stream `cluster2D_subsets` path has a better framework for project
aggregation:

- `src/main/stream/simple_stream_cluster2D_subsets.f90` tracks project files in
  `rec_list` / `chunk_rec` and passes arrays of project filenames into a merge
  utility.
- `src/fileio/simple_projfile_utils.f90::merge_chunk_projfiles` accepts an
  array of project files, allocates output project segments, transfers
  orientations with `transfer_ori`, remaps stack and particle indices, and
  optionally writes the merged project.
- That framework keeps commanders thin and puts reusable project-file merge
  mechanics outside the commander.

The refactor reuses this framework shape: `merge_projects` builds an ordered
list of source project files from `projtab`, then calls a project-file merge
utility in `simple_projfile_utils`.

## Project Field Contract

The merger should determine the populated project data segments for every input
project:

- `os_mic`
- `os_stk`
- `os_ptcl2D`
- `os_ptcl3D`
- `os_cls2D`
- `os_cls3D`
- `os_out`
- `os_optics`

Every segment must be either populated in all inputs or empty in all inputs.
This is the meaning of "the individual fields match" for this phase. A
stack-only project can merge with another stack-only project; a movie-only
project can merge with another movie-only project; a particle project can merge
with another particle project that has the same populated particle/support
segments. A stack-only project should not be silently merged with a
stack-plus-particle project.

Metadata-only segments such as `projinfo`, `jobproc`, and `compenv` are not a
data shape. The merged project should copy them from the first input and update
the output project filename.

Complete source rows should be transferred with `transfer_ori`. The merger
should mutate only fields whose values are local to the source project and
therefore cannot remain verbatim in the merged namespace.

## Particle State Contract

When `os_ptcl2D` is present, `state` is an authoritative input field. The
merger must preserve the `state` value for every copied `ptcl2D` row verbatim
after append and index remapping. It must not reset all particles to active,
infer `ptcl2D` state from class rows, drop deselected rows, or normalize
positive state labels unless a separate explicit pruning/relabeling command is
requested.

If `os_ptcl2D` is present and `os_ptcl3D` is absent, the helper may preserve the
historical stream behavior of creating `os_ptcl3D` from the merged `ptcl2D`
rows and deleting 2D-clustering fields. If both particle segments are present,
both should be transferred and remapped as independent authoritative segments.

## Index Contracts

`stkind` is not a value to preserve verbatim across projects. It is a foreign
key from each particle row into that project's local `os_stk` table. During an
N-project merge, stack rows are concatenated into one output `os_stk` table, so
every copied particle row with stacks present must have `stkind` rewritten to
the corresponding output stack index.

For stack/particle projects, the merger must:

1. Transfer each source stack row.
2. Rewrite stack `fromp` / `top` into the output particle index range when
   those fields are present.
3. Require valid particle `stkind` values when particles are merged with
   stacks.
4. Rewrite each transferred particle `stkind` from the source stack index to
   the output stack index.

For stack-only projects, stack rows are transferred as-is. The absence of
particles should not make the stack field unmergeable, and stack-only merges do
not rewrite `fromp` / `top`.

When class segments are present, class rows and particle `class` assignments
must be offset into the merged class namespace. The helper should not infer
selection from class rows.

## Optics Assignment Contract

`ogid` is a row-level assignment namespace. Values from different source
projects can collide, so the merged project must rewrite positive `ogid` values
into a single output namespace.

The remapping trigger is the presence of row-level `ogid` values on any copied
segment, not the presence of `os_optics`. For each source project, collect the
maximum positive `ogid` namespace used by its rows, allocate a collision-free
output offset, and rewrite every copied row that has an `ogid` field.

If `os_optics` is part of the matching input field set, append/remap those rows
with the same output `ogid` namespace and update `ogname` when appropriate. If
no source project has `os_optics`, no durable optics rows need to be synthesized.

## Core Refactor

Add or finish a project-level CTF-model resolver in `src/main/project`, with
behavior roughly equivalent to:

1. Identify the relevant row for `oritype` and particle index.
2. Read CTF-model parameters from the stage's primary row: `os_mic` for
   micrographs, `os_stk` for stacks and particles.
3. Do not consult `os_optics` to fill missing CTF-model values.
4. Require the returned model to be complete unless CTF is explicitly disabled;
   missing authoritative fields should fail validation rather than be inferred
   from metadata.
5. Read image- or particle-specific defocus from the current particle-row
   locations.
6. Return a complete `ctfparams`.

Update `sp_project%get_ctfparams` to call this resolver if direct duplicated
logic remains. That change covers the major CTF-aware 2D/3D consumers because
polar-FT CTF matrix generation, class averaging, and reconstruction helper
paths already call `get_ctfparams`.

Add validation helpers:

- `validate_project_field_shape`: ensure all inputs have matching populated
  data segments.
- `validate_ctf_model_rows`: ensure every active `mic` or `stk` row that needs
  CTF has complete row-level `smpd`, `kv`, `cs`, `fraca`, `ctf`, and
  phase-plate state.
- `validate_particle_ctf_rows`: ensure particles with CTF enabled through their
  stack have defocus/astigmatism/phase-shift fields required by their CTF mode.
- `validate_ptcl2D_state_rows`: when `ptcl2D` exists, record the state vector
  so the merger can assert that copied rows are unchanged.
- `build_ogid_remap`: collect row-level `ogid` assignments per source project
  and allocate a collision-free output `ogid` namespace independent of
  `os_optics`.
- `validate_common_analysis_smpd`: enforce or report the common sampling
  distance expected by current 2D/3D parameter derivation when stack/mic rows
  carry `smpd`.
- `validate_common_particle_box`: enforce identical particle image dimensions
  when stack rows carry `box`.

## Project Merge Changes

`merge_projects` should remain a thin commander over the reusable N-project
merge utility:

1. Accept only a project file table for N input projects through `projtab`.
   Do not preserve the legacy two-project `projfile` plus `projfile_target`
   interface for this refactor.
2. Require an explicit output project name, `projfile_merged`, so the merged
   project is not written over one of the inputs.
3. Resolve a relative `projfile_merged` against the execution directory without
   calling an existence-based canonicalizer; the output file may not exist yet.
4. Build an ordered `project_fnames(:)` array, analogous to the stream
   `cluster2D_subsets` / `merge_chunk_projfiles` call sites.
5. Read each source through existing `sp_project` segment-info and
   segment-read methods so `merge_projects` can inspect arbitrary populated
   fields without rewriting source project metadata.
6. Call the shared project-file merge helper:
   `merge_selected_project_files(project_fnames, projfile_out, merged_proj, ...)`.
7. Keep the commander responsible only for CLI validation, input normalization,
   writing, and `simple_end`.

The shared helper should:

1. Read all source projects and determine their populated project field shape.
2. Fail if the populated data segments do not match across inputs.
3. Fail if no mergeable data segments exist beyond metadata.
4. Validate row-level CTF-model data where CTF is enabled.
5. Validate identical analysis sampling distance and particle box dimensions
   when those fields exist. Do not run `reextract` in this merge path.
6. Allocate every populated output segment up front.
7. Copy source rows with `transfer_ori`.
8. Remap only merge-local fields: stack `fromp` / `top` when particles are
   present, particle `stkind`, particle/class `class`, and row-level `ogid`.
9. Preserve row-level CTF-model values by transferring complete `os_mic` and
   `os_stk` rows. Do not overwrite `kv`, `cs`, `fraca`, `ctf`, `phaseplate`, or
   `smpd` from global parameters.
10. Preserve particle-level CTF values and exact `os_ptcl2D%state` flags by
    transferring complete particle rows.
11. Remap row-level `ogid` assignments for every copied segment that carries
    `ogid`, independent of whether the source project has `os_optics`.
12. Preserve/remap `os_optics` only when it is part of the matching input field
    set.
13. If a source project has incomplete authoritative row-level CTF fields, fail
    validation even if corresponding values exist in `os_optics`.

Performance policy:

- The merge is allowed to keep all input and output metadata rows in memory;
  this refactor is not optimized for memory-constrained machines.
- The merge must not rewrite movie, micrograph, stack, or particle image data.
- The hot row-transfer/remap loops should be OpenMP-parallel where each
  iteration writes to a precomputed, distinct output row.
- Large validation and namespace-discovery scans should also use OpenMP when
  they can reduce to deterministic scalar results, such as the first invalid
  particle row or the maximum row-level `ogid` in a source project.
- The final project metadata write is still a single output-file serialization
  step through the existing project writer.

This keeps merge behavior predictable and avoids accidental regrouping.

## Optional Optics Metadata

When `os_optics` exists for every source project, it should be kept consistent
with row-level values and remapped with the same `ogid` namespace used for the
rows. When no source project has `os_optics`, the row-level `ogid` assignment
should still be remapped and the project should still be complete and
scientifically valid.

`os_optics` is not a repair source for this refactor. It may be preserved and
remapped as metadata, but it should not be read to fill missing `smpd`, `kv`,
`cs`, `fraca`, CTF flag, or phase-plate state on authoritative rows.

If future STAR export needs an optics table from a merged project that lacks
`os_optics`, the STAR adapter can synthesize a temporary export-only optics
table from row-level CTF-model groups. It should not require or persist
`os_optics` just to export.

## 2D and 3D Analysis

Once `get_ctfparams` uses the resolver consistently, most CTF-aware 2D/3D
analysis should pick up mixed `kv`, `cs`, and `fraca` values without broad
algorithm changes.

Still audit these areas explicitly:

- `src/main/pftc/simple_polarft_ctf.f90`
- `src/main/class/simple_classaverager_restore.f90`
- `src/main/volume` reconstruction helpers that apply CTF
- any direct construction of `ctf(params%smpd, params%kv, params%cs, ...)`
  outside simulation or validation code

Global `params%kv` and `params%cs` should not be used for real project CTF
evaluation. They can remain defaults for simulation, validation, and
homogeneous import.

## Tests

Add project-level and commander-level tests for matching field shapes:

- `merge_projects` accepts an N-project file table and routes through the shared
  project-file merge helper rather than an in-commander append loop.
- Stack-only projects merge when all inputs have `os_stk` and no particle
  segments.
- Movie/micrograph-only projects merge when all inputs have `os_mic` and no
  stack/particle segments.
- Particle projects merge when all populated data segments match.
- Metadata-only projects fail with a clear error that identifies the offending
  `projtab` entry and hints at a nearby data-bearing stage project when one is
  present.
- Projects with mismatched populated data segments fail with a clear field-shape
  error.
- Input projects can have identical analysis `smpd` and compatible boxes but
  different `kv`, `cs`, and/or `fraca` on their `os_stk` rows.
- `merge_projects` preserves row-level mic/stack CTF-model values.
- `merge_projects` preserves the exact `os_ptcl2D%state` vector from each input
  project in the corresponding output row range, including `0` entries and any
  positive labels.
- `merge_projects` preserves particle-level defocus fields and remaps `stkind`
  correctly after append.
- For every output stack row, all output particles in that stack's `fromp` /
  `top` range have `stkind` equal to the output stack index.
- Class rows and particle `class` values are offset when class segments are
  part of the matching input field set.
- `sp_project%get_ctfparams('ptcl2D', iptcl)` returns different `kv`/`cs` for
  particles originating from different projects.
- Cluster2D and refine3D CTF matrix construction call through the resolver.
- The merger works when no source project has `os_optics`.
- The merger remaps row-level `ogid` assignments even when no source project
  has `os_optics`.
- If two source projects both use `ogid=1`, the output rows from those projects
  receive distinct output `ogid` values even when neither source has
  `os_optics`.
- A project with incomplete row-level CTF fields fails validation even if
  `os_optics` contains matching values.
- A mixed-sampling or mixed-box project fails validation for this merge path,
  with instructions to rescale/re-extract before merging.
- The old two-project `projfile` / `projfile_target` interface is not supported
  by this refactor.

## Suggested Migration Plan

1. Add resolver and validation helpers with no behavior change.
2. Route `get_ctfparams` through the resolver without introducing an
   `os_optics` path for missing authoritative fields.
3. Add a reusable project-file merge helper alongside `merge_chunk_projfiles`.
4. Refactor `merge_projects` to build a project filename array from `projtab`
   and call the shared helper.
5. Generalize the helper from a 2D-selected particle-only merge to a matching
   project-field merge.
6. Make stack/particle merge paths require identical `smpd` and box size when
   those fields are present; leave scaling/re-extraction as an explicit
   pre-merge step.
7. Preserve row-level CTF-model fields and exact `os_ptcl2D%state` flags during
   row transfer.
8. Remap local indices and namespaces during row transfer independently of
   whether `os_optics` exists.
9. Add merge tests for stack-only, mic-only, particle, and heterogeneous CTF
   particle projects.
10. Audit direct CTF-field lookups in 2D/3D paths and route them through the
    resolver where needed.
11. Leave heterogeneous import support for a separate future refactor.

## Open Questions

- Should mixed optional `os_optics` presence be supported by explicitly dropping
  optics rows from the merged output, or should matching field shape remain the
  rule?
- Should phase-plate state be required to be stack-level for merged particle
  projects, or should per-particle phase-plate state be allowed?
- Should selected class provenance be retained as optional annotations after
  class/output segments are concatenated?
