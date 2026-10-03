# Project Stack Indexing Policy

This policy defines the indexing contract between SIMPLE project metadata,
particle rows, and physical particle stack files. It is intended to remove
ambiguity around `merge_projects`, `prune_project`, `fix_projfile`, and
downstream stack readers such as 2D and 3D analysis.

## Core Rule

`fromp`, `top`, and `nptcls` in the `os_stk` field refer to the project, not
to the physical stack file.

`indstk` refers to the physical stack file, not to the project.

These are separate index domains and must not be interchanged.

## Index Domains

For projects with particle rows and stack rows:

```text
ptcl row index        row in os_ptcl2D or os_ptcl3D
stkind                row in os_stk for the particle's stack
fromp/top/nptcls      current project particle range/count for an os_stk row
indstk                physical 1-based image index inside the stack file
nptcls_stk            physical image count in the stack file
```

The stack path stored on an `os_stk` row identifies the physical stack file. A
particle image is read by resolving:

```text
particle row -> stkind -> os_stk row -> stack path
particle row -> indstk -> physical image number in that stack
```

## Field Semantics

### `os_stk%fromp`, `os_stk%top`, and `os_stk%nptcls`

These fields describe the particle rows currently present in the project for a
stack row.

They do not describe physical image positions in the backing stack file.

They do not inherently mean "active particles" or "`state > 0` particles." They
count project particle rows assigned to that stack row. If a project still
contains rows with `state = 0`, those rows are part of the project and
therefore are included in `fromp/top/nptcls`.

If a pruning operation removes rows from the project, then removed rows are no
longer counted. In that case `fromp/top/nptcls` describe the post-prune project,
not the original project and not the physical stack file.

Required invariants when particle rows are present:

```text
top - fromp + 1 == nptcls
stack ranges are contiguous in project-row space
stack ranges cover the particle rows in the project segment
particle rows in a stack range have the matching stkind
```

### `os_stk%nptcls_stk`

`nptcls_stk` describes the number of physical images in the backing stack file.

It may be larger than `nptcls` after pruning or selection if the project
metadata points to the original stack file.

It should equal `nptcls` only when the backing stack file itself contains
exactly the particle rows present in the project, for example after writing a
new stack file from those rows.

### `os_ptcl2D/os_ptcl3D%stkind`

`stkind` identifies the `os_stk` row associated with a particle. Because stack
rows are project metadata, `stkind` is a project-local index and may be remapped
when projects are merged or stack rows are renumbered.

### `os_ptcl2D/os_ptcl3D%indstk`

`indstk` identifies the particle's physical image number inside the stack file
referenced by its `stkind`.

Every particle row must carry a positive `indstk`, and every stack row with
particle rows must carry `nptcls_stk`. Both are required from release 4 on:

```text
1 <= indstk <= nptcls_stk
```

There is no run-time fallback. Earlier releases derived a missing `indstk` from
the project rows (`particle_project_row - fromp + 1`) and a missing
`nptcls_stk` from the range count. That mapping is only correct while the stack
file holds exactly one image per project row, and old stream projects carried
`indstk` values that were wrong even when present; the fallback hid that bug. A
project that lacks these fields stops with an error in mapping, pruning,
merging and STAR export. `fix_projfile` (see the Validation Policy) converts
such a project once, where the physical indices can be proven.

## Selection State

Particle `state` is not part of the stack-index definition. It is a selection
or activity flag on a project particle row.

Use `state > 0` when a workflow needs the number of active particles. Use
`fromp/top/nptcls` when a workflow needs the number of particle rows currently
present in the project for a stack.

Do not infer active-particle counts from `nptcls` unless the project operation
has explicitly removed all inactive rows.

## Merge Policy

`merge_projects` concatenates project metadata. It must preserve physical stack
indices.

When merging particle projects:

- Remap `stkind`, because stack rows are renumbered in the merged project.
- Remap `fromp` and `top`, because project particle rows are concatenated.
- Preserve `nptcls` for each stack row.
- Preserve `indstk`. A source stack row without `nptcls_stk`, or a particle
  whose `indstk` is missing or outside `1..nptcls_stk`, stops the merge; run
  `fix_projfile` on that source first.
- Preserve `nptcls_stk`, because the physical stack files are not rewritten.
- Preserve stack paths and row-level stack metadata, including CTF-model fields.
- Remap row-level `ogid` values as project metadata.
- Drop analysis products such as `cls2D`, `cls3D`, and `out` unless a future
  policy explicitly defines how to merge them.

The important merge invariant is:

```text
project row indices may change
stack row indices may change
physical image indices must not change
```

Because `indstk` is carried over unchanged, the merged global particle row is
never used as a stack index.

## Prune Policy

Pruning must distinguish between metadata-only pruning and materializing
pruning.

### Metadata-only prune

If pruning removes particle rows from the project but keeps the original stack
files:

- Remove deleted rows from `os_ptcl2D` and `os_ptcl3D`.
- Update `os_stk%fromp`, `os_stk%top`, and `os_stk%nptcls` to describe the
  post-prune project.
- Remap `stkind` if stack rows are renumbered.
- Preserve `indstk`, because the physical image positions in the original
  stack files have not changed. A retained particle without a valid `indstk`
  stops the prune before any row is moved.
- Preserve `nptcls_stk`, because the physical stack files have not changed.
- Preserve stack paths.

This is the common case that protects downstream readers from using project row
numbers as stack slice numbers.

### Materializing prune

If pruning writes new physical stack files containing only the retained
particles:

- Update stack paths to the new physical files.
- Rewrite `indstk` to the physical image positions in the new stack files,
  normally `1..nptcls`.
- Set `nptcls_stk` to the physical image count in the new stack file.
- Set `fromp/top/nptcls` to the project row ranges.
- The number of particle rows is the same as the physical image count:
  `nptcls = nptcls_stk`

This is the only case where pruning should rewrite `indstk` based on the new
stack-file order.

### Invalid prune state

The following state is invalid:

```text
stack path still points to the original physical stack file
indstk has been recomputed from project-row order
```

That state corrupts the particle-to-image mapping and can make distributed
workers request image indices that do not correspond to the intended particles.

## Reader Policy

Any code that reads particle images from a stack file must resolve the
particle's `stkind` and use `indstk` as the physical image index
(`map_ptcl_ind2stk_ind`). A missing `nptcls_stk`, or an `indstk` that is
missing or outside `1..nptcls_stk`, is an error. Readers must never use a
project particle row as the physical stack-file image index.

## Validation Policy

Project validation should check:

- Stack project ranges are contiguous and cover the particle rows.
- `top - fromp + 1 == nptcls`.
- Particle `stkind` values identify valid stack rows.
- Every stack row with particle rows has `nptcls_stk`.
- Every particle row has `1 <= indstk <= nptcls_stk`.

Project-writing code must write both. Validation must not silently convert
physical `indstk` values into project-row indices.

The `fix_projfile` program brings a project from an earlier release to this
contract and writes `input_name_fixed.simple`:

- It repairs stack ranges from the reference particle segment (`ptcl2D` when
  present), then gives every row of both segments the `stkind` of the stack
  whose range owns it; a stored `stkind` naming another stack is replaced, and a
  row no stack owns is an error.
- A project with particle rows but no stack rows is an error: its particles
  cannot be mapped to images.
- A stack row without `nptcls_stk` gets the image count from its file header.
  Such a row predates the stack-index fix, so its stored `indstk` values are not
  trusted: when the stack holds exactly one image per project row, every
  particle gets `indstk = particle_project_row - fromp + 1`, checked against
  `1..nptcls_stk`; otherwise the physical images cannot be recovered and the
  stack is reported.
- A stack row that already has `nptcls_stk` keeps its `indstk` values; a
  missing or out-of-range one is reported.
- It reports every repair and error, and writes the fixed project only when no
  error remains. Reading and rewriting the project also brings the particle
  records to the current width.

## Test Policy

Tests for merge, prune, and stack readers should include pruned-style projects
where:

```text
nptcls <= nptcls_stk
fromp/top describe current project rows
indstk contains non-contiguous physical image indices
```

For example, a project stack row may have:

```text
fromp = 1
top = 3
nptcls = 3
nptcls_stk = 6
particle indstk values = [1, 4, 6]
```

After merging this project with another, the second project's `fromp/top` and
`stkind` values may be remapped, but its `indstk` values must remain physical
indices into its original stack file unless new physical stack files are
written.

Tests should also include an earlier-release project where `nptcls_stk` is
absent and `indstk` is missing, zero, or wrong, on a stack that holds one image
per project row: `fix_projfile` must take `nptcls_stk` from the stack header and
set `indstk = particle_project_row - fromp + 1`, after which the project maps
cleanly.
