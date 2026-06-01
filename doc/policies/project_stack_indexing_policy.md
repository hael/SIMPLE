# Project Stack Indexing Policy

This policy defines the indexing contract between SIMPLE project metadata,
particle rows, and physical particle stack files. It is intended to remove
ambiguity around `merge_projects`, `prune_project`, and downstream stack readers
such as 2D and 3D analysis.

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

Modern project outputs should write positive `indstk` values for particle rows.
Missing or non-positive `indstk` values are legacy compatibility cases.
For old stream projects, even a positive `indstk` may be unreliable when the
stack row does not carry `nptcls_stk`.

When a legacy fallback is needed, or when a project-rewriting operation cannot
trust the input `indstk`, `fromp/top` define the project rows belonging to the
stack row, while the stack file itself is 1-based. Therefore the fallback
physical stack index is:

```text
indstk = particle_project_row - fromp + 1
```

The fallback must never use the raw project particle row as the stack image
index unless `fromp == 1`. If `state = 0` rows exist, they still occupy project
rows and therefore still occupy positions in this fallback mapping.

Required invariant:

```text
1 <= indstk <= nptcls_stk
```

When `nptcls_stk` is absent in a legacy project and no physical stack inspection
has been done, the fallback physical stack count is:

```text
nptcls_stk = top - fromp + 1
```

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
indices when they are credible and repair legacy old-stream indexing when they
are not.

When merging particle projects:

- Remap `stkind`, because stack rows are renumbered in the merged project.
- Remap `fromp` and `top`, because project particle rows are concatenated.
- Preserve `nptcls` for each stack row.
- Treat input `indstk` as a candidate, not as automatically correct.
- If the source stack row has `nptcls_stk`, preserve positive `indstk` values
  that are within `1..nptcls_stk`.
- If the source stack row lacks `nptcls_stk`, treat the source as legacy and
  derive `indstk` from `source_particle_project_row - source_fromp + 1`, even
  when an `indstk` value is present.
- If `indstk` is missing, non-positive, or out of range, derive it from
  `source_particle_project_row - source_fromp + 1`.
- Preserve `nptcls_stk` when present, because the physical stack files are not
  rewritten.
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

The fallback formula is applied in the source project's project-row domain
before the output `fromp/top` remapping. This prevents the merged global
particle row from being used as a stack index.

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
- Treat input `indstk` as a candidate, not as automatically correct.
- If the source stack row has `nptcls_stk`, preserve positive `indstk` values
  that are within `1..nptcls_stk`, because the physical image positions in the
  original stack files have not changed.
- If the source stack row lacks `nptcls_stk`, treat the source as legacy and
  derive `indstk` from `old_particle_project_row - old_fromp + 1`, even when an
  `indstk` value is present.
- If `indstk` is missing, non-positive, or out of range, derive it from the
  legacy fallback `old_particle_project_row - old_fromp + 1`.
- Preserve `nptcls_stk` when present, because the physical stack files have not
  changed.
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

Any code that reads particle images from a stack file must use `indstk` as the
physical image index when the stack row has `nptcls_stk` and `indstk` is
present, positive, and within range.

If `nptcls_stk` is absent, or if `indstk` is missing, non-positive, or out of
range, readers may use the legacy fallback:

```text
indstk = particle_project_row - fromp + 1
```

The reader must first resolve the particle's `stkind`, read `fromp/top` from
that stack row, and verify that the particle row lies inside that range. Readers
must not use the raw particle row index as the physical stack-file image index
except in the special case where `fromp == 1`.

## Validation Policy

Project validation should check:

- Stack project ranges are contiguous and cover the particle rows.
- `top - fromp + 1 == nptcls`.
- Particle `stkind` values identify valid stack rows.
- New project-writing code should write `nptcls_stk` on stack rows when
  particle rows are present.
- When `nptcls_stk` is missing in a legacy project, validation may inspect the
  physical stack file or use `top - fromp + 1` as the fallback stack count.
- New project-writing code should write positive particle `indstk` values.
- Missing, non-positive, or untrusted particle `indstk` values may be accepted
  through the `particle_project_row - fromp + 1` fallback.
- Particle `indstk` values are not larger than `nptcls_stk` when
  `nptcls_stk` is present.
- If `nptcls_stk` is absent, validation must not assume existing `indstk`
  values are correct; merge/prune-style repair should derive them from
  `fromp/top`.

Validation must not silently convert physical `indstk` values into project-row
indices.

The `validate_projfile` program applies this policy to an input project and
writes `input_name_validated.simple`. It reports warnings, errors, and repairs
encountered during validation, but writes a best-effort repaired project so that
legacy stream outputs can be normalized before downstream 2D or 3D analysis.

## Test Policy

Tests for merge, prune, and stack readers should include pruned-style projects
where:

```text
nptcls < nptcls_stk
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

Tests should also include old-stream-style projects where `nptcls_stk` is
absent and `indstk` is missing, zero, or wrong. Merge and prune should repair
those rows with:

```text
indstk = source_particle_project_row - source_fromp + 1
```
