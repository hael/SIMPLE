# Project Stack Indexing Code Review

Date: 2026-06-03

This review compares the current written stack-indexing policy against the code
paths that create, merge, prune, validate, export, or consume SIMPLE project
particle-stack metadata.

The key policy is `doc/policies/project_stack_indexing_policy.md`.

## Executive Summary

The written policy is coherent. The core reader and the main `merge_projects`
implementation mostly implement it correctly.

The remaining risk is in project writers and rewriters that still perform
partial remaps by hand. Several stream and chunk-merge paths remap project
indices (`fromp`, `top`, `stkind`) without normalizing the physical stack index
(`indstk`) or the physical stack count (`nptcls_stk`). Those paths can produce a
project that looks internally plausible but later asks a stack reader for image
`N + 1` from a physical stack containing only `N` images.

This is consistent with downstream failures such as:

```text
ERROR! index i out of range: 50 / 49; simple_discrete_stack_io.f90; line: 61
```

The physical stack reader is only reporting the symptom. The corruption is
introduced earlier when project metadata and physical stack metadata are allowed
to drift apart.

There is one especially important nuance: after a metadata-only prune, such as
`map_cavgs_selection` followed by `prune=yes`, `fromp/top` no longer describe a
contiguous physical image range in the backing stack. They describe only the
remaining project rows. From that point on, the project must carry valid
`indstk` values. Falling back to `particle_project_row - fromp + 1` is no longer
safe unless the project row range is known to still mirror the physical stack
order.

## Policy Contract

The policy separates two index domains:

- `os_stk%fromp`, `os_stk%top`, and `os_stk%nptcls` are project-row metadata.
- `os_ptcl2D/os_ptcl3D%indstk` is the physical image index inside the stack file.
- `os_stk%nptcls_stk` is the physical image count in the stack file.
- `os_ptcl2D/os_ptcl3D%stkind` is the project-local stack row index.

Required invariants when particle rows and stack rows are present:

```text
top - fromp + 1 == nptcls
stack ranges are contiguous in project-row space
stack ranges cover the particle rows in the project segment
particle rows in a stack range have matching stkind
1 <= indstk <= nptcls_stk
```

When `indstk` cannot be trusted and the original stack file is still being used,
the fallback is:

```text
indstk = particle_project_row - source_fromp + 1
```

That fallback is in the source project's project-row domain. It must not use the
output project row after merge/remap.

The fallback is valid only when the source project row range still corresponds to
the physical stack order. Once a metadata-only prune removes rows while keeping
the original stack file, `fromp/top` becomes a project-only range and must not be
used to reconstruct physical indices for missing or invalid `indstk` values.

`state` is not part of the stack-indexing contract. `state = 0` rows still count
in `fromp/top/nptcls` while they remain present in the project.

## Findings

### P1: metadata-only `prune=yes` makes `fromp/top` unusable as a physical fallback

Files:

- `src/main/commanders/simple/simple_commanders_cavgs.f90`
- `src/main/project/simple_sp_project_ptcl.f90`

Relevant lines:

- `simple_commanders_cavgs.f90:831-833`
- `simple_commanders_cavgs.f90:295-297`
- `simple_commanders_cavgs.f90:409-411`
- `simple_commanders_cavgs.f90:627`
- `simple_sp_project_ptcl.f90:430-532`

Several class-average selection commands call `map_cavgs_selection`, then
optionally call `prune_particles` when `prune=yes`. This is a metadata-only
prune: it removes particle rows from the project and rewrites `fromp/top/nptcls`,
but it does not rewrite the physical particle stacks.

After this operation, the output stack row range is dense in project-row space,
but the physical images may be sparse. For example:

```text
original physical stack images: 1 2 3 4 5 6
selected project rows/images:  1 4 6

after metadata-only prune:
    fromp = 1
    top = 3
    nptcls = 3
    nptcls_stk = 6
    required indstk values = [1, 4, 6]
```

In this state, `iptcl - fromp + 1` would produce `[1, 2, 3]`, which is a valid
project-row sequence but the wrong physical image sequence.

Recommendation:

Any project with `nptcls_stk > nptcls` or `nptcls_stk > top - fromp + 1` should
be treated as physically sparse. In that state:

- valid existing `indstk` is mandatory,
- missing, zero, or out-of-range `indstk` should hard fail,
- readers, merge, validation, and metadata-only prune should not fall back to
  `iptcl - fromp + 1`,
- the only automatic repair is a materializing prune that rewrites the stack and
  resets `indstk = 1..nptcls` and `nptcls_stk = nptcls`.

This should be promoted to the formal policy.

### P1: the canonical reader fallback is too permissive for sparse project ranges

File: `src/main/project/simple_sp_project_ptcl.f90`

Relevant lines: 58-110

`map_ptcl_ind2stk_ind` currently reads `indstk` whenever present, regardless of
whether `nptcls_stk` is present. If `indstk` is missing or non-positive, it falls
back to `iptcl - fromp + 1`:

```fortran
if( ptcl_field%isthere(iptcl, 'indstk') )then
    ind_in_stk = ptcl_field%get_int(iptcl, 'indstk')
endif
if( ind_in_stk < 1 )then
    call set_indstk_from_range
endif
```

This is unsafe after metadata-only pruning. If `fromp/top` is already a
post-prune project range, the fallback returns an in-range but semantically
wrong physical index. It can silently read the wrong particles rather than
throwing an out-of-range error.

It is also more permissive than the written policy for old stream projects:
positive `indstk` is trusted even when `nptcls_stk` is absent, although the
policy says old stream `indstk` may be unreliable without `nptcls_stk`.

Recommendation:

The reader should compute the project range count:

```text
nptcls_project = top - fromp + 1
```

Then:

- if `nptcls_stk` is present and `nptcls_stk > nptcls_project`, require valid
  `indstk` and do not fallback,
- if `nptcls_stk` is present and `nptcls_stk == nptcls_project`, fallback is
  acceptable for missing legacy `indstk`,
- if `nptcls_stk` is absent, do not trust positive `indstk`; use fallback only as
  a legacy unverified mapping and report/repair upstream when possible.

### P1: `merge_chunk_projfiles` bypasses the merge indexing policy

File: `src/fileio/simple_projfile_utils.f90`

Relevant lines: 174-191

The chunk merge implementation loops over source stack ranges, mutates the
source particle `stkind`, transfers the particle row, remaps output `fromp/top`,
and transfers the stack row:

```fortran
fromp = chunks(ic)%os_stk%get_fromp(i)
top   = chunks(ic)%os_stk%get_top(i)
do j = fromp,top
    iptcl_glob = iptcl_glob + 1
    call chunks(ic)%os_ptcl2D%set_stkind(j, istk)
    call merged_proj%os_ptcl2D%transfer_ori(iptcl_glob, chunks(ic)%os_ptcl2D, j)
enddo
top_glob = fromp_glob + top - fromp
call chunks(ic)%os_stk%set(i, 'fromp', fromp_glob)
call chunks(ic)%os_stk%set(i, 'top',   top_glob)
call merged_proj%os_stk%transfer_ori(istk, chunks(ic)%os_stk, i)
```

This does not call the policy-compliant `resolved_particle_indstk` logic used by
`merge_selected_project_files`. It does not treat source `indstk` as a candidate,
does not repair old-stream `indstk`, and does not explicitly verify that
`nptcls_stk` is credible.

This function is used in important stream paths, including:

- `src/main/stream/simple_microchunked2D.f90:1418`
- `src/main/stream/simple_microchunked2D.f90:1822`
- `src/main/stream/simple_stream_p05_sieve_cavgs.f90:405`
- `src/main/commanders/simple/simple_commanders_project_core.f90:389`

Recommendation:

Refactor `resolved_particle_indstk` out of `merge_selected_project_files` into a
shared project-indexing helper. Use it in `merge_chunk_projfiles` for both
`ptcl2D` and any `ptcl3D`-derived output. After remapping, assert:

```text
top - fromp + 1 == nptcls
1 <= indstk <= nptcls_stk
```

### P1: stream project writers do partial project remaps

Several stream-side writers transfer stack and particle rows, then remap only
`fromp/top` and `stkind`. They do not consistently rewrite or validate `nptcls`,
`nptcls_stk`, or `indstk`.

Primary targets:

- `src/main/project/simple_sp_project_core.f90:615`
- `src/main/stream/simple_stream_p04_refpick_extract_new.f90:416`
- `src/main/stream/simple_stream_p04_refpick_extract.f90:471`
- `src/main/stream/simple_stream_p06_pool2D_new.f90:403`
- `src/main/stream/simple_stream_p06_pool2D.f90:349`
- `src/main/stream/simple_stream_pool2D_utils.f90:314`
- `src/main/stream/simple_stream_cluster2D_microchunked.f90:205`
- `src/main/stream/simple_stream_cluster2D_subsets.f90:289`

Typical pattern:

```fortran
call dst%os_stk%transfer_ori(new_stk, src%os_stk, old_stk)
call dst%os_stk%set(new_stk, 'fromp', new_fromp)
call dst%os_stk%set(new_stk, 'top',   new_top)

call dst%os_ptcl2D%transfer_ori(new_ptcl, src%os_ptcl2D, old_ptcl)
call dst%os_ptcl2D%set_stkind(new_ptcl, new_stk)
```

This pattern is incomplete under the policy. Any code that remaps `fromp/top` and
`stkind` must also resolve `indstk` in the source stack domain and must preserve
or establish `nptcls_stk` correctly.

Recommendation:

Introduce a shared routine for copying a stack range from one project into
another. The routine should:

- copy stack metadata and stack path,
- set output `fromp/top/nptcls`,
- preserve `nptcls_stk` if present and credible,
- physically inspect the stack only when necessary and available,
- copy particle rows,
- set output `stkind`,
- resolve output `indstk` from the source row using the policy.

Stream code should use that routine instead of open-coding stack/particle range
copies.

### P1: `validate_projfile` can report a project as repaired without proving physical stack consistency

File: `src/fileio/simple_projfile_utils.f90`

Relevant lines: 1115-1133 and 1138-1176

`repair_stack_nptcls_stk` trusts an existing `nptcls_stk` if it is greater than
or equal to the project stack count:

```fortran
if( nptcls_stk >= stack_counts(istk) .and. nptcls_stk >= 0 )then
    trusted_nptcls_stk(istk) = .true.
```

This is not enough to establish physical correctness. If a stream project writes
`nptcls_stk = 50` but the file has 49 images, validation will trust the metadata
and downstream stack reads can still fail.

More seriously, `repair_particle_segment` expands `nptcls_stk` to cover fallback
`indstk`:

```fortran
if( fallback_indstk > nptcls_stk )then
    call proj%os_stk%set(stkind, 'nptcls_stk', fallback_indstk)
```

That changes a physical count without rewriting or inspecting the physical stack
file. Under the policy, `nptcls_stk` is not a repair convenience; it is the
physical image count.

Recommendation:

Validation should separate metadata repairs from physical verification.

For each stack row with a stack path:

- inspect the physical stack file when possible,
- set `nptcls_stk` from the physical file count,
- if the file cannot be inspected, report that `nptcls_stk` is unverified,
- never expand `nptcls_stk` merely to accommodate a fallback index,
- flag any particle whose resolved `indstk` exceeds verified `nptcls_stk`.

If the intended behavior is best-effort metadata normalization, the program
should say so explicitly and should not claim full conformity when physical stack
inspection was skipped or failed.

### P2: metadata-only `prune_particles` has an unsafe fallback

File: `src/main/project/simple_sp_project_ptcl.f90`

Relevant lines: 489-505

`prune_particles` preserves valid source `indstk` when `nptcls_stk` exists, but
if the source `indstk` is missing or invalid it falls back to:

```fortran
indstk = iptcl - fromp + 1
```

The fallback is policy-correct in principle, but the implementation does not
recheck that the fallback is within `1..nptcls_stk` before writing it to the
output particle rows. More importantly, the fallback is not policy-correct after
a previous metadata-only prune, because the source `fromp/top` range may already
be a project-only sparse range relative to the original physical stack.

Recommendation:

Before using the fallback, determine whether the source stack is physically
sparse:

```text
nptcls_project = top - fromp + 1
nptcls_stk > nptcls_project
```

If the stack is physically sparse, require a valid existing `indstk`. Do not
fall back to `iptcl - fromp + 1`. If the stack is not sparse, fallback is
acceptable, but the result must still satisfy:

```text
1 <= indstk <= nptcls_stk
```

If the input cannot meet these constraints, the project is not safely repairable
as a metadata-only prune. It should either hard fail or require a materializing
prune that rewrites the physical stack.

### P2: `append_project` is an older merge path that bypasses the policy

File: `src/main/project/simple_sp_project_core.f90`

Relevant lines: 263-283

`append_project` appends stack and particle rows, remaps `fromp/top`, and remaps
`stkind`. It does not normalize or validate `indstk`.

This is used by `concatenate_projects`:

- `src/main/commanders/simple/simple_commanders_project_core.f90:357`

Recommendation:

Either route `concatenate_projects` through the policy-compliant
`merge_selected_project_files`, or update `append_project` to use the shared
index resolver.

### P2: STAR stream export leaks physical/project index semantics

File: `src/main/star/simple_starproject_stream.f90`

Relevant lines: 246-248 and 299-301

`get_stkname_and_ind` returns a physical stack index in `ind_in_stk`. The export
code then uses `ind_in_stk` as a row index into `os_ptcl2D`:

```fortran
call spproj%get_stkname_and_ind('ptcl2D', i, stkname, ind_in_stk)
stkind = floor(spproj%os_ptcl2d%get(ind_in_stk, 'stkind'))
```

That should use the project particle row `i`, not the physical stack index
`ind_in_stk`.

Recommendation:

Replace those lookups with:

```fortran
stkind = spproj%os_ptcl2d%get_int(i, 'stkind')
```

This is probably not the direct cause of distributed `abinitio2D` stack-read
errors, but it is a clear violation of the same indexing separation.

### P3: current tests cover the good path, not the risky paths

File: `src/main/project/simple_project_merge_tester.f90`

Relevant lines: 21-157

The current tests verify:

- `merge_selected_project_files` preserves physical `indstk` when `nptcls_stk`
  exists,
- legacy missing `nptcls_stk` falls back to `row - fromp + 1`,
- `state = 0` rows remain present during validation.

Those tests are useful, but they do not exercise the code paths most likely to
produce stream indexing corruption:

- `merge_chunk_projfiles`,
- `projrecords2proj`,
- `simple_stream_p04_refpick_extract(_new)%write_project`,
- `simple_stream_p06_pool2D(_new)%import_sets_into_pool`,
- `simple_stream_pool2D_utils` pool subset project generation,
- `simple_stream_cluster2D_microchunked` chunk generation,
- metadata-only `prune_particles` with bad old-stream `indstk`,
- validation against a real physical stack with fewer images than metadata says.

Recommendation:

Add a small test fixture with a physical stack count mismatch:

```text
project says nptcls = 50
project says nptcls_stk = 50
physical stack file contains 49 images
one particle resolves to indstk = 50
```

`validate_projfile` should flag this as physically incompatible. It should not
write a project that still appears valid and later fails in `dstack_io%read`.

## Code That Looks Policy-Compliant

### Canonical reader

File: `src/main/project/simple_sp_project_ptcl.f90`

Relevant lines: 8-115

`map_ptcl_ind2stk_ind` centralizes the reader policy. That is good, but the
current fallback is too permissive for sparse metadata-only pruned projects. It:

- resolves `stkind`,
- reads `fromp/top`,
- uses `nptcls_stk` when present,
- checks `indstk <= nptcls_stk` when `nptcls_stk` is present,
- falls back to `iptcl - fromp + 1` when `indstk` is missing or non-positive.

That fallback should be disallowed when `nptcls_stk` shows that the project range
is physically sparse.

The downstream physical stack reader then correctly fails if the resolved
physical index exceeds the actual physical file count:

- `src/fileio/simple_discrete_stack_io.f90:48`

### Main `merge_projects`

File: `src/fileio/simple_projfile_utils.f90`

Relevant lines: 247-471 and 758-809

`merge_selected_project_files` is the best current implementation of the policy.
Its `resolved_particle_indstk` helper treats source `indstk` as a candidate,
trusts it only when `nptcls_stk` exists and the value is in range, and otherwise
uses the source-domain fallback:

```text
source_particle_project_row - source_fromp + 1
```

This behavior should be shared rather than kept local to one merge routine.

### Fresh stack imports

File: `src/main/project/simple_sp_project_stk.f90`

Relevant lines: 7-168

Fresh stack imports physically inspect the stack and write:

```text
nptcls = physical count
nptcls_stk = physical count
fromp/top = project range
indstk = 1..physical count
```

This is policy-compliant.

### Materializing `prune_project`

File: `src/main/commanders/simple/simple_commanders_project_ptcl.f90`

Relevant lines: 652-821

The materializing prune path reads through `map_ptcl_ind2stk_ind`, writes new
physical stack files, and then writes:

```text
indstk = 1..ptcl_cnt
nptcls_stk = ptcl_cnt
nptcls = ptcl_cnt
fromp/top = new project range
```

This is the correct behavior when pruning rewrites the physical stack files.

## Recommended Fix Strategy

### 1. Introduce one shared resolver

Add a shared helper with semantics equivalent to the current
`resolved_particle_indstk` in `merge_selected_project_files`.

Inputs:

```text
source os_stk
source particle segment
source particle row
```

Output:

```text
resolved physical indstk
status: preserved, fallback, invalid
```

Rules:

- require valid source `stkind`,
- require source `fromp/top`,
- require the source row to fall inside that source stack range,
- if source stack has credible `nptcls_stk`, preserve positive in-range
  `indstk`,
- otherwise derive `indstk = source_row - source_fromp + 1`,
- never derive from the output project row.

### 2. Introduce one shared stack-range copy routine

Every stream/chunk project writer should use one routine that copies a stack row
and its particle rows together. That routine should update all coupled fields:

```text
os_stk%fromp
os_stk%top
os_stk%nptcls
os_stk%nptcls_stk
os_ptcl2D/os_ptcl3D%stkind
os_ptcl2D/os_ptcl3D%indstk
```

This removes the current repeated open-coded pattern.

### 3. Make validation physically aware

`validate_projfile` should inspect physical stack files when possible. If it
cannot inspect a stack, it should report that explicitly and avoid claiming full
physical conformity.

It should never increase `nptcls_stk` unless a physical stack file was actually
rewritten or physically inspected and found to contain that many images.

### 4. Add tests for the risky paths

Minimum regression tests:

- `merge_chunk_projfiles` with source project `fromp > 1` and stale `indstk`.
- `projrecords2proj` with old-stream-style missing `nptcls_stk`.
- metadata-only `prune_particles` where fallback would exceed physical
  `nptcls_stk`.
- metadata-only `prune_particles` after a prior project-only prune where
  `nptcls_stk > nptcls` and one retained particle has missing `indstk`.
- `validate_projfile` with real physical stack count less than metadata.
- one microchunk or pool project generation path that exercises the shared copy
  routine.

## Notes On `reject_mics`

Changing the default of `reject_mics` is not itself an indexing-policy violation.
It can expose the bug by changing which micrographs or particle rows enter later
stream stages, but `state` and rejection policy should not alter the meaning of
`fromp/top/nptcls` while rows remain present.

If a path removes rows, it must rewrite project ranges. If it keeps rows with
`state = 0`, those rows remain counted by `fromp/top/nptcls`.

The indexing bug is not that particles are rejected. The bug is that some code
paths remap project rows without also resolving physical stack indices.

## Acceptance Criteria

A project emitted by merge, stream, prune, or validation should satisfy:

```text
for every stack row:
    top - fromp + 1 == nptcls
    fromp/top ranges are contiguous and cover the particle segment

for every particle row:
    stkind points to a valid stack row
    particle row lies inside os_stk(stkind)%fromp/top
    indstk is positive
    indstk <= nptcls_stk when nptcls_stk is known

for every stack row with an inspectable stack file:
    nptcls_stk == physical image count in that file
```

And for any code that rewrites project stack ranges:

```text
project row indices may change
stack row indices may change
physical image indices must not change unless the physical stack is rewritten
```
