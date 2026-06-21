# Dual Particle Representative Source Policy

This document defines the policy for projects that carry paired raw and
denoised particle representatives.

The current implementation uses one active particle image source at a time.
Matching/alignment, state assignment, sigma estimation, and reconstruction must
all consume the same representative images. Independent matching and
reconstruction source selectors are not supported.

Related policy documents:

- [refine3D_policy.md](refine3D_policy.md)
- [refine3D_auto_policy.md](refine3D_auto_policy.md)
- [refine3D_multi_policy.md](refine3D_multi_policy.md)
- [project_stack_indexing_policy.md](project_stack_indexing_policy.md)
- [sigma_calculation_policy.md](sigma_calculation_policy.md)
- [nonuniform_filtering_policy.md](nonuniform_filtering_policy.md)

## 1. Scientific Objective

The workflow supports projects with two stored image representations per
logical particle:

1. a raw representative image
2. a denoised representative image

At execution time, a command selects exactly one representation:

```text
ptcl_src = raw | den
```

The selected source is used consistently for:

- particle matching and alignment
- state assignment
- shift and in-plane updates
- matcher-owned sigma estimation
- Cartesian reconstruction image input

This replaces the earlier split-role model. A run must not align using raw
particles while reconstructing from denoised particles, or align using denoised
particles while reconstructing from raw particles.

## 2. Required Project Metadata

Dual representative projects store the denoised stack path as metadata on the
same `os_stk` row as the raw stack.

Required fields:

```text
os_stk%stk       raw stack path
os_stk%stk_den   denoised stack path
```

Example:

```text
os_stk row 1:
  stk        = raw_particles_001.mrcs
  stk_den    = my_arbitrary_denoised_reps_for_001.mrcs
  nptcls_stk = 10000
  ctf        = flip
```

`stk_den` is intentionally not derived from `stk` by suffix. Arbitrary denoised
stack names are valid.

`stk_den` is also intentionally not a separate counted stack row. Adding
denoised stacks as ordinary extra `os_stk` entries would change stack counts and
break existing particle-to-stack indexing assumptions.

The same particle row and physical image index identify the raw and denoised
representatives:

```text
ptcl3D row -> stkind -> os_stk row -> stk     -> raw image path
ptcl3D row -> stkind -> os_stk row -> stk_den -> denoised image path
ptcl3D row -> indstk -> physical image number in either paired stack
```

The denoised stack must have the same physical image order as the raw stack.
The same `indstk` value must identify paired raw and denoised images.

## 3. Import Contract

Import commands must be able to attach denoised representative paths to
`os_stk%stk_den`.

Supported import parameters:

```text
stk_den       denoised stack paired with stk
stktab_den    denoised stack table paired row-by-row with stktab
```

For a single stack:

```text
simple_exec prg=import_particles \
  stk=raw_particles.mrcs \
  stk_den=denoised_particles_any_name.mrcs \
  ctf=flip \
  ...
```

For multiple stacks:

```text
simple_exec prg=import_particles \
  stktab=raw_stacks.txt \
  stktab_den=denoised_stacks.txt \
  ctf=flip \
  ...
```

For `stktab_den`, row `i` must be the denoised stack corresponding to raw
`stktab` row `i`.

## 4. Public Source Selection

Projects with `stk_den` metadata must be able to switch between raw and
denoised representatives for debugging and testing.

The public selector is:

```text
ptcl_src = raw | den
```

Default:

```text
ptcl_src = raw
```

Defaulting to raw preserves existing behavior for projects that do not carry
denoised representatives.

If `ptcl_src=raw`, commands read `os_stk%stk`.

If `ptcl_src=den`, commands read `os_stk%stk_den`.

No command should silently substitute raw images for denoised images, or
denoised images for raw images. If a requested source is unavailable, the
workflow must fail early.

## 5. V1 Scope

V1 support applies to `ptcl3D` particle representatives used by
`reconstruct3D`, `refine3D`, `refine3D_auto`, and `refine3D_multi`.

`ptcl_src=den` is not supported for `ptcl2D` or other 2D applications.
2D workflows must continue to read raw particle stacks.

V1 supports:

- raw/denoised representative pairs
- arbitrary denoised stack names stored in `os_stk%stk_den`
- `reconstruct3D` source switching through `ptcl_src`
- `refine3D` source switching through `ptcl_src`
- `refine3D_auto` source parameter pass-through
- `refine3D_multi` source parameter pass-through
- raw and denoised stacks interpreted as `ctf=flip`

V1 does not support:

- suffix-derived denoised stack paths
- denoised stacks as counted extra `os_stk` rows
- independent denoised orientation tables
- separate denoised CTF metadata
- separate matching and reconstruction image sources in one run

## 6. CTF Convention

V1 requires raw and denoised representatives to be in `ctf=flip` status when
`ptcl_src=den`.

This means:

- raw `os_stk%ctf` must be `flip`
- denoised images are interpreted as already matching the same `ctf=flip`
  convention
- `stk_den` does not provide independent CTF parameters

If a raw stack is `ctf=yes` or `ctf=no`, denoised-source execution must hard
fail in V1.

This restriction can be relaxed only after a separate policy defines how
denoised representatives and CTF correction interact.

## 7. Sigma And ML Regularization

Sigma calculation follows `ptcl_src`.

When `objfun=euclid`, the matcher calculates residual sigma from the same image
source used for matching and reconstruction.

```text
ptcl_src=raw  -> raw images supply matching, sigma estimation, and reconstruction
ptcl_src=den  -> denoised images supply matching, sigma estimation, and reconstruction
```

`ml_reg=yes` may still consume grouped sigma spectra through the existing
`builder%esig` path, but those spectra must be consistent with the selected
particle source for the run.

The previous mixed-source mode, where raw sigma estimates were combined with
denoised reconstruction images inside one refinement execution, is not part of
the current policy.

## 8. Validation Rules

Validation should happen as early as the relevant source information is
available.

When `ptcl_src=den`, fail unless all of the following are true:

- the workflow uses `ptcl3D` particle representatives
- every active particle has a valid `stkind`
- every active particle has or can derive a valid physical `indstk`
- every referenced raw `os_stk` row has `ctf=flip`
- every referenced raw `os_stk` row has a non-empty `stk_den`
- every referenced `stk_den` path exists
- raw and denoised stack dimensions match in X and Y at import or read time
- raw and denoised physical image counts match at import time
- every active particle's `indstk` is within the denoised physical image count

Validation must use the same stack-indexing contract as ordinary particle
reads:

```text
stkind identifies the os_stk row
indstk identifies the physical image number
fromp/top describe project row ranges, not physical image indices
```

## 9. Execution Ownership

### Commanders

Commanders own command shape, defaults, and early validation. They should not
own low-level image-read loops or assembled-volume postprocessing.

`reconstruct3D`, `refine3D`, `refine3D_auto`, and `refine3D_multi` must expose
or preserve `ptcl_src`.

Wrapper commands must pass `ptcl_src` through to their child `refine3D` or
`reconstruct3D` invocations.

### Strategy Layer

`simple_refine3D_strategy.f90` must preserve shared-memory and distributed
parity. `ptcl_src` must be present in worker command lines and must resolve the
same `stk_den` paths on workers as on the master.

The strategy may dispatch validation before scheduling workers. It should not
take over batch image I/O.

### Project Layer

The project layer owns paired-source resolution. The source resolver maps:

```text
ptcl3D particle index + ptcl_src -> stack path + physical image index
```

The helper chooses:

```text
ptcl_src=raw -> os_stk%stk
ptcl_src=den -> os_stk%stk_den
```

This keeps `stk_den` metadata policy out of matcher and reconstruction loops.

### Matcher / Batch I/O Layer

The matcher owns per-batch particle image loading.

For online matching with Cartesian reconstruction enabled, the matcher must read
the selected particle batch from disk exactly once. Reconstruction buffers must
be copied from the already-loaded matcher batch, not filled by a second read
from either the same or a different stack source.

### Reconstruction Layer

`simple_matcher_3Drec.f90` should continue to prepare Fourier planes and grid
particles using the current `ptcl3D` orientation, state, shift, CTF metadata,
and sigma state. It should use `ptcl_src` only to select which stack path the
single reconstruction read consumes.

### Volume Assembly Layer

`volassemble` remains unchanged. It reduces partial reconstructions, restores
even/odd volumes, computes FSCs, runs automasking and nonuniform filtering, and
writes project output metadata from whatever partial reconstruction inputs the
matcher or `reconstruct3D` produced.

## 10. Refine3D Flow

For one matcher batch:

```text
1. sample or receive ptcl3D particle indices for update
2. read selected ptcl_src images through the source-aware stack resolver
3. copy the in-memory batch to reconstruction buffers if reconstruction is active
4. preprocess and polarize the selected images for matching
5. run the selected search strategy or consume ASSIGNMENT.dat
6. update orientation, state, shift, and sigma from the selected-source matcher pass
7. prepare reconstruction Fourier planes from the copied selected-source images
8. grid reconstruction planes using updated ptcl3D metadata
9. write partition-local Cartesian partials
```

There is no second particle-stack read between steps 2 and 7.

## 11. Reconstruct3D Flow

`reconstruct3D` uses `ptcl_src` to select the image source for the entire
reconstruction.

```text
ptcl_src=raw -> read os_stk%stk
ptcl_src=den -> read os_stk%stk_den
```

It still uses the same `ptcl3D` orientation, state, shift, CTF, and even/odd
metadata. Only the physical image source changes.

## 12. Distributed Execution

Distributed and shared-memory runs must have the same scientific behavior.

For distributed workflows:

- master-side validation should check all referenced `stk_den` paths before
  scheduling workers when `ptcl_src=den`
- worker command lines must include `ptcl_src`
- workers must resolve `stk_den` paths with the same project/current
  working-directory rules used for `stk`
- partition-local partials have the same filenames and downstream contracts as
  ordinary partials
- `volassemble` should not need to know whether partials came from raw or
  denoised images

If a worker cannot access a `stk_den` path that passed master validation, the
worker should fail hard with the raw stack path, `stk_den` path, and partition
information in the diagnostic.

## 13. Implementation Shape

The implementation prefers narrow helper APIs over spreading source-selection
conditionals through search code.

Implemented project helper:

```fortran
call spproj%get_ptcl_source_stkname_and_ind(oritype, iptcl, source, &
    stkname, ind_in_stk)
```

`source` accepts:

```text
raw
den
```

Implemented batch I/O helper:

```fortran
call discrete_read_imgbatch_source(params, build, source, nptcls, pinds, &
    batchlims, imgs)
```

The raw default path may use the existing optimized `discrete_read_imgbatch`
helper. Denoised source reads use `discrete_read_imgbatch_source` with
`source='den'`.

## 14. Error Messages

Errors should name the selected source and both stack paths where possible.

Examples:

```text
ptcl_src=den requires oritype=ptcl3D
os_stk row 17 has no stk_den entry
denoised source requires raw stack ctf=flip
missing denoised stack for raw stack raw_particles.mrcs: stk_den=my_denoised_particles.mrcs
raw/denoised stack image counts differ: raw_particles.mrcs has 120000, my_denoised_particles.mrcs has 119998
raw/denoised stack dimensions differ: raw_particles.mrcs is 256x256, my_denoised_particles.mrcs is 384x384
particle indstk is outside denoised stack range
refine3D_auto did not propagate ptcl_src to child refine3D command
```

## 15. Testing Requirements

Minimum tests:

- import stores `stk_den` on the correct `os_stk` row
- `stktab_den` maps row-by-row to `stktab`
- validation passes for matching raw/denoised dimensions and counts
- validation fails when `stk_den` is missing
- validation fails when a `stk_den` file is missing
- validation fails when dimensions differ
- validation fails when image counts differ
- validation fails for non-`ptcl3D` denoised-source use
- validation fails when raw stack CTF is not `flip`
- `ptcl_src=raw` preserves existing `reconstruct3D` and `refine3D` behavior
- `ptcl_src=den` makes `reconstruct3D` read `stk_den`
- `ptcl_src=den` makes `refine3D` use denoised images for matching and
  reconstruction
- `refine3D` online Cartesian reconstruction copies reconstruction buffers from
  the already-loaded matcher batch rather than reading particle stacks twice
- `refine3D_auto` propagates `ptcl_src` to child commands
- `refine3D_multi` propagates `ptcl_src` to child commands
- users can switch between `ptcl_src=raw` and `ptcl_src=den` from the same
  dual-representative project without rewriting stack metadata

Distributed tests should verify that worker command lines preserve `ptcl_src`.

## 16. Open Future Extensions

Future versions may consider:

- support for CTF conventions other than `flip`
- source-provenance labels in output volumes and FSC reports
- side-by-side raw and denoised reconstruction diagnostics

Those extensions should be policy-gated separately. The first implementation
should stay narrow enough that the single-source particle-domain and
volume-domain contracts remain easy to audit.
