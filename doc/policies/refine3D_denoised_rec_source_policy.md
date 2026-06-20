# Dual Particle Representative Source Policy

This document defines the intended policy for projects that carry paired raw
and denoised particle representatives. The main production use case is to keep
raw representatives as the source of particle-domain decisions while denoised
representatives provide the Cartesian reconstruction signal.

This is a design and implementation note for the V1 dual-source particle
workflow. The implemented default remains raw-only; dual-source behavior is
activated explicitly with the source-selection parameters described below.

Related policy documents:

- [refine3D_policy.md](refine3D_policy.md)
- [refine3D_auto_policy.md](refine3D_auto_policy.md)
- [refine3D_multi_policy.md](refine3D_multi_policy.md)
- [project_stack_indexing_policy.md](project_stack_indexing_policy.md)
- [sigma_calculation_policy.md](sigma_calculation_policy.md)
- [nonuniform_filtering_policy.md](nonuniform_filtering_policy.md)

## 1. Scientific Objective

The workflow supports projects with two image sources per logical particle:

1. a raw representative image
2. a denoised representative image

The production refinement mode is:

```text
alignment / assignment / shift update / sigma update = raw representatives
Cartesian reconstruction signal                       = denoised representatives
```

The denoised representatives therefore influence future alignment only through
the assembled volumes they produce at the end of an iteration. They do not
replace raw-image particle-domain decisions in the production mode.

In fixed-point form, one production iteration becomes:

```text
current volumes
  -> reprojection model
  -> raw-particle alignment / state assignment / shift update
  -> raw-particle sigma update
  -> denoised-particle Cartesian partial reconstruction
  -> volassemble and volume postprocessing
  -> next iteration volumes
```

This is not a soft-assignment EM redesign. It remains the existing SIMPLE
probabilistic pre-alignment followed by hard particle assignment workflow, with
explicit image-source selection for particle-domain and reconstruction-domain
work.

## 2. Required Project Metadata

Dual representative projects must store the denoised stack path as metadata on
the same `os_stk` row as the raw stack.

Required fields:

```text
os_stk%stk       raw stack path
os_stk%stk_den   denoised stack path
```

Example:

```text
os_stk row N:
  stk        = raw_micrograph_001_particles.mrcs
  stk_den    = my_arbitrary_denoised_reps_for_001.mrcs
  ctf        = flip
  box        = ...
  smpd       = ...
  nptcls     = project row count for this stack
  nptcls_stk = physical image count in stk
```

`stk_den` is intentionally not derived from `stk` by suffix. Arbitrary denoised
stack names must be supported.

`stk_den` is also intentionally not a separate counted stack row. Adding
denoised stacks as ordinary extra `os_stk` entries would change stack counts and
could confuse stack traversal, partitioning, project merging, pruning, and
sigma grouping.

The same particle row and physical image index identify the raw and denoised
representatives:

```text
ptcl3D row -> stkind -> os_stk row -> stk     -> raw image path
ptcl3D row -> stkind -> os_stk row -> stk_den -> denoised image path
ptcl3D row -> indstk -> physical image index inside both stacks
```

The denoised stack must have the same physical image order as the raw stack.
The same `indstk` value must identify paired raw and denoised images.

## 3. Import Contract

Import commands must be able to attach denoised representative paths to
`os_stk%stk_den`.

Recommended public import parameters:

```text
stk_den       denoised stack paired with stk
stktab_den    denoised stack table paired row-by-row with stktab
```

Single-stack import:

```text
simple_exec import_particles \
  stk=raw_particles.mrcs \
  stk_den=denoised_particles_any_name.mrcs \
  projfile=my_project.simple \
  ctf=flip ...
```

Stack-table import:

```text
simple_exec import_particles \
  stktab=raw_stacks.txt \
  stktab_den=denoised_stacks.txt \
  projfile=my_project.simple \
  ctf=flip ...
```

For `stktab_den`, row `i` must be the denoised stack corresponding to raw
`stktab` row `i`.

The existing generic `stk2` parameter should not become the durable public
contract for this workflow. It is too overloaded across other commands. The
project metadata key and command-line vocabulary should use `stk_den` so the
role is explicit.

## 4. Supported Commands

Dual-source support is required for the commands that reconstruct 3D volumes
from `ptcl3D` representatives:

- `reconstruct3D`
- `refine3D`
- `refine3D_auto`
- `refine3D_multi`

Wrapper commands must pass image-source parameters through to their
`refine3D` or `reconstruct3D` child command lines. They must not silently drop,
reset, or reinterpret them.

If a workflow does not support dual sources, it must reject the source-selection
parameters explicitly rather than silently falling back to raw stacks.

## 5. Public Source Selection

Projects with `stk_den` metadata must be able to switch between raw and
denoised representatives for debugging and testing.

The implementation should expose image-source selection by role, because
refinement has two different image consumers:

```text
match_src = raw | den
rec_src   = raw | den
```

Defaults:

```text
match_src = raw
rec_src   = raw
```

Defaulting both to raw preserves existing behavior for all projects and
commands.

The intended production dual-source mode is:

```text
match_src = raw
rec_src   = den
```

Diagnostic combinations are allowed when the required source exists:

```text
match_src = raw      rec_src = raw       current behavior
match_src = raw      rec_src = den       production dual-source mode
match_src = den      rec_src = den       denoised-only diagnostic
match_src = den      rec_src = raw       cross-source diagnostic
```

`match_src` and `rec_src` are literal selectors. No command should
silently substitute raw images for denoised images, or denoised images for raw
images, inside either phase. If `match_src=den`, alignment, state
assignment, shift updates, and matcher-owned sigma updates use the denoised
representatives. If `rec_src=den`, Cartesian reconstruction uses the
denoised representatives.

Using `match_src=den` changes particle-domain statistics and should
be documented as diagnostic/testing mode unless a later policy promotes it to a
scientific workflow.

For `reconstruct3D`, only `rec_src` is relevant because there is no
matcher pass.

## 6. V1 Scope

V1 support applies to `ptcl3D` particle representatives. The selected commands
may be high-level wrappers, but the underlying image-source contract is still
for `ptcl3D` rows and `os_stk` stack metadata.

V1 supports:

- raw/denoised `ptcl3D` representative pairs
- arbitrary denoised stack names stored in `os_stk%stk_den`
- `reconstruct3D` source switching
- `refine3D` source switching
- `refine3D_auto` source parameter pass-through
- `refine3D_multi` source parameter pass-through
- raw and denoised stacks interpreted as `ctf=flip`
- raw-image sigma estimation for the production dual-source mode

V1 does not support:

- suffix-derived denoised stack paths
- denoised stacks as counted extra `os_stk` rows
- `ptcl2D`, `cls2D`, or `cls3D` dual-source workflows
- independent denoised-image CTF metadata
- separate denoised orientation tables
- denoised sigma estimation as the production ML regularization source

## 7. CTF Convention

V1 requires raw and denoised representatives to be in the `ctf=flip`
convention.

Policy:

- raw stack metadata must report `ctf=flip`
- denoised images are interpreted as already matching the same `ctf=flip`
  convention
- `stk_den` does not provide independent CTF parameters
- CTF parameters continue to be read through the raw project metadata

If a raw stack is `ctf=yes` or `ctf=no`, denoised-source reconstruction or
denoised-source matching must hard fail in V1.

This restriction can be relaxed only after a separate policy defines how
denoised representatives and CTF correction interact.

## 8. Sigma And ML Regularization

Sigma calculation follows the matching source.

When `objfun=euclid`, the matcher must calculate residual sigma from the image
source used for particle-domain matching. In production dual-source mode that
source is raw:

```text
match_src = raw
rec_src   = den
```

When `ml_reg=yes`, reconstruction may consume grouped sigma spectra through the
existing `builder%esig` path. In production dual-source mode those spectra are
raw-derived spectra applied during denoised reconstruction.

The intended production policy is:

```text
alignment image source       = raw
assignment image source      = raw
shift-update image source    = raw
sigma-estimation source      = raw
ML sigma consumption source  = raw sigma, applied during denoised reconstruction
reconstruction image source  = denoised
```

Diagnostic source switching is allowed and intentionally literal:

```text
match_src = den  -> sigma-estimation source = denoised
match_src = raw       -> sigma-estimation source = raw
```

Users who need raw sigma estimates while reconstructing denoised volumes should
run `match_src=raw rec_src=den`. Users who select
`match_src=den` are explicitly choosing denoised particle-domain
statistics for alignment and sigma updates.

Probabilistic modes keep their existing split: probability-table commands
choose assignments, while the following matcher pass owns residual sigma
updates after assignments are applied.

## 9. Validation Rules

Validation should happen as early as the relevant source information is
available.

When any active source selector is `denoised`, fail unless all of the following
are true:

- the workflow uses `ptcl3D` particle representatives
- every active `ptcl3D` particle has a valid `stkind`
- every active `ptcl3D` particle has or can derive a valid physical `indstk`
- every referenced raw `os_stk` row has `ctf=flip`
- every referenced raw `os_stk` row has a non-empty `stk_den`
- every referenced `stk_den` path exists
- raw and denoised stack dimensions match in X and Y at import or read time
- raw and denoised physical image counts match at import time
- every active particle's `indstk` is within the denoised physical image count

For `reconstruct3D`, validation is needed when `rec_src=den`.

For `refine3D`, `refine3D_auto`, and `refine3D_multi`, validation is needed
when either `match_src=den` or `rec_src=den`.

Validation must use the same stack-indexing contract as ordinary particle
reads:

```text
stkind identifies the os_stk row
indstk identifies the physical image number
fromp/top describe project row ranges, not physical image indices
```

The validation path must not infer denoised physical indices from project row
numbers when a valid `indstk` exists.

## 10. Execution Ownership

### Commanders

Commanders own command shape, defaults, and early validation. They should not
own low-level image-read loops or assembled-volume postprocessing.

`reconstruct3D` must apply `rec_src` to its reconstruction image reads.

`refine3D` must apply `match_src` to matcher/probability-table particle
input and `rec_src` to Cartesian reconstruction particle input.

`refine3D_auto` and `refine3D_multi` must pass source selectors through to
their internal `refine3D` or `reconstruct3D` invocations.

### Strategy Layer

`simple_refine3D_strategy.f90` must preserve shared-memory and distributed
parity. New source-selection parameters must be present in worker command lines
and must resolve the same `stk_den` paths on workers as on the master.

The strategy may dispatch validation before scheduling workers. It should not
take over batch image I/O.

### Project Layer

The project layer owns paired-source resolution. A helper should map:

```text
ptcl3D particle index + source role -> stack path + physical image index
```

The helper must choose:

```text
source=raw      -> os_stk%stk
source=denoised -> os_stk%stk_den
```

This keeps `stk_den` metadata policy out of matcher and reconstruction loops.

### Matcher / Batch I/O Layer

The matcher owns the per-batch decision to read the selected matching source
for particle-domain work and the selected reconstruction source for Cartesian
partial reconstruction.

The raw default path should remain unchanged.

### Reconstruction Layer

`simple_matcher_3Drec.f90` should continue to prepare Fourier planes and grid
particles using the current `ptcl3D` orientation, state, shift, CTF metadata,
and sigma state. It should not decide how `stk_den` paths are resolved.

### Volume Assembly Layer

`volassemble` remains unchanged. It reduces partial reconstructions, restores
even/odd volumes, computes FSCs, runs automasking and nonuniform filtering, and
writes project output metadata from whatever partial reconstruction inputs the
matcher or `reconstruct3D` produced.

No additional volume-domain branch is required.

## 11. Refine3D Flow

For one matcher batch in production dual-source mode:

```text
1. sample or receive ptcl3D particle indices for update
2. read raw images through the source-aware stack resolver
3. preprocess and polarize raw images for matching
4. run the selected search strategy or consume ASSIGNMENT.dat
5. update orientation, state, shift, and sigma from the raw matcher pass
6. read paired denoised images into the reconstruction image array
7. prepare reconstruction Fourier planes from denoised images
8. grid reconstruction planes using updated ptcl3D metadata
9. write partition-local Cartesian partials
```

For diagnostic source combinations, steps 2 and 6 use the selected
`match_src` and `rec_src` respectively.

## 12. Reconstruct3D Flow

`reconstruct3D` should use `rec_src` to select the image source for the
entire reconstruction.

```text
rec_src=raw       -> read os_stk%stk
rec_src=den       -> read os_stk%stk_den
```

It should still use the same `ptcl3D` orientation, state, shift, CTF, and
even/odd metadata. Only the physical image source changes.

This is important for debugging: after a dual-representative project has been
created, users must be able to reconstruct from either raw or denoised
representatives without rewriting the project.

## 13. Wrapper Workflow Requirements

### `refine3D_auto`

`refine3D_auto` must expose or preserve the source selectors and pass them to
the child `refine3D` invocations. It must not force a raw-only child command
when the parent requested denoised reconstruction.

If `refine3D_auto` runs terminal or helper reconstruction commands, those
commands must receive the selected `rec_src`.

### `refine3D_multi`

`refine3D_multi` must expose or preserve the source selectors and pass them to
all child `refine3D` or reconstruction invocations. State count and
multi-volume policy must not alter the image-source selection.

Per-state outputs remain ordinary state volumes. Their image-source provenance
comes from the selected reconstruction source, not from a different output
filename convention.

### Other wrappers

Any other wrapper that internally invokes `refine3D` or `reconstruct3D` should
either pass through the source selectors or reject them explicitly.

## 14. Distributed Execution

Distributed and shared-memory runs must have the same scientific behavior.

For distributed workflows:

- master-side validation should check all referenced `stk_den` paths before
  scheduling workers
- worker command lines must include `match_src` and `rec_src`
- workers must resolve `stk_den` paths with the same project/current
  working-directory rules used for `stk`
- partition-local partials have the same filenames and downstream contracts as
  ordinary partials
- `volassemble` should not need to know whether partials came from raw or
  denoised images

If a worker cannot access a `stk_den` path that passed master validation, the
worker should fail hard with the raw stack path, `stk_den` path, and partition
information in the diagnostic.

## 15. Implementation Shape

The implementation prefers narrow helper APIs over spreading source-selection
conditionals through search code.

Implemented project helper:

```fortran
call spproj%get_ptcl_source_stkname_and_ind(oritype, iptcl, source, &
    stkname, ind_in_stk)
```

`source` should accept at least:

```text
raw
denoised
```

The helper should:

- accept only `oritype='ptcl3D'` for dual-source V1
- call the existing stack-index mapping to validate `stkind` and `indstk`
- choose `os_stk%stk` for `source=raw`
- choose `os_stk%stk_den` for `source=denoised`
- return the physical `indstk` unchanged
- avoid mutating `os_stk`

Implemented batch I/O helper:

```fortran
call discrete_read_imgbatch_source(params, build, source, nptcls, pinds, &
    batchlims, imgs)
```

This helper reads into either `build%imgbatch` or a caller-provided
reconstruction image array. The raw default path preserves the existing
optimized behavior.

## 16. Error Messages

Errors should name the source role and both stack paths where possible.

Examples:

```text
source=denoised requires oritype=ptcl3D
os_stk row 17 has no stk_den entry
denoised source requires raw stack ctf=flip
missing denoised stack for raw stack raw_particles.mrcs: stk_den=my_denoised_particles.mrcs
raw/denoised stack image counts differ: raw_particles.mrcs has 120000, my_denoised_particles.mrcs has 119998
raw/denoised stack dimensions differ: raw_particles.mrcs is 256x256, my_denoised_particles.mrcs is 384x384
particle indstk is outside denoised stack range
refine3D_auto did not propagate rec_src to child refine3D command
```

## 17. Testing Requirements

Minimum tests:

- import stores `stk_den` on the correct `os_stk` row
- `stktab_den` maps row-by-row to `stktab`
- validation passes for matching raw/denoised dimensions and counts
- validation fails when `stk_den` is missing
- validation fails when a `stk_den` file is missing
- validation fails when dimensions differ
- validation fails when image counts differ
- validation fails for non-`ptcl3D` dual-source use
- validation fails when raw stack CTF is not `flip`
- `rec_src=raw` preserves existing `reconstruct3D` behavior
- `rec_src=den` makes `reconstruct3D` read `stk_den`
- `refine3D` with `match_src=raw rec_src=den` reads raw
  images for matching and denoised images for reconstruction
- `refine3D` with `match_src=den rec_src=raw` reads denoised
  images for matching and raw images for reconstruction
- `refine3D` with `match_src=den rec_src=den` reads
  denoised images for both matching and reconstruction
- `refine3D_auto` propagates source selectors to child commands
- `refine3D_multi` propagates source selectors to child commands
- users can run all four raw/denoised combinations from the same project
  without rewriting stack metadata
- ML-regularized production dual-source refinement consumes raw-derived sigma
  spectra

Distributed tests should verify that worker command lines preserve
`match_src` and `rec_src`.

## 18. Follow-Up Documentation Updates

Keep the following user-facing documentation aligned with this policy:

- `refine3D_policy.md` to describe supported dual-source behavior
- `refine3D_auto_policy.md` and `refine3D_multi_policy.md` for pass-through
  requirements
- `project_stack_indexing_policy.md` to define `stk_den`
- SIMPLE UI/help text for `import_particles`, `reconstruct3D`, `refine3D`,
  `refine3D_auto`, and `refine3D_multi`
- user-facing examples that demonstrate raw and denoised stack import

This document should remain the detailed policy anchor for the V1
implementation.

## 19. Open Future Extensions

Future versions may consider:

- dual-source support for `cls3D` class-average workflows
- support for CTF conventions other than `flip`
- explicit diagnostic sigma products for `match_src=den`
- source-provenance labels in output volumes and FSC reports
- side-by-side raw and denoised reconstruction diagnostics

Those extensions should be policy-gated separately. The first implementation
should stay narrow enough that the raw particle-domain contract and existing
volume-domain contract remain easy to audit.
