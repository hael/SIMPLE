# Sigma Calculation Policy

## Scope

This document defines Euclidean sigma calculation, persistence, update, and
consumption in SIMPLE. It applies to 2D/3D matching, ab initio workflows,
refine3D, reconstruction, class-average restoration, streaming pools/chunks,
and secondary consumers such as Flex PCA.

The implementation currently supports two persistence contracts selected by
`sigma_store=legacy|canonical`. `legacy` remains the default until the
canonical runtime validation matrix is complete. The scientific calculation
and consumption gates are the same in both modes.

## Scientific gates

Sigma calculation is a Euclidean-objective feature, with one explicit
wrapper-owned transition for external-reference CC pose initialization.

- `objfun=euclid` computes and updates sigma spectra.
- `objfun=cc` normally neither computes nor consumes sigma.
- Internal `cc_emit_sigma=yes` lets an external-reference CC matcher preserve
  an image-power bootstrap and replace committed particles' active shells with
  reference-conditioned residuals. This does not make CC scoring consume sigma.
- Restoration and reconstruction consume sigma only when `ml_reg=yes` and the
  effective objective is Euclidean.

The update gate is therefore

```text
cc_objfun == OBJFUN_EUCLID or cc_emit_sigma == yes
```

while the ML consumption gate remains `params%l_ml_reg`. These gates must not
be conflated.

Probabilistic refinement is not an exception. `prob_align*` consumes the
current grouped model for Euclidean scoring and owns candidate tables and
assignments, but the following matcher owns the durable residual update after
the assignment has been applied. Matcher-side sigma calculation must run in
both probabilistic and non-probabilistic modes.

## Runtime ownership

`builder%esig` (`euclid_sigma2`) owns the in-memory state:

- `sigma2_noise`: grouped spectra expanded onto particle rows for matching and
  ML reconstruction/restoration;
- `sigma2_part`: residual spectra produced for updated particles;
- `sigma2_groups`: the current even/odd grouped model.

Matcher/search code calculates per-particle residuals. Sigma commanders and
parallelization strategies own bootstrap and transaction orchestration.
`simple_sigma2_state` owns canonical validation, identity, reduction, and
state transitions; `simple_sigma2_state_file` owns byte-level canonical I/O.
Reconstruction and class restoration are consumers only.

## Grouping policy

Both current policies remain supported:

- `sigma_est=global`: one spectrum for each even/odd half-set;
- `sigma_est=group`: one spectrum for each `stkind` and half-set.

Grouping is a derived model. Per-particle residual spectra remain the
authoritative update state. Regrouping must reduce those records; it must not
reconstruct unchanged particle records from a group mean.

## Native-grid policy

Persistent spectra use the original particle grid:

```text
native shell range = 1 : fdim(params%box)-1
native identity    = params%box, params%smpd
```

`box_crop` and `smpd_crop` determine only the shells a stage consumes. A crop
change does not invalidate persistent state. A native `box` or `smpd` change
does invalidate it and requires a particle-power rebuild; persistence code
never interpolates records to migrate between native grids.

Reconstruction keeps ownership of its existing Fourier-plane grid adaptation.
Sigma persistence and stage initialization must not add a second interpolation
path.

## Canonical persistence contract

Canonical mode registers one explicit `sigma2_state.bin` path in project
`projinfo`. Runtime code does not scan for a highest iteration or infer state
from worker filenames.

The committed store contains:

1. a versioned header with native grid, particle count, grouping, generation,
   provenance, offsets, state, and integrity metadata;
2. the derived native-grid even/odd grouped model;
3. one native-grid spectrum for each physical particle row;
4. integrity data for records and sections.

Its order-sensitive layout identity is derived from project lineage plus each
ordered particle's normalized stack reference and stack index. It excludes the
particle `state` flag, so logical activation changes do not change row identity.
It also excludes `nparts`, worker number, execution mode, crop, and iteration.

Every committed state must have a valid header and file size, pass integrity
checks, match the project layout and native grid, contain finite positive active
records, and have grouped curves equal to reduction of the particle table.

## Canonical lifecycle

### Bootstrap and recovery

`calc_pspec` initializes a missing or invalid canonical state from masked
particle power. It uses one even/odd-stratified sample of at most 25,000 active
particles, applies the established average-image subtraction and positivity
conditioning, propagates the resulting half-set curves to particle records,
reduces the requested grouping, and commits the first generation.

Workflow entry points validate the registered state against the current native
grid, ordered layout, grouping, and committed state. A mismatch rebuilds from
particle power. Canonical checkpoint continuation reuses a valid committed
generation and rebuilds only when it is missing or stale.

Direct `reconstruct3D`, `bootstrap_rec3D`, and Flex PCA initialize canonical
state when they need Euclidean sigma and the associated particle project has no
valid state. The legacy half-map-difference STAR bootstrap remains confined to
legacy `bootstrap_rec3D`.

### Search update and commit

During 2D/3D matching, residual sigma is calculated after the particle
orientation, shift, class, or state assignment is committed:

- 2D uses `refkind=class`;
- 3D uses `refkind=proj`.

The residual is formed in polar Fourier space after shift, rotation, and CTF
application. Per-shell squared residual energy retains the established
normalization by `2*pftsz`.

Each canonical iteration is transactional:

1. the master creates `sigma2_state.next`, copying the committed generation
   when unchanged rows must survive;
2. each worker writes one exclusive temporary range for its assigned global
   particle rows;
3. the master verifies generation/layout identity and exact non-overlapping
   scheduled coverage, merges the ranges, and reduces active records into the
   grouped model;
4. commit-time scientific and integrity validation runs;
5. the candidate is synced and atomically published over the committed file,
   followed by directory sync.

A failed or incomplete candidate never replaces the previous committed state.
The committed format is independent of partition count and number padding.

Fractional and `update_missing` passes copy untouched records and replace only
scheduled rows. If a sparse selection contains no rows, workers still provide
the exact empty/unchanged transaction coverage required for a valid commit.

The safe implementation always uses exclusive worker range files. Direct
multi-client writes into one candidate are not exposed; they remain a possible
future optimization that requires filesystem-specific concurrency validation.

### Particle-set changes

- Append accepts old records only when the old layout digest matches the new
  layout prefix. It preserves that prefix, initializes the suffix, regroups,
  and commits.
- Logical removal retains the row but excludes `state=0` from reduction and
  ordinary consumption.
- Equal-size reorder or untrusted compaction fails identity validation and
  rebuilds unless the owning operation supplies an explicit row map.
- Canonical chunk and general project concatenation use the known row-wise
  concatenation map, preserve particle records exactly, remap stack groups,
  reduce the new grouping, and register a new output-owned state.
- A mixed canonical/legacy general project merge carries no state path into the
  output; the next Euclidean workflow rebuilds it. Stream chunk aggregation
  rejects mixed persistence contracts.

## Workflow policy

### 2D

`abinitio2D` and its SGD/checkpoint variants bootstrap at the first stage,
commit matcher residuals each iteration, and use the same canonical state for
ML class-average restoration. Shared-memory and distributed paths implement the
same candidate transaction. Canonical final output does not register or create
an iteration STAR.

Independent `abinitio2D_chunks` give each subset project its own canonical
lineage. Stream chunks likewise own isolated states rather than inheriting a
parent project's registration.

The streaming pool owns `pool/sigma2_state.bin`. Because pool membership can
change without a stable old-to-new row map, the safe policy is to rebuild the
exact pool layout with `calc_pspec` before its clustering transaction. Ordinary
append paths outside arbitrary pool selection preserve verified prefixes.

### 3D

Particle `abinitio3D`, `abinitio3D_cavgs`, refine3D, refine3D_auto,
refine3D_states, and classify3D_refs use the common canonical 3D matcher
transaction in shared and distributed execution. Full, fractional,
update-missing, probabilistic, and optional CC residual-emission passes follow
the same ownership rule.

`abinitio3D_cavgs` creates a workflow-local state for its temporary class-average
particle project and never inherits the input project's particle state. Docked
state relabeling changes only logical activity, so its reconstruction may reuse
the committed generation before the following matcher opens the next candidate.

Gridding and PCG reconstruction resolve the registered state through
`load_sigma2_groups`, validate native grid/layout/grouping, and load the grouped
model through the same boundary.

### Secondary consumers

Flex PCA validates or initializes canonical state and uses the common grouped
loader for colored-noise whitening. Nano3D explicitly runs correlation mode and
therefore does not participate in the sigma persistence lifecycle.

## Legacy compatibility and conversion

Legacy mode retains the existing runtime artifacts:

- `sigma2_noise_partN.dat` per-particle partition files;
- `sigma2_it_<iter>.star` grouped handoff files.

Legacy behavior remains available while `sigma_store` defaults to `legacy`, but
canonical workflows never read or write those runtime artifacts.

`sigma2_convert` is the only explicit format boundary:

- `parts_import` validates and imports a complete, exactly covering legacy part
  set into per-particle canonical records;
- `star_import` expands a grouped STAR over the current project membership,
  labels the state as a lossy STAR seed, reduces it again, and commits it;
- `star_export` validates project identity and writes only the current grouped
  section to a user-selected STAR file.

Conversion never runs automatically and refuses to overwrite an existing
canonical output state.

## Consumption policy

Matching consumes `sigma2_noise` whenever the objective is Euclidean.
Restoration/reconstruction consumes it only when `params%l_ml_reg` is true.
The loader expands the selected even/odd global or stack group onto particle
rows. Fourier-plane generation performs the established sigma resampling and
divides both the complex sample and CTF power contribution by sigma.

The principal implementation boundaries are:

- `src/main/simple_euclid_sigma2.f90`: in-memory model and legacy translation;
- `src/main/simple_sigma2_state.f90`: canonical identity and transactions;
- `src/fileio/simple_sigma2_state_file.f90`: canonical binary I/O;
- `src/fileio/simple_sigma2_files.f90`: dual-mode consumer boundary;
- 2D/3D matcher strategies: residual production;
- Euclid and workflow strategies: bootstrap and commit orchestration;
- reconstruction and class-average restoration: gated consumption.

## Cutover

Canonical mode is functionally complete but remains opt-in pending the runtime
matrix in `doc/refactoring_notes/canonical_sigma2_state_refactoring.md`. After
that matrix passes, maintainers may change the default and remove legacy runtime
part/STAR discovery. Both `sigma_est=global` and `sigma_est=group` remain
scientifically supported; removing grouped estimation is not part of this
refactor.
