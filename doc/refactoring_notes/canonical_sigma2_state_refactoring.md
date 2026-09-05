# Canonical Sigma2 State Refactoring

Date: 2026-09-04

Status: canonical opt-in implementation complete across the current 2D, 3D,
streaming, reconstruction, project-merge, conversion, and Flex consumers;
awaiting the consolidated maintainer build and runtime matrix before cutover.

Purpose: single living design, implementation plan, and validation record.

## 1. Decision Summary

Replace partition-local `sigma2_noise_partN.dat` files and
iteration-numbered `sigma2_it_N.star` files with one current
`sigma2_state.bin` per particle-project lineage.

The canonical store contains the latest native-grid residual spectrum for
each global particle row and the latest grouped model derived from those
records. Its identity does not include `nparts`, worker number, execution
mode, or iteration. Only the last committed generation is retained.

Runtime workflows neither read nor write STAR sigma files. A dedicated
`sigma2_convert` command provides explicit backward-compatible import and
export at the format boundary. It is not part of normal execution.

All new or invalid sigma state is initialized from particle power spectra.
Existing particle records survive append operations; only new rows are
initialized. Logical removal uses the project `state` flag. Physical reorder
or compaction requires an explicit row map. A native `box` or `smpd` change
invalidates the store; crop-only stage changes do not.

Distributed updates use a transactional candidate. Worker-exclusive files
are the safe default. Direct writes to disjoint candidate ranges are enabled
only on a filesystem/site configuration that has passed the production
concurrency diagnostic.

Both current scientific policies remain supported:
`sigma_est=global` and `sigma_est=group`.

## 2. Why This Model

The current in-memory distinction is valid:

- per-particle residual spectra are the authoritative update state;
- grouped even/odd spectra are derived models used by Euclidean matching and
  ML reconstruction/restoration.

The current file identities are not scientific requirements. Part files
encode scheduler ranges, and iteration STAR files encode file-selection
history. This makes restart depend on the old `nparts` and can corrupt the
meaning of fractional updates if individual residuals are reconstructed from
a group mean.

For an updated subset `U`, the correct group sum is

```text
sum(i not in U) r_i + sum(i in U) r'_i
```

Replacing unchanged records by the old group mean is generally not
equivalent. The canonical store therefore preserves every `r_i`; regrouping
and repartitioning become independent operations.

## 3. Canonical Store Contract

### Ownership and identity

`builder%esig` remains the runtime owner. A focused sigma-store domain API
owns validation and state transitions; byte-level operations belong in
`src/fileio`. Commanders and strategies orchestrate the lifecycle but never
calculate offsets or parse the format.

The project registers the canonical path explicitly. Workflows must not scan
directories for the highest-numbered sigma file.

The versioned header records:

- magic, format version, scalar kind, native shell bounds;
- particle count, native `box`, native `smpd`;
- grouping policy, group count, and section offsets;
- committed generation, provenance, state, and integrity metadata;
- an order-sensitive particle-layout identity.

The layout identity is a digest of the ordered stable particle keys owned by
the project (lineage plus normalized stack reference and stack index). The
particle `state` flag is excluded so activation changes do not invalidate the
layout. On append, the old rows are accepted only if their prefix digest
matches the committed layout. A same-size reorder therefore cannot be
mistaken for the old layout.

### Stored data and invariants

The file contains:

1. a fixed header;
2. the current native-grid grouped even/odd model;
3. one fixed-size native-grid spectrum per physical particle row;
4. integrity information for the committed generation.

Global particle row is the only record key. `nparts`, worker ranges, and
iteration are absent from the committed format. Inactive rows may remain
stored but are excluded from consolidation and ordinary consumers.

Every committed store must satisfy:

- header and file-size validation;
- layout and native-grid compatibility;
- complete, finite, positive active-particle records;
- grouped curves equal to reduction of the committed particle records under
  the recorded grouping policy;
- integrity validation for all committed sections.

The grouped section is always native-grid. Cropped consumers select the
needed shell interval through the existing grid-adaptation path.

## 4. Lifecycle Policy

### Initialization

A missing or invalid store is rebuilt with `calc_pspec` for every active
particle, followed by grouped reduction. This applies to a native-grid
mismatch and to a reorder for which no trusted map exists.

Canonical workflow bootstrap checks validate the registered header, grid,
layout, grouping, and committed state. Existence of an iteration STAR is not a
criterion in canonical mode.

`bootstrap_rec3D` currently derives a grouped curve from half-map differences
and cannot create the required per-particle state. Therefore:

- a reconstruction command with an associated particle project initializes
  canonical state through `calc_pspec`;
- a reconstruction-only command with no particle table may use the existing
  half-map estimate transiently, but it does not create or replace a canonical
  particle store.

### Updates and grouping

Workers read and write global particle ranges. A full update replaces every
active record; a fractional update replaces only selected records and copies
all others unchanged into the candidate. After the barrier, the master
reduces active particle records into `global` or `group` curves and commits
the complete generation.

### Particle-set changes

- **Append:** verify the old-layout prefix, copy old records, run `calc_pspec`
  only for the new suffix, regroup, and commit.
- **Logical removal:** keep the record but exclude `state=0` rows from
  grouping and consumption.
- **Physical compaction/reorder:** rewrite through the explicit old-to-new map
  produced by the owning project operation. Rebuild if no trusted map exists.
- **Native-grid change:** rebuild all records; persistence migration never
  interpolates shells.

## 5. Distributed Write and Recovery

The committed `sigma2_state.bin` is immutable. Every update targets
`sigma2_state.next`; only the master may create, size, validate, publish, or
discard it.

Protocol:

1. The master removes a stale candidate, creates a new exclusive candidate,
   copies the committed generation when unchanged rows must survive, records
   the scheduled global ranges, flushes, closes, and launches workers.
2. Each worker writes exactly its assigned records, flushes and syncs its
   output, closes it, and only then emits the normal completion sentinel.
3. After the distributed barrier, the master opens fresh handles, verifies
   exact non-overlapping range coverage and checksums, derives the grouped
   section, and validates the complete candidate.
4. The master flushes and syncs the candidate, atomically renames it over the
   committed file, and syncs the containing directory. A failed validation
   leaves the previous committed file untouched.

The implemented range-I/O mechanism is:

- `local` (safe default): each worker writes one exclusive temporary range
  file; the master validates and merges them. This avoids concurrent
  multi-client writes and cache interactions. These files are candidate
  material, not persistent state, so their names may contain worker/range
  information without making the committed state depend on `nparts`.
- `direct` remains a possible later optimization: workers would use GFortran
  unformatted stream `write(pos=...)` on disjoint candidate ranges. It has no
  activation path in this refactor. Adding one requires a production
  concurrency diagnostic for the exact filesystem and mount configuration.

All current workflows therefore use `local`; there is no workflow parameter
that can opt into unvalidated shared-file writes.

The implementation provides explicit low-level support for
`flush -> file fsync -> close`, atomic rename without the delete-first
behavior of `simple_rename`, and directory fsync. Filesystem probing and mount
reporting remain optional conveniences, not substitutes for a future direct-I/O
diagnostic.

## 6. STAR Compatibility Converter

`sigma2_convert` is a dedicated command and module with no dependency from
the runtime sigma-store API. It supports:

- **Exact legacy import:** assemble a canonical store from a complete set of
  legacy part files, using their recorded global ranges, and validate exact
  coverage. An accompanying grouped STAR may be checked against the derived
  model but is not authoritative.
- **STAR-only import:** expand the STAR grouped curve according to the current
  project group membership to seed per-particle records, then derive and
  validate the canonical grouped section. The command labels this conversion
  as lossy and records `provenance=star_group_seed`.
- **STAR export:** write the current grouped section to a user-selected STAR
  path for older tools. Export does not create iteration history or change
  canonical state.

Import requires an explicit project and input path; export requires an
explicit output path. The converter never runs automatically, searches for
iterations, or silently replaces a valid store. If conversion is not
requested, an old project follows the normal power-spectrum initialization
path.

## 7. Rollout and Implementation Plan

The canonical path is introduced as an explicit, typed opt-in and remains off
by default until the validation gate is complete. The selector must use the
normal `parameters` registration and propagation path; it must not be an ad
hoc worker environment flag. After all consumers pass and maintainers approve
cutover, canonical becomes the only runtime path and the temporary selector
and legacy runtime implementation are removed. The converter remains.

1. **Store and safe transactions — complete:** API, header/layout digest,
   local range merge, integrity and recovery, and explicit conversion boundary.
   Direct shared-file writes are deferred as an optional optimization.
2. **Initialization and reduction — complete:** per-particle power spectra,
   prefix-preserving append, blockwise reduction, and bootstrap identity checks.
3. **3D migration — complete for canonical opt-in:** matcher full/fractional
   updates, refine3D variants, external-reference emission, and gridding/PCG
   reconstruction.
4. **2D/restoration migration — complete for canonical opt-in:** cluster2D,
   abinitio2D, SGD/checkpoint paths, probabilistic assignment, and class-average
   restoration.
5. **Streaming/secondary migration — complete for canonical opt-in:** isolated
   chunk/pool lineages, safe dynamic-pool rebuild, project concatenation, Flex,
   cleanup, retention, and project-output behavior. Nano remains correlation-only
   and therefore has no sigma state to migrate.
6. **Cutover — pending maintainer validation:** run the consolidated matrix,
   then decide whether to change the default and remove legacy runtime I/O.
   `sigma2_convert` remains as the compatibility boundary.

### Implemented scope

The canonical opt-in now provides:

- a versioned binary store with fixed offsets, record and section integrity,
  order-sensitive particle-layout identity, deep validation, candidate copying,
  exact local-range coverage, blockwise even/odd group reduction, file and
  directory sync, and atomic publication;
- explicit project registration and a typed
  `sigma_store=legacy|canonical` selector, which remains `legacy` by default;
- shared-memory and distributed transactions for 2D and 3D matchers, including
  fractional updates, `update_missing`, probabilistic modes, and the optional
  CC residual-emission path;
- power-spectrum initialization and committed-state recovery for abinitio2D,
  abinitio2D SGD/checkpoint continuation, particle abinitio3D,
  abinitio3D_cavgs, refine3D variants, direct reconstruct3D,
  bootstrap_rec3D, and Flex PCA;
- canonical loading for class-average restoration plus gridding and PCG
  reconstruction at native and cropped working grids;
- isolated state ownership for stream chunks and changing pools. Normal append
  preserves an identity-verified prefix; an arbitrary dynamically selected pool
  without an explicit old-to-new map is rebuilt from particle power;
- exact row-wise canonical concatenation and regrouping in chunk aggregation and
  `merge_projects`. Mixed canonical/legacy general project merges discard the
  inherited registration so a later workflow rebuilds it rather than consuming
  a stale first-project path;
- an explicit `sigma2_convert` developer command for exact legacy part import,
  lossy grouped-STAR import, and grouped-STAR export, with project identity,
  native-grid, grouping, and exact-range checks;
- canonical UI opt-ins on the owning 2D, 3D, stream, reconstruction, and Flex
  programs. Normal canonical workflows do not create legacy part files or
  iteration-numbered sigma STAR files.

`simple_test_sigma2_state` covers both grouping policies, prefix identity,
candidate preparation, exact range merge, commit publication, corruption and
coverage rejection, and preservation of the prior committed generation.

The safe local-range protocol is the completed production implementation.
Direct multi-client candidate writes are deliberately deferred as an optional
performance project; they are not required for canonical semantics or cutover.

## 8. Acceptance Criteria

Storage and recovery:

- changing `nparts`, `numlen`, or shared/distributed mode reuses the same
  committed state;
- patterned disjoint writes pass repeated multi-process checks on every
  direct-enabled filesystem, including boundaries inside pages/stripes;
- missing, overlapping, truncated, or interrupted ranges reject the candidate
  and preserve the prior file byte-for-byte;
- unknown or failed filesystem validation selects the local merge path.

Scientific equivalence against the pre-refactor baseline:

- full and fractional updates preserve updated and untouched particle rows;
- `sigma_est=global` and `sigma_est=group` reproduce grouped even/odd curves;
- Euclidean scores, assignments, and sigma-weighted reconstruction/restoration
  inputs agree within the existing numerical tolerance;
- probabilistic modes retain matcher-owned residual updates; `objfun=cc` and
  `ml_reg` retain their present consumption gates.

Lifecycle and compatibility:

- append preserves the old prefix and initializes only new rows;
- inactive rows, explicit compaction maps, untrusted reorder, native-grid
  changes, and crop-only stages follow Section 4 exactly;
- an equal-size reordered project is rejected by the layout digest;
- exact legacy import, lossy STAR-only import, and STAR export are tested;
- normal workflows create no part or iteration STAR sigma files;
- every listed 2D, 3D, streaming, and secondary consumer uses the canonical
  API before cutover.

## 9. Validation Record

On 2026-09-04, source-only checks for the first implementation slice passed:

- `git diff --check`;
- repository Fortran source-index generation, including both new modules;
- command/UI registration audit, with only the three pre-existing unrelated
  name mismatches.

The maintainer subsequently reported that `simple_test_sigma2_state` passes.
The first `simple_test_units` run exposed that
`simple_abspath(..., check_exists=.false.)` does not reliably resolve a future
file on this platform. Registration was changed to construct an absolute path
from the current directory when the target does not yet exist, and the
maintainer then reported that `simple_test_units` passes. Record tested
operating systems, filesystems and mount options, scheduler/backend,
direct/local mode, numerical baseline revision, and results here before
cutover.

After adding the base-`refine3D` opt-in, the source-only checks were repeated:

- `git diff --check` passed;
- repository Fortran source-index generation completed and the generated
  module graph remained acyclic;
- the command/UI registration audit again reported only its three pre-existing
  unrelated name mismatches.

On 2026-09-05, after the base-`refine3D` plumbing and the added
update-preparation coverage, the maintainer reran `simple_test_sigma2_state`.
It completed with `SIMPLE_TEST_SIGMA2_STATE NORMAL STOP`.

Before attempting a production 3D run, the validation order was changed to
start with a laptop-scale `abinitio2D` comparison. The `refine3D` UI opt-in was
removed, and the canonical path was connected to shared-memory `abinitio2D`,
its in-memory `cluster2D` iterations, and final class-average regularization.
Source-only checks for this pivot passed: focused `git diff --check`, Fortran
index generation, and an acyclic generated module graph; the UI audit retained
only its three pre-existing unrelated mismatches. A concurrent user edit in
`doc/algorithms/abinitio2d.md` has separate trailing whitespace and was left
untouched. The maintainer subsequently reported that the canonical
shared-memory `abinitio2D` run completed gracefully and that its class averages
look excellent, providing the scientific gate for the distributed slice.

The distributed `cluster2D` path now uses the same canonical transaction as the
shared-memory path: master candidate preparation before worker dispatch,
partition-local range production, and exact-coverage master consolidation and
commit. At that intermediate stage checkpoint resume and streaming-SGD were
still explicitly gated. This
distributed slice passed focused `git diff --check`, Fortran source-index
generation with an acyclic module graph, and the UI audit with only its three
pre-existing unrelated mismatches. The maintainer subsequently reported that
the distributed test passed.

The first `abinitio3D_cavgs` exposure then passed focused `git diff --check`,
Fortran source-index generation with an acyclic module graph, and the UI audit
with the same three pre-existing unrelated mismatches. Compilation and runtime
validation were left to the maintainer, and docked multi-state canonical mode
was still explicitly gated for that first 3D slice.

The maintainer's first distributed `abinitio3D_cavgs` run reached normal
`refine3D` and map-symmetrization completion, then exposed a legacy-only sigma
file check in the standalone gridding reconstruction child. The shared
`load_sigma2_groups` boundary now accepts the project explicitly: canonical
mode resolves the registered committed state, validates native grid, ordered
layout and grouping identity, and loads its grouped curves; legacy discovery
and STAR loading remain unchanged. All gridding and PCG reconstruction callers
were updated through that common boundary. Focused `git diff --check`, Fortran
index generation with an acyclic module graph, and the unchanged UI audit
passed. The maintainer subsequently validated both shared-memory and
distributed `abinitio3D_cavgs` execution with `sigma_store=canonical`, including
the standalone reconstruction consumer that exposed the original failure.

Post-run artifact inspection found one remaining `sigma2_it_98.star`. The final
original-sampling reconstruction was rebuilding its child command line from
selected parameters and omitted `sigma_store`; its missing-STAR check therefore
entered the legacy half-map bootstrap and wrote an iteration STAR before the
regularized reconstruction. The final reconstruction now propagates both
`sigma_store` and `sigma_est`. In canonical mode it validates the registered
store against the original native grid, project layout, grouping, and committed
state. A valid store is consumed directly; an invalid native identity is rebuilt
from the associated particle project with `calc_pspec`, as required by Section
4. The legacy half-map/STAR bootstrap remains unchanged for legacy mode.

Focused `git diff --check`, Fortran source-index generation, and an acyclic
generated module graph passed after this correction. The maintainer subsequently
rebuilt and reran the canonical workflow and confirmed the clean-artifact
result with no `sigma2_it_*.star` created.

The maintainer then validated an independent multi-state canonical
`abinitio3D_cavgs` run. The remaining class-average 3D path was the docked
split. Its old pre-reconstruction `calc_group_sigmas` call exists only to
materialize the legacy iteration STAR from already current particle sigmas.
Canonical state needs no equivalent mutation because split relabelling changes
neither ordered particle identity nor global/stack membership. Canonical docked
mode now skips that legacy barrier and reuses the committed generation for the
split reconstruction; the following matcher creates and commits its normal
candidate transaction. Source-only validation passed, but shared-memory and
distributed docked runtime tests remain outstanding.

The integration pass then removed the remaining artificial workflow gates and
completed the canonical opt-in across the current runtime surface. This added
abinitio2D checkpoint recovery and SGD, refine3D sparse/update-missing and CC
emission transactions, direct reconstruct3D and bootstrap initialization,
stream chunk/pool ownership, prefix-preserving append, exact canonical state
concatenation in both chunk aggregation and `merge_projects`, the explicit
converter, and Flex initialization. Canonical selectors are exposed on the
owning 2D, 3D, stream, reconstruction, and Flex programs while the default
remains legacy until the consolidated runtime matrix passes.

After the integration pass, source-only validation completed successfully:

- focused `git diff --check` for `src/main`, `src/fileio`, and the sigma test;
- repository Fortran source-index generation;
- an acyclic generated module graph;
- command/UI registration audit with only the same three pre-existing unrelated
  name/instance mismatches.

Compilation and runtime execution were not performed by the agent, in
accordance with repository policy. The maintainer will run the consolidated
matrix in Section 10.

## 10. Consolidated Maintainer Test Matrix

1. Build the changed executables and run `simple_test_sigma2_state` plus the
   normal unit suite.
2. Run canonical abinitio2D shared and distributed, then resume a checkpoint at
   a later stage; run the developer SGD variant once.
3. Run canonical particle abinitio3D and the already established
   abinitio3D_cavgs cases, including docked multi-state if available.
4. Run canonical refine3D shared and distributed with a fractional update and
   `update_missing=yes`; cover an external-reference CC initialization path.
5. Run direct canonical reconstruct3D and bootstrap_rec3D from a project with no
   registered state and verify that power-spectrum initialization occurs.
6. Run a small canonical stream through chunk and pool updates, including one
   pool membership change and one append.
7. Exercise `sigma2_convert` exact part import, grouped STAR export/import, and
   refusal to overwrite an existing canonical output.
8. Merge two canonical chunk projects and two canonical general projects;
   verify the merged project owns a new state and runs with a different
   `nparts`. Also merge a mixed canonical/legacy pair and confirm that no stale
   state path is retained.
9. Run canonical Flex PCA once and confirm it either reuses a matching state or
   initializes one from particle power.
10. For every canonical workflow, confirm that no `sigma2_noise_part*.dat` or
    `sigma2_it_*.star` runtime artifacts are produced and that changing
    `nparts` does not change the registered committed state identity.
