# Canonical Sigma2 State Refactoring

Date: 2026-09-04

Status: proposed for maintainer review; no implementation has started.

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

`sigma2_stage_needs_bootstrap` must be changed to validate the canonical
header, grid, and layout; existence of an iteration STAR is no longer a
criterion.

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

Two range-I/O mechanisms share this protocol:

- `local` (safe default): each worker writes one exclusive temporary range
  file; the master validates and merges them. This avoids concurrent
  multi-client writes and cache interactions. These files are candidate
  material, not persistent state, so their names may contain worker/range
  information without making the committed state depend on `nparts`.
- `direct` (validated optimization): workers use GFortran unformatted stream
  `write(pos=...)` on disjoint candidate ranges. Byte positioning alone does
  not guarantee correct behavior on every distributed filesystem. Direct
  mode is allowed only after the production diagnostic has passed on that
  filesystem and mount configuration.

The master records `range_io=local|direct` in the candidate so every worker
obeys one decision. Unknown and unvalidated filesystems use `local`. A
site-level setting may enable `direct` only for a recorded validated mount;
it is not a workflow parameter and does not change scientific semantics.

The implementation needs explicit low-level support for
`flush -> file fsync -> close`, atomic rename without the delete-first
behavior of `simple_rename`, and directory fsync. File fsync and rename
bindings exist; directory-open/fsync support must be added. Filesystem probing
and mount reporting are optional conveniences, not substitutes for the
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

1. **Store and diagnostics:** implement the API, header/layout digest,
   transaction, local merge path, direct-write diagnostic, converter tests,
   and a pre-refactor numerical baseline.
2. **Initialization and reduction:** write particle spectra directly to the
   store; support append; reduce records blockwise into grouped curves; update
   bootstrap validity checks.
3. **3D migration:** migrate matcher range I/O, full/fractional updates,
   refine3D variants, reconstruction and external-reference paths. Remove the
   partition-change workaround from the canonical path.
4. **2D/restoration migration:** migrate cluster2D, abinitio2D,
   probabilistic assignment, sigma-only passes, and restoration while
   preserving existing scientific gates.
5. **Streaming/secondary migration:** give each independent chunk or pool
   lineage its own store; merge through explicit particle maps; migrate Flex,
   nano, cleanup, retention, and project-output consumers.
6. **Cutover:** complete the validation matrix, make canonical state the sole
   runtime contract, remove part/iteration discovery and runtime STAR I/O,
   update `doc/policies/sigma_calculation_policy.md`, and retain
   `sigma2_convert` as the compatibility boundary.

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

No source implementation, compilation, linking, or runtime testing has been
performed. Per repository policy, the user will compile and run the staged
validation suites. Record tested operating systems, filesystems and mount
options, scheduler/backend, direct/local mode, numerical baseline revision,
and results here before cutover.
