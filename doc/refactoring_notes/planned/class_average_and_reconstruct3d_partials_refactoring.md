# Class-Average State and Reconstruct3D Partials Refactoring

Date: 2026-09-05

Revised: 2026-09-11 after completion and runtime validation of the canonical
sigma2 cutover.

Status: planning refreshed; implementation has not started. The sigma2
prerequisite is satisfied, but the pre-implementation gates in Section 9 must
be closed before production code changes begin.

Purpose: single living design, review, and validation record for removing
partition identity from durable class-average carry-over state, retiring one
dead streaming reconciliation path, and removing one invalid reconstruct3D
continuation check. The completed sigma2 work is the implementation template,
but class-average state has different scientific semantics and must not copy
the sigma2 data model mechanically.

Related contracts:

- [Canonical Sigma2 State Refactoring](../completed/canonical_sigma2_state_refactoring.md)
- [Abinitio2D Policy](../../policies/2D/abinitio2D_policy.md)
- [Class-Average Bootstrap Policy](../../policies/2D/class_average_bootstrap_policy.md)
- [Importance Sampling and Fractional Update Policy](../../policies/importance_sampling_fractional_update_policy.md)
- [Refine3D Policy](../../policies/3D/refine3D_policy.md)
- [Reconstruct3D PCG Policy](../../policies/3D/reconstruct3D_pcg_policy.md)

## 1. Decision Summary

The refreshed plan contains three independently reviewable changes, in this
order:

1. **Remove the gridding continuation file-count check.** On `continue=yes`
   with fractional update and the gridding backend, refine3D currently
   requires the previous directory to contain exactly `nparts * 4` partial
   files. Those files are transient and are deleted before every round, so the
   check is not a valid continuation invariant. Reconstruct3D partial semantics
   otherwise remain unchanged.
2. **Retire the dead streaming chunk reconciliation.** Source audit now finds
   no caller of `stream_chunk%read` and no consumer of the partless files made
   by its internal `average_into`. The division by `nparts_chunk` therefore
   does not currently under-weight pool carry-over; it is dead code, not an
   active numerical bug. Remove the unused binding and its averaging code in a
   focused cleanup, after one final repository-wide call-site check.
3. **Introduce canonical class-average carry-over state.** Workers accumulate
   current-only contributions from zero. The assembly owner reduces them in a
   deterministic order, applies the class-local recurrence exactly once, and
   publishes one registered, versioned `cavg_state.bin` per 2D lineage. The
   committed file is independent of `nparts`, worker number, execution mode,
   stage, and iteration.

The earlier proposal for two committed half files is withdrawn. Two renames
cannot form one atomic publication. Even numerator, odd numerator, even
CTF-squared sum, and odd CTF-squared sum must be sections of one committed
container and become visible through one atomic replace.

No public CLI selector will be introduced. The sigma2 rollout showed that a
temporary persistence selector propagated through nested command lines creates
a second configuration surface and can leave unexercised child paths. Rollout
will instead proceed vertically, one owning workflow at a time, with a single
carry-over path inside each migrated workflow.

## 2. Lessons Carried Forward from Canonical Sigma2

The following are requirements, not analogies to be reconsidered during
implementation:

| Sigma2 lesson | Application to class-average state |
| --- | --- |
| Separate authoritative state from derived views. | Persist unregularized Fourier numerator and CTF-squared sufficient statistics. Restored class-average stacks, FRCs, and display products remain derived outputs. |
| Durable identity must be scientific, not scheduler-shaped. | The committed identity excludes `nparts`, `part`, `numlen`, queue mode, stage, and iteration. Worker identity exists only in generation-scoped transient payloads. |
| Register the owner explicitly. | The owning project records `cavg_state`; workflows never discover state by globbing `cavgs_*_partN` or choosing a file by iteration number. |
| Keep the committed generation immutable while work is in flight. | Workers write current-only temporary payloads; only the assembly owner may create, blend, validate, publish, or discard a candidate. |
| Make the whole scientific unit atomic. | All four accumulator sections live in one container published by one atomic rename followed by directory sync. |
| Validate before publication and fail closed. | Header, file size, identity, coverage, finiteness, section integrity, and candidate generation are checked before the old committed file is replaced. |
| Define lag-one visibility at consumers. | The previous committed class state remains visible until the current update has been reduced and all outputs that consume it have succeeded. The candidate is state for the next invocation. |
| Prefer safe local worker payloads. | Workers write exclusive local payloads and the master merges them. Direct concurrent writes to a shared candidate are outside this refactor. |
| Treat append, reorder, removal, and grid changes as first-class lifecycle events. | Carry-over acceptance or rebuild behavior is explicit in Section 5; it is not inferred from matching array dimensions. |
| Cut over vertically and remove the legacy path. | The state-file API is tested first, then shared 2D, distributed 2D, and streaming are migrated. Legacy reads are removed from each workflow when that workflow crosses the boundary. |
| Audit every handoff, not only the main matcher. | Project copies, chunk/pool lineage creation, checkpoint restart, box changes, cleanup, and nested commands are part of the acceptance matrix. |

Two sigma2 design choices must *not* be copied:

- Sigma2 can preserve exact per-particle records across fractional updates.
  Class-average carry-over stores per-class aggregates and therefore cannot
  exactly remove or move an old individual contribution after a particle is
  deactivated or changes class.
- Sigma2 grouped curves are derived from authoritative particle rows. A
  class-average aggregate has no smaller exact source in the proposed store;
  it is itself the state of a stochastic recurrence. Its header identity is a
  synchronization guard, not a claim that the aggregate equals a fresh sum of
  every particle currently assigned to the class.

## 3. Current-System Review Findings

Verified against `master` on 2026-09-11.

### 3.1 Class-average carry-over is partition-shaped

`cavger_init_online` in
`src/main/class/simple_classaverager_restore.f90` reads the previous four files
for the local `part`, scales them by `1 - f(c)` from
`get_class_update_fracs`, accumulates the current worker subset, and rewrites
the same files. `cavger_assemble_sums_from_parts` then sums the worker files.

The identity

```text
sum_p (1 - f(c)) * previous_p(c) = (1 - f(c)) * previous(c)
```

holds only while the old and new worker layouts are identical. Shrinking
`nparts` silently loses prior mass from removed parts. Growing it requires
prior files that do not exist. Cleanup enumerates only the current `nparts`,
so stale files from a wider layout can survive. These are persistence defects,
not requirements of the numerical recurrence.

### 3.2 The class update fraction is aggregate policy

`get_class_update_fracs` calculates, for each active populated class, the
fraction of previously updated particles carrying the maximum current
`sampled` marker. `apply_weights2cavgs` applies one scalar `1 - f(c)` to each
of the four arrays for that class.

This is an intentionally class-local stochastic update. It does not identify
which old Fourier contributions belong to which particles, and it cannot
transport old mass exactly when class assignments or the active set change.
The canonical refactor preserves this recurrence; it does not silently
reinterpret the aggregate as exact current membership.

### 3.3 Streaming `average_into` is dead

`stream_chunk%read` sums each `cavgs_*_partN` and
`ctfsqsums_*_partN` family, divides by `nparts_chunk`, and writes a partless
MRC file. Repository-wide source search finds neither a call to this type-bound
procedure nor a reader of the four resulting partless names. Chunk memoization
passes the converged project path to the pool; pool fractional carry-over is
maintained by the pool's own partition-local files.

Consequently:

- the division does not affect the current pool recurrence;
- changing it from a mean to a sum would only modify dead output and would not
  fix the partition-shaped persistence problem; and
- the right preliminary change is removal, with a final call-site and artifact
  audit, rather than preserving it as a migration bridge.

### 3.4 Missing gridding partials are zero by design

`simple_matcher_3Drec` writes no partial for a state with no particles in a
worker part, and writes both halves for a populated state even when one half
is empty. An absent state/part pair in
`read_gridding_pair_accumulators` therefore represents a zero contribution.
The normal distributed completion barrier remains responsible for worker
success. Explicit-empty markers and a reconstruct3D manifest would broaden
this refactor without correcting an observed contract violation.

### 3.5 The continuation count check is invalid

`simple_refine3D_strategy.f90` requires `nparts * 4` gridding partials in the
previous directory on `continue=yes`, while
`remove_partial_rec_files` deletes those transient files before every round.
The check can only reject a valid continuation or pass because unrelated stale
files happen to match the count. Remove the check; do not make reconstruct3D
partials durable to satisfy it.

### 3.6 PCG already owns its payload contract

The PCG raw `(B, D)` payload is versioned, published through a temporary file,
records its geometry and provenance, and is reduced in ascending order. The
canonical class-average work neither consumes nor changes it. No PCG change is
planned.

### 3.7 Rebuild and commit hooks exist, but are not yet a lifecycle

The offline `make_cavgs` path already calls
`cavger_assemble_sums(.false.)`, which builds all four accumulators from zero.
Shared matcher restoration also starts from zero whenever fractional restore is
disabled. These are suitable numerical primitives for a canonical rebuild, but
identity rejection must be decided by a master before fractional sampling; a
worker cannot discover a bad store and repair only its own range.

The distributed `cavgassemble` owner already reduces part files, restores and
writes class products, then performs a full project write. It is the natural
candidate/commit boundary, subject to the crash-order trace in Section 9.

The streaming pool currently queues `calc_pspec` before `cluster2D` when pool
membership changes. That rebuilds canonical sigma2 for the selected pool
project; it does not rebuild class-average carry-over. Streaming therefore has
no proven full-pool class-state recovery route yet and remains the principal
design blocker.

## 4. Scientific State Contract

### 4.1 Stored meaning

For every class `c`, the committed file contains four unregularized Fourier
accumulators:

```text
previous_numerator_even(c)
previous_numerator_odd(c)
previous_ctfsq_even(c)
previous_ctfsq_odd(c)
```

They are captured before CTF-density correction, ML-prior attachment, inverse
FFT, low-resolution even/odd insertion, gridding correction, or any other
restoration mutation. Restored class averages and FRCs are never accepted as a
substitute for this state.

Workers and the shared-memory matcher accumulate the selected current
particles into zeroed arrays. After deterministic reduction, the assembly
owner computes, independently for each class and half,

```text
candidate_numerator(c) =
    (1 - f(c)) * previous_numerator(c) + current_numerator(c)

candidate_ctfsq(c) =
    (1 - f(c)) * previous_ctfsq(c) + current_ctfsq(c)
```

The current term already contains the realized sampled mass and is not divided
by `f(c)`. Edge behavior remains:

- a bootstrap or forced rebuild ignores previous state and performs a full
  current accumulation;
- `f(c) = 1` replaces the old class state;
- `0 < f(c) < 1` blends once at the owner;
- `f(c) = 0` retains the old class state exactly; and
- no previous or current support follows the existing zero-support recovery.

### 4.2 Approximation boundary

The aggregate recurrence above is the compatibility target. It is not an exact
sum over the project's current class assignments after particles move between
classes. Exact membership transport would require per-particle Fourier
numerator and CTF-squared records, with storage and I/O comparable to four
particle-image banks. That alternative is out of scope unless equivalence
testing shows that preserving the aggregate recurrence is scientifically
unacceptable.

This distinction must appear in the 2D policy documentation and in the state
API comments. Terms such as "exact class sum" or "mass transport" must not be
used for the aggregate store.

## 5. Ownership, Identity, and Lifecycle

### 5.1 Ownership and registration

The class-average/restoration domain remains the numerical owner. A focused
class-average-state API owns validation and state transitions; byte layout and
publication primitives belong in `src/fileio`. Strategies and commanders
orchestrate when to prepare, reduce, restore, and commit, but do not parse the
container or calculate offsets.

The owning project registers a `cavg_state` path using the same relocation
rules established for `sigma2_state`: a bare filename is resolved beside the
project, while an explicit path remains explicit. New chunk and pool projects
must delete inherited registration until they have deliberately created or
adopted their own lineage state.

### 5.2 Compatibility identity

The header identity must include:

- native class-state geometry: `box_crop`, `smpd_crop`, padded array bounds,
  scalar/complex kinds, and class count;
- an order-sensitive layout digest for every physical particle row, using the
  same stable key as sigma2: project lineage plus normalized stack reference
  and `indstk` selected through `stkind`;
- a checkpoint digest, over the committed layout length, of each particle's
  active flag, class assignment, even/odd membership, and whether
  `updatecnt > 0`; and
- a format/recurrence version so a future scientific-policy change cannot
  silently consume an older aggregate.

The checkpoint digest is a fail-closed synchronization token. It records the
project state from which the aggregate generation was produced; it does not
assert that the stochastic aggregate is a fresh exact sum of those assignments.
The numeric `updatecnt` and the `sampled` generation are excluded: sampling
legitimately changes them before carry-over is loaded, while the boolean
`updatecnt > 0` records whether a row can already be represented in the old
aggregate. The digest allows either side of a project/state publication
interruption to be detected on restart without requiring a cross-file atomic
transaction.

The identity must exclude `nparts`, worker number, `numlen`, queue mode, stage,
and iteration. A generation counter is recorded for transaction validation and
diagnostics, not for file discovery.

### 5.3 Lifecycle matrix

| Event | Required behavior |
| --- | --- |
| Same contributing layout and matching checkpoint | Accept the committed generation. |
| Append of particles | Accept only when the old layout and checkpoint are exact prefixes. Preserve the old aggregate; suffix rows have no prior contribution even if sampling has just changed their `updatecnt` from zero. |
| First update of appended particles | Add their current contributions and publish a candidate covering the extended layout and checkpoint. |
| Activation, deactivation, removal, reorder, compaction, or reassignment outside the normal owner-controlled update | Reject carry-over and run a full accumulation, unless the owner can prove an equivalent complete rebuild in that invocation. |
| Normal sampled reassignment inside one update | Validate the input against the pre-update project, run the aggregate recurrence once, and label the candidate with the post-update checkpoint digest. |
| Class-count or class-numbering change | Rebuild; never truncate, pad, or reinterpret old class sections automatically. |
| Crop-grid expansion with unchanged physical extent | Permit only the existing validated Fourier zero-fill/padding transform, applied once by the state owner to a candidate. Never pad each worker's old state. |
| Any other class-state geometry or representation change | Rebuild; persistence migration never interpolates arbitrary Fourier grids. |
| Missing, corrupt, stale, or unregistered state | Rebuild from the full contributing lineage; do not scan for legacy part files. |

The full-rebuild path is a real workflow requirement, not merely an error
message. It must run before fractional sampling mutates `sampled`/`updatecnt`:
the owner can use the existing zero-based `cavger_assemble_sums(.false.)`
capability to seed a valid committed state from the full project, and then run
the requested fractional iteration. Before migration, each owner must prove it
can schedule all required particles and reproduce all four accumulator
sections from zero.

For streaming, identity is computed from the full pool lineage, not from the
temporarily selected worker subproject. If the pool owner cannot make the full
lineage available for validation and rebuild, streaming cutover remains blocked.

## 6. Container and Transaction Contract

### 6.1 Committed file

There is one partless file per lineage:

```text
cavg_state.bin
```

It contains:

1. a fixed, versioned header;
2. even numerator;
3. odd numerator;
4. even CTF-squared sum;
5. odd CTF-squared sum; and
6. integrity metadata for every committed section.

The header records magic, format and recurrence versions, kinds, geometry,
array bounds, class count, section offsets and sizes, generation, both identity
digests, state/provenance flags, and integrity values. Validation checks exact
file size, header consistency, finite payload values, section integrity, and
project compatibility.

The implementation should reuse the generic POSIX publication behavior added
for sigma2 (`flush -> file sync -> close`, atomic replace without delete-first,
then directory sync). The class-average format remains a separate domain
module; it must not depend on sigma-specific headers or particle-spectrum
semantics.

### 6.2 Distributed transient payloads

Workers never open the committed state. Each writes one generation-scoped,
exclusive current-contribution payload containing all four accumulator
sections plus:

- generation and pre-update checkpoint identity;
- worker number and assigned global particle range;
- geometry, class count, and section sizes; and
- payload integrity metadata.

These files may contain `part` because they are transaction material, not
durable scientific state. The master accepts a set only when generation and
identity match and scheduled particle ranges have exact, non-overlapping
coverage. It reduces accepted payloads in ascending part number. A missing,
stale, overlapping, truncated, corrupt, or non-finite payload rejects the
candidate and preserves the committed file.

Direct concurrent writes to disjoint regions of a shared candidate have no
activation path in this refactor. They require the same filesystem-specific
production validation deferred by the sigma2 work.

### 6.3 Publication and recovery

1. Validate the registered committed input against the pre-update project. If
   invalid, select the explicit full-rebuild path.
2. Allocate zeroed current accumulators and a new internal generation. In
   distributed mode, remove only stale transaction files for this lineage and
   launch workers with generation-scoped output names.
3. Require the normal worker barrier, validate exact transient coverage, and
   reduce in ascending part order.
4. Blend the valid previous state once at the assembly owner, or omit it for a
   rebuild. Write the complete candidate container.
5. Validate the candidate deeply. Restore class averages from owned copies of
   the candidate arrays and complete FRC, class document, and project output.
   The project output already carries the registered, stable `cavg_state` path.
6. Publish the candidate over `cavg_state.bin` with one atomic replace and sync
   the directory. This rename is the final commit point.
7. Delete only the transaction files belonging to the completed generation.

The precise order of the project-file write and state rename cannot make two
files atomic. Recovery instead relies on the checkpoint digest: after a crash,
an old project with a new state or a new project with an old state fails
identity validation and enters the full-rebuild path. Publication must never
silently accept whichever file happens to exist.

The previous committed generation stays visible until Step 6. A failure after
the project checkpoint but before publication leaves the previous state
byte-identical and is detected as a digest mismatch on restart. Because state
publication is last, there is no normal path that publishes a new state and
then attempts another project checkpoint write.

## 7. Compatibility and Scope Boundaries

- Legacy `cavgs_*_partN.mrc` and `ctfsqsums_*_partN.mrc` files are never
  discovered automatically. An old or incompatible run performs a full
  accumulation.
- No legacy converter is planned. A converter would need the old `nparts`, a
  complete trustworthy part set, and the matching project checkpoint; there is
  no established live migration need. Add one only if such a need is produced
  before cutover.
- No public persistence selector, reconstruct3D manifest, explicit-empty state
  markers, or PCG generation field is part of this work.
- Iteration class-average stacks, FRCs, JPEGs, and restored outputs retain their
  current naming and retention policies. Only unregularized carry-over state is
  canonicalized.
- Selection remains owned by sampling/probabilistic workflow code.
  Class-average code consumes `sampled` and `updatecnt`; it does not choose or
  mutate the selected subset.

## 8. Rollout Plan

Each slice is separately reviewable and must preserve one production path per
migrated workflow.

1. **Independent cleanups.** Remove the invalid gridding continuation count
   check. Confirm `stream_chunk%read` and its partless products remain unused,
   then remove them. Record source-only validation separately.
2. **State-file foundation.** Add the container, identity, deep validation,
   local worker-payload, candidate, and atomic-publication APIs plus focused
   unit tests. Do not connect a production workflow yet.
3. **Shared-memory 2D.** Migrate shared `cluster2D` and `abinitio2D` through
   the common class-average assembly boundary. Establish full, fractional,
   rebuild, restart, and failure baselines before proceeding.
4. **Distributed 2D.** Convert workers to current-only payloads and the master
   to exact coverage validation, deterministic reduction, one blend, and one
   publication. Change `nparts` across a restart as the principal regression
   test.
5. **Streaming.** Give each chunk and pool an explicitly registered lineage.
   Validate against the full pool project, cover append and pool membership
   changes, and prove the rebuild route. Streaming is not enabled merely
   because batch 2D works.
6. **Cutover and cleanup.** Remove remaining persistent part-file reads,
   preservation, padding, averaging, and cleanup branches. Keep only
   generation-scoped current payloads. Update policies and generated code maps,
   then audit all runtime artifacts and child command lines.

No slice introduces `cavg_store=legacy|canonical` or an equivalent CLI flag.
The state API can be unit-tested directly, and workflow equivalence can be
tested on focused commits without making persistence policy a user parameter.

## 9. Pre-Implementation Gates

All of the following must be answered in this note before any production edit
in Section 8 starts:

1. **Aggregate policy accepted.** Confirm that compatibility with the existing
   class-local stochastic recurrence is the goal, and that per-particle exact
   class contribution storage is not required.
2. **Checkpoint digest timing verified.** The candidate fields are now defined
   in Section 5.2. Trace standard, probabilistic, and streaming sampling to
   confirm that validation can occur before class/state/EO changes, and define
   fixed text/integer serialization without depending on native derived-type
   layout.
3. **Rebuild routes traced.** `cavger_assemble_sums(.false.)` is the existing
   zero-based primitive. For shared `cluster2D`, distributed `cluster2D`,
   `abinitio2D`, a stream chunk, and the dynamic pool, identify the master that
   can invoke it on the full project before fractional sampling and publish the
   resulting seed.
4. **Streaming ownership proven.** Demonstrate that the pool owner can validate
   against its full contributing project rather than the selected subproject,
   and define behavior for removal of already-contributing pool particles.
5. **Commit boundary traced.** Identify the exact project write and class-output
   write sites around which candidate publication occurs, and verify both
   one-sided crash outcomes fail closed through the checkpoint digest.
6. **Numerical baseline recorded.** Capture representative pre-refactor
   committed part sums and restored outputs, including fractional class moves,
   so the first canonical implementation has an independent comparison target.

Items 1, 3, and 4 are scientific/workflow decisions. If any cannot be closed
without changing the recurrence, return to design review rather than beginning
an implementation that only changes file names.

## 10. Validation Matrix

### Reconstruct3D cleanup

- `continue=yes` with fractional gridding refinement succeeds when the previous
  round has removed its transient partials.
- Reconstruction output matches an equivalent non-continued run within the
  existing tolerance.
- PCG and gridding partial read/write rules are otherwise unchanged.

### Dead streaming cleanup

- Repository source and generated-index searches find no call to
  `stream_chunk%read` and no reader of its four partless output names before
  removal.
- Normal chunk memoization and pool import still use the converged project
  handoff.

### State-file unit contract

- Round-trip all four sections and header identities.
- Reject bad magic/version/kind, wrong size or offsets, non-finite data,
  geometry/class mismatch, digest mismatch, stale generation, and corrupt
  sections.
- Reject missing, overlapping, duplicate, wrong-generation, and wrong-identity
  worker payloads; accept exact coverage in ascending reduction order.
- An interrupted or rejected candidate leaves the committed file byte-identical.
- A successful commit publishes all four sections together; no mixed-half state
  can be observed.

### Scientific and workflow equivalence

- Full update and fractional `f(c) = 0`, intermediate, and `1` cases agree with
  the pre-refactor recurrence within a recorded tolerance.
- Shared and distributed results agree within the same tolerance for
  `nparts = 1` and several larger values.
- Changing `nparts`, `numlen`, or execution mode between invocations does not
  change the committed identity or multiply/drop previous state.
- Bootstrap and every rejected-identity case run a complete zero-based rebuild.
- Stage 1, zero-support, newly empty classes, ML-regularized restoration, class
  documents, and FRCs preserve current behavior.
- Append of never-updated particles preserves prior state; their first update
  adds them exactly once.
- Unsupported reorder, active-set change, class-count change, and geometry
  change reject carry-over and rebuild.
- A physical-extent-preserving crop-grid expansion migrates through the one
  owner-side padding route and agrees with the current low-frequency state;
  arbitrary geometry changes rebuild.
- An old/new project-state crash mismatch is detected in both directions.
- Independent batch, chunk, and pool lineages never inherit or share a
  `cavg_state` registration accidentally.
- Streaming covers one append, one membership change, one `nparts_chunk`
  change, and one forced rebuild.

### Cutover audit

- No worker reads committed class-average state.
- No runtime carry-over path discovers `cavgs_*_partN` or
  `ctfsqsums_*_partN` by glob or fixed worker name.
- Normal workflows leave no durable partition-shaped class-average state.
- Project moves/copies resolve the registered state by the documented path
  rules.
- Current policies and generated navigation describe the canonical-only path.

## 11. Acceptance Criteria

- The `nparts * 4` reconstruct3D continuation check is gone, with no broader
  reconstruct3D transaction change.
- The dead `stream_chunk%read` averaging path is removed and is not replaced by
  another chunk-part reconciliation file.
- `cavg_state.bin` is the only durable runtime source of unregularized
  class-average carry-over statistics.
- All four accumulator sections are validated and published atomically as one
  generation.
- Workers produce current-only transaction payloads; previous state is blended
  exactly once by the assembly owner.
- Changing partition count or shared/distributed mode neither invalidates nor
  scales carry-over.
- Project/state checkpoint mismatches fail closed and enter a proven rebuild
  path.
- The aggregate approximation boundary and lifecycle invalidation rules are
  documented in the owning 2D policies.
- Legacy persistent part-file read, preservation, padding, averaging, and
  cleanup paths are absent after cutover.

## 12. Validation Record

- 2026-09-05: initial source review confirmed partition-shaped class-average
  carry-over and the invalid gridding continuation check. Missing gridding
  state/part files were confirmed to represent zero contribution by design, so
  the proposed reconstruct3D transaction layer was dropped. No implementation,
  compilation, or runtime validation was performed.
- 2026-09-11: planning review after canonical sigma2 cutover. The maintainer
  reports the canonical-only sigma2 code runs successfully. The completed
  refactor was reviewed for storage ownership, identity, atomic publication,
  lag-one visibility, safe local worker payloads, lifecycle handling, rollout,
  and cutover lessons. The class-average proposal was corrected from a
  non-atomic file pair to one container, restored integrity/generation checks,
  added explicit project/checkpoint identity and lifecycle rules, prohibited a
  public rollout selector, and elevated rebuild/streaming ownership questions
  to pre-implementation gates.
- 2026-09-11: repository-wide source inspection found no caller of
  `stream_chunk%read` and no consumer of the partless files produced by
  `average_into`; item 2 is now classified as dead-code removal rather than a
  numerical fix. Item 1 remains present at
  `src/main/strategies/parallelization/simple_refine3D_strategy.f90` and remains
  valid. `git diff --check` passes and every relative contract link in this note
  resolves. No production code, compilation, or runtime execution is part of
  this planning change.
