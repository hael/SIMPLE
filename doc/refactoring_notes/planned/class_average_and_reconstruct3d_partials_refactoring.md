# Class-Average State and Reconstruct3D Partials Refactoring

Date: 2026-09-05 (scope reduced after review, same day)

Status: planned, reduced scope; implementation has not started and must not
start until the canonical sigma2 refactoring has landed and been validated.

Purpose: single living design, review, and validation record for removing
partition identity from durable class-average carry-over state, and for two
small fixes to reconstruct3D partial handling. A wider transaction-layer
proposal was reviewed and dropped; section 2 records why so it is not
re-proposed without new evidence.

Related contracts:

- [Canonical Sigma2 State Refactoring](canonical_sigma2_state_refactoring.md)
  (prerequisite; shares files and durability primitives)
- [Abinitio2D Policy](../policies/abinitio2D_policy.md)
- [Importance Sampling and Fractional Update Policy](../policies/importance_sampling_fractional_update_policy.md)
- [Refine3D Policy](../policies/refine3D_policy.md)
- [Reconstruct3D PCG Policy](../policies/reconstruct3D_pcg_policy.md)

## 1. Decision Summary

Three pieces of work, in this order, all after the sigma2 refactoring:

1. **Remove the gridding continuation file-count check.** On `continue=yes`
   with fractional update and the gridding backend, refine3D requires the
   previous directory to contain exactly `nparts * 4` partial files. Those
   partials are deleted before every round, so the check is a false
   invariant. One-line removal.
2. **Audit the streaming chunk part-file averaging.** `stream_chunk%read`
   sums chunk part files and divides by `nparts_chunk`, producing the mean of
   the parts. Restored class averages are unaffected (numerator and
   CTF-squared sum scale equally), but if the pool consumes these files as
   prior mass in a fractional blend, that mass is under-weighted by
   `1/nparts_chunk`. Determine whether it does. If so, fix it as a standalone
   bug ahead of item 3.
3. **Canonical class-average carry-over state, trimmed.** Replace the
   per-worker read-modify-write of `cavgs_*_partN` and `ctfsqsums_*_partN`
   with: workers accumulate the current subset from zero; the assembly owner
   blends once with the class-local realized fraction; one partless file pair
   per lineage, published by atomic rename. No source-update identity, no
   checksums, no converter, no manifest.

Not adopted: the reconstruct3D generation/manifest/completion-record
transaction layer, explicit-empty markers, checksummed worker payloads, and
idempotent commit identities. See section 2.

## 2. Review Findings

Verified against the tree on 2026-09-05.

### 2.1 Class-average carry-over is partition-shaped (confirmed)

`cavger_init_online` (`src/main/class/simple_classaverager_restore.f90`)
reads the previous files for the local `part`, scales them by `1 - f(c)` from
`get_class_update_fracs`, accumulates the current subset, and the matcher
rewrites the same four part files. The master sums parts. The identity

```text
sum_p (1 - f(c)) * prev_p(c) = (1 - f(c)) * prev(c)
```

holds only while `nparts` and the scheduled ranges are unchanged. Shrinking
`nparts` under fractional update silently drops the mass of the lost parts;
growing it makes the new parts fail to read. Cleanup enumerates only the
current `nparts`, so a wider stale layout survives. Within one `abinitio2D`
or `cluster2D` invocation `nparts` is constant, so ordinary batch runs never
trigger this. The exposure is reruns in a populated directory and the
streaming chunk/pool boundary.

### 2.2 Streaming reconciliation rescales mass (confirmed, consequence open)

`average_into` in `src/main/stream/simple_stream_chunk.f90` divides the
summed parts by `nparts_chunk`. Whether the result is ever used as fractional
prior mass in the pool is not established. This is item 2 above and is the
strongest available argument for item 3; it should be settled first.

### 2.3 Missing gridding partials are zero by design (proposal claim rejected)

`simple_matcher_3Drec` writes no partial for a state that has no particles in
a part, and writes both halves for a populated state even when one half is
empty. An absent pair in `read_gridding_pair_accumulators` therefore means
"no particles of this state in this part", and a zero contribution is
correct. `state_has_partials` is, under the same write rule, equivalent to
"some part had particles of this state". The ambiguity with a failed worker
requires a worker to emit its scheduler completion sentinel without writing
its files. The barrier makes that a bug-only path, not an operational one.
Explicit-empty markers, completion records, and generation tokens would
close a gap that the existing contract already closes in practice, at the
cost of touching the scheduler contract, the parameter object, both backends,
and volassemble.

### 2.4 Continuation file-count check is inconsistent (confirmed)

`simple_refine3D_strategy.f90` requires `nparts * 4` gridding partials in the
previous directory on `continue=yes`, and `remove_partial_rec_files` deletes
them before each round. The check can only block valid continuations or pass
by coincidence. Item 1.

### 2.5 PCG already carries identity

The PCG raw `(B, D)` payload is versioned, written through temporary file and
rename, records geometry, partition, particle count, and provenance, and is
reduced in ascending order. Extending it with a master generation was the
only PCG change proposed; without the transaction layer it has no consumer.
No PCG change is planned.

### 2.6 Dropped machinery and why

- Source-update identity and idempotent commit: an interrupted iteration is
  rerun from the previous project state, not resumed mid-commit. Nothing
  replays a commit twice today.
- Checksums on worker payloads: truncated MRC and rho files already fail
  with I/O errors on read; a checksum would catch only bit corruption, which
  is not a failure mode this workflow has seen.
- Legacy part-file converter: no live run needs migration; a rebuild is
  cheaper and safer.
- Reconstruct3D master manifest, worker completion records, quarantine
  policy: see 2.3.
- Particle-count sidecars for gridding partials were considered as a minimal
  "counts close" check and deferred; volassemble already derives state
  population from the project and the write rule in 2.3 keeps the file set
  consistent with it.

## 3. Sequencing and Prerequisites

- The canonical sigma2 refactoring touches `simple_classaverager_restore.f90`,
  the 2D and 3D matchers, and volassemble. Nothing here starts until that
  work is committed and its runtime validation is recorded, so that a
  regression in either can be bisected independently.
- Item 1 and item 2 are independent of each other and of item 3. Do them
  first, each as its own commit with its own check.
- Item 3 reuses the sigma2 atomic-publication primitive (write candidate,
  sync, rename, sync directory). No second rename implementation.

## 4. Invariants Preserved by Item 3

- Particle selection stays with the sampling and probabilistic workflow;
  class-average code consumes `sampled` and `updatecnt` and never changes
  the subset.
- Class restoration stays class-local and uses realized class update
  fractions from `get_class_update_fracs`.
- Even and odd accumulators stay independent until the existing merge point.
- Reduction order stays ascending in part number and is the only summation
  order.
- No additional particle-stack pass is introduced.
- The canonical state is captured before CTF-density correction, ML-prior
  attachment, inverse FFT, low-resolution even/odd insertion, or gridding
  correction mutate the restoration work arrays. It holds reusable
  unregularized sufficient statistics, not restored class averages.
- Iteration class-average stacks remain output artifacts and are never used
  as a substitute for accumulator state.

## 5. Trimmed Canonical Class-Average State (Item 3)

### 5.1 Files

One partless pair per 2D lineage, replacing the four `*_partN` families:

```text
cavg_state_even.bin     even numerator + even CTF-squared sum
cavg_state_odd.bin      odd  numerator + odd  CTF-squared sum
```

Exact naming and whether the two halves share one container are
implementation choices. The identity must not contain `nparts`, `part`,
`numlen`, execution mode, stage, or iteration.

The header records: magic and format version, `box_crop`, `smpd_crop`, array
bounds, class count, and a particle-layout digest (ordered stable particle
keys and even/odd membership). Class assignments are excluded from the
digest; they are expected to change and the recurrence transports mass
between classes through the current assignments. A geometry, class-count, or
layout mismatch rejects carry-over and forces the full-rebuild path.

### 5.2 Recurrence

Workers and the shared-memory path accumulate only the current selected
particles into zeroed arrays. They never read committed state. After
reduction, the assembly owner computes per class and half

```text
new_numerator(c) = (1 - f(c)) * previous_numerator(c) + current_partial_numerator(c)
new_ctfsq(c)     = (1 - f(c)) * previous_ctfsq(c)     + current_partial_ctfsq(c)
```

with `f(c)` from `get_class_update_fracs`. The current partial already
carries the sampled mass and is not divided by `f(c)`. Edge behavior is
unchanged from today:

- stage 1 and any full-rebuild path ignore previous state;
- `f(c) = 1` replaces, `0 < f(c) < 1` blends, `f(c) = 0` retains exactly;
- a class with neither previous nor current support stays explicitly empty
  and follows the existing zero-support recovery.

### 5.3 Publication

1. Reduce current contributions (ascending part order).
2. Blend into owned work arrays.
3. Write and sync the candidate pair without touching the committed pair.
4. Run restoration on a copy of the blended arrays; write class averages,
   FRCs, and project updates.
5. Atomically rename the candidate over the committed pair and sync the
   directory.

If any step fails the committed pair is unchanged and the iteration is rerun
from the previous project state, as today. No commit identity is recorded.

### 5.4 Execution paths

- Shared memory: accumulate one current update, blend, publish, restore.
- Distributed: workers write current-only part files scoped to the
  iteration; the master reduces, blends once, publishes, restores, and deletes
  the current part files. No worker preloads previous state; otherwise every
  worker would add a copy of it and the result would again depend on
  `nparts`.
- Streaming: each chunk and the pool own one lineage and one committed pair.
  The `nparts_chunk` averaging in `stream_chunk%read` is removed once the pool
  reads the committed pair directly (or after item 2 has shown it must be
  replaced by a sum, whichever comes first).

### 5.5 Rollout

1. Shared-memory `cluster2D` and `abinitio2D` first, with a laptop-scale
   `abinitio2D` run as the runtime equivalence gate.
2. Distributed `cluster2D` second, compared against the shared-memory result
   and against the pre-refactor distributed result.
3. Streaming chunk/pool last, as a separate decision with its own baseline,
   because its lineage boundaries are the least exercised.
4. Remove `cavger_readwrite_partial_sums('read')`, the part-file cleanup
   branches, and the streaming averaging only after every consumer is
   migrated.

## 6. Validation

### Item 1

- `continue=yes` with fractional gridding refinement succeeds from a
  directory whose partials were deleted by the previous round, and the
  result matches a run that was not continued.

### Item 2

- A written answer to whether the averaged chunk files enter a fractional
  blend, with the code path named. If yes, a fix and a before/after
  comparison of pool class averages at the first fractional iteration.

### Item 3

- Full update: committed pair and restored outputs agree with the
  pre-refactor accumulator within a recorded tolerance (summation order
  changes, so bit identity is not expected).
- Fractional update: agreement for `f(c) = 0`, intermediate, and `f(c) = 1`.
- Shared-memory and distributed runs agree within the same tolerance.
- `nparts = 1`, several `nparts` values, and an `nparts` change between
  invocations produce the same committed pair for the same selected records.
- Stage 1 does not consume previous state.
- Zero-support and newly empty classes follow the documented behavior.
- ML-regularized class averages and FRCs agree with baseline.
- Geometry, class-count, or layout mismatch rejects carry-over.
- A failed iteration leaves the committed pair byte-identical.
- Independent stream lineages do not share state.
- Source check: no worker reads committed state; no runtime path discovers
  carry-over by `*_partN` glob.

## 7. Acceptance Criteria

- The `nparts * 4` continuation check is gone.
- The streaming averaging question has a recorded answer and, if needed, a
  landed fix.
- The committed class-average pair is the only runtime carry-over source
  for unregularized fractional class-average statistics.
- Changing partition count or execution mode does not invalidate or
  multiply carry-over.
- Legacy part-file read, preservation, padding, averaging, and cleanup
  branches are removed after cutover.
- Reconstruct3D partial handling is otherwise unchanged.

## 8. Validation Record

- 2026-09-05: design review. Claims in 2.1, 2.2, 2.4 confirmed by source
  inspection; the "missing partial becomes zero" hazard (2.3) found to be by
  design and the transaction layer dropped. Scope reduced to items 1 to 3.
  No implementation, compilation, or runtime validation performed.
