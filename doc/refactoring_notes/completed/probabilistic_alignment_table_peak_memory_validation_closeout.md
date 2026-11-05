# Probabilistic-alignment table: validation closeout

Status: closed and archived on 2026-11-04 by the maintainer's decision to retain
the missing historical memory/equivalence validation as explicit limits rather
than require historical reruns for archival. No peak-memory target, scientific
equivalence or runtime pass is newly established by this disposition.

Historical record, 2026-11-02: reconstructed from the [completed refactoring record](probabilistic_alignment_table_peak_memory_refactoring_completed.md),
especially its Validation section. The refactor is implemented; scientific-equivalence and
multi-state peak-memory evidence had not been reconciled. No new execution evidence is supplied here.

## Unrecorded validation retained as limits

The following were the original evidence requests, not completed measurements:

- Measure worker and master peak memory on representative single-state and multi-state cases,
  including partition and OpenMP scaling. Use a comparable baseline and retain commands and logs.
- Compare sampled particles, evaluated references and candidate order, scores, rotations, shifts,
  duplicate updates, final assignments and RNG state for fixed input, seeds and thread count.
- Cover supported dense and neighborhood modes, Euclidean/CC objectives, shifts on/off,
  populated/empty states, and shared/distributed execution. Temporary files may change bytes;
  compare decoded candidate streams and the final assignment contract.
- Exercise the 2D reuse and affected 3D refinement/solve workflows through existing test entries.

The completed record's acceptance targets are: neighborhood worker candidates
`O(R*T + Qbatch)`, global candidates `O(Q)`, dense/state candidates bounded by the batch,
thread storage linear in `T` without a hidden particle multiplier, and one published table per part.

## Historical comparison contract and limits

The old validation matrix names retired modes and assignment policies, including `prob_assign`.
Current posterior assignment has deliberately replaced that behavior. Identify which invariants
still apply and which comparisons must use an archived baseline; do not restore retired APIs or
claim byte identity across an intentional scientific-policy change. Preserve the original matrix
as evidence in the completed record, and record any accepted replacement criteria here.

The original closeout required results with baseline/source revisions and
platform, or an explicit maintainer disposition. The latter was given on
2026-11-04. The requests remain recorded here without claiming their results
or changing any current workflow's numerical acceptance criteria.
