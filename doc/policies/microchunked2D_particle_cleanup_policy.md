# Microchunked 2D Particle Cleanup Policy

## Scope

This document defines particle cleanup policy for the multi-tier streaming 2D workflow implemented by `microchunked2D`.

It covers:

- how particles are deselected during class-average rejection
- how consumed chunks are finalized and marked complete
- how accepted and rejected particle counts are tracked
- how cleaned particle metadata is combined into final outputs
- restart and idempotency behavior based on sentinel files

Primary implementation owner:

- `src/main/stream/simple_microchunked2D.f90`

Related policy:

- `doc/policies/class_average_auto_rejection_policy.md`

## Public policy

Particle cleanup is workflow-managed and happens in deterministic stages:

1. class-average rejection updates class and particle states inside each chunk project
2. consumed chunks are marked complete using sentinel files
3. finalized match chunks are copied to the completion directory
4. completed match chunks can be merged into one combined project

Cleanup is monotonic:

- a particle moved to `state=0` remains deselected for downstream accounting
- chunk lifecycle flags only advance forward (`abinitio2D_complete`, `rejection_complete`, `complete`)

## Ownership

### `reject_cavgs`

Responsibilities:

- apply class rejection and convert class decisions into particle cleanup
- map class states to particle states via `map2ptcls_state`
- set particle `class` and `class_match` to zero where `state=0`
- write chunk metadata and `REJECTION_FINISHED` sentinel
- store per-chunk selected-particle count (`nptcls_selected`)

### `collect_and_reject`

Responsibilities:

- detect completed jobs (`ABINITIO2D_FINISHED`)
- trigger `reject_cavgs` exactly once per chunk lifecycle
- finalize eligible match chunks into completion outputs
- update cumulative accepted/rejected counters for finalized match chunks
- write `COMPLETE` sentinel when chunk finalization is done

### Tier-generation routines

Responsibilities:

- consume upstream chunks only after rejection completes
- mark consumed source chunks complete and write `COMPLETE`

## Cleanup stages

### Stage 1: Intra-chunk particle cleanup

When a chunk is ab initio complete and not yet rejection-complete:

- class states are edited (`state=0` for rejected classes)
- class states are propagated to particles with `map2ptcls_state`
- particle class labels are cleared for deselected particles:
- `class=0`
- `class_match=0`

Artifacts written per chunk:

- `*_selected.mrc`
- `*_rejected.mrc`
- `REJECTION_FINISHED`

### Stage 2: Inter-tier consumption cleanup

Chunk consumption rules:

- pass-1 chunks are consumed into pass-2 only if `rejection_complete=.true.` and `complete=.false.`
- pass-2 chunks are consumed into refchunk under the same gate
- match chunks are generated from rejection-complete pass-2 chunks after refchunk rejection is complete

On successful consumption, source chunks are marked complete:

- in-memory: `complete=.true.`
- filesystem sentinel: `COMPLETE`

### Stage 3: Match finalization cleanup

For each match chunk with `rejection_complete=.true.` and `complete=.false.`:

- count selected particles (`state > 0`) and rejected particles (`state = 0`)
- accumulate totals in workflow counters:
- `n_accepted_ptcls`
- `n_rejected_ptcls`
- copy project metadata to `completedir`
- write `COMPLETE` sentinel and set `complete=.true.`

Policy note:

- cumulative accepted/rejected totals are based on finalized match chunks only
- pass-1/pass-2 selected counts are tracked separately as unconsumed inventory metrics

### Stage 4: Final combined output cleanup

`combine_completed_match_chunks` merges all complete match chunk projects into one combined project.

Cleanup/normalization expectations for combined output:

- preserve only class-assignment-relevant 2D state for downstream use
- refresh class populations from merged particle assignments
- ensure class-to-particle state consistency by remapping after merge

## Sentinel contract and restart behavior

Sentinels are authoritative for lifecycle recovery:

- `ABINITIO2D_FINISHED` means job-level compute finished
- `REJECTION_FINISHED` means cleanup/rejection pass finished
- `COMPLETE` means chunk has been consumed/finalized and should not be reprocessed

Import routines restore chunk lifecycle by probing these files. This provides restart safety and prevents double counting or duplicate cleanup.

## Counters and metrics policy

### Cumulative counters

- `get_n_accepted_ptcls()` returns cumulative accepted particles from finalized match chunks
- `get_n_rejected_ptcls()` returns cumulative rejected particles from finalized match chunks

### Inventory counters

- `get_n_pass_1_non_rejected_ptcls()` sums selected particles in rejection-complete, not-yet-consumed pass-1 chunks
- `get_n_pass_2_non_rejected_ptcls()` sums selected particles in rejection-complete, not-yet-consumed pass-2 chunks

These two counter classes must not be mixed when reporting final completion metrics.

## Invariants

- `state=0` particles must have `class=0` and `class_match=0`
- no chunk may be consumed before rejection completes
- no chunk may be finalized twice once `complete=.true.`
- restart import must reconstruct lifecycle flags from sentinels without data loss

## Failure handling policy

Cleanup/rejection should no-op safely when preconditions are not met:

- chunk not ab initio complete
- chunk already rejection-complete
- missing or empty class-average stack
- class-average count mismatch

In these cases, lifecycle flags must not be advanced.

## Governance and change control

Any change to cleanup semantics must include:

- updates to this policy
- updates to related rejection policy when criteria or propagation logic changes
- regression tests covering lifecycle gates, sentinel writes, and counter accounting
- release-note entry when user-visible acceptance/rejection totals can change
