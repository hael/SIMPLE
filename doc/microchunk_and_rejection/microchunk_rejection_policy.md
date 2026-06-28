# Microchunk Rejection and Cleanup Policy

## Scope

This document describes the current class-average rejection and particle cleanup policy in `simple_microchunked2D`.

The live stream path is implemented by:

- `src/main/stream/simple_microchunked2D.f90`
- `src/main/stream/simple_cluster2D_rejector.f90`

`model_cavgs_rejection` and the shared `src/main/cavg_quality` backend are documented separately. They are not called by the current `simple_microchunked2D` rejection path in this tree.

## Owners

`collect_and_reject` is the workflow owner. It detects finished jobs, calls rejection exactly once per eligible chunk, finalizes match chunks, updates runtime accepted/rejected counters, and marks completed lifecycle stages with sentinel files.

`reject_cavgs` owns per-chunk rejection and cleanup. It reads class averages, applies the rejection engine, writes selected/rejected class-average stacks, propagates class rejection to particle states, records selected-particle counts, and writes `REJECTION_FINISHED`.

`cluster2D_rejector` owns the scalar rejection criteria. It maintains the cumulative per-class rejection mask and exposes final rejected flags/states.

Generation routines consume rejection-complete chunks only after their project files for the next tier have been written successfully.

## Public Policy

Class-average rejection runs after a chunk has finished ab initio 2D classification.

Policy invariants:

- `state=1` means kept; `state=0` means rejected.
- rejection criteria are cumulative within one pass.
- a rejected class maps to deselected particles.
- a deselected particle has `os_ptcl2D%class=0` and `os_ptcl2D%class_match=0`.
- a chunk is not consumed by the next tier until rejection is complete.
- lifecycle flags only advance: `abinitio2D_complete`, `rejection_complete`, `complete`, or `failed`.
- when `microchunked2D%new(..., skip_pass_2=.true.)` is used, the same rejection and sentinel contract applies, but pass-2 import and generation are skipped.

## Lifecycle Sentinels

Sentinels are authoritative for restart:

- `ABINITIO2D_FINISHED`: chunk job completed.
- `REJECTION_FINISHED`: class-average rejection and particle cleanup completed.
- `REJECTION_FAILED`: rejection could not be performed safely.
- `COMPLETE`: chunk has been consumed or finalized and should not be reprocessed.

On import, `simple_microchunked2D` reconstructs chunk flags from these files:

- `failed = REJECTION_FAILED exists`
- `rejection_complete = REJECTION_FINISHED exists or failed`
- `complete = COMPLETE exists or failed`

Chunks that are neither complete nor failed have their command lines regenerated for restart.

## Rejection Eligibility

`collect_and_reject` marks a running chunk as `abinitio2D_complete` when `ABINITIO2D_FINISHED` appears, then calls `reject_cavgs`.

`reject_cavgs` returns without changing the chunk when:

- `abinitio2D_complete` is false;
- `failed` is true;
- `rejection_complete` is true.

If the class-average stack is empty, or the number of class averages does not match the number of `cls2D` rows, `reject_cavgs` writes `REJECTION_FAILED` and `COMPLETE`, sets `failed`, `rejection_complete`, and `complete`, and exits.

## Rejection Criteria

Criteria are applied in this fixed order:

1. population
2. FSC resolution
3. mask geometry
4. local variance

### Population

Reject class `i` when:

```text
pop(i) < ceiling(sum(pop) * threshold_fraction)
```

Equality with the threshold is kept.

Effective threshold fractions:

- pass 1: engine default `0.005`
- pass 2: `DEFAULT_MICRO_P2_POP_THRESH = 0.0035`
- refchunk: `DEFAULT_REF_POP_THRESH = 0.0025`
- match: `DEFAULT_REF_POP_THRESH = 0.0025`

### Resolution

Reject class `i` when:

```text
res(i) > RES_THRESHOLD
```

`RES_THRESHOLD = 40.0` Angstrom. Equality with the threshold is kept.

### Mask Geometry

For each class average, the rejector:

1. edge-normalizes the image;
2. low-pass filters to 30 Angstrom;
3. applies Otsu thresholding;
4. finds connected components;
5. removes connected components whose diameter spans the full image;
6. rejects the class if no valid component remains;
7. rejects the class if any component centroid lies outside the mask radius;
8. rejects the class if the largest component has more than `MASK_THRESHOLD` pixels outside the mask disc.

`MASK_THRESHOLD = 10.0` pixels.

### Local Variance

For each class average, the rejector:

1. edge-normalizes the image;
2. low-pass filters to 10 Angstrom;
3. applies Otsu thresholding;
4. measures local variance inside and outside the foreground mask with a window of 10 pixels.

Classes with both local-variance scores near zero are rejected unconditionally:

```text
abs(score_inside) <= ZERO_SCORE_EPS and abs(score_outside) <= ZERO_SCORE_EPS
```

`ZERO_SCORE_EPS = 1.0e-6`.

Remaining classes are robust-z-scored separately inside and outside the mask, excluding the zero-score classes. A class is rejected when one region is below the strong threshold and the other is below the weak threshold:

```text
(z_inside < strong and z_outside < weak) or
(z_inside < weak   and z_outside < strong)
```

Effective local-variance thresholds:

- pass 1: engine defaults `strong=-0.5`, `weak=-0.1`
- pass 2: `strong=-1.0`, `weak=-1.0`
- refchunk: `strong=-2.0`, `weak=-2.0`
- match: `strong=-2.0`, `weak=-2.0`

## Per-Chunk Cleanup

After rejection, `reject_cavgs`:

1. writes rejected and selected class-average stacks using `_rejected.mrc` and `_selected.mrc` suffixes;
2. writes JPEG contact sheets for non-empty selected/rejected stacks;
3. sets rejected `os_cls2D` states to zero;
4. mirrors the class states to `os_cls3D`;
5. calls `map2ptcls_state`;
6. clears `os_ptcl2D%class` and `os_ptcl2D%class_match` for particles with `state=0`;
7. writes the chunk project file;
8. writes `REJECTION_FINISHED`;
9. stores `chunk%nptcls_selected` from `os_ptcl2D%count_state_gt_zero()`;
10. sets `chunk%rejection_complete=.true.`.

For the reference chunk, `reject_cavgs` also records the full class-average stack path and box size for match-chunk generation.

When `DEBUG=.true.`, `reject_cavgs` also writes a `_deselected` project snapshot for inspection.

## Inter-Tier Consumption

In the normal full ladder, pass-1 chunks are eligible for pass-2 generation only when:

```text
rejection_complete and not complete and not failed
```

After a pass-2 project is written successfully, consumed pass-1 chunks are marked `complete` and receive `COMPLETE`.

Pass-2 chunks are eligible for reference generation under the same gate. Reference generation consumes all eligible pass-2 chunks, but those pass-2 chunks are not marked complete until match chunks are generated from them.

Match chunks are generated only after the reference stack and box are available. Each eligible pass-2 chunk is copied into one match chunk, then the source pass-2 chunk is marked `complete` and receives `COMPLETE`.

After `LAST_IMPORT_TIMEOUT`, if all pass-2 chunks are complete or failed and the reference chunk is complete, remaining rejection-complete pass-1 chunks can be merged into one final match chunk. Those pass-1 chunks are then marked `complete` and receive `COMPLETE`.

When `skip_pass_2` is active, pass-1 chunks feed the reference chunk directly. Match generation then creates one match chunk per remaining rejection-complete pass-1 chunk, and those source pass-1 chunks are marked `complete` only after the match projects are written.

## Match Finalization

For each match chunk with `rejection_complete=.true.` and `complete=.false.`, `collect_and_reject`:

1. counts selected particles as `state > 0`;
2. increments runtime accepted/rejected particle counters;
3. refreshes latest-match JPEG metadata;
4. copies the project file to `completedir`;
5. writes `COMPLETE`;
6. sets `chunk%complete=.true.`.

Runtime counters are updated when match chunks are finalized in the current process. Sentinel files and completed project copies are the durable restart record.

## Combined Output

`combine_completed_match_chunks` merges complete match chunk projects into `microchunks_match_combined.simple` or the caller-supplied combined project path.

The combined output:

- includes only complete match chunks;
- keeps particle class assignment but strips other 2D clustering parameters;
- replaces stale `cls2D`, `cls3D`, and output metadata with reference chunk class metadata;
- attaches the reference class-average stack;
- recomputes class populations from merged particle assignments;
- calls `map2ptcls_state`;
- no-ops when there are no match chunks, no complete match chunks, or the combined file already exists.

## Finished State

The workflow is finished only when:

- at least one pass-1 chunk exists and all pass-1 chunks are complete or failed;
- unless `skip_pass_2` is active, at least one pass-2 chunk exists and all pass-2 chunks are complete or failed;
- the reference chunk is complete and not failed;
- at least one match chunk exists and all match chunks are complete or failed.

## Verification

Engine-level unit coverage lives in:

- `src/main/stream/simple_cluster2D_rejector_tester.f90`

Policy-sensitive checks:

- threshold boundary behavior for population and resolution;
- cumulative rejection semantics;
- mask rejection with no valid connected component;
- local-variance zero-score rejection and robust-z-score thresholds;
- rejection failure on empty or mismatched class-average stacks;
- sentinel reconstruction on restart;
- class-to-particle state propagation;
- match finalization counters and `COMPLETE` writes.
