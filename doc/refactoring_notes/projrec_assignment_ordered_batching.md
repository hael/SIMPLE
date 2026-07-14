# Projection-direction reconstruction: deferred shared-assembly option

## Status

**Deferred; do not implement this refactor without a demonstrated aggregate
memory or scheduling problem.**

`projrec=yes` currently provides numerical equivalence but no established
speedup. The present state-homogeneous implementation already removes the
`nstates` multiplication of reconstruction storage within an individual
distributed worker: each worker creates one current-state reconstructor,
writes its Cartesian partial, and releases it before the next state.

The remaining design in this note is therefore not a performance proposal. Its
single material benefit is a different placement of cubic reconstruction
storage across the distributed job:

```text
current projrec=yes:     one reconstructor per active distributed worker
deferred shared design:  no reconstructors in distributed workers;
                         one reconstructor in shared-memory assembly
```

For `nparts` simultaneously active workers, the deferred design removes up to
`nparts` concurrent worker reconstruction objects and replaces them with one
current-state reconstruction object in the shared-memory assembly phase. It
does not reduce the peak size of that one 3D reconstruction or the alignment
reprojection model.

## Current policy

Keep the current implementation unless measurements show that aggregate worker
memory, scheduler placement, or accelerator/host capacity is the limiting
resource. Do not pursue this work merely to obtain a speedup.

The deferred design's phase boundary would be:

```text
Phase A: distributed registration, assignment, and 2D partial assembly
    current alignment batches and all-state reprojection model
    finalize orientation, shift, state, projection-direction, and even/odd data
    write the worker's raw 2D numerator/CTF² partial sums

Phase B: shared-memory reconstruction
    reduce one state's at-most-nspace 2D class sums across worker parts
    compact-grid and postprocess one current-state reconstructor
```

## Deferred design, if memory pressure justifies it

Workers would still perform particle-level 2D partial-sum assembly in parallel;
that work cannot be deferred because every worker owns a different particle
subset. They would write raw numerator and CTF-squared partial class sums using
the existing distributed 2D partial-sum representation, with four dedicated
files per part (even/odd numerator and even/odd CTF²).

The slice mapping is state major:

```text
slice = (state - 1) * nspace + projection_direction
```

The shared-memory assembly phase would then process states serially:

```text
for state = 1, nstates
    reduce this state's at-most-nspace 2D class sums from all parts
    initialize one reconstructor
    compact-grid the assembled raw 2D sums
    run existing state-local restoration and postprocessing
    release the 2D sums and reconstructor
end for
```

This holds at most `nspace` 2D class sums for a state (with the usual even/odd
and numerator/CTF² companions), not `nstates * nspace` sums. For `nspace <=
5000`, that shared-memory 2D working set is considered acceptable. The worker
must write partial slices incrementally and must not allocate a full
`nstates * nspace` stack to do so.

If compact 3D gridding of the assembled state-local stack proved too slow, the
same assembled `nspace` class sums could instead be supplied to a distributed
3D reconstruction job. That is a scale-out option within the deferred design,
not an expected benefit of the current path.

## Advantages and costs

| Item | Assessment |
| --- | --- |
| Worker memory | Removes all cubic reconstructor objects from distributed workers. |
| Shared-memory assembly | Owns one 3D reconstructor at a time, plus at-most-`nspace` 2D class sums for its state. |
| File count | Four regular partial-sum files per part; no file-per-key expansion. |
| Numerical path | Retains raw 2D numerator/CTF² sums and compact 3D insertion, with no real-space conversion. |
| Runtime | No speed advantage is assumed; extra partial-sum I/O and reduction may make it slower. |
| Engineering scope | Requires complete-key grouping, incremental partial-stack I/O, and a state-sliced 2D reduction/assembly path. |

## Revisit gate

Reopen this work only if profiling shows one of the following:

1. Aggregate worker reconstructors are preventing a viable `nparts` setting.
2. Worker memory limits force a materially smaller box, state count, or job
   concurrency than the shared-memory design would permit.
3. A distributed 3D reconstruction from assembled class sums has an
   independently demonstrated runtime advantage.

Before implementation, quantify worker and shared-memory peak RSS, partial-sum
I/O volume, assembly time, and final-map equivalence against the current
validated `projrec=yes` route.
