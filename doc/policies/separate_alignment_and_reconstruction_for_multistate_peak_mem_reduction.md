# Separate alignment and reconstruction for multi-state peak-memory reduction

## Status

Active and implemented for normal CPU multi-state reconstruction in `refine3D`,
`refine3D_auto`, `refine3D_multi`, `abinitio3D`, and shared reconstruction
entry points. The auto and multi commanders execute their refinement stages
through the common `refine3D` path.

## Policy

Alignment and partial reconstruction are separate phases. Alignment retains the
complete all-state reprojection model; reconstruction begins only after final
assignments are frozen and alignment-only memory is released.

Distributed workers reconstruct one `(state, even_or_odd)` group at a time.
They use one standalone half-map `reconstructor`, write its existing Cartesian
partial, and destroy it before the next group. Workers must not construct the
EO composite `build%eorecvol` during normal CPU reconstruction.

`volassemble` owns `build%eorecvol`: it is the first stage that requires both
half maps concurrently for reduction, FSC, restoration, and postprocessing.

```text
Phase A — alignment
    retain the all-state reprojection model
    determine final orientation, shift, state, and even/odd labels

assignment barrier
    release PFTC/reprojections, matching images, search workspaces, and batches

Phase B — worker partial reconstruction
    group selected particles by (state, even/odd)
    for each populated state
        for even, then odd
            construct one standalone reconstructor
            reread and grid only this state/half particle group
            write the existing half-map partial and rho artifact
            destroy the reconstructor
        end for
    end for
```

An empty half of a populated state still writes an explicit zero partial. This
preserves the paired even/odd partial-file contract consumed by `volassemble`.

## Ownership and lifetime

| Phase | May remain allocated | Must not remain allocated |
| --- | --- | --- |
| Alignment | PFTC/reprojection model, matching images, search and probability workspaces | Partial reconstruction volumes and reconstruction image buffers |
| Barrier | Final orientation/state/half metadata, selected particle indices, CTF/project metadata, symmetry, Euclidean sigma data | PFTC/reprojections, matching images, search strategies, particle PFTs, correlation caches, alignment batches |
| Worker reconstruction | One half-map reconstructor, bounded reconstruction image/Fourier buffers, one `(state, half)` index range | EO composite, another state/half reconstructor, alignment objects |
| Assembly | EO composite and state-local postprocessing objects | Worker alignment state and arrays of state volumes |

The worker memory target is therefore:

```text
max(alignment working set, one state/half reconstruction working set)
```

It must be independent of `nstates` in cubic reconstruction storage.

## Required invariants

1. Reconstruction uses the exact selected particle subset from the completed
   alignment phase.
2. Orientation, shift, state, and even/odd labels are read-only after the
   assignment barrier.
3. Every valid selected particle belongs to exactly one `(state, half)` group.
4. A worker has at most one initialized standalone half-map reconstructor.
5. Each populated state produces both expected half-map partial artifacts.
6. Existing Cartesian partial names, payloads, CTF handling, interpolation,
   symmetry, and downstream assembly behavior are unchanged.
7. Euclidean sigma data may remain available through reconstruction; retaining
   it must not retain PFTC or reprojection allocations.

## Scope limits

This policy does not change alignment search order, objectives, or reprojection
model residency. It does not change the Cartesian worker-partial contract or
the `volassemble` reduction and postprocessing path. Any alternative artifact
or reconstruction-ownership model requires a separate design decision.

The specialized OpenMP-offload kernel may use an EO composite for `nstates=1`.
For `nstates>1`, it must dispatch and return before that object is constructed,
using the normal state/half CPU path.

## Verification

Validate one and multiple states, unbalanced populations, empty states and
empty halves, both particle sources, ML regularization, symmetry, distributed
execution, and normal downstream assembly/postprocessing.

Check particle counts by state and half, Cartesian partial payloads, restored
half maps, FSC/resolution, final map quality, reconstruction runtime, and peak
RSS. Confirm at the assignment barrier that PFTC and alignment image/search
objects are absent before the first half-map reconstructor is created.

Implementation anchors: `simple_strategy3D_matcher`, `simple_matcher_3Drec`,
`simple_reconstructor_openmpoffload`, and `commander_volassemble`.
