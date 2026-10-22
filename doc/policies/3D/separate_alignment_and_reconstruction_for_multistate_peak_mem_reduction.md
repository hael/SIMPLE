# Separate alignment and reconstruction for multi-state peak-memory reduction

## Status

Active and implemented for normal CPU multi-state reconstruction in `refine3D`,
`refine3D_auto`, `refine3D_states`, `solve3D`, and shared reconstruction
entry points. The commanders of `refine3D_auto`, `refine3D_states`,
`refine3D_states_cavgs` and `classify3D_refs` execute their refinement stages
through the common `refine3D` path.

## Policy

Alignment and partial reconstruction are separate phases. Alignment retains the
complete all-state reprojection model; reconstruction begins only after final
assignments are frozen and alignment-only memory is released.

Distributed workers reconstruct one `(state, even_or_odd)` group at a time.
They use one standalone half-map `reconstructor`, write its existing Cartesian
partial, and destroy it before the next group. Workers must not hold both half
maps at once during normal CPU reconstruction.

`volassemble_exec` is the first stage that holds both half maps concurrently
(its even, odd, read and sum accumulators) for reduction, FSC, restoration, and
postprocessing.

```text
Phase A — alignment
    retain the all-state reprojection model
    determine final orientation, shift, state, and even/odd labels

assignment barrier
    release PFTC/reprojections, matching images, search workspaces, and batches

Phase B — worker partial reconstruction
    group selected particles by (state, even/odd): hard labels, or weighted membership
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

Under a state weight set (`m_estimator=flex`: a frozen weight per particle and
state, published by `flex_pca`), the grouping follows the weights instead of
the labels (`group_pinds_by_weights` in `simple_matcher_3Drec`). A particle
joins the group of every state whose weight for it is above zero, in the half
of its even/odd label, and enters that state's reconstructor with its weight.
With `state=` only that state is grouped. A particle is therefore read once per
state it weighs into. The weight scales both the data term and the density term
of the inserted plane; partial names and payload formats are unchanged. The PCG
(preconditioned conjugate gradients) worker groups the same way, with its
membership threshold of 0.01 (`PCG_WEIGHT_THRESHOLD`).

## Ownership and lifetime

| Phase | May remain allocated | Must not remain allocated |
| --- | --- | --- |
| Alignment | PFTC/reprojection model, matching images, search and probability workspaces | Partial reconstruction volumes and reconstruction image buffers |
| Barrier | Final orientation/state/half metadata, selected particle indices, CTF/project metadata, symmetry, Euclidean sigma data | PFTC/reprojections, matching images, search strategies, particle PFTs, correlation caches, alignment batches |
| Worker reconstruction | One half-map reconstructor, bounded reconstruction image/Fourier buffers, one `(state, half)` index range (with its weights under a weight set) | Another state/half reconstructor, alignment objects |
| Assembly | Both half-map accumulators and state-local postprocessing objects | Worker alignment state and arrays of state volumes |
| All three worker phases (alignment, barrier, reconstruction) | The particle buffer: the decoded raw images of the pass's share, one immutable slab filled before the search and read by copy into the bounded batches of the search and the reconstruction, released after the reconstruction. A cropped reconstruction instead keeps one record per particle, the untapered cropped observation the search's batches make from their own transforms, sealed after the search and read by copy into the reconstruction's batches | Anything else derived from those images (prepared, filtered or transformed particles) |

The worker memory target is therefore:

```text
max(alignment working set, one state/half reconstruction working set) + particle buffer
```

The particle buffer is allocated only within its budget: half the task memory a queue system grants
(`job_memory_per_task` of the computing environment), otherwise the node's physical memory divided by twice the number of processes it runs (the parts and their
master, or one process in shared memory). A share above the budget, or an allocation that fails, leaves the
buffer empty, and the search and the reconstruction read their batches from the stacks as before. The
budget, the buffer's size and the outcome are in the search's benchmark file. Plan:
[particle_read_path_policy.md](../particle_read_path_policy.md); the design is in
`doc/implementation_notes/completed/particle_io_few_particle_stacks.md`, section 3.1.

It must be independent of `nstates` in cubic reconstruction storage.

## Required invariants

1. Reconstruction uses the exact selected particle subset from the completed
   alignment phase.
2. Orientation, shift, state, and even/odd labels are read-only after the
   assignment barrier. State weights are read-only throughout the run.
3. Under hard labels, every valid selected particle belongs to exactly one
   `(state, half)` group. Under a state weight set, a valid selected particle
   belongs to one group per state it weighs into: it may enter several states,
   once per state, and in each it is in the half of its even/odd label. A
   particle without a weight above the backend's threshold belongs to no group.
4. A worker has at most one initialized standalone half-map reconstructor.
5. Each populated state produces both expected half-map partial artifacts.
   Under a weight set, a state is populated when it has weighted members.
6. Existing Cartesian partial names, payloads, CTF handling, interpolation,
   symmetry, and downstream assembly behavior are unchanged.
7. Euclidean sigma data may remain available through reconstruction; retaining
   it must not retain PFTC or reprojection allocations.

## Scope limits

This policy does not change alignment search order, objectives, or reprojection
model residency. It does not change the Cartesian worker-partial contract or
the `volassemble` reduction and postprocessing path. Any alternative artifact
or reconstruction-ownership model requires a separate design decision.

The OpenMP-offload kernel (`calc_3Drec_gpu`, built only with
`USE_OPENMP_OFFLOAD`; no caller routes to it) keeps one expanded even/odd
accumulator pair for `nstates=1`. For `nstates>1`, it dispatches to the normal
state/half CPU path before allocating that pair.

## Verification

Validate one and multiple states, unbalanced populations, empty states and
empty halves, both particle sources, ML regularization, symmetry, distributed
execution, hard labels and fractional state weights, and normal downstream
assembly/postprocessing.

Check particle counts (under weights, memberships) by state and half,
Cartesian partial payloads, restored half maps, FSC/resolution, final map
quality, reconstruction runtime, and peak RSS. Confirm at the assignment
barrier that PFTC and alignment image/search objects are absent before the
first half-map reconstructor is created.

Implementation anchors: `simple_strategy3D_matcher`, `simple_matcher_3Drec`,
`simple_state_weight_set`, `simple_reconstructor_openmpoffload`, and volume
assembly (`volassemble_exec` in `simple_volassemble_strategy.f90`, which the
distributed `refine3D` master runs on its own builder and the
`commander_volassemble` wrapper runs for the command-line callers).
