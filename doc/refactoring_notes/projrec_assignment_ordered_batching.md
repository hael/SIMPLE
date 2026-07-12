# Projection-direction reconstruction: assignment-ordered batching

## Summary

`projrec=yes` should use two distinct batching policies with a hard phase
boundary between them:

1. **Alignment batching** exists to perform 3D registration, shift refinement,
   and state assignment efficiently. It may retain the current particle order
   and current matcher batching policy.
2. **Reconstruction batching** is created only after the alignment results are
   final. It reorders the selected particles by their assigned state,
   even/odd half, and discrete projection direction so that each complete raw
   2D Fourier numerator/CTF-squared sum can be written once and released.

The reconstruction path must not inherit the alignment batches. Alignment
batches are arbitrary slices of the selected particle list and therefore mix
states and projection directions. Accumulating compact reconstruction data in
that order either keeps many 2D sums alive or spills duplicate partial records
for the same projection direction. Both outcomes weaken the intended memory
reduction.

This separation must be visible in the top-level matcher as two sequential
loops. Compact reconstruction must not be performed opportunistically inside
the registration loop.

The fundamental reconstruction unit is a completed key:

```text
(state, even_or_odd, projection_direction)
```

A reconstruction batch is a memory-bounded collection of complete keys, not an
arbitrary particle range.

## Scope decision

The immediate memory refactor is deliberately limited to reconstruction
ownership. It does **not** state-window, shard, or otherwise redesign the
all-state reprojection model used by registration. Keeping only one state's
reprojections resident would require decomposing the existing particle-centric
search procedures into prepare, state-window search, candidate merge, shift
refinement, and final-assignment phases. That is a substantially larger search
refactor and is not justified as the first step.

The simpler high-value boundary is:

```text
distributed workers: registration + compact 2D sufficient statistics only
volassemble:          all 3D insertion, restoration, and postprocessing
```

The current `calc_projdir3Drec` implementation does not yet satisfy this
boundary. It calls `init_rec` with volume initialization enabled,
`preprecvols` constructs one `reconstructor_eo` for every populated state, and
the worker inserts its compact 2D sums into `build%eorecvols(state)` before
writing legacy Cartesian partials. This validates the compact interpolation
path, but worker memory still contains `O(nstates * box^3)` reconstruction
storage.

`volassemble`, by contrast, already has the correct outer lifetime: it loops
over states serially and uses singleton current-state reconstruction objects.
The compact-input branch must preserve that property. Its single-state peak
may include several cubic restoration and nonuniform-filter work images, but
those lifetimes must never be multiplied by `nstates`.

## Goals

- Keep the current registration and state-assignment machinery unchanged.
- Freeze the final orientation, shift, state, projection-direction, sampling,
  and even/odd assignments before reconstruction starts.
- Visit every selected particle exactly once during compact reconstruction.
- Assemble each populated projection-direction sum completely before writing
  it.
- Write at most one compact record for each populated key and distributed
  part.
- Release or reuse a completed key's 2D storage immediately after writing.
- Make worker reconstruction memory independent of `nstates` and bounded by a
  small number of 2D accumulators.
- Remove all expanded Cartesian reconstruction volumes from distributed
  `projrec=yes` workers; allocating one worker volume per state is explicitly
  forbidden.
- Preserve the existing raw numerator/CTF-squared convention and the current
  2D and 3D Kaiser-Bessel interpolation kernels.
- Keep `projrec=no` behavior unchanged.

## Non-goals

- This refactor does not change the alignment search order or objective.
- It does not change reprojection-model residency or refactor the 3D search
  procedures into state-windowed searches.
- It does not prioritize removal of the comparatively small
  `nrefs * nthreads` search workspaces; the reprojection bank remains the
  dominant alignment allocation.
- It does not change particle sampling or rebalance particles between
  distributed parts after alignment.
- It does not restore or normalize projection-direction averages in 2D.
- It does not transform compact sums through real space.
- It does not require all `nspace * nstates` planes to be resident at once.
- It does not initially redesign volume postprocessing or nonuniform filtering.

## Current Problem

The integrated `refine3D` matcher currently creates alignment batches from the
selected particle list. These batches are appropriate for:

- particle image reading and preparation
- polar reference matching
- orientation and shift search
- probabilistic or greedy state assignment
- updating the in-memory project orientation field

They are not appropriate for compact reconstruction. A single alignment batch
can contain particles from many states, projection directions, and even/odd
halves. Conversely, particles belonging to one projection-direction sum can be
spread across many alignment batches.

The current `calc_projdir3Drec` path partially compensates by constructing a
state list and scanning projection directions inside each image batch. It still
initializes the normal multi-state reconstruction volumes and ultimately
converts the compact sums into Cartesian partial reconstructions on the worker.
It therefore validates the numerical path but does not yet establish the
desired memory ownership.

## Required Phase Boundary

The matcher must have an explicit transition:

```text
Phase A: registration and assignment
    existing alignment batches
    update orientation/shift/state
    finalize discrete projection labels
    finalize selected particle set and even/odd labels
                         |
                         | assignment freeze / barrier
                         v
Phase B: compact reconstruction
    build assignment-ordered reconstruction index
    process complete projection-direction keys
    append each completed 2D sum once
    release accumulator storage
```

No orientation, shift, state, projection-direction, sampling, or even/odd value
used by reconstruction may change after the barrier. If a later workflow step
needs to alter assignments, it must occur before reconstruction indexing or
force that index to be rebuilt.

## Required Top-level Loop Structure

The integrated matcher should be organized around two explicit loops.

### Loop 1: 3D registration and state assignment

```text
for alignment_batch in current_alignment_batches
    read/build particle images for matching
    run 3D orientation and shift search
    run state assignment or consume probabilistic assignment
    update the project orientation field
end for

finalize assignment metadata
derive final discrete projection-direction labels
write/commit orientation results as required by the workflow
freeze assignments
```

This loop owns all mutation of particle orientation, shift, state, sampling,
and projection-direction metadata. It must not allocate compact reconstruction
accumulators or insert particles into partial reconstructions when
`projrec=yes`.

### Loop 2: assignment-ordered 2D partial-average assembly

```text
build grouped index from the frozen selected particles
open one temporary compact container for this distributed part

for reconstruction_batch in complete_key_batches
    for key in reconstruction_batch
        reset one accumulator for (state, eo, proj)
        for io_chunk in particles_belonging_to_key
            reread particle images for reconstruction
            generate reconstruction Fourier planes
            accumulate the residual in-plane rotations in 2D
        end for
        append the completed numerator/CTF^2 record once
        release or reuse the accumulator
    end for
end for

finalize and atomically publish the compact part container
```

The second loop is read-only with respect to assignment metadata. It is
expected to reread particle images because its order is deliberately different
from the alignment loop. Attempting to retain the alignment-loop particle
images until reconstruction would trade the removed volume memory for a large
particle-image cache and defeat the memory policy.

The two loops must use the same selected particle set. The grouped index is a
permutation of that set after invalid state-zero entries are handled according
to existing reconstruction policy; it is not a new sampling operation.

For distributed execution, the barrier is local to each worker after its
assigned particle subset has completed registration. It does not require
redistributing particles between workers. The same key may occur in several
distributed parts; each part contributes one complete local record for that
key, and `volassemble` combines the part contributions linearly.

## Reconstruction Key and Ordering Policy

The canonical key order should be:

```text
state -> even/odd -> projection_direction
```

The key must include even/odd because the numerator and CTF-squared sums remain
separate until half-map reconstruction and FSC calculation.

Within one key, particle indices may be ordered for storage locality without
changing the numerical meaning. A useful secondary order is:

```text
source stack/file -> particle image index
```

This secondary order can reduce random reads when particles assigned to the
same projection direction come from multiple stacks. It must never split the
logical ownership of the key: all chunks for that key feed the same live 2D
accumulator until the key is complete.

## Building the Reconstruction Index

A comparison sort is unnecessary. The key domain is already discrete and
bounded by `nstates * 2 * nspace`, so a counting-sort/prefix-sum layout is the
preferred implementation.

For the exact selected `pinds` used by the registration phase:

1. Ensure the final closest projection-direction label is available from the
   final orientation field.
2. Count particles for each `(state, eo, proj)` key.
3. Prefix-sum the counts to produce key offsets.
4. Scatter particle indices into one grouped particle-index array.
5. Materialize a compact list of populated keys and their `[first,last]`
   ranges.

Suggested logical data shape:

```text
type projdir_group_index
    keys(:)       ! populated (state, eo, proj) tuples
    offsets(:)    ! start/end boundaries in grouped_pinds
    grouped_pinds(:)
end type
```

The dense counting workspace contains integers only. Its size is
`nstates * 2 * nspace`, which is acceptable compared with Fourier image data.
The actual 2D accumulators must be allocated only for populated keys selected
for the current reconstruction batch.

The index builder must reject or skip, according to existing policy:

- state zero
- state outside `1:nstates`
- missing or invalid projection-direction labels
- invalid even/odd labels
- particles not present in the selected reconstruction subset

It must report counts so the following invariant can be checked:

```text
number indexed == number of valid selected reconstruction particles
```

## Reconstruction Batch Definition

A reconstruction batch is a set of complete keys whose accumulators fit within
an explicit memory budget.

For each key in the batch:

1. Allocate or reset one native-grid 2D Fourier numerator/CTF-squared
   accumulator.
2. Read that key's particles in ordinary I/O chunks. The I/O chunk size is
   independent of the reconstruction batch size.
3. Generate padded particle Fourier planes using the existing reconstruction
   preparation path.
4. Apply the residual in-plane rotation and accumulate with the current 2D
   Kaiser-Bessel machinery.
5. After all particles in the key have been consumed, append one compact record
   to the part container.
6. Immediately reset or release the accumulator before moving to another key.

This distinction is important:

```text
reconstruction batch = complete projection-direction keys held concurrently
I/O chunk             = particle images read together while completing a key
```

A key with more particles than the image I/O limit is processed in multiple
I/O chunks, but it is still written once.

## Memory Policy and Parallelism

The lowest-memory implementation uses one live key accumulator at a time:

```text
peak compact memory ~= one 2D numerator + one 2D CTF-squared plane
```

This should be the correctness-first implementation. It is deterministic and
makes the memory objective easy to verify.

Key-level parallelism can be added later by processing a small number of keys
concurrently. Each lane must own its accumulator and particle Fourier-plane
workspace. The number of lanes must be derived from a reconstruction-memory
budget, not automatically equated to `nthr`.

Conceptually:

```text
nlanes = min(populated_keys_remaining,
             configured_or_derived_memory_budget / bytes_per_key_accumulator)
```

OpenMP workers may cooperate within a key or own different keys, but no two
threads may update the same accumulator without an explicit reduction policy.
The current class-average ownership model, where one thread owns a target
slice, is the preferred pattern.

The worker memory target becomes:

```text
O(number_of_selected_particles) integer index
+ O(nstates * nspace) integer counts/offsets
+ O(reconstruction_key_lanes * box^2) Fourier sums
+ bounded particle image/Fourier I/O buffers
```

There must be no `O(nstates * box^3)` reconstruction allocation on the
`projrec=yes` worker path.

## Compact Output Policy

Use one compact container per distributed part rather than files per key or per
alignment batch:

```text
PROJDIR_PART001.bin
PROJDIR_PART002.bin
...
```

Each completed key is appended once. The container should include:

- magic and format version
- box, sampling, `nspace`, `nstates`, symmetry, and interpolation convention
- distributed part identifier
- a record directory or footer index
- for each record:
  - state
  - even/odd half
  - projection-direction index
  - particle population
  - raw Fourier numerator half-plane
  - raw CTF-squared half-plane

The on-disk representation should use canonical logical half-plane ordering,
not expose FFTW physical addressing. This keeps the format independent of the
in-memory class-average stack layout.

The worker should write to a temporary filename and atomically rename it only
after the directory/footer and all records are complete. The queue completion
sentinel must not be declared before the compact container is closed.

## `volassemble` Consumption Policy

`volassemble` already processes states serially. The compact branch should keep
that policy:

1. Initialize or reset one current-state expanded `reconstructor_eo`.
2. For every distributed part container, seek to records for the current state.
3. Stream one record at a time.
4. Recover the projection orientation from the shared eulspace index.
5. Insert the raw numerator and CTF-squared plane through
   `grid_plane_compact` into the current even or odd reconstruction.
6. After all part records for the state have been inserted, compress once.
7. Continue through the existing density correction, FSC, trailing
   reconstruction, automasking, nonuniform filtering, and output policy.
8. Release state-specific work objects before advancing to the next state.

It is not necessary to assemble a combined `nspace`-sized 2D stack in
`volassemble`. Both 2D accumulation and 3D insertion are linear, so one complete
part-local key record can be inserted directly into the current-state 3D
accumulator.

The compact branch should not allocate the Cartesian `eorecvol_read` buffer and
should not read or write legacy Cartesian partial numerator/rho volumes.

The state loop is also the required allocation boundary. Before advancing to
the next state, the compact branch must release or reset all reconstruction,
restoration, mask, and nonuniform-filter objects whose contents belong only to
the completed state. Arrays of lightweight per-state metadata such as
populations, resolutions, filenames, and update fractions are acceptable;
arrays of allocated 3D images or reconstructors are not.

## Ownership Boundaries

### `simple_strategy3D_matcher`

Owns:

- completion of registration/state assignment
- the assignment-freeze barrier
- the exact selected particle subset
- transition from alignment batching to reconstruction indexing

It must not define the compact file layout or perform volume assembly.

### Reconstruction grouping helper under `strategies/search`

Owns:

- building the `(state, eo, proj)` group index
- reconstruction key iteration
- memory-bounded key batching
- I/O chunk iteration within a complete key

This helper should be independent of individual 3D search strategies.

### `fourier_2d_accumulator`

Owns:

- raw 2D numerator/CTF-squared accumulation
- residual in-plane rotation
- current normalized 2D Kaiser-Bessel interpolation
- export of an un-restored compact Fourier plane

It must not know about distributed parts, state scheduling, or container
filenames.

### Compact-part container helper

Owns:

- format versioning and validation
- append-once record writing
- state/key indexing
- record streaming for assembly
- atomic completion

This is a reconstruction artifact and should be reusable by both integrated
`refine3D` and standalone `reconstruct3D`.

### `volassemble`

Owns:

- selecting Cartesian or compact partial input from `projrec`
- streaming compact records state by state
- 3D Kaiser-Bessel insertion
- all existing reconstruction restoration and postprocessing policy

## Required Invariants

1. The reconstruction subset is exactly the subset selected for the iteration.
2. Assignment fields are immutable after reconstruction indexing starts.
3. Every valid selected particle contributes exactly once.
4. Every populated `(state, eo, proj)` key produces exactly one record per
   distributed part.
5. A record is written only after all local particles for its key have been
   accumulated.
6. Raw numerator and CTF-squared data use identical 2D interpolation weights.
7. No 2D CTF-density correction, normalization, FFT/IFFT round trip, or
   real-space conversion occurs.
8. `volassemble` applies the existing 3D compact-source scaling exactly once.
9. Empty states and empty keys allocate no Fourier sum storage.
10. `projrec=no` retains the legacy Cartesian reconstruction path.

## Refactoring Sequence

### Phase 1: Separate alignment and reconstruction lifecycles

- Make the assignment-freeze point explicit in `simple_strategy3D_matcher`.
- Ensure final projection labels are derived after all orientation/state
  updates.
- Keep `projrec=yes` reconstruction disabled inside alignment batches.
- Call reconstruction image/Fourier-plane setup without initializing
  `build%eorecvols(:)` on the compact path.
- Remove worker calls to `grid_plane_compact`, `write_partial_recs`, and the
  legacy Cartesian partial-volume lifecycle when `projrec=yes`.
- Split reconstruction-buffer cleanup from Cartesian partial writing so the
  particle-image and 2D Fourier preparation machinery remains reusable without
  constructing volumes.

### Phase 2: Assignment-ordered grouping

- Add the counting-sort/prefix-sum group index.
- Add consistency diagnostics for indexed particles and populated keys.
- Unit-test mixed-state, mixed-even/odd, empty-state, and invalid-label cases.

### Phase 3: Complete-key compact accumulation

- Start with one live key accumulator.
- Read large keys through multiple image I/O chunks.
- Verify that each key is exported once.
- Compare exported sums with the current state-wide accumulator path.

### Phase 4: Single part-container output

- Add the versioned compact format and atomic writer.
- Add round-trip and corruption/incompatibility tests.
- Update stale-artifact cleanup and queue completion ordering.

### Phase 5: State-streamed compact `volassemble`

- Add the `projrec=yes` input branch.
- Stream records directly into one current-state reconstruction.
- Eliminate `eorecvol_read` from the compact branch.
- Assert through profiling and lifecycle checks that no prior or future
  state's 3D reconstruction objects remain allocated while assembling the
  current state.
- Preserve existing FSC, trailing, automask, and nonuniform-filter behavior.

### Phase 6: Controlled parallelism

- Measure single-key throughput and peak RSS first.
- Add a small number of complete-key lanes only if worker time is limiting.
- Keep lane count memory-budgeted and deterministic.

## Validation Matrix

Validate against the current working `projrec=yes` implementation for:

- one and multiple states
- uneven state populations and empty states
- C1 and nontrivial point-group symmetry
- even/odd and combined-even/odd reconstruction
- raw and denoised reconstruction sources
- ML regularization
- fractional update and trailing reconstruction
- shared-memory and distributed execution
- ordinary and high-`nspace` neighborhood modes
- nonuniform filtering and automasking

Compare:

- per-key raw Fourier numerator and CTF-squared sums
- current-state 3D numerator and sampling density before restoration
- restored half-map scale and norms
- FSC curves and reported resolutions
- final reconstruction quality
- worker and `volassemble` peak RSS
- compact artifact size and I/O time

Floating-point differences caused by a changed summation order are acceptable,
but tolerances must be defined before replacing the current path.

## Acceptance Criteria

- The alignment phase continues to use its current batching policy.
- Reconstruction grouping occurs only after assignment is final.
- Worker reconstruction peak RSS does not grow with `nstates`.
- The worker never initializes `build%eorecvols(:)` for `projrec=yes`.
- The worker never performs 3D insertion or writes Cartesian partial
  reconstruction volumes for `projrec=yes`.
- Only a memory-bounded number of 2D key accumulators is resident.
- Each populated key is written once per distributed part.
- One compact file is produced per part, avoiding per-key and per-batch file
  clutter.
- `volassemble` remains state-serial and owns only one current-state expanded
  reconstruction.
- `volassemble` owns no array of allocated state volumes; any state-indexed
  arrays contain metadata only.
- Final numerical quality remains equivalent to the validated current
  `projrec=yes` path.
