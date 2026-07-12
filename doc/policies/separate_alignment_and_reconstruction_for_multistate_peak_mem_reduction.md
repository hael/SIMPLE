# Separate alignment and reconstruction for multi-state peak-memory reduction

## Policy status

This is the active memory-lifetime policy for normal CPU multi-state
reconstruction in `refine3D`, `abinitio3D`, and the shared reconstruction
machinery used by `reconstruct3D`. It applies to both `projrec=no` and
`projrec=yes`.

## Implementation status (2026-07-12)

This policy is implemented for the normal CPU paths:

- `simple_strategy3D_matcher` completes all alignment batches without
  constructing or updating reconstruction volumes;
- alignment images, strategies, batch workspaces, and the complete PFTC are
  destroyed at the assignment barrier;
- the now alignment-only strategy3D toolbox is also destroyed at that barrier;
- the matcher verifies that PFTC teardown completed before reconstruction;
- benchmark output for the representative worker records current RSS after
  alignment teardown and after reconstruction, alongside peak RSS;
- `projrec=yes` finalizes discrete projection labels before worker orientation
  metadata is written;
- `calc_3Drec` and `calc_projdir3Drec` group the final selected particles by
  hard state and report the resulting state populations;
- both routines use `init_rec` and one singleton `build%eorecvol` at a time;
- each state partial is written through the existing Cartesian artifact
  contract and the singleton reconstructor is killed before the next state;
- the complete Euclidean sigma object remains available through reconstruction
  for both standard and non-standard consumers.

`build%eorecvols(:)` and its legacy multi-volume reconstruction helpers have
been removed. The OpenMP-offload path uses the singleton for a single state
and dispatches multi-state work to the common state-homogeneous CPU kernel.
Compact 2D part output remains a later refactor.

## Policy

Multi-state execution must separate alignment from partial reconstruction,
then reconstruct one state at a time from state-homogeneous particle batches.

This policy preserves the current worker output contract: every distributed
part still writes one Cartesian even/odd partial reconstruction for each
populated state, and `volassemble` continues to consume those files unchanged.
It does not introduce the compact projection-direction container described in
[`projrec_assignment_ordered_batching.md`](../refactoring_notes/projrec_assignment_ordered_batching.md).

The worker loop becomes:

```text
Phase A: align all selected particles
    use existing alignment batches
    update orientation, shift, state, and projection metadata
    perform no partial-reconstruction allocation or update

freeze assignments
build a state-grouped particle index

Phase B: reconstruct state 1
    allocate one reconstructor_eo
    process only particles assigned to state 1
    write the existing state-1 Cartesian partial
    destroy the reconstructor

Phase B: reconstruct state 2
    reuse the same lifecycle
...
```

Peak worker reconstruction memory therefore contains one state reconstruction,
not one reconstruction for every populated state.

The phase boundary also permits alignment memory to be released before the
first reconstruction object is constructed. In particular, the complete PFTC
object and its all-state reprojection model are alignment-only. The intended
worker peak becomes approximately:

```text
max(alignment working set, one-state reconstruction working set)
```

rather than the two working sets overlapping.

## Motivation

Before this policy was implemented, the `projrec=no` matcher initialized
reconstruction before entering the alignment loop. `init_rec` called
`preprecvols`, which constructed
`build%eorecvols(state)` for every populated state. After each alignment batch,
`maybe_restore_batch` prepared reconstruction Fourier planes and inserted them
into the corresponding state volume.

The corresponding legacy `projrec=yes` path delayed reconstruction until
alignment was complete, but `calc_projdir3Drec` also called `init_rec` with
volume initialization enabled. It accumulated projection-direction sums state
by state, then inserted them into an already allocated array of state
reconstructors.

Both paths therefore retain expanded reconstruction objects with memory
approximately proportional to:

```text
nstates * box^3
```

Neither path requires that lifetime. Final state assignments are available
after alignment, and reconstruction is linear over the particles assigned to
each state.

## Goals

- Complete registration and hard state assignment before any partial-volume
  reconstruction begins.
- Freeze the exact selected particle set and its final reconstruction metadata.
- Group the selected particles by final state.
- Allocate and initialize only one `reconstructor_eo` on a worker at a time.
- Write the existing Cartesian partial-reconstruction artifact for the current
  `(state, distributed part)` before advancing to the next state.
- Apply the same lifetime policy to `projrec=no` and `projrec=yes`.
- Keep `volassemble` and its current Cartesian input contract unchanged.
- Make worker partial-reconstruction memory `O(box^3)`, independent of
  `nstates`.
- Destroy the all-state reprojection model, particle PFTs, correlation
  memoization, and search workspaces at the assignment barrier.
- Construct reconstruction images, Fourier planes, interpolation caches, and
  the current-state reconstructor only after alignment memory has been
  released.
- Minimize simultaneous live objects within Phase B by using state-local and
  I/O-batch-local lifetimes.
- Preserve particle selection, interpolation, CTF, even/odd, symmetry, and
  restoration conventions.

## Non-goals

- This refactor does not change the all-state reprojection model used for
  alignment.
- It does not change any 3D search procedure, objective, or stochastic search
  order.
- It does not reorder alignment batches.
- It does not introduce compact 2D partial files or move
  `grid_plane_compact` into `volassemble`.
- It does not make projection-direction-homogeneous batches. State grouping is
  sufficient for the first memory reduction.
- It does not reduce the single-state peak inside `volassemble` or redesign
  nonuniform filtering.
- It does not change `projrec=no` or `projrec=yes` numerical conventions.

## Required phase boundary

The matcher must expose two sequential phases.

### Phase A: alignment and assignment

```text
for alignment_batch in existing_alignment_batches
    read and prepare matching images
    run the existing 3D search strategy
    update orientation and shift
    update hard state assignment
    update sigma statistics when required
end for

finalize projection-direction labels when required
write orientation results according to existing worker policy
freeze reconstruction metadata
```

This phase must not:

- construct any reconstruction volume through `init_rec`;
- allocate `ptcl_rec_imgs` solely to update partial volumes;
- call `prep_imgs4rec` or any particle-to-volume insertion routine;
- construct or reset any state reconstructor;
- insert any Fourier plane into a reconstruction.

The reprojection model and matching image buffers can be released after the
assignment freeze if they are not needed by reconstruction. Reconstruction
must use the final in-memory orientation field, including final shifts and
state labels.

The assignment barrier must perform an explicit alignment teardown before
Phase B allocation begins. It is not sufficient to rely on final matcher
cleanup after reconstruction.

### Phase B: state-homogeneous partial reconstruction

Build a grouped index over the exact `pinds` selected in Phase A:

```text
state_counts(1:nstates)
state_offsets(1:nstates+1)
grouped_pinds(1:nvalid)
```

A counting/prefix-sum grouping is sufficient:

1. Count valid selected particles by final state.
2. Prefix-sum the counts.
3. Scatter the selected particle indices into `grouped_pinds`.
4. Validate that every eligible selected particle occurs exactly once.

State-zero and invalid-state particles must follow the existing reconstruction
rejection policy and must be reported in the grouping diagnostics.

For each populated state:

```text
initialize one current-state reconstructor_eo
reset its even and odd accumulators

for io_batch in grouped_pinds(state_range)
    reread the reconstruction-source particle images
    prepare reconstruction Fourier planes using final orientations/shifts

    if projrec=no
        grid each particle directly into the current-state reconstructor
    else
        accumulate the current state's projection-direction 2D sums
        grid completed compact sums into the current-state reconstructor
    end if
end for

compress and write the existing Cartesian partial for this state and part
destroy/reset all state-specific reconstruction storage
```

Only after destruction of the current state may the next state reconstructor
be initialized.

## Phase-lifetime policy

Every substantial object should belong to one of three lifetime classes:

1. alignment-only;
2. retained metadata shared across the barrier;
3. reconstruction-only, preferably state-local or I/O-batch-local.

Objects must not remain live merely because they are members of the shared
`builder`. Ownership should follow the last actual consumer.

### Alignment-only objects: release at the barrier

The following objects have no reconstruction consumer and should be destroyed
after final assignments and any alignment-derived sigma output have been
committed:

- `build%pftc`, including:
  - even and odd raw reprojection PFTs;
  - all four memoized reference banks;
  - particle PFTs and their memoized forms;
  - polar CTF matrices;
  - per-thread correlation and shift-search workspaces;
- the global `s3D` strategy allocation through `clean_strategy3D`;
- per-particle strategy objects and specifications;
- probabilistic assignment/search objects after their results are consumed;
- alignment search-order, subspace, and in-plane-search workspaces;
- `ptcl_match_imgs` and `ptcl_match_imgs_pad`;
- the alignment `build%imgbatch` allocation;
- alignment-only batch ranges, counters, and incremental-shift buffers;
- temporary reference projectors or volumes used to prepare the reprojection
  model.

`polarft_calc%kill` already releases both raw and memoized reference banks,
particle banks, CTF matrices, and per-thread heaps. It should therefore be
called at the barrier, not after partial reconstruction as in the current
matcher ordering.

The strategy3D toolbox no longer owns reconstruction storage. It may therefore
be destroyed at the barrier once no remaining alignment consumer needs its
in-plane rotation metadata.

### Metadata retained across the barrier

Reconstruction still requires a comparatively small set of authoritative
metadata:

- `params` and command-line reconstruction policy;
- the selected `pinds` and the new state-grouped index;
- `spproj` and `spproj_field`, including stack/source information, CTF
  parameters, final orientation, shift, state, projection, and even/odd labels;
- point-group symmetry;
- reconstruction masks and sampling/box metadata;
- eulspace for `projrec=yes` projection-direction insertion;
- output naming and distributed part metadata.

The direct `projrec=no` kernel does not intrinsically need the reprojection
eulspace, but retaining it for the first implementation is acceptable because
it is small compared with the PFTC reference banks. A later path-specific
cleanup may release eulspace and subspace maps when the selected reconstruction
kernel does not consume them.

### Euclidean sigma state retained across the barrier

`build%esig` is intentionally retained until reconstruction is complete.
Alignment uses it to calculate and write Euclidean sigma updates, while the
standard ML-regularized reconstruction currently uses:

```text
build%esig%sigma2_noise
build%esig%get_kfromto()
```

The sigma arrays are not a significant contributor to the multi-state memory
problem, and retaining the complete object leaves them available to future or
non-standard reconstruction kernels. This is a deliberate extensibility
policy, not merely an ML-regularization exception.

Retaining `build%esig` must not retain `build%pftc`. PFTC teardown nullifies
its pointer to `sigma2_noise`, while the underlying arrays remain independently
owned by `build%esig`. The sigma object is destroyed only after the final state
reconstruction has been written.

### Reconstruction-only objects: construct late

The following allocations should not exist during alignment:

- the singleton current-state `reconstructor_eo`;
- `ptcl_rec_imgs` or the equivalent reconstruction-source image batch;
- padded reconstruction image heaps;
- `fplane_type` arrays;
- reconstruction FFT/logical-address maps;
- `projrec=yes` 2D numerator and CTF-squared accumulators;
- state-local population/projection lookup workspaces.

Construct the reusable reconstruction I/O buffers only after alignment
teardown and size them from the reconstruction batch policy, not the alignment
batch maximum. This permits a smaller reconstruction batch when memory is more
valuable than I/O throughput.

Interpolation address maps should be memoized only immediately before their
first reconstruction consumer and forgotten when reconstruction is complete.
They must not be constructed as a side effect of alignment preparation and
kept alive across the barrier.

### State-local teardown

At the end of each state:

1. compress and write the current Cartesian partial;
2. close and verify the output artifact;
3. destroy/reset the current `reconstructor_eo`;
4. destroy all `projrec=yes` projection-direction accumulators and maps;
5. release any state-local particle-index or population workspace;
6. reuse, rather than multiply, the bounded reconstruction image/Fourier-plane
   buffers for the next state.

The reusable image and Fourier-plane buffers may remain allocated across
states because their size is independent of `nstates`. State-dependent cubic
or dense projection-direction data may not.

### End-of-reconstruction teardown

After the final state partial is committed:

- release the reusable image batch, padded-image heap, and Fourier planes;
- forget reconstruction FFT maps;
- release the retained Euclidean sigma object;
- release the grouped particle index;
- perform ordinary builder/workflow cleanup and only then signal worker
  completion.

Deallocation may not reduce the process high-water RSS reported by
`get_peak_rss_bytes`, and a system allocator may retain freed pages for reuse.
The correctness requirement is that alignment and reconstruction allocations
are not simultaneously live. Validation should therefore record current RSS
or explicit allocation accounting at the barrier as well as final peak RSS.

## `projrec=no` policy

The Cartesian path should use the existing particle-to-volume numerical
machinery, but its entry point should accept one state-homogeneous particle
range and one current-state reconstructor.

The important behavioral change is I/O lifetime. Today, reconstruction images
can be prepared alongside matching images and inserted immediately after each
alignment batch. After this refactor, particle images are reread during Phase B
in state order. This adds reconstruction I/O and preprocessing, but removes the
state-multiplied volume allocation and prevents reconstruction from observing
intermediate assignments.

The state-homogeneous routine receives one current-state reconstructor. The
caller has already selected one state, and all valid input particles must
belong to it. A mismatched state is an invariant violation.

## `projrec=yes` policy

The projection-direction path should retain its current numerical sequence for
this first refactor:

1. Accumulate raw 2D Fourier numerator and CTF-squared sums for the current
   state.
2. Insert those sums through the current `grid_plane_compact` convention.
3. Write the resulting existing Cartesian partial reconstruction.

The only required memory change is that the insertion target is the singleton
current-state reconstructor. No other state's 2D sums or 3D reconstruction may
be resident.

This is intentionally an intermediate design. A later refactor may write the
compact 2D sums directly and move 3D insertion to `volassemble`, as specified
in [`projrec_assignment_ordered_batching.md`](../refactoring_notes/projrec_assignment_ordered_batching.md).

## Reconstructor ownership

Worker reconstruction should use a singleton object, conceptually:

```text
build%eorecvol
```

or an equivalent local `reconstructor_eo`. The reconstruction code must not
index an array of live reconstructors by state.

The scientific state number remains an argument used for filenames, metadata,
and state validation. It must not imply persistent ownership of a separate
volume object.

## Output and assembly compatibility

The output filenames and payloads remain the current Cartesian partials:

```text
refine3D_partial_rec_fbody(state, part, numlen_part)
```

Each populated state is compressed and written immediately after its worker
reconstruction completes. Empty states produce no partial, following current
policy.

Because the artifact contract is unchanged, `volassemble` requires no new
reader or reduction algorithm. It already loops over states serially and
reduces the corresponding part files into one current-state reconstruction.

Stale-partial cleanup and distributed completion signaling must retain their
current ordering. A worker must not report successful completion until all of
its populated state partials are closed and visible.

## Memory policy

The intended Phase B worker peak is:

```text
one current-state reconstructor_eo
+ bounded particle image/Fourier-plane buffers
+ projrec=yes current-state 2D projection sums
+ O(number of selected particles) integer grouping index
+ O(nstates) counts and offsets
```

Forbidden allocations include:

- one initialized `reconstructor_eo` per populated state;
- retained reconstruction images from the alignment phase;
- multiple state volumes kept alive to simplify final writing;
- simultaneous projection-direction sums for several states.

The total worker lifetime should follow this allocation timeline:

```text
startup metadata
    -> allocate PFTC/reprojections and alignment buffers
    -> align all selected particles
    -> freeze assignments
    -> destroy PFTC/reprojections and alignment buffers
    -> allocate bounded reconstruction buffers
    -> for each state: allocate one reconstructor, write, destroy
    -> destroy reconstruction buffers
```

No arrow may advance while an object owned exclusively by the preceding phase
remains live.

The `projrec=yes` current-state 2D sums may still scale with populated
projection directions within one state. Reducing that to complete-key batching
belongs to the later compact-output refactor.

## Proposed code decomposition

### `simple_strategy3D_matcher`

Owns:

- the Phase A alignment loop;
- the assignment-freeze barrier;
- creation of the state-grouped index;
- iteration over populated states;
- dispatch to the selected state reconstruction kernel.

It should remove `maybe_init_reconstruction` and `maybe_restore_batch` from the
alignment lifecycle when partial reconstruction is requested.

### State-grouping helper under `strategies/search`

Owns:

- counting and prefix-sum grouping by state;
- ranges into `grouped_pinds`;
- validation and population diagnostics.

It must be independent of `projrec` so both reconstruction kernels consume the
same grouping.

### `simple_matcher_3Drec`

Owns state-local reconstruction:

- initialize one reconstructor;
- allocate bounded reconstruction image/Fourier-plane buffers;
- reconstruct one homogeneous state using either direct or projdir insertion;
- write that state's existing Cartesian partial;
- release the state reconstructor and buffers.

The module exposes state-oriented entry points and never implicitly allocates
or writes an array of state reconstructors.

## Migration sequence

### Phase 1: freeze alignment before reconstruction

- Remove reconstruction initialization from before the matcher batch loop.
- Remove Cartesian partial updates from `maybe_restore_batch`.
- Stop constructing `ptcl_rec_imgs` during alignment.
- Complete all orientation, shift, state, sigma, and orientation-output work.
- Add an explicit assignment-freeze point.
- Move `clean_strategy3D`, probabilistic-object cleanup, matching-image cleanup,
  and `build%pftc%kill` to the assignment barrier.
- Retain `build%esig` across the barrier for standard and non-standard
  reconstruction consumers, without retaining `build%pftc`.
- Add barrier diagnostics for current RSS and the existence/allocation state of
  the principal alignment objects.

### Phase 2: build the state grouping

- Add the counting/prefix-sum grouped index over the selected `pinds`.
- Record state populations and rejected/invalid entries.
- Validate exact-once membership.
- Test empty, highly imbalanced, and state-zero cases.

### Phase 3: singleton Cartesian reconstruction

- Add a `projrec=no` state-local reconstruction entry point.
- Construct reconstruction image buffers and Fourier planes only after the
  alignment teardown.
- Allow reconstruction batch size to be selected independently of alignment
  batch size.
- Reread images in state-homogeneous I/O batches.
- Allocate, write, and destroy one state reconstructor at a time.
- Preserve the current Cartesian partial format and scaling.

### Phase 4: singleton projection-direction reconstruction

- Adapt `calc_projdir3Drec` to accept one state range and one reconstructor.
- Allocate projection-direction sums only for the current state.
- Insert through the existing compact gridding convention.
- Write and destroy the current state before advancing.
- Release all current-state projection lookup and accumulator storage before
  initializing the next state.

### Phase 5: remove obsolete array ownership (complete)

- Removed `build%eorecvols(:)` from the strategy3D toolbox.
- Removed the legacy multi-volume reconstruction helpers.
- Converted the OpenMP-offload single-state path to singleton ownership and
  routed its multi-state case to the common state-homogeneous CPU kernel.
- Keep `build%eorecvol` as the explicit current reconstruction object.

## Validation

For both `projrec=no` and `projrec=yes`, compare before and after for:

- one state and multiple states;
- balanced and highly imbalanced state populations;
- empty states and state-zero particles;
- shared-memory and distributed workers;
- raw and denoised particle sources;
- even/odd half assignment;
- ML regularization;
- C1 and nontrivial point-group symmetry;
- trailing reconstruction and fractional updates;
- ordinary and nonuniform `volassemble` processing.

Compare:

- selected particle counts by state and even/odd half;
- pre-compression numerator and sampling-density volumes;
- written Cartesian partial payloads;
- restored half maps and merged volumes;
- FSC curves and reported resolution;
- final reconstruction quality;
- worker peak RSS and reconstruction I/O time.
- current RSS or live-byte accounting immediately before and after the
  assignment teardown;
- allocation-state assertions showing that PFTC reference and memo banks are
  absent before the first state reconstructor is initialized.

Changed particle summation order may produce floating-point differences.
Numerical tolerances must be defined, while particle membership, state
ownership, and interpolation conventions must remain exact.

## Acceptance criteria

- Alignment completes without allocating or updating partial reconstruction
  volumes.
- Reconstruction consumes the final frozen assignments.
- The all-state PFTC/reprojection model is destroyed before any state
  reconstructor is initialized.
- Alignment-only images, PFTs, CTF matrices, memoization, and search workspaces
  are absent during Phase B.
- Reconstruction-only buffers are constructed after the assignment barrier.
- The exact selected reconstruction subset is grouped once by final state.
- At most one initialized state `reconstructor_eo` exists on a worker.
- Worker reconstruction peak does not grow as `nstates * box^3`.
- `projrec=no` and `projrec=yes` both follow the state-homogeneous lifecycle.
- Existing Cartesian partial filenames and payload conventions are preserved.
- The complete Euclidean sigma object remains available throughout
  reconstruction without retaining the PFTC object.
- `volassemble` operates without modification to its partial-reduction
  contract.
- Reconstruction results remain equivalent within predefined floating-point
  tolerances.
