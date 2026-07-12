# Probabilistic alignment table peak-memory refactoring

## Status

Approved and under implementation as of July 2026. The common candidate model,
batch-streamed workers, compact neighborhood merge, and probabilistic matcher
allocation gate are implemented. Scientific equivalence and representative
multi-state peak-RSS validation remain required before this note is marked
complete.

The target is probability-table memory in `abinitio3D`,
`abinitio3D_cavgs`, `refine3D_auto`, `refine3D_multi`, and direct
`refine3D` runs using `prob`, `prob_state`, or `prob_neigh`.

This revision deliberately favors a small common data model over specialized
storage classes for each probability mode.

## Decision summary

The future refactor should introduce only two persistent concepts:

1. one compact record describing an evaluated 3D probability candidate;
2. one particle-oriented candidate store used by all probability modes.

A transient reference-oriented assignment index may be built from that store
when global balancing requires it. It is workspace, not another authoritative
table representation.

`prob`, `prob_state`, and all `prob_neigh` modes should differ in how they
generate candidates, not in the type used to hold candidates.

The main memory changes are:

- make worker candidate storage batch-local;
- make `prob_neigh` storage proportional to evaluated candidates;
- compact dense `prob` candidate records;
- remove repeated particle/state/projection metadata;
- remove full-size work arrays that duplicate the candidate store;
- preserve the complete all-state reprojection model and all current searches.

## Relationship to other memory work

This note concerns probabilistic pre-alignment only. It complements:

- [`separate_alignment_and_reconstruction_for_multistate_peak_mem_reduction.md`](../policies/separate_alignment_and_reconstruction_for_multistate_peak_mem_reduction.md);
- [`projrec_assignment_ordered_batching.md`](projrec_assignment_ordered_batching.md).

It does not change reconstruction ownership or the alignment/reconstruction
barrier.

## Required invariants

### Scientific behavior

- Keep the complete PFTC reprojection model for every state resident until
  probabilistic scoring is complete.
- Do not state-window, shard, regenerate, or reload the reference model during
  search.
- Preserve objective values, CTF handling, noise normalization, in-plane
  sampling, shift seeding, and shift refinement.
- Preserve `prob_assign=legacy|likelihood` semantics.
- Preserve the candidate set and candidate order for
  `prob_neigh_mode=shc|snhc|geom|state`.
- Preserve the current global hard-assignment and state-balancing procedure.
- Keep the final `ASSIGNMENT.dat` contract consumed by
  `simple_strategy3D_prob` unchanged.

### Stochastic behavior

- Preserve random-number call order and per-thread generators.
- Preserve search permutations, `put_last`, early exits, Top-K selection, and
  tie behavior.
- Do not use a second scoring pass to count candidates.

### Parallel behavior

- One OpenMP thread continues to own one particle evaluation at a time.
- Threads use private candidate and search workspaces.
- Threads do not append directly to a shared growable array or shared file
  position.
- Completed particle results are merged in canonical particle order after the
  parallel region or into precomputed non-overlapping offsets.

## Non-goals

- Refactoring the fundamental probability searches.
- Changing particle sampling, update fractions, or state balancing.
- Redesigning PFTC or reducing reprojection-model residency.
- Changing the global assignment algorithm in the first implementation.
- Combining probability-table construction with reconstruction.
- Removing ML sigma data.
- Creating one new candidate/table class for every probability mode.
- Refactoring `simple_eul_prob_tab2D` in the first implementation. It may reuse
  the common structures later.

## Current memory problem

Define:

```text
S      = active states
R      = nspace * S
P      = sampled particles represented by an object
Ppart  = sampled particles in a distributed part
T      = OpenMP threads
Q      = evaluated particle/reference candidates
Qbatch = evaluated candidates in one particle batch
```

The current `ptcl_ref` contains approximately 44 bytes with the default
four-byte component kinds. `eul_prob_tab` stores:

```text
loc_tab(R,P)       approximately 44 R P bytes
state_tab(S,P)     approximately 44 S P bytes
assgn_map(P)       approximately 44 P bytes
```

`state_tab` and `assgn_map` are allocated in roles that do not use them.
`loc_tab` repeats particle, state, and projection metadata in every record.

`eul_prob_tab_neigh` additionally stores:

```text
eval_touched_refs(R,P)    approximately 4 R P bytes
```

It therefore uses roughly `48 R P` bytes for a table that is logically sparse.
It later converts the table to sparse partition files and a sparse assignment
graph, so the dense representation provides no scientific benefit.

The default `refine3D_multi` geometry stage uses `nspace=5000` and
`nspace_sub=500`. It evaluates roughly ten fine projection directions per
state and particle but allocates 5000 slots per state and particle: about one
used slot in 500.

Dense `prob` genuinely evaluates all references. Its global assignment also
allocates `stab_inds(P,R)`, approximately `4 R P` bytes, while `loc_tab` is
still live. The index matrix is part of the current balancing algorithm and is
not an immediate removal target.

## Minimal common data model

### Candidate record

Use one small record for every evaluated 3D probability candidate:

```text
type prob_candidate
    integer :: iref  ! canonical full PFTC reference index
    integer :: inpl
    real    :: dist
    real    :: x, y
    logical :: has_sh
end type
```

This record works for dense `prob`, `prob_state`, and all neighborhood modes.

- Particle identity comes from the owning particle range.
- `iref` is the canonical full reference index
  `(state-1)*nspace+projection`; state, projection, and compact active-reference
  rank are derived through centralized helper routines.
- `icls`, `npeaks`, and `frac` are not candidate fields in 3D.
- The final assignment is converted to the existing `ptcl_ref` only when an
  assignment is selected.

The first implementation should use this readable record rather than a highly
packed set of mode-specific arrays. A later representation change is justified
only if profiling shows this record remains a major peak contributor.

With current four-byte component kinds this record is approximately 24 bytes.
For dense global `prob`, the common records plus the retained `stab_inds`
matrix are therefore approximately `28 R P` bytes, versus approximately
`48 R P` for the current `ptcl_ref` table plus indices.

### Particle-oriented candidate store

Use one flat store for all modes:

```text
type prob_candidate_store
    integer, allocatable        :: pinds(:)
    integer(int64), allocatable :: offsets(:)     ! size P+1
    type(prob_candidate), allocatable :: candidates(:)
    real, allocatable           :: seed_shifts(:,:) ! 2,P
    logical, allocatable        :: seed_has_sh(:)   ! P
end type
```

Candidates for particle `i` occupy:

```text
offsets(i) : offsets(i+1)-1
```

This is the authoritative table representation.

- For dense `prob`, each particle normally owns `R` candidates.
- For `prob_state`, each particle owns one candidate per active state.
- For `prob_neigh`, each particle owns only candidates actually evaluated.

The container does not need mode-specific subclasses. Mode behavior remains in
the existing search/fill procedures.

### Reference-oriented assignment workspace

Global balancing needs candidates grouped by reference. Build a transient index
from the common store:

```text
ref_offsets(R+1)
particle_ids(Q)
candidate_positions(Q)
```

`candidate_positions` points back to the authoritative candidate record. The
workspace may contain sorted distances if the current frontier implementation
requires them, but it must not copy complete candidate records.

The workspace is destroyed after `assgn_map` has been produced.

## Common ownership model

Avoid a constructor that allocates the union of every role. The same candidate
store should support three simple lifetimes:

### Worker

- owns the complete PFTC model;
- owns one particle/PFT batch;
- owns thread-private current-particle candidate buffers;
- owns one completed candidate batch awaiting output;
- does not own a global assignment map or a partition-wide candidate table.

### Global assignment

- owns the merged candidate store for all sampled particles;
- owns the final assignment map;
- temporarily owns the reference-oriented assignment workspace;
- does not own a PFTC model.

### Matcher assignment reader

- owns sampled particle IDs and the final assignment map only.

These can be explicit initialization paths on the existing classes. They do
not require separate data types for dense, state, and neighborhood modes.
Constructor role is expressed by dedicated procedures (`new_worker`,
`new_state`, `new_neigh_global`, and `new_assignment`), not optional ownership
flags such as `with_assignment` or `compact_global`.

## Worker batching policy

Probability workers already batch particle images and particle PFTs. Candidate
storage must use the same boundary.

```text
load complete reprojection model
open partition output

for particle_batch
    prepare particle images/PFTs

    parallel particle scoring
        generate candidates with existing mode-specific search
        store them in one thread-private current-particle buffer
        perform existing candidate refinement
        retain first-touch order and duplicate-update behavior
    end parallel scoring

    append thread-owned chunks carrying particle-local indices
    preserve candidate order within each particle
    release or reuse batch storage
end for

finalize partition output
release PFTC
```

For `prob_neigh`, worker candidate memory becomes:

```text
O(R*T + Qbatch)
```

rather than `O(R*Ppart)`. Direct stochastic modes may legitimately need a
current-particle buffer approaching `R`; that does not justify `R` slots for
every particle in the partition.

For dense `prob`, worker candidate memory becomes `O(R*batchsz)`. Global dense
assignment still needs all candidates, but each candidate is substantially
smaller than `ptcl_ref`.

Per-particle allocatable arrays are not preferred because they create many
small allocations. Use flat thread buffers and one flat batch store.

## Common partition format

Use one chunk format for all 3D probability modes. Do not create
separate dense, state-only, and neighborhood formats.

### File-count policy

Chunks are records inside one append-only stream; they are not separate files.
The hard distributed artifact contract is:

```text
during worker execution:  one open stream per distributed part
after worker completion:  one closed table per distributed part
after global assignment:  one ASSIGNMENT.dat
```

File count is therefore `O(nparts)` and independent of particles, candidates,
states, threads, and particle batches.

The implementation must not create:

- per-batch files;
- per-thread files;
- per-state files;
- per-particle or per-candidate files;
- directory trees that represent sparse arrays as filesystem objects.

Each worker opens its stream once. After every particle batch, the
serial writer appends one large internal chunk to that same open stream. OpenMP
threads only contribute memory buffers and never open their own output files.
The worker patches the final header and closes the stream after the last batch.

After all part streams have been decoded and `ASSIGNMENT.dat` has been
successfully published, the part streams should be removed according to the
normal temporary-artifact cleanup policy. Failed runs may retain them for
diagnosis and restart investigation.

One file per part may itself stress a metadata server at exceptionally large
part counts. If this becomes a measured problem, a later hierarchical reducer
may consolidate groups of completed part streams into a bounded number of
shards. Do not introduce concurrent multi-worker writes to one shared file,
MPI-I/O, HDF5, or another container dependency in the first implementation.

```text
file header:
    nrefs, nptcls
    chunk count and total Q

partition metadata:
    pinds
    seed_nrots, seed shifts, and flags

repeated thread/batch chunk:
    candidate count
    particle-local indices
    prob_candidate records
```

The writer patches final counts, validates byte sizes, and atomically publishes
the file. The reader consumes chunks directly into the common candidate store.

Use 64-bit candidate counts, offsets, and file-position arithmetic. Chunk size
must be bounded by a target byte budget, not a fixed number of particles.

Temporary probability files may change format because they are created and
consumed within one run. They carry no format-version field and support no
reuse of tables from older executions. `ASSIGNMENT.dat` remains unchanged.

## Derived metadata cleanup

In 3D every projection exists for every active state. Replace:

```text
proj_exists(nspace,nstates)
sinds(R)
jinds(R)
```

with small common maps:

```text
active_states(S)
state_to_active_rank(nstates)
```

Central helper procedures convert between:

- canonical full candidate reference;
- compact active-reference rank;
- actual state and projection.

Do not duplicate this arithmetic throughout the probability classes.

## Workspace cleanup

After the common store is established:

- remove dense `dists_refs(R,T)` and select refinement candidates from the
  current particle's candidate list;
- remove `evaluated_ref_ids(R,T)` because candidate IDs are already stored;
- remove `state_eval_dists(R,T)` by scanning the candidate list per state;
- replace `fullref_to_sparse_ref(R)` with the common active-state map;
- allocate in-plane vectors only for objective/mode branches that consume them;
- keep the current stochastic permutation and `ran_tabu` storage unless a
  separate proof shows they can change without RNG drift.

These changes use the common record/store and should not introduce new
mode-specific workspace types.

## Staged implementation

### Phase 0: measurement — implemented separately, validation pending

Add allocation accounting for `R`, `P`, `Ppart`, `T`, `Q`, candidate-count
distribution, and major live byte totals. Include process peak RSS in existing
`l_bench_glob=.true.` output. Distributed execution should continue reporting
one representative part, not one file per part.

### Phase 1: lifetime cleanup — implemented

- Construct the global neighborhood object only after workers finish. It
  currently overlaps the worker object in shared-memory execution.
- Allocate `state_tab` only for state assignment.
- Allocate `assgn_map` only for global assignment and matcher reading.
- Remove derived reference metadata.
- Bound existing I/O buffers by bytes.

This phase does not change candidate representation or file format.

### Phase 2: common candidate record and store — implemented

- Introduce `prob_candidate` and `prob_candidate_store` in one shared module.
- Convert candidate recording/refinement helpers to operate on them.
- Preserve a temporary dense-versus-common comparison path for tests.
- Do not add mode-specific candidate types.

### Phase 3: batch-streamed workers — implemented

- Stream common-store chunks for `prob_neigh`, `prob`, and `prob_state`.
- Remove partition-wide worker tables.
- Validate candidate counts, order, uniqueness, and file sizes per batch.

Implement and validate `geom` first, then exercise the same store with
`state`, `shc`, and `snhc`. This is test sequencing, not separate production
implementations.

### Phase 4: common global merge — neighborhood implemented

- Merge partition chunks into one common candidate store.
- Adapt neighborhood normalization, fallback, and assignment graph creation to
  candidate positions.
- Adapt dense and state assignment to the same record accessors.
- Remove global dense neighborhood and touched-reference matrices.

Dense and state global assignment use the common candidate record in their
existing rectangular access layout. Converting that layout to flat offsets
would not materially reduce its payload and is deferred unless profiling shows
an additional benefit.

### Phase 5: workspace cleanup — partially implemented

- Remove duplicated full-reference/thread arrays.
- Retain `stab_inds(P,R)` for dense global assignment initially.
- Consider changing dense frontier storage only in a later reviewed proposal.

The probabilistic matcher no longer allocates the unused per-thread
`proj_space_shift`, `proj_space_corrs`, `proj_space_inplinds`, or sampling
threshold arrays. Remaining `R*T` scoring workspaces are retained where they
are needed to preserve current Top-K and tie behavior.

### Phase 6: remove validation path and consider 2D reuse

Once exact equivalence is established, remove the old dense neighborhood path
and its temporary comparison switch. Then separately decide whether
`simple_eul_prob_tab2D` should reuse the common record/store.

## Validation

For fixed seeds, inputs, and thread count, compare old and new paths for:

- sampled particle IDs;
- evaluated reference IDs and per-particle candidate order;
- distances, in-plane indices, shifts, and `has_sh`;
- duplicate-update behavior;
- final assignment records;
- final RNG state in instrumented tests.

Required mode coverage:

- `prob_state` and dense `prob`;
- `prob_neigh_mode=geom|state|shc|snhc`;
- `prob_assign=legacy|likelihood`;
- Euclidean and CC objectives where supported;
- shifts on and off;
- single-state, multi-state, and empty-state cases;
- one, two, and production OpenMP thread counts;
- shared-memory and distributed execution;
- `abinitio3D`, `abinitio3D_cavgs`, `refine3D_auto`, and
  `refine3D_multi`.

Temporary partition bytes need not match after format changes. Their decoded
candidate streams must match. The final assignment artifact should be
byte-identical on the same compiler/platform when stochastic state is fixed.

Memory acceptance criteria:

```text
prob_neigh worker candidates    O(R*T + Qbatch), not O(R*Ppart)
prob_neigh global candidates    O(Q), not O(R*P)
dense/state worker candidates   O(candidates per batch), not per partition
thread workspaces               linear in T, with no hidden P multiplier
published table files           exactly one per part, independent of batches
```

## Main risks

- Candidate order changes during thread-buffer merge.
- RNG drift from altered traversal or a second scoring pass.
- Duplicate coarse/fine candidates becoming duplicate records rather than
  updates.
- Sparse candidate lookup increasing runtime.
- 32-bit overflow in candidate counts or file positions.
- Allocator fragmentation from per-particle arrays.
- Hidden dense `loc_tab(iref,iptcl)` dependencies in fallback or balancing.

Mitigations are explicit particle indices in every chunk, preserved
within-particle candidate order, single-pass scoring, common record-update
helpers, candidate-position indices, 64-bit counts/offsets/file positions,
flat buffers, and exact old/new candidate-stream comparison.

## Approved decisions

1. `prob_candidate` is the sole 3D candidate record.
2. Flat particle-oriented storage is used for sparse global assignment; dense
   and state assignment retain rectangular indexing over the same record.
3. All 3D probability modes use one chunked partition format, without a
   version field or backward-compatibility path.
4. Partition artifacts remain exactly one file per distributed part.
5. Exact decoded-candidate and assignment equivalence is the scientific
   acceptance criterion; representative datasets still need to be selected
   for final validation.

## Expected files

Primary:

- `src/main/simple_eul_prob_tab.f90`
- `src/main/simple_eul_prob_tab_neigh.f90`
- `src/main/simple_eul_prob_tab_utils.f90`
- `src/main/commanders/simple/simple_commanders_prob.f90`

Likely one small shared module should own `prob_candidate`,
`prob_candidate_store`, mapping helpers, and chunk I/O. It should not own search
logic or create a parallel hierarchy of mode-specific table classes.

Validation touches the refine3D matcher/orchestration and staged ab initio
controllers. PFTC should remain a consumer through unchanged scoring
interfaces, not a refactoring target.
