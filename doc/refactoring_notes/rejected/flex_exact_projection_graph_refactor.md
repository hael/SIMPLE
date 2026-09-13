# Refactoring Proposal: Exact-Projection Sparse Graphs for FLEX Analysis

**Status:** Proposal for review only  
**Implementation status:** No implementation changes authorized or made  
**Date:** July 28, 2026  
**Scope:** FLEX graph construction, component-aware diffusion embedding, output contracts, distributed execution, and downstream safeguards

## 1. Executive summary

This refactor replaces FLEX's expanding angular candidate search with an exact projection-assignment rule:

\[
\mathcal{C}_i = \{j \ne i : p_j = p_i\},
\]

where \(p_i\) is the projection-direction assignment inherited from the orientation search. For every particle, FLEX will compare residual features only with particles having the identical assignment. It will retain the nearest `k_nn` members of that set, with a proposed default of `k_nn=50`. There will be no cross-projection graph edges and no fallback that expands to neighboring projection directions.

The graph remains globally stored and sparse, but it is intentionally a disconnected union of projection-specific components. Consequently, the current single global diffusion-map eigensolve, global ICM rank selection, global k-medoids, global kernel weights, and immediate 3D preimage reconstruction are no longer mathematically valid. The refactor must therefore perform diffusion embedding independently within each projection assignment and must block the current global state-volume path until a separate projection-aware latent-alignment model is available.

The recommended implementation is staged:

1. Implement exact-assignment graph construction and component-wise embeddings.
2. Produce per-projection coordinates, spectra, and diagnostics.
3. Validate whether same-projection embeddings recover known flexibility on simulation and reproduce across independent halves.
4. Only then implement a graph-independent, projection-aware alignment of the local coordinates for global 3D reconstruction.

This staging tests the central hypothesis - that cross-view comparisons are suppressing FLEX - without introducing unvalidated coordinate stitching or misleading state volumes.

## 2. Decisions fixed by this proposal

The following are design requirements, not open implementation choices:

| Requirement | Decision |
|---|---|
| Neighbor eligibility | A particle may consider only particles with the identical projection-direction assignment |
| Cross-projection graph edges | Prohibited |
| Neighbor retention | Retain all selected nearest `k_nn` intra-assignment neighbors |
| Initial `k_nn` | Change FLEX default from 100 to 50 |
| Candidate population | Examine every other particle in the same assignment by default |
| Candidate cap | Do not introduce in the first implementation; reserve an optional cap only if profiling demonstrates a need |
| Projection grid | Inherited from the orientation search; FLEX must not tune or redefine `nspace` |
| Sparse representation | Preserve sparse k-nearest-neighbor CSR storage |
| Disconnected topology | Treat as intentional and process explicitly |
| State volumes | Do not construct global states from unaligned component coordinates |

## 3. Clarification of terminology

“Identical projection direction” in this proposal means identical discrete projection assignment from the orientation search. It does not mean that the continuous Euler directions are numerically identical.

The assignment is nevertheless the correct comparison unit because the registration and residual construction for particles in a given assignment use a compatible in-plane transformation and the same assigned mean projection. Particles assigned to different projection directions must not be compared in registered-image space.

Two distinct sizes must be kept separate:

- `n_p`: number of active particles assigned to projection component `p`.
- `k_nn`: maximum number of nearest same-assignment neighbors retained for each particle.

The candidate set normally has size `n_p-1`; the stored graph retains at most `k_nn` outgoing neighbors per particle before symmetric union.

## 4. Current behavior that must change

### 4.1 Projection assignments are currently recomputed

`prepare_flex_diffmap_features` currently builds an explicit Euler grid and assigns each active particle to the closest grid direction with `find_closest_proj`. The comment states that stored `proj` values may have been generated using a different `nspace`.

Relevant source:

- `src/main/flex/simple_flex_diffmap_features.f90`, approximately lines 169-182

Required change:

- FLEX must consume the projection assignment associated with the orientation search used for registration.
- FLEX must not silently create a different grid or reassign particles on a new FLEX-selected `nspace`.
- If the workflow still needs `nspace` to reconstruct the orientation-search grid or mean projection stack, it must be inherited as workflow metadata/plumbing and validated against the stored assignments.
- A user-provided FLEX override that changes the orientation-search grid must be rejected.
- If assignment provenance is absent or inconsistent, fail with an explicit diagnostic rather than silently rebinning.

### 4.2 Angular bins currently expand until a particle cap is reached

`find_gated_euclidean_neighbors_rows` currently:

1. sorts projection bins by angular separation;
2. visits bins in that order;
3. evaluates candidate residual distances;
4. stops when `nang_nbrs` candidates have been examined;
5. retains the closest `k_nn`.

Relevant source:

- `src/main/pca/simple_diff_map_graphs.f90`, approximately lines 90-205

Required change:

- Remove angular sorting and the direction-expansion loop from the FLEX graph path.
- For row `i`, iterate only over the members of component `proj_ids(i)`.
- Examine all `n_p-1` eligible particles in phase 1.
- Retain the closest `min(k_nn,n_p-1)` candidates.
- Never expand to another projection assignment when `n_p-1 < k_nn`.

### 4.3 Current defaults encode the expanding search

Current FLEX defaults are:

- `k_nn=100`
- `nang_nbrs=1000`
- `view_balance=yes`
- `dm_alpha=0`

Relevant source:

- `src/main/flex/simple_flex_analysis_strategy.f90`, approximately lines 142-158
- `src/main/ui/simple/simple_ui_denoise.f90`, approximately lines 350-376

Required change:

- Set the FLEX default `k_nn=50`.
- Remove `nang_nbrs` from the active exact-assignment graph interface.
- Do not reinterpret `nang_nbrs` as an intra-component cap; use a separately named parameter if a cap is later justified.
- Remove or disable graph-stage `view_balance` for component-wise embeddings. Within a single component, inverse occupancy is constant and cancels under normalization.
- Retain any angular balancing needed for later 3D reconstruction as a decoder-level concern, not a graph-affinity concern.

## 5. Proposed graph contract

### 5.1 Mathematical definition

For particle `i` with residual feature vector \(x_i\) and projection assignment \(p_i\):

\[
d_{ij}^{2} = \|x_i-x_j\|^{2}
\quad \text{only if} \quad p_i=p_j.
\]

Let:

\[
k_i = \min(k_{\mathrm{nn}}, n_{p_i}-1).
\]

The directed neighbor set \(N_i\) contains the `k_i` smallest eligible distances. All of these directed neighbors are retained. The final undirected graph uses the current symmetric-union rule:

\[
(i,j) \in E
\quad \text{if} \quad
j \in N_i \ \text{or}\ i \in N_j.
\]

The edge kernel is evaluated within each projection component. No edge may satisfy \(p_i \ne p_j\).

### 5.2 Component occupancy policy

Because `nspace` is fixed by orientation search, FLEX cannot guarantee uniform occupancy. The graph must support variable effective neighbor counts.

Recommended policy:

| Component occupancy `n_p` | Effective graph behavior |
|---:|---|
| `n_p = 1` | No graph edge; component is reported and excluded from embedding |
| `n_p = 2` | One mutual edge; report as non-embeddable for diffusion coordinates |
| `3 <= n_p <= k_nn` | Retain all `n_p-1` eligible neighbors; embed with at most `n_p-2` nontrivial modes |
| `n_p > k_nn` | Retain the closest `k_nn` eligible neighbors |

No component may borrow particles from another assignment.

The run summary must report:

- occupied component count;
- minimum, median, mean, and maximum `n_p`;
- components below 3 particles;
- components below `k_nn+1`;
- particle fraction in underfull components;
- effective `k_i` distribution;
- connected-component count within every projection assignment.

### 5.3 Candidate cap

The first implementation should compare all same-assignment candidates. Its computational complexity is:

\[
O\left(d \sum_p n_p^2\right)
\]

for feature dimension `d`, while graph storage remains:

\[
O(N k_{\mathrm{nn}}).
\]

The finely sampled angular grid should keep individual `n_p` small enough for this to be practical.

If profiling later demonstrates that one or more components are too large, add an optional parameter named `max_intra_candidates`, with:

- default `0`, meaning all same-assignment candidates;
- no change to projection eligibility;
- deterministic, row-order-invariant candidate selection;
- mandatory reporting of affected components and sampling fractions;
- validation that results are stable relative to exhaustive same-component search.

A naive “first N rows” cap is prohibited because it recreates acquisition-, frame-, or state-order bias.

## 6. Recommended implementation structure

### 6.1 Generic partitioned graph builder

Add a partition-aware graph API rather than modifying the semantics of the current angularly gated routine:

```fortran
build_partitioned_euclidean_knn_graph(
    features,
    partition_ids,
    k_nn,
    graph,
    component_summary,
    bandwidth_mode,
    bandwidth_tune
)
```

Suggested helper:

```fortran
find_partitioned_euclidean_neighbors_rows(
    features,
    partition_ids,
    k_nn,
    rows,
    nbrs,
    d2s,
    nneighbors,
    ncandidates
)
```

Use the generic name `partition_ids` inside `simple_diff_map_graphs.f90`; the FLEX caller supplies projection assignments.

Do not reuse `find_gated_euclidean_neighbors_rows` with a special value of `nang_nbrs`. Exact partitioning is a different graph contract and should be explicit in names, validation, logs, and tests.

### 6.2 Prefer per-component construction followed by CSR concatenation

The cleanest implementation is:

1. group global feature rows by `partition_id`;
2. build a local Euclidean kNN graph for each group;
3. estimate the kernel scale within that group;
4. map local row indices back to global feature rows;
5. concatenate the component CSR blocks into one global CSR container.

Advantages:

- existing `build_euclidean_knn_graph` already clamps local `k` to `n_p-1`;
- every candidate in the component can be examined exactly;
- bandwidth is estimated on comparable residuals only;
- variable `k_i` is naturally supported;
- block structure is explicit;
- no padded invalid neighbors are needed;
- component-specific diagnostics are easy to produce.

The global `diffmap_graph` can remain the storage container, but its metadata must not imply a single connected diffusion operator. Add or return a component descriptor rather than overloading `graph%k_nn`.

Suggested descriptor:

```fortran
type :: partitioned_graph_summary
    integer :: npartitions
    integer :: nembedded
    integer :: nunderfull
    integer, allocatable :: partition_id(:)
    integer, allocatable :: occupancy(:)
    integer, allocatable :: k_effective(:)
    integer, allocatable :: internal_components(:)
    real,    allocatable :: bandwidth(:)
end type
```

### 6.3 Bandwidth estimation

The current graph estimates one global Ferguson or median bandwidth from all retained distances. That is inappropriate when no cross-projection distance is meaningful.

Required behavior:

- estimate bandwidth independently for each projection component;
- use the requested Ferguson estimator when it is numerically supported;
- fall back to the component's median positive kth-neighbor distance;
- report the estimator, bandwidth, fallback status, and finite-distance count per component;
- never pool residual distances across projection assignments to set a common kernel scale.

## 7. Component-aware diffusion embedding

### 7.1 Why the current global eigensolve cannot be retained

With no bridges, the global graph is block diagonal:

\[
W = \operatorname{diag}(W_1,\ldots,W_P).
\]

The normalized diffusion operator has one trivial eigenvector per connected component. The current `embed_graph` requests one additional eigenpair and removes only one trivial mode. A global call would therefore return projection-component indicator modes rather than conformational coordinates.

Relevant source:

- `src/main/pca/simple_diffusion_maps.f90`, approximately lines 55-128

### 7.2 Required embedding behavior

Add a FLEX-specific component wrapper that:

1. extracts each projection block;
2. verifies its internal connected-component count;
3. calls the existing generic `embed_graph` separately for each connected block;
4. removes the trivial mode independently within every block;
5. performs ICM/eigengap rank selection independently for each block;
6. maps local coordinates back to global particle rows;
7. preserves projection ID and local mode count for every row.

Do not change generic `embed_graph` behavior for unrelated SIMPLE workflows.

For component `p`, the maximum requested nontrivial rank is:

\[
d_p \le \min(\texttt{neigs}, n_p-2).
\]

If the kNN symmetric union fragments a nominal projection component, report that as a graph-quality failure. Do not silently treat the fragments as equivalent conformational coordinates.

### 7.3 Coordinate output contract

The present `flex_embedding_result` assumes one dense global coordinate table with a shared meaning for each column. That contract cannot represent unaligned per-projection modes.

Introduce a separate result type, for example:

```fortran
type :: flex_partitioned_embedding_result
    integer :: nptcls
    integer :: npartitions
    integer :: max_local_modes
    integer, allocatable :: pinds(:)
    integer, allocatable :: partition_id(:)
    integer, allocatable :: local_mode_count(:)
    real(dp), allocatable :: z_local(:,:)
    real(dp), allocatable :: eigvals_local(:,:)
end type
```

Unused coordinate cells must be accompanied by `local_mode_count`; downstream code must never interpret padding as a physical zero.

Recommended text outputs:

`flex_projection_component_summary.txt`

```text
# projection_index occupancy k_effective internal_components selected_modes bandwidth fallback
```

`flex_projection_local_coordinates.txt`

```text
# feature_row raw_particle_index projection_index local_row selected_modes psi1 ... psiM
```

`flex_projection_local_spectra.txt`

```text
# projection_index mode eigenvalue selected
```

## 8. Downstream behavior and safety boundary

### 8.1 Global k-medoids must be disabled

The current preimage path constructs pairwise distances over the global `raw_coords` table and performs one global k-medoids clustering:

- `src/main/flex/simple_flex_diffmap_preimage.f90`, approximately lines 17-54

Distances between coordinates from different projection components have no meaning because:

- eigenvector signs are arbitrary;
- mode order may differ;
- selected rank may differ;
- coordinate scale and nonlinear parameterization may differ;
- different local modes may represent different phenomena.

Therefore, the exact-projection graph mode must not call:

- `select_flex_diffmap_preimages` globally;
- `build_flex_preimage_kernel_weights` globally;
- `reconstruct_flex_diffmap_weighted_states` from unaligned local weights;
- `reconstruct_flex_diffmap_local_linear_states` from unaligned local coordinates.

### 8.2 Local medoids are diagnostic only

Per-component medoids or coordinate quantiles may be written for inspection and simulation truth testing. They must not be pooled into global 3D states until cross-component coordinate identity and sense have been established.

### 8.3 `nano3D` compatibility

`fit_flex_analysis_embedding` currently returns `flex_embedding_result`, which is consumed by `simple_trajectory_chunker` and `nano3D`.

Exact-projection component embeddings must not be silently packed into that global result. Until component alignment is implemented:

- make the partitioned result a distinct API;
- reject requests that require the current global embedding contract;
- leave the existing global path available behind its current mode for controlled A/B comparisons;
- do not change `nano3D` semantics as part of the graph-only phase.

## 9. Future global alignment without graph bridges

No graph bridges are required or proposed. Global 3D reconstruction can instead align local coordinates through the physically compatible projection model.

Let \(u_i^{(p)}\) be the local diffusion coordinates for particle `i` in projection component `p`. Introduce a component-specific latent transform \(R_p\) and shared 3D signed basis volumes \(B_q\):

\[
r_i \approx
C_i P_i
\left(
\sum_q [R_{p_i}u_i^{(p_i)}]_q B_q
\right)
+ \epsilon_i,
\]

where:

- \(r_i\) is the consensus-subtracted particle residual;
- \(P_i\) is the particle projection operator;
- \(C_i\) is its CTF operator;
- \(R_p\) maps a component's local coordinates into a shared latent basis;
- \(B_q\) are jointly estimated signed 3D variability volumes.

For a one-dimensional local coordinate, \(R_p\) initially reduces to a sign and scale. For multiple modes it may include permutation, sign, and a low-dimensional orthogonal transform.

An alternating fit could:

1. normalize local coordinates within each projection component without forcing endpoint occupancy;
2. initialize component transforms deterministically;
3. fit shared basis volumes with the existing projection-aware M-step machinery;
4. infer particle latent coordinates from the shared basis;
5. update each \(R_p\) using a regularized component-wise alignment;
6. repeat under independent-half validation.

This aligns components through a common 3D forward model, not through direct comparisons between incompatibly registered projection images.

The existing `simple_flex_projected_latent_model.f90` contains reusable CTF- and projection-aware M-step/E-step infrastructure, but component-specific transforms and their identifiability constraints would be new work. This should be a separately reviewed second phase.

Do not use histogram matching or independent stretching of every component to `[0,1]`; those operations would force apparent endpoints and recreate a major ManifoldEM false-positive mechanism.

## 10. Distributed execution

Current distributed FLEX:

1. distributes feature preparation;
2. assembles a global feature table on the master;
3. constructs the graph on the master;
4. performs one global embedding;
5. distributes reconstruction.

For the exact-projection graph phase:

- retain distributed feature preparation;
- retain the assembled global feature-row order and mapping;
- group rows by projection assignment on the master;
- parallelize independent component graph construction with OpenMP initially;
- build and embed each component independently;
- write component outputs from the master;
- do not schedule reconstruction until aligned global coordinates exist.

The legacy stage-2 worker path serializes a fixed number of neighbors per row. If retained, it must support variable `k_i` and verify that every serialized neighbor has the same projection assignment. It may be simpler to retire or bypass this unused path for exact-projection mode rather than expand its file contract prematurely.

Shared-memory and distributed runs must produce identical component membership, neighbor identities, distances, spectra within numerical tolerance, and particle-row mappings.

## 11. Parameter and UI changes

Recommended exact-projection interface:

| Parameter | Proposed behavior |
|---|---|
| `k_nn` | Default 50; maximum retained same-assignment neighbors |
| `nang_nbrs` | Not used in exact-projection mode; retain only for the legacy angular-expansion mode during A/B validation |
| `max_intra_candidates` | Deferred; default 0 means exhaustive same-assignment candidates |
| `nspace` | Inherited orientation-search metadata/plumbing; not a FLEX tuning control |
| `view_balance` | Not used during component graph embedding |
| `bandwidth_mode` | Applied independently within each projection component |
| `bandwidth_tune` | Applied independently within each projection component |
| `neigs` | Maximum local modes per component, capped at `n_p-2` |
| `npreimages` | Disabled for unaligned exact-projection component embeddings |
| `preimage_mode` | Disabled until global latent alignment is available |

For transition and A/B testing, introduce a clearly named graph policy:

```text
flex_graph_scope=(angular_expand|exact_projection){exact_projection}
```

The default should change only after the component-aware embedding path and regression tests pass. Output filenames must record the graph policy to prevent accidental comparison of incompatible coordinate semantics.

## 12. File-level refactoring map

| File | Proposed change |
|---|---|
| `src/main/pca/simple_diff_map_graphs.f90` | Add generic partitioned kNN construction; exact same-partition validation; component-wise bandwidth; CSR block concatenation |
| `src/main/pca/simple_diffusion_maps.f90` | Keep generic eigensolver unchanged; optionally expose small helpers for local-to-global index mapping |
| `src/main/flex/simple_flex_diffmap_features.f90` | Consume and validate orientation-search projection assignments; preserve assignment provenance in feature maps |
| `src/main/flex/simple_flex_analysis_strategy.f90` | Add exact-projection policy, default `k_nn=50`, component grouping, per-component embedding/rank selection, diagnostics, and downstream reconstruction guard |
| `src/main/flex/simple_flex_diffmap_preimage.f90` | No phase-1 algorithm change; prevent global calls from partitioned embeddings; optional local diagnostic medoids later |
| `src/main/flex/simple_flex_diffmap_rec3D.f90` | No graph-phase reconstruction; later accept only aligned global latents |
| `src/main/flex/simple_flex_projected_latent_model.f90` | Phase 2 only: add component transforms and alternating projection-aware alignment |
| `src/main/ui/simple/simple_ui_denoise.f90` | Expose graph policy; set exact-mode neighbor description/default; remove `nspace` as a FLEX tuning concept; disable incompatible preimage controls |
| `src/main/nano/simple_trajectory_chunker.f90` | No phase-1 change; reject partitioned embeddings until a reviewed alignment contract exists |
| `src/main/flex/README.md` | Document exact-assignment topology, component coordinate semantics, and reconstruction safety boundary |

## 13. Tests

### 13.1 Unit tests: graph eligibility

1. **No cross-assignment edges**
   - Construct particles with interleaved row order and distinct partition IDs.
   - Assert `partition_id(i) == partition_id(j)` for every CSR edge.

2. **Exhaustive candidate count**
   - With `max_intra_candidates=0`, assert `ncandidates(i)=n_{p_i}-1`.

3. **Nearest-neighbor retention**
   - Compare retained neighbors against brute-force within-component distances.
   - Assert all closest `min(k_nn,n_p-1)` neighbors are present.

4. **Underfull component behavior**
   - Test occupancies 1, 2, 3, 25, 50, and 51 with `k_nn=50`.
   - Verify exclusion/embedding rules and `k_i`.

5. **Row-order invariance**
   - Permute global particle rows.
   - After mapping back to particle IDs, require identical neighbors and distances.

6. **No fallback**
   - Create a component with two particles beside a highly populated angular neighbor.
   - Verify that no neighboring-component particle is added.

### 13.2 Unit tests: component embedding

7. **One trivial mode removed per component**
   - Construct several connected block graphs.
   - Verify that projection indicator vectors do not appear in returned local coordinates.

8. **Internal fragmentation detection**
   - Construct a nominal projection component containing two disconnected clusters.
   - Require an explicit diagnostic or failure according to the reviewed policy.

9. **Rank bounds**
   - Verify selected rank never exceeds `min(neigs,n_p-2)`.

10. **Independent normalization**
    - Rescale the features of one component.
    - Verify that another component's coordinates and spectrum do not change.

### 13.3 Distributed tests

11. **Shared/distributed parity**
    - Compare component membership, retained neighbors, distances, bandwidths, spectra, and coordinates.

12. **Partition-boundary independence**
    - Change qsys row partitions while preserving input particle order.
    - Require identical scientific outputs.

13. **Feature-map provenance**
    - Verify raw particle index, registered stack index, and orientation-search projection assignment survive all part assembly steps.

### 13.4 Scientific validation

14. **IgG per-projection truth recovery**
    - Measure per-component Spearman correlation between local modes and MD frame.
    - Report the best mode, sign-adjusted correlation, occupancy, and half-set reproducibility per projection.

15. **Current-versus-exact graph A/B**
    - Hold features, particles, CTF handling, bandwidth policy, and rank limits fixed.
    - Compare current angular expansion with exact projection components.

16. **Neighbor purity**
    - On simulation, compare the distribution of MD-frame differences among retained neighbors.

17. **Independent-half local embeddings**
    - Fit exact-projection embeddings separately in even and odd halves.
    - Compare local spectral structure and truth recovery without sharing coordinates.

18. **`k_nn` sensitivity**
    - At minimum test 25, 50, and the largest feasible value below component occupancy.
    - The default 50 should be accepted only if it improves truth recovery or stability without fragmenting components.

19. **Order sensitivity**
    - Repeat after particle-row, stack-order, and within-projection permutations.

20. **Performance**
    - Record `sum(n_p^2)`, distance evaluations, peak memory, graph time, embedding time, and largest component.

## 14. Acceptance criteria

Phase 1 is acceptable only when all of the following hold:

1. Every graph edge joins particles with the identical inherited projection assignment.
2. With no candidate cap, every eligible same-assignment particle is examined.
3. The closest `min(k_nn,n_p-1)` neighbors are retained for every particle.
4. No underfull component expands into another projection assignment.
5. Graph results are invariant to particle-row ordering.
6. `nspace` is inherited and validated, not retuned by FLEX.
7. Bandwidth is estimated independently within each component.
8. One trivial eigenvector is removed independently per connected component.
9. Component-indicator eigenvectors are not presented as conformational coordinates.
10. Shared-memory and distributed results agree within numerical tolerance.
11. Exact-projection mode cannot invoke global k-medoids or reconstruct global state volumes from unaligned coordinates.
12. IgG simulation reports per-projection truth recovery and independent-half stability before any claim that flexibility has been modeled.

## 15. Principal risks and mitigations

| Risk | Severity | Mitigation |
|---|---|---|
| Component occupancy is too low for `k_nn=50` | High | Variable `k_i=min(k_nn,n_p-1)`; explicit underfull diagnostics; no cross-projection fallback |
| kNN graph fragments within a projection assignment | High | Component connectivity diagnostics; test lower/higher `k_nn`; do not silently merge fragments |
| Independent local modes cannot be assembled into global states | Critical | Block global preimages; implement separately reviewed projection-aware latent alignment |
| Per-component bandwidths create incomparable eigenvalue scales | Medium | Never compare raw eigenvalue magnitude across components without normalization; report local spectra |
| Stored projection assignments lack grid provenance | High | Validate inherited search metadata and fail on ambiguity |
| Exhaustive same-component comparisons become expensive | Medium | Profile first; add deterministic cap or approximate search only with exhaustive-reference validation |
| Existing callers expect one global dense embedding | High | Introduce a distinct result type and reject incompatible callers |
| Exact gating reduces orientation mixing but does not solve low SNR | Critical scientific risk | Use IgG truth and independent halves to decide whether local denoising is still required |

## 16. Proposed implementation sequence

### Phase 1A: Graph contract

- Add partitioned graph builder.
- Consume inherited projection assignments.
- Set exact-mode default `k_nn=50`.
- Add exhaustive candidate and no-cross-edge tests.
- Add occupancy/connectivity diagnostics.

### Phase 1B: Component embedding

- Embed every valid projection component independently.
- Add ragged/component-aware output type and text artifacts.
- Add local spectra and rank-selection diagnostics.
- Validate shared/distributed parity.

### Phase 1C: Scientific diagnostic

- Run current-versus-exact graph A/B on IgG.
- Measure local truth recovery, neighbor purity, half reproducibility, and `k_nn` sensitivity.
- Decide whether the exact graph resolves the FLEX embedding failure.

### Phase 2: Global latent alignment

- Write a separate design review for component-specific latent transforms and the projection-aware 3D decoder.
- Fit and validate independent halves.
- Re-enable global state or basis-volume generation only after signed-basis and held-out tests pass.

## 17. Review decisions still required

The graph policy itself is fully specified. The following implementation choices should be approved before coding:

1. Whether underfull components with `3 <= n_p <= k_nn` should be embedded using all neighbors, as recommended, or excluded from scientific output.
2. Whether the legacy `angular_expand` mode should remain temporarily for A/B testing or be removed immediately.
3. Whether phase 1 should expose local diagnostic medoids/quantiles, clearly labeled as non-global, or only coordinates and spectra.
4. Whether `flex_partitioned_embedding_result` should live beside `flex_embedding_result` or in a new transport module.
5. Whether global latent alignment is intentionally deferred until the exact-projection IgG diagnostic demonstrates recoverable local conformational signal.

## Final recommendation

Implement phases 1A-1C before modifying the 3D decoder. The exact-projection graph is a clean and scientifically motivated test of the current cross-view mixing hypothesis. Its disconnected structure is not a defect; it is the required consequence of respecting registration compatibility. The corresponding safeguard is to treat its coordinates as local until a separate projection-aware 3D model aligns them.

The first success criterion should therefore be improved and reproducible per-projection recovery of the known IgG coordinate, not the immediate production of visually different global volumes.
