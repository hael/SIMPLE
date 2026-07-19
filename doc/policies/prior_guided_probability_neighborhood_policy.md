# Prior-Guided Probabilistic Reference-Support Policy

This policy defines a late-stage 3D refinement strategy that uses probability-
or likelihood-derived evidence to reduce the reference search space for each
particle. The retained references are a particle-specific support set, not a
geometrically defined neighborhood. The support must remain usable when the
source and consumer stages use different angular resolutions (`nspace`).

The existing 2D rank-list implementation is historical test code only. It is
not required for the 3D design and is not the artifact contract described here.

## 1. Objective

The existing probabilistic workflow condenses a probability table into one hard
assignment:

```text
probability table -> assignment.dat -> matcher -> project orientation update
```

The new artifact is auxiliary:

```text
probability table -> reference-support artifact -> later sparse search
```

The hard assignment remains authoritative. The support artifact may reorder or
restrict later candidate evaluation, but it does not directly update project
orientations, reconstruct volumes, or perform soft-assignment EM.

## 2. Probabilistic model

For particle `i`, let `r = (state, projection)` identify a reference direction
and let `u = (in-plane angle, x shift, y shift)` denote nuisance pose variables.
The current likelihood path returns a noise-normalized Euclidean distance
`d_i(r)` after the configured nuisance-pose handling. It is a negative
log-likelihood scale, not an objective score that needs an ad hoc power:

```text
ell_i(r) = exp(-d_i(r))
```

For an optional reference prior `pi(r)`, the desired reference support is
based on:

```text
q_i(r) ∝ pi(r) * integral[ likelihood_i(r,u) * prior(u|r) du ]
```

If integration over nuisance variables is too expensive, the first
implementation may use the profile likelihood:

```text
q_i(r) ∝ pi(r) * max_u likelihood_i(r,u)
```

The source tables currently use the profile form: the in-plane/shift handling
produces one calibrated distance per reference. A future marginal producer
could replace the likelihood term with an integral over nuisance variables, but
it must retain the same likelihood calibration.

The current 3D probability candidates already contain a reference id, best
in-plane index, shift, and distance, so they support this profile form. A full
joint distribution over reference, in-plane angle, and shifts is not required
for a reference-support policy: the consumer searches nuisance variables
within each retained reference. It becomes necessary only if the policy also
prunes those nuisance dimensions. The support is therefore a calibrated,
truncated reference likelihood measure, not a claim that a full posterior over
SO(3) and all pose variables has been represented.

The implementation must distinguish these two meanings in artifact metadata:
`marginal` or `profile` support. Distances are not probabilities; the artifact
must record the likelihood transform, prior convention, and normalization.

## 3. Geometry-independent support set

For each particle, normalize the valid reference weights over the complete
source candidate set. Use a per-particle minimum shift only for numerical
stability:

```text
dmin = min_r d_i(r)
w_i(r) = exp(-(d_i(r)-dmin)) / sum_r exp(-(d_i(r)-dmin))
```

The shift cancels exactly in normalized weights. With a non-uniform prior,
multiply by `pi(r)` before normalization; do not apply a second power or a
second likelihood calibration.

The preferred support set is the smallest set whose cumulative weight reaches a
configured target, for example `1 - alpha`:

```text
S_i = smallest set with sum_{r in S_i} w_i(r) >= 1 - alpha
```

This set may contain multiple disconnected angular modes. It is therefore not
restricted by geometric proximity: a distant high-probability projection is
retained, while a nearby low-probability region may be omitted.

The consumer should enforce `K_min` and `K_max`, and should fall back to a
broader existing search when the support is flat, under-covered, invalid, or
exceeds `K_max`. Entropy, effective candidate count, top-two margin, retained
mass, and source coverage should be recorded for diagnostics.

The support cardinality should be adaptive. The existing
`detect_peak_thres_fdr` routine provides a suitable first rule: robust
median/MAD score calibration, one-sided lower-tail p-values for distances, and
Benjamini–Hochberg FDR control. For each particle, apply it to the dense
reference score landscape with configured minimum and maximum support sizes.
Retain candidates below the resulting particle-specific threshold, subject to
the retained-mass and `K_min`/`K_max` guards. This preserves multiple separated
low-distance modes and avoids treating a fixed fraction of references as a
posterior neighborhood.

The FDR parameter `PRIOR3D_FDR_Q=0.25` is a false-discovery-rate target for the
statistical thresholding step. It does not mean that 25% of the projection
directions should be retained; the realized support size remains particle-
dependent and is bounded by the support policy.

A sparse table containing only candidates chosen by an earlier geometric or
random search does not by itself define a global posterior support set. The
first 3D producer should therefore use dense `prob` tables, or apply known
proposal-probability corrections when sampling references.

The current 2D implementation is useful as a template for the adaptive count,
but it should not be copied literally: it uses the FDR result to choose how
many lower-distance candidates to retain, then writes only their ids. The 3D
implementation should retain the threshold, scores, coverage, and weights.
The FDR threshold is a support-selection rule, not by itself a calibrated
posterior probability. The artifact must retain both the selected scores and
the normalization/coverage information needed to distinguish a calibrated
credible support set from an adaptive likelihood neighborhood.

## 4. Resolution remapping

An artifact generated at a smaller `nspace` must not be rejected merely because
the consumer uses finer angular sampling. A source projection direction
represents support on the two-dimensional projection-direction surface, not
only one target grid index.

For each source `(state, source projection)` row, construct the source and
consumer Euler projection spaces with the same point-group symmetry and use
the existing symmetry-aware angular distance machinery to map it into the
consumer space:

```text
source projection orientation
        |
        +-- nearest 3 consumer projections on the angular surface
        |
        +-- optional full-to-subspace mask for controlled expansion
```

Selecting three finer projections is a good first remapping stencil because
projection directions form a two-dimensional angular surface. The first
remapping rule should therefore select `N_REMAP_NEIGH = 3` nearest consumer
projections using the equivalent of
`pgrpsym%nearest_proj_neighbors(consumer_eulspace, source_orientation, 3, mask)`.
Distances must be symmetry-aware angular distances, not differences between
Euler-angle components. The three-way expansion is deduplicated across source
rows and states, and each mapped row records its source projection and angular
remapping distance.

Three neighbors should not be assumed to be universally sufficient. A coarse
source projection represents evidence over a local cell, while the three
nearest target points are only a center stencil. The remapper must use the
source-to-target cell map as a coverage diagnostic. If the source cell is wide,
the resolution ratio is large, or the three-point angular span is inadequate,
it should add a bounded number of target projections from the corresponding
`simple_eulspace_neigh_map` mask or fall back to a broader non-prior search.
The current-stage target scores, not equal splitting of the source weight among
three points, determine the final ordering.

The existing `simple_eulspace_neigh_map` machinery remains useful for the
optional expansion path: it maps the finer consumer space to source-like
subspaces and returns a full-resolution mask/list for selected subspaces. The
direct three-neighbor expansion is the default because it bounds candidate
growth; a subspace expansion is a controlled guard when the source support is
too sparse or the resolution ratio is large.

The supported direction is resolution refinement: the source `nspace` must be
less than or equal to the consumer `nspace`. Coarsening is not part of the
workflow and does not need an aggregation rule. The current-stage target
scores determine the final support ordering; the artifact supplies candidate
provenance and a proposal support, not a permanent index map.

The artifact should store source projection orientation data, or enough source
space metadata to reproduce it exactly. Storing only a projection integer is
unsafe if the Euler-space construction changes. The remapping helper must be
cached per `(pgrp, source nspace, consumer nspace)` and tested for symmetry,
duplicate removal, and coverage. The existing symmetry-aware nearest-neighbor
implementation must be audited before reuse: every symmetry branch must
initialize its per-projection distance accumulator before taking minima.

## 5. 3D source and consumer

The initial 3D source is the last dense `refine=prob` stage immediately before
the late `prob_neigh` stages. With the current controller constants this is
stage 5, immediately before `PROB_NEIGH_REFINE_STAGE = 6`:

```text
sample particles
run prob_tab workers
merge dense tables
extract adaptive reference support
perform hard assignment
write assignment.dat and support artifact
```

Extraction occurs after all partition tables have been merged and before the
global probability object is destroyed. It uses exactly the sampled particle
subset selected by `prob_align`; it must not resample particles.

The first consumer is:

```text
prob_neigh_mode=prior
```

It should use the existing sparse PFTC scoring and in-plane/shift refinement
machinery, replacing stochastic or geometry-only reference selection with the
saved support set. The implementation may:

- evaluate the saved fine `(state, projection)` candidates directly;
- map them to coarse subspaces and evaluate those subspaces; or
- use a hybrid of fine candidates plus a controlled coarse expansion.

The first implementation should prefer direct fine-candidate evaluation plus a
small configurable expansion when needed for robustness. Membership in the
support set must remain probability-based, not distance-based.

Fallback must preserve the existing stage policy. If the artifact is absent,
stale, incompatible after remapping, too flat, or insufficiently covered, use
the configured non-prior neighborhood mode rather than failing the refinement.

The controller timing must be explicit:

- the final `prob` stage writes the support artifact only when a later prior
  consumer is requested;
- the following `prob_neigh` stages use `prob_neigh_mode=prior`;
- the source artifact is not regenerated from a `prob_neigh` sparse table; and
- a source-stage iteration is written only after its final dense probability
  table has been merged, not during an earlier iteration of that stage.

For a docked multi-state split, the pre-split artifact belongs to the old
update epoch and cannot be consumed after the split. The controller must insert
or expose a fresh post-split dense `prob` refresh (or an explicitly equivalent
full projection-support producer) before enabling `prob_neigh_mode=prior` in
the new epoch. The intervening `prob_state` stabilization stage may be used as
fallback, but it is not a valid projection-support source unless it produces
projection-level evidence.

## 6. Artifact contract

The support artifact is separate from `assignment.dat`. It should use a
versioned metadata header plus fixed-width, particle-grouped binary rows.

The row file should mirror the efficient 2D direct-access design:

- fixed record length determined by `K_max` and the row schema;
- record address keyed by particle id or a validated particle-slot index;
- one direct-access record per represented particle, with an explicit count;
- separate direct-access units for concurrent OpenMP readers; and
- no text parsing or per-row allocation in the hot consumer path.

A small stream-format metadata sidecar is acceptable and preferred for the
variable header. The direct-access row file remains the performance-critical
artifact. Direct records must carry enough information to reject stale rows,
and sparse particle ids must not be mistaken for valid zero-filled rows.

Header fields:

```text
magic and format version
workflow and source program
source stage/iteration and aggregation window
particle count and represented particle ids
source nstates and source nspace
consumer-remapping version and nearest-neighbor count
box/sampling/cropping metadata when needed
objective and noise/scale metadata
reference prior and orientation-weight convention
support mode: marginal or profile
alpha, K_min, K_max
source coverage/proposal mode
update-epoch identifier
direct-record schema and maximum row capacity
```

Per-particle candidate fields:

```text
particle id
rank and cumulative retained weight
state and projection
raw score or log weight
normalized weight
best in-plane index
x/y shift and has-shift flag
optional coarse subspace id
entropy, effective candidate count, top-two margin
source coverage/touched fraction
source projection orientation and remapping distance when remapped
```

The preferred 3D identity is explicit `(state, projection)`, rather than an
opaque flattened reference id. This makes validation and future projection
remapping safer when `nspace` changes.

## 7. Compatibility and fallback

The consumer must reject or ignore an artifact when:

- represented particle ids do not match the current sampled subset;
- `nstates`, box, sampling, or crop context is incompatible and cannot be
  remapped;
- objective, noise scale, or support mode is incompatible;
- the artifact belongs to an earlier docked multi-state update epoch;
- source coverage is insufficient for a global support interpretation; or
- the artifact is empty, malformed, flat, or otherwise invalid.

An `nspace` change is not itself an invalidation condition. The consumer must
construct the source/consumer projection spaces, apply the resolution-remapping
rule, and validate that every retained source row maps to at least one valid
consumer projection. Rejection is allowed only when the source orientation
metadata, symmetry, or remapping coverage is invalid.

Fallback behavior must be explicit and logged, including the number of
particles that used prior support versus the configured fallback mode, the
realized neighborhood-size distribution, remapping distances, and the number
of source rows expanded.

The fallback must remain bounded. For single-state and
`multivol_mode=independent` workflows, an incompatible prior falls back to the
`prob_neigh_mode=state` coarse-scored neighborhood. For
`multivol_mode=docked`, it falls back to the `prob_neigh_mode=geom` projected
neighborhood. It must never silently expand to a full reference search merely
because the prior artifact is unavailable.

## 8. Ownership and invariants

`simple_commanders_prob.f90` owns sampling, worker dispatch, partition merging,
hard-assignment writing, and support-artifact writing.

`simple_eul_prob_tab.f90` owns dense 3D candidate scoring and support
extraction. `simple_eul_prob_tab_neigh.f90` owns prior-guided sparse candidate
evaluation. Matchers remain responsible for final hard state/projection,
in-plane, shift, and sigma updates.

The following invariants must hold:

- `assignment.dat` remains the hard-assignment handoff;
- support extraction never changes the sampled particle subset;
- support rows are derived only from validated candidate evidence;
- source/consumer resolution changes are handled through explicit angular
  remapping, never by silently reusing projection indices;
- low-probability regions are omitted by the probabilistic support rule, not by
  an implicit geometric-neighborhood assumption;
- multiple separated high-probability modes remain eligible;
- volume assembly and reconstruction consume only the hard updated project;
- prior support never crosses an explicit update-epoch reset; and
- existing `prob` and non-prior `prob_neigh` modes remain unchanged.

## 9. Initial implementation sequence

1. Factor the support record, metadata, and direct-access row I/O into a shared
   3D module. Preserve the 2D fixed-record performance pattern, but add explicit
   header/version/row-count validation.
2. Implement dense 3D profile/marginal score extraction from `eul_prob_tab`.
3. Reuse `detect_peak_thres_fdr` for adaptive per-particle support thresholds,
   then apply retained-mass and `K_min`/`K_max` guards.
4. Implement and test a cached source-to-consumer angular remapper using the
   existing Euler spaces, symmetry distances, three nearest finer projections,
   and optional `simple_eulspace_neigh_map` masks. Fix and unit-test any
   nearest-neighbor accumulator or symmetry-path defects found in that helper.
5. Change the 3D controller so the final `prob` stage writes support and the
   subsequent `prob_neigh` stages consume it; add a fresh post-docked-epoch
   producer before prior consumption.
6. Add `prob_neigh_mode=prior` to parameter validation, UI metadata, commander
   dispatch, and `simple_eul_prob_tab_neigh` candidate construction.
7. Add adaptive fallback logging and tests for coarse-source-to-finer-consumer
   mappings, cell-coverage guards, symmetry, multimodal support, and duplicate
   removal.
8. Validate shared-memory and distributed equivalence, direct-access concurrent
   reads, stale metadata, missing rows, changed resolutions, and broad/flat
   posteriors.

Acceptance requires that a prior-guided run never silently converts an
incomplete candidate sample into a claimed global posterior support set, and
that ambiguous particles fall back to a broader search rather than becoming
overconfidently trapped. Acceptance also requires that a valid prior generated
at one angular resolution remains usable after remapping to another.
