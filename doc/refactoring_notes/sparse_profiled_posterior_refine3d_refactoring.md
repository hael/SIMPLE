# Sparse Profiled Posterior for 3D Neighborhood Refinement

## Status

The static sparse-posterior producer, direct-access reader, finer-grid
subspace mapping, and bounded coarse-neighbor exploration are implemented in
the 3D path. Production validation has shown that a static support is not
sufficient: the posterior-guided stage can lose the correct support as the
reference volume changes. The adaptive posterior refresh lifecycle described
below is now partially implemented: merged current candidate evidence is
written back as the next bounded posterior generation. Generation/update-epoch
metadata and the explicit proposal-mixture regularization remain hardening
work before the implementation can be considered production ready.

This note replaces the former 3D `prob_neigh_mode=prior` implementation and
associated code. It does not change the 2D implementation.

The user-facing mode is `prob_neigh_mode=posterior`. The former `prior` mode,
artifact, remapping workspace, controller branches, and reporting branches have
been removed; no compatibility alias is required.

## Decision summary

The final dense `refine=prob` phase should produce a sparse, per-particle,
profiled posterior over state and projection direction. This first artifact is
a seed posterior. Later neighborhood refinement stages should use it to
generate and prioritize a small candidate set on their own angular grid, then
evaluate the current likelihood on those target candidates and refresh the
posterior for the next iteration.

The artifact is therefore:

```text
final dense prob table
        -> sparse profiled posterior artifact
        -> mapped support plus bounded exploration candidates
        -> current target-grid likelihood evaluation
        -> posterior refresh for the next iteration
        -> existing hard assignment/update path
```

The artifact is not a replacement for `assignment.dat`, does not directly
update orientations, and is not a full joint posterior over in-plane angles
and shifts.

This reuses the useful framework introduced by commit `843c296`:

- the calibrated `dist -> exp(-dist)` likelihood path;
- sparse per-particle candidate bookkeeping;
- adaptive likelihood support selection; and
- direct-access, particle-keyed records for fast later reads.

It does not reuse the experimental posterior assignment literally. That code
draws one assignment from already evaluated candidates; it does not persist a
reusable posterior distribution.

### Runtime evidence motivating adaptive support

The first production test shows the failure mode clearly. At the transition
into posterior-guided refinement, the FSC(0.143) resolution worsened from
approximately `4.158 A` to `6.237 A`, while the mean in-plane displacement
rose above `43 degrees`. Across the subsequent iterations,
`POSTERIOR NEIGH NPEAKS` remained exactly `47.232` and the search fraction
remained exactly `1.889%`. The mapped support was therefore being reused
unchanged while the reference and assignments evolved.

This is not evidence that the subspace mapping is wrong. It demonstrates that
mapping a fixed support is insufficient when later likelihood evaluations move
the particle posterior. The stage-7 segmentation fault in the same run is a
separate robustness issue and must be debugged independently with bounds and
debug instrumentation.

## Scope and non-goals

### In scope

- 3D `refine=prob` as the posterior producer;
- later `refine=prob_neigh` stages as posterior-guided consumers;
- source-to-consumer resolution refinement, including larger `nspace`;
- state-aware sparse support;
- existing single-state, independent multi-volume, and docked workflows;
- serial, distributed, and OpenMP execution;
- direct-access artifact I/O.

### Non-goals for the first implementation

- a full posterior over all in-plane angles, shifts, and projection directions;
- a continuous SO(3) density estimator;
- exact posterior interpolation onto every target-grid point;
- aggregation when moving to coarser `nspace`;
- changes to the existing hard-assignment algorithm;
- changes to the non-posterior `state`, `geom`, `shc`, or `snhc` searches.

The first implementation should be a posterior over discrete `(state,
projection)` candidates with the in-plane and shift pose profiled by the
existing likelihood calculation.

## Probabilistic definition

Let `y_i` denote the observed image of particle `i`. A full pose would contain
the discrete state `s`, projection direction `r`, in-plane angle `phi`, and
shifts `Delta=(x,y)`. In Bayesian notation, the full posterior would be:

```text
p_i(s,r,phi,Delta | y_i)
    proportional to p(y_i | s,r,phi,Delta) * pi(s,r,phi,Delta)
```

where `pi` is the prior over poses. The present implementation does not
construct this full object. It collapses the nuisance pose
`eta=(phi,Delta)` to one representative value: for every discrete `(s,r)`,
the existing likelihood path finds or samples the relevant in-plane/shift
solution and returns its calibrated noise-normalized distance. Writing that
result as `d_i(s,r) = d_i(s,r,hat(eta_i(s,r)))`, the profile likelihood is:

```text
L_i(s,r) = exp(-d_i(s,r))
```

Thus, the word *posterior* means a normalized probability distribution over
the discrete state/projection hypotheses conditional on particle image `y_i`,
after nuisance pose has been profiled. With a uniform state/projection prior,
define:

```text
dmin_i = min_(s,r) d_i(s,r)
w_i(s,r) = exp(-(d_i(s,r) - dmin_i))
p_i(s,r) = w_i(s,r) / sum_(s,r) w_i(s,r)
```

The minimum shift is numerical stabilization only; it must not alter the
normalized probabilities. If a non-uniform reference prior is introduced
later, it belongs explicitly in the numerator before normalization.

This defines a proper normalized discrete posterior over the evaluated
profile-score candidates, provided that the producer evaluates the complete
dense candidate table and normalizes over all valid candidates. It is a
posterior over the reduced hypothesis space `(s,r)` under the calibrated
profile-score model; it is not a full posterior over `(s,r,phi,Delta)` because
the nuisance variables were collapsed rather than integrated out. In
particular, it does not claim uncertainty over the in-plane angle or shifts.

There is one implementation qualification. In strict statistical language,
profiling means choosing the maximizing nuisance value,
`hat(eta)=argmax_eta p(y_i|s,r,eta)`. The current `gen_likelihood_val` path can
select a likelihood-weighted in-plane candidate from the evaluated set rather
than always taking that maximum. In that case `L_i(s,r)` is a calibrated,
stochastic plug-in/profile score. The normalization is still mathematically
well-defined and useful for search-space restriction, but it should not be
overstated as an exact marginal posterior or as an exact maximum-profile
posterior.

The distinction from a marginal posterior is important. A marginal posterior
would be:

```text
p_i(s,r | y_i) = sum_or_integral_(phi,Delta)
                  p_i(s,r,phi,Delta | y_i)
```

whereas the implemented profile posterior is proportional to the likelihood
at one representative nuisance solution. Collapsing nuisance variables this
way is cheaper and is useful for deciding which projection directions deserve
later evaluation, but it can
under-represent a direction whose probability is spread over several
in-plane/shift alternatives. That limitation is deliberate and documented.

The sparse artifact is also not itself a posterior on a changed target grid.
It stores a truncated source-grid posterior that guides target candidate
generation. The consumer remaps its support geometrically, reevaluates the
current likelihood on the target grid, and lets that target likelihood drive
the subsequent assignment. Therefore, the source posterior supplies a
probabilistic search restriction; it does not claim that source probabilities
were exactly interpolated or conserved on the target grid.

### Profiled versus marginal support

The first artifact must identify its support mode as `profile`. The retained
in-plane index, shift, and distance describe the selected/profiled nuisance
pose associated with each projection candidate. The consumer remains
responsible for searching or refining nuisance variables.

Calling this artifact a full pose posterior would be incorrect. A future
`marginal` producer could replace the profile likelihood with a log-sum-exp or
equivalent integration over in-plane/shift alternatives, but that is outside
this refactoring.

### Forward-looking extension: Rao–Blackwellized pose posterior

A useful intermediate step toward a full pose posterior would retain the
current sparse posterior over state and projection direction, while adding a
conditional in-plane posterior only for the retained projection directions.
This avoids expanding the expensive in-plane dimension over the complete
angular search space.

For a retained direction `(s,r)`, let `Phi_i(s,r)` be the bounded set of
evaluated in-plane angles and let `d_i(s,r,phi)` be their calibrated
likelihood distances. The conditional in-plane distribution would be:

```text
q_i(phi | s,r,y_i) =
    exp(-d_i(s,r,phi)) * pi(phi)
    --------------------------------
    sum_(phi' in Phi_i(s,r)) exp(-d_i(s,r,phi')) * pi(phi')
```

The retained joint approximation would then factor as:

```text
q_i(s,r,phi | y_i)
    = p_i(s,r | y_i) * q_i(phi | s,r,y_i)
```

where `p_i(s,r | y_i)` is the existing sparse posterior. This is
Rao–Blackwellized in the practical sense that the conditional in-plane
uncertainty is represented explicitly for important state/projection
hypotheses, while the much larger set of low-probability directions remains
truncated.

Shifts would deliberately remain outside this model. For each retained
`(s,r,phi)`, the existing continuous shift search would continue to provide a
point estimate `Delta_hat_i(s,r,phi)`. The result would therefore be a
posterior over `(state, projection, in-plane angle)` conditional on profiled
shifts, not a full posterior over shifts.

The artifact can preserve the current fast direct-access design by using
bounded variable-length rows: each particle record has fixed maximum storage,
an explicit number of used projection entries, and an explicit number of used
conditional in-plane entries (or flattened pose entries). Unused capacity is
ignored. The metadata and particle-keyed records remain unchanged in spirit,
and one private reader can still be used per concurrent OpenMP worker.

The consumer could use the conditional in-plane distribution to seed or
prioritize target-grid evaluations. It should still reevaluate the target
likelihood after angular remapping; source-grid in-plane probabilities should
not be copied blindly because the in-plane coordinate frame can change under
projection-direction remapping.

This extension would require a likelihood API that returns the calibrated
in-plane score vector for a selected reference, rather than only the one
in-plane candidate currently returned by `gen_likelihood_val`. It would not
require enumerating a shift grid. Validation should check conditional-row
normalization, retained projection support mass, agreement with brute-force
in-plane evaluation on small cases, artifact size bounds, and unchanged
fallback behavior.

## Adaptive posterior lifecycle

The sparse artifact is a posterior for the model and reference volume used by
the producer iteration. It is not a timeless truth about the particle. Later
iterations change the reconstructed reference, noise estimates, and sometimes
the state assignment, so the consumer must treat the previous artifact as a
warm-start proposal rather than as immutable support.

For particle `i` and iteration `t`, let `q_t(h)` be the retained support
weights for hypothesis `h=(state, projection)`, and let `L_t(h)` be the current
calibrated likelihood obtained from the existing `dist -> likelihood`
transformation. The next support is built from:

```text
C_t = support(q_t) union exploration_t

q_(t+1)(h) proportional to
    L_t(h)^beta * ((1-alpha) * q_t(h) + alpha * pi_explore_t(h))
```

followed by normalization and bounded support selection. The previous support
is therefore a proposal/continuity term, not a second copy of the same image
evidence. Naively multiplying likelihoods from successive iterations would
double-count the particle image while the reference model is changing.

### Bounded exploration is mandatory

Updating weights only inside the old support cannot recover a direction that
was previously discarded. Every posterior-guided iteration must therefore
evaluate a small, deterministic exploration set outside the retained support.
The exploration set must be bounded and should be generated in the same coarse
subspace used by the posterior mapping. It may include:

- a small hash- or sequence-based sample from currently unrepresented coarse
  cells;
- the particle's current orientation cell and a bounded local state-aware
  probe; and
- additional cells requested when the support is flat, low-mass, or unstable.

The primary exploration structure is now fixed: for every retained coarse
projection cell, include that cell plus the three nearest coarse projection
cells in the `nspace_sub` angular search space. The builder constructs this
symmetry-aware projection-neighbor graph once. Each selected coarse cell is
then expanded through the existing `full2sub` mapping onto the active finer
grid. Thus a posterior row covering 25 coarse cells maps to approximately four
times the finer-grid multiplicity, rather than to an angular-radius scan of all
`nspace` directions.

If a particle's current projection cell is not covered by that mapped ring, or
if the row is empty, the consumer performs a conditional global probe using
the ordinary mode that would have been selected without
`prob_neigh_mode=posterior` (`state` for the single-state/independent path and
`geom` for docked, geometric, and input-orientation workflows). This probe is
restricted to the affected particles and its candidates are merged with the
mapped posterior candidates. It is a bounded recovery mechanism, not a
full-search fallback for the valid posterior path.

Exploration is part of `prob_neigh_mode=posterior`; it is not a silent switch
to a full search and should not be reported as fallback. The default budget
must remain small enough that the total candidate count scales with retained
support, the coarse-neighbor factor, and the grid-resolution ratio, not with
`nspace`.

Useful refresh triggers include a large increase in the best current distance,
a large change in the selected orientation or in-plane solution, low retained
support mass, a flat support, and failure of exploration candidates to overlap
the previous support. Thresholds should be named policy constants and exposed
in diagnostics rather than embedded as unexplained literals.

### Refresh timing and ownership

The lifecycle is:

```text
seed: final dense refine=prob merge
    -> write generation G0

posterior iteration t:
    read generation Gt
    map retained support and add bounded exploration
    evaluate current likelihoods
    merge distributed candidate tables
    refresh and atomically publish generation G(t+1)
    reconstruct/update the reference volume

next iteration:
    consume generation G(t+1)
```

The refresh should occur after the distributed neighborhood tables have been
merged, because that is the first point where the complete current candidate
evidence is available. It should be published before the next iteration reads
the artifact. The current iteration must continue using one pinned generation;
workers must never observe a file while it is being rewritten.

The current implementation performs this refresh from the merged candidate
store using the calibrated `exp(-dist)` likelihood transform, normalizes over
the evaluated candidate set, retains the bounded top support, and publishes a
temporary direct-access data file plus metadata before promotion. It therefore
updates the support from later evidence without multiplying likelihoods across
iterations. The explicit `(1-alpha) q_t + alpha pi_explore_t` proposal mixture
and generation/update-epoch fields are still follow-on hardening work.

The producer and consumer must remain in the particle-domain side of the
workflow. Volume assembly may change the reference used by the next iteration,
but it should not become responsible for posterior serialization.

### Adaptive artifact metadata

Each published generation should record:

```text
posterior generation and producer iteration
reference/update epoch used for the likelihoods
source grid and coarse-subspace mapping contract
number of posterior-guided and exploration evaluations
retained support mass and entropy summaries
refresh reason and trigger flags
```

The direct-access record format remains particle-keyed and bounded. Publication
should use a temporary generation-specific file followed by an atomic rename,
with metadata published consistently with the binary records. A consumer that
detects an epoch or generation mismatch must not mix rows from different
artifacts.

### Sparse support selection

The artifact stores a bounded support rather than all `nspace` candidates. The
support is selected from the complete dense table using the existing adaptive
peak/FDR logic, then ordered by posterior weight:

```text
K_i = bounded adaptive support size
S_i = top K_i valid candidates
support_mass_i = sum_(r in S_i) p_i(r)
```

`POSTERIOR3D_FDR_Q=0.25` remains a false-discovery-rate target for identifying
candidate peaks. It does not mean that 25% of projection directions are
retained.

The artifact must preserve both the normalized retained weights and
`support_mass`. Since truncation generally makes `support_mass < 1`, the
consumer must not silently treat the retained weights as a complete posterior.
It may either condition on the retained support or use the missing mass to
request bounded exploration/expansion. This is a posterior refresh signal,
not by itself a reason to switch to another neighborhood mode.

Minimum and maximum support sizes should remain enforced. Flat or ambiguous
particles should receive a broader bounded exploration set instead of an
overconfident small support.

## Resolution changes

The source posterior is defined on the source angular grid, but its retained
entries must be stored with Euler coordinates and state identity, not only
source reference indices. A source projection index has no meaning on a target
grid with a different `nspace`.

The builder already constructs a `simple_srchspace_map` between the consumer
fine grid and a coarse subspace. The artifact metadata distinguishes the full
source grid (`source_nspace`) from its source coarse grid
(`source_nspace_sub`). Each retained record stores its full source projection
identity and source Euler coordinates. For a changed grid, the consumer maps
the source Euler direction to the nearest target coarse representative, then
uses the target map:

```text
coarse(t) = full2sub_map(t)
```

For a retained posterior projection `r`, the consumer evaluates exactly the
fine-grid directions assigned to the mapped target coarse cell:

```text
T(r) = { t : full2sub_map(t) = r }
```

When source and target full grids are identical, the target builder map is used
directly. When the target grid is finer, the source Euler coordinates are
mapped once per source projection and cached, so a retained source direction
maps to the appropriate target coarse cell and then to roughly the expected
resolution-ratio number of fine directions. Twenty-five retained directions
therefore produce roughly fifty target candidates in a 2:1 refinement, with
deduplication across states and support entries. The target likelihood remains
authoritative for the final ordering.

This mapping is preferable to expanding an angular-radius ball around every
retained direction: it preserves the source-to-target cell ownership already
used by the builder and makes the candidate count scale with the resolution
ratio rather than with the area of overlapping angular neighborhoods. The
source and target coarse-subspace resolutions need not be identical; the
source Euler coordinates establish the cross-grid correspondence. The
consumer target `nspace` may be finer.

Only refinement to a finer angular grid is supported. Coarsening and posterior
aggregation are explicitly out of scope.

### Mapping ownership and performance

The full-to-subspace map is a property of the source and target Euler spaces,
not of an individual particle. It must be constructed once per refinement run
and materialized as compact per-subspace target-index lists. Workers then only
look up the lists for the particle's retained posterior directions. The
consumer must not independently scan the full target space for every posterior
row.

The mapping cost should scale with the number of distinct retained source
directions times the target coarse-space size, not with:

```text
number of particles * source nspace * target nspace
```

This is the principal performance requirement of the posterior implementation.

## Artifact contract

The existing fixed-record direct-access pattern should be retained. It is
appropriate for random particle reads and avoids text parsing or per-particle
file searches in the consumer hot path.

### Header metadata

The versioned metadata must include:

```text
magic and format version
source refine stage and final source iteration
update-epoch identifier
represented particle ids and count
source nstates, full source nspace, and source nspace_sub
point-group/symmetry identifier
box, crop, sampling, and relevant objective/noise metadata
likelihood transform = profile_exp_minus_dist
support mode = profile
FDR/support parameters and Kmax
mapping version and source-to-subspace map contract
direct-record schema and record capacity
```

### Per-particle record

Each record should contain:

```text
particle id
nsel and support_mass
state and source projection identity
source Euler angles
profile distance and normalized posterior weight
cumulative retained weight/rank
best in-plane index
x/y shift and has-shift flag
optional target-mapping diagnostics
```

The source Euler angles are authoritative for resolution remapping. The source
projection index is retained for provenance and same-grid fast paths; it is not
treated as a target-grid index after an `nspace` change.

The reader must expose a particle-row operation. It must not require scanning
all records to locate one particle.

## Controller timing

The final dense probability iteration is the seed producer point:

```text
final refine=prob iteration
    -> merge dense probability tables
    -> extract and normalize sparse profiled posterior
    -> publish posterior generation G0
    -> write existing hard assignment
```

The seed artifact must be written only after the final dense table has been
merged. It must not be generated from an earlier iteration. Once posterior
neighborhood refinement starts, the merged current `prob_neigh` candidate
table becomes a second producer point:

```text
posterior prob_neigh iteration
    -> merge current sparse candidate/evidence tables
    -> rebuild bounded per-particle posterior rows
    -> publish generation G(t+1) atomically
```

The refresh must use the current calibrated likelihood values, not raw
objective scores and not a cumulative product of likelihoods from previous
iterations. The next iteration consumes the newly published generation.

The controller should create the seed artifact only when a later
`prob_neigh_mode=posterior` stage is actually requested. It should then keep
the refresh path active for each posterior-guided iteration until the stage
ends. The artifact generation and reference/update epoch must be carried into
the child commands or equivalent workflow metadata so a stale artifact cannot
silently cross a stage or state-split boundary.

For docked multi-volume workflows, an explicit posterior-guided request should
override the default `geom` mode. The artifact must carry the update epoch; it
must not cross a docked state-split/update boundary without a new producer.

## Consumer behavior

The consumer should perform the following per-particle operation:

```text
read and pin sparse posterior generation
validate particle/state/epoch coverage
map retained Euler directions to bounded target candidates
add the coarse three-neighbor exploration ring
conditionally probe affected particles with the ordinary search mode
deduplicate target references
evaluate current target likelihood and nuisance pose
run the existing neighborhood assignment/update logic
merge current candidate evidence
publish the next posterior generation
```

The posterior is guidance. It should not bypass the existing matcher,
assignment, shift, or sigma update paths.

The current target likelihood remains authoritative. The previous posterior
weight may order candidates or regularize the refresh, but it must not be copied
onto remapped target directions as if it were an exact target-grid posterior.

The mode must retain multiple disconnected posterior modes. A purely local
geometric expansion around the highest-weight direction is not equivalent to
posterior guidance.

## Fallback policy

Adaptive exploration and fallback are separate mechanisms. A valid mapped
posterior row must remain in `prob_neigh_mode=posterior`; the normal finer-grid
resolution change must never trigger fallback. Exploration is the bounded
recovery mechanism for uncertain or drifting support.

Fallback must remain bounded and must be reported per particle and globally.
It must never silently become a full `nspace` search.

- single-state and `multivol_mode=independent`: fall back to the bounded
  `prob_neigh_mode=state` search;
- docked, geometric, and input-orientation workflows: fall back to the bounded
  `prob_neigh_mode=geom` search;
- missing, stale, malformed, or insufficiently covered rows: use the same
  corresponding fallback for those particles only;
- low retained mass, flat posterior, excessive remapping span, or a large
  evidence drift: request bounded posterior exploration before falling back;
- only an invalid artifact, an invalid row, or exhausted bounded exploration
  may enter the corresponding ordinary-mode fallback.

The convergence report should include:

```text
posterior-guided particles
posterior generations consumed and published
exploration particles and candidate counts
exploration hits outside the previous support
fallback particles by mode/reason
mean/quantiles of nsel
mean/quantiles of support_mass
mean/quantiles of posterior entropy
mapped target candidates before/after deduplication
mapping-distance and coverage diagnostics
```

The dynamic support count must be reported through the existing `npeaks`
convergence field, preserving the useful 2D diagnostic without claiming that
it is a global fixed neighborhood size.

## Parallel and I/O requirements

- retain one private candidate/mapping workspace per OpenMP thread;
- use one direct-access reader unit per concurrent reader, as in the existing
  fast-record pattern;
- avoid allocation inside the per-particle candidate loop;
- deduplicate with private markers or generation stamps rather than clearing
  an `nrefs`-sized logical array for every particle;
- construct mapping/index data before the hot parallel loop;
- keep candidate counts bounded so posterior guidance cannot create a hidden
  full-grid workload; and
- preserve the current distributed partitioning and canonical particle order;
- pin one posterior generation for the duration of an iteration; and
- publish refreshed records through a generation-specific temporary file and
  atomic replacement, never by mutating records visible to active readers.

## Removal and replacement inventory

The implementation removes the current `prior` concept rather than leave
two competing 3D support mechanisms in the code base:

- replace `simple_prob_prior3D.f90` and its `prior3d_writer`/
  `prior3d_reader` types with posterior-support I/O;
- remove `PRIOR3D_*` remapping constants and the worker-local `prior3d_ws`;
- replace all parameter, UI, commander, builder, controller, strategy, and
  convergence checks for `prob_neigh_mode=prior` with `posterior`;
- replace `write_prior3D` with the final dense-`prob` posterior writer;
- replace the current full source-to-target remap initialization in
  `simple_eul_prob_tab_neigh.f90` with builder-owned subspace-map lookups;
- preserve the existing `npeaks` field in `assgn_map` and the convergence
  statistics, but remove any prior-specific naming or separate posterior peak
  variable; and
- keep the existing `state`/`geom` fallback branches as ordinary modes, not as
  hidden variants of the removed `prior` mode.

## Refactoring sequence and implementation status

1. Define the posterior artifact schema and replace the current `prior3D`
   writer/reader with posterior-support I/O. **Complete.**
2. Extract the final dense `refine=prob` table into normalized sparse profile
   posterior rows. **Complete.**
3. Correctly record retained mass, profile metadata, source Euler coordinates,
   and per-particle diagnostics. **Partially complete; adaptive diagnostics
   remain.**
4. Implement the source-to-target subspace mapper using the builder's
   `simple_srchspace_map` outputs and add coverage and duplicate diagnostics.
   **Complete for the current mapped path.**
5. Replace the current worker-local full remap initialization with particle-row
   reads and compact per-subspace target-index lists. **Complete.**
6. Feed mapped candidates into the existing neighborhood search while keeping
   target-grid likelihood evaluation authoritative. **Complete.**
7. Add bounded posterior exploration outside the retained support, using a
   named exploration budget and deterministic per-particle selection. **Complete
   for the coarse three-neighbor ring and conditional ordinary-mode probe.**
8. Rebuild posterior rows from the merged current neighborhood evidence using
   the calibrated likelihood transform. Treat the previous artifact as a
   proposal/continuity term, not as cumulatively re-used evidence. **Current
   evidence refresh is complete; explicit proposal mixing remains.**
9. Add generation/epoch metadata, pinned per-iteration reads, atomic refresh
   publication, and controller timing for the seed and refresh producers.
   **Seed/refresh timing and temporary-file publication are complete; explicit
   generation/epoch metadata remains.**
10. Add convergence reporting for posterior generations, exploration hits,
    support mass, entropy, refresh triggers, and bounded fallback coverage,
    while continuing to report candidate count through `npeaks`.
11. Add bounds/debug validation for the posterior consumer and finer-grid
    candidate store, then validate same-grid, finer-grid, refresh, recovery,
    distributed, OpenMP, and artifact-replacement behavior.
12. Remove the old `prior` implementation, all associated code paths, and the
    old artifact. **Complete; runtime safety validation remains.**

## Acceptance criteria

The refactoring is ready for implementation completion when:

- a dense `refine=prob` row normalizes to a valid source-grid posterior;
- truncation preserves explicit `support_mass` and never masquerades as full
  normalization;
- same-grid posterior-guided search reproduces the intended candidate support;
- finer-grid consumers use Euler/symmetry mapping rather than source indices;
- a retained source direction maps to all target directions assigned to its
  builder subspace cell;
- mapped target candidates are deduplicated without probability inflation;
- multiple separated modes remain eligible;
- a valid finer-grid mapping does not trigger fallback;
- bounded exploration can discover support outside the previous posterior;
- refreshed posterior rows respond to later current-likelihood evidence;
- missing or malformed support falls back only through the bounded
  state/geom rules;
- no worker performs a full source-to-target remap or target-grid scan per
  posterior row;
- direct-access reads are safe under distributed/OpenMP execution;
- the seed is written by the final dense `prob` iteration, and later
  generations are written only after merged current neighborhood evidence;
- readers never observe a partially replaced artifact; and
- convergence output reports the dynamic posterior-guided count through the
  existing `npeaks` statistics, exploration diagnostics, generation updates,
  and fallback coverage.

## Locked design decisions for implementation

- The external mode is `prob_neigh_mode=posterior`; `prior` is removed.
- The source artifact is a sparse profiled posterior, not a one-sample
  assignment and not a full in-plane posterior.
- The target candidate set is the union of fine-grid directions assigned to
  the coarse subspace cells selected by the sparse posterior.
- No `posterior_npeaks` parameter or separate reporting field is introduced.
  The existing `npeaks` value and convergence reporting are reused for the
  dynamic number of posterior-guided target directions.
- The target-stage likelihood determines final candidate ranking. Posterior
  weights may guide candidate ordering, but must not be copied onto target
  directions without a mass-preserving mapping.
- The final dense `prob` artifact is the seed generation; posterior-guided
  iterations refresh later generations from current evidence.
- Previous posterior weights are a tempered proposal/continuity term, not a
  cumulative product of repeated image likelihoods.
- Every posterior-guided iteration includes a bounded exploration budget so
  excluded directions can be recovered without a full search.
- A valid mapped support remains in posterior mode. Missing or malformed
  support uses the existing bounded state/geom fallback rules and never a full
  search.
