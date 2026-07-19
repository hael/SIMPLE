# Sparse Profiled Posterior for 3D Neighborhood Refinement

## Status

Implemented in the 3D path, pending runtime validation on production data.
This note replaces the former 3D `prob_neigh_mode=prior` implementation and
associated code. It does not change the 2D implementation.

The user-facing mode is `prob_neigh_mode=posterior`. The `prior` mode,
artifact, remapping workspace, controller branches, and reporting branches are
to be removed; no compatibility alias is required.

## Decision summary

The final dense `refine=prob` phase should produce a sparse, per-particle,
profiled posterior over state and projection direction. Later neighborhood
refinement stages should use that posterior to generate and prioritize a small
candidate set on their own angular grid, then evaluate the current likelihood
on those target candidates.

The artifact is therefore:

```text
final dense prob table
        -> sparse profiled posterior artifact
        -> target-grid candidate generation
        -> existing neighborhood likelihood evaluation
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
trigger a bounded fallback/expansion.

Minimum and maximum support sizes should remain enforced. Flat or ambiguous
particles should receive a broader bounded search instead of an overconfident
small support.

## Resolution changes

The source posterior is defined on the source angular grid, but its retained
entries must be stored with Euler coordinates and state identity, not only
source reference indices. A source projection index has no meaning on a target
grid with a different `nspace`.

For each retained source direction:

1. construct or access the target Euler space with the same symmetry;
2. locate a bounded set of nearby target projection directions using
   symmetry-aware angular distance;
3. deduplicate target references across source support entries; and
4. evaluate the current target-grid likelihood for the resulting candidates.

The target stencil should not be a fixed three-neighbor remap. For each
particle, use the process-provided angular distance estimate
`DIST BTW BEST ORIS (DEG)` and define:

```text
radius_i = 2 * DIST_BTW_BEST_ORIS_i
```

Select every target projection direction within `radius_i`, using
symmetry-aware angular distance. If the radius contains fewer than three
directions, add the nearest target directions until three are present. The
minimum-three rule is a lower bound, not the normal selection rule. The
convergence `AVG` of `DIST BTW BEST ORIS (DEG)` is only a report; the mapper
must use the corresponding per-particle estimate. If the internal angular
distance routine uses radians, convert the degree radius once at the API
boundary.

The mapped entries are a candidate support and the target likelihood remains
authoritative for the final ordering. If source weights are also used as
target prior factors, the mapping must be mass-preserving:

```text
target_prior(t) = sum_r source_mass(r) * K(t, r)
```

where `K` is a normalized angular kernel or an equivalent cell/barycentric
mapping. The mapping must account for duplicate target references and, where
needed, source/target angular cell area. Equal copying to the nearest three
directions is not a valid probability mapping.

The source posterior must not be used to assign zero probability to all target
directions outside the selected angular-radius neighborhoods. The omitted
source mass and mapping coverage must be retained as diagnostics and used to
request a bounded expansion or fallback when appropriate.

Only refinement to a finer angular grid is supported. Coarsening and posterior
aggregation are explicitly out of scope.

### Mapping ownership and performance

The angular-distance index is a property of the source and target Euler
spaces, not of an individual particle. It must be constructed once per
refinement run, or cached by `(pgrp, source grid, target grid, mapping
version)`. The per-particle radius is then applied as a query against that
index. Workers must not independently scan the full source artifact and
perform a full target-space nearest-neighbor search.

The angular search should cover only the per-particle radius neighborhood. Its
cost must scale with retained support and the selected angular neighborhood,
not with:

```text
number of particles * source nspace * target nspace
```

This is the principal performance requirement that the current prior
implementation does not satisfy.

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
source nstates, nspace, and nspace_sub
point-group/symmetry identifier
box, crop, sampling, and relevant objective/noise metadata
likelihood transform = profile_exp_minus_dist
support mode = profile
FDR/support parameters and Kmax
mapping version and angular-radius rule
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
projection index is retained for provenance and same-grid fast paths only.

The reader must expose a particle-row operation. It must not require scanning
all records to locate one particle.

## Controller timing

The final dense probability iteration is the only producer point:

```text
final refine=prob iteration
    -> merge dense probability tables
    -> extract and normalize sparse profiled posterior
    -> write posterior artifact
    -> write existing hard assignment
```

The artifact must be written only after the final dense table has been merged.
It must not be regenerated from a sparse `prob_neigh` table or from an earlier
iteration.

The following neighborhood stages consume the artifact. The controller should
write it only when a later `prob_neigh_mode=posterior` stage is actually
requested.

For docked multi-volume workflows, an explicit posterior-guided request should
override the default `geom` mode. The artifact must carry the update epoch; it
must not cross a docked state-split/update boundary without a new producer.

## Consumer behavior

The consumer should perform the following per-particle operation:

```text
read sparse posterior row
validate particle/state/epoch coverage
map retained Euler directions to bounded target candidates
deduplicate target references
evaluate current target likelihood and nuisance pose
run the existing neighborhood assignment/update logic
```

The posterior is guidance. It should not bypass the existing matcher,
assignment, shift, or sigma update paths.

The mode must retain multiple disconnected posterior modes. A purely local
geometric expansion around the highest-weight direction is not equivalent to
posterior guidance.

## Fallback policy

Fallback must remain bounded and must be reported per particle and globally.
It must never silently become a full `nspace` search.

- single-state and `multivol_mode=independent`: fall back to the bounded
  `prob_neigh_mode=state` search;
- `multivol_mode=docked`: fall back to the bounded `prob_neigh_mode=geom`
  search;
- missing, stale, malformed, or insufficiently covered rows: use the same
  corresponding fallback for those particles only;
- low retained mass, flat posterior, or excessive remapping span: request a
  bounded expansion before falling back.

The convergence report should include:

```text
posterior-guided particles
fallback particles by mode/reason
mean/quantiles of nsel
mean/quantiles of support_mass
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
- preserve the current distributed partitioning and canonical particle order.

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
  `simple_eul_prob_tab_neigh.f90` with per-particle angular-radius mapping;
- preserve the existing `npeaks` field in `assgn_map` and the convergence
  statistics, but remove any prior-specific naming or separate posterior peak
  variable; and
- keep the existing `state`/`geom` fallback branches as ordinary modes, not as
  hidden variants of the removed `prior` mode.

## Refactoring sequence and implementation status

1. Define the posterior artifact schema and replace the current `prior3D`
   writer/reader with posterior-support I/O.
2. Extract the final dense `refine=prob` table into normalized sparse profile
   posterior rows.
3. Correctly record retained mass, profile metadata, source Euler coordinates,
   and per-particle diagnostics.
4. Implement the per-particle angular-radius mapper using
   `2 * DIST BTW BEST ORIS (DEG)` and a minimum of three target directions.
   Add coverage and duplicate diagnostics.
5. Replace the current worker-local full remap initialization with particle-row
   reads and the angular-radius mapper.
6. Feed mapped candidates into the existing neighborhood search while keeping
   target-grid likelihood evaluation authoritative.
7. Add controller timing, update-epoch validation, fallback reporting, and
   convergence reporting through the existing `npeaks` field and statistics.
8. Remove the old `prior` implementation, all associated code paths, and the
   old artifact after same-grid, finer-grid, and fallback behavior are
   validated. This removal and replacement is now complete; production-data
   validation remains.

## Acceptance criteria

The refactoring is ready for implementation completion when:

- a dense `refine=prob` row normalizes to a valid source-grid posterior;
- truncation preserves explicit `support_mass` and never masquerades as full
  normalization;
- same-grid posterior-guided search reproduces the intended candidate support;
- finer-grid consumers use Euler/symmetry mapping rather than source indices;
- target support contains every direction within twice the per-particle
  `DIST BTW BEST ORIS (DEG)` estimate, with at least three directions;
- mapped target candidates are deduplicated without probability inflation;
- multiple separated modes remain eligible;
- missing or weak support falls back to bounded state/geom search;
- no worker performs a full source-to-target remap per particle or per batch;
- direct-access reads are safe under distributed/OpenMP execution;
- the final `prob` iteration alone writes the artifact; and
- convergence output reports the dynamic posterior-guided count through the
  existing `npeaks` statistics and reports fallback coverage.

## Locked design decisions for implementation

- The external mode is `prob_neigh_mode=posterior`; `prior` is removed.
- The source artifact is a sparse profiled posterior, not a one-sample
  assignment and not a full in-plane posterior.
- The target candidate set is all directions within twice the per-particle
  `DIST BTW BEST ORIS (DEG)` estimate, with a minimum of three directions.
- No `posterior_npeaks` parameter or separate reporting field is introduced.
  The existing `npeaks` value and convergence reporting are reused for the
  dynamic number of posterior-guided target directions.
- The target-stage likelihood determines final candidate ranking. Posterior
  weights may guide candidate ordering, but must not be copied onto target
  directions without a mass-preserving mapping.
- Missing or weak support uses the existing bounded state/geom fallback rules;
  it must not trigger a full search.
