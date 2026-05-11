# Continuous SO(3) Probability Model for PFTC Assignment

**Date:** May 10, 2026  
**Status:** Design note / first-draft implementation plan

## Purpose

The current `eul_prob` classes define a probabilistic model over discrete reference points on SO(3). In practice, they build particle-reference score tables over an Euler-space grid, optionally sparsify the grid through neighborhoods, then perform global stochastic assignment from those evaluated entries.

This note explores a continuous alternative: represent the particle-local orientation uncertainty as a proposal distribution on SO(3), sample a small set of trial orientations from that distribution, extract polar central sections from the current reference volume on the fly, evaluate those sections with the existing PFTC objective, and feed the evaluated trials into a global sparse assignment step.

The immediate goal is not to analytically integrate over SO(3). The recommended first goal is a continuous sparse candidate generator whose outputs are evaluated and assigned using the same objective semantics SIMPLE already trusts.

## Current Code Shape

The relevant existing implementation is concentrated in:

- `src/main/simple_eul_prob_tab.f90`
  - dense probability-table construction over active discrete projection references
  - `fill_tab`
  - `ref_assign`
  - `ref_score_tab`
  - `read_tab_to_glob`
  - `write_assignment`
- `src/main/simple_eul_prob_tab_neigh.f90`
  - sparse neighborhood probability-table construction
  - `fill_tab_neigh`
  - `record_sparse_eval`
  - `write_tab_neigh`
  - `read_sparse_tab_to_glob`
  - `ref_assign_neigh`
- `src/main/simple_eul_prob_tab_utils.f90`
  - distance/correlation conversion
  - bounded and power sampling utilities
  - seeded shift materialization
- `src/main/pftc/simple_polarft_corr.f90`
  - exact PFTC objective evaluation once a reference PFT exists
- `src/main/image/simple_projector_pft_batch.f90`
  - batched extraction of discrete reference polar sections
- `src/main/image/simple_projector_pft.f90`
  - single-orientation polar section extraction adapters
- `src/main/strategies/search/simple_strategy3D_prob.f90`
  - consumption of probability-table assignments during search/update
- `src/main/strategies/search/simple_strategy3D_matcher.f90`
  - partition-level probabilistic assignment use and matcher integration

The dense class stores candidates as `loc_tab(ri,i)`, where `ri` maps to a fixed `(state, iproj)` reference. The neighborhood class keeps this same reference identity but stores only touched entries using `eval_touched_refs`.

A continuous candidate cannot be represented fully by only `ri`, because its orientation may not coincide with any fixed `iproj`. It may still need a nearest `iproj` for legacy bookkeeping, but its actual orientation must be carried explicitly.

## Proposed Model

For particle `i`, define a proposal distribution on SO(3):

```text
q_i(R, s) = mixture of local orientation kernels
```

where:

- `R` is a continuous orientation.
- `s` is the state.
- mixture centers come from previous assignments, coarse discrete candidates, sparse neighborhood representatives, or random exploration.
- mixture widths describe current uncertainty in projection direction and optionally in-plane angle.

The continuous proposal is used to generate a small candidate set:

```text
C_i = { (state, R_sample, proposal_weight) }
```

Each sampled candidate is then evaluated exactly:

1. Extract a polar central section from the current state volume at `R_sample`.
2. Evaluate all in-plane rotations using the existing PFTC objective.
3. Optionally sample or maximize the in-plane rotation using the existing probabilistic in-plane policy.
4. Optionally run existing continuous shift refinement for top candidates.
5. Store the evaluated candidate in a sparse table.
6. Perform global assignment over sparse evaluated candidates.

This makes the continuous distribution a candidate-generation model first. A later implementation can use proposal-density correction to approximate a true posterior over SO(3).

## Proposal Parameterization

Avoid Euler-space Gaussians. They are singular and geometry-dependent.

The recommended first parameterization is a local tangent-plane distribution around one or more orientation centers:

```text
R_sample = exp([delta]x) * R_center
delta    = a*u + b*v
```

where `u` and `v` are tangent directions for the local projection direction. In-plane rotation can remain discrete through the current PFTC all-rotation evaluation.

Reasonable first kernels:

- tangent-plane Gaussian around previous orientation
- tangent-plane Gaussian around top coarse discrete references
- wider heavy-tailed tangent proposal for exploration
- state-specific mixtures built from top candidates per state

Longer-term kernels could include:

- von Mises-Fisher or Kent-like distributions on projection direction
- Bingham-like antipodal kernels where projection-direction sign or mirror symmetry matters
- a factorized `S2 projection direction x S1 in-plane` model
- a full SO(3) kernel if an implementation need emerges

For a first implementation, keep in-plane rotation discrete and let PFTC select or probabilistically sample it.

## Assignment Semantics

There are two possible interpretations.

### Proposal Model First

The continuous distribution only proposes candidates. Final assignment uses evaluated distances exactly as current sparse assignment does.

This is the recommended first path because:

- it preserves SIMPLE's current hard-assignment workflow
- it avoids requiring correct SO(3) quadrature weights immediately
- it is easier to compare directly against `eul_prob_tab_neigh`
- it keeps global sparse assignment conceptually intact

### True Posterior Approximation Later

Candidate weights can later include proposal-density correction:

```text
weight(R) proportional to likelihood(X_i | R, s) * prior(R, s) / q_i(R, s)
```

This would move the code closer to importance sampling on SO(3). It is more statistically complete, but it introduces new responsibilities:

- evaluate or approximate `q_i(R,s)`
- define orientation/state priors
- handle repeated or nearby samples
- normalize candidate weights without biasing states or proposal centers

This should not be the first implementation.

## Sparse Candidate Representation

The current `ptcl_ref` structure stores:

```text
pind, istate, iproj, inpl, dist, x, y, has_sh
```

A continuous candidate needs at least:

```text
pind
istate
nearest_iproj or source_iproj
continuous orientation, preferably e1/e2/e3 or a rotation representation
inpl
dist
x, y, has_sh
proposal metadata, optional
```

There are two implementation options.

### Option A: Sidecar Orientation Arrays

Keep `ptcl_ref` for most existing assignment data and add parallel arrays in a new continuous sparse class:

```text
loc_tab-like ptcl_ref candidate values
candidate_eulers(3, max_candidates, nptcls)
candidate_weights(max_candidates, nptcls)
candidate_source(max_candidates, nptcls)
```

This is less invasive and is preferred for a prototype.

### Option B: New Candidate Type

Define a sibling type, for example:

```fortran
type :: ptcl_cont_ref
    integer :: pind = 0
    integer :: istate = 0
    integer :: iproj_nearest = 0
    integer :: inpl = 0
    real    :: eul(3) = 0.
    real    :: dist = huge(1.0)
    real    :: x = 0.
    real    :: y = 0.
    real    :: proposal_logw = 0.
    logical :: has_sh = .false.
end type
```

This is cleaner long term but touches more I/O and assignment code.

## First-Draft Algorithm

For each particle:

1. Read previous orientation and state from `spproj_field`.
2. Estimate a seed shift using the previous assignment, as current `fill_tab` and `fill_tab_neigh` already do.
3. Build proposal centers:
   - previous orientation
   - top coarse discrete references from a cheap sweep or subspace sweep
   - optional random global orientations for exploration
   - optional state-specific centers
4. For each proposal center:
   - sample `m` tangent-plane perturbations
   - convert each perturbation to a full orientation
   - reject or fold through symmetry if needed
5. For each sampled orientation:
   - extract a temporary polar central section from the current state volume
   - evaluate all in-plane rotations using existing PFTC scoring
   - sample or select in-plane rotation
   - store distance and orientation in a sparse candidate list
6. Refine shifts for the top candidates per particle/state using existing `pftc_shsrch_grad` logic when possible.
7. Write sparse candidate files for partitions.
8. Read all sparse candidate files into a global object.
9. Perform global sparse assignment over evaluated continuous candidates.
10. Write an assignment map carrying the accepted continuous orientation, in-plane rotation, state, shift, and distance.

The first version should only use continuous candidates after a coarse discrete or sparse neighborhood pass. That preserves global coverage while testing whether continuous SO(3) sampling can reduce the number of exact section evaluations.

## Integration Strategy

The least risky implementation path is a sibling class:

```text
simple_eul_prob_tab_cont
```

or, if the first version is explicitly neighborhood-driven:

```text
simple_eul_prob_tab_neigh_cont
```

Recommended ownership:

- `eul_prob_tab_cont` owns sampled candidate lists, candidate I/O, and global assignment.
- A small SO(3) proposal helper owns tangent sampling and mixture-center generation.
- `projector`/PFTC remain numerical engines and should not own probabilistic policy.
- Existing matcher/strategy code should consume the final assignment map without caring whether the candidate came from a grid point or a continuous sample.

The new class can reuse:

- seed-shift handling from `eul_prob_tab`
- sparse touched-candidate concepts from `eul_prob_tab_neigh`
- sparse graph assignment structure from `ref_assign_neigh`
- distance conversion utilities from `simple_eul_prob_tab_utils`

It should avoid reusing:

- `ri` as the sole candidate identity
- `nrefs` as the number of candidate alternatives
- global arrays shaped as `(nrefs,nptcls)` for continuous candidate data

## Normalization and Bias Control

Continuous candidates create a new normalization issue: each particle may have a different number and distribution of sampled orientations.

For the first proposal-model implementation:

1. Convert raw objective to distance exactly as now.
2. Normalize only over candidates evaluated for the same particle.
3. Keep state assignment from being biased by unequal candidate counts across states.
4. If sampling per state is unequal, assign state priors or normalize the best candidates per state before global assignment.

For later importance-sampling semantics:

```text
score = objective_distance - log_prior + log_proposal
```

or an equivalent likelihood-weight form could be introduced. This should be deferred until the proposal-model path is validated.

## First-Draft Implementation Plan

### Draft 0: Continuous Candidate Data Model Only

Purpose: make the assignment pipeline capable of carrying continuous orientations without changing how candidates are generated.

Implementation:

- Add a prototype continuous sparse candidate type or sidecar orientation arrays.
- Generate candidates at exact existing grid orientations.
- Write/read sparse candidate files.
- Assign globally.
- Confirm output assignments match or closely reproduce the discrete sparse path.

Expected result:

- The continuous assignment representation should be able to emulate current sparse grid assignment.

### Draft 1: Tangent Samples Around Previous Orientation

Purpose: test on-the-fly sampled section extraction without needing a coarse proposal mixture.

Implementation:

- For each particle/state, sample a small number of orientations around the previous orientation.
- Extract polar sections on-the-fly.
- Evaluate with PFTC.
- Assign from the sampled candidate list.

Expected result:

- On synthetic data with small pose perturbations, sampled candidates should improve or preserve assignment quality.
- On broad uncertainty cases, performance may be poor; this is expected because the proposal is too local.

### Draft 2: Coarse-Discrete Mixture Proposal

Purpose: combine global discrete coverage with continuous local exploration.

Implementation:

- Use a cheap discrete sweep or existing neighborhood/subspace sweep to identify `K` centers per particle/state.
- Around each center, sample `m` tangent perturbations.
- Keep total exact section evaluations bounded by `K*m`.
- Optionally include previous orientation as a persistent mixture center.

Expected result:

- For the same number of evaluated candidates, continuous samples should outperform fixed neighborhood grid candidates when the true orientation lies between grid points.
- This draft is the first serious comparison against `eul_prob_tab_neigh`.

### Draft 3: Sparse Global Assignment Over Continuous Candidates

Purpose: replace `ref -> particle` sparse graph assumptions with candidate-based sparse assignment.

Implementation:

- Build a candidate graph where each particle has `M_i` evaluated candidates.
- If global balancing by reference is still desired, define candidate groups by nearest projection reference or state.
- Otherwise perform per-particle stochastic assignment with optional state balancing.
- Preserve fallbacks for particles with no valid candidates.

Expected result:

- Assignments should remain robust when candidate sets differ per particle.
- The accepted orientation must be the continuous sampled orientation, not only the nearest grid reference.

### Draft 4: Optional Proposal-Density Correction

Purpose: move from candidate-generation semantics toward true continuous probabilistic inference.

Implementation:

- Track mixture component and proposal log-density for each sample.
- Define a prior over state/orientation.
- Adjust candidate weights by `prior / proposal`.
- Compare assignments with and without correction.

Expected result:

- Correction should matter only when proposal densities are highly nonuniform.
- If it destabilizes assignment, keep it disabled until the proposal model is better calibrated.

## Tests and Validation

### Unit-Level Tests

1. **Tangent sampler geometry**
   - Sample many perturbations around a known orientation.
   - Verify angular displacement distribution matches requested width.
   - Verify no Euler singularity behavior near poles.

2. **Symmetry folding**
   - For cyclic/dihedral/symmetric cases, verify equivalent orientations are handled consistently.
   - Confirm nearest-grid bookkeeping remains stable.

3. **Grid-emulation test**
   - Create continuous candidates exactly at grid orientations.
   - Verify assignment output matches the sparse discrete path within expected stochastic variation.

4. **Sparse candidate I/O**
   - Write and read continuous candidate files.
   - Verify particle ids, states, continuous Euler angles, in-plane indices, shifts, and distances survive round-trip.

5. **Candidate normalization**
   - Test equal and unequal candidate counts per state.
   - Verify assignment is not biased simply because one state has more samples.

### Synthetic Integration Tests

1. **Known orientation, between grid points**
   - Generate particles from known orientations deliberately placed between discrete projection directions.
   - Compare:
     - discrete `eul_prob`
     - sparse neighborhood `eul_prob`
     - continuous proposal samples
   - Measure angular error and objective score.

2. **Fixed evaluation budget**
   - Hold the number of exact PFTC evaluations constant.
   - Compare fixed grid neighborhood candidates versus continuous samples around coarse centers.

3. **Coarse grid plus continuous proposal**
   - Run a coarse global grid with continuous local sampling.
   - Compare against the normal dense grid.
   - Measure runtime, orientation accuracy, and reconstruction quality.

4. **Multimodal ambiguity**
   - Construct or select particles with ambiguous orientations.
   - Verify mixture proposals preserve multiple modes rather than collapsing immediately to previous assignment.

5. **State movement**
   - In a multi-state case, verify continuous proposal assignment can move particles between states without being dominated by sample-count imbalance.

### Performance Measurements

Measure:

- number of on-the-fly section extractions per particle
- cost per extracted continuous section
- cost per PFTC all-in-plane evaluation
- cost of shift refinement per accepted/top candidate
- total assignment-table I/O size
- global sparse assignment time
- accuracy/runtime tradeoff versus `eul_prob_tab_neigh`

The important comparison is:

```text
same or lower evaluation count, better orientation accuracy
```

not simply higher score after adding more candidates.

## Design Risks

1. **Proposal collapse**
   - If samples are drawn only near the previous orientation, the model may fail to escape bad assignments.
   - Mitigation: use mixture centers from coarse discrete candidates and include occasional broad exploration.

2. **Sample-count bias across states**
   - More samples in one state can create unfair assignment opportunities.
   - Mitigation: equalize candidate budgets per state or normalize state candidates carefully.

3. **Loss of global sparse assignment semantics**
   - Existing `ref_assign_neigh` uses active references as graph anchors.
   - Continuous candidates need a candidate-based graph or a well-defined grouping by nearest reference/state.

4. **I/O and assignment-map compatibility**
   - Downstream code expects assignment records that can be converted into `ori` fields.
   - Mitigation: keep nearest `iproj` for compatibility but store and apply actual continuous Euler angles.

5. **Symmetry ambiguity**
   - Continuous sampling near symmetry-equivalent orientations can duplicate candidates.
   - Mitigation: fold or deduplicate samples under the point group where possible.

6. **On-the-fly extraction cost**
   - Continuous candidates cannot rely on precomputed full reference banks.
   - Mitigation: sample small candidate sets and compare against reduced discrete grids, not against zero-cost references.

## Recommended First Milestone

The first useful milestone is:

1. Add a continuous sparse candidate representation that can emulate current sparse grid candidates.
2. Run a grid-emulation assignment test.
3. Add tangent sampling around previous orientation for a synthetic known-orientation case.
4. Compare orientation error before and after continuous sampled candidates.

Only after this works should we add coarse-discrete mixture centers and performance comparisons against `eul_prob_tab_neigh`.

## Long-Term Direction

If validated, the `eul_prob` family could evolve from:

```text
probability table over fixed projection references
```

to:

```text
particle-local SO(3) proposal model -> sparse exact PFTC evaluations -> global assignment
```

This would keep SIMPLE's current strengths:

- exact PFTC evaluation once a polar section exists
- stochastic hard assignment
- compatibility with state and shift refinement workflows
- sparse global assignment for distributed execution

while introducing a true continuum over SO(3) where it matters: candidate generation.
