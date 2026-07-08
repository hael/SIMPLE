# Prior-Guided Probability Neighborhood Policy

This document records the policy for prior-guided probability neighborhoods in
late 2D and 3D refinement stages. The current 2D implementation is deliberately
small: it writes and consumes a per-particle ranked class list. That rank-list
mode is a testable stepping stone toward the richer per-particle top-K posterior
model that should be implemented for 3D.

The goal is to preserve useful probability-table information that is currently
collapsed into one sampled hard assignment. The preserved information should
act as a prior-guided neighborhood model for later search. It should accelerate
late-stage refinement when the posterior is concentrated and improve
robustness when the posterior identifies plausible alternative classes or
projection directions.

The primary near-term objective is to condense probabilistic registration
evidence from mature pre-neighborhood stages into a per-particle ranked
candidate model. With the prior-guided 2D schedule, stage 5 produces
`posterior_topk_stage05.dat` and stage 6 consumes it. In 2D this artifact is a
rank prior over classes, not yet a full posterior. In 3D the target artifact
should carry probabilities or scores, in-plane/shift support, and concentration
diagnostics so the neighborhood size can become adaptive per particle.

Related policies:

- [abinitio2D_policy.md](abinitio2D_policy.md)
- [abinitio3D_policy.md](abinitio3D_policy.md)
- [refine3D_policy.md](refine3D_policy.md)
- [importance_sampling_fractional_update_policy.md](importance_sampling_fractional_update_policy.md)
- [sigma_calculation_policy.md](sigma_calculation_policy.md)

## 1. Core Model

Probabilistic pre-alignment remains a particle-domain pre-step followed by a
hard-assignment matcher update. The new model does not turn SIMPLE into a
soft-assignment volume-integrated EM workflow.

The general target artifact is a per-particle top-K posterior list:

```text
particle i:
    candidate 1, score/probability, in-plane, shift, rank, diagnostics
    candidate 2, score/probability, in-plane, shift, rank, diagnostics
    ...
    candidate K, score/probability, in-plane, shift, rank, diagnostics
```

For staged workflows, this artifact is treated as a condensed evidence model,
not only a single-stage snapshot. The preferred long-term policy is to
aggregate ranking evidence from designated late probabilistic stages while
allowing early basin-finding stages to be excluded.

The existing hard assignment remains authoritative for updating project
orientation records:

```text
assignment.dat -> matcher -> ptcl2D/ptcl3D hard update
```

The prior artifact is an auxiliary search-prior handoff:

```text
prior artifact -> next probabilistic/neighborhood stage -> candidate
ordering, neighborhood selection, and in-plane/shift refinement seeding
```

Because SIMPLE's existing probability-table values are score-derived and
workflow-normalized, this document uses "posterior" as a practical model term.
The current 2D implementation is more modest: it stores only per-particle
ranked class ids. The full posterior artifact, especially for 3D, should record
enough scores, weights, and diagnostics to distinguish a concentrated posterior
from a broad or unreliable one.

Condensation model for selected source stages:

```text
choose source stages S (2D rank-prior implementation: S = {5})
aggregate_score_i(c) = reducer({score_i(c, s) for s in S when available})
rank_i = order candidates by aggregate_score_i (or aggregate_weight_i)
```

The 2D rank-prior implementation uses a single source stage, so no cross-stage
reducer is needed. The 3D posterior implementation can later use a weighted
mean, robust mean, or best-stage score and record the reducer in artifact
metadata.

## 2. Ownership

`simple_commanders_prob.f90` owns orchestration:

- choose and record the outer sampled particle subset
- dispatch probability-table workers
- aggregate table outputs
- write `assignment.dat`
- write the optional prior artifact

`simple_eul_prob_tab.f90`, `simple_eul_prob_tab_neigh.f90`, and
`simple_eul_prob_tab2D.f90` own candidate-level scoring and top-K extraction
from dense or sparse probability tables.

`simple_strategy2D_matcher.f90` and `simple_strategy3D_matcher.f90` own
consuming hard assignment and, in prior-guided modes, consuming the
prior artifact as a candidate ordering or restricted neighborhood. The matcher
still owns final class/projection/state/in-plane/shift updates and Euclidean
sigma updates after assignment.

Class-average restoration, explicit class-average assembly, partial
reconstruction writing, `volassemble`, FSC calculation, automasking, and
nonuniform filtering do not own the prior-guided posterior model. They consume the
same hard updated project state and sampled-update bookkeeping as before.

## 3. Artifact Contract

The prior artifact should be separate from `assignment.dat`.

Reasons:

- `assignment.dat` is the existing single-assignment handoff consumed by
  `strategy2D_prob` and `strategy3D_prob`.
- Existing workflows and tests can ignore the new artifact.
- Fallback behavior is simple when the top-K artifact is absent or invalid.
- The artifact can evolve from the current 2D rank-list format into a richer
  posterior format without changing the hard-assignment binary format.

Suggested file families:

```text
posterior_topk_stage05.dat
posterior_topk.dat
posterior_topk_part*.dat
posterior_topk_neigh.dat
posterior_topk_neigh_part*.dat
```

The exact names should follow the local constants and file-body conventions
used for `ASSIGNMENT_FBODY` and `DIST_FBODY`.

The current 2D rank-prior artifact is intentionally compact:

```text
header:
    magic
    version
    nparticles represented
    nclasses
    K
    source iteration
    start iteration
body:
    pinds(nparticles)
    counts(nparticles)
    class_ids(K, nparticles)
```

This format is enough to test whether a learned class-rank neighborhood can
accelerate or stabilize the final 2D probabilistic stage.

The full top-K posterior row model, targeted first for 3D, is:

```text
pind
rank
candidate id
raw distance or objective score
normalized weight
in-plane index
x shift
y shift
has_shift flag
search fraction or touched fraction
entropy or effective candidate count
source mode/version
source stage window
```

For 2D, the candidate id is `class`.

For 3D, the candidate id is either:

- full reference id, equivalent to `(state - 1) * nspace + proj`
- explicit `(state, proj)` pair
- optional coarse subspace id when the artifact is used only as a subspace
  neighborhood prior

The preferred 3D implementation is explicit `(state, proj)` rows. This keeps
the artifact readable and avoids encoding assumptions when `nspace` changes
between stages.

## 4. Top-K Extraction

Top-K extraction should happen after partition tables have been merged into the
global probability object and before the object is destroyed.

The existing flow is:

```text
prob_align*:
    sample particles
    run prob_tab* workers
    read/merge partition dist tables
    ref_assign/state_assign
    write assignment.dat
```

The prior-guided flow is:

```text
prob_align*:
    sample particles
    run prob_tab* workers
    read/merge partition dist tables
    extract prior rows
    ref_assign/state_assign
    write assignment.dat
    write posterior_topk*.dat
```

Prior extraction must use only the sampled particle subset already selected by
the probabilistic pre-step. It must not resample particles.

The current 2D implementation selects top ranked classes from the sparse
evaluated candidates in stage 5 and stores class ids only. It uses a fixed
fraction of classes as the retained rank list. If no sparse candidates are
available for a particle, it falls back to the available dense class table for
that particle.

The future full posterior implementation should select top-K from all valid
candidates for dense tables. Sparse tables should select top-K only from
evaluated candidates and should record the touched fraction. If sparse coverage
is too small or no valid candidate exists for a particle, the top-K artifact
should either omit that particle or include a clearly marked fallback row.

In the full posterior model, weight normalization should be local to each
particle:

```text
w_i(c) = normalize(score_i(c) over retained candidate set)
```

The retained set may be all candidates for dense tables or all evaluated
candidates for sparse tables. The full posterior artifact should record whether
the weights come from dense or sparse coverage.

When rows are condensed from multiple source stages, metadata should also
record stage support (for example contributing stage count) so consumers can
downweight weakly supported rankings.

## 5. 2D Implementation Policy

### Current 2D Flow

`abinitio2D` uses staged `cluster2D` invocations. The default late stages switch
from sampled SNHC to probabilistic assignment:

- stages 1-3: `refine=snhc_smpl`
- stages 4-5: `refine=prob_snhc`
- stage 6, only when prior-guided refinement is requested: `refine=prob_prior`
- terminal all-particle refresh after sampled staged updates: dense
  `refine=prob`

`prob_align2D` owns outer subset sampling and writes the existing hard
assignment artifact. `prob_tab2D` reproduces the sampled subset and fills
`eul_prob_tab2D%loc_tab(class, particle)`. `cluster2D_exec` later reproduces
the same subset and consumes `assignment.dat` through `strategy2D_prob`.

### Implemented 2D Rank-Prior Model

For every sampled particle represented in stage 5, store:

```text
top K class ids ranked by stage-5 sparse probabilistic score
```

This defines a prior-guided class neighborhood for the final 2D search pass:

```text
N_i = ranked class ids from posterior_topk_stage05.dat
```

The 2D implementation intentionally does not yet store normalized weights,
in-plane rows, shifts, entropy, effective K, margins, or touched-fraction
diagnostics. It is a simplified, testable `refine=prob_prior` mode for asking a
narrow question: does a late per-particle class-rank prior improve or accelerate
the final 2D probabilistic stage?

Implemented mode:

```text
refine=prob_prior
```

The implementation routes through the existing `prob_align2D` and
`cluster2D_exec` ownership split.

### 2D Consumption

The prior-guided 2D consumer should:

1. read `assignment.dat` exactly as today for the hard update
2. read `posterior_topk_stage05.dat`
3. for each particle with a valid prior row, move the ranked class list to the
   front of the sparse search order
4. evaluate exactly the stored prior classes for particles with a prior row
5. refine in-plane/shift using the existing sparse 2D helpers
6. fall back to existing `prob_snhc`-compatible stochastic sparse behavior for
   particles without a valid prior row or when the artifact is missing or
   dimension-incompatible

The current 2D neighborhood size is fixed by the stored rank-list length. A
later 2D extension may use concentration diagnostics from objective/probability
statistics to set `K_i` between configured bounds:

```text
concentrated posterior -> small K_i
broad posterior -> larger K_i
low sparse coverage -> broaden further or fallback
```

The top-K artifact should not directly restore class averages. Class-average
restoration continues to consume the hard assignment after the matcher update.

### 2D Stage Scheduling

For the prior-guided 2D schedule, the first practical implementation schedule is:

- stages 1-3: basin-establishment and stabilization stages; do not require
  prior construction
- stage 4: turn on `refine=prob_snhc`
- stage 5: write `posterior_topk_stage05.dat` from probabilistic registration
  when the user requested prior-guided refinement
- stage 6: consume that prior through new `refine=prob_prior`
  using the stored per-particle class rank lists
- terminal dense all-particle refresh: optionally produce diagnostics, but do
  not require prior-guided consumption

This schedule keeps prior construction late and targeted, while giving the
consumer one final refinement stage to exploit the learned per-particle class
neighborhoods. The controller should expose `PROB_PRIOR_STAGE = 6` so the
consumer stage is named explicitly instead of inferred from the probabilistic
start stage.

### 2D Refine Mode Contract: `refine=prob_prior`

`refine=prob_prior` is a new opt-in 2D refine mode intended to be orthogonal to
existing `prob_snhc` and `prob` behavior.

Mode intent:

- keep the probabilistic pre-alignment ownership split (`prob_align2D` pre-step,
  `cluster2D_exec` matcher/update)
- consume a prior artifact written by stage-5 probabilistic registration
- evaluate stored class rank lists rather than stochastic class proposals for
  particles with a valid prior row

Fallback in `refine=prob_prior`:

- if the prior artifact is missing or incompatible, continue with existing
  `prob_snhc`-compatible behavior
- if a particle is not represented by a valid prior row, use the usual sparse
  stochastic class order for that particle
- flatness, posterior concentration, and sparse-coverage diagnostics are not
  available in the current 2D rank-list artifact

Orthogonality requirements:

- no behavior change for existing refine values (`prob_snhc`, `prob`,
  `snhc_smpl`, ...)
- no format change to `assignment.dat`
- no mandatory new command-line flags for existing workflows
- activation only when `refine=prob_prior` is explicitly requested

### 2D Binary Artifact Policy (Fast I/O)

The 2D prior artifact should use binary handling patterns equivalent to
`assignment.dat` handling in the probabilistic pipeline.

Suggested files:

```text
posterior_topk_stage05.dat
posterior_topk_stage05_part*.dat
```

I/O behavior:

- the 2D rank-prior implementation writes `posterior_topk_stage05.dat` after partition
  rows have been merged into the global probability table
- a future distributed extension can add `posterior_topk_stage05_part*.dat`
  writers if global merge memory or wall-time becomes limiting
- writer/reader use the same compact binary stream I/O style as
  `assignment.dat`/`DIST_FBODY` flows
- include compact header with magic/version and basic run metadata; future
  richer posterior artifacts should extend this enough to reject stale or
  incompatible priors quickly

Performance notes:

- single-pass write, single-pass read, row-major particle grouping
- avoid text parsing and per-row dynamic allocation in hot paths
- preserve partition-local write path to keep distributed scaling behavior

## 6. 3D Implementation Policy

### Current 3D Flow

`abinitio3D` stages are configured by `simple_abinitio_controller.f90`.
Middle stages use `refine=prob`; late stages use `refine=prob_neigh`.

`prob_align` and `prob_align_neigh` own 3D probability pre-alignment. They
merge partition tables into either `eul_prob_tab` or `eul_prob_tab_neigh`,
perform hard assignment, and write `assignment.dat`. `refine3D_exec` then
consumes the hard assignment with `strategy3D_prob`.

`prob_neigh` currently chooses sparse neighborhoods using mode-specific local
rules:

- `state`: score coarse subspace representatives independently per state
- `sum`: score coarse subspace representatives by summed distance across states
- `geom`: use the coarse subspace containing the previous projection
- `shc` or `snhc`: direct stochastic sparse search

### Proposed 3D Top-K Posterior Model

The 3D implementation should not merely copy the simplified 2D rank-list
artifact. The 2D mode is useful because it is small and testable, but 3D needs
the richer posterior form: scores or weights, pose support, compatibility
metadata, and concentration diagnostics.

For every sampled particle, store:

```text
q_i(state, projection), top K state/projection candidates
best sampled/refined in-plane index for each candidate
shift for candidates whose shift was refined
posterior concentration diagnostics
optional coarse subspace id
```

This defines a prior-guided projection/state neighborhood:

```text
N_i = top-K (state, projection) candidates from q_i(state, projection)
```

The first consumer should be a new `prob_neigh_mode`:

```text
prob_neigh_mode=prior
```

`prob_neigh_mode=prior` should be implemented as a prior-guided variant of the
existing `prob_neigh_mode=snhc` sparse search: keep the same sparse-scoring
execution pattern, but replace random/stochastic reference selection with
selection from saved ranked candidates produced by earlier probabilistic
refinement. Unlike the current 2D rank-prior artifact, the 3D artifact should
retain enough posterior evidence to tune the neighborhood width per particle.

The prior rows consumed by `prob_neigh_mode=prior` should come from earlier
probabilistic refinement stages (default condensation window: stages 4-5 with
`refine=prob`), not from `prob_neigh` geometry-driven sparse modes.

In this mode, `simple_eul_prob_tab_neigh` should build the sparse evaluated
candidate set from saved ranked rows. If the saved rows contain fine
projection ids, the implementation can either:

- evaluate exactly those fine `(state, proj)` candidates and refine the best
  in-plane/shift candidates
- map the fine projections to coarse subspaces and evaluate all references in
  those subspaces
- use a hybrid policy: always include top-K fine candidates and also include
  their neighboring coarse subspaces

The hybrid policy is the safest first `prior` mode. It preserves local
robustness around each prior-ranked candidate while keeping the candidate set
smaller than a dense probability table.

### 3D Consumption

The prior-guided 3D consumer should:

1. read `assignment.dat` exactly as today for the hard update
2. read `posterior_topk*.dat` when `prob_neigh_mode=prior`
3. map retained rows to `(state, proj)` and optionally coarse subspaces
4. evaluate the prior-ranked candidate set through the existing PFTC sparse
  scoring path used by `prob_neigh`, replacing stochastic `snhc` candidate
  selection with prior-ranked selection
5. refine in-plane and shift for the best candidates according to existing
   `npeaks_inpl` and shift-search policy
6. fall back to `prob_neigh_mode=snhc` (or user-requested non-prior mode) when
  the artifact is absent, stale, incompatible with the current `nspace`, or
  too flat

Neighborhood size in `prior` mode should be adaptive per particle using
objective/probability diagnostics from condensed rows. Concentrated rankings
can use smaller candidate sets; broad or unreliable rankings should increase
candidate count or trigger fallback.

The fallback default should match the stage policy:

- single-state late `abinitio3D`: fallback to `snhc`
- docked or multi-state late `abinitio3D`: fallback to `sum` or `snhc`
- explicit user override: respect the requested non-prior mode

### 3D Stage Scheduling

For ordinary `abinitio3D`, the first practical schedule is:

- stages 1-3: warmup basin-establishment stages; optional to exclude from
  condensation
- stages 4-5: produce ranked state/projection evidence and write condensed rows
  (default condensation window: stages 4-5)
- stage 6 and later: consume condensed rows through `prob_neigh_mode=prior`
  with adaptive per-particle neighborhood size

For `multivol_mode=independent`, the default schedule stops at stage 5, before
late `prob_neigh`. That mode can still write top-K diagnostics from stage 5,
but prior-guided neighborhood consumption only applies when the user extends
the run past the default independent schedule.

For `multivol_mode=docked`, the split starts a new multi-state update epoch.
Do not carry pre-split top-K state/projection rows into the post-split epoch.
The first post-split stage should either produce a fresh top-K artifact or use
the full-update pooled-state `prob_neigh_mode=state` stabilization behavior.
Later post-split stages can consume prior rows produced within the post-split
epoch through `prob_neigh_mode=prior`.

## 7. Staleness and Compatibility

A prior artifact is valid only for a compatible search context. The current 2D
rank-prior file records a compact header with magic, version, particle count,
class count, retained K, source iteration, and start iteration. That is enough
for the first 2D test mode, but it is not the full compatibility contract.

The full 3D posterior file header should record:

```text
format version
workflow kind: 2D or 3D
source program and refine mode
source iteration/stage
nparticles represented
nclasses for 2D
nstates and nspace for 3D
box/smpd/crop metadata if needed
objective kind
prob_athres
K
whether rows came from dense or sparse tables
```

The full posterior consumer must reject or ignore an artifact when:

- particle ids do not match the current sampled subset
- class count differs in 2D
- state count differs in 3D
- 3D `nspace` differs and no reliable projection/subspace remapping exists
- the artifact was written before a docked split or other explicit update-epoch
  reset
- the artifact is too sparse for the requested mode

For 3D schedules that change `nspace` between stages, the first 3D implementation
should either store enough angular metadata to remap projections or treat the
artifact as invalid across `nspace` changes. The conservative path is to
consume prior rows only when `nspace` is unchanged, then later add a
projection-remapping helper through the Euler-space APIs.

## 8. Posterior Diagnostics

The current 2D rank-prior artifact does not carry posterior diagnostics; it only
stores ranked class ids. This keeps the first implementation easy to validate
but means 2D `prob_prior` cannot distinguish a sharp rank list from a flat or
under-sampled one except through missing/incompatible-row fallback.

Every particle-level top-K group in the full posterior artifact should carry
concentration diagnostics:

- entropy over retained candidates
- effective candidate count
- top-1 minus top-2 margin
- retained probability mass if known
- sparse touched fraction when applicable

Full posterior consumers should use these diagnostics to decide whether to
narrow or broaden search:

- concentrated posterior: evaluate a small prior-guided neighborhood
- broad posterior: increase K or fall back to existing broader search
- low touched fraction: treat weights as a sparse prior, not a complete
  posterior
- invalid or empty posterior: ignore the artifact for that particle

This prevents the prior-guided model from becoming an overconfident trap in
ambiguous or under-sampled regions.

Adaptive neighborhood policy can be implemented with threshold rules over
entropy, top-1/top-2 margin, effective candidate count, and touched fraction.
The first 3D implementation should use simple thresholds and log realized `K_i`
statistics to support later tuning.

## 9. Implementation Steps

Implementation status and next steps:

1. Implemented for 2D: write a compact class-rank prior from
   `eul_prob_tab2D` during the stage-5 `prob_snhc` pass.
2. Implemented for 2D: consume `posterior_topk_stage05.dat` in stage 6 when
   `refine=prob_prior` is requested.
3. Implemented for 2D: preserve `assignment.dat` as the hard-assignment
   handoff and use the prior artifact only as search-order/neighborhood input.
4. Next for 2D: add tests that confirm stage 5 writes the artifact, stage 6
   consumes it, missing/incompatible artifacts fall back, and existing refine
   modes remain unchanged.
5. Next for 3D: add a full top-K posterior data model and binary reader/writer
   near the existing probability-table code.
6. Next for 3D: add dense top-K extraction for `eul_prob_tab`, including scores
   or weights, in-plane/shift support, and concentration diagnostics.
7. Next for 3D: write 3D `posterior_topk.dat` from `prob_align`.
8. Next for 3D: add `prob_neigh_mode=prior` to consume 3D top-K rows when
   compatible, with adaptive per-particle neighborhood sizing.
9. Later: extend sparse top-K extraction for `eul_prob_tab_neigh`, recording
   touched fraction explicitly.
10. Later: add integration tests that compare shared-memory and distributed
   artifact equivalence and validate adaptive neighborhood-size decisions for
   concentrated versus broad posteriors.

2D first implementation detail:

- stage 4 turns on `refine=prob_snhc`
- stage 5 writes `posterior_topk_stage05.dat`
- stage 6 with `refine=prob_prior` reads and consumes that file
- if file validation fails, execution continues with `refine=prob_snhc`

The 2D path is the lower-risk first target because class ids are stable across
late stages and do not require projection remapping. It should be treated as a
rank-prior prototype, not as the final posterior artifact contract.

## 10. Invariants

- The outer sampled particle subset is selected once by `prob_align*`.
- Probability-table workers and matchers reproduce the recorded subset.
- Prior extraction does not resample particles.
- `assignment.dat` remains the hard-assignment handoff.
- The prior artifact is optional and must have a safe fallback path.
- Matchers remain responsible for final particle/class/projection/state,
  in-plane, shift, and sigma updates.
- Class-average restoration and volume assembly do not consume prior rows.
- Prior-guided neighborhoods must not cross an explicit update-epoch reset,
  such as the docked multi-state split.
- Sparse full-posterior rows must record sparse coverage diagnostics.
- 3D full-posterior consumption must validate `nstates` and `nspace`
  compatibility.
- Existing `prob`, `prob_snhc`, and `prob_neigh` behavior must remain
  available without `prob_neigh_mode=prior` artifacts.
- Existing 2D modes must remain unchanged unless `refine=prob_prior` is set.

## 11. Review Checklist

For any implementation of this policy, check:

- Does the new artifact remain separate from `assignment.dat`?
- Is the artifact written after probability-table aggregation and before the
  global probability object is destroyed?
- Does the matcher reproduce the same outer sampled subset before consuming
  prior rows?
- For the 2D rank-prior mode, are missing, incompatible, and unrepresented
  particle rows handled by fallback behavior?
- For the full posterior mode, are absent, stale, flat, sparse, and
  incompatible artifacts handled by fallback behavior?
- Does 2D keep class-average restoration tied to the hard assignment only?
- Does 3D keep volume assembly tied to hard updated project state only?
- Are shared-memory and distributed paths writing equivalent prior artifacts?
- For the full posterior mode, are entropy, margin, effective candidate count,
  and sparse touched fraction available for diagnostics?
- Does docked multi-state refinement avoid carrying pre-split prior
  neighborhoods into the post-split epoch?
- Does `refine=prob_prior` remain orthogonal (no behavior change) for existing
  2D refine modes when the new mode is not selected?
