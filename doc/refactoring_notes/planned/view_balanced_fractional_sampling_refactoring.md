# Canonical view-balanced fractional-update sampling refactoring

Date: 2026-10-04

Status: design direction accepted; implementation not started.

Validation level: static source inspection only. No Fortran source was changed,
compiled, or executed while preparing this plan.

This is the single living design record for this refactor. Update its contract,
review findings, implementation milestones, and validation evidence here rather
than creating parallel specification and plan documents.

## 1. Decision summary

Replace the 3-D fractional-update use of `balance=yes` and `partition=yes`
with one explicit selector:

```text
view_balance_mode=none|cavg|proj
```

The values mean:

| Mode | Sampling-group source | Scientific meaning |
|---|---|---|
| `none` | No groups | Select globally from the lowest `updatecnt` tiers. |
| `cavg` | Hierarchical clustering of aligned selected class averages | Equalize the sample over view groups inferred from class-average images, not over individual 2-D class labels. |
| `proj` | Symmetry-aware bins of current trusted `ptcl3D` projection directions | Equalize the sample over occupied projection-direction bins. |

There is deliberately no mode that gives every raw 2-D class an equal quota.
That policy is not view-balanced when a preferred view produces many class
representatives. If the value `class` is ever preferred over `cavg` for UI
wording, it must still mean class-average HAC view groups; it must never restore
one-group-per-class sampling.

The mode selects how groups are constructed. The existing outer
fractional-update bookkeeping remains authoritative:

- `update_frac` determines the target sample size;
- `sampled` identifies the current sampling round;
- `updatecnt` is persistent selection history;
- probabilistic pre-alignment samples once, and probability-table generation
  plus the matcher reproduce that exact subset;
- trailing reconstruction consumes realized update fractions and never chooses
  a particle subset.

## 2. Problem statement

The current `balance` flag does not define a scientific grouping. It tells
[`sample_ptcls4update3D`](../../../src/main/strategies/search/simple_matcher_smpl_and_lplims.f90)
to read `clssmp.bin` and apply equal quotas to whatever groups the caller wrote.
The same file and sampler therefore represent several incompatible meanings.

### 2.1 Current implementations

| Caller | Current setup | Actual groups |
|---|---|---|
| `solve3D balance=yes partition=no` | `get_class_sample_stats` over selected classes | One group per selected 2-D class. |
| `solve3D partition=yes` | `make_view_partition_class_samples` | Average-linkage HAC groups of aligned selected class averages. |
| `refine3D_states` | Local `make_projdir_class_samples` helper | Occupied bins of a 5,000-direction 3-D projection grid. |
| `classify3D_refs balance=yes partition=no` | `get_class_sample_stats` | One group per selected 2-D class. |
| `classify3D_refs partition=yes` | Read `cls2D.cluster` | Pre-existing cluster labels produced elsewhere; no HAC is performed here. |
| `sample_classes partition=yes` | Read `cls2D.cluster` | Pre-existing cluster labels used for project subdivision, not the `solve3D` HAC contract. |

The `partition` parameter is additionally confusable with scheduler partitions,
distributed worker partitions, even/odd partitioning, and project subdivision.
It is not an acceptable long-term name for view sampling.

### 2.2 Why balancing raw classes is wrong

Let one physical view produce ten well-populated 2-D classes while another view
produces one. Equal per-class quotas allocate roughly ten times as many sampled
particles to the first view until capacities saturate. The sampler balances
class representatives, not projection directions.

The `solve3D partition=yes` implementation addresses this by aligning selected
class averages, forming a view-similarity distance, and applying average-linkage
clustering before quotas are assigned. All similar class representatives then
share one view-group quota. This must become the canonical class-average-based
definition across workflows.

### 2.3 Additional correctness problems exposed by the refactor

The current grouped allocator increments every unsaturated group in one loop
body and stops only after the total reaches or exceeds the target. Except for a
special trim in the sticky-cohort path, a request for `K` particles can therefore
return up to `K + G - 1`, where `G` is the number of active groups.

The current `clssmp.bin` serialization records group indices, populations,
particle indices, and ranking scores, but not:

- which grouping definition produced it;
- the project or active-row layout it belongs to;
- symmetry or projection-grid settings;
- HAC settings;
- whether the artifact is stale after an add-on or continuation changes the
  particle layout.

Both issues are part of this refactor because a canonical mode is not reliable
without an exact quota contract and validated group provenance.

## 3. Public and typed parameter contract

Add the following typed parameter and UI metadata:

```text
view_balance_mode=none|cavg|proj
```

Use an enum/multi UI input, not a binary flag. The `parameters` field is a
normal character value registered in the parse pipeline and validated after
all sources have been resolved.

Mode-specific advanced controls are:

```text
nview_groups=20
view_group_crit=cc|sig|res|hybrid
view_nspace=5000
```

`nview_groups` and `view_group_crit` apply only to `cavg`.
`view_nspace` applies only to `proj` and is deliberately separate from the
stage's search `nspace`: the view-balance grid must not change accidentally
because a search stage changes its candidate density.

Validation rules:

1. `cavg` requires selected `cls2D` images, active `ptcl2D` class assignments,
   `nview_groups >= 1`, and a supported distance criterion.
2. `proj` requires trusted non-virgin `ptcl3D` orientations and
   `view_nspace >= 1`.
3. `none` rejects mode-specific controls rather than silently ignoring them.
4. Full-update mode does not build or read view groups. The selected mode may
   remain recorded for provenance, but it has no effect on the all-active
   subset.
5. Fill-in and final missing-update passes remain coverage policies and do not
   claim view balance.

Remove `balance` and `partition` from the 3-D sampling contract. They must not
be aliases for `view_balance_mode`. Release 4 does not need a backwards-
compatibility command alias.

`balance` may remain temporarily for unrelated project-selection operations,
but it must not be read by 3-D outer sampling. Dead or hidden `partition`
branches in project utilities should be removed or replaced with an
operation-specific selector after checking their callers; they must not be
described as the same view-balancing feature.

## 4. Sampling-group semantics

### 4.1 Common group contract

Both non-`none` modes produce one array of generic sampling groups. Each group
contains:

- a stable group identifier within the artifact;
- the eligible global particle indices;
- the population represented by those indices;
- an optional ranking key used only by `frac_best` or greedy selection;
- enough provenance to verify the artifact before use.

Groups cover every active eligible particle exactly once. Duplicate membership,
missing active rows, out-of-range particle indices, and an empty group set are
hard errors.

The group definition is pooled over conformational states, preserving current
behavior. Introducing state-by-view joint quotas is a separate scientific
change and is not part of this refactor.

### 4.2 `cavg`: class-average HAC groups

Preserve the scientific behavior of
[`simple_view_partition_sampling`](../../../src/main/strategies/search/simple_view_partition_sampling.f90):

1. select `cls2D` rows with `state > 0`;
2. align and compare the selected class averages under the class-average
   clustering parameters (`objfun=cc`, no CTF, `lp=6 A`, `trs=10`);
3. form the requested distance (`cc`, `sig`, `res`, or `hybrid`);
4. apply average-linkage HAC into `min(nview_groups, nselected_classes)` groups;
5. map every active particle into the HAC group of its 2-D class;
6. within a HAC group, preserve rank-percentile interleaving across source
   classes so that `frac_best` does not simply favor the class with the largest
   raw correlation scale;
7. do not write or reinterpret the project's `cls2D.cluster` field.

The class-average images are the evidence from which views are inferred. The
individual class labels are not themselves balancing groups.

### 4.3 `proj`: direct projection-direction groups

Extract the current local `make_projdir_class_samples` implementation from
`exec_refine3D_states` into the common view-group owner:

1. construct a symmetry-aware orientation grid with `view_nspace` directions;
2. assign every active particle's current `ptcl3D` orientation to a projection
   direction through the existing `set_projs` convention;
3. create one group per occupied projection bin;
4. retain particle correlation only as the optional ranking key;
5. treat the project's point group and the grid definition as artifact
   provenance.

This mode is valid only after an orientation-producing step has made the poses
scientifically meaningful. Random or untrusted placeholder orientations must
not be presented as projection-balanced evidence.

### 4.4 Exact quota allocation

For `N` active eligible particles and requested fraction `f`, retain the
existing whole-project target:

```text
K = min(N, max(1, nint(f * N)))
```

The grouped allocator must return exactly `K` distinct particles over the whole
project. Its saturation-aware water filling obeys:

1. no group receives more particles than its current eligible capacity;
2. among unsaturated groups, quotas differ by at most one;
3. the quota sum is exactly `K`;
4. a remainder is rotated or deterministically seeded by the sampling round so
   the same low-numbered groups are not favored every iteration;
5. within each group, lower `updatecnt` tiers are exhausted first;
6. the cutoff tier is sampled uniformly without replacement unless greedy
   selection is explicitly active;
7. `frac_best` restricts the eligible capacity before quotas are finalized;
8. a distributed worker may receive zero members of the already-selected
   global subset and must complete the existing empty transaction.

Rename the group-aware selection API so its name states the generic contract,
for example:

```text
sample4update_class  -> sample4update_groups
class_sample         -> sampling_group
clssmp.bin           -> view_sampling_groups.bin
```

The final names may follow local conventions, but no production API or artifact
should imply that projection bins or HAC groups are raw classes.

## 5. Workflow defaults and lifecycle

| Workflow | Resolved default | Build point and lifecycle |
|---|---|---|
| Normal particle `solve3D` | `cavg` | Build once before the first fractional stage. The selected class averages are fixed for the run. |
| External-volume `solve3D` | `proj` | Run the full CC pose-initialization pass first, then build projection groups before fractional refinement. |
| `solve3D_addon` | Inherit the base manifest's resolved mode | Rebuild groups for the augmented active layout using the inherited settings; never reuse base-run particle indices. |
| `refine3D_states` | `proj` | Build at the beginning of each parent frequency stage from the latest committed orientations; freeze through each child probabilistic transaction. |
| `classify3D_refs` | `proj` | Move group preparation after the full CC pose-initialization pass; then rebuild at parent frequency-stage boundaries. `cavg` remains an explicit alternative when class averages exist. |
| Direct base `refine3D` | `none` | The base command does not invent a workflow-specific scientific default. A non-`none` direct invocation prepares and validates its groups on the master before iteration dispatch. |
| `refine3D_auto` | `none` in this refactor | Changing automatic single-state refinement to `proj` is a separate behavior decision with its own comparison baseline. |

Within one probabilistic iteration the group artifact and chosen outer subset
are immutable:

```text
prepare/validate groups
    -> prob_align selects once and writes sampled/updatecnt
    -> prob_tab reproduces the subset
    -> refine3D matcher reproduces the subset
    -> partial reconstruction uses the matcher-resident images
    -> assembly consumes realized fractions
```

`sticky_class_sampling` is a class-specific name for a cohort eligibility
rule. Rename it to `sticky_view_sampling` or `sticky_sampling_cohort`. The
cohort remains an eligibility mask; it does not define the view groups. When a
projection group set is rebuilt, its capacities must be calculated from the
eligible cohort if sticky sampling is active.

## 6. Versioned group artifact

Replace the untyped matrix-only `clssmp.bin` handoff with a versioned view-group
artifact. It must contain or cover with an integrity digest:

- schema name and version;
- resolved `view_balance_mode`;
- project row count and particle-layout identity;
- active/eligible row identity;
- group count, group populations, particle indices, and ranking keys;
- `pgrp` plus `view_nspace` for `proj`;
- class selection identity, `nview_groups`, `view_group_crit`, and fixed
  alignment settings for `cavg`;
- whether sticky-cohort filtering is active;
- a content checksum or digest.

The artifact is transient workflow state, not a substitute for project
metadata. A missing artifact may be rebuilt by its owning master when the mode
and source evidence are available. A mismatched artifact is never consumed
silently.

Human-readable inspection output remains required:

- resolved mode and construction settings;
- source class/bin count and final group count;
- particles and assigned quota per group;
- maximum/minimum population share and sample share;
- saturation and exact target/realized totals;
- the class-to-HAC-group table and representative stacks for `cavg`.

Rename `VIEW PARTITION` log and file labels to `VIEW BALANCE` or `VIEW SAMPLING`
so they cannot be confused with execution partitions. Preserve the useful
class-average group stacks, but use mode-neutral group names.

## 7. Ownership

### 7.1 UI and parameters

The UI and typed-parameter lifecycle own:

- registration of `view_balance_mode` and mode-specific controls;
- per-program defaults and visibility;
- enum and consistency validation;
- rejection of unsupported evidence/mode combinations.

Touch the normal owners together:

- `src/main/params/simple_parameters.f90`;
- `src/main/params/simple_parameters_parse.f90`;
- `src/main/params/simple_parameters_phases.f90`;
- `src/main/ui/simple/simple_ui_solve3D.f90`;
- `src/main/ui/simple/simple_ui_refine3D.f90`;
- `src/main/ui/simple/simple_ui_heterogeneity.f90`.

### 7.2 Commanders and stage controllers

Workflow commanders choose the resolved mode, decide when the source evidence
is trusted, and request group construction. They do not implement HAC, pose
binning, or quota allocation inline.

The `solve3D` stage controller carries the resolved mode to child `refine3D`
commands. It strips expensive builder-only controls after the master has
created the artifact, but it must not strip the resolved mode: the matcher
dispatch needs to know whether it should consume groups or sample globally.

### 7.3 Search strategy/domain layer

Create or reshape a focused owner under `src/main/strategies/search` for view
group construction and artifact validation. If the implementation owns
persistent allocatable state, model it as an encapsulated type with private
components and `new`/`kill` symmetry. Keep HAC and projection builders behind
one public contract rather than leaving one nested in a commander.

### 7.4 Orientation collection

`simple_oris` owns exact group quota assignment and particle selection. It
continues to update `sampled` and `updatecnt`. It does not align class averages,
construct projection grids, choose workflow defaults, or own the persisted
artifact lifecycle.

### 7.5 Probabilistic and reconstruction paths

No sampling-policy change belongs in `simple_eul_prob_tab*`, the matcher
candidate search, `volassemble`, or trailing-reconstruction blending. Preserve:

- sample once and reproduce;
- outer subset selection versus inner candidate importance sampling;
- matcher single-read particle I/O;
- realized state-local update fractions for trailing reconstruction.

## 8. Manifest, continuation, split, and add-on contracts

The `solve3D` manifest currently records and replays `partition`, `nclust`, and
`clust_crit`. Replace them with the resolved mode and its canonical controls.
Because the scientific meaning and persisted input set change, bump the
manifest schema version. A version-1 manifest is rejected clearly rather than
silently translating ambiguous `partition` semantics.

`solve3D_addon` must:

1. inherit the base solution's resolved mode and builder settings;
2. validate that the new project contains the evidence required by that mode;
3. rebuild groups over the current augmented active layout;
4. record the new group artifact identity in the add-on run evidence;
5. preserve empty distributed-worker transactions.

The docked split checkpoint consumes the same validated group artifact used by
the surrounding run. It must not deserialize a mode-free legacy class-sample
file. Sticky cohort selection remains exact-K and is applied within the
resolved view groups.

Continuation rebuilds a missing or stale transient group artifact from the
committed project at the defined workflow boundary. It does not discover a
random file solely by its conventional filename.

## 9. Staged implementation

### Phase 0: baselines and contract tests

1. Capture current grouped-sampler behavior, including the target overshoot.
2. Construct a deterministic synthetic preferred-view case in which many
   classes represent one view and few classes represent another.
3. Record current `solve3D partition=yes` group tables and sampled counts on a
   small representative project.
4. Record the current `refine3D_states` projection-bin table and sampled counts.

Exit gate: expected old and new behavior is explicit before production logic
changes.

### Phase 1: generic group model and exact allocator

1. Introduce the generic group type and versioned file I/O.
2. Move the group sampler to an exact-K, saturation-aware allocator.
3. Preserve global lowest-`updatecnt`, greedy, `frac_best`, sticky-cohort, and
   distributed-empty behavior.
4. Add artifact validation and focused unit tests.

Exit gate: the generic sampler is independently tested and no workflow yet
depends on ambiguous class-specific names.

### Phase 2: canonical group builders

1. Refactor the existing `solve3D partition=yes` HAC implementation into the
   common `cavg` builder.
2. Extract `make_projdir_class_samples` from `exec_refine3D_states` into the
   common `proj` builder.
3. Give both builders common diagnostics and artifact publication.

Exit gate: both modes emit validated artifacts consumed by the same exact
group sampler.

### Phase 3: parameter and workflow migration

1. Add `view_balance_mode` and its mode-specific controls to parameters and UI.
2. Make normal `solve3D` default to `cavg` and remove raw-class balancing.
3. Make `refine3D_states` explicitly use `proj`.
4. Move `classify3D_refs` group construction after CC pose initialization and
   make it explicitly use `proj` by default.
5. Keep direct `refine3D` and `refine3D_auto` at `none` unless explicitly
   configured otherwise.
6. Rename sticky-cohort controls and update child command construction.

Exit gate: no 3-D workflow dispatches on `params%balance` or interprets
`params%partition` as view sampling.

### Phase 4: split, add-on, and manifest migration

1. Update docked split selection to the generic artifact and exact quota API.
2. Bump the `solve3D` manifest schema and record/replay the canonical mode.
3. Rebuild add-on groups over the augmented particle layout.
4. Update manifest, split-checkpoint, and add-on tests.

Exit gate: all persisted and replayed sampling policy has one unambiguous
meaning.

### Phase 5: cleanup and documentation

1. Remove 3-D uses of `balance`, `partition`, `nclust`, `clust_crit`,
   `class_sample`, `sample4update_class`, and `CLASS_SAMPLING_FILE` where their
   old semantics no longer apply.
2. Remove or rename unrelated hidden `partition` branches in project utilities.
3. Update policies, algorithms, help text, generated UI metadata, code maps,
   and relevant skills.
4. Record runtime validation and remaining follow-up decisions in this note.

Exit gate: repository searches find no stale claim that raw class balancing is
view-balanced or that `partition=yes` has a shared scientific meaning.

## 10. Test plan

Tests follow `doc/policies/test_environment_policy.md`: they are assertion-
based procedures beside the code they test, registered in an existing area
suite. This refactor does not add a CTest entry.

### 10.1 Unit tests

Extend `simple_oris_tester.f90` or add the nearest owning tester to assert:

- exact `K` for group counts that do not divide the target;
- equal water-filled quotas before saturation;
- correct redistribution after small groups saturate;
- deterministic or seeded rotation of remainder quotas;
- lowest-`updatecnt` selection within each group;
- `frac_best` capacity applied before quota allocation;
- sticky-cohort eligibility and exact-K behavior;
- no duplicate particle indices;
- a distributed range with no selected particles returns an empty sample when
  allowed;
- full update and fill-in do not consume view groups.

Add a focused view-group tester for:

- two synthetic class-average view families with unequal numbers of duplicate
  class representatives: `cavg` must create the requested HAC view groups and
  must not assign a quota per raw class;
- known Euler directions mapped to the expected symmetry-aware `proj` bins;
- every active row appearing exactly once;
- rejected stale layout, wrong mode, wrong symmetry/grid, corrupt checksum,
  duplicate membership, and missing membership;
- artifact write/read round trips for both modes.

Expected values must come from explicit synthetic geometry, a small independent
water-filling emulation, or closed-form counts, not from output captured from
the implementation under test.

### 10.2 Workflow and command-shape tests

Use existing workflow/high-level suites to assert:

- normal `solve3D` resolves to `cavg` and writes one validated group artifact;
- a full-sampling `solve3D` bypasses group construction;
- external-volume `solve3D` does not build `proj` groups before CC pose
  initialization;
- `refine3D_states` resolves to `proj` and rebuilds only at parent-stage
  boundaries;
- `classify3D_refs` prepares `proj` groups after CC pose initialization;
- probabilistic pre-alignment, tables, and matcher use the same sampled round;
- shared-memory and distributed runs realize the same global group counts for a
  fixed seed, including empty worker partitions;
- docked split selection is exact-K;
- `solve3D_addon` inherits the base mode but rebuilds membership for the new
  layout;
- a manifest with the old schema or ambiguous legacy keys is rejected clearly.

### 10.3 Scientific regression

The decisive synthetic case contains two true view families:

- view A is represented by many selected 2-D classes;
- view B is represented by one or a few selected classes;
- class populations and the total target are chosen so neither view saturates.

The expected result is:

- the recorded legacy one-class-one-group baseline is biased toward view A;
- `cavg` gives the two HAC view groups quotas differing by at most one;
- `proj` gives the corresponding occupied projection regions water-filled
  quotas, subject only to the declared binning and capacity rules;
- the exact total is `K` in both canonical modes;
- repeated rounds equalize `updatecnt` within each view group.

Maintainer-run representative workflows should compare orientation coverage,
group populations, sampled group shares, convergence, FSC, and final maps
between the baseline and the new policy. A view-biased real dataset is the
important acceptance case; a nearly isotropic dataset is the negative control
and should not regress materially.

## 11. Documentation and generated surfaces

Update at least:

- `doc/policies/importance_sampling_fractional_update_policy.md`;
- `doc/policies/3D/solve3D_policy.md`;
- `doc/policies/3D/solve3D_addon_policy.md`;
- `doc/policies/heterogeneity/refine3D_states_policy.md`;
- `doc/policies/3D/classify3D_refs_policy.md`;
- `doc/algorithms/sampling_and_fractional_updates.md`;
- `doc/algorithms/solve3d.md` and `doc/algorithms/refine3d.md` where they
  describe outer sampling;
- `.github/skills/simple-solve3d-importance-sampling/` and the refine3D
  sampling map;
- UI/NICE metadata generated from the program definitions;
- code-overview indexes after source stabilization.

Use `view group` for the scientific strata and reserve `partition` for actual
execution/data partitioning. Use `outer fractional-update sampling` for the
particle subset and `inner probabilistic candidate sampling` for pose
candidates.

## 12. Risks and mitigations

| Risk | Level | Mitigation |
|---|---:|---|
| HAC cost becomes a startup bottleneck. | Medium | Build `cavg` groups once per workflow, retain the current master-phase thread budget, and report distance-matrix time. |
| Projection groups are built from untrusted or random poses. | High | Validate mode preconditions and order external-reference CC initialization before `proj` construction. |
| Rebuilding `proj` groups changes sticky-cohort quotas unexpectedly. | Medium | Treat the cohort as the eligibility set when capacities are calculated and pin the behavior with exact-K tests. |
| Exact-K allocation changes realized update fractions relative to the overshooting legacy sampler. | Intended but scientifically relevant | Record the baseline, expose target and realized counts in logs, and compare convergence/trailing behavior. |
| Add-on or continuation consumes stale particle indices. | High | Validate layout identity and rebuild transient membership for changed layouts. |
| A child command loses the resolved mode after builder-only keys are stripped. | Medium | Distinguish the persistent resolved mode from construction controls in one helper and test command shape. |
| `cavg` and `proj` are compared as if they were numerically identical. | Medium | Document them as different evidence models sharing one quota contract; validate each against its own geometric expectation. |
| Broad renaming disturbs unrelated `balance` or scheduler partition behavior. | Medium | Restrict the first removal to 3-D sampling paths; audit unrelated callers before operation-specific cleanup. |

## 13. Non-goals

- Changing the nominal `nsample` or `update_frac` schedules.
- Changing probability-table candidate importance sampling or hard assignment.
- Introducing state-by-view joint quotas.
- Replacing fill-in or final missing-update coverage passes.
- Moving subset selection into volume assembly or trailing reconstruction.
- Adding a second particle-stack read during matcher reconstruction.
- Writing HAC labels into `cls2D.cluster` as durable project state.
- Reusing a base-run view-group membership file after an add-on changes the
  particle layout.
- Treating scheduler or distributed partitions as view groups.
- Enabling projection-balanced sampling in `refine3D_auto` without a separate
  scientific-policy decision.

## 14. Acceptance criteria

The refactor is complete when:

1. Every 3-D fractional-update workflow has one explicit resolved
   `view_balance_mode`.
2. No 3-D view-balanced path gives one quota to every raw 2-D class.
3. `cavg` means the canonical aligned-class-average HAC implementation in every
   caller.
4. `proj` means the canonical symmetry-aware current-orientation binning
   implementation in every caller.
5. Grouped sampling returns exactly the requested whole-project target and
   obeys the saturation and fairness invariants.
6. The group artifact is versioned, mode-aware, layout-validated, and rejected
   when stale or corrupt.
7. Probabilistic paths still sample once and reproduce the same subset.
8. `sampled`, `updatecnt`, realized update fractions, fill-in, trailing
   reconstruction, and matcher single-read I/O retain their existing ownership.
9. `solve3D` split and add-on workflows preserve their cohort and manifest
   contracts with the new artifact.
10. Repository documentation and UI no longer describe `balance=yes` or
    `partition=yes` as a shared view-balancing policy.
11. Focused unit and workflow tests pass when run by the maintainer, and the
    view-biased scientific comparison shows the intended correction without a
    material isotropic-data regression.
12. Static validation and all outstanding compilation/runtime checks are
    recorded in this living note without claiming unobserved results.

## 15. Implementation record

No implementation has started.

Preparation checks:

- current branch and dirty worktree were inspected before creating this note;
- current behavior was traced through solve3D, refine3D_states,
  classify3D_refs, the shared matcher sampler, `simple_oris`, split checkpoint,
  and solve3D manifest paths;
- no compilation or runtime tests were run, in accordance with repository
  policy.

Update this section after each implementation phase with the verified source
milestone, static checks, maintainer-run tests, and outstanding work.
