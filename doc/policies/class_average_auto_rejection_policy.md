# Class-Average Auto-Rejection Policy

## Scope

This document defines the policy for automatic rejection of 2D class averages in the streaming microchunked 2D workflow.

It covers:

- the rejection owner and execution point
- the ordered rejection criteria
- per-tier threshold policy
- state propagation and output artifacts
- invariants, compatibility rules, and review guidance

Primary implementation owners:

- rejection engine: `cluster2D_rejector` in `src/main/stream/simple_cluster2D_rejector.f90`
- workflow integration: `reject_cavgs` in `src/main/stream/simple_microchunked2D.f90`

## Public policy

Class-average auto-rejection is workflow-managed and runs after class averages are produced.

Policy statements:

- rejection is mandatory when chunk rejection is executed
- criteria are cumulative: once rejected by any criterion, a class remains rejected
- selection state uses binary semantics: `state=1` kept, `state=0` rejected
- class rejection is propagated to particle-level states and class assignments

## Ownership

### `reject_cavgs` (workflow owner)

Responsibilities:

- load class averages and class metadata
- instantiate `cluster2D_rejector` with class images and mask diameter
- select tier-specific threshold overrides where required
- execute criteria in policy order
- write selected and rejected class-average stacks
- propagate final class states to `os_cls2D`, `os_cls3D`, and mapped particle states
- write completion sentinel and bookkeeping fields

### `cluster2D_rejector` (criteria engine)

Responsibilities:

- hold per-class rejection mask
- apply criterion-specific rejection rules
- expose final `l_rejected` and integer `states`

Non-responsibilities:

- chunk orchestration
- project metadata I/O
- output stack naming policy

## Rejection criteria and order

Criteria are executed in this fixed order:

1. population (`reject_pop`)
2. FSC resolution (`reject_res`)
3. mask leakage/centroid tests (`reject_mask`)
4. local-variance z-score tests (`reject_local_variance`)

Order is currently policy-significant because criteria are cumulative and diagnostics are logged per pass.

## Criterion rules

### 1. Population rejection

Rule:

- reject class `i` if `pop(i) < ceiling(sum(pop) * threshold_fraction)`

Default engine fraction (`cluster2D_rejector`):

- `POP_PERCENT_THRESHOLD = 0.005`

Tier-specific workflow overrides are applied by `reject_cavgs` (see threshold matrix).

### 2. Resolution rejection

Rule:

- reject class `i` if `res(i) > RES_THRESHOLD`

Current threshold:

- `RES_THRESHOLD = 40.0` Angstrom

Boundary behavior:

- `res == 40.0` is kept

### 3. Mask rejection

Rule summary:

- build Otsu-binarized foreground after edge normalization and band-pass
- connected components (CCs) spanning full image diameter are pruned
- reject when any remaining CC centroid lies outside mask radius
- reject when largest CC has more than `MASK_THRESHOLD` pixels outside disc mask
- reject if no valid CC remains after pruning

Current threshold:

- `MASK_THRESHOLD = 10.0` pixels (outside-mask allowance)

Mask geometry:

- radius in pixels is derived from workflow mask diameter and class-average sampling

### 4. Local-variance rejection

Rule summary:

- compute inside/outside local-variance scores per class
- classes with both scores exactly zero are unconditionally rejected
- non-zero classes are robust-z-scored per region (excluding unconditional-zero rejects)
- reject class if one region is below strong threshold and the other below weak threshold (in either ordering)

Engine defaults:

- `LOCVAR_STRONG_THRESH = -0.5`
- `LOCVAR_WEAK_THRESH = -0.1`

Tier-specific workflow overrides are applied by `reject_cavgs` (see threshold matrix).

## Threshold matrix by chunk tier

The workflow applies these effective thresholds:

- REFCHUNK and MATCH:
- `reject_pop` threshold fraction: `DEFAULT_REF_POP_THRESH = 0.0025`
- `reject_local_variance`: `DEFAULT_REF_LOCVAR_STRONG_THRESH = -2.0`, `DEFAULT_REF_LOCVAR_WEAK_THRESH = -2.0`
- PASS 2:
- `reject_pop` threshold fraction: `DEFAULT_MICRO_P2_POP_THRESH = 0.0035`
- `reject_local_variance`: `DEFAULT_MICRO_P2_LOCVAR_STRONG_THRESH = -1.0`, `DEFAULT_MICRO_P2_LOCVAR_WEAK_THRESH = -1.0`
- fallback/default labels:
- `reject_pop` uses engine default `0.005`
- `reject_local_variance` uses engine defaults `(-0.5, -0.1)`

Global (non-tiered) thresholds currently remain:

- resolution: `40.0` Angstrom
- mask leakage allowance: `10.0` pixels

## State propagation and artifacts

After rejection is finalized for a chunk:

- class states are updated (`os_cls2D`, `os_cls3D`)
- particle states are remapped via `map2ptcls_state`
- particle `class` and `class_match` are zeroed for deselected particles

Artifacts:

- rejected class stack: `*_rejected.mrc`
- selected class stack: `*_selected.mrc`
- completion sentinel: `REJECTION_FINISHED`

Debug-mode artifact:

- optional `*_deselected` project metadata snapshot for inspection

## Preconditions and fail-fast behavior

Rejection pass is skipped or aborted in these cases:

- chunk not complete
- chunk already marked rejection-complete
- no class averages available
- class-average count mismatch with expected class count

Engine-level hard precondition:

- number of class metadata rows must match rejector class count

## Invariants

- `l_rejected(i)=.true.` implies final class state `0`
- `l_rejected(i)=.false.` implies final class state `1`
- rejection is monotonic within a pass: no criterion can re-enable a class
- state propagation must preserve consistency between class and particle views

## Verification policy

Unit coverage for engine behavior is maintained in:

- `src/main/stream/simple_cluster2D_rejector_tester.f90`

Minimum expected test categories:

- lifecycle safety (`new`/`kill`)
- population threshold and boundary behavior
- resolution threshold and boundary behavior
- cumulative rejection semantics
- `get_states`/`get_rejected` consistency
- deterministic mask and local-variance edge cases

Integration verification should additionally confirm:

- tier-specific threshold overrides are actually used by label
- state propagation modifies both class and particle spaces as intended
- expected artifacts and sentinel files are emitted once per completed chunk

## Governance and change control

Any change to rejection thresholds, criterion order, or criterion logic should include:

- update to this policy document
- update to in-code threshold comments where present
- corresponding unit/integration test updates
- explicit release-note entry for behavioral impact

Source-of-truth rule:

- executable constants in module declarations are normative
- prose comments describing values are informative and may drift if not updated
