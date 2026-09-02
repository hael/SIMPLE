# State Refinement and Reference-Guided Classification

SIMPLE exposes two multi-state workflows over the same base projection-
matching and reconstruction mechanisms. Their boundary is scientific intent
and reference provenance.

## Conformational State Refinement

`refine3D_states` starts from a same-lineage particle orientation scaffold and
state references. Its `pose_policy` determines the permitted geometry:

- `fixed`: preserve projection direction while optimizing in-plane angle and
  x/y translations and choosing the state stochastically;
- `local`: search the current geometric neighborhood across all states;
- `global`: perform full state-pooled probabilistic pose matching.

`global` is the default. Local angular, in-plane, and shift limits are derived
automatically and may be overridden individually.

Existing state labels/maps, stochastic initialization, FLEX initialization,
and an `abinitio3D` split checkpoint are supported. The default sample target
is 10,000 particles per state, capped at 100,000.

## Reference-Guided 3D Classification

`classify3D_refs` requires a complete external `vol1..volN` set. Because those
references may not share the particle pose history, the workflow first runs
one fixed-reference CC pose-initialization pass:

```text
fixed references at the common lpstart bandwidth + at most 100,000 particles
    -> one broad CC pose/state assignment pass
    -> residual sigma consolidation
    -> data-derived checkpoint maps
    -> Euclidean global probabilistic classification
```

For `classify3D_refs`, `lpstart` therefore controls both the fixed-reference
CC initialization pass and the start of the subsequent common Euclidean
frequency schedule. Its default is 10 Å, and every state uses the same value.

No reconstruction modifies the supplied references during the CC pass. The
checkpoint maps, not the external inputs, become the authoritative working
references for subsequent Euclidean matching.

## Common Frequency Schedule

Both workflows consume `simple_refine3D_stage_plan`. It divides the requested
iteration budget into short blocks and returns a common per-block low-pass,
crop, translation limit, and global iteration range. Every state uses the same
bandwidth at a given block so competitive evidence is comparable.

Each block is a normal base-`refine3D` call:

```text
candidate preparation -> hard pose/state update -> partial half maps
-> volume assembly -> filtered references for the next block
```

## Coverage and Final Maps

Stochastic fractional updates are balanced when suitable metadata exists. A
final missing-update pass assigns active particles with `updatecnt=0` without
replacing the last staged maps. Only then does a fresh all-particle
`reconstruct3D` pass produce authoritative native-sampling state maps, half
maps, FSCs, resolution products, and orthogonal reprojections.

## Ownership

- Workflow orchestration:
  `src/main/commanders/simple/simple_commanders_refine3D.f90`
- CC pose initialization:
  `src/main/simple_external_reference_pose_initialization.f90`
- Frequency planning: `src/main/simple_refine3D_stage_plan.f90`
- Base refinement:
  `src/main/strategies/parallelization/simple_refine3D_strategy.f90`
- Search and assignment:
  `src/main/strategies/search/simple_strategy3D_prob.f90` and
  `src/main/strategies/search/simple_strategy3D_matcher.f90`
- Policies: [refine3D_states_policy.md](../policies/refine3D_states_policy.md)
  and [classify3D_refs_policy.md](../policies/classify3D_refs_policy.md)
