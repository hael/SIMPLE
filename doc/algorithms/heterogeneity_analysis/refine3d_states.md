# Conformational State Refinement (`refine3D_states`)

## Problem

Given particles with a common pose scaffold from consensus refinement, refine
several same-lineage conformational maps and their particle assignments. The
shared estimator, bandwidth schedule, and final reconstruction are described
in the [heterogeneity-analysis overview](README.md).

The input may contain an existing label set with project state maps, the
[`flex_pca`](flex_pca.md) initializer (the default for state-0/1 input), the
stochastic state initializer (`flex=no`), or an `abinitio3D` docked checkpoint.
Supplied reference volumes are not accepted; that input belongs to
[`classify3D_refs`](classify3d_refs.md).

## Algorithm

`pose_policy` controls how much of the pose may move while states are decided:

- `fixed` freezes the projection direction and optimizes only the in-plane
  angle, shift, and state. This is classification given the consensus geometry.
- `local` permits the direction to move within the coarse Voronoi cell of its
  previous value. A bound `alpha` in degrees becomes
  `nspace_sub = min(5000, max(2, 2/(1 - cos alpha)))`, matching the number of
  cells to the solid angle of that cone.
- `global` (the default) performs a full coarse peak search across all
  directions and states in every iteration.

The workflow samples 10,000 particles per state, capped at 100,000 overall,
and balances the sample over projection-direction bins so every view
contributes to the state decision. State and pose assignments then follow the
common likelihood-ratio and probabilistic-assignment rules before each state is
reconstructed from its members.

## Rationale

A shared pose and map lineage makes direct Euclidean state competition valid.
The selectable pose policy separates nearly fixed-pose classification from
local or global conformational refinement, while balanced view sampling keeps
an initially larger state from winning solely through better angular coverage.

## Implementation

- Workflow: `src/main/commanders/simple/simple_commanders_refine3D.f90`.
- Frequency schedule: `src/main/simple_refine3D_stage_plan.f90`.
- State assignment: `src/main/simple_eul_prob_tab.f90` and
  `src/main/strategies/search/simple_strategy3D_prob.f90`.
- Policy: `doc/policies/heterogeneity/refine3D_states_policy.md`.
