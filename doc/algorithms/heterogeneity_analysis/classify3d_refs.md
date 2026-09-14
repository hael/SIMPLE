# Reference-Guided 3D Classification (`classify3D_refs`)

## Problem

Given particles and a complete external set of state volumes, classify the
particles without importing incompatible amplitude scaling, noise statistics,
or pose history from those references. The common multi-state estimator,
bandwidth schedule, and final reconstruction are described in the
[heterogeneity-analysis overview](README.md).

This workflow is for references that may come from another dataset or docking.
When the states and particles already share a refinement lineage, use
[`refine3D_states`](refine3d_states.md).

## Algorithm

Direct Euclidean comparison is unsafe before the external maps and particle
data share a noise model and pose scaffold. The workflow therefore creates
data-derived checkpoint maps first:

```text
external maps, low-pass 15 A, nspace = 2500, at most 100,000 particles
  -> one greedy correlation-objective pass assigning (state, pose)
  -> noise power estimated from that pass's residuals
  -> checkpoint maps reconstructed from the assigned particles
  -> Euclidean probabilistic multi-state refinement from the checkpoints.
```

The correlation objective is scale-free, so the initialization tolerates
references whose amplitude spectrum differs from the data. The external maps
remain unchanged; after the checkpoint, only data-derived maps act as
references. Coverage is accepted only when at least 99 percent of sampled
particles were updated and every state is non-empty.

The subsequent state decisions use the shared whitened Euclidean loss at a
common bandwidth. Balanced fractional updates and the final all-particle
[reconstruction](../reconstruction.md) produce the delivered state and half
maps.

## Rationale

The correlation pre-pass changes the initialization objective, not the final
estimator. It supplies the data-derived noise model and pose scaffold that the
Euclidean likelihood requires while retaining the external maps only as an
initial source of structural hypotheses.

## Implementation

- Workflow: `src/main/commanders/simple/simple_commanders_refine3D.f90`.
- Correlation initialization:
  `src/main/simple_external_reference_pose_initialization.f90`.
- Frequency schedule: `src/main/simple_refine3D_stage_plan.f90`.
- State assignment: `src/main/simple_eul_prob_tab.f90` and
  `src/main/strategies/search/simple_strategy3D_prob.f90`.
- Policy: `doc/policies/3D/classify3D_refs_policy.md`.
