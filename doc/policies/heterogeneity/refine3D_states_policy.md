# Conformational State Refinement Policy

This document defines the public and scientific contract of
`refine3D_states`. It is the canonical same-lineage multi-state workflow; no
historical command alias is registered or routed.

Related policies:

- [refine3D_policy.md](refine3D_policy.md)
- [classify3D_refs_policy.md](classify3D_refs_policy.md)
- [importance_sampling_fractional_update_policy.md](importance_sampling_fractional_update_policy.md)
- [nonuniform_filtering_policy.md](nonuniform_filtering_policy.md)

Primary implementation:

- `src/main/ui/simple/simple_ui_refine3D.f90`
- `src/main/exec/simple_exec_refine3D.f90`
- `src/main/commanders/simple/simple_commanders_refine3D.f90`
- `src/main/simple_refine3D_stage_plan.f90`
- `src/main/strategies/search/simple_strategy3D_prob.f90`

## 1. Scientific Scope

`refine3D_states` refines conformational states from a particle project with a
meaningful 3D orientation scaffold. State maps and poses must have the same
particle/project lineage. Independently derived references belong to
`classify3D_refs` because they require reference-conditioned CC pose
initialization before Euclidean refinement.

The wrapper owns state initialization, pose-policy selection, sampling,
frequency planning, coverage enforcement, and final reconstruction. Base
`refine3D` and its search strategies own candidate scoring and committed
particle updates. Reconstruction and volume modules retain numerical ownership
of state maps, half maps, FSCs, masks, and filtering.

## 2. Input and State Initialization

The project must contain active particles and meaningful 3D orientations. The
workflow may start from:

1. populated multi-state labels plus compatible project state maps;
2. state-0/1 input plus `nstates` and a complete `vol1..volN` set;
3. state-0/1 input plus distributed startup reconstruction;
4. `flex=yes`, which obtains labels and maps from `flex_pca`;
5. an `abinitio3D` split checkpoint.

Partial `vol1..volN` input is rejected. Existing multi-state labels determine
the effective state count; an explicit `nstates` must agree. Every accepted
state must be populated.

## 3. Pose Policy

`pose_policy` is the only public pose-search selector and defaults to
`global`.

| Policy | Scientific meaning | Internal search mapping |
| --- | --- | --- |
| `fixed` | Keep each particle's projection direction; stochastically choose its state while optimizing the in-plane angle and in-plane translations | `refine=prob_state` |
| `local` | Search state, projection direction, in-plane angle, and translations inside the current geometric neighborhood | `refine=prob_neigh`, `prob_neigh_mode=geom` |
| `global` | Search every pose degree of freedom through full state-pooled probabilistic matching | `refine=prob_neigh`, `prob_neigh_mode=state` |

`fixed` does not freeze the complete stored pose record. Its invariant is:

```text
projection direction after update = projection direction before update
```

The committed in-plane angle, x/y translations, state, correlation, and update
accounting come from the selected optimized state candidate. Discarding those
optimized in-plane values would violate the policy.

For `local`, angular, in-plane, and shift bounds are automatic. Advanced
`local_ang_bound`, `local_inpl_bound`, and `local_shift_bound` values override
only the corresponding automatic bound and are rejected for other policies.

The implementation-shaped `multivol_mode` and `prob_neigh_mode` controls are
not public inputs to this workflow; the commander derives them from
`pose_policy`.

## 4. Sampling and Frequency Planning

The automatic per-iteration target is 10,000 particles per state, capped at
100,000. If the active count exceeds the target, the wrapper uses stochastic
fractional updates and projection-balanced class sampling. Otherwise it uses a
full update.

`lpstart` and `lpstop` define one common frequency schedule for all states.
`simple_refine3D_stage_plan` returns short blocks containing the low-pass,
crop, translation limit, and global iteration range. Both
`refine3D_states` and `classify3D_refs` consume this planner. A state-specific
frequency schedule is outside the current contract because it would make
competitive evidence state-dependent.

Every planned frequency block is executed so the workflow reaches `lpstop`;
state-overlap diagnostics do not terminate the march at an earlier bandwidth.

## 5. `abinitio3D` Handoff

For docked multi-state ab initio work, `abinitio3D` owns the single-state
scaffold and split-checkpoint construction. The checkpoint preserves state
labels, maps, sampled/update metadata, the capped cohort, realized update
fraction, and iteration position. Post-split refinement is dispatched once to
`refine3D_states` with `pose_policy=local`.

The split checkpoint is constructed by
`simple_abinitio3D_split_checkpoint`; the old post-split state-refinement loop
is not a second production path.

## 6. Focus Evidence Boundary

A future focus mask may restrict evidence used to discriminate states, but it
must remain orthogonal to `pose_policy`. It must not become the authoritative
reconstruction mask, redefine the stored pose beyond the selected policy, or
replace state automasking and nonuniform filtering. No public focus-mask input
is enabled until the matcher can enforce this separation.

## 7. Completion and Outputs

Before final reconstruction, every active particle must have `updatecnt > 0`.
A missing-update pass fills remaining assignments without intermediate volume
reconstruction. Final state maps are then reconstructed and postprocessed from
all active particles at native project sampling. The workflow writes normal
state volumes, half maps, FSC/resolution records, diagnostic low-pass outputs,
and orthogonal reprojections.

## 8. Validation

User-side validation must cover:

- all three pose policies and the `global` default;
- projection-direction identity plus in-plane/translation updates for `fixed`;
- automatic and overridden local bounds;
- monotonic common frequency marching through `lpstop`;
- stochastic/full sampling and final update coverage;
- shared-memory and distributed execution;
- split-checkpoint handoff from `abinitio3D`;
- native-sampling final maps and expected artifacts.

Compilation and runtime tests are performed by the user. No Linux or BOX
result is recorded as passing without observed output.
