# Conformational State Refinement Policy

This document defines the public and scientific contract of
`refine3D_states`. It is the canonical same-lineage multi-state workflow; no
historical command alias is registered or routed.

Related policies:

- [refine3D policy](../3D/refine3D_policy.md)
- [classify3D references policy](../3D/classify3D_refs_policy.md)
- [importance-sampling and fractional-update policy](../importance_sampling_fractional_update_policy.md)
- [nonuniform filtering policy](../NU/nonuniform_filtering_policy.md)

Primary implementation:

- `src/main/ui/simple/simple_ui_refine3D.f90`
- `src/main/exec/simple_exec_refine3D.f90`
- `src/main/commanders/simple/simple_commanders_refine3D.f90`
- `src/main/simple_refine3D_stage_plan.f90`
- `src/main/strategies/search/probabilistic/simple_strategy3D_prob.f90`

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

1. populated multi-state labels plus compatible project state maps, which
   continue without state initialization;
2. state-0/1 input plus `nstates`, initialized by `flex_pca` by default
   (`flex=yes`): labels and maps derive from the project consensus map under a
   population floor (`min_state_frac`, default 0.1) so that no
   under-populated flex cluster enters the volume refinement;
3. state-0/1 input plus `nstates` and `flex=no`, initialized by the
   distributed stochastic startup reconstruction;
4. a `solve3D` split checkpoint whose state maps are registered in the
   project `out` segment.

`flex` defaults to `yes`. It is skipped automatically when the project already
carries multi-state labels, and an explicit `flex=yes` on such a project is an
error. Flex initialization requires `nstates >= 3`; two-state stochastic
initialization needs `flex=no`.

`vol1..volN` input is rejected: starting state maps must come from the project
lineage. Classification against supplied references belongs to
`classify3D_refs`. Existing multi-state labels determine the effective state
count; an explicit `nstates` must agree. Every accepted state must be
populated.

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
fractional updates over the sampling units of `balance=cavg`: one unit per
selected 2D class, grouped into `nclust` (default 20) groups of similar class
averages; every draw shares its target equally over the groups, then over the
classes of a group, then lowest `updatecnt` first inside a class. Without
selected class averages in the project the wrapper falls back to
`balance=none` (global lowest `updatecnt` tiers). Particles inactive in
`ptcl3D` (state 0, for example those flex initialization dropped) are not part
of any unit. Otherwise it uses a full update. The sampler never uses the state
maps, poses or projection directions to choose particles.

With equal quotas one full visit of the active particles takes
`sweep = max over units of ceil(pop / quota)` draws, computed from the unit
table before the first stage. When states are initialized stochastically under
`local` or `global`, the `prob_state` init phase runs at least `sweep`
iterations and at most the larger of `sweep` and ten. The state-overlap exit
cannot end the phase before the sweep completes, so every active particle
receives an initial state label before `prob_neigh` refinement starts. The
automatic `maxits` is four target updates per particle times `sweep`, between
10 and 50 iterations.

Each `prob_neigh` frequency block refines one particle cohort:
`refine3D_states` sets the internal `cohort_sampling=yes` on the block's command
line, the block draws a new cohort at its first iteration and rescores that
cohort at the others, so state labelling and pose refinement converge on one
particle set per block, and `% PARTICLES UPDATED SO FAR` grows only at block
boundaries. A full visit therefore takes `sweep` blocks. Before the first stage
the wrapper prints the unit table with the expected visits per particle over
the planned blocks and warns when the least-visited unit falls below one visit
or the most-visited exceeds ten times the target; neither `nsample` nor the
frequency march is adjusted, and the terminal missing-update pass labels any
particle the march did not reach. The `prob_state` phase and the missing-update
pass draw every iteration. Holding one cohort for `k` iterations of trailing
reconstruction gives it a cumulative map coefficient `1 − (1−u)^k`; the
equal quota biases the composition of each partial reconstruction towards
about `min(pop, quota)` particles per class, which moves with `nsample`,
`nclust` and the class selection. Under cohorts `updatecnt` counts
iterations, not draws; coverage is read from `% PARTICLES UPDATED SO FAR`.
`sticky_class_sampling` is not an input of this workflow; the `solve3D`
docked handoff still sets it internally to keep the post-split cohort, with its
meaning unchanged.

`lpstart` and `lpstop` define one common frequency schedule for all states.
`simple_refine3D_stage_plan` returns short blocks containing the low-pass,
crop, translation limit, and global iteration range. Both
`refine3D_states` and `classify3D_refs` consume this planner. A state-specific
frequency schedule is outside the current contract because it would make
competitive evidence state-dependent.

Every planned frequency block is executed so the workflow reaches `lpstop`;
state-overlap diagnostics do not terminate the march at an earlier bandwidth.

## 5. `solve3D` Handoff

For docked multi-state `solve3D` work, `solve3D` owns the single-state
scaffold and split-checkpoint construction. The checkpoint preserves state
labels, maps, sampled/update metadata, the capped cohort, realized update
fraction, and iteration position. Post-split refinement is dispatched once to
`refine3D_states` with `pose_policy=local`; the checkpoint state maps reach it
through the project `out` segment, not as `vol1..volN` inputs.

The split checkpoint is constructed by
`simple_solve3D_split_checkpoint`; the old post-split state-refinement loop
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
reconstruction. Final state maps are then produced by the shared ending
`calc_final_rec` (module `simple_final_rec`), the same routine that closes
`solve3D`, `refine3D_auto`, and `classify3D_refs`: committed canonical
sigmas are reused when valid at native sampling, otherwise `bootstrap_rec3D`
rebuilds them, and the shipped maps are Euclidean ML reconstructions from all
active particles at native project sampling. The workflow writes normal state
volumes, half maps, FSC/resolution records, diagnostic low-pass outputs, and
orthogonal reprojections.

## 8. Validation

User-side validation must cover:

- all three pose policies and the `global` default;
- projection-direction identity plus in-plane/translation updates for `fixed`;
- automatic and overridden local bounds;
- monotonic common frequency marching through `lpstop`;
- stochastic/full sampling and final update coverage;
- shared-memory and distributed execution;
- split-checkpoint handoff from `solve3D`;
- native-sampling final maps and expected artifacts.

Compilation and runtime tests are performed by the user. No Linux or BOX
result is recorded as passing without observed output.
