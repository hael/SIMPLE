# Multi-State Refinement and Reference-Guided Classification

## Problem

A dataset contains `S` discrete conformations. Each particle carries an
unknown state label `s_i` in addition to its pose, and each state has its own
map `x_s`. The estimator is the [refine3D](refine3d.md) alternation with the
reference index extended to `(state, direction)` pairs: assign
`(s_i, R_i, theta_i, shift_i)` per particle by comparing against the
projections of all `S` maps, then reconstruct each map from its members.

The difficulty is not the estimator but its initialization and schedule. The
state label is discrete and the maps compete for particles, so an early
imbalance is self-reinforcing: a state with more particles gets a better map,
which wins more particles. Two workflows manage this for two kinds of input.

## Objective and state assignment

Every state is compared with the same whitened Euclidean loss at the same
bandwidth, so the per-particle state decision is a likelihood-ratio test:

```text
s_i = argmin_s  min_{R,theta,shift} L_i(s, R, theta, shift).
```

In the probabilistic table the state is chosen by a deterministic argmin over
the heads of the balanced assignment loop
([sampling](sampling_and_fractional_updates.md)); only `refine=prob_state`
draws the state from the full softmax `exp(-(d_s - d_min))` over states. The
projection within the chosen state is then drawn stochastically as usual.
No shift seed from one state is reused to rank another state.

## Conformational state refinement (`refine3D_states`)

Input: particles with a common pose scaffold (from a consensus refinement)
and either state references, an existing label set, a flex-PCA initializer,
or an `abinitio3D` docked checkpoint.

`pose_policy` fixes how much of the pose may move while states are being
decided:

- `fixed`: the projection direction is frozen at its input value; only the
  in-plane angle, shift, and state change. This is pure classification given
  the consensus geometry.
- `local`: the direction may move within the coarse Voronoi cell of its
  previous value. A user bound `alpha` in degrees is converted to a subspace
  size `nspace_sub = min(5000, max(2, 2/(1 - cos alpha)))`, the number of
  cells whose solid angle matches a cone of half-angle `alpha`.
- `global` (default): full coarse peak search across all directions and all
  states each iteration.

The sampling target is 10 000 particles per state, capped at 100 000 in
total; sampling is balanced over projection-direction bins so that every
view contributes to the state decision.

## Reference-guided classification (`classify3D_refs`)

Input: a complete external set of `S` volumes that may not share the
particles' pose history, for example maps from another dataset or from
docking. Comparing such references directly with the Euclidean objective is
unsafe because the noise model and the pose scaffold both assume references
built from these particles. The workflow therefore first builds
data-derived references:

```text
external maps, low-pass 15 A, nspace = 2500, at most 100 000 particles
  -> one greedy correlation-objective pass assigning (state, pose)
  -> noise power estimated from the residuals of that pass
  -> reconstruction of checkpoint maps from the assigned particles
  -> Euclidean probabilistic multi-state refinement from the checkpoints.
```

The correlation objective is scale-free, so it tolerates references whose
amplitude spectrum does not match the data. The external maps are never
modified; from the checkpoint onward the data-derived maps are the
references. Coverage is verified: at least 99 percent of the sampled
particles must have been updated and every state must be non-empty.

## Frequency schedule

Both workflows split the iteration budget into blocks of three iterations
and march the bandwidth linearly in Fourier index:

```text
k_start = max(5, index(lpstart)),      k_stop = min(box/2 - 2, index(lpstop)),
n_blocks = ceil(n_iterations / 3),
k_b = k_start + (b - 1) (k_stop - k_start) / (n_blocks - 2),
lp_b = max( resolution(k_b), lpstop ),
```

with defaults `lpstart = 10 A`, `lpstop = 6 A`, so the last two blocks run at
`lpstop`. Every state uses the same band in a block, which is what makes the
state likelihoods comparable. The iteration cap without an explicit `maxits`
is `clamp(ceil(4 N_active / N_per_iteration), 10, 50)`, about four expected
updates per particle.

## Coverage and final maps

Fractional updates are balanced over bins. A final pass assigns any active
particle with `updatecnt = 0` without touching the last staged maps. Only then
is a fresh all-particle reconstruction made at native sampling, producing
state maps, half maps, FSCs, and orthogonal reprojections. The
`combine_eo` final iteration of single-state refinement is not used, since a
merged half pair would break the independence needed for per-state FSCs.

## Rationale

- Deciding states with a likelihood ratio at a common bandwidth is the
  discrete analogue of the pose decision; nothing in the estimator changes,
  only the index set.
- Balanced sampling over views and the balanced assignment loop are the two
  mechanisms that keep the competition from collapsing to one state.
- The correlation pre-pass for external references is a change of objective,
  not a change of estimator: it produces the noise model and the pose scaffold
  that the Euclidean objective requires, from the data at hand.

## Implementation

- Workflows: `src/main/commanders/simple/simple_commanders_refine3D.f90`.
- Correlation pose initialization:
  `src/main/simple_external_reference_pose_initialization.f90`.
- Frequency schedule: `src/main/simple_refine3D_stage_plan.f90`.
- State assignment: `src/main/simple_eul_prob_tab.f90`,
  `src/main/strategies/search/simple_strategy3D_prob.f90`.
- Policies: `doc/policies/refine3D_states_policy.md`,
  `doc/policies/classify3D_refs_policy.md`.
