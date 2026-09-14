# Heterogeneity Analysis

Heterogeneity analysis extends the single-map model when one consensus density
cannot explain the particles. SIMPLE provides one continuous model and two
discrete-state workflows, selected by the provenance of the starting model and
particle poses:

- [`flex_pca`](flex_pca.md) estimates continuous variability with a fixed-pose,
  low-rank latent volume model and derives discrete state maps from that latent
  space.
- [`refine3D_states`](refine3d_states.md) refines same-lineage conformational
  states from an existing particle and reference scaffold.
- [`classify3D_refs`](classify3d_refs.md) classifies particles against a complete
  external reference set, then reconstructs data-derived state maps before
  continuing refinement.
- `ptcl3D_state_consensus` is a metadata utility in the same category: it
  combines the per-particle state assignments from a file table of projects
  (for example several `refine3D_states` runs) into one consensus `ptcl3D`
  field. It has no algorithm chapter.

## Shared discrete-state model

For `S` discrete conformations, particle `i` carries a state label `s_i` in
addition to its pose, and each state has a map `x_s`. This is the
[refine3D](../refine3d.md) alternation with the reference index extended to
`(state, direction)` pairs: compare the particle with projections of every map,
commit `(s_i, R_i, theta_i, shift_i)`, then [reconstruct](../reconstruction.md)
each map from its members.

Every state is compared with the same whitened Euclidean loss at the same
bandwidth, so the state decision is a likelihood-ratio test:

```text
s_i = argmin_s  min_{R,theta,shift} L_i(s, R, theta, shift).
```

In the probabilistic table the state is chosen by a deterministic argmin over
the heads of the balanced assignment loop
([sampling](../sampling_and_fractional_updates.md)); only
`refine=prob_state` draws the state from the full softmax
`exp(-(d_s - d_min))`. The projection within the chosen state is then drawn
stochastically as usual. A shift seed from one state is never reused to rank
another state.

## Shared frequency and coverage schedule

The discrete workflows split the iteration budget into blocks of three and
march the bandwidth linearly in Fourier index:

```text
k_start = max(5, index(lpstart)),      k_stop = min(box/2 - 2, index(lpstop)),
n_blocks = ceil(n_iterations / 3),
k_b = k_start + (b - 1) (k_stop - k_start) / (n_blocks - 2),
lp_b = max(resolution(k_b), lpstop).
```

With defaults `lpstart = 10 A` and `lpstop = 6 A`, the last two blocks run at
`lpstop`. A common band keeps state likelihoods comparable. Without an explicit
`maxits`, the cap is
`clamp(ceil(4 N_active / N_per_iteration), 10, 50)`, approximately four
expected updates per particle.

Fractional updates are balanced over projection-direction bins. A final pass
assigns any active particle with `updatecnt = 0` without changing the last
staged maps. A fresh all-particle reconstruction at native sampling then writes
state maps, half maps, FSCs, and orthogonal reprojections. Single-state
`combine_eo` finalization is not used because merging the half pair would break
the independence required for per-state FSCs.

The common bandwidth and balanced assignment keep the competing maps on equal
statistical footing; the workflow-specific initialization determines whether
that competition begins from a same-lineage scaffold or external references.
