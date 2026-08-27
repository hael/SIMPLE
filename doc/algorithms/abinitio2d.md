# Staged Ab Initio 2D Classification

## Problem

Construct stable 2D class averages and assignments from particles without
assuming trusted input references. The workflow must explore broadly at low
resolution, introduce probabilistic assignment only after classes become
meaningful, and finish with complete all-particle metadata.

## Stage controller

`abinitio2D` is a controller around repeated `cluster2D` invocations. It owns
the low-pass schedule, crop and translation policy, reference initialization,
search-mode transitions, outer sample size, and final coverage passes. The
matcher and classaverager retain their ordinary meanings at every stage.

The production schedule uses sampled stochastic-neighborhood search in stages
1 and 2. From stage 3 onward, `refine=prob` uses dense probabilistic candidate
tables. In `prob_snhc`, intermediate stages retain a sparse neighborhood table,
but the last staged call becomes dense `prob` so the previous class remains an
explicit candidate and class-overlap convergence can be measured.

## Initial references and frequency march

Initial class references are generated under the requested initialization
policy and begin at a coarse automatic low-pass limit. Each later stage admits
finer frequencies and may change crop and translation limits. A fixed public
`lp` override begins only when ML-regularized stages are active; the initial
Gaussian-reference stage keeps its automatic coarse limit so initialization
regularization remains meaningful.

The workflow is Cartesian. Polar Fourier transforms accelerate matching, but
there is no public alternative polar workflow branch.

## Outer sample target

The run-local target is

```text
N_target = nsample, default 200000
rho      = min(1, min(N_active,N_target)/N_active).
```

Stage 1 may use a random subset but never restores old class sums
fractionally. Later stages use update-count-biased sampling or reproduce the
subset chosen by probabilistic pre-alignment. When `rho<1`, previous class sums
are restored using class-local realized fractions.

The probabilistic pre-step chooses the outer subset exactly once. Table workers
and the matcher reproduce it; none may draw a second particle sample.

## Candidate sampling

For an active particle, the probability table evaluates an explicit top-K
class/in-plane support. Euclidean distances are treated as negative log
weights:

```text
p_j = exp[-(d_j-d_min)] / sum_l exp[-(d_l-d_min)].
```

Correlation mode uses `d=1-clamp(cc,0,1)` as a monotone pseudo-likelihood.
One candidate is sampled, optionally profile-refined, and written to the
assignment artifact. The matcher then performs the hard particle update. The
class average is not a soft weighted sum over all table candidates.

## Continuous in-plane polish

`inpl_cont=yes` replaces callback-style local angle/shift refinement with a
joint bounded `(sx,sy,theta)` solve after a class/cell decision. Probability
artifacts remain discrete: they store the rounded in-plane cell and its shift.
Only the durable hard assignment owns a fractional angle. Search-control
statistics use discrete cells so sub-grid polish cannot cause premature
convergence. The complete optimizer contract is in
[continuous in-plane refinement](continuous_inplane_refinement_abinitio2D.md).

## Coverage and final products

`fillin=yes` during staged work is a coverage guard, not a missing-only
sampler. After any sampled staged update, `abinitio2D` therefore runs a separate
dense greedy all-particle pass with fractional update disabled. It refreshes
class, angle, and shift for every active particle before final class averaging.

The workflow publishes final merged/even/odd class stacks, FRCs, ranked class
products, assignments, distances, sigma state, and updated project segments.

## Implementation

- Orchestration: `src/main/commanders/simple/simple_commanders_abinitio2D.f90`.
- Stage policy: `src/main/simple_abinitio2D_controller.f90`.
- Iteration execution: `src/main/strategies/parallelization/simple_cluster2D_strategy.f90`.
- Search and restoration: `src/main/strategies/search/simple_strategy2D_matcher.f90`
  and `src/main/class/simple_classaverager.f90`.
- Policy: `doc/policies/abinitio2D_policy.md`.

