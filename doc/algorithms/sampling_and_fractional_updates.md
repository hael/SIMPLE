# Probabilistic Pre-Alignment and Fractional Updates

## Problem

An iteration of 2D or 3D alignment costs `O(N x R x nrots)` for `N`
particles and `R` references. Two different subsampling ideas reduce it, and
they must not be confused:

- the **outer** problem: which particles are updated this iteration;
- the **inner** problem: for an updated particle, which of the `R x nrots`
  candidate poses are evaluated and how one is chosen.

The outer problem is a stochastic-gradient idea: a fraction of the data gives
an unbiased but noisier update of the references. The inner problem is
importance sampling over poses. This chapter covers both, and the rule that
blends a partial update into the running reference estimate.

## Outer sampling

Let `N_active` be the particles with state above zero and let
`update_frac` be the requested fraction. The subset size is

```text
n = min(N_active, max(1, round(update_frac * N_active))).
```

Each particle carries two counters: `updatecnt`, the number of times it has
been updated, and `sampled`, the index of the last iteration that updated it.
The policies for choosing the `n` particles are:

- **random**: uniform without replacement.
- **count-tiered** (the default for fractional 3D and later 2D stages): take
  every particle at the lowest `updatecnt`, then the next tier, until the
  budget is filled; within the tier at the cutoff, draw uniformly. This is
  lowest-count-first rather than a weighted draw, so coverage of the dataset
  is as even as possible.
- **fill-in**: the same, starting explicitly at `updatecnt = 0`.
- **class-balanced** (`balance=yes`): distribute the budget across classes
  (2D) or projection-direction bins (3D) by round-robin increments capped by
  each bin's population, so every bin gets an equal share up to its size;
  then apply the count-tiered rule within each bin. With `frac_best < 1` only
  the best-scoring fraction of each bin is eligible, and `greedy_sampling`
  takes the top of each bin deterministically.
- **reproduce**: reuse exactly the particles with `sampled` equal to the
  current maximum. This is what lets a probabilistic iteration choose the
  subset once and have every subsequent step, table filling, matching, and
  reconstruction, act on the same particles.

A run whose target `nsample` exceeds 90 percent of the active set switches to
full updates and disables the blending below; the gain would not justify the
added noise.

## Inner sampling: the probability table

For a probabilistic iteration, before any hard assignment, a table of losses
is computed. For each active particle `i` and each reference `j` (in 3D,
`j` enumerates `nstates x nspace` projections; in 2D, the active classes):

1. a common shift seed is found by shift-only optimization at the particle's
   previous pose;
2. every reference is scored at that shift over all in-plane rotations, and
   the in-plane rotation is drawn from the truncated softmax over the best
   `K_inpl` rotations;
3. the best `K_ref` references are re-refined in shift.

The truncation `K` is derived from an angular threshold rather than fixed:

```text
K = min(n, max(1, floor(athres * n / 180))),
athres = min(prob_athres, mean angular change of the last iteration),
```

with `prob_athres` defaulting to 10 degrees. As the alignment settles, the
measured angular change shrinks and the support tightens automatically.

Losses are whitened negative log-likelihoods (`d = -log(score)` for the
Euclidean objective, `d = 1 - max(cc, 0)` for correlation). A draw from a
support of `K` candidates uses

```text
w_j = exp[-(d_j - d_min)],     p_j = w_j / sum_l w_l.
```

No temperature is applied: the noise normalization of the objective already
puts `d` in natural log-likelihood units. The support is a local search
distribution, not an approximation to the full posterior over SO(3); the
method commits one candidate and is not a marginalizing EM.

## Balanced assignment

The table is not consumed particle by particle. It is consumed as one global
assignment problem, which is what prevents a few references from absorbing
most particles:

1. for each reference `j`, sort all particles by `d_ij`;
2. each reference keeps a *head*, its best still-unassigned particle;
3. repeat until every particle is assigned: draw one reference from the
   softmax over the head losses (truncated to the best `K` heads), give it its
   head particle, mark the particle taken, advance all heads.

Every reference competes at every step and a reference that wins takes exactly
one particle, so populations equalize implicitly without quotas. The result is
a stochastic greedy bipartite matching rather than an optimal (Hungarian)
one; optimality is not wanted, since the draw is the exploration mechanism.

In multi-state 3D, the state label is assigned first by the same loop with a
deterministic argmin over heads (only `refine=prob_state` samples the state),
and the within-state projection is then drawn stochastically. Neighborhood
variants restrict which references are scored per particle: a stochastic
subset (`shc`, `snhc`), the coarse Voronoi cell containing the previous
projection (`geom`), or the pooled top-`npeaks` coarse cells across states
(`state`). The assignment loop is unchanged.

## Blending partial updates

The requested `update_frac` is a target; restoration uses the realized
fraction, computed from `sampled` and `updatecnt`:

```text
rho_k = #{updated and active in k} / #{active in k},
```

for each 2D class `k`, or each 3D state.

**2D.** Class accumulators are carried forward with weight `1 - rho_k`
([Cluster2D](cluster2d_class_averaging.md)).

**3D trailing reconstruction.** A persistent chain stores, per state and half,
the unregularized Fourier numerator and sampling density at full-dataset mass.
With `f` the realized fraction that produced the current partial sums and `u`
the desired update weight (`u = f` unless overridden), the blend in the
accumulator domain is

```text
A_new = (u/f) A_current + (1 - u) A_previous,
```

applied identically to numerator and density. The factor `u/f` rescales the
partial sums so the current data carry weight `u` while total sampling mass is
conserved: `(u/f)(f D) + (1-u) D = D`. One density division after the blend
therefore restores a correctly normalized map. Blending in the accumulator
domain rather than between restored volumes keeps the FSC computed on the
blended halves honest, because both halves are still ratios of sums.

Because the exponential moving average `A_new` decays previous contributions
geometrically, the reference at iteration `t` is a weighted sum over the last
`~1/u` iterations of partial reconstructions. This is why fractional 3D
refinement converges at all: each partial reconstruction alone would be too
noisy to align against.

## Guards

- Full-update mode selects every active particle and disables trailing.
- Restoration and volume assembly read `sampled` but never choose a subset.
- Final maps and final class averages require an all-particle coverage pass
  when staged fractional updates left any active particle unseen.

## Implementation

- Sampling policies: `src/main/ori/simple_oris_sampling.f90`;
  dispatch in `src/main/strategies/search/simple_matcher_smpl_and_lplims.f90`.
- Probability tables and assignment: `src/main/simple_eul_prob_tab*.f90`.
- Trailing blend: `src/main/commanders/simple/simple_commanders_rec_distr.f90`.
- Orchestration: `src/main/commanders/simple/simple_commanders_prob.f90`.
- Policy: `doc/policies/importance_sampling_fractional_update_policy.md`.
