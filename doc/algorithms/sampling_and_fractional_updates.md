# Sampling and Fractional Updates

## Problem

An iteration of 2D or 3D alignment costs `O(N x R x nrots)` for `N`
particles and `R` references. Two different subsampling ideas reduce that
cost, and they must not be confused:

- the **outer** problem: which particles are updated in this iteration;
- the **inner** problem: for an updated particle, which of the `R x nrots`
  candidate poses are evaluated, and how one of them is chosen.

The outer problem is a stochastic-gradient idea: a fraction of the data gives
a cheaper but noisier update of the references. The inner problem is
importance sampling over poses. This chapter covers both, the balanced
assignment that consumes the inner samples, and the rule that blends a
partial update into the running reference estimate.

## Outer sampling

Let `N_active` be the number of active particles (state above zero) and
`update_frac` the requested fraction. The subset size is

```text
n = min(N_active, max(1, round(update_frac * N_active))).
```

Each particle carries two counters: `updatecnt`, the number of times it has
been updated, and `sampled`, the last sampling round that selected it. The
`n` particles are chosen by one of these policies:

- **random**: uniform without replacement. The first stage of `solve2D` draws
  its subset this way once and keeps it for all of the stage's iterations.
- **count-tiered** (`balance=none`, the default of `refine3D` and of the
  later `solve2D` stages): take every particle at the lowest `updatecnt`,
  then the next tier, until the budget is filled; within the tier at the
  cutoff, draw uniformly. This is lowest-count-first rather than a weighted
  draw, so coverage of the dataset is as even as possible.
- **fill-in** (late 3D stages): the same, but always starting from the
  particles that have never been updated (`updatecnt = 0`).
- **class units** (`balance=class` or `balance=cavg`; `cavg` is the default
  of `solve3D`, `refine3D_states`, and `classify3D_refs`): share the budget
  over sampling units built from the selected 2D classes. Each unit is one
  selected class, with its active particles ordered by their 2D score.

  - Under `class`, every unit gets an equal share, capped at its population.
  - Under `cavg`, the class averages are first clustered into `nclust`
    groups (default 20) by average linkage on their rotation-, shift-, and
    mirror-invariant correlation, so that a preferred view split over many
    classes becomes one group. The budget is then shared equally over the
    groups and, within a group, equally over its units; the remainder of a
    split goes to the least-updated units. The group level keeps the views
    balanced, and the unit level keeps the classes inside a view, which may
    be different conformational states, on an equal footing.

  In both modes the count-tiered rule applies inside each unit. With
  `frac_best < 1`, only the best-scoring fraction of each unit is eligible,
  and `greedy_sampling` takes the top of each unit deterministically.
- **reproduce**: reuse exactly the particles whose `sampled` equals the
  current maximum. This lets a probabilistic iteration choose the subset once
  and have every later step (table filling, matching, and reconstruction) act
  on the same particles.
- **cohort** (`refine3D_states`): the set drawn at the first iteration of a
  frequency block is kept for the whole block. Later iterations reproduce it
  and advance its counters, so that state labels and poses converge on one
  set.
- **missing**: every active particle that has never been updated. This is the
  terminal coverage pass.

Nothing derived from the 3D maps, the poses, or the projection directions
decides which particles are sampled. The maps carry bias, and a sampler
driven by them would feed that bias back; the 2D class averages are evidence
independent of the maps.

With equal quotas, visiting every particle once takes
`max over units of ceil(pop / quota)` iterations rather than
`1 / update_frac`, because a unit larger than its share takes longer. This
sweep length is printed in a table of the units before the first stage, and
nothing is adjusted automatically to meet it.

In `solve3D`, a target `nsample` above 90 percent of the active particles
switches to full updates (`solve2D` switches at 99 percent). Full updates
select every active particle and disable the blending described below: that
close to the full set, subsampling would save little and add noise.

## Inner sampling: the probability table

Before any hard assignment, a probabilistic iteration fills a table of losses
between every particle in the subset and every reference `j`: in 3D the
`nstates x nspace` projections, in 2D the active classes. For each particle:

1. a common shift seed is found by shift-only optimization at the particle's
   previous pose;
2. every reference is scored at that shift over all in-plane rotations, and
   the in-plane rotation is drawn from the best `K_inpl` rotations by the
   softmax defined below;
3. the best `K_ref` references are refined again in shift.

For a list of `n` candidates, the truncation `K` is derived from an angular
threshold rather than fixed:

```text
K      = min(n, max(1, floor(athres * n / 180))),
athres = min(prob_athres, mean angular change of the last iteration),
```

with `prob_athres` defaulting to 10 degrees. As the alignment settles, the
measured angular change shrinks and the candidate list tightens
automatically.

For the Euclidean objective, the loss is the whitened negative
log-likelihood, `d = -log(score)`; for correlation, it is
`d = 1 - max(cc, 0)`. A draw from `K` candidates uses

```text
w_j = exp[-(d_j - d_min)],     p_j = w_j / sum_l w_l.
```

No temperature is applied, because the noise normalization of the objective
already puts `d` in natural log-likelihood units. The `K` candidates are only
the best-scoring poses near the top of the list, used to decide where to
search; they are not meant to represent the probability of every possible
orientation. One candidate is drawn and committed as a hard assignment. A
particle is never spread over several orientations in proportion to their
probabilities, as it is in expectation-maximization refinement.

## Balanced assignment

The table is not consumed particle by particle. It is consumed as one global
assignment problem, which is what prevents a few references from absorbing
most of the particles:

1. for each reference `j`, sort all particles by `d_ij`;
2. each reference keeps a *head*, its best still-unassigned particle;
3. until every particle is assigned, repeat: draw one reference from the
   softmax over the head losses (truncated to the best `K` heads), give it
   its head particle, mark the particle as taken, and advance all heads.

Every reference competes at every step, and a reference that wins takes
exactly one particle, so populations even out without explicit quotas. The
result is a stochastic greedy matching between particles and references, not
the optimal assignment that the Hungarian algorithm would find. Optimality is
not wanted here, because the draw is the exploration mechanism.

In multi-state 3D, the state label is assigned first, by the same loop with a
deterministic choice of the lowest head loss, and the projection within the
state is then drawn
stochastically. Neighborhood variants (`prob_neigh`) restrict which
references are scored for each particle: a stochastic subset (`shc`,
`snhc`), the coarse cell of projection directions that contains the previous
projection (`geom`), or the best `npeaks` coarse cells pooled across states
(`state`); the cells are defined in [refine3D](refine3d.md). The assignment
loop is the same in every case.

## Blending partial updates

The requested `update_frac` is only a target. Restoration uses the realized
fraction of each 2D class or 3D state `k`,

```text
rho_k = n_k / N_k,
```

where `n_k` counts the active members sampled in the current round and `N_k`
the active members that have been updated at least once (`updatecnt > 0`).

**Population rule.** Every carried set records `M_k`, the population its sums
represent. With `u` the applied update weight (`u = rho_k` unless
overridden), the assembly step blends

```text
A_new = (u/rho_k) A_current + (1 - u) (N_k/M_k) A_previous,    M_k <- N_k,
```

so the sampling mass after the blend is `u N_k + (1 - u) N_k = N_k`, the
represented population, whatever joined or left the group. With an unchanged
population (`M_k = N_k`), this reduces to the familiar
`(u/rho_k) A_current + (1 - u) A_previous`. The rule keeps the mass right,
not the membership: old contributions are removed in proportion, not
particle by particle.

**2D.** Class accumulators are carried forward with weight `(N_k - n_k)/M_k`,
the population rule with `u = rho_k`
([Refine2D](refine2d_class_averaging.md)).

**3D trailing reconstruction.** A persistent chain stores, per state and
half, the unregularized Fourier numerator and sampling density at the mass of
the population it represents. With `f` the realized fraction of the state
(`rho_k` above) and `u` the desired update weight (`u = f` unless
overridden), the blend in the accumulator domain is

```text
A_new = (u/f) A_current + (1 - u) (N/M) A_previous,
```

applied identically to the numerator and the density. The factor `u/f`
rescales the partial sums so that the current data carry weight `u` while
the total sampling mass is that of the represented population:
`(u/f)(f D) + (1-u)(N/M) M d = D` for a per-particle density `d` and
`D = N d`. One density division after the blend therefore restores a
correctly normalized map. Blending the accumulators rather than restored
volumes keeps the FSC computed on the blended halves valid, because both
halves are still ratios of sums; finished, restored half maps are never
blended.

When a blend is due but no chain exists yet, the iteration does not blend.
The chain is seeded with the current partial sums scaled by `1/f`, so that it
carries the mass of the full population, and the current sample alone is
restored and used for that iteration. The chain then fills over about `1/f`
iterations. The override of `u` (`ufrac_trec`) applies to single-state runs
only; multi-state runs always use the realized fraction of each state.

Because the blend is an exponential moving average, previous contributions
decay geometrically, and the reference at iteration `t` is a weighted sum of
the partial reconstructions of roughly the last `1/u` iterations. This is
why fractional 3D refinement converges at all: each partial reconstruction
on its own would be too noisy to align against.

## Guards

- Full-update mode selects every active particle and disables trailing.
- Restoration and volume assembly read `sampled` but never choose a subset.
- Final maps and final class averages require coverage when staged
  fractional updates left any active particle unseen: `solve2D` ends with a
  dense all-particle pass, and 3D runs end with the missing-update pass.

## Implementation

- Sampling policies: `src/main/ori/simple_oris_sampling.f90`;
  dispatch in `src/main/strategies/search/simple_matcher_smpl_and_lplims.f90`.
- Sampling units, class-average groups, and the unit table:
  `src/main/strategies/search/simple_view_partition_sampling.f90`.
- Probability tables and assignment:
  `src/main/strategies/search/probabilistic/simple_eul_prob_tab*.f90`.
- Trailing blend: `src/main/commanders/simple/simple_commanders_rec_distr.f90`;
  the distributed PCG chain in
  `src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90`.
- Orchestration: `src/main/commanders/simple/simple_commanders_prob.f90`.
- Policy: `doc/policies/importance_sampling_fractional_update_policy.md`.
