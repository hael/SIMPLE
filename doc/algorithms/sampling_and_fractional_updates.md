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
importance sampling over poses. This chapter covers both, the assignment
that consumes the inner samples, and the rule that blends a partial update
into the running reference estimate.

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
   one in-plane rotation is drawn from the row's tempered posterior over all
   of them (below);
3. the best `K_ref` references are refined again in shift.

For the Euclidean objective, the loss `d` is the whitened negative
log-likelihood, `d = -log(score)`, normalized by the particle's weighted
power: one nat is `1/T_i` of it with

```text
T_i = pi wsqsum_i / (2 pftsz),
```

the particle's noise temperature, several hundred to a thousand on PfCRT
(measured 2026-10-22, doc/implementation_notes/planned/pose_marginalization.md,
sections 3 and 14). A draw from `n` candidates uses

```text
w_j = exp[-beta T_i (d_j - d_min)],     p_j = w_j / sum_l w_l,
```

the tempered posterior p^beta with `beta = prob_beta` (default 0.3, for
small membrane proteins; 1 is the physical posterior; below one hotter). No
truncation: the exploration of the draw is the data's width, not a fixed
angular window. The tables refuse the correlation objective. One candidate
is drawn and committed as a hard assignment. A particle is never spread
over several orientations in proportion to their probabilities, as it is
in expectation-maximization refinement.

The randomness of the in-plane draw is load-bearing: it is what keeps the
direction assignment diverse against the references' memory of where each
particle was last put (the self-reinforcement of the same note, 1-7 nats
on the own half). Scoring every entry at its best grid rotation instead
(`K_inpl = 1`, tried as a switch on 2026-10-23 and removed again) takes the
search noise out as intended but traps the run: 5 of 10 PfCRT restarts
ended near 7 A against 4.0-4.6 A for the other 5 and 4.0 A for the draw
(`PfCRT_posterior/greedy_inpl`). The in-plane angle a particle ends with is
the continuous refinement's in either case; what the draw buys is the
diversity of the direction assignment.

Before the posterior (until 2026-10-24) the in-plane rotation was drawn by
the untempered softmax `exp[-(d_j - d_min)]` over the best `K_inpl`
rotations, the truncation derived from an angular threshold,
`K = min(n, max(1, floor(athres * n / 180)))` with
`athres = min(prob_athres, mean angular change of the last iteration)` and
`prob_athres` 10 degrees: since `T_i` is hundreds, a uniform draw over a
random rotation within `athres` of the best (on PfCRT more than a nat worse
than the best grid pose for 17-27% of the particles; the note, section 14),
which `athres` shrinking with the measured in-plane change annealed as the
alignment settled. The tables' fills no longer use it; the SHC strategies
of the early 2D stages keep their own samplers.

## Assignment

The table is consumed particle by particle (2026-10-23; the only assignment
since 2026-10-24, when the class-balanced frontier was deleted from every
table): each particle's reference is one draw from its own row by the
tempered posterior above, within the state of its closest reference, `T_i`
written by the table worker into the part file. Particles are assigned to
references; no reference is forced to be populated, and no particle takes
an argmin. The stage then ends on the resolution exit or its budget, not on
the pose overlap (the committed pose is a draw and stable from a stage's
first iterations). The 2D class assignment is the same per-particle draw
over the row (`solve2D_policy.md`); the neighbourhood tables follow the rule
on their candidate store (the state of the closest valid candidate, one
draw within it), the current pose among the candidates so that the draw is
an acceptance step and not a forced move.

The temperature is the band's: T_i grows with the number of Fourier
components the stage searches, so at solve3D's stages 1-2 (rings 2-5 on a
137 A box) T is about 33 where it is hundreds from stage 3 on, and the
physical posterior is flat there and a delta later; the ladder is itself a
cooling schedule. In ten streptavidin restarts at beta 0.3, with the current
pose missing from the shc neighbourhood's candidates, stages 1-2 did nothing
and stage 3 solved the structure from scratch (the note, section 19); the
incumbent is a candidate since, and whether the early stages also need a
colder draw is the open measurement. `prob_beta` is the
search's annealing knob, since at the physical temperature the posterior
over directions is a delta once the map has any resolution (n_eff 1.0-1.2)
and the map-pose alternation then stops at the nearest fixed point (section
17 of the note). A refine3D run with `prob_beta_end` anneals geometrically
from `prob_beta` to it over its iterations; solve3D does so on its last
stage, to 1, so the final map is weighted at the physical posterior. The
design and its measurements are in
`doc/implementation_notes/planned/pose_marginalization.md`.

### The balanced frontier, deleted

Until 2026-10-24 the tables consumed the table as one global assignment
problem: for each reference, sort its particles by `d_ij`; each reference
keeps a head, its best still-unassigned particle; until every particle is
assigned, draw one reference from the softmax over the head losses
(truncated to the best `K` heads), give it its head particle and advance
the heads. Every reference competed at every step and a winning reference
took exactly one particle, so populations evened out without quotas: a
stochastic greedy matching between particles and references, which forced
every reference to be populated whether or not any particle's data put it
there. Hans, 2026-10-24: methodologically it makes no sense; particles are
assigned to references, not references to particles. It was deleted from
the dense 3D table, the 2D tables and the neighbourhood table, with its
samplers and reference tests (the note, section 19).

Neighborhood variants (`prob_neigh`) restrict which references are scored
for each particle: a stochastic subset (`shc`, `snhc`), the coarse cell of
projection directions that contains the previous projection (`geom`), or
the best `npeaks` coarse cells pooled across states (`state`); the cells are
defined in [refine3D](refine3d.md). The assignment is the same in every
case.

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
