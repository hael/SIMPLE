# Probabilistic Pre-Alignment and Fractional Updates

## Problem

Reduce the cost of large 2D and 3D iterations without confusing two different
sampling problems: which particles are updated, and which pose/class candidates
are explored for each updated particle.

## Two sampling layers

The **outer sampler** chooses the particle subset for the iteration. It records
the current round in `sampled` and increments cumulative `updatecnt`. Balanced,
random, update-count-biased, class-balanced, state-balanced, and fill-in paths
are alternative policies for this layer.

The **inner sampler** operates only inside that fixed subset. A probabilistic
table evaluates a bounded support of reference, orientation, state, or in-plane
candidates and draws candidates from score-derived weights.

The layers are independent:

```text
outer subset -> probability-table candidates -> one assignment -> hard update.
```

Candidate sampling may not add, remove, or resample particles.

## Sample-once-and-reuse

For a probabilistic iteration:

1. the pre-alignment commander chooses the outer subset;
2. the project records its `sampled` generation;
3. table workers evaluate only that generation;
4. the table outputs merge into one assignment artifact;
5. the matcher reproduces the same subset;
6. hard updates and partial reconstruction/averaging are written.

`sample4update_reprod` is the reproduction contract. Inferring participation
from nominal `update_frac`, or drawing a new sample in a worker, changes the
estimator.

## Importance weights

For evaluated candidate costs `d_j`, SIMPLE uses stabilized exponential
weights over the retained support:

```text
w_j = exp[-(d_j-d_min)],       p_j = w_j / sum_l w_l.
```

The support is intentionally truncated. It is the local search distribution,
not a claim to represent the full posterior over all SO(2) or SO(3). After one
candidate is drawn, shifts and sometimes in-plane angle may be profile/MAP
refined; the result is still one hard particle assignment.

## Realized update fractions

The requested `update_frac` is a target. Restoration uses realized membership
from `sampled`, active labels, and `updatecnt>0`:

- `rho_k` for each 2D class;
- `rho_s` for each 3D state;
- a global fraction only for callers whose object is genuinely global.

This matters when balancing and integer quotas give different classes or states
different realized fractions.

## 2D restoration

For class `k`, prior even/odd numerator and CTF-squared sums receive weight
`1-rho_k`, while current partial sums retain their current-update mass. Empty
updated classes therefore keep their previous estimate; fully updated classes
replace it. Stage 1 of `abinitio2D` intentionally starts without this memory.

## 3D trailing reconstruction

The persistent state/half chain stores unregularized Fourier numerators and
sampling densities at full-dataset mass. Let `f` be the realized fraction that
produced the current partials and `u` the desired map-update weight (`u=f`
unless an explicit single-state override exists). The accumulator recurrence is

```text
A_new = (u/f) A_current + (1-u) A_previous.
```

Thus current mass is `uD`, previous mass is `(1-u)D`, and one density correction
after the blend restores the map. The manifest and all state/half components
form one validated artifact set; invalid provenance causes re-seeding, not a
partial reuse.

## Guards

- Full-update mode selects every active particle and disables trailing.
- Probabilistic artifacts retain rounded in-plane cells; durable fractional
  pose precision is committed by the matcher afterward.
- Restoration and volume assembly consume sampling state but never choose a
  new subset.
- Final scientific maps and final class products require an all-particle
  coverage pass when staged fractional updates left any active particle unseen.

## Implementation

- Sampling helpers: `src/main/strategies/search/simple_matcher_smpl_and_lplims.f90`.
- Bookkeeping: `src/main/ori/simple_oris_sampling.f90` and
  `src/main/ori/simple_oris_getters.f90`.
- Probability tables: `src/main/simple_eul_prob_tab*.f90`.
- Orchestration: `src/main/commanders/simple/simple_commanders_prob.f90`.
- Policy: `doc/policies/importance_sampling_fractional_update_policy.md`.

