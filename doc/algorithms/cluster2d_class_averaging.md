# Cluster2D and CTF-Aware Class Averaging

## Problem

Partition particle images into aligned 2D classes and estimate one
CTF-corrected class average per class without mixing particle-domain assignment
with class-domain restoration.

## Particle pose and objective

For particle `i`, a 2D assignment is `(class, theta, sx, sy)`. Current class
averages are transformed to polar Fourier representations. Search strategies
evaluate a discrete class/in-plane support and refine shifts according to the
selected objective:

- Euclidean mode uses noise-normalized distances and updates the grouped noise
  estimate;
- correlation mode maximizes normalized correlation;
- probabilistic modes consume a precomputed bounded candidate assignment, then
  still commit one hard class and pose.

Search modes differ in how they explore the candidate support—greedy,
stochastic-neighborhood hill climbing, sampled variants, in-plane-only, or
probabilistic—but share the same durable particle fields.

## Iteration

One `cluster2D` iteration:

1. chooses or reproduces the outer particle subset;
2. reads current class references and FRC metadata;
3. prepares Fourier references and the active frequency band;
4. reads particles in bounded batches;
5. searches one assignment per active particle;
6. updates class, in-plane angle, shift, score, and optional sigma state;
7. accumulates even/odd class numerator and CTF-squared sums from the same
   in-memory batch;
8. restores merged class averages and calculates convergence statistics.

The online path reads each batch once for both matching and accumulation.
Distributed workers write partial sums; assembly performs the same scientific
restoration after deterministic reduction.

## Fourier accumulation and restoration

Each aligned particle Fourier plane contributes to its assigned class through
the 2D Kaiser-Bessel interpolation kernel. For class `k`, the accumulator keeps
an even and odd numerator `B_k(q)` and CTF-squared density `D_k(q)`. Conceptually,

```text
B_k(q) = sum_{i in k} w_i(q) conj(CTF_i(q)) Y_i^aligned(q)
D_k(q) = sum_{i in k} w_i(q) |CTF_i(q)|^2.
```

Restoration divides the numerator by the regularized density, applies the
gridding correction and masks, and transforms back to real space. Even/odd
products remain available for FRC; the merged class is their derived
combination.

## FRC and search bandwidth

For each class, Fourier ring correlation compares its even and odd averages.
The stored curves provide 0.5 and 0.143 resolution estimates and drive
class-specific filtering and later alignment-band choices. Degenerate curves
are bounded to the first wavelength or Nyquist rather than emitting undefined
resolution values.

## Fractional restoration

When only a subset is updated, prior class sums are carried forward per class,
not through one global fraction. If `rho_k` is the realized updated fraction of
class `k`, restoration combines

```text
new accumulator = current partial + (1-rho_k) * previous accumulator,
```

with even/odd numerator and CTF-squared sums kept separate. A class with no
updated members retains its previous contribution; a fully updated class
replaces it. Stage 1 of `abinitio2D` deliberately disables this carry-over.

## Implementation

- Workflow and search: `src/main/commanders/simple/simple_commanders_cluster2D.f90`,
  `src/main/strategies/parallelization/simple_cluster2D_strategy.f90`, and
  `src/main/strategies/search/simple_strategy2D_matcher.f90`.
- Class sums and restoration: `src/main/class/simple_classaverager.f90`.
- FRCs: `src/main/class/simple_class_frcs.f90`.
- Policy: `doc/policies/abinitio2D_policy.md` and
  `doc/policies/importance_sampling_fractional_update_policy.md`.

