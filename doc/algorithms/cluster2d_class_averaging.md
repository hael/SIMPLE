# Cluster2D and CTF-Aware Class Averaging

## Problem

Given `N` particle images and a number of classes `K`, find for each particle a
class label `c_i`, an in-plane rotation `theta_i`, and a shift `(sx_i, sy_i)`,
and for each class a CTF-corrected mean image `A_k`, such that each particle
is well explained as a rotated, shifted, CTF-modulated copy of its class mean.

The problem is a joint clustering and registration problem with a
non-convex objective, so it is solved by alternation: fix the class means and
update the assignments (the **search** step), then fix the assignments and
update the class means (the **restoration** step). Everything below is one
round of that alternation. The schedule that wraps rounds into a
coarse-to-fine run is described in [ab initio 2D](abinitio2d.md).

## Model

The image `X_i` is compared with the class mean `A_k` in polar Fourier
coordinates `(phi, k)`: `k` indexes rings (spatial frequency) and `phi` indexes
`nrots` equispaced angles. An in-plane rotation is then a cyclic shift of the
angular index, and a translation is a phase ramp. The forward model is

```text
X_i(phi,k) = CTF_i(phi,k) * [S_s R_theta A_c](phi,k) + n_i(phi,k),
```

with `n_i` Gaussian, independent across Fourier samples, with per-ring
variance `sigma2_i(k)`. The particle is never CTF-corrected; the CTF is
applied to the reference, which is the numerically stable direction.

**Euclidean objective** (the default). The normalized negative log-likelihood of
particle `i` at rotation `r`, shift `s`, class `c` is

```text
L_i(c,r,s) = sum_k w_k sum_phi |X_i - CTF_i [S_s R_r A_c]|^2
             / sum_k w_k sum_phi |X_i|^2,          w_k = k / sigma2_i(k).
```

The ring weight `k` accounts for the number of Fourier samples per ring, and
`1/sigma2` whitens the noise. Dividing by the weighted particle energy makes
`L=1` for an empty reference and `L=0` for a perfect fit, so values are
comparable across particles. The score used for ranking is `exp(-L)`.

**Correlation objective** (`objfun=cc`) is the scale-free alternative

```text
cc_i(c,r,s) = <X_i, CTF_i S_s R_r A_c> / (||X_i|| ||CTF_i S_s R_r A_c||),
```

used when no noise model is available yet.

Both objectives are evaluated for all `nrots` rotations at once: the sum over
`phi` of a product of two angular functions is a circular correlation, so one
inverse FFT along the angular index yields every rotation. Rings outside the
current band `[k_lo, k_hi]` are excluded; `k_hi` is the low-pass limit set by
the stage schedule or, once class FRCs exist, by the median of the three best
class resolutions at FRC 0.143.

**Noise model.** `sigma2_i(k)` is estimated from the residual at the
committed pose,

```text
sigma2_contrib_i(k) = (1/nrots) sum_phi |CTF_i [S R A_c] - X_i|^2,
```

averaged over all particles in the same group (by default one group per
micrograph stack, separately for even and odd halves), and lagged by one
iteration. Before any alignment exists it is bootstrapped from the background
power spectrum outside the particle mask.

## Algorithm: the search step

For each active particle the previous pose is the seed. The class candidates
are visited in a fresh random permutation with the previous class placed last,
so that a random subset of classes is always explored before the old answer is
reconsidered. The search modes differ in how much of that permutation is
visited and how the winner is chosen:

- **greedy**: evaluate every non-empty class, take the argmax over all
  rotations.
- **snhc** (stochastic neighborhood hill climbing): visit at most
  `n_bound = max(2, K (1 - f))` classes; for each, pick uniformly at random a
  rotation whose score beats the previous score and stop at the first class
  where one exists. If none beats it, keep the best seen. The searched
  fraction `1 - f` is annealed:

  ```text
  f = clamp( 0.5 * 0.8^(15 it/extr_lim - 2), 0, 0.5 ),
  ```

  so about half the classes are searched early and all of them once
  `it > extr_lim` (default 15). Iteration 1 is forced greedy.
- **snhc_smpl**: as snhc, but the rotation within each class and the final
  class are *sampled* from the top candidates with a power-transformed uniform
  draw (`rank = 1 - u^p`, `p = 2` during the annealed phase, `p = 4` after).
  Sampling replaces argmax so that near-ties keep exploring.
- **inpl**: class fixed, only rotation and shift are updated.
- **prob** and **prob_snhc**: the class, rotation, and shift are read from a
  globally computed probability table (see
  [sampling](sampling_and_fractional_updates.md)); the search step only
  commits and polishes them.

Shifts are refined with L-BFGS-B on the same objective, bounded by `trs`,
starting from a small coarse box around the seed. Shift search is disabled in
the first two iterations and `trs` is set automatically to `0.07 * mask
radius` (clamped to 5 to 6 pixels) once the search space fraction exceeds 75
percent, so shifts are only trusted after classes have stabilized.

The committed pose is discrete in angle. An optional joint continuous
`(sx, sy, theta)` polish, described in
[continuous in-plane refinement](continuous_inplane_refinement_abinitio2D.md),
can then refine it without changing what was selected.

## Algorithm: the restoration step

Each aligned particle plane is inserted into its class accumulator with a
Kaiser-Bessel kernel on a 2x padded Fourier grid, keeping even and odd halves
separate:

```text
B_k(q) = sum_{i in k} CTF_i(q) X_i^aligned(q) / sigma2_i(q)
D_k(q) = sum_{i in k} |CTF_i(q)|^2         / sigma2_i(q).
```

Without ML regularization the `1/sigma2` factors are absent. `B_k / D_k` is
the least-squares (Wiener-numerator over CTF-squared) estimate of the class
mean; each particle contributes to the class in proportion to how much signal
its CTF transfers at that frequency.

Restoration then:

1. divides `B_k` by `D_k` (bare division inside Nyquist);
2. soft-masks, subtracts the edge mean, and computes the FRC between the even
   and odd class averages;
3. with `ml_reg=yes`, converts the FRC to a signal-to-noise estimate
   `SSNR(k) = FRC/(1-FRC)`, sets `tau2(k) = SSNR(k) sigma2_k(k)` where
   `sigma2_k` is the mean noise power implied by `D_k`, and re-divides with

   ```text
   A_k = B_k / (D_k + 1/(tau tau2(k)))      (tau = 1 by default),
   ```

   which is the MAP estimate under a Gaussian prior with per-ring signal
   power `tau2`. Rings below index 6 are left unregularized because the
   signal there is effectively infinite. Halves with fewer than three
   particles receive an additional unit ridge;
4. replaces the lowest-resolution rings of both halves by the merged average
   up to the ring where the FRC drops below 0.7 (at least ring 4), because the
   halves are indistinguishable there and sharing them stabilizes matching;
5. applies the real-space gridding correction (the reciprocal of the
   Kaiser-Bessel instrument function) after the inverse FFT.

**FRC.** For class `k`,

```text
FRC_k(q) = Re sum_ring E_k O_k* / sqrt(sum |E_k|^2 sum |O_k|^2),
```

gives a per-class resolution at FRC 0.5 and 0.143. It is used to build a
per-class Wiener filter `2 FRC/(1+FRC)` for the next search step (unless ML
regularization has already done the equivalent in the restoration), to set the
global search band, and to rank classes.

## Fractional restoration

When only a subset of particles is updated in an iteration, the previous
accumulators are carried forward class by class. If `rho_k` is the fraction of
class `k`'s active members that were updated,

```text
B_k <- B_k^partial + (1 - rho_k) B_k^previous,
D_k <- D_k^partial + (1 - rho_k) D_k^previous.
```

A class with no updated members keeps its previous estimate intact, a fully
updated class replaces it, and the blend is exact in the sense that the total
sampling density remains that of the full dataset. Using the realized
per-class fraction rather than the requested global one matters because
balanced sampling gives different classes different fractions. Stage 1 of
`abinitio2D` runs without this memory so that random initial references are
overwritten rather than blended.

## Convergence

Per particle, the search records whether the class changed (`mi_class`), the
fraction of the class space visited (`frac`), and the in-plane angular change.
A run has converged when the averages satisfy

| regime | condition |
|---|---|
| full updates | `mi_class > 0.80` and `frac > 90 %` |
| fractional updates or streaming | `mi_class > 0.95` and `frac > 95 %` |
| in-plane only | mean in-plane change below 0.5 degrees |

Requiring a high visited fraction prevents an early stochastic iteration, which
by construction visits little, from looking converged.

## Rationale

- Whitening by `sigma2` and weighting by ring count turns the objective into a
  proper log-likelihood, so scores from different particles are on one scale.
  That is what allows probabilistic sampling and ML regularization to use the
  same numbers.
- Stochastic first-improvement acceptance is a form of simulated annealing:
  it keeps poor early classes from capturing everything by not always jumping to
  the current best, while the annealed search fraction guarantees eventual
  exhaustive search.
- Keeping even and odd sums separate costs nothing and gives an unbiased
  resolution estimate per class, which in turn drives the filter and the
  bandwidth used in the next search step.

## Implementation

- Objective and all-rotation evaluation: `src/main/pftc/simple_polarft_corr.f90`.
- Noise model: `src/main/simple_euclid_sigma2.f90`, `src/main/simple_sigma2_state.f90`.
- Search strategies: `src/main/strategies/search/simple_strategy2D_*.f90`;
  annealing constants in `src/utils/math/simple_decay_funs.f90`.
- Accumulation and restoration: `src/main/class/simple_classaverager*.f90`;
  FRCs in `src/main/class/simple_class_frcs.f90`.
- Convergence: `src/main/simple_convergence.f90`, limits in
  `src/defs/simple_defs_conv.f90`.
- Iteration driver: `src/main/strategies/parallelization/simple_cluster2D_strategy.f90`.
