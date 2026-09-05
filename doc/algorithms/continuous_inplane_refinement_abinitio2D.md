# Continuous In-Plane Refinement

## Problem

The 2D and 3D search steps choose the in-plane rotation from a grid of
`nrots` angles, `nrots = 2 pftsz` where `pftsz ~ pi * mask radius` is the
number of angular samples of the polar Fourier transform. The grid spacing is
therefore about one pixel of arc at the mask edge. Shifts are already
continuous. This chapter describes how the committed pose is polished to a
continuous angle jointly with the shift, and why the polish is kept strictly
separate from the search that selected the pose.

## Model

Fix the particle, the class or projection reference `A`, and the CTF. Let
`c(phi, k)` be the polar Fourier coefficients of the reference and `x(phi, k)`
those of the particle. The whitened Euclidean loss at continuous rotation
`theta` and shift `s` is

```text
L(s, theta) = 1 + [ sum_k w_k ( |A|^2 CTF^2 - 2 Re x_s^* R_theta A ) ] / (2 nrots sum_k w_k |x|^2),
w_k = k / sigma2(k),
```

where `x_s` is the particle with the shift phase ramp applied. Because the
reference and particle are sampled on a ring in angle, the dependence on
`theta` of the cross term is a circular correlation, which is a finite
trigonometric series:

```text
L(s, theta) = c_0(s) + sum_{m=1}^{pftsz-1} 2 Re[ c_m(s) e^{i m theta} ]   (+ Nyquist term),
dL/dtheta   = - sum_m 2 m Im[ c_m(s) e^{i m theta} ],
dL/ds_j     = b_0j(s) + sum_m 2 Re[ c^{(j)}_m(s) e^{i m theta} ].
```

The angular coefficients `c_m` are obtained by one FFT along the ring of the
products `x_s^* A`; the shift derivative coefficients `c^{(j)}_m` come from
the same FFT applied to `x_s^* A` multiplied by the shift-phase argument in
`x` and `y`. The three series are computed in one batched FFT and evaluated
at any real `theta` directly, so the loss and its full gradient at a
fractional angle cost the same as at a grid angle and are exactly smooth in
all three variables. The argument products must be formed before the angular
transform; forming them after would require a convolution in angular
frequency.

For the correlation objective, `cc = N(theta)/sqrt(D(theta) C)` with `N` the
cross series, `D` the shift-independent reference power series, and `C` the
particle power; its derivatives follow by the quotient rule, with the angle
derivative `(N' - N D'/(2D))/sqrt(DC)`.

## Algorithm

With `inpl_cont=yes`, the search runs unchanged until a class (or state and
direction) and an integer rotation index have been committed. Then one
bounded joint solve polishes the pose:

1. seed at the committed integer rotation and its native shift;
2. bound the rotation variable to plus or minus two grid cells around the
   seed and the shifts to `[-trs, trs]`; other cells are never rescanned;
3. run L-BFGS-B on `L(s, theta)` with the analytic gradient (three
   variables, at most `maxits_sh` iterations);
4. accept the result only if it is finite, valid for the objective's range,
   improves on the seed's re-scored value by a relative tolerance of `1e-4`
   (with a 1 percent floor for the correlation objective), and does not sit
   within `1e-3` of the range of any bound. A solution pinned to a bound is a
   sign that the true optimum lies in a neighboring cell the search rejected,
   and is not trusted.

On acceptance the fractional angle is persisted (`e3`), the nearest integer
index is kept for compatibility, and shift and score are taken from the same
joint pose. On rejection the committed discrete pose stands with its
re-scored value.

## What the polish must not change

The joint solve runs exactly once per particle per iteration, after
selection. Candidate scoring, probability-table construction, and shift
seeding all use the discrete route. Two consequences are deliberate:

- **Search trajectory identity.** The sequence of classes and cells visited,
  and the assignment committed, is identical with and without the polish.
  Polishing candidates *before* selection would sharpen the score contrast
  the stochastic sampler sees and collapse exploration on difficult
  particles.
- **Statistics parity.** Every statistic that steers the schedule (in-plane
  angular change, class overlap, visited fraction) is computed from the
  discrete cells on both sides of the comparison. A fractional angle would
  otherwise register as sub-grid motion, read as premature convergence, and
  truncate the annealing that difficult particles need. The fractional angle
  is used only by reconstruction and persisted for the next iteration's seed.

The polish is available for the raw Euclidean and correlation objectives. It
is not applied in the SGD streaming variant or in time-series shift-only
search.

## Rationale

The alignment precision that reconstruction can use is limited by the noise,
not by the grid, once the grid spacing is finer than the angular uncertainty.
The polish therefore mainly removes the systematic quantization of the
angle at the highest frequencies, where one pixel of arc at the mask edge is
a full phase cycle. Separating polish from search keeps the search's
annealing schedule, which was tuned on the discrete route, valid by
construction.

## Implementation

- Joint objective and gradients: `src/main/pftc/simple_polarft_corr.f90`.
- Joint optimizer and acceptance guards: `src/main/pftc/simple_pftc_shsrch_grad.f90`.
- 2D and 3D integration: `src/main/strategies/search/simple_strategy2D_srch.f90`,
  `src/main/strategies/search/simple_strategy3D_srch.f90`.
- Validation harness: `production/tests/simple_test_continuous_inplane_rotation2D*.f90`.
- Design note: `doc/implementation_notes/continuous_inplane_rotation_polar.md`.
