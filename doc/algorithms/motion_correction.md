# Motion Correction

## Problem

A dose-fractionated movie records `n` frames of the same field of view while
the specimen drifts and deforms under the beam. The task is to estimate a
displacement field `u(x, y, t)` and produce an integrated micrograph in which
every frame has been warped by `-u` before summation. The estimate is made in
two stages: a single translation per frame for the whole field of view, then a
smooth space-time polynomial for the residual local motion. A model is only
accepted if it reproduces the measured data within tolerance; otherwise the
simpler model is used.

## Preprocessing

Frames are gain-corrected, optionally EER-fractionated, hot-pixel cured
(values outside 6 standard deviations replaced by a local window average),
optionally truncated by dose, and Fourier-cropped to the smallest FFT-friendly
box that still supports the final alignment resolution `lp_stop`. The central
frame is the fixed reference (frame 1 under the RELION convention).

## Global alignment

**Objective.** With frame transforms `F_j` and per-frame weights `w_j`, the
shift vector `d = (d_1..d_n)` maximizes the leave-one-out correlation of each
frame with the weighted sum of the others,

```text
J(d) = sum_i corr( F_i shifted by d_i ,  sum_{j != i} w_j F_j shifted by d_j ),
```

evaluated in a band `[hp, lp]` with a cosine-tapered band-pass and an
optional B-factor envelope `exp(-B s^2 / 4)` (`B = 50 A^2` by default in
preprocessing). Leaving frame `i` out of its own reference removes the
trivial self-correlation that would otherwise bias `d_i` to zero.

**Discrete phase.** Iterate (3 to 10 times): form the weighted reference,
compute the full cross-correlation map of each frame by FFT, take the integer
peak within `[-trs, trs]^2`, and refine to sub-pixel precision by parabolic
interpolation along each axis, `delta = (a - c) / (2(a + c - 2b))`. The band
tightens from `lp_start = 8 A` to `lp_stop = 5 A` over the iterations, and
the weights are recomputed from the correlations:

```text
cc_i <- max(cc_i, 0);  normalize to [0, 1] if max/min >= 1.5;
w_i = exp(-(1 - cc_i)) / sum_j exp(-(1 - cc_j)).
```

Convergence is an inter-iteration shift RMSD below 0.5 pixels after the band
has reached `lp_stop`. When a movie has fewer than 12 frames, shifts are
additionally constrained to a cubic polynomial in frame index.

**Continuous phase.** The same objective is then optimized per frame in
Fourier space by L-BFGS-B, where a shift is a phase ramp and the gradient is
analytic, within `+/- 5` pixels, for 2 to 5 iterations until the correlation
gain falls below 0.1 percent and the RMSD below 0.1 pixels. The global shifts
are applied to all frames before local estimation.

## Local patch alignment

The field of view is tiled into patches of about 200 A on a side (at least
200 pixels after scaling), `n_x = floor(width / 200 A)`, and the same hybrid
aligner is run on each patch stack independently, yielding a measured
trajectory `d_p(t)` per patch center `(x_p, y_p)`.

## Deformation model

The residual displacement is modeled separately in `x` and `y` as a
polynomial that is cubic in time and quadratic in space, with 18 terms:

```text
u_d(x, y, t) = sum_{k=1}^{18} c_{d,k} phi_k(x, y, t),
phi in { t, t^2, t^3 } (x) { 1, x, x^2, y, y^2, xy },
```

in normalized patch coordinates with `t` measured from the fixed frame, so
that `u = 0` there by construction. The coefficients are the least-squares
solution (by SVD) against the measured patch trajectories. Quadratic in space
captures the dome-like doming of the support film under the beam; cubic in
time captures the fast initial burst and slow later drift. This is the
MotionCor2 model (Zheng et al., Nature Methods 14, 331 (2017)).

**Robust variant.** `patch_refine` first fits on a padded grid, trims the 10
percent of points with the largest residuals, refits, evaluates the fitted
field at the patch centers to seed a second patch alignment, and fits the
final model to that.

## Model acceptance

The fit is judged by the root-mean-square deviation between fitted and
measured patch shifts, per axis, over all frames and patches. Both axes must
be below 4 pixels (5 under the RELION convention). If not, the patch grid is
halved in each direction and the fit retried once; if it still fails, the
local model is discarded and the global translation alone is applied. A
polynomial that cannot reproduce the patch trajectories is more likely to be
fitting patch-alignment failures than real motion.

## Dose weighting and integration

When dose weighting is on, each frame is filtered before warping by the
Grant and Grigorieff critical-exposure model. With cumulative exposure `e_i`
at frame `i` and critical exposure `N_e(k) = 0.245 k^{-1.665} + 2.81`
electrons per square Angstrom at spatial frequency `k` (scaled by 0.8 at
200 kV and 0.64 at 100 kV), the per-frame Fourier weight is
`q_i(k) = q_1 r^{i-1}` with `r = exp(-dose_per_frame / (2 N_e))`, normalized so
that `sum_i q_i^2 = 1` at every frequency. The frames are then warped through
the accepted deformation field, weighted by `w_i`, and summed. A second
uniform-weight, non-dose-weighted sum is produced for CTF estimation, where
the high-frequency Thon rings that dose weighting attenuates are needed.

## Implementation

- Workflow and integration: `src/main/motion/simple_motion_correct.f90`,
  `src/main/motion/simple_motion_correct_iter.f90`.
- Hybrid discrete/continuous aligner: `src/main/motion/simple_motion_align_hybrid.f90`.
- Patch trajectories and deformation model: `src/main/motion/simple_motion_patched.f90`.
- Frame weights: `src/utils/math/simple_stat.f90`; dose weighting:
  `src/main/image/simple_image_ops.f90`.
