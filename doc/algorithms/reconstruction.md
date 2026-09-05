# Cartesian 3D Reconstruction

## Problem

Given particle images with fixed poses, shifts, CTFs, state labels, half-set
labels, and per-shell noise powers, estimate the even and odd 3D volumes of
each state. This is the linear inverse problem of tomography with a
frequency-dependent, sign-changing transfer function. Poses are treated as
known; their estimation is the concern of [refine3D](refine3d.md).

## Forward model

For volume `x` with Fourier transform `F`, particle `i` observes

```text
y_i = C_i S_i G_i F + n_i,       n_i ~ N(0, diag(sigma2_i(k))),
```

where `G_i` gathers the central section at orientation `R_i` from the 3D
Fourier grid with the Kaiser-Bessel (KB) kernel, `S_i` is the shift phase
ramp, and `C_i` is the CTF. The KB kernel is

```text
apod(u) = (1/W) I_0( beta sqrt(1 - (2u/W)^2) ),   |u| <= W/2,
```

with window `W = 3` grid points and `beta` chosen for oversampling factor
`alpha = 2`; the 3D Fourier grid is `2x` the native box so that the kernel's
real-space envelope falls outside the reconstructed volume. Particles are
grouped by state and half before accumulation, and point-group symmetry is
imposed by inserting every symmetry-related orientation of each particle.

## Gridding backend (default)

The weighted least-squares solution of the forward model, ignoring the
coupling between neighboring grid points, is the ratio of two accumulators:

```text
B(q) = sum_i sum_{g in sym} G_{ig}^T [ C_i S_i^* y_i / sigma2_i ]
D(q) = sum_i sum_{g in sym} G_{ig}^T [ |C_i|^2 / sigma2_i ].
```

`B` is the CTF-weighted backprojection of the whitened data, `D` the
CTF-squared sampling density. Each particle plane is deposited into the `3^3`
grid neighborhood of every polar sample with separable KB weights normalized
to unit sum per axis. Without ML regularization the `1/sigma2` factors are
absent. Only the `k <= 0` half is stored; Friedel symmetry supplies the rest.

Restoration:

1. **Density division.** `F(q) = B(q) / D(q)` inside the resolution limit,
   zero beyond. The division is bare: with ML regularization the
   Wiener term has already been added to `D` ([refine3D](refine3d.md)), and
   without it the density of a real dataset is never small enough inside the
   limit to need a floor.
2. **Inverse FFT** to the padded real-space grid.
3. **Gridding correction.** Multiplication by the reciprocal of the KB
   kernel's real-space envelope, computed as the discrete Fourier transform of
   the normalized stencil so that it matches exactly what deposition did.
4. **Crop** to the native box.

The merged map is restored from `B_even + B_odd` and `D_even + D_odd` (with ML
regularization, after each half has received its own prior), not by averaging
restored halves; the FSC is computed on the halves before the gridding
correction, which is a pure amplitude envelope and cancels in the ratio.

## PCG backend (opt-in)

`rec_backend=pcg` solves the same model without the decoupling
approximation. With `K_i = sigma_i^{-1} C_i S_i G_i`, the normal equations are

```text
( sum_i K_i^T K_i + Lambda ) F = sum_i K_i^T sigma_i^{-1} y_i,
```

i.e. `(H + Lambda) F = B`, where `H` couples grid points through the KB
kernel overlap and `Lambda` is an optional prior. Because every `K_i^T K_i`
is a shift-invariant deposit of `|C_i|^2 / sigma2_i` at one orientation, `H`
is a convolution operator on the padded lattice whose kernel is the Fourier
transform of `D` deposited with the *squared* stencil. Applying `H` is then
one padded FFT, a pointwise product, and an inverse FFT, so the cost of an
iteration is independent of the number of particles: the particle loop runs
once, to build `B` and `D`, and never again.

The solve is left-preconditioned conjugate gradients with the reciprocal of
`D` plus a shell-relative floor as preconditioner (an absolute floor would
amplify unsampled modes by orders of magnitude). A soft spherical support
`P` may be imposed as `(P H P) u = P b`, which keeps the operator symmetric
positive semidefinite. Modes outside the reachable Fourier sphere stay zero.
Convergence is reported as the relative residual `||b - HF||/||b||` and the
relative update `||dF||/||F||`; the iteration count is kept short because the
Toeplitz kernel is an approximation of `H` whose positive-curvature region is
limited.

Two solves are made from one accumulation: a base solve producing
unfiltered halves, whose FSC remains the resolution authority, and an ML
replay that adds the FSC-derived `1/tau2` shell prior and warm-starts from the
previous same-half solution. PCG maps receive no gridding correction and no
second density division. The backend currently excludes trailing
reconstruction, projection-direction reconstruction, and conical
regularization; those requests are rejected rather than silently rerouted.

## Weighted and sparse variants

When particle weights lie in `[0, 1]` (flex PCA state maps), the density
becomes small and erratic in poorly occupied regions. A shell-relative floor
`D(q) >= mean_shell(D) / 1000` is applied before division in that path only.

## Complexity

Accumulation is `O(N_particles x N_sym x N_polar x 27)`; restoration is one
3D FFT on the padded grid. Distributed execution reduces per-partition
accumulators by summation, which is exact.

## Rationale

- Depositing whitened, CTF-multiplied data and dividing by CTF-squared
  density is the least-squares estimator for a frequency-dependent transfer
  function: where the CTF is near zero for one defocus, particles at other
  defoci fill the gap in `D`, and no particle is ever divided by its own CTF.
- Kaiser-Bessel interpolation has near-optimal energy concentration for a
  given window, and its instrument function is known analytically, so the
  interpolation blur is removed exactly rather than approximately.
- Even/odd halves reconstructed from disjoint particles are the only way to
  obtain a resolution estimate whose noise terms are independent.

## Implementation

- Gridding accumulator and restoration: `src/main/volume/simple_reconstructor.f90`;
  kernel in `src/main/interp/simple_kbinterpol.f90`, envelope in
  `src/main/interp/simple_gridding.f90`.
- Plane preparation (CTF, shift, whitening): `src/main/image/simple_image_ctf.f90`.
- Even/odd pair restoration: `src/main/commanders/simple/simple_commanders_rec_distr.f90`.
- Half-map diagnostics (FSC, cFAR): `src/main/volume/simple_halfmap_diagnostics.f90`.
- PCG operator and solver: `src/main/volume/simple_reconstructor_pcg.f90`;
  orchestration in `src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90`.
- Policies: `doc/policies/KB_Interpolation_Policy.md`,
  `doc/policies/reconstruct3D_pcg_policy.md`.
