# Cartesian 3D Reconstruction

## Problem

Estimate dense even and odd 3D volumes from fixed particle poses, shifts,
state labels, CTFs, half assignments, and optional noise weights. Pose search
is upstream: reconstruction treats all of those quantities as immutable input.

## Shared forward model

Let `x` be a real 3D volume and `F` its Fourier transform. For particle `i`,
`G_i` gathers the oriented Fourier plane with the Kaiser-Bessel (KB)
interpolator, `S_i` applies the shift phase, `C_i` is the CTF, and `N_i` is the
diagonal noise covariance. The observation model is

```text
y_i = C_i S_i G_i F x + noise.
```

Particles are separated by state and even/odd half before accumulation.
Point-group symmetry is applied by inserting or gathering every symmetry-related
orientation under the same data model.

## Gridding backend

The default backend backprojects the CTF-weighted particle planes to a 3D
Fourier numerator and accumulates a CTF-squared sampling density:

```text
B(q) = sum_i G_i^dagger conj(C_i S_i) N_i^-1 y_i
D(q) = sum_i G_i^dagger |C_i|^2 N_i^-1.
```

Restoration folds Fourier symmetry, regularizes/divides by the density, applies
the discrete KB gridding correction, and inverse-transforms the dense map.
The output uses the data-quotient amplitude convention; no legacy box-size
multiplier or divider remains.

## PCG backend

`rec_backend=pcg` is an opt-in fixed-pose weighted least-squares solve. With

```text
K_i = N_i^-1/2 C_i S_i G_i F,
```

it solves

```text
(sum_i K_i^dagger K_i + Lambda) x =
 sum_i K_i^dagger N_i^-1/2 y_i.
```

The normal operator contains only `|C_i|^2/sigma2_i`; the unit-modulus shift
phase cancels there but remains complex in the right-hand side. A soft
spherical support `P` changes the system to `(P H P)u=P b`, keeping the
constrained operator symmetric positive semidefinite.

## Fused accumulation and operator construction

Particles are read in bounded batches. One KB traversal per particle updates
the complex right-hand side `B` and real density precursor `D`; PCG iterations
then perform no image I/O.

The production kernel operator converts `D` into a padded Toeplitz/Gram kernel.
One padded FFT, pointwise multiplication, and inverse FFT apply the approximate
normal operator, so iteration cost is independent of particle count. The exact
matrix-free particle loop is a small-fixture oracle, not a real-data fallback.

Both the gather and deposition KB envelopes use exact discrete stencil
transforms and are removed on the correct side of the operator. The unknown is
native-size, while Fourier operations use the same 2x padded lattice needed for
non-wrapping Gram convolution. Scatter coloring handles interior windows in
parallel and wrapped rim windows serially to retain deterministic summation.

## Preconditioner and solve

The reciprocal preconditioner is the sampling density plus a shell-relative
floor. An absolute floor would over-amplify poorly sampled modes by many orders
of magnitude. Modes outside the reachable Fourier sphere remain zero.

Left-preconditioned conjugate gradients report the true residual
`||b-Hx||/||b||` and the relative update `||dx||/||x||`. Non-positive curvature,
non-finite state, invalid preconditioners, and zero RHS are hard failures.
Production kernel solves are deliberately short because over-iterating the
approximate operator can leave its useful positive-curvature region.

## Base and ML maps

Refinement PCG performs two solves from one raw accumulation. The base solve
produces `_unfil` halves; their FSC/cFAR remains the resolution authority. An
ML replay adds FSC/SSNR shell precision and produces the standard halves,
warm-starting from the previous same-half ML map when compatible, otherwise
from the base solution. Priors are never folded into raw `B` or `D`.

After either backend restores dense halves, volume assembly owns the common FSC,
merged map, nonuniform filtering, mask, project, and postprocessing path. PCG
maps must not receive gridding correction or a second density division.

## Current PCG exclusions

Production PCG rejects fractional/trailing reconstruction, projection-direction
reconstruction, conical-FSC regularization, and matrix-free workflow execution.
It must not silently fall back to another estimator.

## Implementation

- Gridding: `src/main/volume/simple_reconstructor.f90`; even/odd pair I/O and
  restoration are contained in `restore_state_from_parts` in
  `src/main/commanders/simple/simple_commanders_rec_distr.f90`.
- Shared half-map diagnostics (FSC, cFAR, resolution, report):
  `src/main/volume/simple_halfmap_diagnostics.f90`.
- PCG operator and solver: `src/main/volume/simple_reconstructor_pcg.f90`.
- Backend orchestration: `src/main/strategies/parallelization/simple_rec3D_strategy.f90`
  and `src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90`.
- Policy: `doc/policies/reconstruct3D_pcg_policy.md`.

