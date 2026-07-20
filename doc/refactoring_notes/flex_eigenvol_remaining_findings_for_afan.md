# `flex_eigenvol`: Remaining Review Findings for Afan

**Date:** July 13, 2026  
**Scope:** Follow-up to the theoretical and implementation review against
`doc/for_developers/ideas/fourier_domain_latent_eigenvolumes.md`.

## Corrections Included in the Current Patch

Two blocking inconsistencies from the review have been corrected.

### Observation and Operator Representation

`gen_fplane4rec` now has an opt-in `observation_model` path. Normal
reconstruction callers retain their existing representation:

```text
cmplx_plane = C* y / sigma^2
ctfsq_plane = |C|^2 / sigma^2
```

For `flex_eigenvol`, the plane contract is now:

```text
cmplx_plane    = shifted y / sigma
transfer_plane = C / sigma
ctfsq_plane    = |C / sigma|^2
```

Here `sigma = sqrt(sigma^2)`. The particle origin shift is applied to the
observation once and is not included again in the forward model. Projection
uses `transfer_plane * P(B)`, and latent M-step insertion explicitly applies
`conjg(transfer_plane)` to the residual. The E-step and M-step therefore use a
single weighted least-squares objective in the same observation space.

### Latent Coordinate Consistency

The post-E-step calls that centered, rotated, and rescaled only the posterior
means `z` have been removed. Those transforms did not update the basis,
posterior covariance, or latent prior variances. Latents, posterior moments,
mode variances, basis volumes, coordinate tables, and trajectory maps now stay
in one coordinate system throughout the fit and at output.

## Remaining Findings

### 1. Canonicalize the Basis Before Calling the Outputs Eigenvolumes

The current outputs are valid projected factor loadings, but they are not yet a
canonical eigenbasis. The basis is not orthogonalized under a sampling-weighted
Fourier inner product, rotated to diagonal latent covariance, ordered by
explained signal variance, or sign-aligned between runs.

A final canonicalization should apply one invertible transform consistently to
the basis and particle coordinates so that `B*z` is unchanged. Posterior
covariances must receive the same coordinate transform if canonicalization is
performed before the final output stage. Closely spaced modes will still need
subspace rather than component-wise comparison.

Until this exists, logs and documentation should describe the components as
linear modes or factor loadings rather than implying uniquely ordered
covariance eigenvectors.

**Fix implemented:**

The final E-step now accumulates the observation-induced basis Gram matrix
before latent-prior precision is added. After that E-step, a terminal
canonicalization orthonormalizes the basis in this sampling-, CTF-, and
noise-weighted Fourier metric, diagonalizes the fitted signal covariance, sorts
components by decreasing signal eigenvalue, and applies a deterministic sign
convention.

The same invertible transform is applied to the Fourier basis, particle latent
coordinates, posterior covariance matrices, and latent prior variances. The
product `B*z`, and therefore the fitted particle model and residuals, is
unchanged. Trajectory volumes use the square root of the canonical signal
eigenvalue as their one-sigma displacement, and the eigenvalues and explained
fractions are written alongside the coordinate table.

Canonicalization is deliberately the last model operation. There is no M-step
after it. Numerically rank-deficient basis metrics fail with a recommendation
to reduce `neigs` rather than silently regularizing an unidentifiable mode.
Closely spaced modes still require subspace rather than component-wise
comparison.

### 2. Enable ML Weighting and Use the Average-Map Nonuniform Filter

The plane machinery now supports the correct `1/sigma` whitening contract, but
the public and private `flex_eigenvol` commanders force `ml_reg=no`. Current
runs therefore use the image normalization performed by preprocessing and an
implicit unit-noise likelihood.

Recommended work:

- allow `ml_reg=yes` and validate the sigma spectra used by the projected model;
- report the weighted objective per observation rather than raw plane energy;
- reuse the residual shell-power estimate produced by refinement for whitening
  and residual diagnostics rather than estimating a separate spectrum inside
  `flex_eigenvol`;
- add an objective/convergence criterion instead of relying only on a fixed
  iteration count;
- derive one nonuniform filter from the supplied average map and apply that
  fixed filter consistently to every fitted mode and trajectory volume;
- use the average-map filter support, weighted reprojection residuals, and
  synthetic recovery tests as the filtering and validation policy.

The mode fit should not estimate a separate data-dependent filter for each
component. Sharing the average-map nonuniform filter keeps the fitted modes in
the spatial-frequency support already justified by the consensus reconstruction
and avoids giving weak modes artificially detailed appearance.

### 3. Document and Test the Approximate Normal Operator

The E-step uses explicit projected basis planes, but the M-step uses the
existing gridding-style voxel-local density approximation. Interpolation
cross-voxel terms from the exact `P^* P` operator are not retained, and masking
and low-pass filtering are applied after the solve. This is a reasonable
generalized-EM/reconstruction approximation, but not an exact PPCA M-step.

Required numerical tests:

- an adjoint test for projection and residual backprojection;
- a homogeneous simulation with nontrivial CTFs and shifts that recovers zero
  latent signal;
- a one-mode simulation that recovers coordinates and the mode up to scale and
  sign;
- shared-memory versus distributed sufficient-statistics equivalence;
- objective tracking that distinguishes expected regularization changes from
  numerical regressions.

The current module-local projected-model test covers only M-step statistics
file I/O and is not registered as a standalone CTest.

### 4. Harden Distributed E-Step Reduction and Cleanup

The distributed E-step reducer currently trusts worker row indices and does not
check for out-of-range, duplicate, missing, or particle-index-mismatched rows.
Mode moments are divided by total `nptcls` when a worker record may be absent or
invalid. E-step state, temporary basis maps, and E-step part files also remain
after successful reduction.

Recommended changes:

- validate the expected `ncomp`, row interval, particle index, iteration, box,
  sampling, and shell range in every part header;
- maintain a `seen(row)` array and fail on duplicate or missing rows;
- divide aggregate moments by the validated record count;
- delete state, part, and temporary basis artifacts after successful reduction,
  while preserving them on failure for diagnosis.

### 5. Improve Initialization and Output Diagnostics

Latents are initialized with deterministic sinusoids over particle row order.
Particle order can correlate with acquisition time, micrograph, or preprocessing
order, so this can seed nuisance modes in a non-convex fit.

Useful alternatives are low-resolution random modes with a fixed seed,
average-map residual modes, or known state-difference modes. Output tables
should also include state, score, and preferably a normalized weighted residual
diagnostic.

## Suggested Order of Work

1. Add synthetic observation/operator and adjoint tests.
2. Enable and validate shell-dependent ML noise whitening.
3. Apply the nonuniform filter derived from the average map.
4. Add final basis canonicalization with a model-preserving coordinate transform.
5. Harden distributed metadata/reduction and clean temporary artifacts.
6. Improve initialization and user-facing diagnostics.

The architecture is otherwise aligned with the idea note: it uses a common 3D
Fourier basis, real particle latents, pose-aware projection, joint latent normal
equations, coupled component updates, physical Fourier-volume constraints, and
a standalone workflow independent of `refine3D`.
