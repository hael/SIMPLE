# `flex_eigenvol` Policy: Fourier-Domain Latent Eigenvolumes

This document records the implemented scientific and software contract for
SIMPLE's projection-aware linear 3D variability workflow. Forward-looking
developments are identified separately and are not part of the production
contract.

## 1. Scope

This policy applies to:

- the `flex_eigenvol` command
- its shared- and distributed-memory execution paths
- the projected latent-model kernels used by those paths
- in-process consumers of `projected_latent_fit_result`, including latent
  trajectory chunking

The current workflow is standard iterative PCA expressed as alternating E- and
M-steps. It is not probabilistic PCA, despite the historical `ppca_inmem` name
that motivated the simplified implementation.

The words **must**, **should**, and **may** distinguish required behavior,
recommended practice, and optional extensions.

## 2. Production Policy

The implemented baseline is:

- a supplied, fixed mean volume
- fixed `ptcl3D` poses and shifts
- real particle latent coordinates
- complex Fourier eigenvolumes
- projection-aware, CTF-aware alternating least-squares updates
- Hermitian Fourier inner products and Hermitian normal-equation solves
- no external sigma array, fitted noise variance, posterior covariance, or
  per-mode uncertainty term
- at most 20 fitted eigenvolumes
- final metric orthonormalization, descending empirical-variance ordering, and
  deterministic signs
- ICM-based rank selection for downstream latent-feature denoising

The shared- and distributed-memory paths must implement the same model and
produce the same kinds of final outputs. Distribution may change how
sufficient statistics are accumulated and transported, but not the scientific
objective.

## 3. Implementation Status

| Capability | Status | Policy meaning |
| --- | --- | --- |
| Projection-aware Fourier EM-PCA | Implemented | Production baseline |
| Hermitian E-step formulation | Implemented | Must be retained |
| Coupled eigenvolume M-step | Implemented | Production baseline |
| Shared-memory execution | Implemented | `nparts=1` path |
| Distributed-memory execution | Implemented | `nparts>1` worker/reduction path |
| Observation-metric orthonormalization | Implemented | Required final canonicalization |
| Descending empirical eigenvalue order | Implemented | Required output order |
| Deterministic component signs | Implemented | Required for reproducible outputs |
| Maximum of 20 fitted modes | Implemented | Current safety and cost limit |
| ICM rank selection in latent trajectory chunking | Implemented | Current downstream denoising policy |
| External or model-estimated sigma | Not implemented by policy | Must not be reintroduced incidentally |
| Posterior latent covariance | Not implemented by policy | Research item requiring a likelihood |
| Spatial structural-variance map | Not implemented | Forward-facing diagnostic output |
| Posterior/predictive uncertainty map | Not implemented | Forward-facing probabilistic extension |
| Project-FSC shrinkage and mode-specific FSC | Not implemented | High-priority cross-validated regularization development |
| Joint pose or shift refinement | Not implemented | Outside the current model |
| ICA or nonlinear latent dynamics | Not implemented | Optional post-PCA research direction |

## 4. Scientific Model

For particle `i`, the volume model is

```text
V_i(k) = mu(k) + sum_q z_iq B_q(k)
```

where:

- `mu` is the supplied consensus volume
- `B_q` is eigenvolume `q`
- `z_iq` is the real coordinate of particle `i` along mode `q`

The particle observation model is

```text
y_i = P_i mu + sum_q z_iq P_i B_q + residual_i
```

`P_i` is particle-specific. It incorporates the fixed 3D pose, origin shift,
CTF, Fourier sampling, and the existing reconstruction conventions. Particle
Fourier planes are prepared through `gen_fplane4rec` with the observation-model
transfer stored, but without passing an external sigma array.

The mean is fixed throughout the fit. `flex_eigenvol` estimates the basis and
particle coordinates; it does not re-estimate the consensus map.

### 4.1 Real latents and complex Fourier volumes

The eigenvolumes have complex Fourier coefficients because they represent real
3D volumes in Fourier space. The physical mixing coefficients are nevertheless
real. The E-step therefore constructs complex Hermitian products and takes
their real symmetric representation for the real latent solve.

Complex Fourier arithmetic must not be mistaken for complex-valued particle
coordinates.

### 4.2 Meaning of EM

The E/M terminology describes alternating conditional least-squares updates:

1. infer point estimates of `z` with the basis fixed
2. update the basis with `z` fixed

There is no posterior expectation, posterior covariance, latent-prior
precision, or noise-variance update in the production objective. Gaussian
sampling is used only to initialize the latent matrix. The initial columns are
centered, orthogonalized, and scaled before the first M-step; that
initialization does not impose a Gaussian prior during fitting.

## 5. E-Step Contract

For each particle, first subtract the projected mean:

```text
r_i = y_i - P_i mu
```

Then construct

```text
G_i(q,r) = Re <P_i B_q, P_i B_r>_H
h_i(q)   = Re <P_i B_q, r_i>_H
```

where `<.,.>_H` is the Hermitian Fourier-plane inner product. The particle
coordinates solve

```text
G_i z_i = h_i
```

using `hermitian_solve`. If the solve reports failure, the current defensive
behavior is to set that latent row to zero.

The E-step also records:

- residual energy after subtracting the mean only
- residual energy after subtracting the fitted modes
- the average observed basis Gram matrix used by final canonicalization

The Hermitian construction is a required invariant. It must not be replaced by
a generic non-Hermitian or SVD formulation without a separately reviewed
scientific reason and equivalent model-preservation tests.

## 6. M-Step Contract

With `z` fixed, the M-step accumulates projection/backprojection sufficient
statistics for all modes and all mode pairs. It solves the coupled small system
at each expanded Fourier-grid location and then finalizes each reconstructed
basis volume using SIMPLE's reconstruction machinery.

The implementation includes a small scale-dependent ridge in the coupled
Fourier solve for numerical stability. This ridge is not a noise estimate, a
latent prior, or posterior regularization.

The basis may be regularized by the configured mask and low-pass limit. No
half-set FSC shrinkage is currently applied.

## 7. Iteration and Final Refit

`flex_eigenvol` performs the configured number of alternating M- and E-steps.
It then performs one final M-step and one final E-step before canonicalization.
The current workflow uses a fixed iteration count; it does not stop early from
an objective or convergence tolerance.

The fitted mode count is clamped to the range `1..20`. Asking for more than 20
modes produces a warning and fits 20. A rank-deficient final observation metric
is a hard failure and should normally be addressed by reducing `neigs`.

## 8. Canonicalization and Eigenvalue Semantics

Alternating matrix factorization is invariant under an invertible change of
latent coordinates:

```text
B z = (B T) (T^-1 z)
```

Raw EM-PCA components therefore have no unique scale, order, rotation, or
sign. Every completed fit must be canonicalized before output.

The implementation uses:

- `M`, the average particle-observed basis Gram matrix from the final E-step
- `S_z`, the empirical covariance of the final latent coordinates

It finds a transform `T` such that the transformed basis and coordinates obey

```text
T^T M T                       = I
Cov(T^-1 z)                   = diag(lambda_1, ..., lambda_Q)
lambda_1 >= ... >= lambda_Q   >= 0
```

The transform is applied to the Fourier eigenvolumes and its inverse to every
particle coordinate, preserving each fitted model volume. Matrix invariants
are checked at runtime. Component signs are fixed deterministically using the
largest-magnitude particle coordinate, with project particle index as the tie
breaker.

The reported eigenvalues are empirical latent variances in a basis that is
orthonormal under the average observation-induced metric. They are not noise
variances or posterior variances. Explained fractions are relative only to the
fitted modes, not to total particle-image variance.

For nondegenerate modes, trajectory amplitudes use `sqrt(lambda_q)`, so the
standard trajectory samples represent `mu + a sqrt(lambda_q) B_q` for `a` in
`-3..3`. The current writer uses a unit step for a numerically zero eigenvalue;
such a mode should not be interpreted as supported variability.

## 9. Shared- and Distributed-Memory Contract

### 9.1 Shared memory

When `nparts=1`, the strategy accumulates M-step statistics and performs
E-step inference within the main process, with OpenMP parallelism inside the
numerical kernels.

### 9.2 Distributed memory

When `nparts>1`, the master dispatches dedicated M-step and E-step workers:

- M-step workers write reducible basis right-hand sides and coupled
  cross-density statistics
- E-step workers write particle rows, point latent estimates, residual
  energies, and basis-metric contributions
- the master reduces these parts, solves and finalizes the common basis, and
  performs the same final canonicalization as the shared path

The distributed state passed between iterations contains particle indices and
point latent coordinates. It does not contain sigma values, posterior
covariances, or mode-variance priors.

Internal state and part-file versions are implementation details and must be
bumped whenever their binary payload changes. A change to one execution path
must be mirrored in the other or explicitly rejected as unsupported.

## 10. Inputs and Supported Operating Envelope

The implemented workflow requires:

- `oritype=ptcl3D`
- a supplied `vol1` mean map, normally from `abinitio3D` or `refine3D`
- exactly one `ptcl3D` state
- active particles with fixed poses and shifts
- a requested mode count between 1 and 20 after clamping

The particle source may follow the existing raw or alternate particle-source
path. Cropping, sampling, CTF, shift, mask, and low-pass conventions must be
inherited from the established builder, project, image, and reconstruction
objects rather than reimplemented inside the PCA kernel.

## 11. Outputs

For an output prefix derived from `outvol`, the command writes:

- eigenvolume 1 as `outvol`
- later eigenvolumes as `<prefix>_002.mrc`, `<prefix>_003.mrc`, and so on
- `<prefix>_zcoords.txt`, containing project particle index, canonical latent
  coordinates, and fitted residual energy
- `<prefix>_eigenvalues.txt`, containing eigenvalue, within-model explained
  fraction, and cumulative fraction
- `<prefix>_modeNNN_m3.mrc` through `m1`, `z0`, and `p1` through `p3`,
  containing mean-plus-mode trajectory maps
- corresponding `_diff.mrc` maps containing trajectory-minus-mean differences

The direct trajectory maps are interpretation aids. They do not replace
particle-subset or kernel-weighted reconstructions when a final high-resolution
state map is required.

In-process callers may receive `projected_latent_fit_result`, which owns the
particle indices, canonical `z`, mean-only and fitted residual energies, and
canonical eigenvalues. It does not own posterior uncertainty fields.

## 12. ICM Feature-Selection Policy

Fitting and downstream feature selection are separate operations. The PCA fit
estimates the requested number of modes, up to 20. It does not use ICM to alter
the fit or suppress eigenvolume files.

Downstream consumers that need denoised latent features should use
`select_spectral_rank_icm` on the descending canonical eigenvalues. The
implemented trajectory chunker:

1. selects a retained prefix with ICM, with at least one mode retained
2. standardizes each selected latent coordinate across particles
3. weights selected modes by normalized canonical eigenvalue
4. assigns zero weight to modes beyond the selected rank

Posterior uncertainty weighting is not part of this policy. New latent-space
consumers should reuse the common ICM selector rather than inventing independent
eigenvalue-gap heuristics.

## 13. Software Ownership

The implementation follows SIMPLE's layer boundaries:

- `src/main/pca/simple_projected_latent_model.f90` owns projected EM-PCA
  kernels, sufficient-statistics I/O, Hermitian inference, and canonicalization
- `src/main/pca/simple_projected_latent_result.f90` owns the reusable fit-result
  value object
- `src/main/strategies/parallelization/simple_flex_eigenvol_strategy.f90` owns
  shared/distributed orchestration, worker dispatch, and public outputs
- `src/main/pca/simple_diff_map_denoise.f90` owns the reusable ICM spectral-rank
  selector
- `src/main/nano/simple_trajectory_chunker.f90` consumes the fit and applies
  ICM-selected latent features to time-constrained chunking
- commanders and UI modules own workflow sequencing, defaults, and user-facing
  parameter registration

The numerical model must not acquire hidden builder, parameter, sigma, or PFTC
globals. Execution context must remain caller-owned and explicitly threaded
through the established interfaces.

## 14. Validation and Change Control

The dedicated projected-latent test currently covers:

- binary round-trip and reduction of M-step sufficient statistics
- canonical transform inversion
- observation-metric orthonormality
- diagonal empirical covariance
- descending eigenvalue order
- E-step part-file payload round-trip used by distributed reduction

Runtime canonicalization repeats the key matrix-invariant checks on real fits.

Changes should also be checked with a representative project in both
`nparts=1` and `nparts>1` modes. The comparison should include residual
reduction, eigenvalue order, subspace agreement, and reconstructed trajectory
maps. Exact component signs should be stable after deterministic anchoring;
small floating-point differences from reduction order are expected.

The following changes require an explicit policy update rather than a quiet
implementation edit:

- replacing the Hermitian formulation
- adding external or fitted sigma terms
- adding a latent prior or posterior covariance to the EM updates
- changing the canonical metric or eigenvalue meaning
- using uncertainty estimates for rank selection or denoising
- increasing the 20-mode production cap

## 15. Known Limitations

- The model is linear and uses fixed poses, shifts, CTF parameters, and mean.
- Pose, shift, CTF, or mean-map error can appear as apparent eigenvolumes.
- Preferred orientation can make the observed basis metric anisotropic or
  rank deficient.
- The mask and low-pass limit are the current basis regularizers; the existing
  project FSC is not yet consumed for eigenvolume shrinkage, and there is no
  mode-specific eigenvolume FSC.
- A fixed iteration count provides no formal convergence certificate.
- ICM selects features for downstream consumers but does not prove that a mode
  represents biological variability.
- Direct mean-plus-mode trajectories may overstate high-resolution detail.
- The workflow reports scalar canonical eigenvalues, not spatial variance or
  uncertainty maps.

## 16. Forward-Facing Developments

Everything in this section is outside the current production contract. These
items are worth considering, but must not be described as implemented until
code and validation exist.

### 16.1 Spatial structural-variance maps

The current canonical model already supplies the ingredients for a population
model-variance map. Let `b(x)` be the vector of real-space eigenvolume values at
voxel `x` and `S_z` the empirical latent covariance. Then

```text
v_structural(x) = b(x)^T S_z b(x)
```

In canonical coordinates, `S_z` is diagonal, giving

```text
v_structural(x) = sum_q lambda_q b_q(x)^2
```

This map would describe spatial variability represented by the fitted linear
model. It would not measure posterior uncertainty and would not by itself prove
that the variability is biological. A useful implementation should optionally
report per-mode contributions as well as their total.

### 16.2 Posterior and predictive uncertainty maps

A probabilistic extension could estimate a particle-specific latent posterior
covariance `C_i = Cov(z_i | y_i)` under a validated likelihood and propagate it
as

```text
v_posterior_i(x) = b(x)^T C_i b(x)
```

This quantity describes uncertainty about particle `i` under the assumed
model. It is not population heterogeneity. Uncertainty in the mean,
eigenvolumes, poses, CTF, and likelihood parameters is a separate contribution
and must either be reported separately or included in a clearly named total
predictive-variance map.

Any posterior extension must retain the Hermitian Fourier observation model.
Canonical rotations must be applied consistently to every covariance so that
propagated maps remain invariant. A pointwise Fourier variance must not simply
be inverse-transformed; real-space variance must be formed by covariance
propagation through the Hermitian-symmetric real-space basis.

Before posterior or predictive maps influence feature selection, trajectory
weighting, or denoising, require:

- a calibrated likelihood/noise model rather than an unverified sigma array
- agreement between shared- and distributed-memory implementations
- homogeneous simulations that should yield no structural signal
- simulations with known heterogeneous signal and known uncertainty
- half-set calibration and reproducibility checks by spatial-frequency shell
- evidence that structural variance, posterior uncertainty, pose error, and
  CTF/model mismatch are not being conflated

Until then, such maps should be diagnostic research outputs only. They must not
replace empirical eigenvalue ordering or ICM feature selection.

### 16.3 Cross-validated FSC shrinkage and mode validation

Any FSC available through the project is already derived from independent
even/odd reconstructions. It is therefore an existing two-fold
cross-validation estimate of shell-wise signal reliability. The first FSC
regularization implementation should consume the appropriate project FSC and
convert it to smooth shell-wise reliability weights for eigenvolume shrinkage.
It should not create a third holdout or a new particle split merely to validate
that shrinkage.

The project FSC can provide a cross-validated frequency-support envelope for
all fitted eigenvolumes: unsupported shells are attenuated toward zero, while
reproducible shells are retained according to an established FSC-to-weight
mapping. This is preferable to treating a nominal FSC threshold only as a hard
low-pass cutoff. In this use, the shrinkage is itself driven by independent
half-set validation already present in the project.

The selected FSC must come from the normal project/FSC ownership path and retain
its state, sampling, and shell conventions. `flex_eigenvol` should not
recalculate an FSC from a merged map or invent a separate full-data reliability
curve.

The project FSC validates frequency support of the reconstructed signal, but it
does not identify or validate individual variability modes. A later refinement
may fit even/odd eigenvolume models and calculate mode- or subspace-specific FSC
to determine whether a particular mode is reproducible. Because modes can
rotate or exchange order within a nearly degenerate subspace, that comparison
must match subspaces before comparing individual signed eigenvolumes.

### 16.4 Convergence and reproducibility diagnostics

Useful additions include objective tracking, relative residual reduction,
subspace-change criteria, controlled random seeding, and optional early
stopping. These diagnostics must work identically at the policy level for
shared and distributed execution even if floating-point reduction order
differs.

### 16.5 Reconstruction-backed state maps

For final interpretation, particle subsets or kernel weights selected from the
latent embedding should be reconstructed through the normal even/odd
reconstruction path. These maps would complement, not silently replace, the
direct linear trajectory outputs.

### 16.6 Optional rotations and richer latent models

ICA may be useful as a model-preserving rotation after a stable PCA embedding.
Any rotation must transform the basis and particle coordinates together and
must preserve reconstruction and validation semantics. Nonlinear trajectories,
multi-state models, or joint pose/latent refinement are separate model families
and should not be added as flags inside the current linear EM-PCA kernel.

## 17. Decision Record

The simplified implementation deliberately mirrors the useful behavior of the
legacy `ppca_inmem` routine: Gaussian random initialization followed by
alternating standard PCA updates. The historical name was misleading because
the working core did not require posterior covariance or iterative noise-model
estimation.

The present workflow keeps that simple numerical policy while retaining the
modern projection-aware Fourier observation model, coupled reconstruction
updates, and Hermitian formulation. Future probabilistic work should build on
this stable baseline as a separately validated extension, not complicate the
baseline before its additional outputs are demonstrably useful.
