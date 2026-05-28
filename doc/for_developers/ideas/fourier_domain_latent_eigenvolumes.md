# Fourier-Domain Latent Eigenvolumes for 3D Variability

**Date:** May 27, 2026  
**Status:** Design note / first-draft implementation plan

## Purpose

This note sketches a SIMPLE-native route for continuous 3D heterogeneity
analysis using Fourier-domain eigenvolumes.

After `abinitio3D` or an early `refine3D` stage, SIMPLE has the ingredients
needed for a projection-aware latent model:

- a consensus 3D map or half maps
- per-particle 3D orientation
- per-particle origin shift
- per-particle CTF parameters
- particle Fourier data in the matching/reconstruction shell range

The central idea is to model structural variability by a small set of complex
Fourier eigenvolumes, rather than by applying PCA directly to raw particle
Fourier vectors. Particle images are pose-dependent 2D samples of a 3D object,
so the latent basis must live in a common 3D Fourier coordinate system.

The long-term target is a continuous latent variable model:

```text
V_i(k) = mu(k) + z_i1 B_1(k) + z_i2 B_2(k) + ... + z_iQ B_Q(k)
```

where `mu` is the consensus Fourier volume, `B_j` are Fourier-domain
eigenvolumes, and `z_ij` are particle latent coordinates. The observation model
uses SIMPLE's existing pose, shift, CTF, and noise-weighting machinery to map
the eigenvolumes to each particle image.

## Current Code Shape

The relevant existing implementation is concentrated in:

- `src/main/pca/simple_complex_ppca.f90`
  - complex/Hermitian PPCA
  - streaming sufficient-statistics accumulation
  - MMSE shrinkage in a complex eigenbasis
  - currently assumes fixed feature coordinates across observations
- `src/main/strategies/search/simple_polarft_lines_ppca_stream.f90`
  - streams prepared particle polar Fourier lines into `complex_ppca`
  - useful scaling pattern, but line indices are not global 3D coordinates
- `src/main/pftc/simple_polarft_access.f90`
  - `get_ptcl_pft`
  - `get_ptcl_line`
  - `get_coord`
  - `get_polar_coords`
- `src/main/image/simple_image_ctf.f90`
  - `gen_fplane4rec`
  - applies origin-shift phase, CTF, and optional sigma weighting
  - outputs `fplane_type` with complex Fourier plane plus CTF squared plane
- `src/main/strategies/search/simple_matcher_3Drec.f90`
  - `prep_imgs4rec`
  - `update_rec`
  - maps prepared particle Fourier planes into the 3D reconstruction lattice
- `src/main/volume/simple_reconstructor.f90`
  - `insert_plane_oversamp`
  - `insert_plane_oversamp_opt`
  - `project_polar`
  - already contains projection and backprojection-like operations between
    2D Fourier planes and 3D Cartesian Fourier volumes
- `src/main/pftc/simple_polarft_core.f90`
  - `vol_pad2ref_pfts`
  - `vol_pad2ref_pfts_opt`
  - extracts polar central sections from 3D volumes for matcher references

The important caveat is that `complex_ppca` is currently a fixed-coordinate
PPCA model. It is appropriate when every observation vector represents the same
feature coordinates. Raw particle PFT entries do not satisfy that requirement:
the same `(irot,k)` index corresponds to different 3D Fourier coordinates for
different particle orientations.

Therefore, the heterogeneity engine should reuse the complex PPCA ideas, but
it must wrap them in projection and backprojection operators.

## Proposed Observation Model

Let particle `i` have orientation `R_i`, shift `t_i`, CTF `C_i`, and latent
coordinate vector `z_i`.

For a 2D Fourier sample `q` in the particle plane:

```text
y_i(q) = H_i(q) * [ mu(R_i q) + sum_a z_i(a) B_a(R_i q) ] + eps_i(q)
```

where:

- `mu(k)` is the consensus Fourier volume
- `B_a(k)` is eigenvolume `a` in the same Fourier grid
- `H_i(q)` includes CTF, shift phase, optional phase-flip convention, and
  optional ML/sigma weighting
- `eps_i(q)` is noise, preferably shell-dependent or particle/shell dependent

In the reconstruction code path, much of `H_i` is already applied when
`gen_fplane4rec` builds `fplane%cmplx_plane` and `fplane%ctfsq_plane`. A first
implementation should reuse those prepared planes rather than inventing another
CTF/shift convention.

For physical interpretability, the recommended first model is:

```text
complex Fourier observations
complex Fourier eigenvolume coefficients
real particle latent coordinates z_i
```

Complex latent coordinates should be avoided initially. They can absorb phase,
shift, and pose errors in ways that are hard to distinguish from structural
variability.

## Why Not PCA on Particle PFTs

It is tempting to build a particle vector by concatenating all polar Fourier
samples and then run `complex_ppca` directly. That is not a 3D variability
model.

For particle `i`, a polar sample corresponds to:

```text
k_3d = R_i q
```

For particle `j`, the same vector index generally corresponds to:

```text
k_3d = R_j q
```

Those are not the same Fourier voxel. A fixed-coordinate PCA over raw PFT
vectors would mix:

- structural variability
- orientation distribution
- CTF differences
- shift or pose errors
- missing-wedge and sampling-density effects

The correct common coordinate system is the 3D Fourier volume, not the raw
particle vector.

## Core Operations

A projection-aware latent eigenvolume engine needs four operations.

### Project Eigenvolumes to a Particle

Given a particle orientation and shift/CTF convention, project `mu` and every
`B_a` onto that particle's Fourier plane:

```text
P_i0(q) = H_i(q) * mu(R_i q)
P_ia(q) = H_i(q) * B_a(R_i q)
```

This is analogous to `project_polar` and `vol_pad2ref_pfts`, but the target is
the actual particle orientation rather than a fixed reference grid.

### Infer Particle Latents

For each particle, solve a small weighted least-squares problem:

```text
r_i(q) = y_i(q) - P_i0(q)
r_i(q) ~= sum_a z_i(a) P_ia(q)
```

The normal equations are only `Q x Q` per particle:

```text
A_i(a,b) = sum_q w_i(q) real( conjg(P_ia(q)) * P_ib(q) )
b_i(a)   = sum_q w_i(q) real( conjg(P_ia(q)) * r_i(q) )
z_i      = solve(A_i, b_i)
```

This step is cheap once the basis projections are available.

### Backproject Weighted Residuals to Eigenvolumes

Given current `z_i`, update each eigenvolume by gridding residual contributions
back into the common 3D Fourier lattice:

```text
res_i(q) = y_i(q) - H_i(q) * [ mu(R_i q) + sum_a z_i(a) B_a(R_i q) ]
```

For each eigenvolume `a`, grid a contribution proportional to:

```text
z_i(a) * conjg(H_i(q)) * res_i(q)
```

The associated density/normalizer should account for:

```text
z_i(a)^2 * |H_i(q)|^2
```

Cross terms between latent dimensions can be handled either by block normal
equations per voxel or by alternating updates. A first implementation can avoid
full per-voxel `Q x Q` storage by updating one eigenvolume at a time or by
keeping `Q` very small.

### Enforce Fourier-Volume Constraints

Each eigenvolume should satisfy the same physical constraints as a real-space
volume:

- Friedel symmetry
- shell-limited support
- optional masking or soft real-space support after inverse FFT
- half-set validation
- scale/sign conventions for stable interpretation

Eigenvolume signs are arbitrary, and eigenvolume order can swap when spectra
are close. Diagnostics should not assume stable signs between runs without
alignment.

## Relation to `simple_complex_ppca`

`simple_complex_ppca` is a useful seed, but not a finished heterogeneity
engine.

What can be reused directly:

- Hermitian covariance thinking
- complex mean/eigenbasis representation
- streaming sufficient-statistics style
- shrinkage/noise-aware interpretation
- careful complex arithmetic and eigenvalue ordering

What must change:

- features cannot be raw particle-vector indices
- observations are incomplete projections of a 3D object
- particle latent coordinates should be real
- CTF/shift/noise operators must be part of the model
- the basis should be stored as 3D Fourier volumes

In other words, the right next object is not a subclass of the current
fixed-vector `complex_ppca`. It should be a projection-aware sibling, for
example:

```text
simple_fourier_latent_eigenvolumes
```

or:

```text
simple_complex_eigenvol_ppca
```

## First-Draft Implementation Plan

### Draft 0: Residual Eigenvolume Smoke Test

Purpose: verify that SIMPLE can accumulate and visualize meaningful 3D
Fourier residual variability with fixed poses.

Implementation:

1. Start from a project after `abinitio3D` or a low-resolution `refine3D`.
2. Reconstruct or read a consensus `mu`.
3. Reuse `prep_imgs4rec` style preprocessing to produce `fplane_type` objects.
4. For each particle, subtract the projection of `mu` from the prepared Fourier
   plane.
5. Grid residual planes into diagnostic volumes by simple particle subsets:
   - random halfsets
   - high/low score groups
   - existing state/class labels
6. Write residual maps and half-map correlations.

Expected result:

- If the residual volume signal is dominated by noise, do not proceed to a
  high-rank model.
- If known compositional or conformational differences appear in residual maps,
  the projection-aware path is viable.

### Draft 1: One-Dimensional Latent Eigenvolume

Purpose: fit one continuous eigenvolume and one scalar latent coordinate per
particle.

Suggested alternating scheme:

1. Initialize `B_1`.
   - use difference between two existing states
   - use residual maps from Draft 0
   - or use a low-resolution random volume with half-set constraints
2. Infer `z_i` for each particle by weighted least squares.
3. Reconstruct `B_1` from residuals weighted by `z_i`.
4. Normalize `B_1` and rescale `z_i`.
5. Repeat for a small number of iterations.

Keep all work low resolution at first. The aim is to test the algebra and
storage model, not to maximize final map resolution.

### Draft 2: Low-Rank Model With Small Q

Purpose: extend from one latent dimension to `Q=2..8`.

Implementation details:

- Infer `z_i(1:Q)` with per-particle `Q x Q` normal equations.
- Update eigenvolumes either:
  - jointly with per-voxel small normal systems, or
  - sequentially with residual deflation
- Orthogonalize eigenvolumes under a sampling-density weighted Fourier inner
  product.
- Track eigenvalue-like variance explained by each dimension.

The first joint implementation can prioritize clarity over maximum speed.
Avoid introducing a large dense covariance matrix over 3D voxels.

### Draft 3: Half-Set Validation and Regularization

Purpose: prevent the latent model from fitting noise.

Implementation:

- Fit even and odd eigenvolume models separately or update from half-specific
  particles.
- Compare eigenvolume half maps by shell.
- Apply FSC-like shrinkage to eigenvolumes.
- Report per-dimension half-set reproducibility.
- Reject or downweight dimensions whose half-set agreement is poor.

This is the analogue of the FSC regularization used for homogeneous maps. It is
also where SIMPLE can avoid a major failure mode of flexible latent models:
plausible but non-reproducible motion.

### Draft 4: Particle Embedding and Map Generation

Purpose: turn the model into an interpretable heterogeneity tool.

Outputs:

- particle latent table:

```text
particle_index z1 z2 ... zQ state eo score optional_local_resolution
```

- eigenvolume stack:

```text
mu.mrc
eigvol_001.mrc
eigvol_002.mrc
...
```

- trajectory maps:

```text
mu + alpha * B_a
```

- clustered or kernel-regressed maps from particles near selected latent
  coordinates

Trajectory maps are useful for inspection, but final high-resolution state maps
should preferably be reconstructed from particle subsets or kernel weights, not
only by direct linear reprojection of low-SNR eigenvolumes.

### Draft 5: ICA or Disentangling Layer

Purpose: separate statistically independent motions after a stable latent
embedding exists.

Recommended first implementation:

1. Fit a low-rank Fourier eigenvolume model.
2. Whiten particle coordinates `z`.
3. Run an in-house ICA rotation on the `N x Q` coordinate matrix.
4. Rotate the eigenvolume basis by the same transform.
5. Validate rotated components by half-set reproducibility and subset
   reconstructions.

ICA should be treated as a coordinate-system choice on top of the latent model,
not as the primary reconstruction engine.

## Storage and API Sketch

A first module could expose an object like:

```fortran
type :: fourier_latent_eigenvol_model
    integer :: q
    integer :: kfromto(2)
    type(image) :: mean_vol
    type(image), allocatable :: eigvols(:)
    real, allocatable :: z(:,:)       ! q x nptcls
    real, allocatable :: dim_power(:)
contains
    procedure :: new
    procedure :: init_from_mean
    procedure :: infer_latents
    procedure :: update_eigenvolumes
    procedure :: normalize_basis
    procedure :: write_model
    procedure :: kill
end type
```

Implementation ownership should stay aligned with current SIMPLE boundaries:

- `src/main/pca/`
  - model object and low-rank latent algebra
- `src/main/strategies/search/`
  - particle batch preparation and pose-aware projection/backprojection glue
- `src/main/volume/`
  - reusable Fourier volume projection/backprojection helpers, if needed
- `src/main/commanders/`
  - experimental workflow command
- `src/main/ui/`
  - UI registration only after the prototype is useful

Avoid making `refine3D` depend on the new model initially. A standalone
analysis command is safer.

## Numerical Policy

Recommended defaults for a prototype:

- low-pass limited input, for example early-stage `abinitio3D` resolution
- small `Q`, preferably `1..4`
- fixed poses and shifts
- no pose refinement inside the latent fit
- use existing CTF and shift conventions from `gen_fplane4rec`
- fit on active particles only
- half-set diagnostics from the start
- write enough intermediate maps to diagnose failure

Potential failure modes:

- pose errors appear as eigenvolumes
- residual CTF errors appear as shell/ring components
- shift errors appear as phase-ramp components
- preferred orientation creates anisotropic eigenvolumes
- too many dimensions fit noise
- direct trajectory maps overstate high-resolution detail

## Open Questions

- Should the first implementation operate on Cartesian `fplane_type` planes
  from reconstruction, or on PFTC polar sections with explicit 3D coordinate
  lookup?
- Should eigenvolume updates use the existing reconstructor gridding machinery
  directly, or a new accumulator specialized for low-rank residual updates?
- How much of `complex_ppca` should be refactored into reusable complex
  Hermitian eigensolver utilities?
- Should latent coordinates be stored in the project file, an external table,
  or both?
- What is the first target workflow: post-`abinitio3D` analysis, post-`refine3D`
  analysis, or an experimental refinement side path?

## Recommended First Step

Build a standalone experimental command that:

1. reads a project with fixed 3D poses
2. reads or reconstructs a consensus map
3. prepares particle Fourier planes with existing reconstruction code
4. fits `Q=1` at low resolution
5. writes `z_i`, `mu`, `eigvol_001`, and half-set diagnostics

This gives SIMPLE an in-house test bed for Fourier-domain latent embedding
without committing the main refinement workflow to the model.
