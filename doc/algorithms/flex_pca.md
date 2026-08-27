# Projection-Aware Flex PCA

## Problem

Infer a low-dimensional description of continuous structural variability from
fixed-pose particle images, then convert selected locations in that latent
space into weighted 3D state initializers for multi-state refinement.

The production program is `flex_pca`. The earlier diffusion-map
`flex_analysis` workflow was removed because it did not produce usable states;
generative per-particle volumes remain design work, not this algorithm.

## Observation model

Let the conformational volume be approximated by

```text
V(z_i) = V_mean + sum_a z_ia U_a,
```

where `U_a` are 3D covariance eigenvolumes and `z_i` is particle `i`'s latent
coordinate. The particle observes only a CTF-modulated projection at its fixed
pose, so ordinary PCA on image pixels would mix viewing geometry with
conformation.

## Mean and covariance probes

The workflow first estimates or loads a consensus mean volume. Residual
particle planes are formed against projections of that mean under the same CTF,
shift, and KB interpolation conventions used for reconstruction.

Rather than materializing the full 3D covariance, SIMPLE selects informative
Fourier covariance columns. Even and odd particle halves accumulate column
numerators and sampling operators independently. Regularized halfset merging
and generalization diagnostics produce real 3D column representatives.

Those representatives are orthonormalized and define a reduced subspace. A
packed reduced covariance solve yields eigenvalues and a low-rank real
eigenbasis. This is a projected-covariance estimator: the viewing and CTF
operators appear in the solve rather than being treated as missing image pixels.

## Latent embedding

For each particle, projections of the mean and eigenvolumes define a small
linear inverse problem. A contrast-aware MAP solve estimates `z_i` and its
precision while accounting for the particle's pose, CTF, and available noise
model. Embedding artifacts are versioned and may be reused only when their
particle identity and model provenance match.

## State targets and weights

State centers are selected in the reliable latent subspace. The production
driver supports deterministic target strategies and bandwidth selection; for a
center `t_s`, kernel weights decrease with latent distance while widening as
needed to satisfy a minimum effective sample size. A particle may contribute
fractionally to several state maps at this stage.

Each state's even and odd Fourier accumulators are reconstructed in one bounded
particle pass. Because kernel weights in `[0,1]` make sampling density sparse and
irregular, restoration applies a shell-relative density floor before division.
FSC-compatible half products and merged state volumes are written.

The winning state weight also becomes a hard `ptcl3D/state` label in the run's
private project; unassigned particles remain state 0. This label and the state
volumes can initialize `refine3D_multi`, where ordinary hard-assignment
multi-state refinement takes over.

## State merging

An optional two-gate agglomerative merge removes over-provisioned states only
when both view-coverage and map-similarity criteria permit it. The resulting
count becomes the effective `nstates`; merging is not inferred from latent
distance alone.

## Implementation

- Driver and embedding cache: `src/main/flex/simple_flex_pca_model.f90`.
- Covariance columns, reduced solve, and embedding:
  `src/main/flex/simple_flex_pca_columns.f90`.
- Weighted reconstruction: `src/main/flex/simple_flex_pca_rec3D.f90`.
- State merge: `src/main/flex/simple_flex_pca_merge.f90`.
- Distributed artifacts: `src/main/flex/simple_flex_pca_distr.f90` and
  `src/main/flex/simple_flex_pca_parts.f90`.
- Subsystem overview: `src/main/flex/README.md`.

