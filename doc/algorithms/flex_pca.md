# Projection-Aware Flex PCA

## Problem

Given particles with fixed poses from a consensus refinement, estimate a
low-dimensional description of continuous structural variability and convert
it into discrete state maps that can seed
[multi-state refinement](heterogeneous_refinement.md).

Ordinary PCA on the images would mix viewing geometry with conformation: two
identical molecules seen from different directions differ more in pixels than
two conformations seen from the same direction. The variability has to be
modeled in the volume domain and observed through the projection operator.

## Model

The conformational volume of particle `i` is a linear combination of a mean
and `r` basis volumes with latent coordinate `z_i` in `R^r`:

```text
V(z_i) = V_mean + sum_q z_iq U_q,      z_i ~ N(0, Gamma),  Gamma diagonal.
```

The particle observes a whitened, CTF-modulated, shifted central section of
that volume with an unknown per-particle contrast `a_i`:

```text
y_i = a_i C_i S_i P_i (V_mean + U z_i) + n_i,      n_i ~ N(0, sigma2 I),
```

in the band-limited Fourier domain, after per-shell whitening. This is
probabilistic PCA with a per-observation linear operator `a_i C_i S_i P_i`, and
it is fitted by expectation-maximization. The projection and CTF operators
appear inside the E-step rather than being treated as missing pixels.

## Algorithm

**Mean.** A consensus map `V_mean` is reconstructed from all particles (or
loaded). Residual planes `r_i = y_i - a_i C_i S_i P_i V_mean` are formed under
the same CTF, shift, and Kaiser-Bessel conventions as
[reconstruction](reconstruction.md).

**Initial basis.** The starting subspace is data-free: the lowest-frequency
Fourier lattice points admitted by the band, chosen greedily with a minimum
separation so they are not neighbors, each realized as a band-limited, masked,
deapodized cosine/sine pair, then orthonormalized. Two calibration scalars
are estimated from the data before iterating: the whitened noise level
`sigma2` from high-frequency shells and an initial prior variance `Gamma^0`,
deliberately over-estimated so the first-iteration prior is weak. Starting
from geometry rather than from a data moment makes two runs comparable and
avoids seeding the EM inside a noise subspace.

**E-step.** For each particle, with `G_i` the `r x r` Gram matrix of the
projected basis at that particle's pose and `b_i` the projected-basis inner
products with the residual,

```text
A_i        = (a_i^2 / sigma2) G_i + Gamma^{-1},
E[z_i]     = A_i^{-1} (a_i / sigma2) b_i,
E[z_i z_i'] = E[z_i] E[z_i]' + A_i^{-1}.
```

The contrast `a_i` is fitted inside a bracket `[0.2, 3]`, in closed form or
by a grid.

**M-step.** The basis volumes are updated by weighted backprojection of
the residuals, coupled per Fourier voxel through the posterior second
moments,

```text
Y_q = sum_i E[z_iq] backproject_i(r_i),   solved against sum_i E[z_i z_i'] (x) D_i,
```

then re-orthonormalized into the next basis, and `Gamma` is set from the
posterior second moment. Even and odd particle halves maintain separate
half-bases so that agreement between them can be measured.

**Convergence.** The iteration stops when the principal angles between
successive subspaces exceed a cosine of 0.97, or when the number of
reproducible dimensions between the even and odd half-bases has plateaued
(tolerance 0.02 for three iterations). The reproducible dimension count, not
the requested rank, is what the downstream steps trust.

**Embedding.** With the converged basis, each particle's latent coordinate is
the MAP solution of the same E-step,

```text
z_i = A_i^{-1} (a_i b_i - a_i^2 c_i) / sigma2,
```

with its posterior precision `Pi_i = A_i` retained as a per-particle
uncertainty. Embeddings are cached and reused only when particle identity and
model provenance match.

## From latent coordinates to states

**Targets.** State centers `t_s` are placed in the reliable subspace, with
components standardized by their variance. The default places `nstates`
centers at equal-occupancy quantiles along a chosen component (or along a
reliability-ordered path), taking each slice's mean over all components as
the target; k-means and diffusion k-center are alternatives.

**Kernel weights.** Particle `i` contributes to state `s` with an
Epanechnikov weight in the particle's own posterior metric,

```text
d_i(s) = (z_i - t_s)' Pi_i (z_i - t_s),
w_i(s) = max(0, 1 - d_i(s) / h_s^2),
```

so an uncertain particle is spread over neighboring states rather than
assigned sharply. The bandwidth `h_s` is set from the sorted distances so
that the effective sample size `(sum w)^2 / sum w^2` reaches a minimum
(default at least 20, and by default the count needed for a stable map), and
grown by 30 percent steps if the support is still too small.

**State maps.** Each state's even and odd accumulators are reconstructed in
one weighted pass through the gridding reconstructor. Because weights in
`[0, 1]` make the sampling density sparse and irregular, a shell-relative
density floor is applied before division. FSC-compatible half maps and merged
maps are written, and the maximum-weight state becomes each particle's hard
label (unweighted particles stay at state 0), so the initializer can be judged
with an ordinary multi-state reconstruction.

**Merging.** Over-provisioned states are merged only when two independent
gates agree: a view-coverage gate (a state whose viewing-direction second
moment is an outlier relative to the others, chi-square with 5 degrees of
freedom on the effect size `chi2/n_eff`, robustly compared by median and MAD)
and a map-similarity gate (a disattenuated FSC-type ratio between the states'
deviations from the ensemble mean above 0.98 by default). Latent distance
alone never triggers a merge, since it is a property of the embedding, not
of the maps.

## Rationale

- Keeping the projection operator in the likelihood is what separates
  conformational variance from viewing variance; the Gram matrix `G_i`
  encodes exactly how much of each basis volume the particle's view can see.
- The EM alternation is the same structure as refinement, with the discrete
  pose replaced by a continuous latent coordinate and the volume replaced by a
  basis. Even/odd half-bases give it the same reproducibility test.
- Posterior-metric kernel weighting turns the embedding's uncertainty into
  soft state membership, and the effective-sample-size floor guarantees that
  each state map has enough particles to be a usable reference.

## Implementation

- EM fit and initializer: `src/main/flex/simple_flex_pca_em*.f90`.
- Driver, targets, kernel weights, merging: `src/main/flex/simple_flex_pca_model.f90`,
  `src/main/flex/simple_flex_pca_merge.f90`.
- Weighted reconstruction: `src/main/flex/simple_flex_pca_rec3D.f90`.
- Projection and backprojection operators:
  `src/main/flex/simple_flex_projected_latent_model.f90`,
  `src/main/flex/simple_flex_reconstructor_latent_ops.f90`.
- Subsystem overview: `src/main/flex/README.md`.
