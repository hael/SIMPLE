# Projection-Aware Flex PCA

## Problem

Given particles with fixed poses from a consensus refinement, estimate a
low-dimensional description of continuous structural variability and convert
it into discrete state maps that can seed
[multi-state refinement](refine3d_states.md).

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

**Mean.** A consensus map `V_mean` is loaded from `vol1`, or taken from the
project's registered consensus map when `vol1` is omitted. Residual planes
`r_i = y_i - a_i C_i S_i P_i V_mean` are formed under the same CTF, shift, and
Kaiser-Bessel conventions as [reconstruction](../reconstruction.md).

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

The contrast `a_i` is fitted in closed form against the projected mean,
`a_i = <T mu, y_i> / <T mu, T mu>`, clamped to `[0.1, 5]`.

**M-step.** The basis volumes are updated by weighted backprojection of
the residuals, coupled per Fourier voxel through the posterior second
moments,

```text
Y_q = sum_i E[z_iq] backproject_i(r_i),   solved against sum_i E[z_i z_i'] (x) D_i,
```

then re-orthonormalized into the next basis, and `Gamma` is set from the
posterior second moment. By default (`rec_backend=pcg`) the coupled solve is
a support-constrained PCG solve of `maxits_pcg=4` iterations, warm-started
from the per-voxel (gridding) solution; `rec_backend=gridding` keeps the
per-voxel solution. Even and odd particle halves maintain separate
half-bases so that agreement between them can be measured.

**Iterations.** Two fits run on disjoint particle halves, each for
`n_probe_iters` EM iterations (default 4). Their bases are frame-aligned and
merged with one joint solve on the summed statistics, and one joint EM
iteration over all particles follows from the merged basis. The only early
stop applies to a rank-1 fit: it ends once the mean principal-angle cosine
between successive bases reaches 0.999999. With two or more components the
non-reproducing tail dominates that mean, so the fit runs to the bound.
Downstream steps trust each component's reproducibility (its match cosine
between two independent half-fits), not the requested rank.

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
components standardized by their variance. The default
(`state_placement=kcenter`) places `nstates` centers by greedy farthest-point
k-center on a reliability-weighted diffusion map of all retained components,
which covers a continuous reaction coordinate and branched compositional
states with the same constants. The centers then initialize a
tied-covariance Gaussian mixture whose responsibilities become the state
weights. `state_placement=equal_occ` instead cuts a reliability-ordered path
through the latent space into `nstates` slices of equal particle count,
taking each slice's mean over all components as the target, so every state
gets the same occupancy; use it when the latent clusters are not well
separated. The same path is the fallback when the diffusion map cannot be
built. `state_axis > 0` cuts equal-occupancy slices along that single
component, and `state_axis < 0` places centers along a density-spread path.
Equal-occupancy placements skip the mixture refit, which would pull the
means onto the dominant mode, and weight particles with the kernel below,
measured along the path.

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

**Population floor.** With `min_state_frac > 0` every delivered state must
hold at least that fraction of the embedded particles. Targets are placed on
the retained particles round by round: clusters below the floor leave the
placement mass together with the particles outside every kernel support, and
the provisioned count is raised by the deficit, until `nstates` clusters
qualify. The most populated qualifying clusters are delivered, surplus
qualifying clusters attach to the nearest delivered target, and the dropped
particles receive a uniformly random delivered label. The delivered maps are
ordinary reconstructions of the labelled particles, and the floor cannot be
combined with the merge or with external targets.

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

- EM fit and initializer: `src/main/flex/fit/simple_flex_probe_fit*.f90`,
  `simple_flex_pca_basis.f90` and `simple_flex_pca_posterior.f90`.
- Driver, targets, kernel weights and merging: `src/main/flex/fit/simple_flex_pca_fit_driver.f90`
  and `src/main/flex/states/simple_flex_pca_{targets,weights,merge}.f90`.
- Weighted reconstruction: `src/main/flex/states/simple_flex_pca_rec3D.f90`.
- Projection and backprojection operators:
  `src/main/flex/fit/simple_flex_reconstructor_latent_ops.f90`.
- Subsystem overview: `src/main/flex/README.md`.
