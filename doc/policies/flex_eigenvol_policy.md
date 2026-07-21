# `flex_eigenvol` Policy: Sparse Diffusion-Map 3D Variability

This policy records the implemented scientific and software contract for
`flex_eigenvol`. Forward-facing developments are separated explicitly from
the production behavior.

## 1. Production model

`flex_eigenvol` is a diffusion-map workflow. It does not fit PCA, probabilistic
PCA, alternating E/M updates, a Gaussian latent prior, an external sigma array,
or a model noise variance.

The implemented input assumptions are:

- one supplied and fixed mean volume (`vol1`);
- one active `ptcl3D` state;
- fixed 3D orientations, in-plane rotations, shifts, and CTF parameters;
- `ctf=yes` or already phase-flipped (`ctf=flip`) particle images;
- one low-pass limit `lp` used for both graph features and 3D mode
  reconstruction.

The workflow returns a continuous low-dimensional coordinate system and up to
20 corresponding 3D mode volumes. A mode volume is a linear reconstruction of
one selected diffusion coordinate; it is not a PCA eigenvector and its
diffusion eigenvalue is not a particle-image variance.

## 2. Canonical registered frame

Each active particle is normalized and transformed once by the inverse of its
stored in-plane shift and in-plane angle. After this transformation it is
directly comparable with the discrete mean-volume reprojection identified by
its `proj` index; no rotation or shift search is performed during graph
construction.

The registered particles are phase-flipped when the input project is
`ctf=yes`. Already phase-flipped particles are not flipped again. The durable
outputs are:

- `flex_registered_particles.mrcs`;
- `flex_registered_particles.simple`;
- `flex_registered_particle_map.txt`.

The registered project retains particle CTF parameters, records the stack as
`ctf=flip`, maps particles to their registered stack entries, and stores zero
in-plane angle and shift. Parameters obtained later in this frame must be
mapped back to raw particles using the established
`denoise_project`/`map_params_from_den` composition convention.

## 3. Graph feature and sparse candidate policy

For particle `i`, the graph feature is the masked, band-limited registered
residual

```text
r_i = phaseflip(register(y_i)) - |CTF_i| register(P_proj(i) mean)
```

serialized in Cartesian image space. The comparison is therefore
CTF-consistent in the phase-flipped observation model.

Particles are binned by the finite discrete projection index. For every
source direction, projection bins are traversed in angular order. At most
`nang_nbrs` candidate particles are compared and only the `k_nn` closest
registered-residual features are retained. The default public parameters are:

```text
k_nn      = 10
nang_nbrs = 1000
```

`nang_nbrs` is a candidate cap, not the number of retained graph edges. It is
never represented as a dense `N x nang_nbrs` table. The retained directed
neighbors are symmetrized into the existing CSR `diffmap_graph` and normalized
with the existing diffusion-map policy.

Angular proximity gates which comparisons are meaningful; registered image
distance decides which gated particles are neighbors.

## 4. Embedding and feature selection

The normalized sparse graph is embedded by the existing sparse diffusion-map
eigensolver. The trivial stationary vector is omitted. Nontrivial eigenpairs
are returned in descending eigenvalue order, and the coordinate convention of
`embed_graph` is retained.

ICM selects a prefix of the spectrum with a minimum rank of one and a maximum
of `min(neigs,20)`. Zero selected modes is invalid. No scientific assumption is
made about a preferred ordering beyond the descending order required by the
current prefix-based ICM selector.

The coordinate and spectrum outputs are:

- `flex_diffmap_coordinates.txt`;
- `flex_diffmap_spectrum.txt`;
- `flex_diffmap_graph.txt`.

## 5. Hermitian 3D reconstruction

The Hermitian Fourier reconstruction formulation is required and must be
retained. After subtracting the fixed mean projection, diffusion mode `q`
accumulates particle `i` with

```text
numerator scale = z_iq
density scale   = z_iq**2
```

through `insert_plane_oversamp_multi_scaled`. The numerator is complex and
Hermitian; the density is real and uses the particle CTF-squared sampling
density. Modes are independent, so reconstruction cost and memory are linear
in the selected rank rather than quadratic in the number of modes.

For phase-flipped input, `gen_fplane4rec` applies `abs(CTF)` to the numerator
and `CTF**2` to the density. The resulting maps are CTF-density-corrected
Hermitian least-squares reconstructions. They must not be described as
unit-transfer backprojections or as statistically Wiener-restored maps.

The configured mask and `lp` are applied during finalization. There is no
post-reconstruction Gram rotation, PCA orthogonalization, trajectory-volume
fan-out, duplicate difference-map output, or sigma-derived scaling.

## 6. Shared- and distributed-memory execution

The shared path performs registered feature preparation, sparse graph
construction, embedding, ICM selection, and Hermitian reconstruction in the
main process, using the established OpenMP kernels.

The distributed path has the same scientific result. Registration, graph
construction, embedding, and ICM remain master-owned because the centered
sparse eigensolve is currently global. Reconstruction particles are divided
among `flex_eigenvol_reconstruct` workers. Each worker writes complex
Hermitian numerators and real density statistics for its particle interval;
the master sums those statistics and performs the same density correction,
masking, filtering, and volume output as the shared path.

The former `flex_eigenvol_mstep` and `flex_eigenvol_estep` programs and their
iterative PCA state are not part of this policy.

## 7. Output and logging policy

Standard output is stage-oriented and concise. It reports:

- particle count, `k_nn`, `nang_nbrs`, `lp`, and maximum rank;
- registered-feature size and preparation time;
- candidate-count range/mean, graph edge count, and graph time;
- selected ICM rank and convergence summary;
- reconstruction progress and stage time;
- registered stack/project names and total time.

The requested `outvol` is mode 1. Further selected modes use the same stem with
three-digit suffixes. Only selected mode volumes are written.

## 8. Implemented versus forward-facing work

| Capability | Status |
| --- | --- |
| Registered Cartesian particle stack and project | Implemented |
| Discrete mean reprojection by `proj` | Implemented |
| Phase-flipped, CTF-matched registered residual | Implemented |
| Angular candidate cap and sparse residual kNN graph | Implemented |
| Sparse diffusion embedding | Implemented |
| ICM prefix rank selection, rank at least one | Implemented |
| Hermitian CTF-density 3D mode reconstruction | Implemented |
| Shared-memory execution | Implemented |
| Distributed reconstruction-statistics workers | Implemented |
| Iterative PCA/PPCA and sigma estimation | Removed from flex workflow |
| Distributed feature and graph preparation | Forward-facing optimization |
| Kernel- or neighborhood-weighted state reconstructions along the diffusion manifold | Forward-facing development |
| Posterior covariance/model-variance map | Forward-facing development |
| Project-FSC shrinkage of mode volumes | Forward-facing development |

Kernel- or neighborhood-weighted state reconstructions along the diffusion
manifold are a valuable later development but are outside this first
replacement.

A spatial variance map derived from posterior covariance or model variance
would also be valuable, but it requires a clearly specified likelihood and
calibrated uncertainty model. It must not be inferred from diffusion
eigenvalues alone.

Project FSCs are already two-fold cross-validated. A future shrinkage method
should use those independently validated FSC estimates directly—preferably
mode-specific FSCs when available—to shrink Fourier coefficients continuously.
This is distinct from applying a hard or soft FSC cutoff as an output filter.
