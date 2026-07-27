# `flex_analysis` policy

This document is the contract of the code that is executed today.  Design
work, validation still to be done, and unused worker paths belong in
`doc/refactoring_notes/flex_analysis_diffusion_map_refactoring.md`.

## 1. Scope and fixed inputs

`flex_analysis` constructs a sparse diffusion-map embedding and representative
3D pre-images from a fixed mean volume and a fixed, selected `ptcl3D` set.  It
does not refine poses, shifts, CTF parameters, the mean volume, a noise model,
or latent coordinates.

The workflow requires:

- `vol1`, a fixed mean volume;
- `oritype=ptcl3D` and exactly one input state;
- an explicit `nspace`; it is never inferred from the highest occupied
  projection index;
- at least three selected particles and at least two requested representatives;
- uniform selected-particle CTF mode (`ctf=yes` or `ctf=flip`), valid unique
  stack slots, and matching `ptcl2D`/`ptcl3D` record counts.

`nspace` defines the discrete projection grid.  Every selected `ptcl3D%proj`
must be in `1:nspace`; an out-of-range value is a hard error, never remapped or
clamped.  `ptcl2D` is checked for structural consistency but its transformation
fields are not used or changed by this workflow.

The flex-specific defaults currently installed by
`simple_flex_analysis_strategy::apply_defaults` are:

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `neigs` | 15 | diffusion-eigenpair scan limit |
| `icm` | `yes` | select a spectral prefix; `no` retains the scan result |
| `k_nn` | 100 | retained residual-image neighbours per particle |
| `nang_nbrs` | 1000 | orientation-gated candidate-particle cap |
| `lp` | 6 Å | graph-feature low-pass limit |
| `bandwidth_mode` | `ferguson` | graph-kernel bandwidth selector |
| `bandwidth_tune` | 1 | Ferguson bandwidth multiplier |
| `view_balance` | `yes` | fixed correction for uneven projection-bin occupancy |
| `npreimages` | 8 | representative state volumes |
| `preimage_mode` | `linear` | local-linear pre-image estimator |
| `preimage_ndim` | 2 | cap on local-linear coordinate dimension |

`k_nn` is clamped to at least one and `nang_nbrs` to at least `k_nn`.  The
candidate cap is a particle count, not a direction count and not the number of
stored graph edges.

## 2. Registered image contract

Feature preparation applies the inverse stored **`ptcl3D`** in-plane shift and
angle once, using the project-aware Cartesian transform path.  The resulting
image has zero in-plane shift and angle and is in the frame of the cached
zero-shift, zero-in-plane-angle mean reprojection at its integer `proj` index.
No registration search occurs during graph construction or reconstruction.

For `ctf=yes` input, the image is phase-flipped once.  For `ctf=flip` input it
is not flipped again.  The copied registered project records `ctf=flip`, keeps
the particle CTF metadata, zeroes the registered `ptcl3D` in-plane transform,
and keeps the projection direction.  Astigmatism orientation is rotated into
that canonical frame.  The original project is not edited.

The durable products are:

- `flex_registered_particles.mrcs` in shared memory, or
  `flex_registered_particles_partNNN.mrcs` in distributed registration;
- `flex_registered_particles.simple`;
- registered-particle map files.

Stacks are written incrementally with a process-level `stack_io` lifetime; the
code does not open an MRC stack for every image.  Mean reprojections are made
once on the `nspace` grid.  In a distributed registration run, workers read the
master-written mean-projection stack rather than repeating reprojection.

All currently selected pre-image paths reconstruct from
`flex_registered_particles.simple`.  They do not map coefficients or poses back
to the native project and therefore do not introduce a second interpolation.

## 3. Sparse graph and embedding

The graph feature is the masked, low-pass registered residual

```text
registered phase-flipped particle - abs(canonical CTF) * mean reprojection(proj).
```

The graph-feature low-pass and mask are not the reconstruction transfer model.
The reconstruction path uses the established CTF-aware plane preparation.
If reusable upstream sigma support is found, flex enables the existing
sigma-scaled preprocessing; flex does not estimate a new sigma or a latent
noise model.

Candidate particles are gathered by angularly ordered projection-direction
bins, up to `nang_nbrs`, and the closest `k_nn` feature distances are retained.
The directed lists are symmetrized into the shared CSR diffusion graph.  No
dense `N x nang_nbrs` candidate table is materialized.

The Gaussian kernel is `exp(-d2 / eps)`.  `median` uses the median positive
k-th-neighbour squared distance; `ferguson` uses the existing log-bandwidth
scan with a median fallback. Flex fixes the Coifman–Lafon exponent internally
to `dm_alpha=0`; `dm_alpha` is not a `flex_analysis` input.

`view_balance=yes` corrects only the known acquisition nuisance: unequal
occupancy of the existing `nspace` projection bins. For occupied bin `p` with
`n_p` selected particles, every particle receives importance weight

```text
omega_i = N / (N_occupied_bins * n_p).
```

Consequently every occupied projection bin has equal total graph measure and
the particle weights have unit mean. There is no correction strength, density
floor, or additional bandwidth to tune. The weighted reversible operator is
formed through its symmetric conjugate; `view_balance=no` recovers the restored
unweighted `dm_alpha=0` graph exactly. Conformational occupancy and feature
support within each view are deliberately retained.

The existing sparse eigensolver omits the stationary vector.  With `icm=yes`,
ICM chooses a nonempty spectral prefix.  With `icm=no`, every nontrivial
eigenpair returned by the requested `neigs` scan is retained.  The files
`flex_diffmap_coordinates.txt`, `flex_diffmap_spectrum.txt`, and
`flex_diffmap_graph.txt` record the result; the graph summary records the
fixed internal normalization and the effective `view_balance` switch.

Representative descriptors are k-medoids obtained through
`cluster_dmat(...,'kmed',...)` on the raw diffusion coordinates returned by
`embed_graph` (`lambda_q psi_q`).  They are manifold locations, not hard
single-particle reconstructions.

## 4. Implemented pre-image estimator

For every medoid, flex builds a soft particle kernel on raw diffusion
coordinates.  Its initial bandwidth is the median positive squared distance
to that medoid.  Bandwidth inflation, when required for effective sample size,
is independent per state.  The kernel is density-corrected by the current
multi-state proxy `q(i)` and normalized to unit column sum.

The default `preimage_mode=linear` uses a local linear, per-Fourier-cell
weighted least-squares fit.  For state target `z_c`, it uses the first
`d=min(preimage_ndim,nmodes)` raw diffusion coordinates and

```text
u_i = z_i - z_c
g_i = [1, u_i]
A = sum_i w_i g_i g_i^T CTF_i^2
b = sum_i w_i g_i D_i
V(z_c) = intercept(A^-1 b).
```

`D_i` and `CTF_i` are supplied by the projected-model reconstruction plane
path.  The coupled per-voxel solve has an observability floor and a small
relative ridge.  The retained intercept is grid-corrected and written as the
state volume; the fitted tangent terms are not output volumes.  This is local
regression at each target, not a global eigenvolume synthesis.

`preimage_mode=constant` remains an explicit compatibility option.  It performs
direct kernel-weighted reconstruction, using the same particle weight in the
Fourier numerator and CTF/sampling-density denominator.

Both implemented modes look for exactly one compatible state-1 FSC in the
source project's `out` segment.  When found, the existing
`fsc2optlp_sub(..., merged=.false.)` filter is applied to every written state.
If metadata, the FSC file, or its shell count is unsuitable, filtering is
skipped and a log message states why.  This is an output policy currently in
the code; it is not a per-state half-map FSC estimate.

## 5. Execution policy

`flex_analysis` uses the standard abstract strategy lifecycle
(`initialize`, `execute`, `finalize_run`, `cleanup`) and factory-selected
shared-memory, worker, and master extensions.  The public program name is
always `flex_analysis`; queue workers are invoked through SIMPLE's normal
private executable dispatch with that same program name.  There is no separate
reconstruction-only program or commander.

Shared-memory execution performs registration, graph construction, embedding,
medoid selection, and the selected pre-image reconstruction in one process.

For `nparts>1`, the master currently distributes **only registered feature
preparation**.  It writes normal text particle assignments, launches stage-1
workers, and builds `flex_registered_particles.simple` over the resulting part
stacks.  It then reads the assembled feature table and performs graph
construction, embedding, medoid selection, and the active local-linear or
constant pre-image reconstruction in the master process.

The source contains worker handlers for graph rows and a legacy
modal-residual-statistics protocol.  The master does not schedule those stages
in the current workflow.  They are not a distributed graph or distributed
local-linear reconstruction implementation and are intentionally documented
as planned work, not as production behavior.

## 6. Diagnostics and outputs

`flex_preimage_identity` is a reconstruction-contract diagnostic for a
registered project.  It compares standard `reconstruct3D` preparation with the
projected-model observation/transfer planes, then checks a constant-latent
residual reconstruction against the ordinary reconstruction.  It disables
masking, output low-pass filtering, and ML regularization for that test.

`flex_preimage_basis_ab` is a separate diagnostic for the legacy global
modal-residual implementation.  It compares raw graph eigenfunctions with a
non-centering whitened basis and inverse-transformed targets.  It is not the
production local-linear state estimator.

State one is written to `outvol` (default `flex_state_001.mrc`); additional
states use the same stem with three-digit suffixes.  The medoid/assignment and
soft-weight information is written to `flex_diffmap_preimages.txt` and
`flex_registered_particle_preimage_map.txt`.

## 7. Exclusions from this policy

The following are not implemented production claims: native-image pre-image
reconstruction or mapping, distributed graph calculation, distributed
local-linear reconstruction, per-state half maps/FSC, posterior variance maps,
frequency-dependent pre-image bandwidths.  Their required contracts and
validation are maintained in the implementation note.
