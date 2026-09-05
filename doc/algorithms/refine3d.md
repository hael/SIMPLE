# Iterative 3D Refinement

## Problem

Given particle images and one or more starting volumes, estimate for each
particle a projection direction, an in-plane rotation, a shift, and (for
several volumes) a state label, and estimate the volumes that best explain the
images under those poses. As in 2D this is solved by alternation: reproject
the current volumes and assign poses, then reconstruct new volumes from the
assigned poses. This chapter is one round of that alternation. Coarse-to-fine
schedules that wrap it are in [ab initio 3D](abinitio3d.md) and
[heterogeneous refinement](heterogeneous_refinement.md).

## Model

By the projection-slice theorem, the 2D Fourier transform of a projection of
volume `x` at rotation `R` is the central section of the 3D transform of `x`
perpendicular to the viewing axis. A reference for direction `R` is therefore
obtained by gathering the 3D Fourier transform on the polar grid of that
plane with a Kaiser-Bessel kernel; no real-space projection is ever formed.
The observation model and objectives are then exactly those of
[Cluster2D](cluster2d_class_averaging.md), with the reference index running
over `nspace` directions per state instead of `K` classes:

```text
X_i = CTF_i S_s R_theta P(R_j) x_state + n_i,
L_i(j, state, theta, s) = sum_k (k/sigma2_i(k)) |X_i - CTF_i [S R P x]|^2 / sum_k (k/sigma2_i(k)) |X_i|^2.
```

Projection directions are placed quasi-uniformly on the sphere by the spiral
construction

```text
h_k = -1 + 2(k-1)/(n-1),   theta_k = acos(h_k),
psi_k = psi_{k-1} + 3.6 / sqrt(n (1 - h_k^2))   (mod 2 pi),
```

restricted to the asymmetric unit of the point group. The angular spacing is
roughly `sqrt(4 pi / n)`; `nspace = 2500` gives about 4 degrees, `20000` about
1.4 degrees.

## Algorithm

One iteration:

1. **Reference preparation.** Read the current even and odd volumes, center
   them (single state, C-symmetric groups only), soft-mask at the mask radius,
   filter, and gather polar central sections for every direction. The filter is
   the FSC-derived Wiener filter `2 FSC/(1+FSC)` unless ML regularization or
   nonuniform filtering has already produced a regularized map, in which case
   no further filter is applied. At low resolution, out to 20 A but no
   further than the shell where the FSC drops below 0.95, the two halves are
   replaced by their average for matching purposes only, since they carry no
   independent information there.
   An envelope mask is never multiplied into the matching reference: removing
   density that is present in the images destroys pose discrimination under
   the Euclidean objective.
2. **Candidate preparation** (probabilistic modes). Build the probability
   table and run the balanced assignment described in
   [sampling](sampling_and_fractional_updates.md).
3. **Search and commit.** For each active particle choose one direction,
   in-plane rotation, shift, and state:
   - `shc`: visit references in random order with the previous one last, take
     the best in-plane rotation for each, and stop at the first reference that
     beats the previous score (first-improvement hill climbing).
   - `snhc_smpl`: as above but bounded to a fraction of references that grows
     with a cosine anneal from 5 percent at the first iteration to 70 percent
     at the last, with direction and rotation drawn by power sampling rather
     than argmax.
   - `neigh`: search a coarse subspace of `nspace_sub` (default 500)
     directions exhaustively, then all fine directions within `athres`
     (default 10 degrees) of the coarse peaks and of the previous direction.
   - `prob`, `prob_neigh`, `prob_state`: read the assignment from the table
     and commit it.
   Shifts are then refined by L-BFGS-B within `trs`, and optionally the
   committed `(sx, sy, theta)` is polished jointly with continuous angle
   ([continuous in-plane refinement](continuous_inplane_refinement_abinitio2D.md)).
   Multi-state search never reuses a shift seed from one state to rank
   another.
4. **Noise update.** Accumulate per-shell residual power at the committed
   poses and reduce it per group for the next iteration.
5. **Partial reconstruction.** Re-read the updated particles and insert their
   Fourier planes into per-state, per-half accumulators under the committed
   poses ([reconstruction](reconstruction.md)).
6. **Volume assembly.** Reduce the partial accumulators, blend with the
   trailing chain if fractional updates are on, restore even and odd maps,
   compute the FSC, regularize, and produce the filtered and masked references
   for the next iteration.
7. **Bandwidth update.** Set the next iteration's low-pass limit from the FSC
   of the best-resolved state at the criterion `lplim_crit` (default 0.143),
   never finer than `lpstop` or the cropped Nyquist.

Shared-memory and distributed execution differ only in where the partial sums
are reduced; the sequence above is the same.

## Resolution estimate

The FSC between the even and odd maps is corrected for the mask by phase
randomization: the unmasked FSC gives the first shell below 0.8, both halves
are phase-randomized beyond it, and for shells two beyond that onset the
reported curve is

```text
FSC = (FSC_masked - FSC_randomized) / (1 - FSC_randomized),
```

clamped to `[0, 1]`. This removes the correlation that the mask itself
induces between the halves. Resolution is quoted at FSC 0.143 and 0.5. A
directional diagnostic, the conical FSC area ratio (cFAR), computes the FSC
inside 256 cones of half-angle 20 degrees and reports the ratio of the
smallest to the largest shell-weighted area under the curve: 1 is isotropic,
values toward 0 indicate preferred orientation.

## ML regularization

With `ml_reg=yes` (Euclidean objective only), the restored half map is the
MAP estimate under a Gaussian signal prior whose per-shell power `tau2` is
inferred from the FSC:

```text
SSNR(k) = FSC(k) / (1 - FSC(k)),      sigma2(k) = 1 / mean_shell(D),
tau2(k) = SSNR(k) sigma2(k),
x_hat(q) = B(q) / [ D(q) + 1/(tau tau2(k)) ],
```

where `B` and `D` are the noise-weighted numerator and CTF-squared density of
[reconstruction](reconstruction.md), and `tau` (default 1) is a fudge factor.
Shells below index 6 are not regularized. Because `1/tau2` grows as the FSC
falls, the estimator shrinks poorly determined high-frequency coefficients
toward zero in proportion to their expected signal, which is the same
operation as the Wiener filter but applied before the density division rather
than after, where it also regularizes the division itself.

## Automasking and nonuniform filtering

Late in a refinement, the reference volume for iteration `n` may be
multiplied by an envelope computed at iteration `n-1`, either a density
automask (low-pass at `amsklp`, Otsu threshold, largest connected component,
one-voxel dilation, 6-voxel cosine edge) or the NU-evidence envelope
([NU-evidence envelope masking](nu_evidence_envelope_mask.md)), and its
bandwidth may be made spatially varying ([nonuniform
filtering](nonuniform_filtering.md)). These operate on the reference only; the
reconstructed halves and the FSC inputs are unchanged.

## Convergence

Per particle the search records `mi_proj` (1 if the direction moved less than
2 degrees), the visited fraction `frac`, and the running mean angular change.
A single-state run has converged when the average `mi_proj` exceeds 0.99 and
the visited fraction exceeds 99 percent. Multi-state runs use the joint
state-and-direction overlap (0.96). Shift search is switched on once the
visited fraction exceeds 75 percent, with `trs = clamp(0.07 msk, 5, 6)`
pixels. A `combine_eo=yes` run ends with one full-update iteration at a
tighter band with the halves merged.

## Rationale

- Working entirely in Fourier space makes reprojection a gather and
  reconstruction a scatter with the same kernel, so the two halves of the
  alternation are adjoint operations on the same representation.
- First-improvement stochastic search with a random visiting order is cheap
  and explores; exhaustive neighborhood search late is accurate. The
  probabilistic mode combines both by decoupling scoring from assignment.
- The FSC drives everything downstream (filter, prior, bandwidth), so it is
  computed on masked, phase-randomization-corrected halves that were
  reconstructed independently, and never from references that have been
  blended for matching.

## Implementation

- Commander: `src/main/commanders/simple/simple_commanders_refine3D.f90`.
- Iteration strategy: `src/main/strategies/parallelization/simple_refine3D_strategy.f90`.
- Reference preparation: `src/main/strategies/search/simple_matcher_refvol_utils.f90`;
  polar gather in `src/main/image/simple_projector_pft_batch.f90`.
- Search strategies: `src/main/strategies/search/simple_strategy3D_*.f90`.
- Volume assembly, FSC, regularization:
  `src/main/commanders/simple/simple_commanders_rec_distr.f90`,
  `src/main/volume/simple_halfmap_diagnostics.f90`, `src/utils/filter/simple_fsc.f90`.
- Convergence: `src/main/simple_convergence.f90`.
- Policy: `doc/policies/refine3D_policy.md`.
