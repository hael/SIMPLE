# Continuous Projection-Direction Exploration in PFTC

**Date:** May 10, 2026  
**Status:** Design note / first-draft implementation plan

## Purpose

SIMPLE already has a strong continuum formulation for origin-shift registration in polar Fourier space. The current 3D projection-direction search, however, remains fundamentally discrete: reference polar central sections are extracted from the 3D reference volume over an Euler-space grid, and each section is compared against the particle PFT while all in-plane rotations are evaluated efficiently.

This note sketches a path for local continuous refinement of projection directions while preserving the existing PFTC advantage: once a polar central section has been extracted, the in-plane rotational correlations are still obtained cheaply by the current FFT-based reduction. The intended use is not to replace the global discrete search immediately, but to accelerate and sharpen it by applying derivative-driven local refinement to a small set of promising discrete projection directions.

## Current Code Shape

The core projection-to-PFT extraction path is:

- `src/main/image/simple_projector_pft_batch.f90`
  - `fproject_polar_batch`
  - `fproject_polar_batch_mirr`
- `src/main/image/simple_projector_pft.f90`
  - `fproject_polar`
  - `fproject_polar_oversamp`
- `src/main/image/simple_projector.f90`
  - `interp_fcomp`
  - `interp_fcomp_oversamp`
- `src/main/pftc/simple_polarft_corr.f90`
  - `gen_objfun_vals`
  - `gen_corrs`
  - `gen_euclids`
  - `gen_corr_grad_for_rot_8`
- `src/main/pftc/simple_pftc_shsrch_grad.f90`
  - existing L-BFGS shift optimizer with discrete in-plane-angle callback
- `src/main/strategies/search/simple_strategy3D_srch.f90`
  - shift/in-plane refinement hooks after discrete projection search
- `src/main/strategies/search/simple_strategy3D_greedy.f90`
  - exhaustive projection-direction scoring over discrete references

The important detail is that extracting a reference polar central section from a volume currently reduces to:

```fortran
px     = real(polar_x(irot,kloc))
py     = real(polar_y(irot,kloc))
loc(1) = px*e_rotmat(1,1) + py*e_rotmat(2,1)
loc(2) = px*e_rotmat(1,2) + py*e_rotmat(2,2)
loc(3) = px*e_rotmat(1,3) + py*e_rotmat(2,3)
refs_state(irot,kloc,iproj) = vol_pad%interp_fcomp_oversamp(loc)
```

That means a continuous projection-direction derivative can be formulated as a derivative of the 3D Fourier interpolation location, rather than as a derivative of the polar correlation framework itself.

## Proposed Formulation

Represent a local projection-direction perturbation by two tangent-plane coordinates around the current projection direction. Let the current orientation matrix be `R0`, and let `a` and `b` be small local tangent coordinates. A perturbed orientation can be represented as:

```text
R(a,b) = exp([a*u + b*v]x) * R0
```

or by an equivalent local tilt basis compatible with SIMPLE's `ori`/Euler conventions. The key requirement is that the local parameterization avoids Euler singularities and updates on the sphere.

For each polar sample `(p,k)`:

```text
q(p,k)      = [polar_x(p,k), polar_y(p,k), 0]
loc(p,k)   = q(p,k) * R(a,b)
P(p,k)     = V_interp(loc(p,k))
dP/da      = grad_V_interp(loc) dot dloc/da
dP/db      = grad_V_interp(loc) dot dloc/db
```

where `P` is the extracted reference PFT section and `V_interp` is the existing Kaiser-Bessel interpolation from `projector%cmat_exp`.

The existing in-plane rotational reduction is retained. For each tangent direction, the derivative PFT section can be run through correlation-gradient kernels analogous to the existing shifted-reference gradient machinery.

For normalized cross-correlation, at a fixed in-plane rotation:

```text
C = N / sqrt(Sref * Sptcl)

dC = dN / sqrt(Sref * Sptcl)
     - 0.5 * C * dSref / Sref
```

where:

- `dN` is the cross term between the particle PFT and the derivative reference PFT.
- `dSref` is the derivative of the CTF-weighted reference square sum.
- `Sptcl` is fixed with respect to projection direction.

For the Euclidean objective, the derivative follows from:

```text
D = sum_k,p w_k * |CTF * P_rot - X|^2
dD = 2 * sum_k,p w_k * real((CTF * dP_rot) * conjg(CTF * P_rot - X))
```

The current code already contains closely related shift-gradient implementations in `gen_corr_cc_grad_for_rot_8` and `gen_euclid_grad_for_rot_8`. Those routines differentiate the shift phase factor. The new work would differentiate the extracted reference section instead.

## Search Policy

The recommended first policy is local refinement of shortlisted discrete candidates:

1. Run the current discrete projection-direction search.
2. Keep the top `npeaks` or `npeaks_inpl` projection candidates.
3. For each candidate, run a local tangent-plane optimizer over two projection-direction variables.
4. Keep the in-plane angle discrete.
5. During optimization, periodically reselect the best in-plane angle, mirroring the shift optimizer's `opt_angle` callback pattern.
6. Accept the continuous-refined candidate only if it improves the stored score.

This avoids asking a local optimizer to solve the global spherical search problem. The discrete grid remains responsible for global coverage; the continuous step sharpens local peaks and may eventually allow a coarser projection grid.

## First-Draft Implementation Plan

### Draft 0: Derivative-Free Tangent Probe

Purpose: establish whether the objective is locally smooth and whether continuous direction refinement improves useful particles.

Implementation:

- Add a small experimental local search helper near the 3D strategy layer, not as a public workflow contract.
- Starting from a top discrete projection direction, generate a tangent basis and evaluate a small stencil:
  - center
  - positive/negative tangent axis 1
  - positive/negative tangent axis 2
  - optionally diagonals
- For each perturbed direction:
  - extract a temporary polar central section from the current reference volume
  - load it into a temporary PFTC reference slot or local scratch route
  - call existing `pftc%gen_objfun_vals`
  - keep the best in-plane rotation and score
- Use this only on top candidates after the normal discrete search.

This draft intentionally avoids analytic derivatives. It tests the landscape and integration points before adding interpolation-gradient code.

Expected result:

- Improvement should be strongest when the projection grid is coarser than the effective angular precision of the data.
- The best perturbation should usually be small relative to the local grid spacing.
- If the objective is noisy or discontinuous, do not proceed to full derivatives until the source is understood.

### Draft 1: PFTC Gradient With Finite-Difference Section Derivatives

Purpose: test the PFTC objective-gradient algebra without yet differentiating Kaiser-Bessel interpolation.

Implementation:

- Build `P0`, `P_plus_a`, `P_minus_a`, `P_plus_b`, `P_minus_b`.
- Approximate:

```text
dP/da = (P_plus_a - P_minus_a) / (2*eps)
dP/db = (P_plus_b - P_minus_b) / (2*eps)
```

- Add a private PFTC routine that evaluates:
  - score at fixed `irot`
  - gradient with respect to the two supplied derivative sections
- Keep the section extraction finite-difference based.

This is the best place to verify the normalized CC and Euclidean derivative expressions against finite differences of the final objective.

Expected result:

- Directional derivatives from the new PFTC routine should match finite differences of `gen_corr_for_rot_8` or `gen_objfun_vals` for fixed `irot`.
- The derivative should become less stable near in-plane rotation ties; this is expected and should be handled by reselecting `irot` between optimizer steps.

### Draft 2: Analytic Kaiser-Bessel Interpolation Gradient

Purpose: remove the most expensive finite-difference section extraction.

Implementation:

- Extend `simple_projector` with a private or experimental interpolation-gradient API:

```fortran
pure subroutine interp_fcomp_oversamp_grad(self, loc, comp, grad_comp)
    class(projector), intent(in) :: self
    real,             intent(in) :: loc(3)
    complex,          intent(out) :: comp
    complex,          intent(out) :: grad_comp(3)
end subroutine
```

- Implement the derivative of the separable KB interpolation:

```text
V(loc) = sum_i,j,k wx_i(x) * wy_j(y) * wz_k(z) * C(i,j,k)

dV/dx = sum_i,j,k dwx_i/dx * wy_j * wz_k * C(i,j,k)
dV/dy = sum_i,j,k wx_i * dwy_j/dy * wz_k * C(i,j,k)
dV/dz = sum_i,j,k wx_i * wy_j * dwz_k/dz * C(i,j,k)
```

- Add either an analytic derivative for `kbinterpol%apod` or a carefully bounded local finite-difference derivative of the 1D KB weight as an intermediate step.
- For oversampled interpolation, account for the padded coordinate scaling:
  - `locpd = loc * OSMPL_PAD_FAC`
  - derivatives with respect to original logical coordinates require the pad factor.
  - existing amplitude scaling `OSMPL_PAD_FAC**3` must also be preserved.

Expected result:

- `interp_fcomp_oversamp_grad` should match finite differences of `interp_fcomp_oversamp` for random locations inside the valid Fourier support.
- The analytic extraction of `dP/da` and `dP/db` should match finite-difference section derivatives.

### Draft 3: Optimizer Object for Projection Direction

Purpose: integrate local continuous refinement in the same style as shift search.

Implementation:

- Add a new optimizer object analogous to `pftc_shsrch_grad`, for example:

```text
simple_pftc_projdir_srch_grad
```

- Keep ownership in the strategy/search layer:
  - the builder still owns reference volumes/PFTC context
  - the search object owns optimizer state
  - the strategy decides when to call the continuous refinement
- The optimizer variables are two tangent-plane coordinates.
- The optimizer callback:
  - generates the perturbed orientation
  - extracts `P`, `dP/da`, `dP/db`
  - evaluates all in-plane rotations or evaluates the current `irot`
  - updates the best in-plane rotation through a callback, as `pftc_shsrch_grad` already does
  - returns score and gradient

Expected result:

- Continuous refinement should be available behind an experimental parameter or compile-time/local branch first.
- It should operate only on shortlisted peaks, never inside the full reference loop in the first implementation.

## Tests and Validation

### Unit-Level Numerical Tests

1. **KB interpolation gradient test**
   - Location: near existing numerical tests, e.g. `src/main/commanders/test/simple_commanders_test_numerics.f90`.
   - Compare `interp_fcomp_oversamp_grad` against finite differences of `interp_fcomp_oversamp`.
   - Test random interior points and points near typical polar-section radii.
   - Avoid window-boundary discontinuities in the first test; add boundary stress tests later.

2. **Section derivative test**
   - Build a small synthetic volume.
   - Extract `P(theta)`.
   - Compare analytic `dP/da`, `dP/db` against finite-difference section extraction.
   - Check both low and high Fourier shells separately.

3. **PFTC objective derivative test at fixed in-plane rotation**
   - Use a stored particle PFT and a synthetic reference section.
   - Compare analytic direction derivatives of the objective against finite differences.
   - Test both `OBJFUN_CC` and `OBJFUN_EUCLID` if both are expected to be supported.

4. **In-plane argmax stability test**
   - Construct cases where the best in-plane rotation is unambiguous.
   - Construct cases near a tie and verify the optimizer does not assume smoothness across the tie.
   - Confirm callback reselection of `irot` works as intended.

### Integration Tests

1. **Single-particle known-orientation recovery**
   - Project a known volume at a known non-grid projection direction.
   - Start from the nearest discrete direction.
   - Verify local continuous refinement moves toward the true direction.
   - Measure angular error before and after refinement.

2. **Grid-coarsening experiment**
   - Run a small controlled refinement with the standard projection grid.
   - Run a coarser grid plus continuous local refinement.
   - Compare:
     - final correlations
     - angular errors against known truth if synthetic
     - runtime
     - number of full section extractions

3. **No-regression discrete path**
   - With continuous projection refinement disabled, results must match the current discrete path.
   - This protects the normal `refine3D` and `abinitio3D` workflows.

4. **Symmetry and mirror checks**
   - Verify behavior with `fproject_polar_batch_mirr`.
   - Confirm local tangent updates respect the projection-direction equivalences used by the current symmetry handling.

### Performance Measurements

Measure separately:

- cost of extracting one central section
- cost of extracting section plus two derivative sections
- cost of all-rotation PFTC objective evaluation
- cost per optimizer step
- total added cost per particle for top-1, top-3, and top-5 candidate refinement

The first useful performance target is not full replacement of discrete search. The first target is:

```text
coarser discrete grid + local continuous refinement
```

matching or improving angular accuracy at lower total runtime.

## Design Risks

1. **Objective nonsmoothness from discrete in-plane argmax**
   - The score is smooth for fixed `irot`, but the max over `irot` is piecewise smooth.
   - Mitigation: optimize fixed `irot` within each step, periodically reselect `irot`, and accept only actual score improvements.

2. **Interpolation window boundary effects**
   - KB interpolation with `nint(loc)` changes the active interpolation window at half-grid boundaries.
   - Mitigation: tests should first avoid boundaries, then characterize boundary behavior. Optimizers may still work because these events are sparse, but gradients will not be perfectly smooth there.

3. **Euler-angle singularities**
   - Direct Euler derivatives are fragile.
   - Mitigation: use a tangent-plane local parameterization on the sphere.

4. **Symmetry bookkeeping**
   - Continuous perturbations must preserve the same projection-direction/state semantics as the discrete search.
   - Mitigation: keep the first implementation local to post-discrete peak refinement and map accepted orientations back through existing `ori` assignment logic.

5. **Memory and scratch pressure**
   - Derivative sections add scratch arrays of roughly `2 * pftsz * nk` complex values per active candidate/thread.
   - Mitigation: allocate per-thread scratch once and reuse it, following the PFTC heap workspace style.

## Recommended First Milestone

The first milestone should be deliberately modest:

1. Implement a derivative-free tangent stencil for the best discrete projection candidate only.
2. Run a synthetic known-orientation test.
3. Confirm that local continuous refinement improves angular error and/or correlation.
4. Only then add derivative-section plumbing.

This gives an early go/no-go signal with minimal impact on the existing PFTC code. If the stencil improves known-orientation recovery, the derivative path is worth pursuing. If it does not, the issue is likely objective landscape, reference quality, symmetry ambiguity, or local parameterization rather than a missing analytic derivative.

## Long-Term Direction

If the derivative path proves stable, the full search could evolve toward:

- a coarser global projection grid
- local tangent refinement of top peaks
- existing discrete in-plane rotational search
- existing continuous shift refinement

That would preserve SIMPLE's PFTC strengths while adding continuous exploration only where it is most likely to pay off.
