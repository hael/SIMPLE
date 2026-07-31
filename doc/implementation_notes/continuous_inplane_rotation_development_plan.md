# Continuous In-Plane Rotation Refinement

## Summary

Implement continuous in-plane angle refinement for `objfun=euclid` in gated stages:

1. Establish a parabolic sub-grid baseline.
2. Add exact continuous angular interpolation from the Euclidean residual’s angular Fourier coefficients.
3. Validate the method and produce the required four-way recovery table.
4. Only if justified, extend L-BFGS-B to jointly refine `(sx, sy, theta)`.
5. Propagate the continuous angle through alignment metadata and class-average assembly.

The existing integer `irot` interface remains compatible throughout the transition.

## Progress

- [x] Phase 0: baseline, Nyquist behavior, and test harness.
  - Confirmed the current Euclidean residual accumulation uses `p=1:pftsz` and omits the `pftsz+1` Nyquist bin.
  - Added the focused `simple_test_euclid_route_identity` CTest target.
  - Confirmed the legacy Euclidean score, raw FFT loss, and scalar gradient route are checked for every discrete rotation by the harness.
  - Compile/link verification passed with the Windows MSYS2 UCRT64 workflow.
  - Runtime validation passed on Oracle Linux 8.10 with 288 rotations: legacy-score error `2.97e-08`, scalar-gradient error `1.82e-07`, tolerance `8.50e-04`.
- [ ] Phase 1: parabolic sub-grid interpolation.
  - Implemented a Euclidean-only three-point parabolic peak offset with periodic neighbors and a safeguarded fallback.
  - Preserved discrete selection for non-Euclidean objectives and retained the integer `cur_inpl_idx` compatibility output.
  - Added the continuous grid-index angle output used for shift-frame rotation and optional downstream propagation.
  - Added a focused synthetic parabola check to `simple_test_euclid_route_identity`.
  - Windows UCRT64 compile/link verification passed; Oracle Linux runtime validation is pending.
- [ ] Phase 2: continuous angular refinement.
- [ ] Stage 1 validation: route identity, gradients, aliasing, and synthetic recovery.
- [ ] Phase 3: joint `(sx, sy, theta)` refinement.
- [ ] Phase 4: downstream metadata and workflow wiring.

## Implementation Changes

### 1. Baseline and numerical contract

- Preserve all current unrelated worktree changes in `simple_pftc_shsrch_grad.f90`, `simple_polarft_corr.f90`, and the current strategy and SGD test files.
- Add Stage 0 parabolic interpolation around the discrete Euclidean residual minimum.
- Use grid-index units internally:
  - integer `j` maps to `angtab(j)`;
  - continuous `theta` uses the same index coordinate;
  - physical angle is `(theta - 1) * get_dang()`.
- Preserve the current rotation convention `REF(phi - theta)` and the corresponding negative rotation derivative.
- Resolve the `pftsz` versus `pftsz+1` Nyquist behavior before interpolation:
  - add a route-identity diagnostic;
  - preserve the current `gen_euclids` bin inclusion semantics;
  - do not silently change the existing grid objective.

### 2. Continuous angular residual API

Extend `simple_polarft_calc.f90` and `simple_polarft_corr.f90` with Euclidean-only operations:

- `gen_euclid_angular_coeffs`
  - compute the same weighted residual spectrum currently produced by `gen_euclids`;
  - retain and return the angular Fourier coefficients;
  - return the best discrete index for initialization.
- `eval_euclid_resid_at_angle`
  - evaluate `r(theta)`, `r'(theta)`, and `r''(theta)` without another FFT;
  - apply interior Fourier-bin weight `2`;
  - apply Nyquist weight `1` when the Nyquist bin is present.
- Add a continuous-angle Euclidean objective/gradient entry point for Stage 2, returning the raw residual and gradients with respect to `(sx, sy, theta)`.
- Keep the coefficient route numerically aligned with the existing FFT precision and normalization.

### 3. Stage 1 angle refinement

Update `simple_pftc_shsrch_grad.f90`:

- Add `cur_inpl_ang` alongside `cur_inpl_idx`.
- Replace the discrete angle callback’s final `maxloc` selection with:
  - discrete initialization from the coefficient minimum;
  - two or three safeguarded Newton iterations on `r'(theta)=0`;
  - clamp movement to ±0.5 grid steps;
  - fall back to the discrete index when `r''(theta) <= 0` or the candidate is non-finite.
- Continue setting `cur_inpl_idx = modulo(nint(cur_inpl_ang)-1,nrots)+1` so existing integer consumers remain valid.
- Add an optional real `theta` output to `grad_shsrch_minimize`.
  - On success, return the refined continuous angle.
  - On failure, retain the existing `irot=0` failure signal and return `theta=0`.
- Use the continuous angle for the shift rotate-back whenever the optional output is available.
- Keep `cc`, hybrid, denoised, direct-only, and streaming paths on their current discrete-angle behavior.

### 4. Stage 2 joint refinement

After the Stage 0–1 validation table is complete:

- Extend the non-direct Euclidean optimizer from two to three variables:
  - `vec = [sx, sy, theta]`;
  - `opt_spec` dimension `3`;
  - third limit bounded to `theta_stage1 ± 2` grid steps.
- Remove the angle callback from this joint path.
- Keep `direct_only` strictly two-dimensional.
- Implement the direct-sum continuous Euclidean gradient first because it reuses the existing residual-gradient structure:
  - memoize the band-limited angular reference derivative;
  - compute the closed-form derivative of the shift phase from the perpendicular polar coordinates;
  - add the rotation component to the existing residual loop.
- Assert the square-box assumption required by the perpendicular phase derivative.
- Preserve periodic angle evaluation even when the local optimizer window crosses the first or last grid index.
- Enforce monotonicity: the Stage 2 accepted residual must not exceed the Stage 1 residual beyond optimizer tolerance.

### 5. Downstream wiring

Only after the numerical stages pass independently:

- Pass the optional continuous angle through 2D and 3D search strategies.
- Replace final `get_rot(irot)` conversions with the returned continuous angle where available.
- Preserve integer indices for search tables, mirroring, and legacy APIs.
- Store continuous degrees in alignment metadata and orientation records.
- Verify class-average assembly and restoration use the continuous angle consistently.
- Do not modify out-of-plane rotation or the separate 3D Cartesian interpolation design.
- Keep `cc` and hybrid/denoised support as later follow-on work.

## Test Plan

Add focused tests under the existing production test framework, preferably alongside the polar PFTC tests.

1. **Route identity**
   - Compare coefficient-route Euclidean values against `gen_euclid_grad_for_rot_8` for every grid rotation.
   - Require agreement to single-precision round-off.
   - Explicitly cover Nyquist weighting, conjugation, index wrapping, and normalization.

2. **Gradient verification**
   - Compare analytic `dtheta`, `dsx`, and `dsy` against central differences at multiple off-grid angles and shifts.
   - Verify second derivatives used by Newton refinement.
   - Test angles near the periodic boundary.

3. **Stage 0 and Stage 1 regression**
   - Compare grid-only, parabolic, and continuous-angle minima on synthetic data.
   - Confirm Stage 1 never returns a worse residual than Stage 0.

4. **Aliasing experiment**
   - Generate a known real-space rotation.
   - Run once with a well-low-passed reference and once with a hard-edged/full-band reference.
   - Report the angular error numerically; do not reduce this to pass/fail.

5. **Synthetic recovery**
   - Recover known continuous `(theta, sx, sy)` under realistic CTF and per-particle noise.
   - Produce the required table:
     - grid-only;
     - Stage 0 parabola;
     - Stage 1 continuous angle;
     - Stage 2 joint refinement.
   - Report RMS error separately for angle and both shifts.

6. **Monotonicity**
   - Assert Stage 2 residual ≤ Stage 1 residual for every particle within optimizer tolerance.

7. **Zero-width compatibility**
   - Clamp the Stage 2 angular window to zero.
   - Require results to match the existing discrete-angle behavior.

8. **Workflow regression**
   - Run the existing polar FFT, class-search, 2D search, 3D search, and SGD tests affected by the API.
   - Confirm non-Euclidean objectives and direct-only paths remain unchanged.

## Acceptance and Rollout

- Stage 0–1 is the first mergeable deliverable.
- Stage 2 cannot begin until route identity, gradient, aliasing, and synthetic recovery tests pass and the four-column recovery table is recorded.
- Stage 3 cannot begin until Stage 2 passes monotonicity and zero-width regression.
- Profile FFT counts per outer iteration and optimizer convergence; do not judge performance using isolated inner-loop flop counts.
- On real data, evaluate ring-wise or radially resolved half-set correlation and cluster iteration count; global FRC alone is insufficient.
- Keep `cc`, hybrid, and denoised objectives explicitly out of the first implementation.

## Assumptions

- The full staged roadmap is desired, with Stage 0–1 as the gated first deliverable.
- Continuous refinement is enabled only for `objfun=euclid` and the non-direct refinement path.
- The current dirty worktree changes belong to the user and must not be overwritten.
- Integer `irot` remains the compatibility index; continuous angle is an additive output and metadata enhancement.
- The original plan’s sign convention and local ±2-grid-step Stage 2 window are authoritative.
