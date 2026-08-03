# Continuous In-Plane Rotation Refinement

## Summary

Implement continuous in-plane angle refinement for `objfun=euclid` in gated stages:

1. Establish a parabolic sub-grid baseline.
2. Add exact continuous angular interpolation from the Euclidean residual’s angular Fourier coefficients.
3. Validate the method and produce the required four-way recovery table.
4. Only if justified, extend L-BFGS-B to jointly refine `(sx, sy, theta)`.
5. Propagate the continuous angle through alignment metadata and class-average assembly.

The feature is opt-in: `inpl_refine=no` is the default, and continuous refinement is
enabled only by `inpl_refine=yes` in the eligible classical Euclidean stages.

The existing integer `irot` interface remains compatible throughout the transition.

## Progress

- [x] Phase 0: baseline, Nyquist behavior, and test harness.
  - Confirmed the current Euclidean residual accumulation uses `p=1:pftsz` and omits the `pftsz+1` Nyquist bin.
  - Added the focused `simple_test_euclid_route_identity` executable; direct SIMPLE invocation with explicit arguments is the authoritative runtime check.
  - Confirmed the legacy Euclidean score, raw FFT loss, and scalar gradient route are checked for every discrete rotation by the harness.
  - Compile/link verification passed with the Windows MSYS2 UCRT64 workflow.
  - Runtime validation passed on Oracle Linux 8.10 with 288 rotations: legacy-score error `2.97e-08`, scalar-gradient error `1.82e-07`, tolerance `8.50e-04`.
- [x] Phase 1: parabolic sub-grid interpolation.
  - Implemented a Euclidean-only three-point parabolic peak offset with periodic neighbors and a safeguarded fallback.
  - Preserved discrete selection for non-Euclidean objectives and retained the integer `cur_inpl_idx` compatibility output.
  - Added the continuous grid-index angle output used for shift-frame rotation and optional downstream propagation.
  - Added a focused synthetic parabola check to `simple_test_euclid_route_identity`.
  - Windows UCRT64 compile/link verification passed.
  - Oracle Linux 8.10 runtime validation passed with 288 rotations: legacy-score error `2.91e-08`, scalar-gradient error `3.60e-07`, tolerance `1.11e-03`; the focused parabolic check also passed.
  - Revalidated the updated `HEAD` on Oracle Linux 8.10 with two independent runs: legacy-score errors `2.87e-08` and `2.91e-08`, scalar-gradient errors `2.70e-07` and `3.70e-07`; both completed with `NORMAL STOP`.
- [x] Phase 2 implementation: continuous angular refinement (classical Euclidean likelihood path only).
  - The current implementation changes only the classical polar Euclidean refinement path and its focused route-identity test.
  - No streaming-SGD implementation or test path is being extended; treat SGD as out of scope for this effort.
  - Windows MSYS2 UCRT64 compile/link verification completed before the build was stopped.
- [x] Phase 2 focused validation: coefficient route identity and angular derivatives.
  - Oracle Linux 8.10 runtime passed with 288 rotations and `NORMAL STOP`.
  - Legacy-score error `2.96171144e-08`; scalar-gradient error `2.26158719e-07`; tolerance `1.08889884e-03`.
  - Angular grid error `1.33313786e-07`; periodic error `9.36750677e-17`; first-derivative error `1.32675199e-08`; second-derivative error `9.03031514e-09`; angular tolerance `1.08889884e-03`.
- [x] Phase 2 extended validation: aliasing experiment, synthetic recovery, and the required four-way recovery table.
  - Added `simple_test_euclid_stage1_validation`, an Oracle-runnable production test for the two-band aliasing comparison and deterministic synthetic recovery report.
  - The report includes grid-only, parabolic, and continuous rows; the Phase 3 joint row is explicitly reported as `NOT_IMPLEMENTED` until Phase 3 exists.
  - The harness now evaluates zero-shift, nonzero-shift, and near-periodic-boundary truths, and repeats the low-pass/full-band comparison with a hard-edged near-Nyquist fixture.
  - All three recovery rows use the same fixed-angle classical Euclidean direct-shift minimizer. The harness reports expected and recovered shift vectors plus acceptance flags.
  - For the polar `shift_ptcl` fixture, the expected candidate is `R(truth_angle) * applied_shift`; this is distinct from the `-R(angle) * applied_shift` corrective convention used by the image-space `rtsq` fixture in `simple_test_sgd_base_suite`.
  - Oracle Linux 8.10 validation passed with `NORMAL STOP`.
  - Three truths were covered: zero shift/zero angle, nonzero shift at `37°`, and nonzero shift near the periodic boundary at `359.375°`.
  - All three stages accepted the common fixed-angle shift refinement. Shift RMS over cases was `0.00704` pixels; the worst individual case was `0.01039` pixels.
  - Continuous angle RMS was `0.3041°`, improving over parabolic `0.3139°` and grid-only `0.4621°`; continuous was better than parabolic for each case.
  - The harder hard-mask near-Nyquist comparison produced low-pass RMS `1.8403e-04°` and full-band RMS `4.1441e-05°`, with no observed aliasing degradation.
  - The Stage 1 report now includes the Phase 3 joint row and four-way recovery summary.
- [x] Phase 3: joint `(sx, sy, theta)` refinement.
  - Added a classical Euclidean continuous-angle gradient API using the angular Fourier residual coefficients and differentiated shift cross terms.
  - Added a three-variable L-BFGS-B joint refinement entry point with a periodic angle window and monotonic acceptance relative to the Stage 1 starting point.
  - Extended the Stage 1 Oracle test to report the Phase 3 row and central-difference gradient error for all three truths.
  - Oracle Linux 8.10 targeted build and runtime validation passed with `NORMAL STOP`.
  - All three joint solutions were accepted; the gradient checks were finite with maximum errors of `9.20e-03`, `2.48e-03`, and `9.40e-03` for the three fixtures.
  - The joint angle RMS was `0.2756°`, improving over continuous-angle RMS `0.3041°`; Windows recompilation was intentionally not repeated on the laptop.
- [x] Phase 4A: classical Euclidean 2D metadata and workflow wiring.
  - Added a gated production gateway in the classical Euclidean 2D search strategy: the discrete candidate is refined by Stage 1, then passed to the Phase 3 three-variable L-BFGS-B path.
  - The public `inpl_refine` key defaults to `no`; the gateway requires `inpl_refine=yes` and is disabled for streaming SGD, time-series shift-only search, probabilistic-table generation, and non-Euclidean objectives.
  - Selected orientations receive the continuous `e3` value while legacy integer `inpl` indices remain populated for search tables and compatibility.
  - A real Oracle Linux classical Euclidean `abinitio2D` workflow completed five iterations with `objfun=euclid` and `sgd=no`, including repeated `SIMPLE_CLUSTER2D NORMAL STOP` and final `SIMPLE_ABINITIO2D NORMAL STOP` markers.
  - Focused route-identity, Stage 1, and SGD regression tests also passed on Oracle Linux; conventional CTest is not used as the SIMPLE acceptance gate because the project test workflow is direct `simple_test_*` execution.
- [ ] Phase 4B: classical Euclidean 3D metadata and workflow wiring.
  - The 3D gateway implementation is present but intentionally not accepted or committed until the 2D pathway is fully reviewed.
  - 3D workflow validation and any required 3D corrections remain pending.
- [ ] 2D-only completion: production metadata and final regression.
  - Added `simple_test_euclid_2d_metadata` to verify persisted continuous `e3`, compatible integer `inpl`, and class-average metadata after a real `abinitio2D` run.
  - The check reads the final project metadata; it does not depend on temporary `algndoc_*.simple` files, which the normal workflow removes during cleanup.
  - When the workflow is run with `mkdir=yes`, launch from the repository `build/` directory and validate the copied project inside a build-local execution directory such as `build/1_abinitio2D`; do not place generated artifacts at the repository root.
  - Before the opt-in switch was added, Oracle compilation and execution passed: `ACTIVE=200`, `OFFGRID=199`, `INVALID=0`, `NONFINITE=0`, `ROTATIONS=288`, and `CLASS_AVERAGES=3`. The default-off and explicit opt-in contracts now require rerunning this regression.
  - Pending Oracle gates: classical checkpoint stop/resume, public `abinitio2D` execution with the non-Euclidean `cc` objective (the internal cluster2D path is not a public `simple_exec prg`), and a fresh real-data Euclidean run with final-project metadata inspection.
  - Checkpoint creation and resume now pass in `build/1_abinitio2D`: stage 1 reached `last_iter=5`, persisted the required reference stacks, stage 2 resumed from iteration 5, completed through iteration 10, generated final class averages/rankings, and ended with `SIMPLE_ABINITIO2D NORMAL STOP`. The resume must run from the build-local execution directory because `mkdir=no` preserves relative reference paths.
  - Restart metadata validation also passed: `ACTIVE=200`, `OFFGRID=200`, `INVALID=0`, `NONFINITE=0`, `ROTATIONS=288`, and `CLASS_AVERAGES=3`.
  - The public `abinitio2D` `objfun=cc` alternate path completed clustering and class-average generation but exposed an existing finalization bug: sigma2 metadata was registered unconditionally although `cc` correctly produces no `sigma2_it_*.star` file. Finalization now skips this optional project entry when the file is absent; Oracle rebuild and rerun remain pending.
  - The metadata test now checks both contracts: `inpl_refine=yes` must produce at least one off-grid continuous `e3`, while `inpl_refine=no` must produce none. The `cc` workflow is tested with the key enabled to confirm the objective gate still keeps continuous refinement off.
  - The first post-push default-off Oracle metadata run exposed stale/off-grid `e3` propagation despite the integer `inpl` result. The orientation writer now reconstructs the grid `e3` from `inpl` unless the explicit continuous gateway marks a valid angle; Oracle rebuild and rerun remain required.
  - After these gates, restore and validate the deferred 3D extension.

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

- Scope this refinement to the classical SIMPLE Euclidean likelihood path. Do not add new SGD integration while that path is being reconsidered.

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
   - Run the existing polar FFT, class-search, 2D search, and 3D search tests affected by the API.
   - Confirm non-Euclidean objectives and direct-only paths remain unchanged.
   - Run public `abinitio2D` with `objfun=euclid inpl_refine=no`, with `objfun=euclid inpl_refine=yes`, and with `objfun=cc inpl_refine=yes`; only the explicit Euclidean opt-in run may persist off-grid continuous `e3` values.

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
- Continuous refinement is opt-in through `inpl_refine=yes`; the default is `inpl_refine=no`.
- The runtime gateway additionally requires a classical, non-SGD, non-time-series, non-probabilistic gated stage.
- The current dirty worktree changes belong to the user and must not be overwritten.
- Integer `irot` remains the compatibility index; continuous angle is an additive output and metadata enhancement.
- The original plan’s sign convention and local ±2-grid-step Stage 2 window are authoritative.
