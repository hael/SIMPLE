# Continuous 3-D pose polishing: Hans clarification

**Status:** AUTHORIZED EXPLORATORY DIAGNOSTIC - NO INTEGRATION CONTRACT
**Date:** 2026-08-19
**Source:** Mazharul's report of a discussion with Hans Elmlund

**Original proposal:**
[continuous_3D_refinement_on_pcg_operator.md](continuous_3D_refinement_on_pcg_operator.md)
**Decommissioned SPEC:**
[continuous_3D_pose_end_polishing_spec.md](continuous_3D_pose_end_polishing_spec.md)
**Decommissioned PLAN:**
[continuous_3D_pose_end_polishing_plan.md](continuous_3D_pose_end_polishing_plan.md)

## Reported requirement

1. The future continuous five-parameter pose polisher is intended to replace
   `inpl_cont` in `refine3D`.
2. Reuse the established `inpl_cont` workflow logic and position, not its 2-D
   PFTC numerical objective.
3. Evaluate each trial pose from a 3-D reference volume through the Cartesian
   central-section forward model developed with the PCG operator.
4. Optimize all five local pose parameters together:

   $$
   q=(\omega_x,\omega_y,\omega_z,\delta t_x,\delta t_y).
   $$

   The three tangent coordinates include the in-plane rotation and both
   out-of-plane rotations. The remaining two coordinates are image shifts.
5. PCG does not optimize these five parameters. It supplies the relevant
   Cartesian volume/forward-operator design; the local pose update uses the
   five-parameter Jacobian and Levenberg--Marquardt.
6. A later `refine3D` integration will need to retain or rematerialize the
   prepared 3-D Cartesian reference efficiently. The current matching path can
   reduce the reference to stored polar Fourier sections for its workers, so
   the required 3-D workspace lifetime must be designed later.
7. Do not implement `refine3D`, UI, parameter, persistence, reconstruction, or
   distributed integration now.

## Immediate experiment requested by Hans

Implement one focused numerical test only:

1. Create or load one known 3-D volume.
2. Generate a clean Cartesian Fourier plane at a known rotation and shift with
   the same PCG `forward_plane` model used by the evaluator.
3. Perturb the starting rotation and shift by controlled known amounts.
4. Run the existing five-parameter Jacobian and bounded LM solve.
5. Sweep increasingly difficult perturbations. The working interpretation is
   to test angular errors through 15 degrees and image-shift errors through at
   least 3 pixels. The exact scenario table is an implementation choice unless
   Hans specifies one.
6. Record where recovery succeeds and where it fails. For every trial, report:
   - injected rotation and shift error;
   - final rotation and shift error;
   - initial and final objective;
   - attempted and accepted LM iterations;
   - terminal status and step-bound or stencil-switch evidence.

This is a capture-range experiment. It should characterize how far the local
five-parameter model can move back to the known pose. It must not add a
production activation.

The sweep is exploratory. Hans expects unexpected behavior because this method
has not been implemented or tested before. There is no predetermined acceptable
final-error threshold for the larger perturbations. The working scientific
hypothesis is that recovery error remains similar through approximately 15
degrees. Report each isolated result before changing another parameter family;
do not hide a failed recovery behind an aggregate pass/fail result.

## Evidence and visualization artifacts

Store all results from one run in one timestamped directory. Include:

- the fixed input 3-D volume;
- the observation at the known pose;
- the prediction at each perturbed starting pose;
- the prediction at each final recovered pose;
- initial and final residual or difference images;
- a machine-readable table and a Markdown summary of all pose and LM metrics;
- the configured Fourier shell range and its nominal low-pass resolution.

The local test recovers a pose, not a 3-D volume. A single 2-D observation
cannot determine a new 3-D reconstruction. Therefore the useful visual
"recovery" artifacts are the final reprojected plane and its difference from
the known observation. Reconstructing another 3-D volume would require many
projections and is outside this first experiment.

## Existing test relationship

`case=pose_recovery` in
`production/tests/simple_continuous_3D_pcg_refinement_recovery_test.f90`
already exercises one clean, matched-operator joint rotation-and-shift recovery
point. It does not yet provide the requested perturbation sweep or a capture-
range table. Its comment calls the observation independently generated, but
the source currently generates it through the same PCG Fourier workspace; that
description must not be used as evidence of an independent forward generator.

## What the experiment can establish

- Correct local behavior of the executed Cartesian Jacobian and LM driver.
- Exact-pose stationarity under the matched numerical model.
- Empirical rotation/shift capture range for the selected volume, shell range,
  bounds, and LM policy.
- Whether failure is caused by iteration limits, bounds, unreliable curvature,
  or convergence to another local solution.

It cannot establish production `refine3D` behavior, FSC improvement, robustness
to noise or CTF variation, reference-volume ownership, or a safe replacement of
`inpl_cont`.

## Initial one-at-a-time sweep

Use one focused test harness, but change only one parameter family at a time:

1. **Exact-pose control:** use the true rotation and shift. This must remain
   stationary.
2. **Shift-only:** perturb one image-shift coordinate while keeping the truth
   rotation and the other shift coordinate fixed. Repeat for each shift
   coordinate and magnitude.
3. **Rotation-only:** perturb one tangent-space rotation coordinate while
   keeping the truth shift and the other rotation coordinates fixed. Repeat for
   each rotation coordinate and angular magnitude through approximately 15
   degrees.

Record every result, including unexpected convergence, divergence, bounds,
stencil switches, and alternate local solutions. Do not add joint
rotation-and-shift perturbations to the initial experiment. Add them only after
the isolated cases behave well enough to justify testing the coupled 5-by-5 LM
system.

The same Cartesian plane dimensions and Fourier sampling are used in every
scenario. Values such as 1 or 3 pixels describe the starting shift error; they
do not describe different Cartesian plane sizes.

The first experiment is currently understood to be clean and matched to the
PCG Cartesian operator. CTF, sigma weighting, noise, independent projection,
and production frequency policy are later experiments unless Hans says they
belong in this initial capture-range measurement.

The existing PFTC-based SPEC and PLAN are decommissioned and must not authorize
integration. This note is the bounded record for the explicitly authorized
exploratory diagnostic. Any later `refine3D` integration requires a new SPEC
and PLAN after the results are reviewed with Hans.
