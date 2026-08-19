# Continuous 3D pose-capture experiment summary

**Date:** 2026-08-19
**Status:** Isolated numerical experiment complete; production integration is not authorized
**Primary evidence:** `continuous_3D_pose_capture_20260819_145957`

## Executive result

The five-parameter analytic Jacobian and Levenberg--Marquardt (LM) implementation
pass direct finite-difference tests. The clean matched-operator experiments also
show that LM can recover many known rotation and shift perturbations with high
accuracy. However, the capture basin depends strongly on the volume, the sign and
size of the displacement, and the optimization route.

The most important result is that a large decrease in the Cartesian-plane image
objective does not prove pose recovery. For one Gaussian fixture, direct joint LM
reduced the objective by approximately 99.7 percent but converged to approximately
43.1 degrees of rotation error. The same initial pose recovered when the shift was
optimized before the joint five-parameter solve. Two other deterministic volumes
did not reproduce the 43-degree endpoint under direct joint LM.

The successful shift-first route has two stages:

1. Keep rotation fixed and optimize only the two shifts $(t_x,t_y)$.
2. Starting from the corrected shift, optimize all five parameters together:
   $(\omega_x,\omega_y,\omega_z,t_x,t_y)$.

The two-stage path is:

$$
(R_0,t_0)
\xrightarrow{\text{shift only}}
(R_0,t_1)
\xrightarrow{\text{joint five-parameter LM}}
(R_1,t_2).
$$

This is evidence of a local, morphology-dependent rotation/shift coupling problem.
It is not evidence of floating-point overflow, a universal 43-degree attractor,
or an elementary sign error in the analytic derivatives.

## Decisions to discuss with Hans

`P0` decisions block the integration design. `P1` decisions must be fixed in the
new SPEC before `pose_cont` can be enabled. The recommendations are proposals, not
an integration contract.

### 1. [P0] Objective and final acceptance

**Decision:** The clarification assigns the five-parameter update to the PCG-derived
Cartesian central-section objective. Should final acceptance also use that objective,
or should `refine3D` rescore the proposal with its executed PFTC raw-Euclidean
objective before replacing the discrete pose?

**Options:** Cartesian objective for both update and acceptance; or Cartesian LM
proposal followed by PFTC acceptance. Using PFTC as the LM objective conflicts with
the current clarification.

**Evidence:** Cartesian LM reduced one objective by 99.7 percent but ended 43.1
degrees from truth. The experiment did not test whether PFTC would reject it.

**Recommendation:** Confirm Cartesian LM as the proposal generator and ask Hans to
choose the final acceptance objective.

### 2. [P0] Optimization sequence

**Decision:** Should production use direct five-parameter LM or shift-first followed
by joint LM?

**Options:** Direct joint LM; shift-first then joint LM; or an adaptive fallback.
Do not use rotation-first.

**Evidence:** Direct joint LM recovered 12 of 15 multi-volume boundary cases.
Shift-first recovered all 15. Rotation-first failed by as much as 116 degrees and
11 pixels.

**Recommendation:** Use shift-first followed by joint LM in the first prototype and
retain direct joint LM only as a comparison route.

### 3. [P1] Workflow position and 3-D reference ownership

**Decision:** Does `pose_cont` refine only the selected discrete winner or every
candidate, and how is the Cartesian reference retained for that work?

**Options:** Refine the final winner or every candidate; retain one immutable
state-and-half reference or rematerialize it per worker/batch. A separate pure
continuous mode can support simulated-data tests.

**Evidence:** Hans said to replace the logical `inpl_cont` action and noted that the
current matcher may discard the 3-D representation after making polar sections.

**Recommendation:** Refine the final selected winner against one immutable
state-and-half reference. Confirm this position and the pure continuous test mode
with Hans before designing storage and distributed ownership.

### 4. [P1] Movement bounds, acceptance, and rollback

**Decision:** How far may a proposal move from its discrete seed, and what happens
when any acceptance condition fails?

**Options:** A fixed hard guard, a soft prior, or a trust region derived from the
current angular/shift sampling. Rejection can either restore the exact discrete seed
or retain a partially accepted pose.

**Evidence:** A seed-centered 15-degree/5-pixel guard contained every endpoint but
did not ensure recovery and sometimes blocked a useful trajectory. Monotone
Cartesian improvement alone also accepted a wrong pose.

**Recommendation:** Derive bounds from the discrete sampling scale and require exact
rollback after final rejection. Do not treat 15 degrees/5 pixels as a production
range.

### 5. [P1] Frequency and interpolation policy

**Decision:** Which frequencies drive the shift-first and joint stages, and how
should LM respond to interpolation-stencil switches?

**Options:** Use the exact current `kfromto` throughout; use lower frequencies for
shift-first and current `kfromto` for joint LM; or create a separate hand-set range.
For stencil switches, record telemetry, reduce the step, or reject it.

**Evidence:** The experiment used fixed shells 2--12, not the production NU-driven
limit. Derivatives passed inside fixed cells, but successful and failed large moves
both crossed many cells.

**Recommendation:** Do not create an independent production frequency range. Ask
whether shift-first may use a lower subset of current `kfromto`; retain bounded steps
and switch telemetry without an unvalidated switch-count threshold.

### 6. [P1] Validation required before enablement

**Decision:** Which gates must pass before `pose_cont=yes` is considered usable?

**Options:** Synthetic recovery with an independent forward generator; noisy/CTF simulation;
shared/distributed production smoke tests; symmetry, persistence, and restart tests;
and a matched workflow comparison using FSC or other held-out evidence.

**Evidence:** The completed tests establish matched-operator derivative correctness
and capture behavior only. They do not establish production or real-data benefit.

**Recommendation:** Require all implementation-correctness gates above. Ask Hans
which dataset-level result is required for scientific acceptance.

## Scope

The experiment uses a clean observation generated by the same PCG Cartesian
Fourier-plane model used by the optimizer. It evaluates three rotation tangent
coordinates and two image shifts.

It does not test:

- production `refine3D` integration;
- the executed PFTC matcher objective;
- noisy or experimental particles;
- CTF variation across particles;
- an independent forward generator;
- point-group symmetry or global gauge handling;
- FSC or real-data improvement.

No production workflow calls the new diagnostic options. The optional LM arguments
preserve the previous behavior when omitted.

## Tests completed

### 1. Componentwise derivative validation

The `rotation_gradient` case compares each of the five analytic residual-Jacobian
columns with centered finite differences of the executed prediction. It also
compares each objective-gradient component independently. The checks are repeated
at two nonstationary poses and with CTF/noise transfer weighting.

Observed worst errors were:

- unweighted residual-Jacobian column: `7.273379e-4`;
- unweighted objective-gradient component: `4.656480e-4`;
- weighted residual-Jacobian column: `5.899270e-4`;
- weighted objective-gradient component: `1.076386e-4`.

All were below the respective `0.015` column and `0.03` component limits. The
focused case passed, and the complete ten-case mother suite passed 10/10 on Oracle
Linux.

This test rules out the basic derivative-sign and compensating-component errors
that a single directional derivative test could miss. It does not remove the
known mathematical boundary at interpolation-stencil switches.

### 2. Exact-pose stationarity

The exact pose successfully produced `finite_no_improvement`, with zero accepted steps, zero
rotation error, and zero shift error. This confirms that the matched observation is
stationary under the tested numerical objective.

### 3. Separate rotation and shift sweeps

The final 138-trial matrix changes one coordinate at a time and then exercises
joint perturbations.

- Separate rotations about all three axes, both signs, recovered through 15
  degrees.
- Separate shifts recovered reliably through approximately 5 to 6 pixels in this
  24-pixel fixture, depending on axis and sign.
- Several shifts at 6.5 pixels and above entered coupled rotation/shift minima.
- Failure was asymmetric. For example, some positive 5.5- and 6-pixel shifts
  recovered while the corresponding negative x shifts did not.

The separate sweep therefore does not define one universal radial capture limit.

### 4. Joint boundary, sign, and multi-axis matrix

The principal positive x-rotation/x-shift results were:

| Initial rotation | Initial shift | Result |
| ---: | ---: | :--- |
| 2 degrees | 0.5 pixel | Recovered |
| 5 degrees | 1 pixel | Recovered |
| 10 degrees | 3 pixels | Recovered |
| 10 degrees | 4 pixels | Recovered |
| 10 degrees | 5 pixels | Converged to about 43.1 degrees and 0.746 pixel error |
| 12.5 degrees | 4 pixels | Recovered |
| 12.5 degrees | 5 pixels | Converged to about 43.1 degrees and 0.746 pixel error |
| 15 degrees | 4 pixels | Recovered |
| 15 degrees | 5 pixels | Converged to about 43.1 degrees and 0.746 pixel error |

The three mixed-sign cases and three simultaneous multi-axis cases in this matrix
recovered. The negative 15-degree/positive 5-pixel case also recovered. The failure
is therefore not determined by perturbation magnitude alone.

### 5. Multi-volume route diagnostic

Five boundary seeds were tested on three deterministic asymmetric volumes with
four routes:

1. direct joint LM;
2. joint LM with a cumulative 15-degree/5-pixel guard around the input seed;
3. shift-only LM followed by joint LM;
4. rotation-only LM followed by joint LM.

The results were:

- **Direct joint LM:** Failed in three positive 5-pixel Gaussian cases, but
  recovered all five cases for `shifted_mixture` and `permuted_texture`.
- **Guarded joint:** Kept every endpoint within the stated seed-centered guard,
  but did not guarantee recovery. It sometimes trapped an incorrect solution or
  blocked a useful coupled trajectory near the guard boundary.
- **Shift then joint:** Recovered all 15 volume/seed combinations to small truth
  error. This was the most robust tested route.
- **Rotation then joint:** Was unsafe. It produced errors up to approximately 116
  degrees and 11 pixels.

The guard is therefore a containment mechanism, not an acceptance or recovery
criterion.

### 6. Objective-path diagnostic

For every mechanism case, the straight right-tangent/linear-shift path from the
seed to truth had no sampled uphill segment. The endpoint objective at truth was
approximately `1e-15` of the starting objective. The paths to the direct-joint
endpoints were also monotone.

The Gaussian failure is not caused by an objective barrier between the seed and
truth. The local LM direction enters another descending valley because the
rotation and shift columns are coupled. The incorrect endpoint retains about
`0.003` of the starting objective, while the true pose reaches numerical zero.

### 7. Default-path regression

The runner executes `case=pose_recovery` without the new optional trajectory,
parameter-mask, or cumulative-guard arguments. It passed. This checks that the LM
API extensions do not change the established default call behavior.

### 8. Evidence-package validation

The final package contains 138 capture rows, 60 mechanism route endpoints,
accepted-step five-vector trajectories, 1,230 objective-path samples, MRC volumes
and observations, logs, source snapshots, and SHA-256 manifests.

The Fortran capture, mechanism, and default-regression cases all returned zero.
The first analysis used tolerances below the observed floating-point reconstruction
error. The corrected analyzer uses a `1e-6`-degree truth-endpoint tolerance and a
combined absolute/relative endpoint-objective tolerance. Reanalysis reports:

```text
POSE_CAPTURE_ANALYSIS: PASS
Package integrity: PASS
Scientific recovery warnings: 29
```

The 29 warnings are scientific capture warnings. They are not missing evidence or
numerical-integrity failures.

## Files changed for this experiment

### Numerical API

- `src/main/volume/simple_reconstructor_pcg.f90:1437` extends `refine_pose_lm`
  with optional accepted-pose trajectories, an optional five-coordinate activity
  mask, and an optional seed-centered cumulative guard.
- `src/main/volume/simple_reconstructor_pcg.f90:1645` applies the active-coordinate
  mask to the scaled gradient and Hessian.

### Fortran tests

- `production/tests/simple_continuous_3D_pcg_refinement_rotation_test.f90:205`
  adds componentwise residual-Jacobian and objective-gradient finite differences.
- `production/tests/simple_continuous_3D_pcg_refinement_pose_capture_test.f90:58`
  implements the 138-trial capture matrix and writes visualization/evidence files.
- `production/tests/simple_continuous_3D_pcg_refinement_pose_mechanism_test.f90:34`
  implements the three-volume route and objective-path diagnostic.
- `production/tests/simple_test_continuous_3D_pcg_refinement.f90:212` registers the
  two isolated child cases.

### Evidence tools

- `production/tests/continuous_3D_pcg_pose_capture/run_oracle_pose_capture.sh`
  runs and packages the three Oracle cases.
- `production/tests/continuous_3D_pcg_pose_capture/analyze_pose_capture.py:381`
  validates the capture matrix, mechanism endpoints, trajectories, objective paths,
  and required MRC artifacts.
- `production/tests/continuous_3D_pcg_pose_capture/reanalyze_oracle_pose_capture.sh`
  reanalyzes complete evidence without rerunning Fortran.
- `production/tests/continuous_3D_pcg_pose_capture/README.md` defines the experiment,
  evidence layout, and interpretation limits.

### Scientific documentation

- `doc/implementation_notes/continuous_3D_pose_polishing_hans_clarification_2026-08-19.md`
  records the clarified isolated-test request and decommissions the earlier
  production-integration direction.
- This report records the completed evidence and discussion questions.

The modified SPEC, PLAN, handoff, generated code indexes, and `.codex/AGENTS.md`
are broader documentation or repository-maintenance changes. They should be
reviewed and staged deliberately rather than assumed to belong to the isolated
diagnostic commit.

## Development status

The immediate experiment requested in the Hans clarification is complete. There
is no remaining numerical test gate in that request.

Production development is still blocked by the scope of the clarification itself.
It authorizes the isolated experiment and explicitly says not to implement
`refine3D`, UI, parameter, persistence, reconstruction, or distributed integration.
It also requires a new SPEC and PLAN after Hans reviews the results. The next safe
action is therefore to discuss the decisions above, then write the new integration
contract. The decommissioned SPEC and PLAN cannot authorize implementation.

## Recommended next sequence

1. Review these six prioritized decisions with Hans.
2. Record his decisions in a new concise SPEC for `pose_cont`.
3. Mark that SPEC FINAL before integration work.
4. Write a focused PLAN that traces the existing `inpl_cont` caller path and the
   required Cartesian reference lifetime.
5. Implement one phase at a time with the feature disabled by default.
