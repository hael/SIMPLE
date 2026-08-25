# Continuous 3-D pose refinement: consolidated validation evidence

**Status:** COMPLETED EVIDENCE RECORD
**Consolidated:** 2026-08-25
**Primary isolated package:** `continuous_3D_pose_capture_20260819_145957`
**Primary historical production package:** `continuous_3D_pose_truth_diagnostic_20260818_114404`
**History and handoff:** [continuous_3D_pose_end_polishing_history_and_handoff.md](continuous_3D_pose_end_polishing_history_and_handoff.md)

## Evidence boundary

This note records observed results. It separates numerical integrity, matched-model capture behavior, historical production-path behavior, and untested scientific claims.

The isolated capture observations were generated with the same Cartesian forward model used by the optimizer. They are matched-operator tests and therefore do not independently validate the forward generator. The historical production truth matrix used reconstructed references and failed its scientific gate. Neither package establishes real-data benefit.

## Numerical derivative evidence

The `rotation_gradient` case compared each of the five analytic residual-Jacobian columns with centered finite differences of the executed prediction. It independently compared each objective-gradient component, repeated the checks at two nonstationary poses, and included CTF/noise transfer weighting.

Observed worst errors were:

| Check | Worst error | Gate |
| --- | ---: | ---: |
| Unweighted residual-Jacobian column | `7.273379e-4` | `0.015` |
| Unweighted objective-gradient component | `4.656480e-4` | `0.03` |
| Weighted residual-Jacobian column | `5.899270e-4` | `0.015` |
| Weighted objective-gradient component | `1.076386e-4` | `0.03` |

The focused case and ten-case mother suite passed on Oracle Linux. These results rule out the earlier failure mode in which only finiteness or one combined directional derivative was tested. They do not validate derivatives at interpolation-stencil switches, where the executed piecewise interpolation model is not differentiable.

## Exact-pose stationarity

For the clean matched Cartesian observation, the exact pose returned `finite_no_improvement` with zero accepted steps, zero rotation error, and zero shift error. This verifies stationarity under the same model used to generate the observation.

It does not establish stationarity against an independently generated observation or a reconstructed reference.

## Separate perturbation sweeps

- Rotations about all three axes, with both signs, recovered through 15 degrees in the tested fixture.
- Separate shifts recovered reliably through approximately 5--6 pixels in the 24-pixel fixture, depending on axis and sign.
- Several shifts at 6.5 pixels and above entered coupled rotation/shift minima.
- Failures were asymmetric; magnitude alone did not determine recovery.

These results characterize one fixture and one numerical policy. They do not define a universal production capture radius.

## Joint pose boundary

The principal positive x-rotation/x-shift results were:

| Initial rotation | Initial shift | Final classification |
| ---: | ---: | --- |
| 2 degrees | 0.5 pixel | Recovered |
| 5 degrees | 1 pixel | Recovered |
| 10 degrees | 3 pixels | Recovered |
| 10 degrees | 4 pixels | Recovered |
| 10 degrees | 5 pixels | Wrong local minimum near 43.1 degrees and 0.746 pixel error |
| 12.5 degrees | 4 pixels | Recovered |
| 12.5 degrees | 5 pixels | Wrong local minimum near 43.1 degrees and 0.746 pixel error |
| 15 degrees | 4 pixels | Recovered |
| 15 degrees | 5 pixels | Wrong local minimum near 43.1 degrees and 0.746 pixel error |

Mixed-sign and simultaneous multi-axis cases in the final matrix recovered. A negative 15-degree/positive 5-pixel case recovered. The failure was therefore not an overflow, universal 43-degree attractor, or simple radial boundary.

## Multi-volume and route evidence

Five boundary seeds were tested on three deterministic asymmetric volumes using four routes:

1. direct joint LM;
2. direct joint LM with a cumulative 15-degree/5-pixel seed guard;
3. shift-only LM followed by joint LM; and
4. rotation-only LM followed by joint LM.

Results:

- Direct joint LM failed in three positive 5-pixel Gaussian cases but recovered all five cases for the other two volumes.
- The cumulative guard contained endpoints but did not guarantee recovery and sometimes blocked a useful trajectory.
- Shift-first then joint recovered all 15 volume/seed combinations to small truth error.
- Rotation-first then joint was unsafe, reaching errors of about 116 degrees and 11 pixels.

This establishes route-dependent local coupling. It does not authorize shift-first as production policy; the correct sequence remains a later contract decision.

## Objective-path evidence

The sampled straight path from each mechanism seed to truth was monotone, and the truth objective was approximately `1e-15` of the starting objective. Paths toward the wrong direct-joint endpoints were also monotone.

In the principal failure, LM reduced the objective by about 99.7 percent while ending about 43.1 degrees from truth. The wrong endpoint retained about `0.003` of the starting objective, while truth reached numerical zero.

Therefore an accepted reduction in the local particle objective is necessary for the implemented LM step but insufficient as a scientific pose-accuracy criterion.

## Evidence-package integrity

The final isolated package contained:

- 138 capture rows;
- 60 mechanism route endpoints;
- accepted-step five-vector trajectories;
- 1,230 objective-path samples;
- three deterministic volumes and their observations;
- starting/final predictions and residual images;
- logs, source snapshots, and SHA-256 manifests; and
- the unchanged default `pose_recovery` regression.

The Fortran cases returned zero. Reanalysis reported:

```text
POSE_CAPTURE_ANALYSIS: PASS
Package integrity: PASS
Scientific recovery warnings: 29
```

The warnings are observed capture failures, not missing evidence or harness failures.

## Historical production truth matrix

The earlier `reconstruct3D`-owned experiment completed 16 truth-controlled arms. Every enabled arm reported `REQUESTED 64 ACTIVE 64 BATCH 500` and took approximately 68.5--70.2 seconds. Parallel execution and harness integrity passed.

The scientific gate failed:

- clean exact rotation drift reached `0.0038695` rad;
- clean perturbed rotation improved by only about 17 percent;
- noisy exact FSC declined;
- noisy perturbed truth-average FSC declined by `0.00630986`; and
- nominal full/FSC05/FSC0143 arms were identical because the pose shell range was not changed by the `lp` label.

This result retired the post-`reconstruct3D` activation. It does not prove that every Cartesian local-pose formulation is scientifically invalid.

## Operator diagnostics retained as historical evidence

Focused diagnostics separated several effects:

- fast Kaiser--Bessel and normalized-stencil derivatives;
- inverse-envelope handling within the declared PCG operator;
- simulator padded projection, Fourier round trip, real-space crop, and native FFT;
- reconstructed-reference bias across iteration/support/lambda choices; and
- threaded equivalence and batching.

The matched-window diagnostic demonstrated that applying the simulator's padded inverse FFT, real-space crop, and native FFT consistently to prediction and Jacobian can reproduce that simulator and recover a local pose. Later review established that this finite-box transform is a simulator-specific model boundary, not automatically a missing operation in the declared production PCG model.

## Source locations at consolidation

- `src/main/volume/simple_reconstructor_pcg.f90:1188` — weighted pose normal terms.
- `src/main/volume/simple_reconstructor_pcg.f90:1244` — pose objective and gradient.
- `src/main/volume/simple_reconstructor_pcg.f90:1418` — bounded five-parameter LM.
- `src/main/volume/simple_reconstructor_pcg.f90:1679` — Cartesian forward plane.
- `src/main/volume/simple_reconstructor_pcg.f90:1766` — particle CTF/sigma transfer.
- `src/main/volume/simple_reconstructor_pcg.f90:1835` — observation whitening.
- `production/tests/simple_continuous_3D_pcg_refinement_rotation_test.f90` — componentwise derivative checks.
- `production/tests/simple_continuous_3D_pcg_refinement_pose_capture_test.f90` — capture matrix.
- `production/tests/simple_continuous_3D_pcg_refinement_pose_mechanism_test.f90` — route and objective-path diagnostics.
- `production/tests/continuous_3D_pcg_pose_capture/` — Oracle package runner and analyzer.

Line anchors describe commit `95e083817` and must be refreshed after source changes.

## Validation still missing

1. Independent accumulation of the complete shell-weighted objective.
2. Independent construction of $J^H J$ and $J^H r$.
3. Independent dense solution of one damped five-by-five LM proposal.
4. Explicit verification that rejected trials restore the complete input pose.
5. A CTF/sigma matrix covering constant, varying, zero-transfer, and invalid-variance cases.
6. A structurally separate Cartesian forward oracle rather than observations generated by the evaluator under test.
7. Repetition across multiple box sizes, shell ranges, and asymmetric volumes under one predeclared tolerance contract.

These are the active numerical-validation gaps in the new draft SPEC and PLAN.
