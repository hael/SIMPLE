# Continuous 3-D pose refinement: consolidated history and handoff

**Status:** COMPLETED HISTORICAL RECORD
**Consolidated:** 2026-08-25
**Pre-consolidation snapshot:** `95e083817`
**Current living development record:** [continuous_3D_pose_end_polishing.md](../continuous_3D_pose_end_polishing.md)
**Validation evidence:** [continuous_3D_pose_end_polishing_validation_evidence.md](continuous_3D_pose_end_polishing_validation_evidence.md)
**Scientific review:** [continuous_3D_pose_end_polishing_scientific_review.md](continuous_3D_pose_end_polishing_scientific_review.md)

## Purpose

This note consolidates the completed design history, research handoff, and two Hans clarifications that previously occupied seven top-level implementation notes. It is historical evidence, not an active implementation contract. The complete predecessor texts remain recoverable from Git at the pre-consolidation snapshot.

## Consolidated predecessor documents

- `continuous_3D_pose_capture_experiment_summary_2026-08-19.md`
- `continuous_3D_pose_end_polishing_handoff.md`
- `continuous_3D_pose_end_polishing_Literature Review and Recommended Path.md`
- `continuous_3D_pose_polishing_hans_clarification_2026-08-19.md`
- `continuous_3D_pose_polishing_hans_clarification_2026-08-20.md`
- the decommissioned `continuous_3D_pose_end_polishing_spec.md`
- the decommissioned `continuous_3D_pose_end_polishing_plan.md`

## Executive classification

The research produced a five-parameter local pose optimizer with three tangent-space rotation coordinates and two image shifts. The optimizer uses an analytic Jacobian and Levenberg--Marquardt (LM). It is not PCG. PCG is a reconstruction solver for the 3-D volume when particle poses are fixed.

The current numerical objective is a shell-variance-weighted Cartesian Fourier-plane Euclidean residual. The existing PFTC matcher objective is not a direct substitute because polar and Cartesian sampling have different integration measures. A PFTC formulation would require its own correctly weighted derivation, including the radial/Jacobian factor; the current Cartesian formulation does not need that polar factor.

The isolated derivative and pose-capture experiments are complete. They support the local analytic derivatives and expose a non-convex, volume-dependent capture basin. They do not authorize normal `refine3D` integration, replacement of `inpl_cont`, or automatic production acceptance.

## Timeline and decisions

### 1. Original reconstruction-coupled proposal

The frozen proposal [continuous_3D_refinement_on_pcg_operator.md](../continuous_3D_refinement_on_pcg_operator.md) began from the Cartesian Fourier workspace developed around `reconstructor_pcg`. It proposed alternating volume reconstruction and local pose updates.

This architectural origin led to the temporary name “PCG pose polishing.” That name is retired. The LM pose solve and the PCG volume solve have different unknowns, equations, and scientific roles.

### 2. Experimental `reconstruct3D` activation

An opt-in post-reconstruction route was implemented under `pcg_pose_polish=yes`, disabled by default. It reconstructed even/odd maps, polished particle poses against fixed same-half references, persisted accepted poses, and reconstructed again.

Operator, derivative, symmetry/gauge, persistence, parallel, and regression diagnostics passed their recorded Oracle gates. A 64-thread truth-controlled matrix completed all 16 arms with `REQUESTED 64 ACTIVE 64 BATCH 500`, but the scientific gate failed: exact clean poses drifted, perturbed rotations recovered weakly, and noisy FSC declined. Frequency-labelled arms were identical because the pose workspace did not receive distinct shell limits.

The route was not part of normal SIMPLE workflows. `abinitio3D` and `refine3D_auto` did not propagate the option, and broad propagation through every `reconstruct3D` caller would not define one coherent polishing stage.

### 3. Production activation removed

On 2026-08-19 Hans redirected the work away from `reconstruct3D` and toward the logical position occupied by continuous pose work in `refine3D`. Commit `1d32b830a` removed the user-visible `reconstruct3D` activation while retaining low-level numerical research routines and focused diagnostics.

The first replacement SPEC and PLAN proposed a PFTC-based `pose_cont` route. They were decommissioned before implementation because Hans clarified that the immediate task was implementation validation, not integration policy.

### 4. Hans clarification of 2026-08-19

The reported requirement was to test the existing five-parameter Cartesian central-section Jacobian and bounded LM in isolation:

1. use a known 3-D volume;
2. generate a clean observation at a known pose;
3. inject controlled rotation or shift errors one family at a time;
4. run the five-parameter local solve;
5. record capture success, failure, bounds, stencil behavior, and terminal state; and
6. produce one timestamped evidence directory with tables, projections, residuals, volumes, logs, and provenance.

No UI, `refine3D`, persistence, reconstruction, distributed, or workflow integration was authorized by that clarification.

### 5. Isolated pose-capture experiment completed

The final experiment contained 138 capture trials plus 60 route endpoints over three deterministic asymmetric volumes. Componentwise residual-Jacobian and objective-gradient finite differences passed. Separate rotations recovered through 15 degrees; separate shifts recovered reliably through about 5--6 pixels in the selected small fixture. Joint recovery was morphology- and route-dependent.

The most important failure lowered the Cartesian objective by about 99.7 percent but ended about 43.1 degrees from the truth pose. Shift-only optimization followed by joint LM recovered all 15 tested volume/seed combinations, while rotation-first staging was unsafe. These results demonstrate local non-convexity and rotation/shift coupling; they do not establish a universal capture radius or a production staging policy.

### 6. Hans clarification of 2026-08-20

Hans corrected the scientific focus:

- do not call the objective “raw Euclidean”; it depends on the variance vector over Fourier shells;
- do not substitute the existing PFTC objective directly because polar coordinates introduce a different integration measure;
- validate the implementation before choosing production policy;
- after complete numerical validation, consider a standalone continuous refinement mode in `refine3D` for simulated-data testing; and
- defer replacement, acceptance, frequency, bounds, rollback, and workflow policy until that standalone mode is understood.

## Settled current interpretation

### The five-parameter unknown

For an input pose $(R_0,t_0)$, the local coordinates are

$$
q=(\omega_x,\omega_y,\omega_z,t_x,t_y),
\qquad
R=R_0\exp([\omega]_\times).
$$

The rotations are local tangent coordinates; the shifts are image-plane pixels.

### The current objective

The intended Cartesian numerical fixture is

$$
\Phi(q)=\frac12\sum_{\mathbf{k}\in\Omega}
\frac{|M(q;\mathbf{k})-Y(\mathbf{k})|^2}{\sigma^2_{s(\mathbf{k})}},
$$

where the prediction uses the same particle-specific CTF, shift-phase, shell-variance, interpolation, and Fourier normalization conventions declared for the Cartesian operator.

The existing `build_transfer` and `whiten_observation` machinery supplies reusable CTF and sigma semantics. Reuse does not remove the need to verify that the pose prediction and derivatives apply those semantics in the same operation order.

### The optimizer

LM proposes a damped local step from the analytic residual Jacobian:

$$
(J^H J+\mu D)\delta=-J^H r.
$$

Accepted-step monotonicity is an internal numerical condition. It is not proof that the final pose is physically correct.

## Retired contracts

### Retired `reconstruct3D` contract

The old contract tied activation to `rec_backend=pcg` and a post-reconstruction same-half map. It is retired because pose refinement belongs to refinement logic, the normal workflows did not activate it, and its truth-controlled scientific gate failed.

### Retired PFTC `pose_cont` contract

The old replacement SPEC required an exact five-parameter PFTC evaluator, normal SHC integration, and a pure one-pass mode. It is retired because it selected an objective and workflow policy before the Cartesian numerical implementation was independently validated, and because the PFTC polar measure cannot be substituted without a new derivation.

The useful architectural ideas retained for later consideration are default-off activation, a standalone simulated-data mode, immutable reference ownership, complete-pose rollback, symmetry-aware persistence, and shared/distributed equivalence. None is authorized as current implementation work.

## Current boundary

The next active slice is numerical implementation validation only. It must independently verify the weighted objective, the five Jacobian columns, $J^H J$, $J^H r$, the damped LM proposal, rejection rollback, CTF/sigma edge cases, and a structurally separate forward oracle over multiple configurations.

Only after those gates pass should a new SPEC decide whether to add a standalone continuous `refine3D` mode. Normal refinement policy remains outside the present contract.

## Historical source anchors

These anchors identify the implementation as it stood when the pose-capture experiment was completed; line numbers must be refreshed before new work:

- `src/main/volume/simple_reconstructor_pcg.f90`: Cartesian forward plane, weighted pose terms, and `refine_pose_lm`;
- `production/tests/simple_continuous_3D_pcg_refinement_rotation_test.f90`: componentwise derivative checks;
- `production/tests/simple_continuous_3D_pcg_refinement_pose_capture_test.f90`: 138-trial capture matrix;
- `production/tests/simple_continuous_3D_pcg_refinement_pose_mechanism_test.f90`: multi-volume route and objective-path diagnostics;
- `production/tests/continuous_3D_pcg_pose_capture/`: Oracle evidence runner, analyzer, and package contract; and
- `production/tests/continuous_3D_pcg_pose_validation/`: historical operator and truth-matrix diagnostics.

## Handoff questions deferred to a later integration SPEC

1. Is a standalone continuous mode only a validation application, or may it later become a production refinement stage?
2. If normal refinement eventually uses continuous pose updates, which discrete or probabilistic seed supplies the local basin?
3. What reference, frequency, pose-prior, symmetry, and half-set contract is authoritative?
4. What external or held-out signal can reject a lower particle objective that worsens truth pose or FSC?
5. What memory, batching, and distributed lifetime owns the immutable Cartesian reference?

These are valid research questions, but they do not block the narrower numerical-validation SPEC except where that SPEC explicitly identifies them.
