# Continuous five-parameter pose refinement numerical-validation PLAN

**Status:** DRAFT — NON-AUTHORITATIVE
**Date:** 2026-08-25
**SPEC:** [continuous_3D_pose_end_polishing_spec.md](continuous_3D_pose_end_polishing_spec.md) (`DRAFT`)
**Completed history:** [completed/continuous_3D_pose_end_polishing_history_and_handoff.md](completed/continuous_3D_pose_end_polishing_history_and_handoff.md)
**Completed evidence:** [completed/continuous_3D_pose_end_polishing_validation_evidence.md](completed/continuous_3D_pose_end_polishing_validation_evidence.md)

## Plan gate

This PLAN may be reviewed while the SPEC is DRAFT, but it does not authorize code changes. Resolve the SPEC blockers, mark the SPEC FINAL, reconcile this PLAN with that contract, and mark the PLAN FINAL before implementation.

## Implementation goal

Independently validate the shell-weighted Cartesian five-parameter pose evaluator and LM driver already present in the research code. The work ends with numerical evidence and a handoff. It does not add a production caller.

## Existing components

- `src/main/volume/simple_reconstructor_pcg.f90:1188` — fused pose normal terms.
- `src/main/volume/simple_reconstructor_pcg.f90:1244` — objective and gradient.
- `src/main/volume/simple_reconstructor_pcg.f90:1418` — bounded five-parameter LM.
- `src/main/volume/simple_reconstructor_pcg.f90:1679` — Cartesian forward plane.
- `src/main/volume/simple_reconstructor_pcg.f90:1766` — CTF/sigma transfer.
- `src/main/volume/simple_reconstructor_pcg.f90:1835` — observation whitening.
- `production/tests/simple_continuous_3D_pcg_refinement_rotation_test.f90` — componentwise derivative checks.
- `production/tests/simple_continuous_3D_pcg_refinement_pose_capture_test.f90` — capture matrix.
- `production/tests/simple_continuous_3D_pcg_refinement_pose_mechanism_test.f90` — route and objective-path diagnostics.
- `production/tests/continuous_3D_pcg_pose_capture/` — Oracle evidence package.

Refresh all line anchors against the implementation commit before changing code.

## Phase 0 — documentation consolidation

**Status:** COMPLETE

- Consolidated the completed experiment, handoff, literature review, clarifications, and retired contract under `doc/implementation_notes/completed/`.
- Replaced the decommissioned top-level SPEC and PLAN with this linked DRAFT pair.
- Updated historical validation runners so deleted note paths no longer break evidence packaging.

**Gate:** every predecessor is recoverable from Git at `95e083817`, active documents link to the archive, and repository searches find no broken predecessor path.

## Phase 1 — finalize the numerical contract

**Status:** BLOCKED BY DRAFT SPEC

1. Fix the objective, normal-term, and LM-step tolerance policy before results are generated.
2. Define the independent forward oracle's model, Fourier normalization, interpolation, CTF, sigma, shell range, and fixed-cell boundary.
3. Confirm that the current slice ends before standalone `refine3D` integration.
4. Mark the SPEC FINAL, reconcile this PLAN, and mark it FINAL.

**Gate:** no SPEC blocker remains and the acceptance criteria are testable without post-result decisions.

## Phase 2 — independent weighted-objective oracle

**Status:** NOT STARTED

Create a test-only objective accumulation that does not call the fused objective routine. For every selected Fourier sample, explicitly form the prediction, residual, shell variance, and scalar contribution.

Cover:

- constant variance;
- smoothly increasing variance by shell;
- strongly varying finite positive variance; and
- more than one nonstationary rotation/shift pose.

Record the independent and fused objective, absolute difference, relative difference, accumulated sample count, and worst shell contribution.

**Gate:** every valid case satisfies the finalized tolerance and no comparison is hidden behind an aggregate pass.

## Phase 3 — independent normal-equation and LM oracle

**Status:** NOT STARTED

Independently materialize the residual and five Jacobian columns, then accumulate

$$
g=J^H r,
\qquad
H=J^H J.
$$

Compare every gradient component and symmetric Hessian entry with the fused implementation. Reproduce scaling, diagonal flooring, damping, the five-by-five solve, predicted reduction, recomputed objective, gain ratio, and accept/reject decision in a separate dense reference calculation.

Include at least:

- one accepted step;
- one finite no-improvement step;
- one rejected trial; and
- one weak or singular system.

**Gate:** all components and the proposed step satisfy the finalized tolerance; every rejected outcome preserves the complete seed pose.

## Phase 4 — CTF and sigma matrix

**Status:** NOT STARTED

Repeat the objective, derivative, and normal-equation checks with:

- ordinary nonzero CTF;
- attenuated shells;
- exact CTF zeros;
- constant sigma;
- smoothly and strongly varying positive sigma; and
- zero, negative, NaN, infinity, or structurally invalid sigma inputs.

The test must distinguish a physically zero CTF contribution from an invalid variance.

**Gate:** valid cases remain finite and match the independent oracles; invalid variance inputs fail through the declared contract.

## Phase 5 — structurally separate forward oracle

**Status:** NOT STARTED

Implement a deliberately slow test-only Cartesian central-section generator that does not call `forward_plane`, `pose_normal_terms`, or their interpolation helper. State its exact normalization and coordinate conventions.

Compare:

- complex forward values;
- each of the five local derivative columns;
- fixed-cell finite differences; and
- explicitly identified stencil-switch cases.

If both the direct Cartesian and finite-box simulator models are retained, report them as different named models. Do not average or silently convert one into the other.

**Gate:** agreement passes away from stencil boundaries, switch cases are reported separately, and no oracle observation is generated by the evaluator under test.

## Phase 6 — configuration robustness

**Status:** NOT STARTED

Run the essential Phase 2--5 checks for at least:

- two box sizes;
- two Fourier shell ranges;
- three deterministic asymmetric volumes;
- multiple rotation axes and shift signs; and
- exact and nonstationary poses.

Reuse the existing capture volumes where appropriate, but keep the independent oracle structurally separate.

**Gate:** the same finalized tolerance policy passes every supported configuration without result-driven adjustment.

## Phase 7 — Oracle Linux verification and evidence package

**Status:** NOT STARTED

After source review and static checks:

1. the user compiles on Oracle Linux;
2. run each focused numerical case;
3. run the complete `simple_test_continuous_3D_pcg_refinement` mother suite;
4. package logs, machine-readable comparisons, source snapshots, configuration, and checksums in one timestamped directory; and
5. analyze every required comparison rather than accepting process exit alone.

**Gate:** the SPEC acceptance criteria pass, the mother suite has no failure, package integrity passes, and no production activation exists.

## Phase 8 — handoff only

**Status:** NOT STARTED

Update this PLAN with observed commands, package path, metrics, decisions, and residual limitations. If the numerical gate passes, open a new SPEC for a standalone continuous `refine3D` validation mode.

Do not add standalone mode code under this PLAN.

## Files expected to change after approval

The exact file set must be confirmed during Phase 1. Expected candidates are:

- `production/tests/simple_continuous_3D_pcg_refinement_rotation_test.f90` or a new dedicated oracle child;
- `production/tests/simple_continuous_3D_pcg_refinement_helpers.f90` for genuinely shared fixtures only;
- `production/tests/simple_test_continuous_3D_pcg_refinement.f90` for case registration;
- `production/CMakeLists.txt` if a new child source is added;
- `production/tests/continuous_3D_pcg_pose_capture/` or a new validation-only evidence runner; and
- this PLAN for phase evidence.

Production UI, parameters, commanders, strategies, and project persistence are not expected to change.

## Test strategy

| SPEC criterion | Verification |
| --- | --- |
| Weighted objective | Independent scalar accumulation across sigma profiles |
| $J^H r$ and $J^H J$ | Componentwise independent matrix/vector comparison |
| LM proposal | Separate dense five-by-five solve and gain-ratio calculation |
| Rollback | Exact before/after pose comparison for every rejection class |
| CTF/sigma | Valid and invalid matrix with explicit terminal outcomes |
| Forward model | Structurally separate generator and fixed-cell finite differences |
| Robustness | Multi-box, multi-range, multi-volume matrix |
| Regression | Existing focused cases plus complete mother suite |
| No integration | Source search for user-facing keys and production callers |

## Review roles

- **Planning:** reconcile the final SPEC with current source ownership and identify the smallest test-only file set.
- **Implementation:** add independent oracle code without changing the fused implementation merely to make comparisons pass.
- **Review:** independently inspect objective signs, units, conjugation, CTF/sigma order, matrix symmetry, damping, rollback, and oracle independence.
- **Verification:** perform static checks locally and authoritative compile/runtime validation on Oracle Linux after explicit user direction.

## Risks

- The oracle can accidentally reuse the implementation and recreate the inverse crime.
- A test can compare scaled terms while missing an error in the unscaled normal equations.
- Complex conjugation or factor-of-two conventions can make a numerically consistent but incorrectly interpreted gradient.
- CTF zeros can be confused with invalid sigma values.
- A tolerance derived after seeing results can make the gate meaningless.
- Stencil-switch cases can be misreported as derivative failures or silently excluded.
- Historical `pcg` names can cause the LM to be described as a PCG algorithm.

## Completion gate

This PLAN is complete only when every finalized SPEC criterion has retained Oracle evidence, review findings are resolved, the worktree contains no production activation, and the handoff clearly states what remains unproven.

Passing this PLAN validates an implementation component. It does not establish scientific benefit or authorize production pose polishing.
