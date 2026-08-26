# Continuous five-parameter pose refinement numerical-validation PLAN

**Status:** DRAFT — NON-AUTHORITATIVE  
**Date:** 2026-08-25  
**SPEC:** [continuous_3D_pose_end_polishing_spec.md](continuous_3D_pose_end_polishing_spec.md) (`DRAFT`)  
**Completed history:** [completed/continuous_3D_pose_end_polishing_history_and_handoff.md](completed/continuous_3D_pose_end_polishing_history_and_handoff.md)  
**Completed evidence:** [completed/continuous_3D_pose_end_polishing_validation_evidence.md](completed/continuous_3D_pose_end_polishing_validation_evidence.md)

## Plan gate

This PLAN may be reviewed while the SPEC is DRAFT, but it does not authorize code changes. Resolve the SPEC blockers, mark the SPEC FINAL, reconcile this PLAN with that contract, and mark the PLAN FINAL before implementation.

## Implementation goal

Move the shell-weighted Cartesian five-parameter pose evaluator and LM driver out of the PCG reconstruction object into a dedicated particle-pose numerical module, extract only proven state-independent overlap into a neutral shared module, and independently validate the resulting implementation. The work ends with numerical evidence and a handoff. It does not add a production caller.

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

These locations describe the current research implementation, not its final ownership. The pose objective, normal terms, and LM driver must move to a particle-pose module. The Phase 1 audit will decide whether lower-level Fourier operations remain reconstruction-owned or move to a neutral shared module.

## Phase 0 — documentation consolidation

**Status:** COMPLETE

- Consolidated the completed experiment, handoff, literature review, clarifications, and retired contract under `doc/implementation_notes/completed/`.
- Replaced the decommissioned top-level SPEC and PLAN with this linked DRAFT pair.
- Updated historical validation runners so deleted note paths no longer break evidence packaging.

**Gate:** every predecessor is recoverable from Git at `95e083817`, active documents link to the archive, and repository searches find no broken predecessor path.

## Phase 1 — finalize the numerical contract

**Status:** BLOCKED BY DRAFT SPEC

1. Define the non-gating calibration fixtures, comparison scales, safety-factor rule, and point at which the objective, normal-term, and LM-step tolerances become frozen.
2. Define the independent forward oracle's model, Fourier normalization, interpolation, CTF, sigma, shell range, and fixed-cell boundary.
3. Produce a routine-by-routine ownership map for the current pose objective, LM driver, Fourier gather/interpolation, transfer construction, observation weighting, workspace data, and lifecycle methods.
4. Select a dedicated particle-pose module boundary and identify only genuinely state-independent operations for a neutral shared module. The particle module must not depend on `reconstructor_pcg`.
5. Keep this slice closed before standalone `refine3D` integration. The later SPEC will make legacy `inpl_cont` and the five-parameter polisher mutually exclusive user choices in the same logical workflow position.
6. Mark the SPEC FINAL, reconcile this PLAN, and mark it FINAL.

**Gate:** no SPEC blocker remains and the acceptance criteria are testable without post-result decisions.

## Phase 2 — separate numerical ownership

**Status:** NOT STARTED

Create the dedicated particle-pose numerical module defined in Phase 1. Move the five-parameter objective, normal terms, bounded LM solve, acceptance/terminal accounting, and rollback into that owner without intentional mathematical changes.

For every overlapping lower-level operation:

- keep it in `reconstructor_pcg` when it depends on PCG solver state, volume reconstruction policy, or reconstruction lifecycle;
- move it to a neutral shared module only when its inputs, outputs, normalization, precision, and failure contract are independent of both owners; and
- update both callers to use the shared operation rather than maintaining copied implementations.

The dependency direction must be:

```text
neutral Fourier numerics <- reconstructor_pcg
neutral Fourier numerics <- particle-pose module
```

The particle-pose module must not import or contain a `reconstructor_pcg` object.

**Gate:** source review confirms the dependency boundary, existing focused component tests still exercise the moved implementation, and lightweight checks find no production caller or user-visible option.

## Phase 3 — numerical calibration and independent weighted-objective oracle

**Status:** NOT STARTED

Create a test-only objective accumulation that does not call the fused objective routine. For every selected Fourier sample, explicitly form the prediction, residual, shell variance, and scalar contribution.

Run named, non-gating calibration fixtures first. Record the numerical differences and relevant scales for the objective, each gradient component, every Hessian entry, and the scaled LM step. Use those measurements together with the implementation precision and a documented safety-factor rule to set a combined absolute-plus-relative comparison of the form

$$
|a-b| \leq \mathrm{atol}+\mathrm{rtol}\max(|a|,|b|).
$$

Freeze the resulting metric definitions and tolerances before running distinct formal acceptance fixtures. Calibration results characterize numerical scale; they do not themselves count as acceptance evidence.

Cover:

- constant variance;
- smoothly increasing variance by shell;
- strongly varying finite positive variance; and
- more than one nonstationary rotation/shift pose.

Record the independent and fused objective, absolute difference, relative difference, accumulated sample count, and worst shell contribution. Preserve the calibration records, tolerance derivation, and frozen values in the evidence package.

**Gate:** the calibration and acceptance fixture lists are disjoint, the tolerance policy is frozen before acceptance begins, every valid acceptance case satisfies it, and no comparison is hidden behind an aggregate pass.

## Phase 4 — independent normal-equation and LM oracle

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

## Phase 5 — CTF and sigma matrix

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

## Phase 6 — independent forward-oracle hierarchy

**Status:** NOT STARTED

Implement four explicitly different forward paths:

1. retain the production fast Kaiser--Bessel gather as the implementation under test;
2. implement a deliberately simple test-only Cartesian central-section interpolator that does not call `forward_plane`, `pose_normal_terms`, or the production stencil helper;
3. implement a brute-force DFT from a real-space array directly to the required rotated Fourier coordinates, with no FFT, Fourier-grid plane extraction, or Kaiser--Bessel interpolation; and
4. retain a separately named finite-box projection-to-FFT path for external discrete-model comparison.

Hard-limit the brute-force DFT to declared tiny test boxes, initially $8^3$ and $16^3$. Match the executed SIMPLE convention: zero-based $x_j=j-N/2$, spatial origin at one-based index $N/2+1$, frequency $([h,k,0]R)/N$, negative-exponential forward transform, and native $1/N^3$ scale after the padded path's `padf**3` correction. Confirm shift phase and CTF order before evaluating comparisons.

Accumulate the normalized direct spatial-frequency derivative $-(2\pi i/N^3)\sum_{\mathbf{x}}\mathbf{x}v(\mathbf{x})\exp(-2\pi i\boldsymbol{\xi}\cdot\mathbf{x})$ with the DFT value. Chain it through the declared right-tangent rotation convention to obtain three rotation columns without differentiating the production interpolation path. Derive the two shift columns separately from the declared Fourier phase.

Add analytic fixtures for:

- a centered unit delta;
- an off-center unit delta;
- two separated points for phase superposition;
- a three-or-more-point non-collinear unequal-amplitude set for rotational derivatives;
- a Gaussian with its sampled/truncated versus continuous interpretation stated; and
- a sphere used only for radial amplitude and normalization, not rotational recovery.

Current asymmetric Gaussian-mixture and mechanism volumes remain useful robustness fixtures, but they do not count as these analytic Fourier identities. In particular, the shared Gaussian fixture uses one-based center $N/2+0.5$, half a voxel from SIMPLE's FFT origin at $N/2+1$, and therefore is not a zero-phase centered object.

Compare and report separately:

- analytic values against their closed-form phase or radial expectations;
- brute-force DFT values against the production gather as a complete interpolation-error measurement;
- slow Cartesian interpolation against production as an interpolation-implementation check;
- finite-box projection against the other paths as a distinct discrete imaging model;
- each of the five local derivative columns;
- fixed-cell finite differences; and
- explicitly identified stencil-switch cases.

Do not apply one tolerance to all paths. Freeze a named scale and tolerance for each relationship after calibration and before formal acceptance. Do not average or silently convert the direct DFT, reference interpolation, and finite-box results.

**Gate:** analytic signs, origins, and scales pass; the direct-DFT interpolation error remains within its frozen contract; the reference interpolator agrees with production within its separately frozen contract away from stencil boundaries; finite-box differences are reported under their distinct model; switch cases are reported separately; and no oracle observation is generated by the evaluator under test.

## Phase 7 — configuration robustness

**Status:** NOT STARTED

Run the essential Phase 3--6 checks for at least:

- two box sizes;
- two Fourier shell ranges;
- three deterministic asymmetric volumes;
- multiple rotation axes and shift signs; and
- exact and nonstationary poses.

Reuse the existing capture volumes where appropriate, but keep the independent oracle structurally separate.

**Gate:** the same finalized tolerance policy passes every supported configuration without result-driven adjustment.

## Phase 8 — Oracle Linux verification and evidence package

**Status:** NOT STARTED

After source review and static checks:

1. the user compiles on Oracle Linux;
2. run each focused numerical case;
3. run the complete `simple_test_continuous_3D_pcg_refinement` mother suite;
4. package logs, machine-readable comparisons, source snapshots, configuration, and checksums in one timestamped directory; and
5. analyze every required comparison rather than accepting process exit alone.

**Gate:** the SPEC acceptance criteria pass, the mother suite has no failure, package integrity passes, and no production activation exists.

## Phase 9 — handoff only

**Status:** NOT STARTED

Update this PLAN with observed commands, package path, metrics, decisions, and residual limitations. If the numerical gate passes, open a new SPEC for a standalone continuous `refine3D` validation mode. That SPEC will follow the logical `inpl_cont` integration position and define the legacy `inpl_cont` route and the five-parameter polisher as mutually exclusive choices.

Do not add standalone mode code under this PLAN.

## Files expected to change after approval

The exact file set must be confirmed during Phase 1. Expected candidates are:

- a new dedicated particle-pose numerical module, with a pose/LM name rather than a PCG name;
- an optional neutral Cartesian Fourier-sampling module if the ownership audit proves shared state-independent operations;
- `src/main/volume/simple_reconstructor_pcg.f90` to remove particle-owned routines and call any approved shared operations without changing reconstruction behavior;
- `production/tests/simple_continuous_3D_pcg_refinement_rotation_test.f90` or a new dedicated oracle child;
- a dedicated tiny-box direct-DFT/analytic-object test module rather than embedding the oracle in a general helper;
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
| Analytic Fourier objects | Closed-form delta/point phases and qualified Gaussian/sphere checks |
| Brute-force DFT | Tiny-box direct sums at rotated continuous frequencies, with no FFT or interpolation |
| Reference interpolation | Structurally separate slow gather and fixed-cell finite differences |
| Finite-box model | Separately named projection-to-FFT comparison without forced equivalence |
| Robustness | Multi-box, multi-range, multi-volume matrix |
| Ownership | Dependency inspection: pose module and reconstructor depend only on the neutral shared layer, never on each other |
| Pose regression | Existing focused cases plus complete mother suite |
| Reconstruction regression | Existing PCG numerical and reconstruction-contract cases after shared extraction |
| No integration | Source search for user-facing keys and production callers |

## Review roles

- **Planning:** reconcile the final SPEC with current source ownership and approve the particle, reconstruction, and neutral shared boundaries routine by routine.
- **Implementation:** separate numerical ownership without intentional mathematical changes, then add independent oracle code without changing the implementation merely to make comparisons pass.
- **Review:** independently inspect objective signs, units, conjugation, CTF/sigma order, matrix symmetry, damping, rollback, and oracle independence.
- **Verification:** perform static checks locally and authoritative compile/runtime validation on Oracle Linux after explicit user direction.

## Risks

- The oracle can accidentally reuse the implementation and recreate the inverse crime.
- A continuous Gaussian or sphere formula can be incorrectly treated as the exact DFT of a sampled, truncated voxel array.
- A symmetric fixture can pass amplitude tests while providing no information about rotational derivatives.
- A test can compare scaled terms while missing an error in the unscaled normal equations.
- Complex conjugation or factor-of-two conventions can make a numerically consistent but incorrectly interpreted gradient.
- CTF zeros can be confused with invalid sigma values.
- Reusing calibration fixtures as acceptance fixtures, or changing a frozen tolerance after seeing acceptance failures, can make the gate meaningless.
- Stencil-switch cases can be misreported as derivative failures or silently excluded.
- Historical `pcg` names can cause the LM to be described as a PCG algorithm.
- Moving routines before identifying hidden workspace dependencies can change normalization or lifecycle behavior.
- Copying overlapping routines instead of extracting a neutral implementation can let reconstruction and particle conventions drift.

## Completion gate

This PLAN is complete only when every finalized SPEC criterion has retained Oracle evidence, review findings are resolved, particle-pose numerics no longer depend on `reconstructor_pcg`, ordinary reconstruction remains within its existing contract, the worktree contains no production activation, and the handoff clearly states what remains unproven.

Passing this PLAN validates an implementation component. It does not establish scientific benefit or authorize production pose polishing.
