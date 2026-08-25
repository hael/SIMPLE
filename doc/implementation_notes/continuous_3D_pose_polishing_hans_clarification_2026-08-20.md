# Continuous 3-D pose polishing: Hans clarification

**Status:** ACTIVE DIRECTION — NUMERICAL VALIDATION ONLY
**Date:** 2026-08-20
**Source:** Email response from Hans Elmlund

**Previous clarification:**
[continuous_3D_pose_polishing_hans_clarification_2026-08-19.md](continuous_3D_pose_polishing_hans_clarification_2026-08-19.md)

**Completed experiment:**
[continuous_3D_pose_capture_experiment_summary_2026-08-19.md](continuous_3D_pose_capture_experiment_summary_2026-08-19.md)

## Hans's response

> I am not sure what you mean by "raw" in this context. We always let the
> Euclidean objective depend on sigma (vector of variances per Fourier shell).
> The PFTC objective cannot be used directly, since it is polar coordinates and
> you will have a Jacobian scaling factor.
>
> You should not plan to add anything at this stage other than a test that
> verifies that your implementation is correct. You are thinking about policy
> now. You should think about implementation validation. When everything is
> tested, I would add a continuous standalone refinement mode to `refine3D`.
> That would be the best way of testing this implementation before we make any
> higher-level policy decisions. Reading up on the theory behind single-particle
> reconstruction and relating that to how things are implemented in SIMPLE also
> seems like an urgent task for you.

## Settled interpretation

### 1. The objective is sigma-weighted Euclidean

Do not describe the objective as "raw Euclidean." In this context, SIMPLE uses
a variance for each Fourier shell. Let $s(\mathbf{k})$ identify the shell of
Cartesian Fourier sample $\mathbf{k}$, and let $M(q;\mathbf{k})$ be the predicted
CTF- and shift-modulated Fourier value at pose
$q=(\omega_x,\omega_y,\omega_z,t_x,t_y)$. The intended objective is

$$
\Phi(q)
=\frac{1}{2}\sum_{\mathbf{k}\in\Omega}
\frac{\left|M(q;\mathbf{k})-Y(\mathbf{k})\right|^2}
{\sigma^2_{s(\mathbf{k})}}.
$$

Equivalently, the model and observation are whitened by
$1/\sqrt{\sigma^2_s}$ before their Euclidean distance is calculated.

The current PCG workspace constructs the CTF/sigma transfer in
[`build_transfer`](../../src/main/volume/simple_reconstructor_pcg.f90#L1781)
and whitens the observation in
[`whiten_observation`](../../src/main/volume/simple_reconstructor_pcg.f90#L1853).

### 2. The PFTC objective is not a direct substitute

PFTC samples a polar Fourier representation. A direct sum over polar samples
does not have the same measure as a sum over a uniform Cartesian plane. The
polar-coordinate formulation introduces a Jacobian or radial weighting factor.

Therefore, do not replace the Cartesian pose objective with the current PFTC
objective and do not plan PFTC rescoring at this stage. This is an
implementation-validation task, not an objective-selection policy task.

### 3. Production policy is deferred

Do not decide the following now:

- whether continuous refinement replaces `inpl_cont`;
- when normal `refine3D` invokes it;
- direct joint versus staged production operation;
- production movement bounds;
- production acceptance or rollback policy;
- FSC or held-out acceptance criteria.

These decisions come only after the numerical implementation and a standalone
continuous refinement mode have been validated.

### 4. The next integration target is standalone validation

After the numerical implementation is fully tested, the next integration step
is a standalone continuous refinement mode in `refine3D`. Its purpose is to test
the implementation through the real SIMPLE data path. It is not yet a policy for
normal refinement or a replacement for discrete search.

No standalone-mode implementation is authorized by this note. It will require a
new concise SPEC and PLAN after the numerical validation gate passes.

## Validation already completed

The current evidence includes:

- finite differences for every rotation and shift Jacobian column;
- independent finite differences for every objective-gradient component;
- a nonconstant shell-variance vector and CTF/sigma-weighted derivative checks;
- exact-pose stationarity;
- controlled rotation-only, shift-only, joint, mixed-sign, and multi-axis trials;
- three deterministic asymmetric volumes;
- LM trajectory, bounds, rollback, and terminal-state evidence;
- the default `pose_recovery` regression; and
- a ten-case Oracle Linux mother-suite pass.

The shell-dependent weighted derivative test begins in
[`simple_continuous_3D_pcg_refinement_rotation_test.f90`](../../production/tests/simple_continuous_3D_pcg_refinement_rotation_test.f90#L137).

This evidence strongly supports the local analytic derivatives. It does not yet
independently validate the complete weighted objective, normal equations, and LM
step.

## Proposed validation-only work plan

Do not add UI, production workflow activation, persistence policy, or distributed
integration during this plan. Hans did not enumerate the tests below; they are the
proposed response to his request for implementation validation.

### Test 1 — weighted-objective oracle

Independently accumulate

$$
\frac{1}{2}\sum_{\mathbf{k}}
\frac{|M(q;\mathbf{k})-Y(\mathbf{k})|^2}{\sigma^2_{s(\mathbf{k})}}
$$

and compare it with the objective returned by the fused implementation. Cover a
constant variance and multiple nonconstant positive shell-variance profiles.

**Gate:** the independent and fused objectives agree within a predeclared
floating-point tolerance.

### Test 2 — Gauss--Newton and LM oracle

Independently construct the five Jacobian columns, residual vector, $J^H J$, and
$J^H r$. Compare them with the gradient and approximate Hessian used by the
implementation. Then compare one damped LM proposal with an explicit dense
five-by-five reference solve.

**Gate:** gradient, approximate Hessian, and proposed step agree componentwise;
rejected trials leave the complete input pose unchanged.

### Test 3 — sigma and CTF matrix

Repeat the weighted objective, derivative, and normal-equation checks with:

- constant variance;
- smoothly increasing variance by shell;
- strongly varying but finite positive variance;
- attenuated shells and CTF zeros; and
- invalid variance inputs, which must be rejected explicitly.

**Gate:** all valid profiles remain finite and agree with their independent
oracles; invalid profiles fail through the declared numerical contract.

### Test 4 — independent forward oracle

Implement a deliberately slow and structurally separate reference calculation
of the same Cartesian central-section model. Do not generate the oracle
observation through the evaluator under test.

**Gate:** forward values and local pose derivatives agree away from interpolation
stencil boundaries. Stencil-boundary cases are reported separately and are not
misrepresented as differentiable points.

### Test 5 — configuration robustness

Repeat the essential oracle checks for more than one box size, Fourier range, and
asymmetric volume. This tests implementation generality; it does not define a
production capture range.

**Gate:** every supported configuration passes the same mathematical contract
without changing tolerances after results are observed.

## Theory-to-code audit

In parallel with the tests, prepare a concise scientific note that maps these
topics to current SIMPLE source:

1. Fourier-slice forward model;
2. CTF multiplication and translation phase;
3. shell-dependent Gaussian noise and whitening;
4. Cartesian versus polar Fourier integration measures;
5. tangent-space rotation coordinates and pixel-shift units;
6. Kaiser--Bessel interpolation and fixed-cell derivative limits;
7. reference-volume preparation and Fourier workspace lifetime; and
8. half-set ownership and the limits of same-data objective improvement.

Use primary literature for the mathematical model and current SIMPLE source for
implementation claims. Keep mathematical facts, source observations, and untested
hypotheses explicitly separate.

## Gate to the next phase

The next phase begins only when Tests 1--5 pass on Oracle Linux and the
theory-to-code audit contains no unresolved operator-definition mismatch.

At that point:

1. write a new concise SPEC for a standalone continuous `refine3D` mode;
2. write its linked implementation PLAN;
3. keep the mode disabled by default;
4. initialize it from supplied poses without normal discrete candidate search;
5. validate exact and perturbed simulated poses through the real SIMPLE path; and
6. defer higher-level refinement policy until that standalone mode is understood.
