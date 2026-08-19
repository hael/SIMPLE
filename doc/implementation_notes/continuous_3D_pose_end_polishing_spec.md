# Continuous 3-D pose refinement in `refine3D` SPEC

**Status:** DECOMMISSIONED (2026-08-19)

> This PFTC-based integration contract is no longer authoritative. Hans
> clarified that the immediate work is an isolated PCG Cartesian pose-capture
> experiment with no `refine3D` integration. The document is retained only as
> design history. See
> [continuous_3D_pose_polishing_hans_clarification_2026-08-19.md](continuous_3D_pose_polishing_hans_clarification_2026-08-19.md).

**PLAN:** [continuous_3D_pose_end_polishing_plan.md](continuous_3D_pose_end_polishing_plan.md)
**Frozen research proposal:** [continuous_3D_refinement_on_pcg_operator.md](continuous_3D_refinement_on_pcg_operator.md)

## Request and context

Hans Elmlund directed this work to the `refine3D` matching workflow. It must
not be presented as PCG and must not run from `reconstruct3D`. The feature must
follow the local-polishing position of `inpl_cont`, remain protected by
`pose_cont=yes|no`, and include a pure continuous mode for simulated-data
testing.

The earlier default-off `reconstruct3D` activation has been removed. Its
scientific evidence is historical input, not the contract for this feature.

## Initial state

`refine3D` performs a discrete or probabilistic pose search. With
`inpl_cont=yes`, it can then refine the selected in-plane angle and two shifts
while keeping the projection direction fixed.

The executed matcher uses PFTC reference planes created from preprocessed
even/odd reference volumes. The reprojection-model header owns the matching
shell range. The matcher objective supplies the CTF, shift phase, noise
weights, normalization, state, and half-set conventions.

## Required outcome

Add an opt-in, local five-parameter refinement

$$
q=(\omega_x,\omega_y,\omega_z,\delta t_x,\delta t_y),
\qquad
R=R_0\exp([\omega]_\times),
\qquad
t=t_0+\delta t.
$$

The first version supports two `refine3D` routes:

1. **Normal SHC route:** refine the selected SHC pose before that iteration
   reconstructs.
2. **Pure route:** `refine=pose_cont pose_cont=yes` starts from stored poses,
   performs no discrete search, and performs one pose pass followed by one
   reconstruction.

## Public contract

- `pose_cont=yes|no` is accepted by `refine3D`; the default is `no`.
- `refine=pose_cont` selects the pure route.
- The initial supported normal route is `refine=shc`.
- The initial objective is `objfun=euclid` with one state.
- The initial routes require normal even/odd ownership and `volrec=yes`.
- The reconstruction backend remains independent of pose refinement.
- `refine3D_auto`, `refine3D_multi`, `abinitio3D`, neighborhood search, and
  probabilistic search do not expose the feature in the first version.
- Unsupported combinations stop with a clear error before matching starts.

When both polishers are enabled in the normal SHC route, the order is:

1. discrete SHC selection;
2. existing `inpl_cont` refinement;
3. five-parameter `pose_cont` refinement;
4. reconstruction from the accepted pose.

The pure route does not run a separate `inpl_cont` solve.

## Scientific and numerical contract

The pose objective must be the executed PFTC raw-Euclidean matcher objective,
not the former PCG Cartesian-plane objective.

- Retain an immutable even/odd reference volume after the same mask, filter,
  centering, low-resolution parity blending, padding, and FFT preparation that
  creates the PFTC reference bank.
- Use the exact `kfromto` stored in the reprojection-model header.
- Use the matcher particle PFT, CTF, Fourier shift phase, sigma weighting,
  normalization, state, symmetry, and even/odd ownership.
- Evaluate reference values and rotation derivatives through the executed
  normalized fast Kaiser--Bessel interpolation path.
- Differentiate the executed polynomial and normalized stencil. Do not
  substitute the ideal Kaiser--Bessel derivative.
- Use right tangent-space rotation increments and pixel-unit image shifts.
- Accept a pose only after a finite, fully recomputed objective reduction.
- Restore the complete input pose after an invalid, singular, bounded-out, or
  non-improving result.
- Keep the accepted orientation in the current symmetry representative and
  refresh its nearest discrete `proj` metadata for restart.
- Use the existing matcher OpenMP particle loop. Shared and distributed
  workers must use the same scientific model.

The first optimizer is per-particle Gauss--Newton/Levenberg--Marquardt. It uses
at most eight iterations, a rotation scale and per-step rotation bound of
$1/\mathrm{msk\_crop}$ radians, a one-pixel per-step shift bound, a total
rotation bound equal to the current angular grid spacing, and a total shift
bound equal to `trs`. A trial step requires a positive reduction and an LM
gain ratio of at least `0.25`.

## Scope

- Exact five-parameter PFTC value, Jacobian, and local optimizer.
- Normal SHC integration and one-pass pure continuous refinement.
- Pose persistence through the normal `ptcl3D` project path.
- Shared and distributed execution equivalence.
- Focused numerical tests and matched simulated beta-gal validation.

## Out of scope

- Changing PCG, gridding, or ordinary `reconstruct3D` behavior.
- Replacing the discrete search.
- Multiple states, CTF refinement, or a new frequency schedule.
- CC or hybrid continuous-pose objectives.
- Automatic propagation to other refinement or ab-initio workflows.
- Production recommendation before the scientific gates pass.

## Acceptance criteria

1. The five analytic derivative columns agree with fixed-cell finite
   differences of the executed PFTC model to relative error at most `1e-2`.
2. At a discrete pose, the new evaluator reproduces the existing PFTC
   Euclidean objective within `1e-4`.
3. Exact synthetic poses remain within `1e-4` radians and `1e-3` pixels.
4. Known local perturbations reduce rotation and shift errors by at least
   50 percent and finish below `8e-4` radians and `2e-3` pixels.
5. Accepted objectives decrease monotonically; every step and total update
   obeys its bound; all terminal outcomes balance the attempted particles.
6. Weak, singular, invalid, rejected, symmetry-equivalent, and stencil-switch
   cases have defined, tested outcomes.
7. `pose_cont=no` preserves the established path and does not create or load
   the continuous-pose reference workspace.
8. The normal and pure routes persist accepted poses, and reconstruction
   consumes those poses.
9. Shared and distributed routes produce equivalent poses and terminal
   accounting.
10. A frozen simulated beta-gal matrix tests all four `inpl_cont`/`pose_cont`
    combinations from identical inputs. Exact-pose stationarity, perturbation
    recovery, metadata integrity, FSC area, cFAR, and runtime are reported.
    FSC-area decline greater than `0.01` fails the scientific gate.
11. Focused tests and the mother suite compile and pass on Oracle Linux.

## Constraints

- Preserve the frozen original proposal.
- Preserve ordinary `refine3D`, `inpl_cont`, and reconstruction behavior when
  `pose_cont` is omitted or `no`.
- Keep numerical algorithms in their owning projector/PFTC modules and
  workflow sequencing in matcher/search modules.
- Oracle Linux supplies compilation and runtime evidence.

## Review decisions awaiting confirmation

Before this SPEC becomes FINAL or FINAL (FROZEN), Hans and the requester must
confirm:

1. both the normal SHC route and the pure one-pass route;
2. the initial Euclidean, one-state, SHC-only boundary;
3. the order `discrete -> inpl_cont -> pose_cont` when both options are on;
4. the stated numerical and scientific acceptance gates.

Implementation must not start while this SPEC remains IN REVIEW.
