# Continuous 3-D pose refinement in `refine3D` PLAN

**Status:** DECOMMISSIONED (2026-08-19)

> This PFTC-based integration plan is no longer authoritative. No phase in this
> document authorizes implementation. The current work is limited to the
> isolated experiment in
> [continuous_3D_pose_polishing_hans_clarification_2026-08-19.md](continuous_3D_pose_polishing_hans_clarification_2026-08-19.md).

**Former SPEC:** [continuous_3D_pose_end_polishing_spec.md](continuous_3D_pose_end_polishing_spec.md) (`DECOMMISSIONED`)
**Frozen research proposal:** [continuous_3D_refinement_on_pcg_operator.md](continuous_3D_refinement_on_pcg_operator.md)

## Plan gate

This PLAN records the agreed technical design. It does not authorize
implementation while the SPEC is IN REVIEW. After Hans and the requester
confirm the review decisions, mark the SPEC `FINAL` or `FINAL (FROZEN)`,
reconcile this PLAN with that contract, and only then start Phase 1.

The former `reconstruct3D` activation is retired. Commit `1d32b830a` removed
its UI, parameter, validation, commander, and distributed-job activation. The
low-level research diagnostics can remain until the PFTC replacements pass,
but the new feature must not depend on the PCG reconstruction objective.

## Implementation goal

Add default-off five-parameter local pose refinement to `refine3D` using its
exact PFTC Euclidean objective and frequency range. Support:

- a normal SHC route after candidate selection; and
- a pure one-pass `refine=pose_cont` route from stored poses.

When `inpl_cont=yes` and `pose_cont=yes` are both active in SHC, run the
existing in-plane solve first and use its accepted pose as the seed for the
five-parameter solve.

## Phase 0 — finalize the contract

**Status:** IN REVIEW

1. Review the concise SPEC with Hans and the requester.
2. Resolve only changes that affect scope, interfaces, numerical policy, or
   acceptance criteria.
3. Mark the SPEC `FINAL` or `FINAL (FROZEN)` after approval.
4. Reconcile this PLAN and mark it `FINAL` before implementation.

**Gate:** the SPEC contains no unresolved blocking decision and both documents
describe the same two routes and acceptance criteria.

## Phase 1 — policy and scaffold

**Status:** NOT STARTED

Add the user-visible and typed policy without adding pose numerics.

- Add `pose_cont=yes|no`, default `no`, to `refine3D` only.
- Add `pose_cont` to the `refine` choices for the pure route.
- Accept only:
  - `refine=shc pose_cont=yes objfun=euclid nstates=1`; or
  - `refine=pose_cont pose_cont=yes objfun=euclid nstates=1 maxits=1`.
- Require `volrec=yes` and the normal even/odd particle partition.
- Do not expose the option through automatic, multi-state, probabilistic,
  neighborhood, or ab-initio workflows.
- Keep `pose_cont=no` free of new allocations, files, calls, and metadata.

Create a mother test with separate policy, numerical, route, recovery, and
metadata children. Put only shared fixtures and assertions in helpers.

**Gate:** source checks pass; Oracle compilation passes; default-off and every
supported or rejected command shape behave as declared.

## Phase 2 — immutable refine3D reference workspace

**Status:** NOT STARTED

Extend the reprojection-model materialization path so the full-pose evaluator
uses the same reference source as the PFTC reference bank.

1. Capture each single-state even/odd reference after existing masking,
   filtering, centering, low-resolution parity blending, and padding
   preparation.
2. Write a versioned immutable reference artifact with parity, state,
   dimensions, sampling, symmetry, and `kfromto` provenance.
3. Load the artifact once per matcher worker when `pose_cont=yes`.
4. Build read-only padded Fourier projectors shared by the worker's OpenMP
   particle threads.
5. Reject a stale or incompatible artifact before processing particles.

The pure and SHC routes use the same workspace. The workspace does not depend
on `rec_backend`.

**Gate:** shared and distributed workers report identical artifact provenance;
the discrete reference sampled from the workspace matches the stored PFTC
reference within `1e-4`.

## Phase 3 — exact five-parameter PFTC evaluator

**Status:** NOT STARTED

Add a single-pose projector API beside the executed optimized polar projector.
For each polar sample, evaluate the normalized fast-KB gather and its spatial
gradient from the immutable Fourier volume.

For

$$
R(\omega)=R_0\exp([\omega]_\times),
$$

form the three rotation columns through the derivative of the rotated sample
coordinate. Form the two shift columns through the existing Fourier shift
phase. Apply the same CTF, sigma weighting, shell weights, normalization, and
operation order as the PFTC raw-Euclidean objective.

The evaluator returns:

- the fully recomputed scalar objective;
- the five-component gradient;
- the symmetric Gauss--Newton matrix;
- minimum stencil margin and stencil-switch telemetry; and
- a finite/identifiable status.

**Gate:** fixed-cell finite differences pass at multiple poses and orientations
with relative error at most `1e-2`; the discrete-pose objective identity is at
most `1e-4`; switch-crossing measurements are reported separately.

## Phase 4 — bounded local LM

**Status:** NOT STARTED

Implement a neutral per-particle five-parameter LM driver. Do not call the PCG
Fourier workspace.

- Use right rotation increments and pixel shifts.
- Use at most eight iterations.
- Scale and bound each rotation step by `1/msk_crop` radians.
- Bound each shift step by one pixel.
- Bound total rotation from the seed by the current angular grid spacing.
- Bound total shift from the seed by `trs`.
- Use initial damping `1e-3` and a diagonal floor based on the executed
  Gauss--Newton matrix.
- Accept only a positive recomputed reduction with gain ratio at least `0.25`.
- Decrease damping after a strong accepted step and increase it after a
  rejected or invalid step.
- Restore the complete seed after no improvement, weak identifiability,
  invalid numerics, or bound failure.

Use explicit terminal states: improved, unchanged, unreliable, bound-rejected,
invalid, and iteration-limit.

**Gate:** exact poses remain stationary; known joint perturbations recover;
weak and singular cases remain unchanged; accepted objectives are monotone;
serial and threaded results agree.

## Phase 5 — normal SHC route

**Status:** NOT STARTED

Insert the local solve after SHC has selected its final state, projection,
in-plane cell, and shift.

- With `inpl_cont=no`, seed the five-parameter solve from the discrete result.
- With `inpl_cont=yes`, run the existing in-plane solve first and seed from its
  accepted result.
- Commit all three Euler angles and two shifts only after full-pose acceptance.
- Preserve state and parity.
- Refresh `proj` to the closest current projection for restart metadata.
- Let the existing partial-reconstruction path consume the final pose.

**Gate:** all four `inpl_cont`/`pose_cont` combinations follow their declared
routes; terminal totals balance active particles; default-off output retains
the established metadata contract.

## Phase 6 — pure one-pass route

**Status:** NOT STARTED

Add a `strategy3D_pose_cont` route that seeds directly from each stored
particle pose.

- Do not run discrete, probabilistic, neighborhood, or in-plane-grid search.
- Do not run a separate `inpl_cont` solve.
- Keep state and parity fixed.
- Perform one local pose pass and one normal reconstruction.
- Reject `maxits` values other than `1` in this initial mode.
- Persist accepted poses and verify that reconstruction consumes them.

**Gate:** route telemetry proves zero discrete candidates were evaluated;
exact stored poses remain stationary; injected perturbations recover; the
reconstructed project passes metadata and restart validation.

## Phase 7 — Oracle and scientific validation

**Status:** NOT STARTED

Run compilation and tests on Oracle Linux only after all source checks pass.

### Numerical and regression package

Run the focused projector, objective, LM, SHC-route, pure-route, metadata, and
parallel-equivalence cases. Then run the complete mother suite.

### Frozen simulated beta-gal matrix

Fork every arm from one immutable project and fixed random state:

| Arm | `inpl_cont` | `pose_cont` | Route |
|---|---:|---:|---|
| Baseline | no | no | discrete SHC only |
| Existing | yes | no | SHC plus in-plane polish |
| New | no | yes | SHC plus five-parameter polish |
| Combined | yes | yes | SHC, in-plane, then five-parameter polish |

Also run pure-mode exact-pose and perturbed-pose controls from the same
completed checkpoint.

Collect:

- rotation and shift error after symmetry minimization and global-gauge
  correction;
- exact-pose drift and perturbed-pose recovery;
- objective before and after;
- FSC area, cFAR, FSC resolutions, and shell curves;
- terminal outcome totals and step bounds;
- shared/distributed pose agreement;
- metadata and reconstruction-consumption checks; and
- runtime and active OpenMP worker count.

**Scientific gate:** exact-pose and recovery thresholds from the SPEC pass,
accepted objectives remain monotone, FSC-area decline is not greater than
`0.01`, and no metadata, ownership, or shared/distributed discrepancy occurs.

Passing this gate validates the first `refine3D` research slice. It does not
authorize production recommendation or propagation to other workflows.

## Deferred propagation

After Phase 7, prepare a separate SPEC/PLAN amendment before adding any of:

- CC or hybrid objectives;
- probabilistic, neighborhood, or greedy search modes;
- repeated pure pose/reconstruction alternation;
- multiple states;
- `refine3D_auto`, `refine3D_multi`, or `abinitio3D`; or
- a dataset-level acceptance or rollback policy.

## Verification and commit policy

- Keep source inspection, compilation, runtime, and scientific acceptance as
  separate evidence levels.
- Do not weaken a threshold to make a failing experiment pass.
- Record each completed phase, important file and current line references,
  commands, results, and outstanding gates in this PLAN.
- Use focused commits only after the phase's declared Oracle gate passes.
- Preserve unrelated dirty worktree changes.
