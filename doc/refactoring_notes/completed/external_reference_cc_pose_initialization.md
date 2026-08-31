# Initialize Poses Against External References Before Euclidean Refinement

## Status

Implemented in the working tree, 2026-08-31. The shared external-reference CC
pose-initialization service and its three callers are present; compilation and
runtime validation remain user-owned. This note is a focused companion to
`refine3d_state_and_reference_workflow_refactor.md`. That document owns the
application boundaries and approved names; this document owns only the policy
for introducing reference volumes that may not share the particles' alignment
history.

## Approved Decision

When reference volumes are independent of the particle alignment, SIMPLE does
not begin with noise-normalized Euclidean matching. It first runs a
reference-conditioned pose-initialization stage with cross-correlation as one pass
over the existing capped pose initialization cohort of at most 100,000 particles:

1. normalize, validate, and low-pass filter every supplied reference to a
   common 15 A limit;
2. keep those references fixed;
3. select the normal pose initialization cohort, capped at 100,000 particles, and
   match that cohort competitively with `objfun=cc` in exactly one pass
   (`volrec=no`, `trail_rec=no`, no fractional-update scheme);
4. accumulate Euclidean residual sigma contributions during that CC pass after
   each hard pose/state assignment;
5. reconstruct checkpoint state maps from the committed assignments and
   consolidate the accumulated sigma statistics;
6. switch to `objfun=euclid` and start the ordinary frequency-marched
   refinement.

The pose initialization invariant is defined over the selected cohort, not the full
active population when more than 100,000 particles are available:

```text
one CC pose initialization pass  =>  each particle in the capped cohort matched once
                              against the fixed reference set
```

No fractional pose initialization schedule, pose initialization coverage counter, or
missing-particle sweep is required. The 100,000-particle cap is the existing
scientific sampling policy, not an update fraction to be completed over later
passes. This is the SIMPLE analogue of RELION's first-iteration
cross-correlation policy, adapted to SIMPLE's existing capped initial-mapping
cohort.

References from the same project lineage, accompanied by a trusted particle
pose scaffold, may start directly with Euclidean refinement. The policy is
therefore based on reference/pose provenance, not merely on whether `vol1` was
supplied.

## Why This Is Needed

### External references are not initially on the data model's scale

An independently produced map may differ from the particles' working model in
normalization, amplitude scale, filtering, masking, and reconstruction
convention. A Euclidean residual interprets those differences as evidence.
Cross-correlation is the safer pose initialization objective because the initial task
is to find pose and state correspondence before estimating meaningful residual
noise statistics.

RELION exposes the same distinction through `do_firstiter_cc`: its source
describes the first iteration as using the signal cross-product instead of the
Gaussian objective. Its refinement tutorial also treats first-iteration
greyscale correction and initial low-pass filtering as explicit reference
initialization concerns:

- <https://github.com/3dem/relion/blob/master/src/ml_optimiser.h>
- <https://relion.readthedocs.io/en/release-5.0/SPA_tutorial/Refine3D.html>

### Sigma statistics should be accumulated during pose initialization

The CC pass establishes each cohort particle's pose and state only as that
particle is matched, so a separate fixed-pose sigma pass cannot assume that a
complete trusted pose set existed before pose initialization. After each hard CC
assignment, the matcher should therefore accumulate the corresponding
Euclidean residual contribution during the same particle read. CC does not
consume those sigmas; it only supplies the committed assignment at which the
bootstrap residual is evaluated.

Reference normalization and the common 15 A low-pass limit reduce the initial
scale mismatch. The resulting sigmas are pose initialization bootstrap statistics for
the first Euclidean stage and are subsequently updated by ordinary Euclidean
refinement. This is preferable to the generic particle power-spectrum
bootstrap (`calc_pspec`), which sets sigma2 to total particle power, signal
included, and ignores the initialized model entirely.

### Populated orientations are not evidence of pose initialization

Base `refine3D` randomizes virgin orientations during initialization. After
that operation, the orientation fields are structurally populated even though
no particle has yet been compared with the input references. State labels may
likewise be seeded before matching.

In the single-pass design this concern is handled structurally: the wrapper
owns the pose initialization stage and its capped cohort, so completion of that stage
— not any inherited metadata — is the evidence that the cohort was initialized.
Random Euler angles and random initial state labels are search seeds only.

### Early subset reconstruction changes the question for later particles

A fractional pose initialization could map a subset of particles and then replace the
supplied references with reconstructions from that subset. Particles
encountered later would no longer be classified against the same models as
particles encountered early. The single capped pass with `volrec=no` and
`trail_rec=no` removes this failure mode: the supplied references remain
immutable for the entire pass, and the first data-derived reconstruction uses
the complete assigned pose initialization cohort.

## Implemented SIMPLE Behavior

The working-tree implementation centralizes the policy:

| Path | Implemented policy |
| --- | --- |
| Base `refine3D` | Remains mechanism-only and does not infer reference provenance. |
| `refine3D_auto` with `vol1` | Exposes `ref_pose_init=cc|none`; `none` warns prominently that the external map is trusted without CC initialization. |
| `refine3D_states` | Starts directly with Euclidean matching because its contract requires a same-lineage pose scaffold and references. |
| `classify3D_refs` | Requires a complete external reference set and always performs the shared CC pose-initialization stage. |
| `abinitio3D` with input volumes | Dispatches the shared CC pose-initialization service before entering its Euclidean stage ladder. |

Some startup and final reconstruction commands set `objfun=cc`. Those are
reconstruction-policy settings and do not constitute CC particle-to-reference
matching.

Relevant current sources are:

- `src/main/params/simple_parameters.f90` for the default objective;
- `src/main/strategies/parallelization/simple_refine3D_strategy.f90` for
  random orientation initialization and objective-dependent sigma setup;
- `src/main/strategies/search/simple_strategy3D_matcher.f90` for the
  matcher-side sigma guards and the existing fixed-pose `refine=sigma` mode;
- `src/main/simple_external_reference_pose_initialization.f90` for the shared
  fixed-reference service;
- `src/main/commanders/simple/simple_commanders_refine3D.f90` for the
  `refine3D_auto`, `refine3D_states`, and `classify3D_refs` callers;
- `src/main/commanders/simple/simple_commanders_abinitio.f90` for the
  `abinitio3D` input-volume route;
- the canonical state- and reference-workflow policy documents for wrapper
  and coverage contracts.

## Reference Provenance Is the Application Boundary

### Trusted-pose initialization

Use this path when:

- particle orientations are meaningful;
- the references and orientations come from the same reconstruction lineage;
- the next task is fixed-pose or bounded-local state refinement.

The target `refine3D_states` workflow belongs here. It may start directly with
`objfun=euclid` because the initial residuals are evaluated around an existing
pose scaffold and compatible maps.

Supplying a same-lineage volume through `vol1` does not turn the run into
external-reference pose initialization.

### External-reference pose initialization

Use this path when:

- references may have been produced independently of the particles;
- existing particle poses were not established against those references;
- broad pose and state correspondence is part of the scientific task.

The target `classify3D_refs` workflow belongs here. Its normal contract should
always include the single capped-cohort CC pose initialization stage before Euclidean
refinement.

### Ambiguous callers

Generic applications such as `refine3D_auto` should not infer provenance from
the presence or absence of numeric orientation fields. `refine3D_auto` should
expose the explicit `ref_pose_init=cc|none` policy. Untrusted external
references remain valid inputs to that mode, but the command must emit a
prominent log and metadata warning that their provenance is untrusted and
that CC pose initialization is being used before Euclidean refinement. Base
`refine3D` should not silently guess the policy; its caller has the scientific
context needed to select the phase.

## Pose-Initialization Stage Contract

### 1. Validate and prepare the supplied references

Before matching:

- require a complete `vol1..volN` set;
- validate box size, sampling, symmetry, mask assumptions, handedness where it
  can be checked, and state count;
- normalize references without using particle assignments that do not yet
  exist;
- apply one common 15 A low-pass limit;
- preserve immutable prepared copies for the complete pose initialization stage.

Every state must be evaluated at the same effective bandwidth. Per-state
frequency adaptation would bias competitive assignment and is out of scope.

### 2. Run one capped-cohort CC pass against the fixed references

The wrapper launches a single `refine3D` iteration with:

- `objfun=cc`;
- the existing pose initialization cohort selection, capped at 100,000 particles;
- no fractional-update schedule or follow-up coverage pass;
- broad angular, in-plane, and translational candidate preparation;
- all supplied states evaluated competitively;
- one hard pose and state assignment committed per particle;
- `volrec=no` and `trail_rec=no`, so no partial or trailing reconstruction can
  modify the supplied references;
- Euclidean residual sigma contributions accumulated after each committed CC
  assignment, even though CC itself does not consume them.

The ordinary matcher continues to own candidate scoring and pose/state
updates. The wrapper owns the fixed-reference stage and the objective
transition. Random Euler angles and random initial state labels may be used
only as search seeds.

For very large datasets, the cost of this stage is one pass over at most
100,000 particles at 15 A, read through the downscaled particle cache. The cap
is final for pose initialization: the workflow must not reinterpret it as the first
fraction of a multi-pass coverage scheme.

### 3. Establish a data-derived checkpoint

After the capped pass:

- validate that every particle in the pose initialization cohort has a valid state and
  pose;
- report state populations and reject empty states;
- reconstruct all state maps from the complete assigned pose initialization cohort;
- produce the normal even/odd, FSC, and reference artifacts;
- record the transition as a restartable phase boundary.

This full reconstruction places the working references on SIMPLE's particle
and reconstruction scale. It must not be a trailing blend with the independent
input maps.

### 4. Consolidate Euclidean statistics from the CC pass

The pose initialization matcher writes per-partition sigma contributions calculated
at each particle's newly committed pose/state assignment. After matching:

- consolidate those contributions with the ordinary group-sigma service so
  the sigma star required by the first Euclidean iteration is in place;
- use global residual statistics for this capped-cohort transition, because
  the cohort is not required to cover every acquisition group;
- do not run a separate fixed-pose sigma pass;
- do not run the generic `calc_pspec` power-spectrum bootstrap on this path.

The underlying residual computation (`gen_sigma_contrib` in
`simple_polarft_corr.f90`) is objective-agnostic. The required implementation
work is to permit the CC pose initialization matcher to emit these Euclidean residual
contributions without making CC depend on sigma values, then preserve
iteration-number continuity so the first Euclidean iteration resolves the
consolidated star.

### 5. Switch to Euclidean refinement

- reset iteration/update-epoch fields deliberately for the refinement phase
  while preserving the initialized pose and state checkpoint;
- switch the effective matching objective to `objfun=euclid`;
- begin the shared frequency-marching plan;
- permit ordinary reconstruction, trailing updates, masking, and
  postprocessing from this point onward.

Because the CC pass has already emitted and consolidated sigma contributions,
the first Euclidean stage takes the existing "reuse existing particle sigma
files" initialization branch; no power-spectrum bootstrap is needed there.
The CC-to-Euclidean transition is a wrapper-controlled stage boundary, not a
hidden objective change inside the matcher.

## Workflow Consequences

### `classify3D_refs`

The single capped-cohort CC pose initialization is mandatory for the
independent-reference workflow. There is no compatibility alias or opt-out in
`classify3D_refs`.

### `refine3D_states`

No pose initialization stage is required for valid same-lineage inputs with trusted
orientations. If independently produced references are supplied, the input
violates the workflow contract and should be rejected or redirected to
`classify3D_refs`, not silently treated as same-lineage maps.

### `refine3D_auto`

`refine3D_auto` exposes `ref_pose_init=cc|none` for external-reference starts.
Untrusted references are accepted rather than rejected, but the command makes
their use conspicuous in logs. The warning identifies that the reference
provenance is untrusted, that the reference may bias pose initialization, and
that the run is entering CC pose initialization before Euclidean refinement.

### `abinitio3D`

`abinitio3D` continues to accept input volumes, and that route dispatches the
same shared pose-initialization service rather than owning another implementation.

## Implementation Ownership

| Responsibility | Owner |
| --- | --- |
| Select trusted versus pose initialization policy | Public workflow commander |
| Validate and prepare external references | Independent-reference wrapper using existing volume helpers |
| Select the capped cohort, force the single pose initialization pass, and own the stage boundary | Independent-reference wrapper/refinement strategy |
| Generate candidates and score with CC | Base search and matcher layers |
| Keep pose initialization references immutable | Wrapper command construction and artifact policy (`volrec=no`, `trail_rec=no`) |
| Accumulate Euclidean residual sigma contributions after each CC assignment | Matcher residual path during the pose initialization read |
| Reconstruct the pose initialization-cohort checkpoint | Existing reconstruction and assembly commands |
| Consolidate group sigmas | Existing `calc_group_sigmas` service |
| Run subsequent frequency-marched refinement | Base `refine3D` under the owning wrapper |

No new matching formula, fractional scheme, or pose initialization coverage
bookkeeping is proposed. The change is orchestration: choose the existing CC
objective for one capped-cohort pass against fixed 15 A references, accumulate
residual sigma contributions after assignment, reconstruct, consolidate, and
switch to existing Euclidean refinement.

## Failure Policy

- Do not honor caller-supplied fractional-update settings inside pose initialization;
  select at most the existing 100,000-particle cohort and process it once.
- Do not replace or blend the fixed external references at any point during
  the pose initialization pass.
- Write sigma contributions during the CC pass only after a hard assignment;
  do not make CC consume them, run a second sigma pass, or run `calc_pspec`.
- Do not continue into Euclidean refinement without the checkpoint
  reconstruction and consolidated sigma statistics in place.
- Do not allow empty states at the transition without an explicit state
  recovery policy.
- Do not carry an external input map into trailing reconstruction as if it
  were a full-data accumulator.
- Preserve the prepared input references and the pose initialization checkpoint for
  diagnosis and restart. The restartable boundary is: checkpoint volumes plus
  the consolidated sigma star produced from the CC-pass contributions.

## Validation Criteria

### Static and policy validation

- Every caller selects trusted-pose or external-reference pose initialization behavior
  explicitly.
- The pose initialization command line has `objfun=cc`, `volrec=no`, and
  `trail_rec=no`; the selected cohort is capped at 100,000 and is not advanced
  through fractional passes.
- Every external state is prepared at the same 15 A low-pass limit.
- Sigma contributions are generated during pose initialization only after assignment,
  consolidated once, and read by the first Euclidean iteration through the
  existing reuse branch.
- `refine3D_auto` exposes pose initialization and emits the required warning for
  untrusted external references.
- `abinitio3D` input volumes dispatch the shared pose initialization implementation.
- Shared-memory and distributed execution implement the same scientific
  contract.
- Removed command names have no compatibility route.

### Coverage and artifact validation

- The pose initialization iteration processes every selected cohort particle exactly
  once and never selects more than 100,000 particles.
- The hashes of prepared external references remain unchanged throughout the
  pose initialization pass.
- No partial or trailing reconstruction replaces the references before the
  checkpoint.
- The checkpoint reconstruction uses the complete initialized cohort and produces
  valid per-state even/odd maps and FSC metadata.
- Sigma accumulation does not modify the initialized pose or state assignment.
- Restart before the pose initialization pass and at the checkpoint preserves phase,
  fixed references, assignments, and consolidated sigmas without falsely
  advancing the objective.

### Scientific validation

Compare legacy Euclidean initialization with CC pose initialization using:

- same-lineage references and known good poses, where the new path should not
  be needed;
- independent references rescaled by different multiplicative amplitudes;
- references with different low-pass and masking histories;
- virgin and deliberately perturbed orientations;
- balanced and strongly unequal state populations.

Measure orientation recovery, state-assignment accuracy or consistency, state
collapse, sensitivity to reference scaling, initialized map quality, and the
stability of the first Euclidean stages. Additionally compare sigma
initialization accumulated during CC against the legacy `calc_pspec` bootstrap
on the first Euclidean stages.

Compilation and runtime tests are intentionally deferred to user-side
validation in accordance with repository policy. Linux and BOX results must
not be claimed unless their output is observed.

## Resolved Review Decisions

Cyril approved and refined the pose initialization contract on 2026-08-31:

1. `classify3D_refs` always initializes poses with CC when references and
   particles have independent provenance.
2. Pose initialization uses the existing cohort cap of 100,000 particles, processes
   that cohort once, and does not use a fractional scheme.
3. Exactly one CC pass is performed before checkpoint reconstruction.
4. Euclidean residual sigma contributions are accumulated during the CC pass
   after each pose/state assignment; there is no separate fixed-pose pass and
   no particle-power-spectrum bootstrap.
5. `refine3D_auto` exposes `ref_pose_init=cc|none`. Untrusted external references are
   accepted, with a prominent warning in logs.
6. `abinitio3D` continues to accept input volumes and dispatches the same
   shared pose initialization implementation.
7. Every external state uses a common 15 A low-pass limit during pose initialization.

## Completion Criteria

The design is complete when:

- external-reference and trusted-pose initialization are distinct public
  contracts;
- independent-reference matching uses one fixed-reference CC pass over the
  capped pose initialization cohort of at most 100,000 particles;
- all external states use the same 15 A pose initialization limit;
- the first data-derived reference is reconstructed from the complete assigned
  pose initialization cohort;
- Euclidean residual sigma contributions are accumulated during the CC pass
  and consolidated before the first Euclidean iteration;
- `refine3D_auto` exposes `ref_pose_init=cc|none` with prominent untrusted-reference
  warnings, and `abinitio3D` input volumes reuse the same implementation;
- the objective transition is restartable and visible in logs and metadata;
- user-side scientific validation demonstrates improved robustness without
  degrading trusted same-lineage refinement.
