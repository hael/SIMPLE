# Acquire Poses Against External References Before Euclidean Refinement

## Status

Proposal for discussion with Cyril, 2026-08-30.

No implementation changes have been made. This note is a focused companion to
`refine3d_state_and_reference_workflow_refactor.md`. That document owns the
application boundaries and proposed names; this document owns only the policy
for introducing reference volumes that may not share the particles' alignment
history.

## Proposed Decision

When reference volumes are independent of the particle alignment, SIMPLE
should not begin with noise-normalized Euclidean matching. It should first run
a reference-acquisition stage with cross-correlation, implemented as a single
full pass over the active particle set:

1. normalize, validate, and low-pass filter the supplied references;
2. keep those references fixed;
3. match all active particles competitively with `objfun=cc` in one pass with
   the update fraction forced to 1.0 (`volrec=no`, `trail_rec=no`, no sigma
   bookkeeping);
4. reconstruct complete state maps from the acquired assignments;
5. initialize the Euclidean noise statistics with a dedicated fixed-pose sigma
   pass against those data-derived maps;
6. switch to `objfun=euclid` and start the ordinary frequency-marched
   refinement.

Because the acquisition stage processes every active particle in one pass,
complete coverage holds by construction:

```text
one full-pass CC iteration  =>  every active particle matched once
                                against the fixed reference set
```

No separate acquisition coverage counter, gate, or missing-particle sweep is
required. This is the SIMPLE analogue of RELION's first-iteration
cross-correlation policy, and — like RELION — it establishes the invariant by
processing the complete particle set rather than by auditing fractional
updates.

References from the same project lineage, accompanied by a trusted particle
pose scaffold, may start directly with Euclidean refinement. The policy is
therefore based on reference/pose provenance, not merely on whether `vol1` was
supplied.

## Why This Is Needed

### External references are not initially on the data model's scale

An independently produced map may differ from the particles' working model in
normalization, amplitude scale, filtering, masking, and reconstruction
convention. A Euclidean residual interprets those differences as evidence.
Cross-correlation is the safer acquisition objective because the initial task
is to find pose and state correspondence before estimating meaningful residual
noise statistics.

RELION exposes the same distinction through `do_firstiter_cc`: its source
describes the first iteration as using the signal cross-product instead of the
Gaussian objective. Its refinement tutorial also treats first-iteration
greyscale correction and initial low-pass filtering as explicit reference
initialization concerns:

- <https://github.com/3dem/relion/blob/master/src/ml_optimiser.h>
- <https://relion.readthedocs.io/en/release-5.0/SPA_tutorial/Refine3D.html>

### Sigma statistics must come from maps on the data scale

The same scale argument applies to noise estimation. Sigma values computed
from residuals against the external references would inherit the very scale
mismatch that motivates CC acquisition. The Euclidean noise statistics should
therefore be estimated only after the checkpoint reconstruction, from
residuals between particles and the data-derived maps at the acquired poses.
This replaces both alternatives considered earlier:

- computing sigmas simultaneously with the CC pass (contaminated by the
  external references' scale);
- the generic particle power-spectrum bootstrap (`calc_pspec`), which sets
  sigma2 to total particle power, signal included, and ignores the acquired
  model entirely.

### Populated orientations are not evidence of acquisition

Base `refine3D` randomizes virgin orientations during initialization. After
that operation, the orientation fields are structurally populated even though
no particle has yet been compared with the input references. State labels may
likewise be seeded before matching.

In the single-full-pass design this concern is handled structurally: the
wrapper owns the acquisition stage and forces a full pass, so completion of
that stage — not any inherited metadata — is the evidence of acquisition.
Random Euler angles and random initial state labels are search seeds only.

### Early subset reconstruction changes the question for later particles

A fractional acquisition could map a subset of particles and then replace the
supplied references with reconstructions from that subset. Particles
encountered later would no longer be classified against the same models as
particles encountered early. The single full pass with `volrec=no` and
`trail_rec=no` removes this failure mode: the supplied references remain
immutable for the entire pass, and the first data-derived reconstruction uses
the complete assigned population.

## Current SIMPLE Behavior

The current policy is fragmented:

| Path | Current effective policy |
| --- | --- |
| Base `refine3D` | Defaults to `objfun=euclid`; virgin orientations are randomized; there is no objective transition tied to acquisition. |
| `refine3D_auto` with `vol1` | Normalizes the input volume but hard-sets Euclidean refinement. It does not perform an external-reference acquisition pass. |
| `refine3D_multi` | Hard-sets Euclidean matching. This is consistent with its intended same-lineage, prior-orientation contract. |
| `refine3D_het` | Hard-sets Euclidean matching. A fractional/trailing start may map only an initial subset before replacing the supplied references. |
| `abinitio3D` with input volumes | Uses Euclidean matching from the start. The current external-volume control flow also assumes prior pose information rather than providing a clean pose-acquisition route. |

Some startup and final reconstruction commands set `objfun=cc`. Those are
reconstruction-policy settings and do not constitute CC particle-to-reference
matching.

Relevant current sources are:

- `src/main/params/simple_parameters.f90` for the default objective;
- `src/main/strategies/parallelization/simple_refine3D_strategy.f90` for
  random orientation initialization and objective-dependent sigma setup;
- `src/main/strategies/search/simple_strategy3D_matcher.f90` for the
  matcher-side sigma guards and the existing fixed-pose `refine=sigma` mode;
- `src/main/commanders/simple/simple_commanders_refine3D.f90` for the
  `refine3D_auto`, `refine3D_multi`, and `refine3D_het` policies;
- `src/main/commanders/simple/simple_commanders_abinitio.f90` for the
  `abinitio3D` input-volume route;
- `doc/policies/refine3D_het_policy.md` and
  `doc/policies/refine3D_multi_policy.md` for wrapper and coverage contracts.

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
external-reference acquisition.

### External-reference pose acquisition

Use this path when:

- references may have been produced independently of the particles;
- existing particle poses were not established against those references;
- broad pose and state correspondence is part of the scientific task.

The target `classify3D_refs` workflow belongs here. Its normal contract should
include the single-full-pass CC acquisition stage before Euclidean refinement.

### Ambiguous callers

Generic applications such as `refine3D_auto` should not infer provenance from
the presence or absence of numeric orientation fields. They should either:

1. require an explicit initialization policy such as `trusted` or
   `acquire_cc`; or
2. reject an external-reference start when the required pose provenance is not
   expressed by the workflow contract.

The exact parameter spelling remains for review. Base `refine3D` should not
silently guess the policy; its caller has the scientific context needed to
select the phase.

## Acquisition-Stage Contract

### 1. Validate and prepare the supplied references

Before matching:

- require a complete `vol1..volN` set;
- validate box size, sampling, symmetry, mask assumptions, handedness where it
  can be checked, and state count;
- normalize references without using particle assignments that do not yet
  exist;
- apply one common conservative low-pass limit;
- preserve immutable prepared copies for the complete acquisition stage.

Every state must be evaluated at the same effective bandwidth. Per-state
frequency adaptation would bias competitive assignment and is out of scope.

### 2. Run one full-pass CC iteration against the fixed references

The wrapper launches a single `refine3D` iteration with:

- `objfun=cc`;
- the update fraction forced to 1.0 for this stage, stripping or overriding
  any fractional-update or sampling-size settings supplied by the caller;
- broad angular, in-plane, and translational candidate preparation;
- all supplied states evaluated competitively;
- one hard pose and state assignment committed per particle;
- `volrec=no` and `trail_rec=no`, so no partial or trailing reconstruction can
  modify the supplied references;
- no sigma computation or sigma file output, because CC does not consume noise
  statistics and residuals against external references must not seed them.

The ordinary matcher continues to own candidate scoring and pose/state
updates. The wrapper owns the fixed-reference stage and the objective
transition. Random Euler angles and random initial state labels may be used
only as search seeds.

For very large datasets, the cost of this stage is one pass over the active
particle set at a conservative low-pass limit, read through the downscaled
particle cache. If a fractional acquisition variant is ever required, it can
be built later on the existing `update_missing` machinery; it is deliberately
out of scope for the first implementation.

### 3. Establish a data-derived checkpoint

After the full pass:

- validate that every active particle has a valid state and pose;
- report state populations and reject empty states;
- reconstruct all state maps from all active assigned particles;
- produce the normal even/odd, FSC, and reference artifacts;
- record the transition as a restartable phase boundary.

This full reconstruction places the working references on SIMPLE's particle
and reconstruction scale. It must not be a trailing blend with the independent
input maps.

### 4. Initialize Euclidean statistics with a fixed-pose sigma pass

Using the acquired poses and the complete data-derived references:

- run one fixed-pose sigma pass (the existing `refine=sigma` matcher mode) in
  which every active particle's residual against the checkpoint maps is
  evaluated at its acquired pose, without searching;
- write the per-partition sigma files and consolidate them with the ordinary
  group-sigma step, so the star file the first Euclidean iteration reads is in
  place;
- do not run the generic `calc_pspec` power-spectrum bootstrap on this path;
  the sigma pass supersedes it.

This is one extra read pass over the particle set. It buys sigma statistics
that are estimated on the data scale, against the model the Euclidean stages
will actually refine — which is the point of deferring them.

Implementation notes: the fixed-pose sigma path exists in
`simple_strategy3D_matcher.f90` (`refine=sigma`), and the underlying residual
computation (`gen_sigma_contrib` in `simple_polarft_corr.f90`) is
objective-agnostic. What is missing is a fresh-initialization branch in the
sigma preparation (`prep_sigmas_objfun`), which currently assumes existing
sigma files, and iteration-number continuity so that the sigma star written by
this pass is the one the first Euclidean iteration resolves.

### 5. Switch to Euclidean refinement

- reset iteration/update-epoch fields deliberately for the refinement phase
  while preserving the acquired pose and state checkpoint;
- switch the effective matching objective to `objfun=euclid`;
- begin the shared frequency-marching plan;
- permit ordinary reconstruction, trailing updates, masking, and
  postprocessing from this point onward.

Because the sigma files already exist, the first Euclidean stage takes the
existing "reuse existing particle sigma files" initialization branch; no new
bootstrap logic is needed there. The CC-to-Euclidean transition is a
wrapper-controlled stage boundary, not a hidden objective change inside the
matcher.

## Workflow Consequences

### `classify3D_refs`

The single-full-pass CC acquisition should become part of the target
scientific contract for the new independent-reference workflow. The legacy
`refine3D_het` alias must retain its current behavior during migration unless
the user explicitly opts into the new acquisition policy.

Initial implementation must therefore be protected by an explicit default-off
key. After numerical and workflow validation, the canonical new command may
make the acquisition contract mandatory while the legacy alias remains
behavior-compatible for the agreed transition period.

### `refine3D_states`

No acquisition stage is required for valid same-lineage inputs with trusted
orientations. If independently produced references are supplied, the input
violates the workflow contract and should be rejected or redirected to
`classify3D_refs`, not silently treated as same-lineage maps.

### `refine3D_auto`

The command needs an explicit decision:

- retain a strict refinement contract and require trusted orientations with an
  external starting volume; or
- expose the same opt-in CC acquisition stage for untrusted starts.

The first option gives the application a clearer boundary. The second is more
flexible but makes provenance and search policy part of its public interface.

### `abinitio3D`

Input volumes should not create another implementation of external-reference
acquisition. If that use case remains supported, `abinitio3D` should dispatch
the shared acquisition service or hand off to `classify3D_refs`. The existing
control flow should be audited separately before it is used as a foundation
for the new policy.

## Implementation Ownership

| Responsibility | Owner |
| --- | --- |
| Select trusted versus acquisition policy | Public workflow commander |
| Validate and prepare external references | Independent-reference wrapper using existing volume helpers |
| Force the full-pass acquisition iteration and own the stage boundary | Independent-reference wrapper/refinement strategy |
| Generate candidates and score with CC | Base search and matcher layers |
| Keep acquisition references immutable | Wrapper command construction and artifact policy (`volrec=no`, `trail_rec=no`) |
| Reconstruct the complete acquisition checkpoint | Existing reconstruction and assembly commands |
| Initialize Euclidean sigma statistics at fixed poses | Existing `refine=sigma` matcher mode plus a fresh-initialization branch in sigma preparation |
| Consolidate group sigmas | Existing `calc_group_sigmas` service |
| Run subsequent frequency-marched refinement | Base `refine3D` under the owning wrapper |

No new matching formula and no new coverage bookkeeping are proposed. The
change is orchestration: choose the existing CC objective for one full-pass
stage against fixed references, reconstruct, estimate sigmas at fixed poses on
the data scale, and switch to existing Euclidean refinement.

## Compatibility and Failure Policy

- Do not change the default behavior of existing public commands without an
  explicit opt-in during initial implementation.
- Do not honor caller-supplied fractional-update settings inside the
  acquisition stage; the stage is a full pass by contract.
- Do not replace or blend the fixed external references at any point during
  the acquisition pass.
- Do not write sigma files during the CC pass, and do not run `calc_pspec` on
  this path; sigma initialization belongs to the fixed-pose pass against the
  checkpoint maps.
- Do not continue into Euclidean refinement without the checkpoint
  reconstruction and consolidated sigma statistics in place.
- Do not allow empty states at the transition without an explicit state
  recovery policy.
- Do not carry an external input map into trailing reconstruction as if it
  were a full-data accumulator.
- Preserve the prepared input references and the acquisition checkpoint for
  diagnosis and restart. The restartable boundary is: checkpoint volumes plus
  the consolidated sigma star.

## Validation Criteria

### Static and policy validation

- Every caller selects trusted-pose or external-acquisition behavior
  explicitly.
- The acquisition command line has `objfun=cc`, `volrec=no`, and
  `trail_rec=no`, and its effective update fraction is 1.0 regardless of
  caller-supplied sampling settings.
- No sigma files are produced during acquisition; sigma files appear only from
  the fixed-pose pass, and the first Euclidean iteration reads them via the
  existing reuse branch rather than bootstrapping.
- Shared-memory and distributed execution implement the same scientific
  contract.
- Existing aliases retain legacy behavior when the opt-in key is disabled.

### Coverage and artifact validation

- The acquisition iteration processes every active particle exactly once.
- The hashes of prepared external references remain unchanged throughout the
  acquisition pass.
- No partial or trailing reconstruction replaces the references before the
  checkpoint.
- The checkpoint reconstruction uses all active assigned particles and
  produces valid per-state even/odd maps and FSC metadata.
- The fixed-pose sigma pass does not modify any pose or state assignment.
- Restart before the acquisition pass, at the checkpoint, and after the sigma
  pass preserves phase, fixed references, and assignments without falsely
  advancing the objective.

### Scientific validation

Compare legacy Euclidean initialization with CC acquisition using:

- same-lineage references and known good poses, where the new path should not
  be needed;
- independent references rescaled by different multiplicative amplitudes;
- references with different low-pass and masking histories;
- virgin and deliberately perturbed orientations;
- balanced and strongly unequal state populations.

Measure orientation recovery, state-assignment accuracy or consistency, state
collapse, sensitivity to reference scaling, acquired map quality, and the
stability of the first Euclidean stages. Additionally compare the fixed-pose
sigma initialization against the legacy `calc_pspec` bootstrap on the first
Euclidean stages.

Compilation and runtime tests are intentionally deferred to user-side
validation in accordance with repository policy. Linux and BOX results must
not be claimed unless their output is observed.

## Questions for Cyril

1. Do you agree that `classify3D_refs` should always acquire poses with CC when
   references and particles have independent provenance?
2. Do you agree the acquisition stage should be one full pass over all active
   particles with the references completely fixed, rather than a fractional
   scheme with coverage bookkeeping?
3. Is a single CC pass sufficient before the checkpoint reconstruction, or is
   there a scientific case for iterating the fixed-reference matching more
   than once?
4. Do you agree that Euclidean sigmas should be initialized from residuals
   against the data-derived checkpoint maps at fixed poses, rather than during
   the CC pass or from particle power spectra?
5. Should `refine3D_auto` reject untrusted external-reference starts, or
   expose an explicit acquisition mode?
6. Should `abinitio3D` continue accepting input volumes, and if so, should
   that route dispatch this same acquisition stage rather than own another
   version?
7. What conservative initial low-pass rule should be shared by all external
   states during acquisition?

## Completion Criteria

The proposal is complete when:

- external-reference and trusted-pose initialization are distinct public
  contracts;
- independent-reference matching uses one full-pass fixed-reference CC stage;
- the first data-derived reference is reconstructed from the complete assigned
  population;
- Euclidean noise statistics are initialized by a fixed-pose sigma pass
  against that checkpoint, on the data scale;
- the objective transition is restartable and visible in logs and metadata;
- legacy behavior remains available during the agreed migration period;
- user-side scientific validation demonstrates improved robustness without
  degrading trusted same-lineage refinement.
