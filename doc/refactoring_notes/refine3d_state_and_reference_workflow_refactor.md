# Refactor Multi-State 3D Workflows Around Scientific Intent

## Status

Proposal for discussion with Cyril, 2026-08-30. Hans has agreed to the
workflow boundary, naming direction, frequency-marching consolidation,
bounded-pose state sorting, and `abinitio3D` handoff described here.

No implementation changes have been made. This is the single living design,
review, implementation, and validation record for this refactor. Review
findings and implementation decisions should be folded into this document
rather than recorded in parallel specifications.

## Executive Summary

The public names `refine3D_multi` and `refine3D_het` do not distinguish the
scientific jobs performed by the workflows. Both are multi-state workflows and
both address structural heterogeneity. Their durable distinction is the
relationship between particles, input orientations, and reference volumes:

- **Conformational State Refinement** starts from a particle project with an
  existing orientation scaffold. State assignments may change while particle
  poses are fixed or refined in bounded neighborhoods. State references belong
  to the same project lineage and may be initialized from an existing split,
  stochastic labels, or FLEX.
- **Reference-guided 3D Classification** competitively matches a particle set
  against supplied references that may have been produced independently of
  the particles. It does not assume a shared pose or state history.

The proposed public commands are:

| Current command | Proposed command | Display name |
| --- | --- | --- |
| `refine3D_multi` | `refine3D_states` | Conformational State Refinement |
| `refine3D_het` | `classify3D_refs` | Reference-guided 3D Classification |

`refine3D_states` will retain stochastic, FLEX, and existing-state
initialization; add explicit fixed and bounded-local pose policies; and consume
a shared frequency-marching plan. `abinitio3D` will stop owning post-split
multi-state refinement and will hand its split checkpoint to
`refine3D_states`.

`classify3D_refs` will retain independent-reference competitive matching. The
fact that both workflows use base `refine3D`, frequency limits, probabilistic
candidate preparation, hard state assignment, and final reconstruction does
not make them the same application. Those are shared mechanisms below two
different public workflow contracts.

## Motivation

### The current names are not user contracts

`multi` and `het` describe properties of the data rather than what the user is
trying to do. The current UI compounds the ambiguity: the two commands have
effectively identical summaries and help text, with only the display names
changing between “Multi-state” and “Heterogeneous.” A user cannot infer:

- whether existing particle orientations are required;
- whether references must belong to the same project;
- whether orientations remain fixed, move locally, or are searched broadly;
- whether the intended result is conformational state sorting or reference-
  guided classification.

The application names and descriptions should answer those questions before a
user inspects advanced parameters.

### The implementation already contains part of the desired boundary

The current `refine3D_multi` policy is a docked, prior-orientation workflow. It
supports:

- `input_oris_refine`, which updates state and pose from an existing
  orientation scaffold;
- `input_oris_fixed`, which changes state/correlation/accounting fields while
  preserving Euler angles, projection direction, in-plane angle, and shifts;
- state initialization from labels, volumes, stochastic assignment, or
  `flex_pca`;
- a state-initialization stage followed by probabilistic-neighborhood
  refinement;
- final all-particle reconstruction at native sampling.

The current `refine3D_het` policy already forces
`multivol_mode=independent`. It supplies explicit references, initializes or
recovers state labels, runs competitive matching through base `refine3D`,
marches matching frequency in short blocks, and reconstructs final state maps.

The refactor should preserve these useful mechanisms while replacing their
ambiguous public framing.

## Decisions Fixed by This Proposal

| Topic | Decision |
| --- | --- |
| Primary application boundary | Same-project orientation scaffold versus independent particle/reference provenance |
| Same-project workflow | `refine3D_states`, Conformational State Refinement |
| Independent-reference workflow | `classify3D_refs`, Reference-guided 3D Classification |
| State initialization | Retain existing-state, stochastic, explicit-volume, and FLEX routes |
| State sorting pose freedom | Expose fixed and bounded-local policies explicitly |
| Focused state sorting | Treat focus evidence separately from pose freedom |
| Frequency marching | Shared planning mechanism consumed by both workflows where appropriate |
| `abinitio3D` boundary | `abinitio3D` owns the ab initio pose scaffold and split checkpoint; `refine3D_states` owns post-split state refinement |
| Base numerical ownership | Matching remains in search/matcher modules; volume assembly remains in `volassemble` and volume-domain modules |
| Compatibility | Preserve old command names and mode spellings as hidden transitional aliases; never silently reinterpret them |

## Target Public Workflow Contracts

### `refine3D_states`: Conformational State Refinement

#### Scientific purpose

Use an existing particle-orientation scaffold to initialize, separate, and
refine conformational states. Typical entry points are:

- continuation immediately after an `abinitio3D` state split;
- additional early conformational exploration from a consensus reconstruction;
- further subdivision of an already reconstructed state;
- focused state sorting while approximately preserving the established pose;
- refinement of existing multi-state labels and maps.

#### Input contract

The input project must contain meaningful 3D orientations for active
particles. State references must belong to the same particle/project lineage,
whether they are supplied explicitly, recovered from project output, or
constructed during state initialization.

The workflow may start from:

1. existing populated state labels and compatible state volumes;
2. a one-state project plus stochastic state initialization;
3. a one-state project plus `flex_pca` initialization;
4. a complete explicit set of same-lineage state volumes;
5. an `abinitio3D` split checkpoint containing state labels, volumes, pose
   anchors, and sampling/update-epoch metadata.

Partial `vol1..volN` input remains invalid. All active states must be populated
before iterative state refinement begins.

#### Pose policy

Replace the implementation-shaped `input_oris_*` public choice with a
scientific `pose_policy`:

| `pose_policy` | State updates | Pose updates | Intended use |
| --- | --- | --- | --- |
| `fixed` | Allowed | Euler angles, projection direction, in-plane angle, and shifts remain unchanged | Pure state assignment and focused state sorting |
| `local` | Allowed | Bounded around the input pose anchor | Default conformational state refinement |

`fixed` is the target spelling for the current `input_oris_fixed` contract.
`local` is the target state-swapping mode requested in this proposal: each
particle keeps one input-pose anchor, all states are evaluated from the same
bounded geometric neighborhood, and a state change must not create an
unrelated pose solution.

The local policy needs explicit bounds for:

- projection-direction deviation;
- in-plane-angle deviation;
- translational deviation.

The exact parameter spellings and automatic defaults should be fixed during
implementation review. The contract is more important than the names: bounds
are measured from the input anchor, are applied consistently across states,
and are checked after every hard assignment. The initial implementation must
be opt-in until bounded-search parity and failure behavior are validated.

A public `free` pose policy is deliberately out of scope for the first
refactor. Broad matching of particles to independently derived references
belongs to `classify3D_refs`. If later evidence supports a free-pose mode for
same-lineage states, it can be added explicitly without weakening the meaning
of `local`.

#### Focused state sorting

Focused evidence and pose freedom are orthogonal controls:

- `pose_policy` controls which geometries may be tested and committed;
- an optional focus mask controls which image/reference evidence contributes
  to state discrimination.

The focus mask must not silently become the final reconstruction mask, redefine
the global particle orientation, or replace the ordinary per-state automask
and nonuniform-filtering policy. Final authoritative maps remain ordinary
all-particle reconstructions unless a separate, explicit reconstruction option
is introduced.

Focused state sorting must be valid with `pose_policy=fixed` and
`pose_policy=local`.

#### Initialization and stage policy

Retain the current initialization capabilities:

- existing state labels and volumes;
- stochastic state assignment;
- FLEX PCA state labels and volumes;
- explicit same-lineage state volumes.

After initialization, state refinement uses the ordinary base `refine3D`
services. The wrapper owns stage planning, state-overlap convergence policy,
update coverage, and final reconstruction. It must not implement candidate
scoring, particle reconstruction, or volume postprocessing itself.

The target workflow accepts `lpstart` and `lpstop` and consumes the shared
frequency-stage plan described below. Frequency marching replaces the current
fixed two-stage low-pass preparation as the general post-initialization path.

### `classify3D_refs`: Reference-guided 3D Classification

#### Scientific purpose

Competitively match particles against a complete set of supplied 3D references
when the particle set and references may be independent. The workflow is
appropriate for model-based classification, reference transfer, and testing a
particle population against externally derived structural hypotheses.

#### Input contract

- `nstates >= 2` is required.
- A complete, compatible `vol1..volN` reference set is required whenever the
  project does not already contain an unambiguous compatible reference set.
- The workflow must not assume that project poses or state labels were derived
  from these references.
- Reference provenance, sampling, box, symmetry, and scale compatibility must
  be validated explicitly.
- Partial reference input remains invalid.

The workflow may initialize missing state labels and virgin orientations, but
that initialization is preparation for independent competitive matching, not
evidence that meaningful states already exist.

#### Search and output contract

The workflow may use broad/global candidate preparation, probabilistic
neighborhood search, and frequency marching. Candidate preparation remains
probabilistic where the selected base mode is probabilistic; the committed
particle update remains a hard state and pose assignment.

Final state reconstruction and postprocessing remain part of the workflow so
that classification produces inspectable maps. The name emphasizes the user's
primary scientific task rather than implying that supplied references and
particles already form one refinement lineage.

If review concludes that iterative improvement of the supplied references is
the primary public promise rather than particle classification, the alternative
name is `refine3D_refs`, with display name “Independent-reference 3D
Refinement.” No other workflow boundary changes under that naming choice.

## Shared Frequency-Marching Contract

Frequency marching is a mechanism, not an application boundary. The current
frequency-block construction inside `exec_refine3D_het` should move to neutral
ownership and be reusable by `refine3D_states`, `classify3D_refs`, and future
callers that need the same policy.

The shared planner should receive explicit inputs and return a plan; it should
not execute a commander or mutate project state. Its contract includes:

- starting and stopping low-pass limits;
- native box and sampling;
- total iteration budget and block length;
- monotonic Fourier-index progression;
- endpoint clamping against crop/Nyquist limits;
- per-block low-pass, translation, crop, iteration range, and LP-set state;
- a stable global iteration horizon for annealing and regularization.

Initial policy should use one common matching-frequency schedule across all
states. State-specific frequency adaptation would change competitive evidence
between states and is a separate scientific feature requiring its own design
and validation.

Recommended ownership is a small workflow-neutral planning module, for example
`src/main/simple_refine3D_stage_plan.f90`, using the existing low-pass/crop
helpers. Commanders remain responsible for consuming the returned plan and
calling base `refine3D`. The planner must not absorb matcher, reconstruction,
filtering, or scheduler behavior.

The present three-iteration HET block size is an implementation default, not a
required universal constant. The plan should represent block length explicitly
so each public workflow can select or inherit policy without copying loops.

## `abinitio3D` Handoff

### Target boundary

`abinitio3D` remains the workhorse for producing an initial 3D orientation
scaffold. In docked multi-state use it owns:

1. the single-state ab initio stages;
2. preparation of the pre-split pose coverage;
3. construction of a valid split checkpoint;
4. selection of the starting frequency and remaining iteration horizon;
5. dispatch to `refine3D_states`.

`refine3D_states` owns:

1. post-split state initialization validation;
2. state-only or bounded-local state swapping;
3. post-split frequency marching;
4. state-overlap convergence and active-particle update coverage;
5. final native-sampling reconstruction and state products.

The existing `abinitio3D` post-split sequence must not remain as a second
implementation of state refinement after the handoff is accepted.

### Required checkpoint contract

The handoff is more than a list of volumes. The current docked split creates a
new multi-state update epoch and carries a persistent sampling cohort. A valid
checkpoint must preserve or explicitly transfer:

- requested and effective state count;
- populated state labels;
- input-pose anchors;
- state volumes and required half-map/reference artifacts;
- current low-pass/crop/translation plan;
- `sampled` current-round markers;
- `updatecnt` persistent update history;
- persistent cohort membership currently represented by `sampled > 0`;
- the intent currently carried by `sticky_class_sampling`;
- realized update fractions needed by trailing reconstruction;
- global iteration numbering and restart provenance.

`refine3D_states` must consume these fields without resampling the cohort or
resetting the update epoch accidentally. `sampled == max(sampled)` must
continue to identify the current update and `sampled > 0` the persistent
cohort until the cohort policy is deliberately retired.

### Migration sequence

For parity, the first implementation should retain the existing split-
checkpoint construction in `abinitio3D` and replace only the subsequent
post-split loop with a `refine3D_states` call.

After parity is established, review whether the stochastic split construction
can move entirely into `refine3D_states`. That second move is desirable only if
the ab initio-specific cohort preparation can be represented as an explicit
input contract rather than hidden caller state.

## Ownership and Source Map

### UI and public command contract

Update:

- `src/main/ui/simple/simple_ui_refine3D.f90`;
- UI program visibility and parameter-visibility policy records;
- generated NICE/CLI metadata through its owning generator rather than by
  editing generated output;
- algorithm and user-facing workflow documentation.

The two new UI summaries must state their input relationship and pose policy
directly. They must not reuse the same generic “multi-state refinement” help
text.

### Execution routing

Update `src/main/exec/simple_exec_refine3D.f90` with explicit routes for the
new commands. Transitional routes for `refine3D_multi` and `refine3D_het`
should normalize to the legacy-compatible contracts before dispatch.

### Commander orchestration

Refactor `src/main/commanders/simple/simple_commanders_refine3D.f90` so that:

- state-refinement initialization, pose policy, stage planning, update
  coverage, and final reconstruction remain in the state workflow commander;
- independent-reference validation, matching orchestration, update coverage,
  and final reconstruction remain in the reference-classification commander;
- duplicated frequency-plan construction is replaced by the shared planner;
- neither wrapper implements candidate scoring or volume postprocessing.

Internal type and procedure names should follow the public names once the
compatibility routes exist. Avoid preserving `multi` and `het` as the primary
internal vocabulary after the contracts are renamed.

### Base refinement and search

Base `refine3D`, probabilistic table generation, and the matcher continue to
own:

- candidate generation and scoring;
- hard state/pose/shift updates;
- sampled-subset reproduction between probabilistic preparation and matching;
- partition-local partial reconstructions;
- single-read particle batch reuse when reconstruction is active.

Bounded-local pose support belongs in the owning search/candidate modules and
parameter validation, not in the wrapper commander. The commander selects the
policy and passes explicit bounds.

### Volume domain

`volassemble` and volume-domain helpers continue to own:

- reduction and restoration of partial reconstructions;
- even/odd handling and trailing-reconstruction weighting;
- FSC and resolution artifacts;
- state-specific automasking and nonuniform filtering;
- final reference products used by the next stage.

Focused classification evidence must not move these responsibilities into the
particle-domain matcher.

### `abinitio3D`

Update:

- `src/main/simple_abinitio_controller.f90`;
- `src/main/commanders/simple/simple_commanders_abinitio.f90`;
- `doc/policies/abinitio3D_policy.md`;
- related fractional-update and trailing-reconstruction policy notes.

The controller should build and dispatch one explicit handoff command rather
than continuing to emit post-split base-`refine3D` stages itself.

## Compatibility and Migration

### Command aliases

For a transitional period:

- `refine3D_multi` remains accepted as a hidden alias;
- `refine3D_het` remains accepted as a hidden alias;
- logs identify the new canonical command and the compatibility mapping;
- restart manifests created with an old command remain readable;
- old aliases do not appear as separate standard UI applications.

### Mode mapping

Legacy mappings must be explicit:

| Legacy command/mode | Transitional target |
| --- | --- |
| `refine3D_multi multivol_mode=input_oris_fixed` | `refine3D_states pose_policy=fixed` |
| `refine3D_multi multivol_mode=input_oris_refine` | Existing pose-refinement behavior behind a compatibility path until bounded-local parity is demonstrated |
| `refine3D_het` | `classify3D_refs` with independent-reference policy |

Do not silently map `input_oris_refine` to a tighter bounded-local search in
legacy runs. Introduce and validate `pose_policy=local` explicitly, then decide
whether it becomes the new-command default. Old commands must retain old
behavior unless the user opts into the new policy.

### Artifact compatibility

Preserve existing state volume, half-map, FSC, assignment, partial
reconstruction, cache, and final reconstruction names in the first
implementation. Renaming public commands does not justify changing scientific
artifact contracts in the same step.

If execution-directory or restart naming embeds the program name, provide an
explicit old-to-new lookup during migration and test continuation in both
directions supported by policy.

## Implementation Plan

### Stage 1: Freeze public and scientific contracts

1. Review this note with Cyril and settle the canonical name of the
   independent-reference workflow.
2. Fix the `pose_policy=local` bounds and focus-mask evidence contract.
3. Document the exact `abinitio3D` split-checkpoint fields and artifacts.
4. Add no new numerical behavior in this stage.

### Stage 2: Introduce canonical names and compatibility routes

1. Register `refine3D_states` and `classify3D_refs` in the UI and router.
2. Give each command distinct summaries, help, required inputs, and standard
   parameters.
3. Preserve old commands as hidden aliases.
4. Rename commander types/procedures after routes are stable.
5. Update policy and algorithm documentation without changing search behavior.

### Stage 3: Extract the shared frequency-stage planner

1. Freeze current HET frequency-march behavior with planner-level tests.
2. Extract plan construction into neutral ownership.
3. Make `classify3D_refs` consume the extracted plan with artifact and
   iteration parity.
4. Make `refine3D_states` consume a frequency plan beginning at its selected
   `lpstart` and ending at `lpstop`.
5. Keep one common schedule across states.

### Stage 4: Add explicit state pose policies

1. Route `pose_policy=fixed` through the existing fixed-orientation behavior.
2. Add opt-in bounded-local candidate generation around a common input anchor.
3. Apply angular, in-plane, and shift bounds consistently across all states.
4. Add focused state-scoring support without changing final reconstruction
   masking or global pose ownership.
5. Preserve probabilistic preparation followed by hard assignment.

### Stage 5: Replace the `abinitio3D` post-split loop

1. Materialize an explicit split checkpoint using the current docked policy.
2. Dispatch `refine3D_states` with the checkpoint, remaining stage plan, and
   update-epoch metadata.
3. Compare the handoff path with the current in-controller post-split path.
4. Remove the duplicate post-split loop only after equivalence gates pass.
5. Retain the pre-split ab initio and checkpoint-construction behavior.

### Stage 6: Complete migration

1. Update tutorials, examples, policy cross-links, generated UI reviews, and
   algorithm descriptions.
2. Mark old commands deprecated in logs and release notes.
3. Remove aliases only in a separately approved compatibility-breaking
   release.
4. Consider moving stochastic split construction from `abinitio3D` into
   `refine3D_states` only as a reviewed follow-up.

## Validation Plan

### Static and interface validation

- New commands are registered once and route to the intended commanders.
- Old aliases remain routable but are not displayed as standard applications.
- UI summaries distinguish same-lineage state refinement from independent-
  reference classification.
- Every new command-line key is registered and validated before parsing.
- UI defaults and commander defaults agree.
- Generated metadata and policy-review tables are regenerated from source.
- Markdown links, code fences, trailing whitespace, and `git diff --check`
  pass.

### Shared frequency planner

- One iteration produces one valid block.
- Non-divisible iteration counts preserve the final remainder block.
- The schedule is monotonic from `lpstart` toward `lpstop`.
- Fourier-index, crop, and Nyquist clamps match current policy.
- Global iteration numbers remain monotonic across blocks and restarts.
- `classify3D_refs` reproduces the current HET block schedule for identical
  inputs.

### State refinement

- Existing populated state projects skip unnecessary state initialization.
- Stochastic and FLEX initialization produce complete populated state sets.
- Partial `vol1..volN` input is rejected.
- `pose_policy=fixed` changes no Euler angle, projection direction, in-plane
  angle, or shift while allowing state changes.
- `pose_policy=local` never exceeds its angular, in-plane, or shift bounds,
  including after state swaps.
- Every state sees the same candidate geometry around each particle anchor.
- Focused state sorting affects classification evidence but not final map-mask
  policy or the stored global pose outside permitted bounds.
- Every active particle has `updatecnt > 0` before final reconstruction.
- Final maps are reconstructed at native project sampling.

### Independent-reference classification

- A particle project and references with independent provenance are accepted
  when box, sampling, symmetry, and scale contracts are satisfied.
- Incomplete or incompatible reference sets fail before matching.
- Virgin orientations and missing state labels follow the documented
  initialization path.
- Competitive matching updates state and pose without assuming prior
  same-lineage assignments.
- Final state maps and their FSC/postprocessing artifacts are produced.

### `abinitio3D` equivalence gate

From the same pre-split project, random seed, and split-stage configuration,
compare the old docked continuation with the new handoff:

- split cohort membership;
- `sampled` and `updatecnt` evolution;
- state populations after label permutation matching;
- pose and shift deltas;
- per-stage update fractions;
- frequency and iteration schedule;
- state volumes, half maps, FSCs, and resolution metadata;
- trailing-reconstruction behavior;
- final project fields and native-sampling outputs;
- restart and continuation artifacts.

State labels are exchangeable, so comparisons must match states by map or
assignment similarity before reporting differences.

### Execution parity and performance

- Shared-memory and distributed routes implement the same scientific
  contracts.
- Probabilistic preparation, matcher execution, and reconstruction reproduce
  the same sampled subset.
- Matching plus reconstruction retains the single particle-stack read per
  active batch.
- The refactor does not introduce a second reconstruction pass inside the
  matcher or duplicate volume assembly in wrappers.
- Peak memory and wall time remain within an agreed tolerance of the existing
  workflows for equivalent policy.

Compilation and runtime tests are intentionally deferred to user-side
validation in accordance with repository policy. No Linux or BOX result may be
claimed without observed output.

## Risks and Review Focus

### Cohort and update-epoch loss at the handoff

This is the highest-risk item. Clearing or reinterpreting `sampled`,
`updatecnt`, or sticky cohort eligibility would change which particles are
seen after the split and how trailing reconstructions are weighted. The
checkpoint must make this policy explicit.

### Pose drift masquerading as state separation

If each state searches a different broad pose neighborhood, the workflow may
explain orientation error as conformational variation. `pose_policy=local`
therefore requires one common particle anchor and explicit post-update bounds.

### Focus-mask leakage

A focus mask used for state evidence must not silently become a reference,
FSC, automask, NU-filter, or final reconstruction mask. Those policies have
different scientific meanings and owners.

### Frequency-dependent state bias

Allowing different states to match at different effective bandwidths can bias
competitive assignment toward the reference with more permissive evidence.
The first implementation uses a common schedule across states.

### Rename without migration

Scripts, restart manifests, UI policy tables, and documentation currently name
`refine3D_multi` and `refine3D_het`. The aliases and artifact-preserving first
stage are required to keep the naming improvement from becoming an avoidable
workflow break.

### Classification versus refinement terminology

`classify3D_refs` still reconstructs and postprocesses maps. Review should
confirm that classification is the primary user intent. If users select the
application principally to iteratively improve externally supplied references,
`refine3D_refs` is the more honest name.

## Questions for Cyril's Review

1. Is **Reference-guided 3D Classification** the right primary user intent, or
   should the independent-reference command be named `refine3D_refs`?
2. Should `pose_policy=local` become the default for the new
   `refine3D_states` command after validation, while legacy aliases retain the
   old behavior?
3. Should local angular, in-plane, and shift bounds be explicit user values,
   automatically derived from sampling, or automatic with advanced overrides?
4. Is the proposed focus-mask contract correct: classification evidence only,
   with authoritative reconstruction and global pose policy remaining
   separate?
5. Is one common frequency schedule across states sufficient initially, with
   any state-specific adaptation deferred?
6. For the first handoff, should `abinitio3D` retain construction of the split
   checkpoint exactly as today, with only post-split refinement delegated?

## Completion Criteria

This refactor is complete when:

- users see two purpose-specific applications with distinct input and pose
  contracts;
- state refinement supports existing, stochastic, and FLEX initialization;
- fixed and bounded-local state sorting are explicit and validated;
- focus evidence is independent of pose and reconstruction masking policy;
- both workflows consume shared frequency-stage planning rather than copied
  loops;
- `abinitio3D` delegates post-split work to `refine3D_states` through a tested
  checkpoint contract;
- old aliases preserve legacy behavior through the agreed migration window;
- source, policies, algorithm descriptions, UI metadata, and restart behavior
  agree with the new names;
- user-side build and runtime validation, including the `abinitio3D`
  equivalence matrix, have been observed and recorded here.
