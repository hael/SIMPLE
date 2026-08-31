# Refactor Multi-State 3D Workflows Around Scientific Intent

## Status

Implementation in progress, 2026-08-31. Cyril and Hans agreed to the workflow
boundary, canonical names, frequency-marching consolidation, three scientific
pose policies, focus-evidence boundary, and `abinitio3D` handoff described
here. The canonical commands, pose policies, shared frequency planner,
external-reference pose initialization path, reusable split-checkpoint builder, and
`abinitio3D` handoff are implemented in the working tree. Compilation and
runtime validation remain user-owned. This is the single living design,
implementation, review, and validation record for the refactor.

Static review fixes applied 2026-08-31:

- `classify3D_refs` pose initialization passed the full fixed-size
  `params%vols(99)` where an `nstates`-sized checkpoint array was required; it
  now uses a local `nstates`-sized array with copy-back (this also removes an
  intent(out) aliasing of `params%vols` against `params`).
- `refine3D_states` input-volume validation accepts any cubic downscaled
  sampling covering the native physical extent (base `refine3D` rescales
  references to the stage crop), so the `abinitio3D` split checkpoint volumes
  reconstructed at the abinitio ladder crop pass validation.
- The split-checkpoint routine now assigns the module-level `update_frac`
  explicitly before emitting the split-stage cline, and the caller receives the
  post-split fraction through a local rather than aliasing the module variable
  with an intent(out) dummy.
- The `refine3D_states` frequency march passes a real per-block `minits`
  (1 once the march-wide minimum is satisfied) instead of `minits == maxits`,
  so blocks can stop early on the state-overlap target while the march still
  completes the low-pass schedule.
- The abinitio3D handoff deletes the inherited `overlap` key: `refine3D_states`
  owns the state-overlap convergence policy.
- `local_*_bound` validation rejects negative values other than the `-1.`
  automatic sentinel; `cc_emit_sigma=yes` with `objfun=euclid` throws.
- The downscaled particle cache is now a 2D-only feature: `cache`/`cache_dir`
  removed from all 3D UI programs, cache consumption removed from the 3D
  matcher/reconstruction/prob paths, and the refine3D strategies and
  `abinitio3D` throw on `cache=yes`.
- `refine3D_states` and `classify3D_refs` no longer expose `autoscale`: staged
  frequency marching always uses planner-driven downscaling (a user-supplied
  `autoscale` is ignored with a notice). `refine3D`/`refine3D_auto` and the 2D
  workflows keep their existing autoscale behavior.

Known open parity deviations to weigh at the equivalence gate: the
`classify3D_refs` planner varies crop and `trslim` per block where the old HET
schedule held them constant; the handoff skips the old split-stage
`prob_state` block (pose_policy=local starts the geom march directly); and
`refine3D_states` rebalances sampling with projection-direction bins where the
old post-split loop used 2D-class bins.

## Executive Summary

The removed names `refine3D_multi` and `refine3D_het` did not distinguish the
scientific jobs performed by the workflows. Both are multi-state workflows and
both address structural heterogeneity. Their durable distinction is the
relationship between particles, input orientations, and reference volumes:

- **Conformational State Refinement** starts from a particle project with an
  existing orientation scaffold. State assignments may change while pose
  search is restricted to in-plane degrees of freedom, constrained to a local
  geometric neighborhood, or allowed to match globally. State references
  belong to the same project lineage and may be initialized from an existing
  split, stochastic labels, or FLEX.
- **Reference-guided 3D Classification** competitively matches a particle set
  against supplied references that may have been produced independently of
  the particles. It does not assume a shared pose or state history.

The approved public commands are:

| Removed implementation-shaped command | Canonical command | Display name |
| --- | --- | --- |
| `refine3D_multi` | `refine3D_states` | Conformational State Refinement |
| `refine3D_het` | `classify3D_refs` | Reference-guided 3D Classification |

`refine3D_states` will retain stochastic, FLEX, and existing-state
initialization; expose explicit `fixed`, `local`, and `global` pose policies;
and consume a shared frequency-marching plan. `abinitio3D` will stop owning
post-split multi-state refinement and will hand its split checkpoint to
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

The former `refine3D_multi` policy was a docked, prior-orientation workflow. It
supports:

- `input_oris_refine`, which updates state and pose from an existing
  orientation scaffold;
- `input_oris_fixed`, which keeps the existing projection direction, optimizes
  in-plane rotation and translation for each state candidate, and samples the
  state assignment stochastically;
- state initialization from labels, volumes, stochastic assignment, or
  `flex_pca`;
- a state-initialization stage followed by probabilistic-neighborhood
  refinement;
- final all-particle reconstruction at native sampling.

The former `refine3D_het` policy forced
`multivol_mode=independent`. It supplies explicit references, initializes or
recovers state labels, runs competitive matching through base `refine3D`,
marches matching frequency in short blocks, and reconstructs final state maps.

The refactor should preserve these useful mechanisms while replacing their
ambiguous public framing.

## Approved Design Decisions

| Topic | Decision |
| --- | --- |
| Primary application boundary | Same-project orientation scaffold versus independent particle/reference provenance |
| Same-project workflow | `refine3D_states`, Conformational State Refinement |
| Independent-reference workflow | `classify3D_refs`, Reference-guided 3D Classification |
| State initialization | Retain existing-state, stochastic, explicit-volume, and FLEX routes |
| State sorting pose freedom | Expose `fixed`, `local`, and `global` policies explicitly; default the new command to `global` |
| Local search limits | Derive angular, in-plane, and shift bounds automatically from sampling; permit explicit user overrides |
| Focused state sorting | Treat focus evidence separately from pose freedom |
| Frequency marching | Shared planning mechanism consumed by both workflows where appropriate |
| `abinitio3D` boundary | `abinitio3D` owns the ab initio pose scaffold and split checkpoint; `refine3D_states` owns post-split state refinement |
| Base numerical ownership | Matching remains in search/matcher modules; volume assembly remains in `volassemble` and volume-domain modules |
| Migration policy | Make a clean break: expose and route only the canonical commands; test environments adopt the new contract |

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
scientific `pose_policy` whose three values describe the permitted search:

| `pose_policy` | State search | Pose search | Intended use | Base search mapping |
| --- | --- | --- | --- | --- |
| `fixed` | Stochastic search over states | Keep the projection direction fixed; optimize the in-plane angle and x/y translations | Pose-anchored state sorting and focused classification | State-only probabilistic preparation (`refine=prob_state`) with in-plane optimization enabled |
| `local` | Search across states | Optimize projection direction within the geometric neighborhood plus in-plane degrees of freedom | Neighborhood-constrained state refinement | `refine=prob_neigh`, `prob_neigh_mode=geom` |
| `global` | Search across states | Optimize all pose degrees of freedom through full probabilistic matching over state-pooled candidate neighborhoods | Wide conformational exploration and recovery from an uncertain pose scaffold | `refine=prob_neigh`, `prob_neigh_mode=state` |

`fixed` names the state-only probability-table behavior directly. The path takes each particle's
existing projection direction, evaluates every active state at that projection,
optimizes in-plane rotation and translation for the state candidates, and
samples the state choice stochastically. The rename introduces no new
numerical behavior for this policy.

`local` is the existing geometric-neighborhood scientific policy. Every active
state is evaluated in the subspace containing the particle's current
projection, with no coarse state-pooled peak search. The refactor should reuse
`prob_neigh_mode=geom`. Angular, in-plane, and shift bounds are derived
automatically from the active sampling and search geometry, with advanced
user overrides for deliberate departures from the automatic policy. The
overrides configure the existing geometric-neighborhood search; they do not
create a parallel neighborhood implementation.

`global` is the full probabilistic-matching policy. Coarse representatives are
scored independently per state, selected neighborhoods are pooled, and the
same pooled projection search space is evaluated for every active state. It
optimizes all pose degrees of freedom, maps to `prob_neigh_mode=state`, and is
the agreed default for the new
`refine3D_states` command. Here `global` describes pose-search breadth; it does
not imply independent reference provenance. Independent references still
belong to `classify3D_refs`.

#### Focused state sorting

Focused evidence and pose freedom are orthogonal controls:

- `pose_policy` controls which geometries may be tested and committed;
- an optional focus mask controls which image/reference evidence contributes
  to state discrimination.

The focus mask must not silently become the final reconstruction mask, redefine
the stored particle pose, or replace the ordinary per-state automask
and nonuniform-filtering policy. Final authoritative maps remain ordinary
all-particle reconstructions unless a separate, explicit reconstruction option
is introduced.

Focused state sorting must be valid with all three pose policies. The focus
mask changes the evidence used to discriminate states, not the meaning of
`fixed`, `local`, or `global`.

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

The approved public name is `classify3D_refs`: classification is the primary
user intent even though the workflow reconstructs and postprocesses the maps
needed to inspect its hard assignments.

## Shared Frequency-Marching Contract

Frequency marching is a mechanism, not an application boundary. The current
frequency-block construction is now owned by neutral
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
3. construction of a valid split checkpoint through a separate reusable
   routine;
4. selection of the starting frequency and remaining iteration horizon;
5. dispatch to `refine3D_states`.

`refine3D_states` owns:

1. post-split state initialization validation;
2. state refinement under the selected `fixed`, `local`, or `global` pose
   policy;
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

Checkpoint construction must be extracted as a separate reusable routine with
explicit inputs and outputs. `abinitio3D` remains its first caller and retains
the current split policy, but the routine must not depend on hidden controller
state so future workflows can construct the same validated handoff.

### Migration sequence

For parity, the first implementation should extract the existing split-
checkpoint construction unchanged into the reusable routine, call it from
`abinitio3D`, and replace only the subsequent post-split loop with a
`refine3D_states` call.

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

`src/main/exec/simple_exec_refine3D.f90` routes only the canonical commands.
There are no compatibility aliases or alternate defaults.

### Commander orchestration

Refactor `src/main/commanders/simple/simple_commanders_refine3D.f90` so that:

- state-refinement initialization, pose policy, stage planning, update
  coverage, and final reconstruction remain in the state workflow commander;
- independent-reference validation, matching orchestration, update coverage,
  and final reconstruction remain in the reference-classification commander;
- duplicated frequency-plan construction is replaced by the shared planner;
- neither wrapper implements candidate scoring or volume postprocessing.

Internal type and procedure names follow the canonical public names; `multi`
and `het` are not retained as production workflow vocabulary.

### Base refinement and search

Base `refine3D`, probabilistic table generation, and the matcher continue to
own:

- candidate generation and scoring;
- hard state/pose/shift updates;
- sampled-subset reproduction between probabilistic preparation and matching;
- partition-local partial reconstructions;
- single-read particle batch reuse when reconstruction is active.

Pose-policy mapping belongs in the owning search/candidate modules and
parameter validation, not in the wrapper commander. The commander selects the
policy; search ownership interprets `fixed` as state-plus-in-plane search,
`local` as `prob_neigh_mode=geom`, and `global` as
`prob_neigh_mode=state`.

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

The controller should call the reusable split-checkpoint routine and dispatch
one explicit handoff command rather than continuing to emit post-split
base-`refine3D` stages itself.

## Clean-Break Migration

The refactor intentionally does not preserve old command routes, hidden
aliases, historical defaults, or public `multivol_mode`/
`prob_neigh_mode` combinations. Test environments and scripts must adopt:

- `refine3D_states pose_policy=fixed|local|global`;
- `classify3D_refs` for external-reference classification.

`pose_policy=fixed` keeps the projection direction and commits optimized
in-plane angle and translations together with the stochastic state choice. It
does not freeze the full stored pose record.

### Artifact compatibility

Preserve existing state volume, half-map, FSC, assignment, partial
reconstruction, cache, and final reconstruction names in the first
implementation. Renaming public commands does not justify changing scientific
artifact contracts in the same step.

Restart and execution-directory validation targets the canonical names only.

## Implementation Plan

### Stage 1: Freeze public and scientific contracts

1. Record the approved canonical name `classify3D_refs` and the three-value
   `pose_policy` contract.
2. Record `global` as the command default and automatic, user-overridable
   local bounds.
3. Document the exact `abinitio3D` split-checkpoint fields and artifacts.
4. Add no new numerical behavior in this stage.

### Stage 2: Introduce canonical names

1. Register `refine3D_states` and `classify3D_refs` in the UI and router.
2. Give each command distinct summaries, help, required inputs, and standard
   parameters.
3. Remove the old public routes.
4. Rename commander types/procedures to the canonical vocabulary.
5. Update policy and algorithm documentation.

### Stage 3: Extract the shared frequency-stage planner

1. Freeze current HET frequency-march behavior with planner-level tests.
2. Extract plan construction into neutral ownership.
3. Make `classify3D_refs` consume the extracted plan with artifact and
   iteration parity.
4. Make `refine3D_states` consume a frequency plan beginning at its selected
   `lpstart` and ending at `lpstop`.
5. Keep one common schedule across states.

### Stage 4: Add explicit state pose policies

1. Expose the existing `input_oris_fixed` behavior as `pose_policy=fixed`:
   stochastic state search with projection direction fixed and only in-plane
   degrees of freedom optimized.
2. Implement `pose_policy=local` through `prob_neigh_mode=geom`, deriving its
   angular, in-plane, and shift bounds automatically while accepting explicit
   user overrides.
3. Implement `pose_policy=global` through `prob_neigh_mode=state` and make it
   the new-command default.
4. Add focused state-scoring support without changing final reconstruction
   masking or stored-pose ownership.
5. Preserve probabilistic preparation followed by hard assignment.

### Stage 5: Replace the `abinitio3D` post-split loop

1. Extract the current docked split-checkpoint construction into a separate,
   reusable routine without changing its policy.
2. Dispatch `refine3D_states` with the checkpoint, remaining stage plan, and
   update-epoch metadata.
3. Compare the handoff path with the current in-controller post-split path.
4. Remove the duplicate post-split loop only after equivalence gates pass.
5. Retain the pre-split ab initio and checkpoint-construction behavior.

### Stage 6: Complete documentation and validation

1. Update tutorials, examples, policy cross-links, generated UI reviews, and
   algorithm descriptions.
2. Remove stale references to the old commands from active policy documents.
3. Consider moving stochastic split construction from `abinitio3D` into
   `refine3D_states` only as a reviewed follow-up.

## Validation Plan

### Static and interface validation

- New commands are registered once and route to the intended commanders.
- Removed command names are neither registered nor routed.
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
- `pose_policy=fixed` preserves the projection direction, optimizes only the
  in-plane angle and translations, and performs stochastic state search.
- `pose_policy=local` uses `prob_neigh_mode=geom`, evaluates every active state
  in the geometric neighborhood containing the current projection, and never
  invokes coarse state-pooled peak selection.
- Automatic local angular, in-plane, and shift bounds reproduce the intended
  sampling-derived neighborhood; each explicit override is honored and
  validated without changing the other automatic limits.
- `pose_policy=global` uses `prob_neigh_mode=state`, performs full
  probabilistic matching over all pose degrees of freedom, and evaluates the
  same pooled candidate geometry for every active state.
- The `refine3D_states` default is `pose_policy=global`.
- `pose_policy=fixed` preserves projection direction while committing the
  optimized in-plane angle, translations, and stochastic state choice.
- Focused state sorting affects classification evidence but not final map-mask
  policy or the stored pose beyond what the selected pose policy permits.
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

### Pose-search breadth masquerading as state separation

Broad search can explain orientation error as conformational variation, while
overly narrow search can lock in a poor scaffold. The three policies must
therefore remain visible scientific choices: `fixed` anchors projection
direction, `local` uses the current geometric neighborhood, and `global` uses
state-pooled probabilistic matching. Within either neighborhood mode, every
state must be compared over the same candidate geometry so state identity does
not alter pose-search opportunity.

### Focus-mask leakage

A focus mask used for state evidence must not silently become a reference,
FSC, automask, NU-filter, or final reconstruction mask. Those policies have
different scientific meanings and owners.

### Frequency-dependent state bias

Allowing different states to match at different effective bandwidths can bias
competitive assignment toward the reference with more permissive evidence.
The first implementation uses a common schedule across states.

### Clean-break adoption

Scripts, restart manifests, UI policy tables, and documentation must adopt the
canonical names and explicit pose policy. Validation must catch stale command
names rather than silently selecting a historical contract.

### Classification name obscuring map production

`classify3D_refs` still reconstructs and postprocesses maps. Its UI summary and
documentation must state this output contract clearly while keeping particle
classification, rather than iterative improvement of independent references,
as the primary scientific intent.

## Resolved Review Decisions

Cyril approved all review items on 2026-08-31:

1. The independent-reference command is `classify3D_refs`, displayed as
   **Reference-guided 3D Classification**.
2. `refine3D_states` defaults to the wide search named
   `pose_policy=global`; there are no legacy aliases.
3. `pose_policy=local` reuses the existing geometric neighborhood policy
   (`prob_neigh_mode=geom`); angular, in-plane, and shift bounds are automatic
   and may be explicitly overridden by the user.
4. A focus mask affects classification evidence only; authoritative
   reconstruction, map masking, and pose policy remain separate.
5. All states initially share one common frequency schedule; state-specific
   adaptation is deferred.
6. The first `abinitio3D` handoff retains current split-checkpoint behavior,
   extracts its construction into a separate reusable routine, and delegates
   only the post-split refinement.

## Completion Criteria

This refactor is complete when:

- users see two purpose-specific applications with distinct input and pose
  contracts;
- state refinement supports existing, stochastic, and FLEX initialization;
- `fixed`, `local`, and `global` state sorting are explicit and validated,
  with `global` the new-command default;
- focus evidence is independent of pose and reconstruction masking policy;
- both workflows consume shared frequency-stage planning rather than copied
  loops;
- `abinitio3D` delegates post-split work to `refine3D_states` through a tested
  checkpoint contract constructed by a separate reusable routine;
- source, policies, algorithm descriptions, UI metadata, and restart behavior
  agree with the new names;
- user-side build and runtime validation, including the `abinitio3D`
  equivalence matrix, have been observed and recorded here.
