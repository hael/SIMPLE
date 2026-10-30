# refine3D_local workflow and architecture cleanup

## Contract and status

Contract: FINAL for the existing local-continuation workflow and the
behavior-preserving corrections requested from Hans's review of `0c4801f`.
The maintainer requested the remaining findings on 2026-10-28. Follow-up source
changes are uncommitted on the current checkout. The maintainer reviews and
accepts the diff; Hans's review recommendations do not imply acceptance of this
implementation.

Execution: the architecture cleanup passed the reviewed Linux fast gate,
focused NU handoff and full integration gate on 2026-10-29. The subsequent
iteration-reference correction passed the focused Linux staging suite,
fast gate and full integration test recorded below. These passes do not
establish recovery on the Streptavidin benchmark.

### Required outcome

Use the existing `refine3D refine=cont` matcher for single-state continuation
from saved poses and genuine half-maps. Remove avoidable startup particle
passes by staging compatible references, committed sigmas and accumulator
chains. Keep the existing iteration strategy and shared final reconstruction.
The public behavior and restrictions are defined once in the
[local-refinement policy](../../policies/3D/refine3D_local_policy.md).

The architecture cleanup must meet these requirements:

1. The commander owns defaults, CLI validation, the single parameter constructor,
   workflow sampling/grid/low-pass policy, child-program calls and final
   reconstruction. State-specific work-project selection stays here because it
   executes the selection child.
2. The startup controller owns the independent source preflight before parameter
   initialization, ordinary project loading/publication and reference/sigma/chain
   staging. Load the selected work project once for the staging operations.
3. The startup strategy owns file destinations, copying, state remapping and
   cleanup. Provenance owns record publication, identity validation and marker
   names. A source path is not overwritten with its staged destination.
4. Provenance has named private identity fields and a reusable validated record.
   Preflight validates reference constituents once; project matching and startup
   staging consume that snapshot. Preserve record version 2 and its serialized
   field order; this is not a new persistent format or a generic cache.
5. Chain compatibility and project discovery each have one owner. Preserve each
   caller's different rejection, logging and cleanup behavior.
6. Portable tests live beside their owners, have isolated fixtures and execute
   without the Linux rejection dispatcher. Keep all meaningful assertions.
7. UI constructors add their command-specific controls explicitly. Their helper
   contains only common inputs. Preserve controls, choices and defaults; correct
   registration help to describe `refine=greedy`.

No new program, CLI parameter, numerical threshold, CMake target, per-source
compiler flag or change to NU filtering is part of this cleanup. Preserve the
existing optimizer, final reconstruction and scientific acceptance criteria.
No main-repository commit or push is authorized.

### Preserved behavior and decisions

- Local selects continuous pose refinement. Auto retains discrete refinement and
  optional `pose_cont=yes` polishing. Base `refine3D refine=cont` stays available.
- Local defaults to `autoscale=no`; auto defaults to `yes`. A deliberate crop and
  explicit `trs` are honored. Native sampling comes from project metadata.
- The sample target remains 25,000. Automatic budgeting targets four updates per
  active particle, with the existing cap and a three-iteration minimum. Explicit
  `maxits` also sets `minits` and can override that minimum. Rotation and shift
  bounds apply around each iteration's seed pose.
- `state=X` creates a single-state work project; frozen `m_estimator=flex` uses
  the positive applied-weight rows. Selection can invalidate sigma/chain reuse.
- Source preflight requires the merged map and genuine even/odd files, compatible
  dimensions/sampling/extent, valid stored poses, binary half labels and selected
  particles in both half-sets. It runs before parameters can synthesize halves.
- A valid producer record binds its constituents to the source project/state,
  ordered rows, halves, population and settings. Independent explicit external
  references retain their caller-supplied alignment/half contract. Same-lineage
  row or half changes cannot use that exception. These checks do not establish
  independent reconstruction or prove a common coordinate frame.
- An absent record permits a logged legacy base pair. Incomplete, mixed or
  modified recorded sets fail preflight. Unverified derived products use the
  existing logged base-half fallback; missing FSC is rebuilt from half-maps.
- Compatible NU reuse carries recorded filtered maps, local resolution, masks and
  matching LP. Unrecorded destination masks are removed. Source files and
  same-file copies are protected. Missing derivatives do not trigger a startup
  particle reconstruction.
- Sigma reuse requires the canonical store and producer identity to match native
  geometry, ordered observations, halves, population and noise settings. Sigma
  estimates may lag poses. Validate both the source and the staged copy. A rejected
  or absent store uses the existing image-power bootstrap owner; no unconditional
  startup power/reconstruction pass is added.
- A chain must match its backend, rows, halves, population, weights, saved poses,
  constituent identities and allowed grid/extent. Do not manufacture sufficient
  statistics from completed half-maps. Copy compatible artifacts with the manifest
  last; otherwise iteration 1 is a counted full update that seeds the chain.
- `refine3D_startup` represents ordinary startup, chain reuse or chain seeding.
  It shares the temporary first-update transition across execution modes and
  restores the requested fraction/flags after seeding. `trail_seed` remains an
  internal transport/reconstruction handshake.
- Preserve inline reconstruction from the matcher's particle buffer, volume-owned
  assembly/filtering and post-assembly sigma commit. Rebuild Cartesian matching
  references each iteration. Local does not become fixed-reference fitting.
- Call `calc_final_rec` once. Its native-grid sigma checks and cropped fallback
  retain their own policy and cold-solve budget. A fixed-pose `refine=sigma`
  fallback is not a pose search, but can prepare a polar model. Do not relabel an
  iteration/partial map as a final all-particle reconstruction.

### Decision sources and rationale

The preceding rules define the current contract. These sources explain their
origin and purpose:

- **H19:** Hans (`hael`) committed the original proposal as `aa7feb246` on
  2026-10-19: *planned refine3D_local workflow for refine=cont pose polishing*.
  Its recommendations were proposal choices, not evidence of implemented code.
- **I27:** Mazharul committed the first implementation deliverable as
  `0c4801f61` on 2026-10-27. That revision retains the implementation records,
  corrections and supplied Linux evidence from 2026-10-23 onward.
- **R27:** Hans's supplied *Architectural and Style Review: Commit 0c4801f*,
  dated 2026-10-27. Findings 1-8 concern ownership, duplication, state, test
  placement and UI structure. The review did not revalidate numerical results.

Read either historical plan without changing the checkout:

```bash
git show aa7feb246:doc/refactoring_notes/planned/refine3d_local_workflow.md
git show 0c4801f61:doc/refactoring_notes/planned/refine3d_local_workflow.md
```

| Decision | Reason | Source |
| --- | --- | --- |
| Extract an explicit local command and reuse `refine=cont` | The existing optimizer already refines saved poses. The old auto route still recalculated power spectra and reconstructed startup maps. Extraction removes avoidable startup work without creating another optimizer. | H19, review findings 1-2 and public contract |
| Exclude global registration and CC pose acquisition | This workflow continues an existing alignment. A new pose-acquisition pass changes that contract. Retiring `pose_cont=only` gives the continuation workflow an explicit entry point. | H19, public contract, original lines 85 and 99 |
| Validate original halves before `params%new` | Parameter setup can copy a merged map into missing half-map files. Validation afterward cannot establish that the original input pair existed. Geometry and recorded identities do not prove independent reconstruction or physical alignment. | H19, review finding 4 and reference handoff |
| Keep separate half references; use `lpstop` for a resolution cap | Explicit `lp` and `nonuniform_lpset` can activate merged-reference matching in the shared loader. Rejecting them preserves the half-specific contract. The later review correction enforced this existing requirement. | H19, public contract, original line 112; I27, review corrections |
| Default local to `autoscale=no` | Native sampling avoids automatically entering the crop-to-native final sigma-bootstrap path. Explicit scaling remains available with its finalization cost. The maintainer reconfirmed this choice after withdrawing the temporary edit to `yes`. | H19, original line 107; maintainer confirmation during this review |
| Keep native geometry project-owned and validate deliberate crops | Fourier cropping can reuse compatible finer maps. Padding a coarse map cannot recover missing frequencies. Native metadata, physical extent and the initial matching band must remain consistent. | H19, reference handoff |
| Validate coherent reference bundles and committed sigma stores | Filenames, timestamps and byte counts cannot establish a common generation. The canonical sigma layout check alone does not verify every half/population/noise setting. Reference and sigma generations can legitimately lag poses; a changed `nparts` alone does not invalidate a store. Never copy pending transactions or worker fragments as committed state. | H19, reference handoff; I27, provenance corrections |
| Recover derived products through volume owners or use the stated base-half fallback | Missing FSC/NU/mask products do not justify a particle pass. Recovery must use the inputs required by the chosen estimator. Do not undo an ML prior or substitute masking for a PCG constrained solve. Missing required `_unfil` inputs need a clear fallback or failure. | H19, reference handoff |
| Reuse complete accumulator chains or count a full first iteration | Finished maps omit the Fourier sums/densities or PCG statistics needed for fractional reconstruction. The first full iteration both refines poses and seeds those statistics. A fresh invocation is not interrupted-run resume; forcing `continue=yes` does not implement this policy. | H19, review finding 7 and iteration policy |
| Preserve the existing iteration pipeline and sampling budget | Inline accumulation avoids another particle read. Assembly owns filtering; sigma commit follows assembly. Keeping the inherited sample/iteration budgets isolates startup changes from optimizer tuning. Stale sampling generations must not count toward new convergence. | H19, iteration policy and public contract |
| Keep the shared final reconstruction | An iteration map can be cropped, fractional or trailed, and can use an older sigma generation. A final map uses all eligible rows, native sampling and its own solver/postprocessing policy. Skipping it requires explicit equivalence across those inputs and outputs. | H19, review finding 8 and finalization policy |
| Narrow startup ownership and share repeated rules | Repeated project reads, implicit host state and duplicated compatibility/discovery rules obscure responsibility and can diverge. The maintainer limited extraction to startup, following `refine3D_states`; defaults, child calls and final reconstruction remain in the commander. | R27, findings 1-6; maintainer's scope correction |
| Move portable checks beside their owners and keep UI differences explicit | A Linux-only rejection harness must not suppress portable provenance, staging and parser checks. Common UI inputs stay shared; constructors describe only their command-specific controls. | R27, findings 7-8 |

### Acceptance contract

Existing assertions protect missing/invalid input rejection before directory
creation, source preservation, same-file safety, constituent mutation, ordered
rows, half/population/observation identity, lagged sigma/reference poses, exact
chain poses, mask remapping/removal, matching LP and fractional transitions.
The analytic impulse fixture has independently known FSC +1/-1 in every shell.
No tolerance or acceptance floor changes are authorized.

H19 also defines broader acceptance: no startup particle pass for a complete
compatible bundle/store; volume-only recovery of missing derivatives; exactly
one necessary sigma bootstrap; simulation-truth pose checks; fractional-update
parity; and native/cropped, weighted-state and backend finalization. The current
test table records coverage, not a claim that this whole matrix has passed.

Any performance comparison must use the same project and settings as the old
continuous-only route. Evaluate the native-grid default separately. Account for
startup particle passes, end-to-end time, peak memory and preparation/matching/
assembly/sigma costs. Removing calls alone does not establish a speedup.

The shared manifest predicate must accept equal identity and same-extent grid
padding, and reject changed population, state layout, weights, component sizes,
a larger prior grid or different extent. Its existing 1% extent tolerance is
preserved. Assembly still discards invalid chains; local startup still refuses
reuse. PCG's additional raw-header validation stays with its existing owner.

Project discovery retains the previous-numbered-directory preference, including
an empty latest directory, and leaves unique-project enforcement/diagnostics to
the caller. Discovery does not create run directories or mutate the command line.

Source-level UI comparison must preserve all 32 auto and 31 local input
specifications, apart from corrected registration prose. Generated display
defaults must remain identical. Linux compilation, affected unit suites, the
fast gate and the focused positive NU-handoff test are required. Static checks
alone establish neither runtime preservation nor reconstruction quality.

## Implementation plan and ownership

The scoped extraction follows the existing `refine3D_states` separation between
commander sequencing and a focused controller. It reuses `parameters`, `cmdline`,
`sp_project`, the two existing refinement strategies and `calc_final_rec`.

| Owner | Responsibility | R27 finding |
| --- | --- | --- |
| `simple_commanders_refine3D` | Defaults/validation, parameter construction, state-selection child, sampling/grid/budget/low-pass decisions, refinement children and final reconstruction | 1, scoped by the maintainer |
| `simple_refine3D_workflow_controller` | Preflight, source path capture, validated reference snapshot, ordinary project loading and startup publication | 1-2 |
| `simple_refine3D_strategy` | Reference/sigma/chain staging and the shared state-row selector; no startup project-file rereads | 3, 6 |
| `simple_refine3D_startup` | Ordinary/reuse/seed state and effective first-update sampling | 5 |
| `simple_refine3D_provenance` | Named identity fields, validated record lifecycle, publication, marker names and compatibility validation | 3-4 |
| `simple_trail_chain_manifest` | Shared population/state/weights/component/grid predicate; no caller cleanup policy | 6 |
| `simple_syslib:resolve_project_files` | Common read-only project discovery for parameters and preflight | 6 |
| Owning tester modules and existing unit areas | Portable isolated fixtures; Linux guard limited to rejection children; see the test table below | 7 |
| `simple_ui_refine3D` | Common workflow inputs plus explicit auto/local constructor controls | 8 |

The validated reference object is created before `params%new`, while source paths
still refer to the input directory. It retains absolute constituent paths through
directory creation. Project matching and staging reuse that object. Sigma and
chain copies still validate both source and destination stores where required.
No whole-`parameters` copy or new build target is needed.

The implementation steps are represented by review findings 1-8: scoped startup
ownership and one loaded work project; provenance/staging separation; named
validated provenance; explicit shared startup state; common rules; test ownership
and portability; explicit UI surfaces. Naming and formatting follow those owners.
The chronology below preserves delivery and review milestones. The decision
table preserves their rationale. Detailed failure logs and intermediate handoffs
remain in I27; they are not current acceptance evidence.

### Test ownership and Linux validation

| Area / selector | Guarantee |
| --- | --- |
| `unit_core suite=syslib` | Project discovery, uniqueness candidates and previous-run precedence |
| `unit_core suite=parameters` | Real parameter construction preserves full/focused high-level suite selection |
| `unit_reconstruction suite=refine3d_provenance` | Identity mutation, constituent integrity, reference rules and record lifecycle |
| `unit_reconstruction suite=trailing_chain_identity` | Manifest round trip, common compatibility predicate and PCG raw-header identity |
| `unit_cart_align3D suite=local_startup` | Linux-only child-process rejection, diagnostic and source-preservation harness |
| `unit_cart_align3D suite=local_staging` | Portable reference/NU/FSC/chain staging, stale-mask cleanup and in-place safety |
| `unit_cart_align3D suite=startup_state` | Requested 0.25 becomes 1.0 only for a seed iteration, then returns with original flags/budget; ordinary/reused chains retain 0.25 |
| `unit_ui` | Registered public surfaces and generated defaults |
| `cont_refine3D_1jxy suite=nu_handoff` | Assembly-produced NU bundle, first matching band, continuous work, fixed halves and source identities |
| full `cont_refine3D_1jxy` | Base CC/Euclid/distributed continuous routes, local continuation and discrete-plus-continuous workflows |

New tester modules use the existing source glob and area-suite registration.
Reconfigure the existing Linux CMake build with its cached options before the
incremental build so those files enter the build. CMake source files and the
15-process fast-gate budget are unchanged. Run the affected unit areas/fast gate
before the focused NU-handoff. Reverse-order unit execution is useful for checking
the moved fixtures. Windows/MSYS builds and SIMPLE execution are prohibited.

Broader PCG/distributed workflow, crop/finalization and matched real-data
measurements remain separate validation. A successful handoff does not prove
map quality, global pose recovery or an end-to-end speed improvement.

## Implementation and evidence record

### Implementation and review milestones

| Date / phase | Deliverable or decision | Source and evidence limit |
| --- | --- | --- |
| 2026-10-19: original proposal | Hans proposed a separate continuous continuation workflow, validated startup reuse, compatible fractional chains and shared final reconstruction. | H19; source review and proposed acceptance criteria, with no implementation or runtime claim |
| 2026-10-23: workflow boundary | Added the command boundary and preflight before parameter-side half synthesis. Migrated `pose_cont=only` callers while keeping auto polishing and base continuous refinement. | I27, preparatory/preflight/boundary records |
| 2026-10-26: startup and continuity | Added saved-reference/sigma reuse, chain staging or full-first-update seeding, provenance and deliberate-crop handling. Same-file copy failures led to source-preservation checks. Temporary focused test stages were restored after focused Linux debugging. | I27, reuse/continuity, provenance, review-correction and gate-restoration records; tests apply to their recorded snapshots |
| 2026-10-27: first committed delivery | Published `refine3D_local` and its tests, policy and implementation record, including the correction that removes destination masks absent from the verified source set. | I27; historical Linux results are summarized below; broader scientific acceptance was not established |
| 2026-10-27: architectural review | Hans reviewed the committed delivery and raised eight maintainability findings plus naming, formatting and documentation corrections. | R27; source review, with no independent numerical revalidation |
| 2026-10-27 to 2026-10-28: response to review | The maintainer narrowed controller extraction. Implemented explicit startup state, one loaded work project, strategy-owned staging, a provenance class, shared rules, portable test ownership and explicit UI controls. Integrated upstream changes and preserved both sets of test registrations. | Current uncommitted cleanup; finding-to-owner mapping is above; recommendations do not imply Hans accepted this implementation |
| 2026-10-29: compiled and tested current delivery | Corrected the polymorphic filename constructor and project argument intent. The maintainer reported a successful Linux build; the supplied fast-gate and focused NU-handoff logs passed review. | Latest evidence below; these results do not validate every backend, finalization route or scientific claim |
| 2026-10-29: full integration reviewed | The supplied full gate passed before and after the iteration-reference correction; results are separated below. Manual review continues. | Current delivery remains uncommitted; final acceptance is not yet recorded |

The first committed delivery and the current review-response delivery are
different source states. Their evidence must remain separate. The current
deliverable is the startup/ownership cleanup, its tests and these documents.

Current source implements findings 1-8 and the related naming/style corrections.
No new numerical algorithm, public setting or persistent format was introduced.
The portable tests are registered separately and the Linux guard covers only the
rejection dispatcher. Existing assertions were moved with their subjects;
compatibility/discovery and provenance-lifecycle checks were added for the shared
APIs. Current Linux evidence is recorded below, separately from historical runs.

Static checks observed for the cleanup: descriptor/source-unit check, test
registry, FLEX dependency DAG and diff whitespace checks pass. A source comparison
retained all 70 original startup assertion statements and found no module-use
cycles across 731 modules. The UI comparison preserved 32 auto and 31 local input
specifications. Display-default generation was
identical before/after (SHA-256
`9b2258e2be3256b19e12cff180c3dae37c31410a6352fc603696c0d3f76dcaaa`).
These checks do not replace the Linux gates above.

### Linux validation reported and reviewed on 2026-10-29

The maintainer reported a successful Linux rebuild after two bounded source
corrections: construct the sigma filename array from a concrete string local,
and propagate the chain stager's `project` argument as `intent(inout)` through
the controller. No numerical settings or acceptance thresholds changed.

| Check | Observed evidence | Scope |
| --- | --- | --- |
| `scripts/run_fast_gate.sh build` | Supplied `run_fast_gate.log`: 15/15 entries passed, 0 failed, 41.54 s wall time | Includes the affected core, reconstruction, Cartesian and UI unit suites |
| `cont_refine3D_1jxy suite=nu_handoff` | Supplied `nu_handoff.log`: focused integration PASS, normal stop, 419.9 s | Committed sigma reuse, verified NU references without import fallback, first matching LP 3.22759 A, continuous work, unchanged half labels and source constituents |

The NU log reports Git stamp `c97106c4`; that stamp does not identify the
uncommitted patch. Attribution to the corrected checkout uses the maintainer's
rebuild report and supplied logs. The supplied full `cont_refine3D_1jxy` log
also passed (`suite=all`, 3652.5 s), including shared/distributed agreement,
fractional and in-place continuation, polishing and NU handoff. These passes
predate the iteration-reference correction below. They are regression evidence,
not matched performance or map-quality evidence; broader PCG coverage remains
unestablished.

### Historical Linux evidence: earlier source only

| Earlier milestone | Observed evidence | Limit |
| --- | --- | --- |
| Workflow-boundary phase, reviewed 2026-10-26 | Recorded fast gate 15/15 in 36.55 s, UI 894/894, Cartesian 808/808, continuous integration PASS in 2104.8 s and FLEX integration PASS in 2081.9 s | Boundary/state-work-project evidence only; startup reuse and the later architecture cleanup were not yet validated |
| Pre-stale-mask focused startup | 364/364 assertions in 4.30 s; fast gate 15/15 in 49.41 s | Superseded by the later stale-mask correction |
| Focused positive NU handoff | PASS in 949.0 s | Predates the stale-mask correction and current cleanup |
| Post-stale-mask correction, 2026-10-27 | User-reported compilation; supplied logs: startup 369/369 in 7.87 s, UI 927/927, fast gate 15/15 in 74.48 s | Predates the current controller/provenance/test/UI cleanup; fast-gate budget warning was not a test failure |

Historical 35-minute/full and 27-minute/focused runs must not be attributed to
later source changes or used as matched performance evidence.

### Iteration reference transition correction

Approved contract: after successful assembly, adopt its merged/even/odd triplet
in both execution modes; retain inputs for states without replacements. The
master's half-map paths previously remained at startup names even after the
command-line merged reference advanced. One strategy-owned adoption routine now
updates parameters, command-line state and distributed worker descriptions
together. Existing PCG reference refresh behavior is preserved.

The local-staging regression uses known, distinct startup/raw/NU impulses and
checks both halves through the production reference reader across two updates.
This protects against stale startup paths and raw fallback when NU maps exist.
Optimizer settings, budgets, assembly products and project ownership are unchanged.

| Correction check | Reviewed Linux result |
| --- | --- |
| `local_staging` | 91/91 assertions, 0.06 s |
| Fast gate | 15/15 entries, 41.51 s |
| Focused NU handoff | PASS, 508.6 s |
| Full `cont_refine3D_1jxy` | `suite=all`, PASS, 4136.7 s; includes in-place continuation and NU handoff |

Raw-half fallback occurs at startup; later local iterations consume NU products.
These results supersede pending statements for this correction. Source association
uses the maintainer's rebuild report because the logs retain stamp `c97106c4`.

### NU and cFAR findings: closed as workflow blockers

The maintainer deferred further NU/metric research: neither observation establishes
a refine3D_local defect or reconstruction quality. Final cFAR=0 was reproduced
from retained maps and follows the existing first-supported-shell scoring rule:
two directions have one sample at k=2 with FSC -0.55494 and therefore zero area.
Producer cFAR 0.8331227 matches its native 0.8331. Empty NU evidence already occurs
before local refinement. The replay found zero seeds at threshold 9.6440, before
ICM or component removal. With automsk=yes the density envelope remains active;
automsk=nu would require density fallback. Threshold suitability remains unresolved.
No NU mathematics, metric scoring or production threshold was changed.

The diagnostic replay FAILED its declared 1e-5 saved-product comparison limit;
geometry, finite values, matching LP and source-preservation checks passed. The
small output mismatch remains unexplained and is not a passing equivalence result.
The maintainer approved removing only its nu_envelope route, settings, routine and
UI entries; removal and static registry/descriptor/whitespace checks are complete.
Workflow regression tests, findings, logs and downloaded maps were retained.

Detailed counters, errors, hashes and retained replay log are in the local
[evidence record](../../../codex_logs/cfar_nu_handoff_20261029_134846/RESULT.md).
Earlier focused evidence (817.4 s) reported producer/local/final cFAR
0.8331/0.9173/0.0000 and FSC=0.143 approximately 3.228/3.120/3.120 A; it is
historical diagnostic evidence, not a matched quality comparison.


Upstream integration: rebased onto origin/master at 47b34ad4d. Retained upstream
half-map sharpening and top-k changes; resolved the test-help conflict by keeping
both half_map_sharpening and refine3d_provenance in registration order. Registry,
descriptor, FLEX DAG and staged whitespace checks passed. Earlier Linux results
precede this integration; the combined source has not been rebuilt or run here.
The pre-rebase recovery snapshot and stash remain available until acceptance.
