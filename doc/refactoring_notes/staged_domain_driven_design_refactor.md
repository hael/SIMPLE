# Staged Domain-Driven Design Refactor

**Status:** Revised proposal; review and X design-case-study actions incorporated

**Implementation status:** No production refactoring has started

**Date:** 2026-08-17

**Review incorporated:**
`doc/refactoring_notes/staged_domain_driven_design_refactor_review.md`

**Design case study incorporated:** X domain-driven architecture and source
review, 2026-08-17

**Delivery model:** Small, validated commits directly on `master`; no isolated
long-lived implementation branch

## 1. Executive decision

SIMPLE should move toward an explicit domain model, but it must not be rewritten
as a new object hierarchy or converted to a new project format in one operation.
The safe approach is a staged, compatibility-preserving migration around the
existing `sp_project`, `ori`/`oris`, commander, strategy, and numerical types.

Every stage must be independently releasable on `master`. A stage may add
types, adapters, validators, or an opt-in execution path, but it must not leave
the default workflow partially migrated or dependent on a later commit. The
accepted legacy behavior and project representation remain authoritative until
the replacement path has passed explicit parity gates.

This is a correctness refactor first. Architectural cleanup, removal of old
APIs, persistence redesign, and new user-visible behavior are later decisions.

## 2. Non-negotiable delivery rules

1. Work is integrated through focused commits on `master`. Do not accumulate
   the refactor on an isolated long-lived branch.
2. `master` must remain usable after every commit. No commit may rely on a
   follow-up commit to restore compilation, project compatibility, or workflow
   behavior.
3. Preserve the current `.simple` representation during the initial migration.
   A persistence-format change requires its own later proposal, versioning,
   migration tooling, backward-read tests, and explicit approval.
4. Preserve existing command-line defaults and numerical behavior. A genuinely
   new behavior or alternate execution path must use an explicit opt-in key
   whose default is disabled.
5. Preserve the legacy path's observable contract, not necessarily its source
   implementation. Before replacement, find or add an independent
   characterization test or accepted baseline that passes on the legacy code.
   The typed implementation may then replace the legacy path in the same
   focused commit, provided the unchanged characterization gate passes. A dual
   production path, opt-in switch, and separate removal commit are not required
   for a behavior-preserving refactor.
6. Each commit must include or extend the tests needed to prove its local
   contract. Test-only preparation commits should precede risky production
   changes when an adequate baseline does not already exist.
7. Use exact comparisons for identifiers, states, indices, strings, dimensions,
   flags, and persisted field presence. Use explicitly justified tolerances only
   for floating-point numerical results.
8. Preserve shared-memory, distributed, streaming, restart, and project-merge
   behavior whenever a touched contract is consumed by those paths.
9. Do not claim a platform, compiler, BOX, distributed, or real-data result
   without observed output from that environment.
10. Agents perform lightweight source and syntax checks unless compilation is
    explicitly requested. The user or CI performs the required compilation and
    runtime validation before the corresponding production commit is accepted.
11. Every opt-in key must follow the normal
    `cmdline -> UI metadata -> parameters -> commander/strategy` path, including
    default-off behavior and validation. Generated UI or argument sources must
    be regenerated through their owning generator and must not be treated as
    independent handwritten authorities.
12. Domain access contracts are owned by the domain or application layer that
    consumes them. `sp_project`, `ori`/`oris`, file formats, command parsing,
    memory mapping, and external-tool representations implement or adapt to
    those contracts; they do not define the domain API.
13. When an optimized representation or algorithm replaces a simpler semantic
    implementation, protect the observable contract with an independent
    oracle. Prefer a deliberately small executable reference where its ongoing
    value justifies the maintenance cost; otherwise use characterization
    fixtures, analytical expectations, or accepted snapshots. This rule does
    not require retaining the replaced production implementation.
14. Construct scientific policy and immutable use-case configuration once at
    the commander/application boundary. Domain services must not inspect method
    names or unrelated fields from the global `parameters` object inside hot
    operations.
15. Performance work must remain self-verifying. A benchmark that exercises an
    optimized, parallel, distributed, persisted, or mapped path must also check
    its accepted result against the applicable independent oracle; timing alone
    is not an acceptance gate.

## 3. Why the current model is difficult to migrate

SIMPLE already contains important domain behavior:

- commanders represent application commands and high-level workflows;
- strategies represent algorithm and execution policy;
- `image`, `ctf`, motion, picking, orientation, reconstruction, and related
  types contain scientific behavior;
- `ctfparams` and several other small structures already resemble value
  objects.

The main problem is the information model. `sp_project` is simultaneously a
project aggregate, persistence gateway, artifact catalog, current analysis
state, and cross-record relationship store. Its movie/micrograph, stack,
particle, class, output, and optics segments are all represented by generic
`oris` collections. Individual `ori` records contain a fixed particle parameter
array plus numeric and character hash tables. Domain distinctions therefore
depend on segment selection, row position, string keys, and caller knowledge.

Examples include:

- movies and integrated micrographs sharing `os_mic` rows and being
  distinguished by `movie`, `intg`, and `imgkind` keys;
- stack and physical-image identity being encoded by `stkind`, `indstk`,
  `fromp`, and `top`;
- CTF parameters being resolved from a combination of stack and particle rows;
- `ptcl2D` and `ptcl3D` representing two mutable current-state projections of
  the same particle population;
- analysis procedures overwriting, transferring, or deleting orientation,
  class, state, and score fields in place.

The first migration goal is therefore typed access to existing information, not
replacement of the numerical algorithms.

## 4. Proposed domain language

The low-level names that become public in M2 are not provisional: M0 must
approve and freeze the identifier, quantity, geometry, and image-location type
names before those modules are committed. The broader aggregate and analysis
names below remain provisional until the owning workflow is migrated. This
allows the ubiquitous language to improve without renaming foundational public
APIs after M2.

| Concept | Proposed role |
|---|---|
| `Microscope` | Entity representing an identified instrument when instrument identity matters |
| `OpticsConfiguration` | Value object for the acquisition/CTF model settings used by a dataset |
| `Movie` | Acquired-data entity with stable identity and image location |
| `Micrograph` | Derived-data entity with explicit lineage to a movie or imported integration |
| `ParticleObservation` | Extracted image observation; do not imply identity of the physical molecule |
| `ParticleSet` | Aggregate and batch-processing boundary for particle observations |
| `DensityMap` | Volumetric density artifact; replaces one ambiguous use of `Model` |
| `ReferenceSet` | Set of 2D or 3D references used by search or refinement |
| `AtomicModel` | Atomic-coordinate model, distinct from a density map or learned model |
| `AnalysisRun` | Identified execution record with inputs, configuration, status, and outputs |
| `MotionCorrectionResult` | Motion field and derived micrograph produced from a movie |
| `CTFEstimate` | Per-micrograph or per-particle defocus, astigmatism, and phase estimate |
| `PickSet` | Identified set of positions produced by a picking run |
| `Classification2DResult` | 2D alignments, assignments, scores, and class-average artifacts |
| `Refinement3DResult` | 3D poses, state assignments, scores, and reconstruction artifacts |
| `HeterogeneityResult` | Latent coordinates, state interpretation, and related artifacts |

`search`, `mask`, and `filter` are normally domain services or operations over
explicit specifications. `compute` is too broad to be useful domain language.
The numerical implementation should remain in its owning subsystem.

## 5. Initial bounded contexts

The boundaries below are intended to prevent a replacement `sp_project` god
object. They do not require separate libraries or processes.

### 5.1 Acquisition

Owns microscope identity, acquisition/optics configuration, movies, imported
micrographs, dose/frame metadata, and source image locations.

### 5.2 Preprocessing

Owns motion-correction runs and results, CTF estimation, picking, extraction,
and lineage from movies through micrographs and pick sets to particle sets.

### 5.3 Particle analysis 2D

Owns 2D alignment, classification, class assignment, class averages, class
selection, and the corresponding current-state materialization.

### 5.4 Structure analysis 3D

Owns orientation search, 3D refinement, reconstruction, half-set semantics,
state assignment, masks, maps, FSC artifacts, and restart contracts.

### 5.5 Heterogeneity analysis

Owns latent representations, continuous/discrete heterogeneity results, state
interpretation, and links to reconstruction outputs.

### 5.6 Workspace and provenance

Owns project persistence, artifact locations, analysis-run metadata, execution
environment, job status, import/export, and compatibility adapters. A project
is a workspace/catalog and persistence boundary; it should not become the owner
of every scientific invariant.

Opt-in controls are infrastructure exposed through the owning workflow. Their
metadata, defaults, parsing, and validation remain in the established UI,
`cmdline`, and `parameters` pipeline; domain objects must not parse command-line
keys or create a parallel configuration path.

### 5.7 Dependency direction and context exclusions

Bounded contexts need both an ownership statement and a negative dependency
rule. The initial exclusions are:

| Context | Must not know about |
|---|---|
| Acquisition | Motion/CTF algorithms, particle alignment, reconstruction, project segment names, or file serialization details |
| Preprocessing | 2D/3D ranking formulas, heterogeneous-state interpretation, scheduler syntax, or report layout |
| Particle analysis 2D | Micrograph file parsing, CTF persistence keys, 3D reconstruction state, or project binary layout |
| Structure analysis 3D | Acquisition file parsing, 2D class-selection implementation, UI metadata, or scheduler syntax |
| Heterogeneity analysis | Project segment layout, acquisition parsing, execution-host details, or unrelated reconstruction policies |
| Workspace and provenance | Scientific scoring, search, filtering, masking, CTF-fitting, or reconstruction formulas |

Dependencies point inward through narrow contracts. Domain modules may depend
on lower-level numerical utilities and their owning scientific modules, but
they must not import project persistence, command parsing, UI, commander, or
external-format modules. Application orchestration may compose domain and
infrastructure implementations without moving scientific formulas into the
commander.

## 6. Value-object priorities

Start with small types that enforce invariants and clarify units:

- typed identifiers: `movie_id`, `micrograph_id`, `stack_id`, `particle_id`,
  `particle_set_id`, `analysis_run_id`, and `optics_group_id`;
- physical quantities: sampling distance, resolution, voltage, spherical
  aberration, dose, defocus, phase shift, and amplitude-contrast fraction;
- geometry: image dimensions, box size, pixel coordinate, physical coordinate,
  Euler orientation, and in-plane/3D shifts;
- storage references: image path plus physical image index and dimensions;
- scientific specifications: mask definition, filter specification, search
  limits, and symmetry.

Constructors must validate ranges and canonical units. Components should be
private where practical. Do not add type wrappers that merely rename a scalar
without enforcing a unit, invariant, identity domain, or meaningful operation.

The current `ctfparams` bundle should eventually be separated into:

- `ctf_model`: sampling distance, voltage, spherical aberration, and amplitude
  contrast;
- `ctf_estimate`: defocus X/Y, astigmatism angle, and numerical phase;
- `ctf_application`: whether the data are raw, CTF-free, or already
  phase-flipped;
- `ctf_fitting_policy`: runtime fitting choices, including phase fitting.

`ctfparams` remains the compatibility DTO until all relevant callers have been
migrated. `ctf_fitting_policy` is runtime-sourced from validated command-line
and parameter state. It is not project metadata and a project adapter must not
invent, persist, or recover it from mic, stack, particle, or optics rows.
Project round-trip parity covers the model, estimate, and application fields,
while fitting-policy parity is tested at the command/configuration boundary.

## 7. HPC and memory constraint

DDD does not require one allocated polymorphic object per particle. SIMPLE may
process millions of particle observations, so the domain layer must preserve
batch access, compact storage, OpenMP-friendly iteration, and structure-of-
arrays layouts in particle-scale hot loops.

`ParticleSet` should be the aggregate and computational boundary. A lightweight
`ParticleView` may expose one logical record without allocating or copying its
image and metadata. Repositories and adapters may return typed columns or batch
views. Typed identifier and quantity objects are API-boundary types. They must
be lowered to compact intrinsic arrays before hot loops and must not become
arrays of derived types in particle-scale kernels unless a measured, reviewed
exception proves no material memory, vectorization, cache-locality, or
throughput regression. The migration must establish the reference memory and
throughput baselines in M0/M1 and compare against them before a particle
representation is adopted broadly.

### 7.1 Semantic and packed representations

The domain model and the high-performance representation are two views of the
same information. A semantic view expresses identifiers, quantities,
relationships, policies, and invariants. A packed view may use structure-of-
arrays columns, CSR offsets, compact tags, stable integer identifiers, and
caller-owned buffers. The packed representation is not a competing domain
model and must not leak its offsets, storage tags, or ownership rules into
callers.

Contracts used by hot paths must support bulk operations, typed columns, or
batch views where necessary. They must not force one polymorphic call, object
construction, allocation, or formatted lookup per particle. If a semantic
contract is too narrow for an optimized implementation, extend the
domain-owned contract with an explicit bulk capability rather than making the
domain service depend on a concrete persistence adapter.

### 7.2 Characterization and reference oracles

Before replacing a legacy operation, first run its characterization test or
baseline against the unchanged implementation. The replacement may remove the
legacy source immediately when the same test exercises the new implementation
and passes in the focused replacement commit.

For each migrated operation, provide an independent oracle appropriate to its
cost and risk:

1. an existing characterization test, deterministic accepted snapshot,
   analytical expectation, or simple implementation that favors clarity and
   complete semantics;
2. deterministic fixtures that exercise normal, boundary, duplicate, empty,
   sparse, and invalid inputs;
3. exact comparison for discrete results and governed numerical comparison for
   floating-point outputs; and
4. persistence and parallel variants where those are part of the production
   path.

Likely reference-oracle candidates include CTF authority resolution, project
identity/index mapping, sparse-stack lookup, field-transfer operations,
particle selection/remapping, orientation transformations, and bounded search
or ranking. A small reference implementation may remain as test-only code when
it provides unique coverage, but this rule does not require retaining the old
production path, a second full reconstruction engine, or duplicate
implementations whose maintenance cost exceeds their diagnostic value.

OpenMP and packed numerical algorithms may remain in their owning scientific
modules. The dependency rule is that semantic entities and value objects do
not manage threads, file layouts, or execution hosts; it is not a requirement
to misclassify parallel numerical implementation as generic infrastructure.

## 8. Compatibility invariants

The following invariants hold until a separately approved change says
otherwise:

### 8.1 Project and persistence

- Existing `.simple` files remain readable.
- A read followed by a write does not lose or reinterpret known fields.
- Untouched segments and unknown extension keys survive adapter-mediated
  updates.
- Segment row order, row counts, field presence, and project-local foreign keys
  remain unchanged unless the invoked legacy operation already changes them.
- Reading a project does not mutate it on disk.
- Partial segment writes preserve all other segments.
- Any later additive artifact format has explicit magic, schema/version,
  endianness and numeric-kind metadata where relevant, provenance, and
  fail-closed validation before payloads are trusted. Such an artifact does not
  replace `.simple` persistence without a separately approved migration.

### 8.2 Acquisition and preprocessing

- Movie and micrograph paths, dimensions, frame counts, sampling distance,
  state, optics assignment, and lineage are preserved.
- Motion correction and CTF estimation retain their current default policies
  and numerical behavior.
- Picking positions, box interpretation, outside-box policy, and extracted
  physical image order are preserved.

### 8.3 Particle indexing

- `stkind` remains the project-local stack foreign key.
- `indstk` remains the physical image index in the backing stack.
- `fromp`/`top` remain project-row ranges, not substitutes for physical indices
  after a metadata-only prune.
- `ptcl2D` and `ptcl3D` row populations and source mappings remain consistent,
  but their analysis fields are not assumed to be mirrors.
- Merge, prune, split, import, extraction, and restart retain the indexing
  policy in `doc/policies/project_stack_indexing_policy.md`.

Before M7, the project contract must contain a reviewed field-authority
manifest. At minimum it must classify fields as follows:

| Field class | Authority and synchronization rule |
|---|---|
| Source identity and physical indexing | `os_stk` plus the particle segment being consumed; `stkind` and `indstk` must resolve the same source in both particle segments when both exist |
| Shared acquisition/CTF provenance | Stack authority for model fields and particle-segment authority for per-particle estimates; duplicated `ptcl2D`/`ptcl3D` values must agree after extraction or an explicit transfer |
| 2D analysis | `ptcl2D` owns class assignment, 2D alignment, and 2D search/score state |
| 3D analysis | `ptcl3D` owns 3D pose, projection assignment, 3D search/score state, and any structural-state interpretation explicitly assigned to the 3D workflow |
| Selection and even/odd state | Authority is workflow-specific and must be declared before migration; equality between particle segments must not be inferred without that declaration |
| Transfer-derived fields | May cross segments only through named legacy transfer operations with an explicit copied, transformed, cleared, and preserved-field contract |

The completed manifest must name individual persisted keys, including dynamic
keys, rather than relying only on these categories. Adapter validation must use
the manifest for the selected workflow and must not enforce global equality of
all `ptcl2D` and `ptcl3D` fields.

### 8.4 Analysis

- Active/rejected state, even/odd assignment, class/state labels, orientations,
  shifts, correlations/objectives, sampling counters, and update counters retain
  their existing meaning.
- Shared and distributed paths produce equivalent project state within the
  currently accepted numerical tolerances.
- Restarted execution consumes the same persisted state as uninterrupted
  execution.
- Artifact names and project `os_out` references remain stable unless an
  explicit user-visible change is approved.

Semantic comparison is mode-specific. A persistence or adapter round trip must
preserve every persisted `sampled` and `updatecnt` value exactly. A restart
comparison must instead follow the owning workflow's checkpoint contract:

- counters declared restart-invariant are compared exactly at the same logical
  boundary;
- counters that legitimately reset or advance are compared through a named
  normalization or invariant, such as selected-cohort membership, nonnegative
  update tiers, per-state update fractions, or equivalent completed work;
- the raw counter arrays are retained in diagnostic snapshots even when the
  pass/fail comparison uses normalized values;
- no counter may be silently excluded merely because a restart changes it.

Every workflow baseline manifest must classify each restart-sensitive field as
`exact`, `tolerance`, `normalized`, `invariant`, or `diagnostic_only`, and give
the policy reference and rationale. Unclassified fields fail the comparator.

## 9. Project lifecycle states to validate

“State validation” has three meanings in this refactor, and all three are
required:

1. every intermediate refactoring commit leaves the program in a valid state;
2. every supported persisted project lifecycle state satisfies its schema and
   relationship invariants;
3. every scientific workflow state retains its accepted meaning and numerical
   behavior.

The project validator and fixtures should cover at least the following states:

| State | Required validation |
|---|---|
| New/empty project | Metadata segments are valid; data segments may be empty; write/read is stable |
| Imported movies | Movie path, image kind, frame count, dimensions, sampling, CTF model, state, and optional optics identity are valid |
| Imported micrographs | Integration path and dimensions are valid; movie lineage may be absent only for directly imported integrations |
| Motion-corrected | Source movie and integration coexist correctly; motion artifacts and dose/sampling metadata are consistent |
| CTF-estimated | Required CTF-model fields and estimate fields are complete and canonical for active records |
| Picked | Box/pick source exists logically; coordinate and particle-count contracts hold |
| Extracted | Stack ranges, `nptcls`, `nptcls_stk`, `stkind`, `indstk`, `ptcl2D`, and `ptcl3D` mappings are valid |
| 2D initialized | Particle selection/even-odd state is valid before class assignment |
| 2D classified | Class assignments, alignments, class rows, class-average indices, populations, and selected/rejected states agree |
| 3D initialized | Reference/map inputs, half sets, particle states, and initial orientation contract are valid |
| 3D refined | Pose, shift, state, projection index, objective, sampling/update counters, and reconstruction artifacts agree |
| Multi-state/heterogeneous | State labels, per-state populations, per-state artifacts, and particle-to-result relations are complete |
| Pruned/selected | Project rows and physical stack indices obey the sparse-stack policy; selection does not silently renumber physical images |
| Merged/split | Colliding project-local identifiers and foreign keys are remapped exactly once; row data and CTF/optics identity are preserved |
| Persisted/reloaded | Semantic snapshot before write equals the snapshot after read |
| Restarted | Restarted workflow state equals an accepted uninterrupted checkpoint under the workflow's explicit exact/normalized/invariant field contract |

The validator must report the first offending segment, row, field, and expected
contract. It must not silently repair invalid projects unless the invoked
command is explicitly a repair/migration command.

## 10. Test architecture

### 10.1 Existing routines to reuse and strengthen

- `src/utils/simple_test_utils.f90` provides the shared assertion counter and
  typed scalar assertions. New tests should use it instead of private ad hoc
  pass/fail flags where practical.
- `src/main/ori/simple_ori_tester.f90`,
  `production/tests/simple_test_ori.f90`, and
  `production/tests/simple_test_oris.f90` cover orientation value behavior and
  serialization.
- `production/tests/simple_test_sp_project.f90` already covers project
  write/read comparison, partial reads, alignment-document merging, phase
  preservation, and non-mutating reads. It is the immediate home for compatible
  project adapter round trips until dedicated tests are warranted.
- `production/tests/simple_test_binoris_io.f90` and
  `production/tests/simple_test_inside_write.f90` exercise binary orientation
  and in-place segment I/O boundaries.
- `production/tests/simple_test_ctf.f90` is the numerical oracle for the CTF
  model and phase convention.
- `production/tests/simple_test_mini_stream.f90` provides an end-to-end
  acquisition/preprocessing/selection/particle-analysis path when its external
  fixtures are available.
- The deterministic snapshot/comparison pattern in
  `simple_test_continuous_inplane_refine3D` should be reused for workflow
  parity. Generate-only mode is not parity; comparison against an independently
  accepted baseline is required.

Existing tests that only print values or use a mutable `test_passed` flag
without asserting the important result should be strengthened before they are
used as a stage gate.

### 10.2 New common testing routines

Develop these as test-support modules before broad production migration:

#### `project_state_fixture`

Creates deterministic, small, synthetic projects representing the lifecycle
states in section 9. Fixtures must include positive, rejected, even/odd,
multi-stack, mixed-CTF, sparse-after-prune, and multi-state records. It must also
include two source projects whose identifiers deliberately collide. For the
current project merge this includes at least stack and optics identifiers. The
fixture must separately include colliding class, state, and artifact labels
whose owning operation declares them preserved rather than remapped, so the
oracle detects both missed/double remaps and inappropriate remapping. The
fixture must not call the new adapter to construct the expected legacy state;
otherwise the implementation would manufacture its own oracle.

#### `validate_project_state`

Checks structural and semantic invariants for a declared lifecycle state. It
should support validation of a complete project and of an intentionally loaded
segment subset. Production use may be proposed later, but the first version is
a test oracle independent of the adapters being tested.

#### `snapshot_project_semantics`

Writes a deterministic semantic snapshot containing:

- populated segment names and row counts;
- ordered field presence and values for every row;
- resolved movie/micrograph lineage;
- stack/particle project and physical indices;
- CTF model and estimate values;
- state, class, pose, shift, score, even/odd, sampled, and update counters,
  including both raw and normalized restart-sensitive representations;
- output artifact kind, state, path, dimensions, and sampling.

Volatile fields such as absolute execution directory, timestamps, hostnames,
and process identifiers must either be normalized explicitly or compared by a
separate provenance contract. They must not be silently omitted without being
listed.

#### `compare_project_semantics`

Compares two semantic snapshots. Integer, logical, string, identity, field-
presence, and row-order contracts are exact. Floating-point tolerances are
field-specific and documented. The comparison prints the first mismatch and a
bounded summary of additional mismatches.

The caller must select a comparison mode:

- `representation`: exact persisted-field parity for adapter, write/read, and
  partial-segment round trips, including exact `sampled` and `updatecnt`;
- `workflow`: accepted same-boundary semantic and numerical parity;
- `restart`: the field classifications and normalizations declared by the
  owning workflow's restart contract;
- `distributed`: the field classifications and tolerances declared by the
  owning workflow's shared/distributed contract.

The comparator fails when a present field has no rule in the selected manifest.
Normalization is explicit and versioned; it is not an ad hoc list of ignored
fields.

#### `compare_numerical_artifacts`

Compares replacement scientific outputs with an accepted legacy
characterization or other independent oracle using artifact-specific metrics.
The legacy and replacement implementations may be run sequentially across the
replacement commit; they need not coexist in production:

- image, stack, and volume dimensions, sampling, image count, and logical
  indexing are exact;
- bitwise checksums are required only for outputs declared deterministic on the
  tested platform and execution mode;
- image and volume data report maximum absolute error, RMS error, relative L2
  error, and real-space correlation; volume comparisons may additionally
  require shell-wise FSC;
- particle stacks report per-image metrics plus aggregate worst-case and
  percentile summaries, without allowing population averages to hide a bad
  image;
- FSC/FRC and other curves require exact axes and lengths, with maximum and RMS
  ordinate errors;
- artifact sets, kinds, state associations, and required filenames are compared
  semantically before numeric contents.

The initial tolerance registry will be the file-backed table
`production/tests/domain_refactor_numerical_tolerances.tsv`. Each row names the
workflow, artifact kind, comparison mode/platform scope, metric, absolute and
relative bounds, minimum correlation where applicable, baseline provenance,
and rationale. Missing rows fail; test code must not supply permissive fallback
tolerances. The approval owner is the maintainer of the workflow/domain that
produces the artifact, with final acceptance by the SIMPLE project maintainer.
A tolerance change is a reviewed contract change in the same commit as its
rationale and independent baseline evidence; it must not be derived solely
from the new path's observed error.

#### `capture_performance_baseline`

Captures wall time, CPU time where available, peak resident memory, particle or
image throughput, thread count, input dimensions/population, compiler/build
identity, host/platform identity, and execution mode for a fixed workload. M0
selects the representative workloads and platform owners; M1 provides the
repeatable harness and records accepted unchanged-`master` baselines. Timings
use documented warm-up and repeated-run statistics. A later stage may not claim
"no material regression" without comparing the same workload and environment,
or documenting why a controlled replacement baseline is required.

#### `assert_adapter_characterization`

For a given legacy row or dataset:

1. construct the fixture independently of the adapter under test;
2. run the legacy implementation and capture explicit semantic expectations;
3. review and observe the characterization test passing before replacement;
4. replace the legacy implementation and run the unchanged expectations
   against the typed adapter;
5. write through the typed adapter into a copy and compare it with the accepted
   semantic and persisted expectations; and
6. serialize, reload, and compare again.

The legacy and typed implementations may be compared directly while both happen
to exist in a developer worktree, but this is optional. The committed test uses
independent expected results and remains valid after the legacy source is
removed.

#### `assert_oracle_equivalence`

Runs the same deterministic fixture through a simple semantic implementation
and the production implementation, then compares the complete observable
result. The comparison includes ordering and tie-break rules, duplicates,
empty results, explicit failure outcomes, and all provenance needed to
interpret a result. For batch or indexed implementations it also verifies that
bulk access changes how data are obtained, not which domain results exist.

Where the production path supports serial, threaded, distributed, persisted,
or memory-mapped execution, the helper treats those as additional
implementations of the same contract. Performance measurements are recorded
only after equivalence succeeds.

Use this helper only when an executable reference earns its maintenance cost.
Otherwise the same gate compares the production result with the accepted
characterization snapshot or analytical expectation.

### 10.3 Required test classes per stage

Each production migration stage must select all applicable classes:

1. value-object constructor, canonicalization, equality, and invalid-input
   tests;
2. legacy-to-domain and domain-to-legacy mapping tests;
3. complete project and partial-segment round trips;
4. lifecycle-state validation, including negative/corrupt fixtures;
5. replacement-path semantic parity against a characterization test or accepted
   baseline captured from the legacy implementation before replacement;
6. simple-reference versus packed/indexed/optimized contract equivalence where
   an executable reference provides sufficient additional value;
7. serial versus threaded and, where applicable, memory-loaded versus mapped
   determinism or governed numerical equivalence;
8. artifact existence, dimensions, sampling, and content checksum/statistical
   comparison through `compare_numerical_artifacts` and the reviewed tolerance
   registry;
9. uninterrupted versus restart parity;
10. shared versus distributed parity;
11. stream/chunk/merge/prune compatibility;
12. representative particle-count memory and throughput comparison against the
    accepted M0/M1 baseline.

### 10.4 Distributed and external-fixture accountability

Distributed, streaming, BOX, and real-data gates must not disappear as silent
skips. For each affected stage:

- use a deterministic local synthetic master/worker or partition/merge fixture
  where the contract can be exercised without a scheduler;
- record external fixture requirements and the named platform owner in the
  stage validation manifest;
- report unavailable platform tests as explicit outstanding gates, not passes;
- do not merge a replacement or close a stage that changes that execution path
  until the required target-platform output has been observed.

A stage that does not touch an optional execution path may mark its gate
`not_applicable` with a reviewed scope justification. It may not mark it
`passed` because the fixture was unavailable.

## 11. Staged commits on `master`

The identifiers below describe independently committable milestones, not one
large feature branch. If a milestone grows beyond a reviewable diff, split it
further while preserving the same entry and exit gates.

### M0: freeze contracts and establish baselines

Changes:

- add this note;
- approve and freeze the foundational public names required by M2 while keeping
  unmigrated aggregate and analysis terminology provisional;
- record an adopt/adapt/reject decision for the X-derived patterns and approve
  the initial bounded-context `must not know about` rules;
- inventory the production readers/writers for the first bounded context;
- identify the first domain-owned contracts and the existing legacy APIs that
  will implement them, without adding production abstractions yet;
- identify or add the characterization gates for the first-slice legacy
  operations and decide where an executable reference provides enough
  additional value to retain;
- capture deterministic accepted semantic and numerical baselines before
  changing production paths;
- select fixed performance workloads and capture initial unchanged-`master`
  wall-time, peak-memory, and throughput results on named environments;
- identify required external fixtures and platform owners.

Exit gate:

- baseline generation and baseline comparison are distinct;
- accepted snapshots come from unchanged `master` behavior;
- every foundational M2 type name is either approved or explicitly deferred;
- initial dependency exclusions, contract ownership, characterization gates,
  and reference-oracle decisions are recorded;
- numerical artifacts have proposed metrics and an identified tolerance owner;
- performance workloads, environments, and owners are recorded even if M1 is
  needed to automate repeatable capture;
- no production behavior changes.

### M1: add project lifecycle fixtures and validators

Changes:

- add `project_state_fixture`, `validate_project_state`, semantic snapshot, and
  comparison support;
- add `compare_numerical_artifacts` and the reviewed, file-backed
  `production/tests/domain_refactor_numerical_tolerances.tsv` registry;
- add the repeatable performance-baseline harness and record accepted
  unchanged-`master` references for the M0 workloads;
- add common characterization/reference-equivalence support that can compare
  accepted semantic results with packed/indexed, serial, and parallel
  implementations without using the implementation under test to construct
  expected results;
- make representative performance workloads self-verifying against their
  semantic or accepted legacy result before recording timings;
- add comparison manifests that classify restart-sensitive and distributed-
  sensitive fields;
- strengthen existing project tests where they currently print without
  asserting;
- cover empty, acquisition, extracted, classified, refined, sparse-pruned,
  merge-collision, merged, and round-tripped states incrementally.

Exit gate:

- validator accepts known-valid fixtures and rejects deliberate single-fault
  fixtures with precise diagnostics;
- semantic comparison fails on unclassified fields;
- numerical comparison fails on missing tolerance rules and detects deliberate
  image, volume, stack-image, curve, and artifact-set perturbations;
- performance references include reproducible workload and environment
  metadata;
- a deliberate characterization/production mismatch and a serial/parallel
  mismatch are detected by the common comparison support;
- unchanged legacy project round trips pass;
- no production path consumes the validator yet.

### M2: introduce foundational value objects

Changes:

- add leaf modules for typed identifiers, quantities, geometry, and image
  locations;
- add constructors and unit tests;
- add a lightweight dependency check preventing those leaf modules from
  importing project persistence, UI, parameters, builder, commander, strategy,
  or external-format modules; consider a separate Fortran module directory or
  leaf library only if it can be introduced without duplicating the complete
  SIMPLE build;
- do not route production workflows through them yet.

Exit gate:

- invalid values are rejected;
- canonical units and conversions are tested;
- committed public type names match the M0 terminology decision;
- no dependency from leaf domain modules to `parameters`, `builder`,
  `sp_project`, commanders, strategies, UI, or external file-format modules;
- the dependency check fails on a deliberately planted forbidden import.

### M3: split and adapt the CTF domain contract

Changes:

- introduce typed CTF model, estimate, application, and fitting-policy types;
- define the read-only CTF resolution contract in the domain layer;
- add lossless conversions between the split domain contract and `ctfparams`,
  plus lossless project-row mappings for only the persisted subset;
- characterize direct legacy CTF resolution, then implement the contract with a
  read-only `sp_project` compatibility adapter over the existing authority
  rules; direct legacy resolution may be replaced in the same focused commit
  when the characterization gate exercises the replacement;
- source `ctf_fitting_policy` only from validated runtime configuration and
  exclude it from project row conversion and project round-trip claims;
- retain `ctfparams` as a compatibility DTO while unmigrated callers still
  require it; the existing production resolver need not remain after its
  characterized replacement is complete.

Exit gate:

- exact field and phase-convention parity for mic, stack, `ptcl2D`, `ptcl3D`,
  and the synthetic no-CTF `cls3D` case;
- unsupported oritypes exercise and preserve the legacy error path rather than
  returning a default domain object;
- mixed stack-level CTF models and per-particle estimates pass;
- contract-based resolution matches the independently captured legacy
  characterization for every supported source type and failure case;
- runtime fitting-policy configuration parity passes independently of project
  persistence, and no test claims that fitting policy round-trips through a
  project;
- numerical `simple_test_ctf` behavior is unchanged.

### M4: add typed acquisition views

Changes:

- add `Movie` and `Micrograph` views and stable typed identities over `os_mic`;
- define domain-owned acquisition read contracts, including bulk views where a
  consumer would otherwise perform repeated string-key lookups;
- represent direct-import and movie-derived micrographs explicitly in the view;
- add read and write adapters while preserving current rows and keys.

Exit gate:

- adapter characterization passes for movie-only, micrograph-only, and
  processed movie rows;
- unknown keys survive write-through;
- `.simple` semantic snapshots are unchanged.

### M5: migrate one preprocessing read path

Changes:

- select one narrow CTF or motion-preparation consumer;
- translate `parameters` into a use-case-specific immutable configuration;
- select and construct concrete scientific policies once at the commander or
  strategy setup boundary;
- consume the typed acquisition/CTF view;
- route any opt-in through UI metadata, generated argument sources, `cmdline`,
  `parameters`, validation, and the owning commander/strategy;
- replace the selected legacy consumer once its characterization gate passes;
  do not add a dual production path solely for the refactor.

Exit gate:

- command-line defaults and observable behavior remain unchanged; because this
  is a behavior-preserving replacement rather than a new feature, no opt-in key
  is required;
- the replacement produces the semantic project snapshots and accepted
  numerical outputs captured from the legacy implementation;
- the migrated operation does not read unrelated `parameters` fields or branch
  on method-name strings below its setup boundary;
- relevant in-memory, distributed, and restart checks pass.

### M6: migrate preprocessing writes and lineage

Changes:

- route motion, CTF, picking, and extraction updates through typed repositories
  one operation at a time;
- encode lineage without changing the current project format initially;
- centralize project-state validation at workflow boundaries in diagnostic or
  explicitly enabled mode before considering default enforcement.

Exit gate:

- every post-operation lifecycle state passes validation;
- typed writes match the independently captured legacy serialization contract;
- repositories expose domain operations and batch contracts rather than raw
  `ori` keys, segment identifiers, or storage offsets;
- extraction indexing, CTF propagation, optics propagation, and artifact paths
  remain unchanged.

### M7: introduce `ParticleSet` and typed indexing

Changes:

- add batch/columnar particle-set views over stack, `ptcl2D`, and `ptcl3D`
  segments;
- define a domain-owned `ParticleSet`/particle-batch read contract and implement
  it first with a compatibility view over the existing `oris` columns;
- use an independent clarity-first index/source oracle for deterministic
  fixtures when its diagnostic value justifies retaining it; otherwise compare
  the batch/columnar implementation with reviewed legacy characterization
  fixtures;
- centralize typed stack and physical-image identity;
- implement the reviewed per-field authority manifest rather than treating
  `ptcl2D` and `ptcl3D` as mirrors;
- add fixtures for every named 2D-to-3D and 3D-to-2D transfer operation,
  asserting the copied, transformed, cleared, preserved, and source-unchanged
  fields;
- migrate one read-only batch consumer before any write path.

Exit gate:

- dense, multi-stack, colliding-local-ID merge, merged, sparse-pruned, and
  transfer-semantics fixtures pass;
- 2D-owned and 3D-owned fields may differ without a false validator failure,
  while shared/index fields obey their declared synchronization rules;
- memory and throughput remain inside the reviewed bounds relative to the
  accepted M0/M1 baselines;
- direct and adapter indexing resolve the same physical images;
- the batch/columnar implementation matches the independent oracle for
  identities, order, selections, source mappings, and explicit failures;
- hot traversal performs no per-particle allocation, formatted key lookup, or
  polymorphic dispatch.

### M8: migrate 2D and 3D analysis one workflow at a time

Changes:

- introduce typed input configurations and result batches;
- define the operation's domain contract and establish an independent
  characterization or small reference oracle, especially for indexing,
  selection, transformation, search, ranking, and materialization boundaries;
- preserve commanders as application orchestration and strategies as execution
  policy;
- keep numerical search, filtering, masking, and reconstruction in their
  owning modules;
- materialize results back into the existing project state through adapters.

Exit gate for each workflow:

- initialization, active/rejected, classified/refined, multi-state, persisted,
  and restarted states validate;
- accepted default-path snapshot and numerical parity pass;
- oracle-versus-optimized and serial-versus-parallel equivalence pass for every
  applicable operation;
- shared/distributed parity passes where supported;
- the next workflow is not started until the current workflow is complete on
  `master`.

### M9: add analysis-run provenance additively

Changes:

- define `AnalysisRun` and typed artifact references;
- model `AnalysisRun` as an application-lifetime object that owns run-specific
  immutable configuration, policy selection, worker workspaces, diagnostics,
  and artifact/provenance coordination; it is not the aggregate root of all
  scientific state and must not replace bounded-context ownership;
- add provenance without replacing current-state segments;
- make any new persistent fields versioned and opt-in initially.

Exit gate:

- projects written without the opt-in remain semantically identical;
- old readers either ignore additive metadata safely or the metadata is stored
  outside the old schema;
- an accepted pre-M9 reader reads an opt-in-written project, preserves unknown
  additive fields through a read/write round trip where the old format supports
  them, and continues the supported legacy workflow without data loss;
- restart and artifact lookup work with and without provenance enabled;
- resource construction, worker ownership, and release order are explicit and
  covered by lifecycle tests.

### M10: encapsulate legacy internals only after migration

Changes:

- inventory remaining direct `os_*%get/set/isthere` access;
- inventory concrete persistence-adapter imports from migrated domain modules
  and eliminate them in favor of the approved contracts;
- remove or privatize direct access in small ownership-aligned commits after
  its characterization gate exists; the replacement and removal may occur in
  the same focused commit.

Exit gate:

- no supported consumer depends on the removed API;
- migrated domain modules pass the dependency-direction check;
- all lifecycle, parity, platform, and performance gates remain green;
- project-format redesign, if still desired, begins as a separate proposal.

## 12. Per-commit acceptance checklist

Before each production commit to `master`:

- confirm the current branch is `master` and inspect the dirty worktree;
- preserve unrelated user changes;
- identify the exact domain and persistence contracts touched;
- confirm contract ownership and run the dependency-direction check when a
  domain module is added or changed;
- add or update the narrowest test before or with the implementation;
- run the existing characterization gate against the legacy implementation, or
  add one before replacement; add an executable reference only where its
  diagnostic value justifies the maintenance cost;
- verify that expected results are not constructed through the implementation
  under test;
- run lightweight syntax/editor diagnostics;
- have the user or CI compile and run the required test set;
- inspect observed output and retain the relevant parity summary;
- verify that numerical tolerances and restart/distributed field manifests cover
  every compared output and field;
- compare performance-sensitive changes with the accepted workload baseline;
- require semantic/oracle equivalence to pass before accepting benchmark
  timing from an optimized, parallel, distributed, persisted, or mapped path;
- run `git diff --check`, review the complete diff, and stage only intended
  files;
- update the matching branch-context note with the verified milestone and
  outstanding platform/data tests;
- create one focused commit with no unrelated cleanup;
- verify `master` is left in a releasable state;
- push to `origin/master` only under the repository's explicit push workflow
  and never force-push.

A failed gate stops the stage. Do not commit a known regression with the intent
to repair it in the next stage.

## 13. Rollback policy

Every migration commit must be mechanically revertible without data migration.
This requires preserving the legacy persisted representation and avoiding
destructive on-read upgrades; it does not require preserving the legacy source
path. If a behavior-preserving typed replacement fails after reaching
`master`, revert the focused commit; do not patch project files in place as an
emergency compatibility measure.

Before an opt-in path is allowed to persist additive fields, its rollback test
must retain an executable built from the accepted pre-change `master` commit,
write a project with the opt-in implementation, and then use the old executable
to read, round-trip, and continue the supported legacy workflow. Unknown fields
must survive when the legacy storage machinery supports unknown keys; if they
cannot survive, the additive data must be stored outside the legacy project
schema or the stage must not proceed. This test is required again before M9
persists additive fields or changes their default behavior.

A behavior-preserving refactor may replace and remove its legacy source path in
one focused commit after the legacy characterization gate has been observed and
the replacement passes the unchanged gate. The recovery points are the
independent characterization evidence and mechanical commit reversion. A
separate opt-in/default-switch sequence is required only for genuinely new
behavior, incompatible persistence, or a change whose risk assessment
explicitly calls for a temporary dual path.

## 14. First recommended implementation slice

Begin with project lifecycle validation and CTF/optics typing, not with particle
analysis or persistence replacement.

Reasons:

- `ctfparams` is already a compact compatibility DTO;
- the numerical `ctf` type and tests provide an independent scientific oracle;
- CTF identity is currently duplicated across mic, stack, particle, and optics
  metadata;
- the existing microscope-parameter policy already documents authority and
  unresolved gaps;
- the slice exercises value objects, project adapters, mixed metadata
  resolution, persistence round trips, and workflow configuration without
  changing the core numerical algorithms.

The first production adoption should remain narrow and read-only. Typed writes
follow only after the lifecycle validator and semantic comparator can prove the
replacement is lossless. Each characterized legacy read or write path may be
replaced and removed in the same focused commit; no parallel production path is
required solely to perform the refactor.

## 15. Explicit non-goals for the initial program

- no wholesale rewrite of commanders or strategies;
- no mechanical entity inheritance hierarchy merely because concepts are
  related in the workflow;
- no requirement for one type per module, one library per bounded context, or a
  repository-wide build split; structural boundaries are introduced only when
  they enforce a valuable dependency rule at acceptable build cost;
- no object-per-particle heap model;
- no new `.simple` format;
- no replacement of optimized numerical kernels with domain wrappers;
- no simultaneous migration of preprocessing, 2D, 3D, stream, and FLEX paths;
- no analysis-event sourcing before typed current-state access is stable;
- no removal of dynamic extension keys until every supported producer and
  consumer has been inventoried;
- no claim that adding derived types alone constitutes DDD.

## 16. Completion criteria

The staged DDD migration is successful when:

- workflows express their inputs and results using stable domain language;
- units, identifiers, and cross-record relationships are checked by types and
  validators rather than caller convention;
- commanders no longer pass the entire global parameter and project surface to
  every scientific operation;
- scientific policies and immutable use-case configuration are selected once
  at the application boundary rather than rediscovered from strings in domain
  operations;
- migrated storage access occurs through domain-owned contracts with explicit
  batch capabilities, and domain modules do not import concrete project or
  external-format adapters;
- direct string-key mutation is confined to compatibility repositories;
- particle operations remain batched and performance-appropriate;
- applicable optimized implementations remain continuously checked against an
  independent characterization, analytical, snapshot, or executable-reference
  oracle;
- performance and determinism workloads verify results before reporting timing;
- every supported project lifecycle state has an executable validator and
  deterministic fixture or accepted real-data baseline;
- the current `.simple` representation can be changed or retained as an
  independent infrastructure decision;
- all changes arrived through small, validated, releasable commits on `master`.

## 17. Review disposition

The 2026-08-17 review is incorporated as follows:

| Review item | Disposition in this revision |
|---|---|
| Restart-sensitive `sampled`/`updatecnt` comparison | Added comparison modes, raw diagnostic retention, and mandatory per-workflow exact/normalized/invariant classification |
| Undefined numerical parity harness | Added `compare_numerical_artifacts`, required metrics, a file-backed tolerance registry, and approval ownership |
| Performance baseline captured too late | Added workload selection and initial capture to M0, automation/reference capture to M1, and baseline-relative M7 gates |
| Incomplete CTF parity matrix | Added synthetic `cls3D` no-CTF behavior and unsupported-oritype error parity to M3 |
| Runtime-only CTF fitting policy | Explicitly excluded it from project persistence and project round-trip claims |
| `ptcl2D`/`ptcl3D` authority ambiguity | Added a field-authority manifest, non-mirror semantics, and transfer-operation fixtures to M7 |
| Terminology freeze timing | M0 now freezes only the foundational public names required by M2; broader aggregate names remain provisional |
| Optional distributed/stream fixtures | Added synthetic fallback, named platform ownership, explicit outstanding gates, and a prohibition on silent-pass skips |
| M9 rollback after additive writes | Added old-executable read/round-trip/continuation validation before additive persistence or default changes |
| Missing UI/params route | Required all opt-ins to follow UI metadata, generated sources, `cmdline`, `parameters`, and owning workflow validation |
| Hot-loop AoS risk | Made compact intrinsic arrays mandatory in particle hot loops unless a benchmarked exception is approved |
| Merge collision coverage | Added colliding project-local identifier fixtures and exact-once remap checks |

## 18. X design-case-study disposition

The 2026-08-17 review of X was used as a design case study, not as code to copy
or as evidence that SIMPLE can be reorganized at the same cost. X is much
smaller and began with domain boundaries; SIMPLE must preserve mature project,
distributed, streaming, restart, and numerical contracts throughout migration.

| X pattern | SIMPLE disposition |
|---|---|
| Value objects with private representation, validated construction, canonical units, and meaningful operations | Adopt incrementally for identifiers, quantities, geometry, CTF contracts, locations, and specifications |
| Semantic objects over packed SoA/CSR storage | Adopt as the required particle-scale representation rule; typed boundaries lower to intrinsic bulk storage before hot loops |
| Domain-owned storage contract with slow and packed implementations | Adapt: define domain-owned contracts with explicit bulk capabilities, retain `sp_project`/`oris` adapters, and avoid concrete-adapter dependencies in migrated domain services |
| Complete slow implementation as an oracle for optimized/indexed behavior | Adapt: retain as test-only code only where its diagnostic value justifies maintenance; otherwise use independent characterization fixtures or analytical expectations |
| Policies built once and passed as immutable configuration | Adopt after the normal SIMPLE UI, generated metadata, `cmdline`, and `parameters` validation path |
| Per-worker preallocated workspace over shared read-only data | Adopt in workflow-specific run contexts and existing numerical owners; do not create an object per particle |
| Compile-enforced domain/infrastructure boundary | Adapt cautiously: begin with dependency checks and a possible leaf target; do not split the complete SIMPLE library as an opening move |
| Self-verifying serial/parallel/persisted benchmark | Adopt: correctness and determinism gates precede timing acceptance |
| Versioned, fail-closed packed artifact | Defer as a separate infrastructure decision; preserve `.simple` initially and require explicit compatibility approval for new formats |
| Shallow domain-specific inheritance | Do not generalize into a Movie/Micrograph/Particle hierarchy; prefer composition, services, result types, and lineage |
| Small-package command specification and one-type-per-module organization | Command specification: see `x_execution_layer_adoption.md`, an independent program that predates this refactor. One-type-per-module: do not transplant |

The X review also exposed a caution that this plan makes explicit: an abstract
semantic contract can become decorative if optimized consumers require methods
available only on a concrete packed type. SIMPLE contracts must therefore be
designed around the actual bulk access needed by production algorithms while
keeping storage layout and persistence ownership outside the domain API.

The X command-object, CLI-specification, and distributed-execution
mechanisms are adopted through a separate, independent program:
`doc/refactoring_notes/x_execution_layer_adoption.md`. It does not depend on
any milestone here and should land first.
