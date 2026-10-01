# FLEX PCA architecture and concrete refactoring plan

Date: 2026-09-18. Revised 2026-09-30.

Status: implementation plan; implementation not started.

Scope: the 3-D `flex_pca` implementation on `master` only. Feature branches that
carry their own FLEX changes are out of scope and are not reconciled here.

Audit baseline: `master` commit `6642aaefa856056bc08261800da582b7a3a0202f`
(2026-09-30). The original audit was taken at `1a68f25b6`; the FLEX commits between
the two (state-weights store, projfile-only input, concurrent even/odd PCG state
solves, per-state delivery, the CTest heterogeneity areas, warning cleanup, and the
removal of every flex GPU/CUDA path in `6642aaefa`) are folded into the measurements
below. FLEX is CPU-only on `master`; nothing in this plan provides for a device path.

This revision supersedes and replaces
`flex_pca_architecture_audit_and_refactoring_plan_2026_09_17.md`, which is deleted.

Validation level: static source inspection only. No compilation, test binary,
workflow, numerical-equivalence run, benchmark, or scientific experiment was run.

Line counts are physical `wc -l` counts, including comments and blank lines. Current
counts are exact. Future file sizes and patch sizes are planning estimates with an
expected uncertainty of about 25%.

This is the single living design record for the refactor. Update it as phases land
rather than creating companion architecture documents.

## 1. Motivation

FLEX has grown into a substantial subsystem. Its commander and
shared/master/worker strategy are useful boundaries, but below them workflow,
distributed transport, scientific state, reconstruction policy, configuration,
persistence, and numerical backends are interleaved.

Ordinary changes therefore cross unrelated lifecycles. The E-step mixes observation
preparation, Cartesian/polar statistics, posterior inference, and residual formation. The M-step mixes EM policy with gridding/PCG storage and transport.
State reconstruction mixes distribution, accumulation, solve, FSC/filtering, filenames,
and project delivery.

The refactor introduces the missing application and backend layers while preserving
the existing estimators, formats, and execution behavior. It is not a numerical
rewrite, performance project, or change to scientific defaults.

## 2. Preservation contract

Every core phase must preserve:

1. The public `flex_pca` command and current 3-D inputs and outputs.
2. One commander-owned `parameters`/builder lifecycle. A resolved FLEX configuration
   may be derived from `parameters`, but it is not another parser or source of defaults.
3. The same scientific phases in shared-memory, distributed-master, and worker runs.
4. Worker partials followed by master-owned reduction and publication.
5. Existing part names, magic values, versions, field order, byte layout, and atomic
   publication until an explicit schema migration.
6. The separate meanings of `rec_backend` and `rec_states_backend`.
7. Existing Cartesian/polar, gridding/PCG, paired-fit, cache/resume, sigma,
   support, filtering, convergence, and project-delivery behavior.
8. The `maxits_pcg=0` gridding/identity gate.

Ownership changes must not be combined with numerical cleanup, loop reordering,
precision changes, optimization, or artifact-format changes.

## 3. Measured current implementation

The complete `src/main/flex` tree is 22,781 lines. Including the FLEX commander and
parallelization strategy gives 23,566 lines. Since the original audit the CUDA
kernels, the GPU façade, and the device E-step and state-accumulation paths were
removed, which took about 6,700 lines out of the numerical backend, 740 out of
`em_iter`, and 300 out of `rec3D`. The tree also gained two owners the first draft did
not know about: the 690-line state-weights store (`simple_flex_weights_state`) and the
two tester modules (341 lines).

| Non-overlapping area | Exact LOC | Architectural treatment |
|---|---:|---|
| Commander and topology strategy | 785 | Keep; narrow their dependencies |
| `model` and `rounds` | 2,130 | Replace the model driver with application services; narrow rounds |
| EM parent and all EM submodules | 9,137 | Split state, E-step, posterior, M-step, transport, and drivers |
| State reconstruction, gridding plus PCG | 1,137 | Split service, backend, transport, and delivery |
| Coupled PCG and latent operators | 3,626 | Keep numerics; add narrow adapters first |
| Other FLEX scientific/support modules, weights store, testers, README | 6,751 | Retain their existing numerical ownership |

The ownership-heavy refactoring surface is therefore **13,189 exact lines**: commander,
strategy, model, rounds, EM family, and state reconstruction. The 3,626-line numerical
backend surface is not part of the initial rewrite.

### 3.1 Control and application concentration

| Current source | Exact LOC | Concentration |
|---|---:|---|
| [`simple_commanders_flex_pca`](../../src/main/commanders/simple/simple_commanders_flex_pca.f90#L40-L364) | 366 | 35-line entry; defaults, sigma bootstrap, and project-derived geometry occupy the rest |
| [`simple_flex_pca_strategy`](../../src/main/strategies/parallelization/simple_flex_pca_strategy.f90#L26-L417) | 419 | Factory, three roles, qsys scheduling, partitions, and worker lifecycle |
| [`simple_flex_pca_rounds`](../../src/main/flex/simple_flex_pca_rounds.f90#L17-L225) | 227 | Executor contract, stage IDs, half policy, schemas, names, and global part directory |
| [`simple_flex_pca_model`](../../src/main/flex/simple_flex_pca_model.f90#L82-L1901) | 1,903 | 534-line application driver, worker dispatch, caches, state inference, delivery, tests, and cleanup |

The commander entry is already thin. The missing boundary is below the strategy: the
534-line [`run_flex_pca`](../../src/main/flex/simple_flex_pca_model.f90#L82-L615)
is the application service but is named and structured as a scientific model.

`rounds` is a valid distributed port, but it currently combines four concerns:
execution, scientific stage identifiers, artifact schema identifiers, and mutable path
state. Role/scheduling calls occur 54 times across ten FLEX files; the target is to
confine them to the application, fit engine, and state-reconstruction service.

### 3.2 Fit concentration

The most coupled fit surface is 5,763 lines:

| Current source | Exact LOC | Main issue |
|---|---:|---|
| [`simple_flex_pca_em`](../../src/main/flex/simple_flex_pca_em.f90#L145-L1078) | 1,078 | 125-line `probe_fit_t`; 781-line parent interface |
| [`simple_flex_pca_em_state`](../../src/main/flex/simple_flex_pca_em_state.f90#L12-L186) | 188 | One lifecycle frees every unrelated state family |
| [`simple_flex_pca_em_iter`](../../src/main/flex/simple_flex_pca_em_iter.f90#L23-L1027) | 1,027 | 412-line single-fit driver, stage config, iteration begin, plus paired/master/worker drivers |
| [`simple_flex_pca_em_estep`](../../src/main/flex/simple_flex_pca_em_estep.f90#L25-L1408) | 1,408 | E-step work through line 665; 743-line part codec after it |
| [`simple_flex_pca_em_mstep`](../../src/main/flex/simple_flex_pca_em_mstep.f90#L28-L491) | 491 | One 450-line gridding/PCG finalizer |
| [`simple_flex_pca_em_solve`](../../src/main/flex/simple_flex_pca_em_solve.f90#L40-L527) | 529 | Shared posterior/MCFA algebra without an owning service |
| [`simple_flex_pca_em_pairmerge`](../../src/main/flex/simple_flex_pca_em_pairmerge.f90#L56-L1038) | 1,042 | Fit-frame merge, backend snapshot algebra, and delivery |

`probe_fit_t` currently owns identity, the scientific model, convergence history,
mixture state, polar banks, per-iteration work, distributed payloads, gridding state,
PCG state, cross-FSC diagnostics, and merge snapshots. File-level submodules have not
created lifetime boundaries because they all mutate the same object.

### 3.3 State-reconstruction concentration

[`simple_flex_pca_rec3D`](../../src/main/flex/simple_flex_pca_rec3D.f90#L32-L648)
is 650 lines. Its 371-line main routine selects a backend, starts rounds, performs
weighted accumulation, reduces worker parts, finalizes maps, filters, and writes output.
[`simple_flex_pca_rec3D_pcg`](../../src/main/flex/simple_flex_pca_rec3D_pcg.f90#L55-L460)
adds 487 lines. Since the original audit its 406-line PCG reconstruction routine has
been reorganised into contained procedures (`deliver_state`, `new_state_operator`,
`accumulate_state_half`, `solve_state_half`, `report_state_half`,
`append_state_table`), solves even and odd halves concurrently on a capped thread
budget, and writes each state as soon as its halves are solved. The delivery tail is
now `deliver_state` at lines 205–279. This is the right internal shape for Phase 7 and
should be extracted as-is, not re-flattened.

The current reference 3-D workflow confirms the desired ownership direction:

- [`simple_strategy3D_matcher`](../../src/main/strategies/search/simple_strategy3D_matcher.f90#L104-L243)
  makes particle-side phase order explicit.
- [`simple_matcher_3Drec`](../../src/main/strategies/search/simple_matcher_3Drec.f90#L23-L107)
  produces partition-local reconstruction artifacts.
- [`simple_rec3D_strategy`](../../src/main/strategies/parallelization/simple_rec3D_strategy.f90#L34-L166)
  owns lifecycle, topology, and backend selection.
- [`volassemble`](../../src/main/commanders/simple/simple_commanders_rec_distr.f90#L871-L1235)
  owns master assembly, diagnostics, products, and metadata.

FLEX should copy the producer-to-assembler ownership split, not the inheritance shape or
the size of those files.

### 3.4 Environment-knob surface

**Correction (Phase 0b/1 measurement, 2026-09-30).** The figure of 50 keys below counted
literal occurrences of `SIMPLE_COV_*` / `SIMPLE_FLEX_*` names, including keys that only
appear in comments and log text. The actual read set on master 6642aaefa was 25 keys
(the `getenv`/`env_*` call sites). After Phase 1 every one of them is read exactly once,
in `flex_run_settings%new` (`simple_flex_pca_run_types`); no other flex module reads
the environment (the SIMPLE-wide `SIMPLE_PTCL_CACHE_DIR` in `plane_cache` is not a
flex knob). `doc/policies/flex_pca_policy.md`, referenced below and in Phase 0b, does
not exist; the remaining keys are documented in the settings type itself.

The tree (as first measured) reads **50 distinct `SIMPLE_COV_*` / `SIMPLE_FLEX_*`
environment keys**, from 21 files, all at the point of use. The GPU removal retired the thirteen `SIMPLE_COV_GPU_*`
and `SIMPLE_FLEX_GPU_*` keys; `simple_flex_pca_em_env` now resolves none of the rest.

| File | Distinct keys | Character |
|---|---:|---|
| `em_iter` | 13 | formulation selection (polar, paired), probe policy, deflation, pairing rule |
| `merge` | 8 | state-merge thresholds and link policy |
| `model` | 6 | composition, deconvolution, `d_tilde`, cross-FSC regularisation |
| `em_mstep`, `em` | 5 each | deflation, Wiener M-step, basis/probe caps |
| `pcg`, `em_pairmerge`, `em_fit` | 4 each | PCG lambda/warm start, deflation, calibration caps |
| 13 further files | 1–3 each | GMM/k-means choice, polar exactness, residency, contrast, state filtering, min-state |

Five keys are read from three or more files (`SIMPLE_COV_XFSC_REG`,
`SIMPLE_COV_MOD4_PAIRING`, `SIMPLE_COV_EM_DEFLATE`, `SIMPLE_COV_DTILDE`,
`SIMPLE_COV_PAIRED_MERGE`), so the same decision is currently re-resolved by
independent readers. The polar keys (`SIMPLE_COV_POLAR_*`) select a formulation inside
the iteration loop, which Section 4 rule 6 forbids.

This surface is the largest single dependency the "resolved settings" value of
Section 6.1 has to absorb. Classifying and retiring it is a behaviour decision, not a
mechanical move, so it is its own phase (Phase 0b) and must precede Phase 1.

### 3.5 Existing test coverage on `master`

FLEX tests run under CTest through the `heterogeneity` area
(`unit_heterogeneity`, `lib_heterogeneity`):

| Test | Exercises | Not exercised |
|---|---|---|
| `simple_flex_pca_tester` (unit) | embedding-cache I/O, auto settings, population-floor placement, kernel bandwidth, state weights | any E-step, M-step, probe part, or delivered map |
| `simple_flex_pca_tester` (lib) | deconvolution of 20,000 particles at realistic sizes | shared vs distributed, paired fits, resume |
| `simple_flex_pcg_tester` | coupled PCG operator check, optional solve and sweep | gridding backend, state PCG delivery |

None of these fixes an input and captures E-step sufficient statistics, M-step
payloads, probe-part bytes, state-weight tables, or delivered maps. Every
"reproduces its ..." gate in Sections 8 and 9 is therefore currently unverifiable.
Phase 0 has to create that evidence; it cannot reuse existing tests.

### 3.6 Conformance to SIMPLE object conventions

Measured against `.claude/skills/simple-modern-fortran` on the current baseline.

Conforms: `use ..., only:` on 212 of 233 imports (the rest are the core API aggregator
and test utilities); `flex_pca_rounds` is an abstract type with deferred procedures,
extended by the three strategies and allocated by a factory; `flex_pcg_t` is a full
object with 41 type-bound procedures and `new`/`kill`; the commander reads `cmdline`
only before the parameters object exists; the EM family uses the parent-plus-submodule
layout; no multi-line `THROW_HARD` forms.

Departs:

| Departure | Measured | Convention |
|---|---|---|
| Central object is a record | `probe_fit_t`: 171 components, 0 type-bound procedures; `new`/`kill` are free procedures; 15 submodules mutate it through dummy arguments | behaviour bound to the type (`reconstructor` 34 bindings, `oris` 216) |
| Types with allocatable state and no lifecycle | `polar_grid_t`, `crossfsc_record`, `crossfsc_file`, `xfsc_ctx_t`, `compose_src_t` are built and freed by free `*_build`/`*_setup`/`*_kill`/`*_teardown` procedures | `new`/`kill` symmetry on the type |
| Module-level singletons | resident plane store (10 module variables), plane cache (8), PCG envelope images (2 `save`d images plus flags), model crop globals and UMAP arrays (4) | state on a created-and-killed object |
| Optional dummies every caller supplies | 34 optional arguments across 20 procedures are passed at every call site (the paired-fit entry passes all five merge outputs on its single call; latent deconvolution has nine such arguments) | explicit dependencies; `optional` only for genuinely optional inputs |
| Parameters inside numerical leaves | `simple_flex_pca_pcg` takes `parameters`/`builder` in 11 signatures and reads 8 fields; `weights` and `latent_ops` do the same on a smaller scale | orchestration values stay out of domain modules |

Value types without allocatable state (`flex_pcg_outcome_t`, `half_solve_rec`) are
records by design and are not departures.

## 4. Target dependency structure

```text
commander
  -> topology strategy: shared | qsys master | worker
     -> flex_pca_application
        |-- resolved run settings + run session
        |-- project gateway + artifact catalog
        |-- rounds / stage executor
        |-- fit engine
        |    |-- E-step backend -> sufficient statistics
        |    |-- shared posterior inference
        |    `-- M-step backend: gridding | coupled PCG
        |-- embedding and state-inference services
        |-- state-reconstruction service
        |    |-- backend: gridding | PCG
        |    `-- common 3-D delivery
        `-- numerical leaf modules
```

Dependencies point downward only:

1. The application sequences scientific phases.
2. The strategy and rounds object own process role, partitioning, launch, and wait.
3. A service may request a stage; numerical backends may not query process role.
4. The project gateway is the only owner of project persistence.
5. The artifact catalog owns identity, paths, and provenance; producer-specific codecs
   own magic values, versions, shapes, and byte layout.
6. Backend selection and fallback occur at stage or iteration boundaries, never inside
   a particle hot loop.
7. Numerical leaves receive prepared values, not `cmdline`, qsys objects, filenames,
   project mutation access, environment policy, or the whole application session.
8. A dummy argument is `optional` only when some valid caller omits it. Every one of the
   34 always-supplied optionals in Section 3.6 becomes a required argument in the phase
   that moves its procedure; new interfaces do not introduce optionals for outputs.

## 5. Concrete source organization

Keep `src/main/flex` flat, matching the existing subsystem convention. Prefixes provide
the grouping without a repository-wide directory move.

Target LOC below describes the intended steady-state size, not an acceptance criterion.

| Target owner | Target file | Estimated LOC | Consolidates or abstracts |
|---|---|---:|---|
| Application | `simple_flex_pca_application.f90` | 300–450 | Master/shared phase order and worker-stage dispatch from `model` |
| Run state | `simple_flex_pca_run_types.f90` | 350–500 | Read-only derived settings, run-owned scientific values, resources, and one `kill` path |
| Project boundary | `simple_flex_pca_project_gateway.f90` | 250–350 | Selection validation, EO repair, labels, project FSC reads, and final segment writes |
| Artifact identity | `simple_flex_pca_artifacts.f90` | 250–350 | Canonical names, run/part roots, provenance, and atomic publication policy |
| Embedding persistence | `simple_flex_pca_embedding_io.f90` | 220–320 | Raw embedding cache and deconvolution block I/O from `model` |
| Fit façade | `simple_flex_pca_fit_engine.f90` | 650–850 | One single/paired/worker iteration engine over `fits(:)` |
| Fit state | `simple_flex_pca_fit_types.f90` | 400–550 | Composed fit specification, model, history, E/M sessions, diagnostics, and iteration workspace |
| E-step | `simple_flex_pca_estep.f90` | 600–800 | Cartesian/polar lifecycle and one batch statistics contract |
| Posterior | `simple_flex_pca_posterior.f90` | 500–550 | ECM/MCFA solve, contrast, moments, likelihood, and Gamma accounting |
| M-step | `simple_flex_pca_mstep.f90` | 500–650 | Backend façade plus gridding/PCG accumulation, reduction, ridge, solve, and snapshot lifecycle |
| Fit transport | `simple_flex_pca_probe_parts.f90` | 700–850 | Typed probe payload and current v5/legacy codecs |
| State inference | `simple_flex_pca_state_service.f90` | 400–550 | Deconvolution, placement, occupancy pruning, CV, merge decision, and reconstruction request |
| State reconstruction | `simple_flex_pca_state_reconstruction.f90` | 180–230 | Backend selection, executor interaction, and one use-case order |
| State gridding | `simple_flex_pca_state_gridding.f90` | 400–480 | Weighted accumulation, worker partial, master reduction, map finalization |
| State PCG | `simple_flex_pca_state_pcg.f90` | 270–320 | Existing state PCG accumulation, reduction, and solve |
| State transport | `simple_flex_pca_state_parts.f90` | 120–160 | Weight table and gridding/PCG part names, codecs, and provenance |
| 3-D delivery | `simple_flex_pca_delivery_3d.f90` | 300–450 | State FSC/filter/mask, filenames, map publication, tables, manifest, and project registration |

The existing commander and strategy remain. `simple_flex_pca_rounds.f90` remains the
distributed port but loses artifact schemas, path globals, and scientific split policy.
`simple_flex_pca_model.f90` is retired after a temporary compatibility wrapper.
`simple_flex_pca_em.f90` becomes a small compatibility façade while callers migrate.

The following existing `em_*` submodules survive as they are, with only their
orchestration and I/O call sites redirected. They are **not** extraction targets and
must not be decomposed in Phases 3–5:

- `simple_flex_pca_em_basis.f90` (341), `simple_flex_pca_em_mean.f90` (577),
  `simple_flex_pca_em_embed.f90` (595), `simple_flex_pca_em_compose.f90` (450),
  `simple_flex_pca_em_crossfsc.f90` (338), `simple_flex_pca_em_polar.f90` (392).
- `simple_flex_pca_em_env.f90` (250) is absorbed into the resolved-settings factory.

The following remain numerical leaves in the core refactor:

- `simple_flex_pca_pcg.f90`: coupled covariance PCG.
- `simple_flex_reconstructor_latent_ops.f90`: projection, backprojection, and coupled
  gridding kernels.
- `simple_flex_pca_polar.f90`: polar numerical kernels.
- `simple_flex_pca_weights`, `targets`, `gmm`, `deconv`, and `merge`: state mathematics.
- Existing mean, basis, embedding, composition, cross-FSC, and pair-merge numerical
  routines, after their orchestration and I/O dependencies move outward.

## 6. Exact consolidations and abstractions

### 6.1 Application and lifetime

The application replaces the current 527-line master/shared driver and 39-line worker
dispatcher. It has six explicit operations:

```text
prepare
obtain_embedding
infer_states
reconstruct_states
publish
execute_worker_stage
```

`flex_pca_run_types` contains two different values:

- `flex_run_settings`: an immutable projection of the already-parsed `parameters`, plus
  facts not retained there such as whether a key was explicit and resolved legacy
  environment switches. It cannot set defaults or mutate `parameters`.
- `flex_run_session`: owns selection, mean/basis, embedding, state values, cache/plane
  lifecycle, and teardown. It replaces application-level module globals and the long
  local allocation list in `run_flex_pca`.

### 6.2 Rounds and artifacts

Keep the existing abstract rounds object, but replace the primitive optional argument
bundle with a typed request:

```text
flex_stage_request
  stage
  fit selection
  iteration / iteration budget
  fit count
  payload kind
```

Move stage identifiers into the stage protocol, the mod-4 rule into fit selection
policy, and paths into the artifact catalog. `rounds` then owns only role, `nparts`,
partition planning, launch, and wait.

The artifact catalog does not become a giant serializer. These codecs stay separate:

- probe parts;
- embedding sufficient-statistics parts;
- mean-scale state;
- state weights and reconstruction parts;
- the per-state weights store (`simple_flex_weights_state`, 690 lines, with its
  `simple_flex_weights_file` reader in `src/fileio`): persistence plus placement
  policy today, split into a codec and a state-service call in Phase 2;
- embedding cache.

Each codec preserves its current schema and validates a typed payload before I/O.

### 6.3 Fit state

Replace the flat `probe_fit_t` with composition:

```text
flex_fit
  spec          identity, selection, artifact keys, resolved policy
  model         mean, basis, eigenvalues, rank, sigma
  history       convergence and persistent MCFA parameters
  estep         selected backend and stage/iteration resources
  mstep         selected system and paired-merge snapshot
  diagnostics   cross-FSC and timings
  iter          per-iteration workspace
```

The distributed payload is not fit state. `flex_probe_part` is a separate value passed
to its codec. This removes the current 20-plus-argument write/read/fold calls and makes
the lifetime of each allocation explicit.

Each composed piece is a SIMPLE object, not a smaller record: it has `new` and `kill`
bound to the type, and the procedures that today take the flat `probe_fit_t` as a
dummy argument become type-bound procedures on the piece that owns the state they
mutate. `flex_fit` itself binds `new`, `kill`, and the per-iteration `begin`/`finish`
pair, and its `kill` calls the pieces' `kill` in reverse order of construction. The
same rule applies to the existing lifecycle-free types: `polar_grid_t`, the cross-FSC
record and file, `xfsc_ctx_t`, and `compose_src_t` gain `new`/`kill` bindings when
their owner module is touched, and their free `*_build`/`*_kill` procedures become
those bindings. The module-level singletons (resident plane store, plane cache, PCG
envelope images, model crop globals) become fields of `flex_run_session` or of the
E-step/M-step session that owns their lifetime.

### 6.4 E-step backend

The E-step backend is a first-class boundary. It is a composed, tagged façade rather
than a Cartesian × polar inheritance pair, so that a further formulation (a nuisance
block, a 2-D class-frame E-step) is one more implementation of the same contract.

```text
select_backend
begin_stage
begin_iteration
form_statistics_batch
form_residual_batch
end_iteration
kill
```

All formulations produce one iteration-owned batch value containing the common
sufficient statistics `G`, `b`, `c`, `e_mm`, and `myv`. It is allocated once and reused;
there is no derived-type allocation per particle. Backend selection happens once, at
`begin_stage`, never inside the batch loop.

### 6.5 Posterior inference

One posterior service consumes sufficient statistics and owns:

- ECM versus MCFA selection;
- fitted contrast;
- latent mean and covariance;
- `E[zz']`, likelihood, and Gamma accounting;
- scaling of moments passed to the M-step.

This consolidates `fit_estep_solve_stats` with the posterior tails repeated in the
Cartesian and polar formers and in the paired pass. Projection code cannot choose a
posterior model.

### 6.6 M-step backend

The M-step contract is:

```text
begin_iteration
accumulate_batch
reduce_local
export_part / import_part
apply_ridge
solve_halves
snapshot_for_pair_merge
kill_iteration
```

The gridding implementation continues to call the existing coupled insertion and solve
kernels. The PCG implementation continues to use `flex_pcg_t`. The abstraction owns
their storage and lifecycle; it does not merge their algebra.

Single, paired, and worker execution use the same E-step, posterior, and M-step services.
The paired layer retains only particle-to-fit ownership, two-fit convergence, cross-fit
diagnostics, and final frame merge.

### 6.7 State reconstruction and delivery

The state-reconstruction service selects gridding or PCG and asks the executor for a
stage. Each backend implements:

```text
begin
accumulate_local_or_write_part
fold_parts
finalize_maps
kill
```

It returns a typed map bundle containing combined/even/odd maps and backend diagnostics.
The common 3-D delivery owner then performs FSC, filtering/masking, naming, publication,
and project updates.

This consolidates the gridding delivery tail at the end of
`reconstruct_flex_weighted_states` and the PCG delivery tail (`rec3D_pcg:205–279`,
`deliver_state`). Current backend differences remain explicit delivery-policy
values; they are not silently normalized.

## 7. Current-to-target move map

| Measured current span | Lines | Destination | Treatment |
|---|---:|---|---|
| `model:82–659` | 578 | application + fit/state services | Split phase order from worker dispatch; preserve call order |
| `model:663–833` | 171 | run settings/session + project gateway | Move validation, memory report, and sigma preparation |
| `model:839–1149` | 311 | state service + embedding I/O | Separate deconvolution from persistence |
| `model:1156–1532` | 377 | state service + weights store + delivery/project gateway | Separate state decisions from publication; `write_flex_weights_store` goes with the weights codec |
| `model:1572–1871` | 300 | state service + delivery | Move population-floor placement and readouts |
| `rounds:17–225` | 209 | rounds + stage protocol + artifacts | Keep executor; split policy/schema/path state |
| `em:145–269` plus `em_state:12–186` | 300 | fit types | Replace flat state and monolithic teardown with composed lifetimes |
| `em:298–1078` | 781 | narrow owner modules | Replace one parent-wide interface with owner-local interfaces |
| `em_iter:23–434` | 412 | fit engine | Decompose the single-fit driver |
| `em_iter:435–637` | 203 | fit engine + run settings | Stage config and iteration begin |
| `em_iter:638–1027` | 390 | fit engine | Paired, master, and worker drivers |
| `em_estep:25–665` | 641 | E-step + posterior + fit engine | Preserve formers, solve, reduction, and paired pass |
| `em_estep:666–1408` | 743 | probe parts | Move codec unchanged before simplifying it |
| `em_mstep:28–491` | 464 | M-step | Split policy from gridding/PCG state |
| `em_solve:1–529` | 529 | posterior | Retain numerical routines; consolidate callers |
| `rec3D:32–402` | 371 | state service + gridding backend + delivery | Separate execution, accumulation, and delivery |
| `rec3D:406–648` | 243 | state parts + project gateway + gridding backend | Separate transport, project FSC, and reconstructor init |
| `rec3D_pcg:55–460` | 406 | state PCG + state parts + delivery | Keep the contained-procedure split; `deliver_state` (205–279) moves to 3-D delivery, the rest to the PCG backend |


Moving a range does not authorize rewriting it. The first commit for each range should
be a mechanical extraction with identical code and an unchanged call path.

## 8. Implementation sequence and evidence

Line budgets have been removed from this section. Section 9 already states that no
phase passes because a file became shorter, and the earlier moved/new/removed
estimates measured nothing that a gate depends on. The rough expectations remain:
roughly 6,000–8,000 relocated lines, net production growth between flat and +1,200
lines, and a Git diff of 15,000–20,000 added-plus-deleted lines because relocated code
counts on both sides.

Each phase instead names the evidence it must produce. "Fixture" means a captured
reference output for a fixed input, stored under the `heterogeneity` CTest area, that
later phases compare against bit-for-bit (integer/byte payloads) or to a stated
tolerance (floating-point statistics).

| Phase | Implementation | Evidence required to pass |
|---:|---|---|
| 0a | Parity harness | Fixtures F1–F8 below exist and pass on the unmodified baseline; the harness is registered in `unit_heterogeneity` (fast subset) and `lib_heterogeneity` (full) |
| 0b | Environment-knob inventory and retirement | Every one of the 50 keys in Section 3.4 is classified as default-path, experiment-only, or dead; dead keys and their code paths are removed; the remaining keys are listed in `doc/policies/flex_pca_policy.md`; F1–F8 unchanged |
| 1 | Run settings/session, stage request, artifact catalog | Every remaining key from 0b is read exactly once, in the settings factory; no format or phase-order change; F1–F8 unchanged |
| 2 | Application, project gateway, state service, weights codec, 3-D delivery behind compatibility calls | `model` no longer owns phase order; project writes have one owner; F5–F8 unchanged |
| 3 | Composed fit state and typed probe payload/codec | Lifecycle test proves iteration cleanup cannot free persistent state; every composed piece has `new`/`kill` bindings and its mutators are type-bound; F3 probe-part bytes identical |
| 4 | Posterior service and E-step backend/session | F1 and F2 statistics and posterior moments identical for the Cartesian and polar paths |
| 5 | One fit engine for single, paired, and worker fits | F4 rounds, convergence trace, fit products, and paired-merge inputs identical |
| 6 | M-step backend contract over gridding and PCG storage | F4 basis, part payloads, ridge, and stopping behaviour identical per backend |
| 7 | State reconstruction split: service, gridding/PCG backends, parts, delivery | F6–F8 maps, FSC/filter products, filenames, and project metadata identical per backend; only the master publishes |
| 8 | Remove compatibility façades and stale documentation | Mechanical checks in Section 9 pass; no dead entry point; README and policy describe only implemented owners |

Fixtures:

| Id | Configuration | Captured reference |
|---|---|---|
| F1 | shared-memory, Cartesian, gridding M-step, single fit, 2 iterations | E-step sufficient statistics per batch (`G`, `b`, `c`, `e_mm`, `myv`), posterior moments, contrast, likelihood |
| F2 | as F1 with polar E-step | same as F1 |
| F3 | distributed master + 2 workers, F1 settings | probe-part bytes from each worker, master fold result, embedding sufficient-statistics parts |
| F4 | paired fit, gridding and PCG M-step (two captures) | per-iteration basis, ridge, convergence trace, paired-merge snapshot |
| F5 | resume from embedding cache | embedding cache bytes, deconvolution block, state labels |
| F6 | state reconstruction, gridding backend, shared | combined/even/odd maps, FSC, filtered maps, weight tables, manifest |
| F7 | state reconstruction, PCG backend, shared | as F6 plus per-state PCG table |
| F8 | state reconstruction, distributed | as F6 plus per-part state weights and gridding parts |

Capturing F3, F4, and F8 needs cluster time; budget it before starting Phase 0a.

### 8.1 Optional late physical splits

After the new contracts are stable, two large numerical leaf files can be split
mechanically without changing their public behavior:

| Leaf | Exact current LOC | Proposed physical split | Additional moved LOC | New scaffolding |
|---|---:|---|---:|---:|
| `simple_flex_reconstructor_latent_ops` | 1,011 | projection, backprojection, coupled grid solve, observation preparation | 950–1,020 | 60–110 |
| `simple_flex_pca_pcg` | 2,615 | geometry, accumulation, operator, solve, support, test support submodules | 2,250–2,450 | 250–400 |

This optional phase adds 3,200–3,470 moved lines and 310–510 interface lines.

## 9. Validation and completion rules

Compilation and runtime gates are user-run unless separately authorized.

| Gate | Pass condition |
|---|---|
| Static ownership | Numerical leaves contain no qsys/role checks, project mutation, final filenames, or environment-policy selection |
| Lifecycle | Each stateful owner has symmetric construction/teardown; success, early return, fallback, and worker exit have one cleanup path |
| E-step | Cartesian and polar paths preserve common sufficient statistics, posterior moments, contrast, and likelihood |
| M-step | Gridding and PCG each preserve accumulation shapes, part payloads, ridge, solve, paired snapshot, and stopping behavior |
| Distribution | Shared and distributed executions preserve phase order, worker inputs, reductions, artifacts, and master-only publication |
| State reconstruction | Gridding and PCG each preserve combined/even/odd maps, FSC/filtering, support, filenames, and project metadata |
| Artifacts and restart | Core phases retain exact formats; wrong shape/version/provenance fails before folding or adoption |
| Documentation | Policy and README describe only implemented owners; this note records completed phases and outstanding gates |

Mechanical checks at completion:

- `cmdline` appears only in commander/bootstrap and the application boundary.
- Environment reads appear only in the resolved-settings factory, and every key read
  is listed in the policy document.
- No environment key retired in Phase 0b is read anywhere.
- Project writes appear only in the project gateway.
- Final output names appear only in artifact/delivery owners.
- Only the application, fit engine, and state-reconstruction service depend on rounds.
- Backend selectors do not appear in particle loops.
- No derived type in `src/main/flex` holds allocatable state without `new`/`kill`
  bindings; no module-scope variable in the tree is mutable (parameters only).
- No `optional` dummy in the tree is supplied by every caller (checked by the script
  that produced Section 3.6).

No phase passes because a file became shorter. It passes only when the claimed ownership
movement and all applicable behavior-preservation gates have evidence.

## 10. Non-goals

- No big-bang rewrite.
- No requirement that gridding and PCG produce the same numerical result.
- No common superclass for scalar reconstruct3D PCG and coupled covariance PCG.
- No topology × formulation × solver subclass matrix.
- No speculative generic-dimensional reconstructor.
- No scheduler, cache, filename, or artifact-schema redesign during extraction.
- No silent reinterpretation of old artifacts.
- No performance or scientific-quality claim without separate matched evidence.

As phases land, update this note and `src/main/flex/README.md`. Do not describe
proposed interfaces as implemented before they exist.

## 11. Implementation status (branch `flex-pca-refactor`, from master 6642aaefa, 2026-09-30)

| Phase | Commit | Status |
|---:|---|---|
| 0a | -- | **Not done.** The fixtures F1--F8 need cluster runs; every "identical" gate below is therefore unverified until the harness is captured on `master` and replayed on this branch. |
| 0b | `bfc5e4cb0` | Dead procedures (17), dead constants (13) and dead flags removed. The knob inventory is the corrected Section 3.4. |
| 1 | `58c129dba`, `306d5930c` | `flex_run_settings` (the one environment reader), `flex_run_session`, `flex_stage_request`, the artifact catalog; `rounds` is executor-only. |
| 2 | `5a396f852`, `39f3ba079` | `flex_pca_application` (prepare / obtain_embedding / infer_states / reconstruct_states / publish, worker entry), project gateway, state service, embedding I/O, 3-D delivery; `simple_flex_pca_model` deleted. |
| 3 | `551539364`, `8018e00e5` | `flex_fit` composed of seven pieces with `kill` each; `flex_probe_part` typed payload, borrowed/restored without copies; v12 and v5 codecs over the payload. |
| 4 | `c8e4f28ab`, `f41bb3375` | `simple_flex_pca_fit_types`, `simple_flex_pca_posterior`; the E-step contract bound to `flex_probe_fit` (`estep_begin_stage`, `estep_bank_prepare`, `estep_batch_begin`, `estep_particle`, `mstep_insert_batch`); `em_state` and `em_solve` deleted. |
| 5 | `a0c5fa125` | `fit_engine_iterate`: one master loop for the single, paired and worker fits over `fit_estep_pass` (k-way merged read list); `probe_subspace_paired`, `paired_estep_pass` and the whole-batch former gone. |
| 6 | `943aa763a` | `simple_flex_pca_mstep`: `flex_fit_mstep` owns gridding/PCG storage behind `begin_iteration` / `accumulate_batch` / `reduce_local` / `apply_ridge` / `solve_halves` / `snapshot_for_pair_merge` / `kill_iteration`. |
| 7 | `05f932fd4` | `flex_states_backend` (abstract; gridding and PCG extensions), `flex_state_delivery` with per-backend `flex_state_delivery_policy`, `simple_flex_pca_state_parts`; `rec3D_pcg` deleted. |
| 8 | `0ede236bb`, `4f12493db` | Never-supplied optionals removed (`fprefix` of `build_covariance_eigenbasis`, `fsc_projfile`, `targets_in`/`zmetric`), always-supplied optionals made required (the paired merge products, `floor_rho`, `apply_ctf_amp`, `mskrad`, `ang_osamp`, `labels`, `proj_out`/`tproj_out`), unused dummies dropped, file-only publics of `simple_flex_pca_em` made private, README rewritten. |
| 9 | `f6006c5d4`, `1cec97bd2`, `efb2f18ac`, this commit | The EM regrouped the hybrid way: `simple_flex_probe_fit` (one parent, the type and its bindings, four submodules `_estep` / `_update` / `_engine` / `_crossfsc`, as `polarft_calc` is split) plus plain modules in dependency order `basis` < `embed` < `pairmerge` < `fit_driver`; `em_env` dissolved into `util`. 13 files became 9; every procedure body moved verbatim (checked per procedure against the previous commit); seeded toy A/B byte-identical to master. |
| 10 | `1214d4eff` .. `cdc006607` | Argument contracts and duplicates. Library calls replace exact flex duplicates (gcd, locate_2, euler2m). Four value records (`simple_flex_pca_records`: flex_selection, flex_fit_model, flex_latent, flex_state_set) compose the session and are what the services take: 23 procedures that took 8 to 24 loose arguments now take 2 to 11; the three posterior kernels take the fit and a thread index with verbatim bodies inside an associate block. `obtain_embedding` is a 61-line selector (resume / composition / paired fit) over one embedding tail; the paired merge's axis weighting and polish live in the paired driver; the unreachable single-fit branch and `build_covariance_eigenbasis` are gone. Shared-memory, PCG M-step, PCG states and gridding M-step arms byte-identical to master; the distributed arm agrees to the last ulp as before. |
| 11 | this commit | Folder layout: `src/main/flex/{run,fit,states}` plus the root for helpers and tests (a pure `git mv`; module names unchanged; the source glob is recursive). The file links in Sections 3 and 7 above are the pre-refactor measurement and keep their original paths. |

Deviations from the plan text worth knowing:

- The E-step backend is not a separate session object: the contract is bound to
  `flex_probe_fit` (Section 6.4's names), with the formulation chosen at
  `estep_begin_stage`. The single and paired drivers share `fit_estep_pass`.
- The single fit iterates against the session's mean by reference
  (`flex_mean_ref`), the paired fits against their own scaled copies; the engine
  never copies a reconstructor.
- The PCG state backend folds a state's raw parts when its maps are finalized
  (one operator pair resident), so `fold_parts` is a no-op for it; the gridding
  backend folds every state's parts up front, as before.
- Delivery differences between the backends are the declared
  `flex_state_delivery_policy` values (mask, project-FSC fallback, the run
  settings' eofilt/filt switches for PCG only, the log tag); nothing was
  normalised. Log lines keyed on the backend or on the fit count were kept.
- Flex helpers that look like library duplicates but are not drop-in and were kept:
  `corr_dp` and `median_of_sorted` (the library's `pearsn_serial_8` and `median` return
  single precision), the 8th-order Butterworth low-pass of the state delivery (the
  library's `butterworth` is the analytic polynomial magnitude), `mem_available_gb` (the
  library has disk space, not memory), the UMAP explicit-state RNG (the library generator
  is global state), the tester's `gauss` (single precision on the global generator), and
  the three SPD solver families inside flex (posterior: scaled, ridged, retrying Cholesky;
  deconv and targets: plain Cholesky), which differ numerically and would change results.
- Build note: seventeen new sources (`run_types`, `stages`, `artifacts`,
  `embedding_io`, `state_service`, `delivery_3d`, `project_gateway`,
  `application`, `fit_types`, `posterior`, `mstep`, `state_parts`,
  `states_backend`, `states_gridding`, `states_pcg`, `state_delivery`) and five
  deletions (`model`, `em_state`, `em_solve`, `rec3D_pcg`, `em_env` helpers)
  mean every existing build directory needs its CMake configure re-run before
  `make` (the source glob is evaluated at configure time).
