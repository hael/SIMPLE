# FLEX PCA architectural audit and staged refactoring plan

Date: 2026-09-17

Status: proposed architecture; implementation not started.

Audit baseline: clean `flex-pca-memspeed` commit
`8d030b8e1ad92291cc5a100a2ce0fe8a2becb587`.

Validation level: static source inspection only. No compilation, test binary,
workflow, numerical-equivalence run, benchmark, or scientific experiment was run.

This is the single living design record for this refactor. Update it as phases land
rather than creating companion architecture documents.

## 1. Motivation

FLEX has grown into a substantial subsystem. Its commander and
shared/master/worker strategy remain useful boundaries, but below them workflow,
distributed transport, scientific state, reconstruction policy, configuration,
persistence, and numerical backends are interleaved.

As a result, changes to gridding, PCG, fitting, reconstruction, or geometry require
understanding several unrelated lifecycles and execution paths at once. This refactor
aims to make those responsibilities explicit, reduce coupling, and preserve current
scientific and distributed behavior. It is not a numerical rewrite, performance
project, branch-reconciliation plan, or change to scientific defaults.

## 2. Preservation contract

The first refactoring phases must preserve:

1. The public `flex_pca` command and current 3-D inputs and outputs.
2. A command-owned `parameters`/builder lifecycle, with explicit child-command
   lifecycles rather than parallel configuration authorities.
3. The same scientific phases in shared-memory and distributed runs.
4. Worker partials followed by master-owned reduction and publication.
5. Existing part names, formats, versions, field order, and atomic publication until
   an explicit schema migration.
6. The separate meanings of `rec_backend` and `rec_states_backend`.
7. Existing gridding/PCG algebra, normalization, support, stopping, filtering, paired
   fitting, cache/resume, sigma, CPU/GPU, and project-delivery behavior.
8. The `maxits_pcg=0` gridding/identity gate.

Do not mix ownership changes with file consolidation, renaming, optimization, or
scientific changes.

## 3. Current architecture and findings

```text
UI / executable
  -> commander: defaults, sigma bootstrap, lifecycle
  -> topology strategy: shared / master / worker
  -> simple_flex_pca_model::run_flex_pca
       -> validation, cache, mean, fit, embedding, states
       -> reconstruction, merge, project delivery, cleanup
       -> EM / gridding / PCG / CPU / GPU implementations

rounds, parameters, builder, project I/O, environment switches,
artifact naming, and mutable run state cross-cut the lower layers
```

The small commander lifecycle and strategy factory are sound and should remain
([`commander:40-74`](../../src/main/commanders/simple/simple_commanders_flex_pca.f90#L40-L74),
[`strategy:26-36,112-130`](../../src/main/strategies/parallelization/simple_flex_pca_strategy.f90#L26-L36)).

| Finding | Source evidence | Required boundary |
|---|---|---|
| **1. The application workflow is hidden in `simple_flex_pca_model`.** The 509-line driver owns validation, selection, cache/resume, fitting, embedding, state construction, reconstruction, project delivery, and teardown. | [`model:81-589`](../../src/main/flex/simple_flex_pca_model.f90#L81-L589) | Add a composed `flex_pca_application` behind the commander. It owns phase order and run invariants; “model” is reserved for scientific values. |
| **2. Distribution is abstracted but not encapsulated.** Domain modules query roles, launch rounds, size buffers, serialize, and reduce. `rounds` also mixes scheduling, stage IDs, formats, paths, and split policy. | [`rounds:17-95,147-225`](../../src/main/flex/simple_flex_pca_rounds.f90#L17-L95), [`EM iteration:489-529`](../../src/main/flex/simple_flex_pca_em_iter.f90#L489-L529), [`state reconstruction:93-175`](../../src/main/flex/simple_flex_pca_rec3D.f90#L93-L175) | Introduce a stage executor. Scientific phases compute typed local results; shared and qsys adapters decide dispatch, transport, reduction, and publication. |
| **3. Three-dimensional geometry is an end-to-end contract.** Validation, EM interfaces, projection code, artifacts, reconstruction, diagnostics, and project delivery expose concrete 3-D types. | [`model validation/delivery:658-699,1411-1448`](../../src/main/flex/simple_flex_pca_model.f90#L658-L699), [`EM API:388-473`](../../src/main/flex/simple_flex_pca_em.f90#L388-L473), [`latent ops:65-210,463-705`](../../src/main/flex/simple_flex_reconstructor_latent_ops.f90#L65-L210) | Separate prepared-observation formation from the common latent-statistics tail. Keep concrete 3-D geometry and delivery in explicit owners; do not introduce a speculative N-dimensional reconstructor. |
| **4. Backend policy is scattered.** The E-step mixes Cartesian/polar statistics, CPU/CUDA paths, preparation, resident banks, fallback gates, posterior inference, and M-step residual formation. The M-step gridding/PCG choice affects allocation, accumulation, payloads, solve, and paired merge. State reconstruction separately combines backend selection, distribution, CPU/GPU execution, filtering, filenames, and project writes. | [`E-step selection:266-388,590-729`](../../src/main/flex/simple_flex_pca_em_iter.f90#L266-L388), [`E-step formers:188-267`](../../src/main/flex/simple_flex_pca_em_estep.f90#L188-L267), [`M-step:98-127`](../../src/main/flex/simple_flex_pca_em_mstep.f90#L98-L127), [`rec3D:35-603`](../../src/main/flex/simple_flex_pca_rec3D.f90#L35-L603) | Use separate E-step, M-step, and state-reconstruction backend contracts. Keep posterior inference shared behind a common sufficient-statistics value, and put common FSC/filter/output/project work in one 3-D delivery service. CPU/CUDA is an implementation choice inside the E-step backend. |
| **5. `probe_fit_t` combines unrelated state and lifetimes.** It owns identity, model, convergence, mixture, polar, iteration scratch, distributed payloads, cross-FSC, merge, and PCG state. Submodules split files, not dependencies. | [`probe_fit_t:135-275`](../../src/main/flex/simple_flex_pca_em.f90#L135-L275), [`lifecycle:12-186`](../../src/main/flex/simple_flex_pca_em_state.f90#L12-L186) | Split immutable fit specification, persistent model/convergence state, per-iteration workspace, backend handle, transport payload, and diagnostics. Prefer composition. |
| **6. Configuration, project state, artifacts, and run lifetime have multiple owners.** UI and commander defaults drift; deep modules read environment policy, mutate parameters/project state, and retain module globals. | [`UI:504-559`](../../src/main/ui/simple/simple_ui_denoise.f90#L504-L559), [`commander:78-130,133-331`](../../src/main/commanders/simple/simple_commanders_flex_pca.f90#L78-L130), [`project/FSC I/O:628-719`](../../src/main/flex/simple_flex_pca_rec3D.f90#L628-L719) | Resolve one immutable `flex_run_config`; add a run-owned context, artifact catalog/codecs, and one application-level project gateway. Numerical code receives narrow prepared values rather than the whole builder. |

The requested 3-D reference code supports this ownership direction:

- [`simple_strategy3D_matcher`](../../src/main/strategies/search/simple_strategy3D_matcher.f90#L52-L215)
  exposes phase transitions and teardown.
- [`simple_matcher_3Drec`](../../src/main/strategies/search/simple_matcher_3Drec.f90#L23-L107)
  produces particle-side partial reconstructions.
- [`simple_rec3D_strategy`](../../src/main/strategies/parallelization/simple_rec3D_strategy.f90#L34-L103)
  owns lifecycle, topology, and backend selection.
- [`volassemble`](../../src/main/commanders/simple/simple_commanders_rec_distr.f90#L843-L1055)
  owns master-side assembly, diagnostics, output, and metadata.

Copy the producer/assembler split and ownership direction, not the size of those files
or their implementation details.

## 4. Target architecture

```text
commander
  -> flex_pca_application / run session
     |-- resolved config and run context
     |-- stage executor: shared | qsys master/worker
     |-- covariance-fit service
     |    |-- E-step
     |    |    |-- observation preparation
     |    |    |-- statistics backend: Cartesian | polar; CPU | CUDA
     |    |    `-- shared posterior inference
     |    `-- M-step system: gridding | coupled PCG
     |-- embedding and state-inference services
     |-- state reconstruction
     |    |-- backend: gridding | PCG
     |    `-- 3-D assembler and delivery
     |-- artifact catalog / codecs
     |-- project gateway
     `-- focused numerical kernels and domain values
```

| Component | Owns | Excludes |
|---|---|---|
| Commander | command normalization and application lifecycle | numerics and reduction |
| `flex_run_config` / context | resolved intent and run-scoped resources | phase sequencing |
| Application | phase order, use-case invariants, final success/failure | numerical loops |
| Stage executor | topology, partitioning, dispatch, completion, reduction | scientific interpretation |
| E-step backend | stage/iteration setup and prepared observations to sufficient statistics | posterior-model policy, project mutation, and M-step accumulation |
| Posterior inference | sufficient statistics to latent moments, contrast, and likelihood | projection, device, and transport selection |
| M-step backend | accumulation, codec, regularization, solve, paired transport | state placement and map delivery |
| State backend | weighted accumulation and backend-specific finalization | filenames and project mutation |
| 3-D delivery | FSC, common filtering/masking, products and project records | backend accumulation internals |
| Artifact/project services | identity, provenance, atomic I/O, final metadata | scheduling and numerical kernels |

Dependency rules:

1. Only the application sequences phases.
2. Only the executor knows process role or qsys topology.
3. Backend selection and fallback occur at stage or iteration boundaries, outside
   particle hot loops.
4. Numerical kernels do not receive `cmdline`, qsys objects, artifact names, project
   mutation access, environment policy, or the whole builder.
5. Codecs stay close to their typed payload; the executor decides when they are used.
6. Final artifacts and project mutation are master/application-owned.
7. Share only demonstrated PCG lattice mechanics; do not merge scalar reconstruct3D
   PCG and coupled FLEX PCG into one solver abstraction.

## 5. Staged migration

Each phase must preserve the current 3-D contract before the next begins.

| Phase | Change | Exit gate |
|---:|---|---|
| 0 | Freeze the baseline: inventory keys/defaults, project mutations, products, cache identities, part schemas, representative shared/distributed commands, and comparison tolerances. | The pinned behavior and artifacts can be reconstructed unambiguously. |
| 1 | Resolve `flex_run_config`; introduce run-owned resources plus concrete project/artifact ownership without changing formats. | Shared and worker roles receive the same resolved decisions; lifecycle ownership is explicit. |
| 2 | Extract named phases from `run_flex_pca` as wrappers around existing routines. | The top-level application is readable and preserves call order, outputs, and algorithms. |
| 3 | Move role checks, scheduling, waiting, serialization, and reduction behind the stage executor, one phase at a time. | Scientific phase services contain no role queries; shared/distributed results remain equivalent at the established tolerance. |
| 4 | Wrap gridding and PCG state reconstruction; extract one common 3-D assembler/delivery service. | Each backend reproduces its own baseline and output metadata; only the master publishes. |
| 5 | Extract M-step backend state and split persistent fit state from iteration workspace. | Single/paired paths preserve payloads, rank, convergence, and model output; iteration cleanup cannot free persistent state. |
| 6 | Extract the E-step backend lifecycle; move builder/project/cache/device preparation out of projection kernels and route Cartesian/polar CPU/CUDA paths through one sufficient-statistics value and shared posterior inference. | Backends accept prepared observations and narrow geometry arguments; each path preserves its baseline, fallback rules, and CPU/CUDA parity gates. |
| Post-core | Formalize artifact schemas and provenance. Add geometry, dimensions, backend, and version only with an explicit legacy-read or rejection policy. | Wrong or stale artifacts fail before folding; atomic publication remains intact. |

Any shared PCG lattice extraction should be a separate follow-up after these boundaries
stabilize. Share box/padding, KB window, limits, wrapping, and envelope construction
only; keep solver algebra and support policy separate.

## 6. Validation

Compilation and runtime gates are user-run unless separately authorized.

| Gate | Pass condition |
|---|---|
| Static boundaries | Numerical modules no longer import qsys/cmdline/project mutation/environment policy outside documented adapters; backend and topology selectors are localized. |
| Lifecycle | Extracted stateful types have clear construction/teardown; success, early return, and failure paths have one owner. |
| CPU/CUDA builds | Existing CPU and available CUDA configurations compile and link without new warnings or fallback changes. |
| Shared/distributed gridding | Fixed-input statistics, models, maps, filenames, and project fields match the baseline at the declared tolerance; only master publishes final state. |
| Shared/distributed PCG | Kernel/RHS shapes, support, ridge, start/stop policy, reduced products, and metadata match the PCG baseline; missing/mismatched parts fail closed. |
| Backend invariants | Each E-step path preserves its statistics and fallback behavior; gridding and PCG preserve their own behavior; cross-backend equality is not required; `maxits_pcg=0` is unchanged. |
| Paired and CPU/GPU paths | Half payloads, frame merge, delivered rank, CUDA parity, and CPU fallback remain within their existing contracts. |
| Artifacts and restart | Formats remain unchanged through the core refactor; after an explicit schema phase, supported legacy schemas are documented and cache/probe/basis/state handoffs reject provenance mismatches. |

No phase passes because line count fell. It passes only when the claimed dependency
movement and applicable behavior-preservation gates have evidence.

## 7. Non-goals and document maintenance

- No big-bang rewrite, numerical cleanup, precision change, loop reordering, or
  optimization mixed with extraction.
- No requirement that gridding and PCG produce the same result.
- No common superclass for scalar and coupled PCG solvers.
- No topology x geometry x solver x CPU/GPU subclass matrix.
- No generic N-dimensional reconstructor or other workflow design in this plan.
- No scheduler, cache, or filename redesign during behavior-preserving phases.
- No silent reinterpretation of old artifacts.
- No performance or scientific-quality claim without separate matched evidence.

As implementation lands, record completed phases and outstanding user-run gates here.
Update `doc/policies/flex_pca_policy.md` when operational ownership changes and
`src/main/flex/README.md` when live module boundaries change. Do not describe proposed
interfaces as implemented before they exist.
