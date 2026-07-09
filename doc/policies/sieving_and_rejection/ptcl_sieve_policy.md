# Particle Sieve Policy

This document defines the behavioral policy for staged particle sieving in
SIMPLE, implemented by
[../../src/main/sieve/simple_ptcl_sieve.f90](../../src/main/sieve/simple_ptcl_sieve.f90).

The policy captures lifecycle, tiering, chunk state transitions, rejection,
completion, and change guardrails for `ptcl_sieve`.

## 1. Scope

`ptcl_sieve` is responsible for staged 2D chunk orchestration:

1. coarse chunk generation from imported project records;
2. optional fine chunk generation from coarse outputs;
3. queue submission and completion polling;
4. class-average rejection and compatibility filtering;
5. per-tier completion accounting and final chunk combination;
6. exposing latest CAVG visualization metadata for stream UI.

It is not responsible for stream file watching/import logic (owned by stream
commanders) or low-level queue backend internals (owned by `qsys_env`).

## 2. Public Contract

Public API surface (type-bound methods on `ptcl_sieve`):

1. lifecycle: `new`, `kill`, `set_final_ingestion`
2. orchestration: `cycle`, `submit`, `collect_and_reject`
3. generation: `generate_chunks_coarse`, `generate_chunks_fine`
4. integration: `combine_completed_chunks`
5. status/query: `get_*` family (`get_finished`, counters, latest JPEG payload)
6. restart recovery: `import_existing_chunks_coarse`, `import_existing_chunks_fine`

Policy note: callers should use `cycle` and query methods as the normal
contract. Lower-level generation/submission calls are exposed for controlled
workflow composition and testing.

## 3. Tier Model

The sieve has two tiers:

1. coarse (pass 1): broad chunking and first rejection
2. fine (pass 2): refined chunking/rejection (optional)

Mode controls:

- `coarse_only=.true.`: execution and terminal accounting stop at coarse tier.
- `coarse_only=.false.`: fine tier is enabled when chunks are produced.
- `pre_chunked=.true.`: coarse chunks are imported from pre-existing project
  files rather than partitioned from a record list.

## 4. Chunk State Policy

Each chunk tracks:

- identity and paths: `id`, `folder`, `projfile`
- counts: `nptcls`, `nptcls_selected`
- lifecycle flags:
  - `abinitio2D_running`
  - `abinitio2D_complete`
  - `rejection_complete`
  - `complete`
  - `failed`

Sentinel files define state transitions:

- `ABINITIO2D_FINISHED` -> `abinitio2D_complete`
- `REJECTION_FINISHED`  -> `rejection_complete`
- `COMPLETE`            -> terminalized chunk
- `REJECTION_FAILED`    -> failed terminal chunk

Import policy from previous runs:

1. sentinel files are authoritative for recovered state;
2. incomplete non-failed chunks must regenerate command lines so they can be
   resubmitted;
3. missing chunk project files are warning-and-skip, not hard stop.

## 5. Cycle Policy

`cycle(project_list)` must execute in this order:

1. `collect_and_reject`
2. `generate_chunks_coarse(project_list)`
3. `generate_chunks_fine()` when not coarse-only
4. `submit`

This order is policy-significant and must not be rearranged without explicit
contract updates, because downstream tier eligibility and counters depend on it.

## 6. Generation Policy

### 6.1 Coarse generation

Coarse chunk creation is driven by particle thresholds and record inclusion:

1. consume non-included project records;
2. build chunk projects under `chunks_coarse`;
3. mark consumed records included;
4. emit/rewrite `imported_projects.txt` from currently included records.

In `pre_chunked` mode, coarse projects are copied from provided per-record
project files instead of repartitioning records.

### 6.2 Fine generation

Fine chunk generation is a merge/promote stage from eligible coarse outputs.
Only coarse chunks that passed rejection and are not already terminalized are
eligible inputs.

## 7. Submission and Scheduling Policy

Submission policy:

1. enforce `nparallel` running limit;
2. prioritize fine chunks over coarse chunks;
3. skip failed, running, or already completed chunks;
4. submit asynchronously via queue environment;
5. restore original working directory after submission pass.

Queue partition override policy:

- `SIMPLE_CHUNK_PARTITION` may override per-chunk partition metadata where
  implemented in generation/merge helpers.

## 8. Rejection Policy

`reject_cavgs` is tier-aware and applies two filters:

1. hard quality rejection (`evaluate_cavg_quality_hard_reject`);
2. compatibility model filtering (`class_compatibility`) for the respective
   tier model (`coarse_compatibility_model` or `fine_compatibility_model`).

Rejection outputs and artifacts:

1. project state is mapped through `map_cavgs_selection` and persisted;
2. selected and rejected class-average stacks/JPEGs are written;
3. `REJECTION_FINISHED` sentinel is emitted on completion;
4. chunk selected-count is updated from particle states.

Compatibility observability policy:

- rejection must log tier metrics (`a/b/c`, deltas, validity flags,
  convergence).
- convergence events should be logged explicitly when reached.

## 9. Completion and Finished Semantics

Chunk completion accounting:

- in coarse-only mode, coarse chunks can finalize accepted/rejected counts
  directly.
- in two-tier mode, terminal accepted/rejected counters are accumulated from
  fine chunks; coarse chunks act as feeders unless no fine tier exists.

`get_finished` contract:

1. requires at least one coarse chunk;
2. all coarse chunks must be complete or failed;
3. if `coarse_only`, this is terminal;
4. if not coarse-only and no fine chunks exist, coarse completion is terminal;
5. otherwise all fine chunks must be complete or failed.

## 10. Combination Policy

`combine_completed_chunks` merges eligible terminal chunk projects into one
project in `completedir`.

Eligibility:

- coarse-only: rejection-complete, non-failed coarse chunks
- two-tier: complete, non-failed fine chunks

No-op policy:

- do nothing when no eligible chunks exist;
- do nothing when target combined project already exists.

## 11. Counter and Query Policy

Counters must remain monotonic and query-safe:

- `get_n_accepted_ptcls`, `get_n_rejected_ptcls`,
  `get_n_accepted_micrographs` are cumulative terminal counters.
- `get_n_pass_1_non_rejected_ptcls` and `get_n_pass_2_non_rejected_ptcls`
  reflect non-terminal per-tier selected counts.
- `get_latest` must return `.false.` safely when latest payload is incomplete
  or uninitialized.

## 12. Failure Handling

Hard-fail conditions include structural inconsistencies (for example missing
required source files in pre-chunked input). Recovery-friendly paths should
prefer warning-and-skip for recoverable restart artifacts (for example missing
one imported chunk project among many).

## 13. Test Policy

Policy-level tests for `ptcl_sieve` must cover:

1. lifecycle defaults and idempotent reset;
2. import recovery from sentinel files;
3. tier counters and running-count semantics;
4. `get_finished` behavior across coarse-only and two-tier modes;
5. empty-cycle behavior on empty record lists;
6. latest-payload safe false-path.

Reference tester module:
[../../src/main/sieve/simple_ptcl_sieve_tester.f90](../../src/main/sieve/simple_ptcl_sieve_tester.f90).

## 14. Change Checklist

When changing `ptcl_sieve` behavior:

1. preserve cycle ordering unless policy is updated;
2. preserve fine-before-coarse submission priority;
3. preserve sentinel-driven state recovery contracts;
4. keep `get_finished` semantics backward compatible;
5. keep counter meaning stable (`terminal cumulative` vs `tier snapshot`);
6. update tester coverage for behavior changes;
7. update this policy document in the same change.
