# Particle Sieve Policy

This document defines the behavioral policy for staged particle sieving in
SIMPLE, implemented by
[simple_ptcl_sieve.f90](../../../src/main/sieve/simple_ptcl_sieve.f90).

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

1. lifecycle: `new`, `kill`, `set_final_ingestion`, `unset_final_ingestion`, `cancel`
2. orchestration: `cycle`, `submit`, `collect_and_reject`, `hand_off_final_set`
3. generation: `generate_chunks_coarse`, `generate_chunks_fine`
4. integration: `combine_completed_chunks`
5. status/query: `get_*` family (`get_finished`, counters, latest JPEG payload)
6. restart recovery: `import_existing_chunks_coarse`, `import_existing_chunks_fine`

Constructor policy (`new`):

- `new(params, settings, completedir, pre_chunked, optics_dir)`. `settings`
  (`ptcl_sieve_settings`, decision 12) holds everything the sieve reads; `params` gives the
  sieve's queue environment only. The streaming stages fill `settings` themselves (p05 with the
  mask diameter of `moldiam.txt`, p03 with the box default, coarse only, 4 chunks of 16 threads);
  `sieve_settings(params)` fills it from a parameters object (the standalone `sieve` commander,
  the tests).
- `optics_dir=<dir>` (optional) is where stream optics assignment publishes its optics maps;
  chunks handed to `completedir` then carry the newest map's groups (section 9).
- `settings%single_pass` enables coarse-only terminal semantics.
- `settings%lpstart` reaches the chunks' `solve2D` only when it was given (`sieve_lpstart`);
  otherwise their command lines carry none and `solve2D` derives the starting limit from the mask
  diameter (`min(max(mskdiam/12, 15), 20)` Å). Parameter validation lifts an unset `lpstart`
  to Nyquist (`fny`), so "given" means above Nyquist, not above 0.
- `settings%use_model` enables learned class-average rejection in both tiers.
- `settings%refs` pre-seeds coarse/fine compatibility models when the file exists; a missing
  file is warning-and-skip (non-fatal).
- tier overrides (0: the default): the population thresholds (`pop_coarse`, `pop_fine`, from
  `nptcls_coarse` and `nptcls_fine`), `lpstop_coarse`, `lpstop_fine`, `box_coarse`, `box_fine`,
  `nsample_coarse`, `nsample_fine`, `ncls_coarse` and `ncls_fine`.
- `settings%nthr` is the chunk jobs' thread count, also given to the sieve's queue environment
  (`qsys_nthr`), and `settings%partition` their queue partition (empty: the environment's).
- The chunks are classified and scored with `settings%mskdiam`, or the box's disc when it is 0
  (decision 7). Before 5 October 2026 scoring always used the box's disc; the sieve preset's
  weights were fitted with it, so p05's selections change and need a validation run before
  release.

Policy note: callers should use `cycle` and query methods as the normal
contract. Lower-level generation/submission calls are exposed for controlled
workflow composition and testing.

## 3. Tier Model

The sieve has two tiers:

1. coarse (pass 1): broad chunking and first rejection
2. fine (pass 2): refined chunking/rejection (optional)

Mode controls:

- `single_pass=yes`: execution and terminal accounting stop at coarse tier.
- `single_pass=no`: fine tier is enabled when chunks are produced.
- `pre_chunked=.true.`: coarse chunks are imported from pre-existing project
  files rather than partitioned from a record list.

## 4. Chunk State Policy

Each chunk tracks:

- identity and paths: `id`, `folder`, `projfile`
- counts: `nptcls`, `nptcls_selected`, `nattempts` (submissions of its 2D job)
- lifecycle flags:
  - `solve2D_running`
  - `solve2D_complete`
  - `rejection_complete`
  - `complete`
  - `failed`

Sentinel files define state transitions:

- `SOLVE2D_FINISHED` -> `solve2D_complete`
- `REJECTION_FINISHED`  -> `rejection_complete`
- `COMPLETE`            -> terminalized chunk
- `REJECTION_FAILED`    -> failed terminal chunk (written with `COMPLETE`, section 12)
- `FINAL_INGESTION`     -> a coarse chunk final ingestion staged (section 6.3):
  `solve2D_complete` and `rejection_complete`
- `EXIT_CODE_solve2D`   -> the exit status the chunk job's script writes; with no
  `SOLVE2D_FINISHED`, the job failed
- `EXIT_CODE_solve2D.job` -> the job's record of itself (pid, host, scheduler and id), written
  when it starts; with no exit status, the job may still run
- `chunk_input.simple`  -> the chunk's project as made, kept at its first submission for a retry
- `<folder>_attempt1`   -> the folder of a failed first attempt: the chunk has used its retry
- `chunk_mics.txt`      -> in a coarse chunk: the micrographs it took (`chunked_mics.txt` format)
- `consumed_coarse.txt` -> in a fine chunk: the ids of the coarse chunks it merged

A chunk exists once its project is written, which is done by temporary and rename. What it took
is written to its folder first (`chunk_mics.txt`, `consumed_coarse.txt`, and `FINAL_INGESTION` for
a staged chunk); what follows its project (the coarse chunks' `COMPLETE` after a fine chunk,
`chunked_mics.txt` after a coarse chunk) is redone by a restart. A retried or renewed chunk keeps
its lists. Chunk ids continue past the highest chunk folder on disk (`chunk_coarse_<id>`,
`chunk_fine_<id>`), so a folder whose making was cut short is never reused.

Import policy from previous runs:

1. sentinel files are authoritative for recovered state;
2. incomplete non-failed chunks must regenerate command lines so they can be
   resubmitted;
3. an incomplete chunk whose job recorded itself and wrote no exit status may
   still run (a crashed stage leaves it): the job is cancelled, the folder set
   aside to `<folder>_unfinished<k>`, and the chunk made afresh from
   `chunk_input.simple`, so a resubmitted job never shares its folder;
4. missing chunk project files are warning-and-skip, not hard stop; every chunk folder on disk
   is looked at, in id order, also past a missing one;
5. the coarse chunks a fine chunk that exists merged (its `consumed_coarse.txt`) are marked
   complete, also when the crash came before their `COMPLETE`, so they are not merged again;
6. the counters (section 11) are rebuilt from the imported chunks as the run counted them;
7. the particle-sieving stage rebuilds `chunked_mics.txt` from the `chunk_mics.txt` of the coarse
   chunks that exist before it reads it (`rebuild_chunked_mics`), so a crash between a chunk's
   project and `chunked_mics.txt` neither imports its micrographs twice nor loses them. While a
   chunk that exists has no list (made before lists were kept, or `pre_chunked`), the file is kept
   as it is.

`cancel` cancels the 2D jobs of the running chunks; the stage calls it when it
stops (p03 and p05). A cancelled chunk is resubmitted after a restart.

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
4. rewrite `chunked_mics.txt` (by temporary and rename): one
   `<project file> <micrograph index>` line per included record, in record
   order. p05 reads it on restart (`read_chunked_mics`) to import every set it
   chunked from again with those micrographs marked included. The included
   records stay a prefix of the record list, which slicing by id relies on.

In `pre_chunked` mode, coarse projects are copied from provided per-record
project files instead of repartitioning records.

### 6.2 Fine generation

Fine chunk generation is a merge/promote stage from eligible coarse outputs.
Only coarse chunks that passed rejection and are not already terminalized are
eligible inputs.

### 6.3 Final ingestion

Final ingestion flushes the particles left below the coarse threshold once no
more input is expected. It is by design that they skip the coarse pass.

1. Trigger: the driving stage calls `set_final_ingestion` once its upstream
   will hand on nothing more. Upstream says so with a marker in its folder
   (`STREAM_IDLE` or `STREAM_FINISHED`, `simple_defs_stream`; `upstream_done`):
   preprocessing (p01) writes `STREAM_IDLE` once no new movie has come for
   `MOVIES_IDLE_TIME_S` (15 minutes) and its sets are done (the movies short of a set of five go
   as a partial set after `PARTIAL_SET_QUIET_S`, 3 minutes, without a new movie), and reference
   picking (p04) once preprocessing is idle or stopped and its own sets are
   done; each writes `STREAM_FINISHED` when it stops. A marker counts once a
   watch made a settle time after it was first seen has found nothing new, so
   every project handed on before it has been taken.
   - The particle-sieving stage (p05) sets final ingestion on p04's marker. It
     withdraws it (`unset_final_ingestion`) when the marker goes or a new set
     arrives (retraction: preprocessing has new movies).
   - The initial analysis (p03) sets it once its "all" set is picked and every
     extraction of it is collected, at `NMICS_PLAN(2)` (500) picked
     micrographs, or with fewer once p01's marker counts and every imported
     micrograph is picked. It then cycles the sieve once more, so the leftover
     chunk exists before cycle 2 asks `get_finished`. It never withdraws it.
2. Staging: `generate_chunks_coarse` puts every record not yet included into
   one more coarse chunk. It flags that chunk `solve2D_complete` and
   `rejection_complete` without running 2D or rejection, and writes
   `FINAL_INGESTION` in its folder, so a restart restores the flags (section 4)
   instead of submitting the chunk for a coarse 2D run.
3. Two-tier mode: the staged chunk feeds the fine tier with all its particles,
   so they get the fine 2D and rejection only. Once final ingestion is set and
   every coarse chunk is complete or failed, the last fine chunk is flushed
   below the fine threshold. No chunk is marked final: fine chunks end in any
   order, so a flag on the last one made could reach the 2D pool before an
   earlier chunk still running, and the pool would publish its final result
   (and 3D start its final run) without that chunk's particles. The end of the
   intake is the empty final set of step 5.
4. Coarse-only mode (`single_pass=yes`, the initial analysis): the staged
   chunk is handed off as it is, with no screening. The initial analysis
   classifies and selects over the combined set in its cycle 2.
5. The final signal is the only one, and always sent: once final ingestion is
   set and every particle is through (no record left to chunk, every coarse
   and fine chunk complete or failed), `hand_off_final_set` hands on an empty
   project with `sieve_final=yes` (`sieve_final_c<ncoarse>_f<nfine>.simple`),
   once per chunk count. It is written in the cycle that hands off the last
   chunk (`collect_and_reject` runs before it), so it follows every chunk's
   hand-off and costs no extra cycle; the pool takes final sets after every
   other (`pool2D_policy.md`, section 3). It also covers nothing pending at
   final ingestion, a last chunk that failed, and a coarse chunk finished with
   nothing selected after the last fine chunk.
6. Retraction reaches the pool: the pool takes back its final note when a set
   with particles and without the flag arrives after it (`pool2D_policy.md`).

## 7. Submission and Scheduling Policy

Submission policy:

1. enforce `nparallel` running limit;
2. prioritize fine chunks over coarse chunks;
3. skip failed, running, or already completed chunks;
4. submit asynchronously via queue environment;
5. restore original working directory after submission pass.

Queue partition override policy:

- `SIMPLE_STREAM_CHUNK_PARTITION` (`simple_defs_environment`) overrides the queue
  partition of the coarse chunks and the merged fine chunks; it is the variable the
  particle-sieving stage's queue reads. The former name `SIMPLE_CHUNK_PARTITION` is no
  longer read.

Worker server: the sieve's queue reuses the persistent-worker server of the process that drives
it. Reusing a server never turns its warm-up cooldown (autoscaling) off; a streaming queue turns
it on (`qsys_env`).

## 8. Rejection Policy

`reject_cavgs` is tier-aware and applies two filters:

1. hard quality rejection (`evaluate_cavg_quality_hard_reject`);
2. compatibility model filtering (`class_compatibility`) for the respective
   tier model (`coarse_compatibility_model` or `fine_compatibility_model`).

Quality model policy (both tiers):

- both tiers load the same preset, `sieve` (`CAVG_QUALITY_MODEL_SIEVE_DEFAULT`), and apply the
  `sieve` hard-gate context. The fine tier does not switch to the `chunk` preset
  (`chunk100mics`) or the `chunk` context, although its merged chunks are as large as the ones
  `model_cavgs_rejection.md` describes for that context. The fine tier is the sieve's second
  pass, not the stream's chunk 2D: its class averages come from the sieve's own 2D runs, so the
  sieve scores them as it scores the coarse tier. This is a decision (stream fix plan,
  decision 8, confirmed 5 October 2026), not a gap.
- as a result, the fine tier applies none of the shared `chunk`/`pool` gates (the `res > 40 Å`
  gate included), nor the `chunk` local-variance and band-pass floors. Flush chunks of a few
  hundred particles are scored the same way.
- when `use_model=yes`, rejection runs model-based quality scoring
  (`evaluate_cavg_quality`) after hard rejection.
- when `use_model=no`, rejection uses hard rejection only.
- giving a tier another preset or context is a policy change. It needs a validation run, and the
  preset's spec must record which tier's tables it was fitted on.

Rejection outputs and artifacts:

1. project state is mapped through `map_cavgs_selection` and persisted;
2. selected and rejected class-average stacks/JPEGs are written;
3. a rejected class records `rejection_reason` as `<tier>_reject:<reason>` with underscores for
   blanks (`coarse_reject:low_population`), since the orientation reader splits values at blanks;
4. an all-class JPEG (`*_all_reasons.jpg`) is written with reason-coded
  borders plus a sidecar key file (`*_all_reasons.jpg.key.txt`);
5. `REJECTION_FINISHED` sentinel is emitted on completion;
6. chunk selected-count is updated from particle states.

Cleanup retention policy (`cleanup_chunk`):

1. cleanup runs after rejection completes;
2. keep lifecycle sentinels used by restart/import recovery:
  `SOLVE2D_FINISHED`, `REJECTION_FINISHED`, `COMPLETE`,
  `REJECTION_FAILED`, `FINAL_INGESTION`;
3. keep chunk project metadata file and `frcs.bin`;
4. keep selected/rejected JPEG renderings;
5. keep all-reasons reason-overlay JPEG and its sidecar key file;
6. keep latest iteration JPEG for the chunk;
7. keep final iteration stacks for all three stack variants when present:
  whole stack (non-`_even`/`_odd`), `_even`, and `_odd`;
8. keep the highest-rank sigma STAR candidate (`sigma*.star`, preferring
  `_iterNNN` when available).

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

Hand-off policy:

- a terminal chunk (a rejection-complete coarse chunk in coarse-only mode, otherwise a
  rejection-complete fine chunk) is handed to the next stage as a copy of its project in
  `completedir`, before its `COMPLETE` sentinel is written;
- with an `optics_dir`, the copy gets the groups of the newest optics map in that directory,
  applied by import index to its micrographs, stacks and particles, with the map's optics
  segment; without a map yet, or without an `optics_dir`, it is an exact copy;
- the chunk's own project is never given optics groups: fine chunks merge coarse chunks, and
  the merge offsets each source's group ids, so the groups are applied only on the hand-off copy;
- the copy is written as `<stem>.tmp` in `completedir` and renamed to `<stem>.simple`, so the
  stage watching `completedir` for `*.simple` never reads a partial project and needs no long
  settle time. Keep the temporary name off the `.simple` suffix and in `completedir` (a rename
  is atomic only within one file system).

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
  `get_n_accepted_micrographs` are cumulative terminal counters: the finalised chunks of the
  terminal tier, and in two-tier mode the coarse chunks finalised with nothing selected (all
  their particles rejected), which no fine chunk counts.
- every counter is cumulative across restarts: `new` rebuilds them from the imported chunks.
- `get_n_pass_1_non_rejected_ptcls` and `get_n_pass_2_non_rejected_ptcls`
  reflect non-terminal per-tier selected counts.
- `get_n_coarse_accepted_ptcls` and `get_n_coarse_rejected_ptcls` are
  cumulative coarse-tier rejection results.
- `get_n_fine_accepted_ptcls` is the cumulative fine-tier accepted count.
- `get_latest` must return `.false.` safely when latest payload is incomplete
  or uninitialized.

## 12. Failure Handling

Hard-fail conditions include structural inconsistencies (for example missing
required source files in pre-chunked input). Recovery-friendly paths should
prefer warning-and-skip for recoverable restart artifacts (for example missing
one imported chunk project among many).

A chunk fails, and the stage does not stop, when its 2D job exits without
`SOLVE2D_FINISHED` (its script wrote `EXIT_CODE_solve2D`), or when its 2D run left
no class averages rejection can score (no output stack, fewer images than
classes, or class rows that do not match) (`fail_chunk`):

1. **The first failure retries it.** The chunk's folder is moved to
   `<folder>_attempt1` and kept, and the chunk starts again in a fresh folder from
   `chunk_input.simple`, with its command line made again.
2. **The second failure drops it.** It is written `REJECTION_FAILED` and
   `COMPLETE`; none of its particles is handed on; it counts as complete for the
   final flush and for `get_finished`. The sieve logs a warning, counts the chunk
   and its particles (`get_n_failed_chunks`, `get_n_failed_ptcls`), and p05 shows
   both in its GUI status.
3. A chunk without `chunk_input.simple` (made before this rule) is dropped at its
   first failure.
4. **A job that vanishes** before its script writes an exit status (walltime, a
   lost node) is found by its job's liveness check (`qsys_async_job`; the chunk's
   job runs through it): gone at two checks `JOB_LIVENESS_S` (5 minutes) apart, it
   gets a lost status and the chunk takes the failure path above (follow-up plan,
   decision 29).

## 13. Test Policy

Policy-level tests for `ptcl_sieve` must cover:

1. lifecycle defaults and idempotent reset;
2. import recovery from sentinel files, including a chunk final ingestion staged;
3. tier counters and running-count semantics;
4. `get_finished` behavior across coarse-only and two-tier modes;
5. empty-cycle behavior on empty record lists;
6. latest-payload safe false-path;
7. the hand-off copy in `completedir`, with the newest optics map's groups when an
   `optics_dir` is given, and no `.tmp` left behind;
8. a failed chunk retried once in a fresh folder, and dropped with its particles
   counted when it fails again.

Reference tester module:
[simple_ptcl_sieve_tester.f90](../../../src/main/sieve/simple_ptcl_sieve_tester.f90).

## 14. Change Checklist

When changing `ptcl_sieve` behavior:

1. preserve cycle ordering unless policy is updated;
2. preserve fine-before-coarse submission priority;
3. preserve sentinel-driven state recovery contracts;
4. keep `get_finished` semantics backward compatible;
5. keep counter meaning stable (`terminal cumulative` vs `tier snapshot`);
6. preserve cleanup artifact retention semantics (including sentinel files and
  whole/even/odd final iteration stack retention) unless policy is explicitly
  revised;
7. update tester coverage for behavior changes;
8. update this policy document in the same change.
