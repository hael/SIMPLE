# Stream policies

Policies for the streaming pipeline: the master (p00) and the stages p01 to p07, their commanders
in `src/main/commanders/stream`, the stage types in `src/main/stream/stages`, the master's parts in
`src/main/stream/master`, the modules the stages share in `src/main/stream/shared`, the 2D pool
and chunk layer in `src/main/stream/pool2D`, and the watcher and utilities in `src/main/stream`.

Each policy states the contract the code keeps today and what a change must preserve. Its
"Known gaps" section names the open defects and the proposals of
`doc/refactoring_notes/completed/stream_area_review_2026-10-02.md`; a proposal becomes policy only when it is
decided and the policy is updated with the code.

| Policy | Read before changing |
|---|---|
| [reference_generation_policy.md](reference_generation_policy.md) | the initial analysis (p03), the picking references, their hand-off to reference picking (p04) |
| [stream_3D_ingestion_policy.md](stream_3D_ingestion_policy.md) | the pool's publications for 3D (p06), their import by multistate 3D (p07) |
| [pool2D_policy.md](pool2D_policy.md) | the 2D pool's schedule, sampling, resolution, dimensions, snapshots (p06 and `simple_stream_pool2D_utils`, `simple_stream_refine2D_utils`) |
| [ipc_policy.md](ipc_policy.md) | the stage-to-master pipes, the GUI metadata that crosses them, the master's heartbeat and the GUI's answers |
| [restart_policy.md](restart_policy.md) | what each stage does when it is started again, by the GUI or in an existing run folder |

Related policies elsewhere: the particle sieve (p03, p05) in
`doc/policies/sieving_and_rejection/ptcl_sieve_policy.md`; class-average quality in
`doc/policies/sieving_and_rejection/model_cavgs_rejection.md`.

Paths in these documents are relative to the repository root.

## Resources

The master gives every stage its threads and parts, and those of the jobs the stages run, from
one table (`src/main/stream/master/simple_stream_master_resources.f90`, decision 14 of the
5 October fix plan). Each value has a named default; a stage's environment variables override it,
read once by the master, which passes the values on the stage command lines and logs the table at
start. The stages hold no resource literals of their own. A stage run on its own, without the
master, takes its commander's defaults and reads none of these variables.

| Stage | Variables (`SIMPLE_STREAM_<SET>_…`) | What they set | Defaults (threads / parts) |
|---|---|---|---|
| p01 preprocessing | `PREPROC_NTHR`, `PREPROC_NPARTS` | the preprocessing jobs | 4 / 16 |
| p02 optics assignment | none | the stage | 1 |
| p03 initial analysis | `REFGEN_NTHR`, `REFGEN_NPARTS` | its 2D jobs and its sieve's chunks (`nthr2D`), its 3D job and reprojection (`nthr3D_pickrefs`, which the user's value overrides), the parts of its 2D and 3D jobs (`nparts`); its own process (picking) keeps 32 threads, its sieve runs 4 chunks (`nchunks`); each of its jobs claims as many threads on a persistent worker as its largest job uses (`worker_nthr`, the larger of `nthr2D` and `nthr3D_pickrefs`) | 16 / 1 |
| p04 reference picking | `PICK_NTHR`, `PICK_NPARTS` | the picking jobs | 8 / 8 |
| p05 particle sieving | `CHUNK_NTHR`, `CHUNK_NPARTS` | the threads of each chunk job and the chunks run at once (`nchunks`) | 16 / 4 |
| p06 pool 2D | `POOL_NTHR`, `POOL_NPARTS` | the pool's iterations | 8 / 6 |
| p07 multistate 3D | `SOLVE3D_NTHR`, `SOLVE3D_NPARTS` | the 3D jobs (`nthr3D`, `nparts3D`); its own process keeps 8 threads | 8 / 8 |

Each stage reads its own `SIMPLE_STREAM_<SET>_PARTITION` for its queue environment: `PREPROC`
(p01), `REFGEN` (p03 and its sieve's chunks), `PICK` (p04), `CHUNK` (p05's sieve), `POOL` (p06),
`SOLVE3D` (p07). A variable that is not a positive integer is ignored with a warning.
`solve2D_chunks` keeps its own `CHUNK` reads (decision 10). p03's extractions run with 4 threads
each (`EXTRACT_NTHR`), several at once.

With a persistent-worker queue (a `*_worker` queue name), the master starts the one worker server
(16 workers of 16 threads, `MASTER_NCUNITS` and `MASTER_QSYS_NTHR`) and gives its address to every
stage except p02 (`worker_server`). A stage starts no workers of its own: its queues' thread counts
(`worker_nthr`, or the queue's `qsys_nthr`) are only what each of its jobs claims on the master's
workers, and different queues of one stage may claim different counts (p03 and its sieve). The
master passes its threads per worker with the address (`worker_server_nthr`), and the stages pass
both on to the jobs that run queues of their own. A queue that claims more stops when it is made
(`qsys_env`): no worker could serve it, and a high-priority job waiting there would hold back
every normal-priority job on the server.
