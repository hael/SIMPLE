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
