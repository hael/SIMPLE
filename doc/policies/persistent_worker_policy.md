# Persistent Worker Policy

This document records durable workflow contracts for SIMPLE's persistent-worker
queue-system overlay. It is policy, not a line-by-line implementation map.

This document is the single policy authority for the persistent-worker backend.

## 1. Core Model

The persistent-worker backend is a dispatch overlay on the existing qsys
system.

Standard qsys backends still generate and run SIMPLE bash scripts. The
persistent-worker overlay changes how those scripts are dispatched:

1. `qsys_env` detects a worker qsys name such as `local_worker` or
   `slurm_worker`
2. the underlying base backend launches long-lived `simple_persistent_worker`
   processes
3. generated SIMPLE job scripts are queued to a TCP worker server instead of
   being submitted directly to the base backend
4. workers poll the server by heartbeat and execute script paths in local
   pthread slots
5. the existing SIMPLE scheduler still observes completion through filesystem
   sentinel files

The persistent-worker backend does not change scientific workflow semantics.
It only changes the transport used for generated distributed scripts.

## 2. Ownership

`simple_qsys_env.f90` owns worker-overlay activation, base-backend creation,
dispatch-backend creation, worker-server reuse checks, worker-process launch,
and qsys environment teardown.

`simple_qsys_factory.f90` owns construction of the internal
`persistent_worker` qsys backend.

`simple_qsys_persistent_worker.f90` owns the qsys backend type used for script
dispatch after the persistent pool is running. It is a qsys submission object,
not the worker-server owner.

`simple_qsys_ctrl.f90` owns script generation, scheduler loops, filesystem
completion polling, and dispatching generated script paths to the persistent
worker server.

`simple_persistent_worker_server.f90` owns the runtime singleton, TCP listener,
worker registry, task queues, worker heartbeats, task selection, and
termination replies.

`simple_persistent_worker_message_*.f90` owns the fixed wire-message types and
their serialisation contracts.

`production/simple_persistent_worker.f90` owns the compute-node worker process:
argument parsing, heartbeat loop, local thread-slot management, script
execution, and worker-process cleanup.

Top-level executables such as `simple_exec`, `simple_stream`, and `single_exec`
own final worker-server shutdown when they launched a persistent pool.

## 3. Backend Selection

A qsys name containing `_worker` activates the persistent-worker overlay. The
local qsys name is stripped before constructing the base backend. The dispatch
backend is then constructed with the internal factory name `persistent_worker`.

Examples:

- `local_worker` launches workers through `local`
- `slurm_worker` launches workers through `slurm`
- `lsf_worker`, `pbs_worker`, and `sge_worker` follow the same overlay pattern

`qsys_env` keeps the two roles separate:

- `base_qsys` and `base_qscripts` launch persistent worker processes through
  the underlying scheduler
- `dispatch_qsys` and `qscripts` queue generated SIMPLE scripts to the running
  persistent-worker server

Do not collapse these roles back into one ambiguous qsys object.

## 4. Runtime and Configuration

The module-level `persistent_worker` runtime lives in
`simple_persistent_worker_server.f90`. It records:

- the owning `persistent_worker_server` pointer
- the launch backend name
- the number of persistent worker processes launched
- the thread capacity per worker process

Worker concurrency is resolved in `qsys_env%new()`:

```text
nthr_per_worker = params%nthr
if worker_nthr > 0, use worker_nthr
if qsys_nthr is present, use qsys_nthr

n_workers = max(1, params%ncunits)
if workers > 0, use workers
```

`worker_nthr` controls advertised thread capacity per worker process.
`workers` controls how many persistent worker processes are launched. The
existing `qsys_ctrl` scheduler still controls in-flight generated scripts using
its normal `ncunits` bookkeeping and filesystem completion files.

An existing worker server may be reused only when it is running and compatible:

- same launch backend
- existing `nthr_per_worker` is at least the requested value
- existing worker count is at least the requested value

Otherwise the code throws rather than silently reusing the wrong pool.

## 5. Worker Server

`persistent_worker_server%new(n_workers, nthr_workers)` allocates listener
state, allocates shared worker data, initialises the listener mutex, starts the
TCP listener (skipped with `client_only`), records the bound port and host IP
list, and connects the dispatch client used by `queue_task`. `n_workers` is the
number of worker slots; `nthr_workers` is the thread capacity per worker.
Verbose logging is the compile-time `DEBUG` parameter; nothing propagates it.

`persistent_worker_data` is the state shared between the caller thread and the
listener pthread. It holds:

- the worker-slot count, active-worker count, and per-worker thread capacity
- queue-pressure and warm-up/scale-down state
- the termination flag and the count of TERMINATE replies sent since it was set

All mutable `persistent_worker_data` fields are protected by the listener
mutex. The worker heartbeat registry, the high, normal, and low priority task
queues, and `job_count` are local to the listener thread and need no mutex.
`queue_task` only sends the task to the listener over the dispatch client
connection; the listener increments its local `job_count` and assigns
`task%job_id` when it enqueues the task.

The listener thread must not perform blocking network replies while holding the
mutex. It copies the selected task to local storage, releases the mutex, and
then sends exactly one reply.

Task selection order is high, normal, low. A task is dispatchable only when the
worker heartbeat reports enough free thread capacity for `task%nthr`.

`qsys_ctrl` sets `task%priority` from `worker_priority` (`high` selects the
high queue, anything else the normal queue); several stream stages and their
jobs set `worker_priority=high`.
Nothing can enqueue into the low queue: the task message carries only a logical
priority flag.

## 6. Message Protocol

The wire protocol enumerates five fixed message types:

| Message | Value | Direction | Purpose |
| --- | --- | --- | --- |
| `WORKER_TERMINATE_MSG` | 1 | server to worker | shut down; the worker cancels running tasks. The same value is the kill sentinel `ipc_tcp_socket_server%kill` sends to the listener |
| `WORKER_HEARTBEAT_MSG` | 2 | worker to server | liveness and thread-load report |
| `WORKER_TASK_MSG` | 3 | none | set by `task%new()` only; never sent on the wire |
| `WORKER_NEW_TASK_MSG` | 4 | `queue_task` to listener, listener to worker | task submission, and the dispatched task (same record, type unchanged) |
| `WORKER_STATUS_MSG` | 5 | listener to worker or `queue_task` | idle reply to a heartbeat with no task; `queue_task` acknowledgement (idle = queued, error = queue full); error reply to unknown types |

The listener sends TERMINATE in reply to a heartbeat during server shutdown,
on cooldown scale-down, and when the `worker_id` is out of range or its UID
does not match the registered one.

Messages travel as their raw memory image. Every concrete message type
overrides `serialise` with its own copy of the transfer body; do not send an
extended type through the base `serialise`.

`TCP_BUFSZ = 1460` is the shared socket buffer size and must be large enough for
every message type.

## 7. Script Dispatch and Completion

Generated SIMPLE bash scripts remain the unit of distributed work.

When the active qsys object is `qsys_persistent_worker`,
`qsys_ctrl` calls `dispatch_task_to_persistent_worker`, creates a
`qsys_persistent_worker_message_task`, fills `nthr` and `script_path`, sets
`priority` from `worker_priority`, and queues the task.

If the worker server is missing or the queue is full, the current dispatch path
logs the failure. It does not currently propagate enqueue failure back to the
scheduler as a failed submission.

The existing scheduler remains the source of truth for completion:

- batch scheduling watches `JOB_FINISHED*`
- streaming scheduling watches `EXIT_CODE_JOB_*` and done/fail stacks

The worker server is a dispatch transport. It is not currently a completion
authority.

In every priority queue, the listener marks a dispatched task inactive
(`job_id = 0`) immediately, and queue compaction then removes it. This is the
current ownership-transfer model: after dispatch, filesystem completion files
decide whether the job finished.

## 8. Worker Process

`simple_persistent_worker` accepts only key-value command-line arguments:

- `nthr`
- `port`
- `server`
- `worker_id`

`port` and `server` are mandatory. Missing or invalid values terminate the
worker with `stop 1`. `nthr < 1` is forced to one thread slot.

Each worker computes a stable process identity:

```text
worker_uid = HOSTNAME_PID
```

The UID is included in every heartbeat. If the server sees a different UID for
an occupied `worker_id`, it sends `WORKER_TERMINATE_MSG` so a restarted process
does not inherit stale identity.

The heartbeat loop sends:

- `worker_id`
- `worker_uid`
- epoch heartbeat time
- total local thread slots
- currently used thread slots

The worker accepts at most one task per heartbeat reply. Each accepted task is
run as:

```text
bash <script_path>
```

The script does not need executable permissions because `bash` is explicit.

On `WORKER_TERMINATE_MSG`, or when the server cannot be reached, the worker
leaves the heartbeat loop and its cleanup cancels the pthreads of tasks that are
still running; it does not wait for them to finish. `WORKER_STATUS_MSG`
replies are ignored.

## 9. Thread Slots

Each worker slot has its own mutex and local state:

- task message
- pthread handle
- `started`
- `complete`
- `nthr`
- `exit_code`

`nthr == 0` is the authoritative idle signal. `get_nthr_used` sums the `nthr`
fields while locking one slot at a time.

A slot must be claimed before `pthread_create` by writing the task, setting
`nthr`, clearing `complete`, and setting `started` under the slot mutex. If a
previous pthread handle still exists, it is joined before the slot is reused.

On normal task exit, the worker thread records the script exit code, sets
`complete = .true.`, and sets `nthr = 0` under the slot mutex.

On worker-process cleanup, incomplete started threads are cancelled, started
threads are joined, slot state is reset, and each slot mutex is destroyed.

## 10. Current Limits

These are current implementation facts, not desired end-state contracts:

- workers do not send completion acknowledgements back to the server
- task timing fields in `qsys_persistent_worker_message_task` exist, but the
  production path does not currently maintain them end-to-end
- enqueue failure is logged but not returned to `qsys_ctrl` as a failed
  submission
- `qsys_env%kill()` does not kill the persistent-worker server; final shutdown
  is done by top-level executables
- top-level shutdown kills the server but does not deallocate and nullify the
  module-level server pointer
- `simple_persistent_worker` contains `is_safe_script_path`, but the task path
  currently bypasses it with `safe_path = .true.`
- `qsys_env` strips `_worker` in the local backend-selection variable; code that
  reads `qdescr('qsys_name')` still needs to be checked for `_worker` names

Do not write future policy as though these limitations have already been
removed.

## 11. Invariants

- Persistent workers execute generated scripts; they do not understand SIMPLE
  task graphs or scientific algorithms.
- The base qsys backend launches persistent worker processes.
- The `persistent_worker` qsys backend dispatches generated scripts to the TCP
  worker server.
- Worker-server task dispatch must preserve existing filesystem completion
  semantics until a real completion-acknowledgement protocol replaces them.
- The listener mutex is initialised by `persistent_worker_server%new` and
  destroyed by `persistent_worker_server%kill` after the listener thread has
  been joined.
- `job_count` increment and task `job_id` assignment stay in the listener
  thread, which owns the task queues; they need no mutex.
- Blocking network replies stay outside the server mutex.
- Concrete wire messages override `serialise` with their own copy of the
  transfer body; extended types are never sent through the base `serialise`.
- Worker slots use `nthr == 0` as the idle signal.
- Worker task execution remains explicit `bash <script_path>`.
- Production use must not claim script-path validation is active while
  `safe_path = .true.` remains in the worker.

## 12. Review Checklist

For persistent-worker changes, check:

- Does `_worker` still select a base backend and a separate dispatch backend?
- Are worker-server reuse checks still preventing incompatible pool reuse?
- Does `qsys_ctrl` still have a valid completion path if worker dispatch
  succeeds?
- Can enqueue failure be observed by the scheduler, or is it still only logged?
- Are task-queue ownership rules explicit for every priority path in use?
- Are server mutex boundaries free of blocking socket writes?
- Are worker thread slots claimed before pthread creation and joined before
  reuse or destruction?
- Are top-level executables still shutting down any launched worker server?
- Is the document still honest about script-path validation and completion
  acknowledgement state?
