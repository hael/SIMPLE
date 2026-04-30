# simple_persistent_worker_server Policy

## Status

This document is a proposal / design record for
`src/utils/qsys/worker/simple_persistent_worker_server.f90`.
It should be promoted to `doc/policies/` once the worker backend has been
exercised in production and the interfaces are considered stable.

See also:
- [simple_qsys_worker_policy.md](simple_qsys_worker_policy.md) — qsys overlay and
  environment-level policy
- [simple_persistent_worker_policy.md](simple_persistent_worker_policy.md) — worker process policy

## Scope

This document defines the architectural policy for the
`simple_persistent_worker_server` module.

It covers:

- the role and responsibilities of `persistent_worker_server`
- the two-type design (`persistent_worker_server` + `persistent_worker_data`)
- the mutex ownership model
- `new()` and `kill()` contracts
- `queue_task()` contract and job-id assignment
- the listener pthread contract (`worker_listener_thread`)
- task selection and the critical section boundary
- `pack_task_queue` compaction
- `set_debug` and the debug logging flag
- the worker-uid identity check
- named constants and their design rationale
- invariants that all future changes must preserve

---

## Core design rule

`persistent_worker_server` is a single-server, multi-worker coordination hub.
Its only job is to:

1. Accept TCP connections from `simple_persistent_worker` processes (heartbeats).
2. Reply with a task, a TERMINATE command, or an idle status.
3. Expose a `queue_task()` call so the SIMPLE control layer can enqueue work.

The server has no knowledge of what the scripts do.  It routes tasks; it does
not execute them.

---

## Two-type design

The module exposes two types with distinct threading roles:

| Type | Accessed by | Role |
|---|---|---|
| `persistent_worker_server` | caller thread only | lifecycle, queueing, configuration |
| `persistent_worker_data` | **both** caller and listener pthread | shared state under mutex |

`persistent_worker_data` is heap-allocated and accessed exclusively through a
`pointer` field on `persistent_worker_server`.  Its address is passed to the
listener thread as a C `void*` via `listener_args%data_ptr`.

The split must be preserved.  No field that is written by the listener thread
may live directly on `persistent_worker_server`.

---

## Mutex ownership

The mutex (`listener_args%mutex`) is owned by `persistent_worker_server`.

Ownership rules:

1. Initialised in `new()` with `c_pthread_mutex_init`.
2. Destroyed in `kill()` with `c_pthread_mutex_destroy`, and only after
   `ipc_socket%kill()` has returned (guaranteeing the listener pthread has
   been joined).
3. `ipc_tcp_socket` borrows the mutex by pointer; it must never init or
   destroy it.
4. No other procedure may init or destroy the mutex.

The mutex protects all mutable fields of `persistent_worker_data`:

- `n_workers`
- `l_terminate`
- `l_debug`
- `workers(:)` (heartbeat registry)
- `tasks_priority_high(:)`, `tasks_priority_norm(:)`, `tasks_priority_low(:)`

`job_count` on `persistent_worker_server` is also updated under the mutex (inside
`queue_task`) because it feeds `task%job_id` which is written into the shared
queue.

---

## `new()` contract

`new(nthr_workers)` must:

1. Check `associated(worker_data)` and `associated(listener_args)`.
   If either is already associated, log a no-op message and return without
   reinitialising.
2. Allocate `listener_args` and `worker_data`.
3. Call `c_pthread_mutex_init` on `listener_args%mutex`.
   If it fails, deallocate both, nullify both pointers, and return.
4. Copy `self%l_debug` → `worker_data%l_debug` **before** spawning the thread.
5. Set `listener_args%data_ptr = c_loc(worker_data)`.
6. Call `ipc_socket%init_server(c_funloc(worker_listener_thread), c_loc(listener_args))`.
7. Query and store `port` and `host_ips` from the socket.

The mutex init and the `l_debug` copy must happen before step 6.  Rationale:
the listener pthread may run immediately after `init_server` returns.

---

## `kill()` contract

`kill()` must proceed in this exact order:

1. Reset `nthr_workers = 0`, `job_count = 0`, call `host_ips%kill()`.
2. If `listener_args` or `worker_data` is not associated, the listener was
   never started.  Call `ipc_socket%kill()`, set `port = 0`, return.
3. Lock the mutex; set `worker_data%l_terminate = .true.`; unlock.
4. Call `c_usleep(KILL_WAIT_TIME_US)` (2 s) to allow workers to receive
   TERMINATE before the socket is closed.
5. Call `ipc_socket%kill()` — this joins the listener pthread.
6. Set `self%port = 0` — only after the thread has been joined so that
   `is_running()` remains `.true.` until shutdown is complete.
7. Call `c_pthread_mutex_destroy(listener_args%mutex)`.
8. Deallocate and nullify `worker_data`.
9. Deallocate and nullify `listener_args`.

Steps 5–9 must not be reordered.  In particular:
- `port = 0` must not be written before step 5.
- `pthread_mutex_destroy` must not be called before step 5.

---

## `queue_task()` contract

`queue_task(task, priority)` enqueues one task and returns `.true.` on
success.

Rules:

1. Guard: if `worker_data` is not associated, log and return `.false.`.
2. Lock the mutex.
3. Increment `job_count`; assign `task%job_id = job_count`.
   Both operations must be inside the same critical section.
4. Find the first slot with `job_id == 0` in the selected queue:
   - `priority == 'high'` → `tasks_priority_high`
   - `priority == 'low'`  → `tasks_priority_low`
   - anything else        → `tasks_priority_norm` (including `'norm'`)
5. Write the task into that slot; set `queued = .true.`.
6. If no free slot is found, log a queue-full warning; set `queued = .false.`.
7. Unlock the mutex.
8. Return `queued`.

The caller must not set `task%job_id` before calling `queue_task`; it is
assigned inside the critical section to guarantee uniqueness.

A `queued = .false.` return is not fatal.  The caller (`simple_qsys_ctrl`)
may retry or log and continue.

---

## Listener pthread contract (`worker_listener_thread`)

`worker_listener_thread` is a `bind(c)` subroutine spawned by
`ipc_tcp_socket`.  It runs in a dedicated OS thread.

### Startup

On entry, the thread must:

1. Recover Fortran pointers from the `c_ptr` argument.
2. If `status` cannot be associated (null `data_ptr`), signal `args%ready = 1`
   under the mutex and return immediately.  Do not proceed to the accept loop.
3. Signal readiness: lock → `args%ready = 1` → unlock.

### Accept loop

One iteration per incoming TCP connection.  The loop exits when:

- `recv_msg` returns `ok = .false.` (socket closed or error), or
- a `WORKER_TERMINATE_MSG` is received (the server itself is shutting down).

### Message dispatch

| Message type | Action |
|---|---|
| `WORKER_TERMINATE_MSG` | Log, echo TERMINATE back to caller, exit loop |
| `WORKER_HEARTBEAT_MSG` | Update worker registry, select task, reply (see below) |
| invalid length (≤ 0 or > `TCP_BUFSZ`) | Reply ERROR status; continue |
| unknown type | Log unknown type; fall through to idle reply |

### Critical section boundary

For heartbeat messages, the critical section (mutex held) must:

- Read `status%l_terminate`
- Update `status%workers(worker_id)` (bounds-checked; see worker-uid rules)
- Scan queues high → norm → low for the first dispatchable task
- Set `task%submitted = .true.` on the selected task
- Call `pack_task_queue` on all three queues (skip if `send_terminate`)
- Copy the selected task to a **local** `dispatch_task`

The mutex must be **released before** any call to `repl_msg`.  Rationale:
`repl_msg` is a blocking network write; holding the mutex during it would
block `queue_task` for an unbounded time.

### Task selection

Selection scans queues in priority order: high, norm, low.
A task is eligible if both:

- `task%job_id > 0` — slot is occupied
- `task%submitted == .false.` — not yet dispatched

A task is dispatchable if additionally:

- `heartbeat_msg%nthr_total - heartbeat_msg%nthr_used >= task%nthr`

A task that is eligible but not dispatchable (insufficient free threads) is
skipped with a debug log line.  The scan continues to the next task in the
same queue, then moves to lower-priority queues.  Only one task is dispatched
per heartbeat.

### Replies

After the critical section, the reply is chosen in order:

1. `send_terminate == .true.` → send `WORKER_TERMINATE_MSG`
2. `has_task == .true.` → send `WORKER_TASK_MSG`
3. Neither → send `WORKER_STATUS_MSG` with `status = WORKER_STATUS_IDLE`

Exactly one reply is sent per accepted connection.

---

## Worker-uid identity check

When a heartbeat arrives for a `worker_id` that already has a non-empty
`worker_uid`, the server checks whether the new heartbeat carries the same
`worker_uid`:

- **Same**: update the registry normally.
- **Different**: log a warning, set `send_terminate = .true.`.

Rationale: a changed `worker_uid` at the same slot means the original worker
process crashed or was restarted.  Sending TERMINATE forces the new process to
re-establish its identity cleanly rather than inheriting stale state.

---

## `pack_task_queue` compaction

`pack_task_queue(queue)` compacts one priority queue in-place.

A task slot is **active** if `job_id > 0 .and. end_time == 0`.  The `end_time`
field is reserved for completion confirmation; until it is set, the task is
considered pending or submitted-but-unconfirmed.  Retaining unconfirmed tasks
allows the server to detect worker crashes (the task will remain visible in
the queue indefinitely until the caller explicitly clears it or the server is
restarted).

`pack_task_queue` must:

1. Shift all active entries to the front of the array, preserving order.
2. Zero-fill (reset to `qsys_persistent_worker_message_task()`) all trailing slots.

`pack_task_queue` must only be called inside the critical section and only
when `send_terminate == .false.`.

---

## Debug logging flag

`l_debug` exists on both `persistent_worker_server` (caller-thread view) and
`persistent_worker_data` (listener-thread view).  They must be kept in sync.

Synchronisation rules:

- `new()` copies `self%l_debug → worker_data%l_debug` before spawning the
  thread.
- `set_debug(l_debug)` updates both: sets `self%l_debug`, then acquires the
  mutex and updates `worker_data%l_debug`.
- `set_debug` is safe to call before or after `new()`.

**Always-on log messages** (regardless of `l_debug`):

- errors and mutex failures
- TERMINATE sent/received
- task dispatch (job id, worker id)
- queue-full rejections
- out-of-range `worker_id`
- changed `worker_uid` warning

**Debug-only log messages** (`l_debug == .true.` required):

- per-heartbeat log lines (worker id, time, thread load)
- per-task "skipping — insufficient threads" lines
- per-heartbeat "no tasks available" idle replies

Do not add unconditional high-frequency log lines to the accept loop.
Unbuffered Fortran I/O in the hot path stalls the accept loop.

---

## Named constants

| Constant | Value | Rationale |
|---|---|---|
| `WORKER_STATUS_IDLE` | 0 | status reply code: server has no pending task |
| `WORKER_STATUS_ERROR` | 2 | status reply code: malformed or invalid message |
| `MAX_WORKERS` | 100 | maximum `worker_id` accepted; array sizing |
| `TASK_QUEUE_SIZE` | 1000 | slots per priority level; array sizing |
| `TCP_BUFSZ` | 1460 | Ethernet MTU (1500) − IP header (20) − TCP header (20) |
| `KILL_WAIT_TIME_US` | 2 000 000 | 2 s grace period before socket close in `kill()` |

`TCP_BUFSZ` is exported (`public`) so `simple_persistent_worker` can declare its receive
buffer to the same size.  All message structs must fit within `TCP_BUFSZ`
bytes.  If a new message type would exceed it, update `TCP_BUFSZ` and this
policy together.

---

## Invariants

The following invariants must be preserved by all future changes:

1. `persistent_worker_data` is the only type whose fields are written by the
   listener pthread.  No field written by the listener thread may live
   directly on `persistent_worker_server`.
2. The mutex is init in `new()`, destroyed in `kill()` after thread join.
   No other procedure may init or destroy it.
3. `job_count` increment and `task%job_id` assignment must be inside the
   same critical section in `queue_task`.
4. `self%port = 0` must be written only after `ipc_socket%kill()` returns
   (i.e., after the listener pthread has been joined).
5. `repl_msg` (blocking network write) must never be called while the
   mutex is held.
6. `pack_task_queue` must only be called inside the critical section and
   only when `send_terminate == .false.`.
7. `worker_data%l_debug` must be set before `ipc_socket%init_server` in
   `new()`, and updated under the mutex in `set_debug`.
8. A changed `worker_uid` at an occupied slot always triggers TERMINATE,
   not a registry update.
9. All message structs must fit within `TCP_BUFSZ` bytes.
10. `pack_task_queue` retains submitted-but-unconfirmed tasks
    (`end_time == 0`) so that unacknowledged dispatches survive compaction.
