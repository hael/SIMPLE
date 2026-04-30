# simple_persistent_worker Policy

## Status

This document is a proposal / design record for the `simple_persistent_worker` process.
It should be promoted to `doc/policies/` once the worker backend has been
exercised in production and the interfaces are considered stable.

See also [simple_qsys_worker_policy.md](simple_qsys_worker_policy.md) for the
server-side and qsys-overlay policy.

## Scope

This document defines the architectural policy for `production/simple_persistent_worker.f90`.

It covers:

- the role and responsibilities of `simple_persistent_worker`
- the command-line interface contract
- the heartbeat loop and server communication
- thread-slot management and the `worker_thread_args` type
- the task execution contract and script-path validation
- the `worker_uid` identity mechanism
- logging conventions
- cleanup and teardown ordering
- invariants that all future changes must preserve

---

## Core design rule

`simple_persistent_worker` is a stateless compute-node agent.  Its only job is to:

1. Advertise available thread capacity to the server via periodic heartbeats.
2. Accept tasks in reply and execute them in isolated pthreads.
3. Exit cleanly when the server sends `WORKER_TERMINATE_MSG`.

`simple_persistent_worker` has no knowledge of the job queue, the task graph, or the
scientific algorithms.  It runs scripts exactly as handed to it.  All
orchestration decisions live in `simple_persistent_worker_server`.

---

## Command-line interface

All arguments are `key=value` pairs read via `command_argument_count` /
`get_command_argument`.  There are no positional arguments.

| Key | Type | Default | Mandatory | Meaning |
|---|---|---|---|---|
| `nthr` | integer | 1 | no | Number of parallel task slots |
| `port` | integer | — | **yes** | TCP port of the central server |
| `server` | string | — | **yes** | Comma-separated list of server IP addresses |
| `worker_id` | integer | 0 | no | Unique ID assigned by `qsys_env` at launch |

Validation rules:

- `nthr < 1` is forced to `1` with a warning; it is not a fatal error.
- Missing or invalid `port` (≤ 0) is fatal: the process calls `stop 1`.
- Missing or empty `server` list is fatal: the process calls `stop 1`.
- `worker_id = 0` is a valid sentinel (server will assign a slot by heartbeat
  if the worker has not been pre-registered).

The argument parser must not reject unknown keys silently.  Future keys can be
added without changing the parser for existing keys.

---

## Worker UID

Each `simple_persistent_worker` process builds a unique identity string at startup:

```
<HOSTNAME>_<PID>
```

This string (`worker_uid`) is:

- stable for the entire lifetime of the process
- included in every `WORKER_HEARTBEAT_MSG` sent to the server
- printed to stdout at startup for log correlation

The server uses `worker_uid` to detect cases where a slot is occupied by a
new process (e.g. after a crash and re-launch) and logs a warning when the
UID changes for a given `worker_id`.

Construction rules:

- `HOSTNAME` is read from the `HOSTNAME` environment variable.
  If that variable is absent or empty, use the string `'unknown'`.
- `PID` is obtained via `c_getpid()` (C interop).
- The two components are joined with `'_'`.
- The combined string must fit in 256 characters.  Truncation is acceptable
  if the hostname is unusually long; the uniqueness guarantee is then
  best-effort.

---

## Heartbeat loop

The heartbeat loop is the main execution loop of `simple_persistent_worker`.  It must:

1. Populate a fresh `qsys_persistent_worker_message_heartbeat` on each iteration:
   - `worker_id`       — from command-line argument
   - `worker_uid`      — computed at startup (constant for process lifetime)
   - `heartbeat_time`  — current epoch time via `c_time(0)`
   - `nthr_total`      — `nthr` from command-line argument (constant)
   - `nthr_used`       — result of `get_nthr_used()` at the time of send
2. Serialise and send the heartbeat; block for a reply with a timeout of
   `HEARTBEAT_TIMEOUT_MS = 5000 ms` and up to `HEARTBEAT_MAX_RETRY = 5`
   retries.
3. On send failure: log and exit the loop (treat as server lost).
4. Decode the reply by reading the first `integer` of the buffer as the
   message type, then dispatch:

| Reply type | Action |
|---|---|
| `WORKER_TERMINATE_MSG` | Set `l_terminate = .true.`; exit loop |
| `WORKER_TASK_MSG` | Validate script path; call `start_worker_thread` |
| `WORKER_STATUS_MSG` | Informational; no action required |
| any other value | Log unknown type; continue |

5. Sleep `POLL_TIME_US = 200 000 µs` (0.2 s) between iterations to avoid
   busy-polling when the server replies with STATUS (no tasks available).

The heartbeat loop must not call `stop` or `error stop`.  All exit paths go
through the loop condition or the TERMINATE case, which allows the cleanup
block to run.

---

## Thread-slot management

### `worker_thread_args`

Each task slot is represented by one `worker_thread_args` element:

| Field | Type | Meaning |
|---|---|---|
| `task_msg` | `qsys_persistent_worker_message_task` | task currently in this slot |
| `mutex` | `c_pthread_mutex_t` | guards all mutable fields in this struct |
| `thread` | `c_pthread_t` | pthread handle |
| `started` | logical | true after successful `pthread_create` |
| `complete` | logical | set by the thread body on exit |
| `nthr` | integer | > 0 while slot is occupied; 0 when idle |
| `exit_code` | integer | script exit status after completion |

Each slot has its own mutex.  The main thread and the task thread communicate
only via mutex-protected reads and writes to these fields.  There is no
global lock on the slot array.

### Slot lifecycle

```
idle:      nthr == 0, started == .false.
running:   nthr > 0,  started == .true.,  complete == .false.
done:      nthr == 0, started == .true.,  complete == .true.   ← must be joined before reuse
```

`start_worker_thread` must check for the `done` state and join the old
pthread handle before reusing the slot.  The join must happen outside the
mutex lock to avoid a deadlock with the thread body, which holds the same
lock when writing its final state.

### Slot claim ordering

The slot must be claimed (all fields set, `started = .true.`) while the mutex
is held, *before* `pthread_create` is called.  Rationale: if the mutex were
released before `pthread_create`, a concurrent call to `start_worker_thread`
could observe `nthr == 0` and steal the same slot.

If `pthread_create` fails, the slot must be reset to idle under the mutex:
`started = .false.`, `complete = .true.`, `nthr = 0`.

### `get_nthr_used`

`get_nthr_used` sums `threads(i)%nthr` across all slots, locking each slot
individually.  It must never hold more than one slot mutex simultaneously.

---

## Task execution contract

### Script-path validation (`is_safe_script_path`)

Before executing a task the worker must validate `task_msg%script_path`:

1. The trimmed path must be non-empty.
2. The path must be absolute (first character is `/`).
3. Every character must be in the allowed set:
   `a–z A–Z 0–9 / . _ -`
4. The file must exist (`inquire(file=…, exist=…)`).

A path that fails any check must be rejected with a warning log line and the
task dropped.  The server will re-queue on the next heartbeat when the slot
becomes free again.

**The `safe_path = .true.` bypass currently in the code is for development
testing only and must be replaced with a real call to `is_safe_script_path`
before production use.**

### Execution

Accepted scripts are run with:

```fortran
call execute_command_line('bash ' // trim(script_path), exitstat=exitstat)
```

The thread body sets `exit_code = exitstat`, `complete = .true.`, and
`nthr = 0` (releasing the slot) under the slot mutex before returning.

No other execution mechanism (system(), popen, …) should be used.  The
`bash` prefix must remain explicit: scripts are not assumed to be executable.

### Exit code handling

The exit code is recorded in `thread_args%exit_code` but is not currently
forwarded to the server.  A future enhancement may include it in a STATUS
reply.  Until then it must still be recorded correctly so it is available for
post-mortem inspection.

---

## Logging conventions

All log output goes to `write(*,…)` in the current implementation.  When the
process is submitted via a qsys backend that redirects stdout to a file (the
`outfile` argument in `exec_simple_prg_in_queue_async`), that file serves as
the worker's log.

Log-message prefix conventions:

| Context | Prefix |
|---|---|
| Main program | `'Worker: '` |
| `start_worker_thread` | `'Worker: '` |
| `worker_task_thread` body | `'Worker thread: '` |

All log lines must identify the relevant `job_id` or slot index when
applicable.  Integers must not be formatted as reals.

Future work: migrate to `write(logfhandle,'(A,...)')` with `int2str()` to
match the rest of the SIMPLE logging convention.

---

## Cleanup and teardown ordering

When the heartbeat loop exits (for any reason), cleanup must proceed in this
order:

1. For each slot `i = 1 … nthr`:
   a. Lock `threads(i)%mutex`.
   b. Read `threads(i)%started`.
   c. Unlock `threads(i)%mutex`.
   d. If `started`, call `pthread_join(threads(i)%thread, …)`.
   e. Under the mutex: reset `started = .false.`, `complete = .true.`,
      `nthr = 0`.
   f. Call `pthread_mutex_destroy(threads(i)%mutex)`.
2. `deallocate(threads)`.
3. Close `logfhandle` if it is not `output_unit` and is open.
4. Print git version and timing.

The mutex must not be destroyed before `pthread_join` completes for that
slot.  `pthread_mutex_destroy` must be called exactly once per slot.

---

## Invariants

The following invariants must be preserved by all future changes:

1. `simple_persistent_worker` has no knowledge of the job queue or task graph.  It
   executes scripts; the server decides what to dispatch.
2. `port` and `server` are mandatory arguments.  The process must call
   `stop 1` if either is absent or invalid.
3. `worker_uid` is computed once at startup and never changes for the
   lifetime of the process.
4. Each slot has its own mutex.  The main thread never holds more than one
   slot mutex simultaneously.
5. A slot is claimed (fields set, `started = .true.`) under the mutex,
   before `pthread_create`.
6. `pthread_join` must be called before the slot's `pthread_mutex_destroy`.
7. The `safe_path = .true.` bypass must not be used in production.
   `is_safe_script_path` must be called before `execute_command_line`.
8. The heartbeat loop must not call `stop`.  All exit paths go through the
   loop condition or the TERMINATE case.
9. `nthr = 0` in a slot is the authoritative signal that the slot is idle.
   No other field may be used as a free/busy indicator.
10. `execute_command_line('bash ' // trim(script_path), …)` is the only
    permitted execution mechanism.  The `bash` prefix must not be removed.
