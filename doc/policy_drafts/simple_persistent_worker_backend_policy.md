# SIMPLE_QSYS Worker Backend Policy

## Status

This document is a proposal / design record, not a description of fully
stabilised behavior.  It should be promoted to `doc/policies/` once the
worker backend has been exercised in production and the interfaces are
considered stable.

## Scope

This document defines the architectural policy for the `_worker`-suffixed
qsys backend in SIMPLE.

It covers:

- the meaning of the `_worker` suffix in a qsys name
- the relationship between the base qsys backends and the worker overlay
- the TCP worker server and its lifecycle
- the message protocol between server and worker processes
- thread-safety obligations
- the `nthr_worker` configuration parameter
- the debug logging flag
- invariants that all future changes must preserve

This document does not cover the internal implementation details of the
base qsys backends (local, SLURM, LSF, PBS, SGE), the IPC socket layer,
or the probabilistic alignment machinery.

---

## Core design rule

The `_worker` suffix is a thin overlay, not a replacement backend.

- A qsys name without the suffix (`local`, `slurm`, …) continues to use the
  standard shell-script/batch-script submission path.
- Appending `_worker` to any supported backend name opts into the TCP worker
  backend.  The suffix is stripped before the factory resolves the base
  backend, so the underlying submission mechanism (local shell, SLURM sbatch,
  …) is still used — but only to launch the persistent `simple_persistent_worker`
  processes.  Once those workers are running, all task dispatch goes through
  the TCP server.

In other words:

- the base qsys backend launches workers once at startup
- the worker qsys backend dispatches tasks to those workers at runtime

That separation must be preserved.  Neither layer should absorb
responsibilities that belong to the other.

---

## Architecture overview

```
simple_qsys_env
  │
  ├─ base qsys (local / slurm / …)   ← used only to launch simple_persistent_worker processes
  │
  ├─ simple_persistent_worker_server        ← TCP accept-loop, heartbeat, task queue
  │     └─ listener thread (pthread)
  │
  └─ simple_qsys_worker (qsys_base)  ← submit_script → submit_task_to_worker_server
        └─ simple_qsys_ctrl           ← script generation + scheduling
```

The `worker_server` module-level pointer in `simple_qsys_worker` is
allocated and managed by `simple_qsys_env`, not by `simple_qsys_worker`
itself.  `simple_qsys_worker` accesses the pointer by USE association.

---

## The `_worker` suffix

### Detection and stripping

`simple_qsys_env%new()` reads the resolved qsys name into a local `string`
variable and calls `has_substr('_worker')` on it.  If the suffix is present:

1. `self%use_workers` is set `.true.`
2. The suffix is stripped from the local variable with `substr_remove`
3. The stripped name is passed to the base qsys factory

This must happen before the first call to `qsys_fac%new()`.  Do not move
or reorder that block.

### Naming convention

Any base backend name may receive the suffix:

| User-supplied name | Base backend used to launch workers |
|---|---|
| `local_worker` | `local` |
| `slurm_worker` | `slurm` |
| `lsf_worker` | `lsf` |
| `pbs_worker` | `pbs` |
| `sge_worker` | `sge` |

No special-casing per backend is needed in `simple_qsys_env`.

---

## Worker server lifecycle

### Startup

`simple_qsys_env%new()` starts the server lazily: it is created only if
`worker_server` is not yet associated.  This means a second call to `new()`
in the same session (e.g. re-initialisation after a stream checkpoint) reuses
the running server.

The startup sequence, in order, must be:

1. Allocate `worker_server`
2. Call `worker_server%new(nthr_workers)` — binds the TCP socket, spawns
   the listener thread, stores the ephemeral port
3. Query `worker_server%get_port()` and `worker_server%get_host_ips()`
4. Call `self%start_workers(n_workers, nthr_workers, port, host_ips)` —
   submits one `simple_persistent_worker` process per computing unit via the base
   qsys backend

No task dispatch may occur before step 4 completes.

### Teardown

`worker_server` teardown is the caller's responsibility (currently the top-
level SIMPLE program).  `simple_qsys_env%kill()` does not deallocate
`worker_server`.  Rationale: the server must outlive individual `qsys_env`
instances within the same session.

To shut down cleanly:

1. Send a `WORKER_TERMINATE_MSG` to each connected worker (handled inside
   `persistent_worker_server%kill()`)
2. Join the listener thread
3. Destroy the mutex
4. Close and zero `self%port`

That order must not be changed.  In particular, the mutex must still exist
when the thread-join happens, and `self%port = 0` must only be written after
`ipc_socket%kill()` returns.

---

## Message protocol

Four fixed-size C-interoperable message types are defined in
`simple_persistent_worker_message_types.f90`:

| Constant | Value | Direction | Purpose |
|---|---|---|---|
| `WORKER_TERMINATE_MSG` | 1 | server → worker | orderly shutdown |
| `WORKER_HEARTBEAT_MSG` | 2 | worker → server | liveness + load |
| `WORKER_TASK_MSG`      | 3 | server → worker | task dispatch |
| `WORKER_STATUS_MSG`    | 4 | server → worker | idle / error reply |

### Serialisation contract

Every concrete message type must override the `serialise` procedure from
`qsys_persistent_worker_message_base`.  The override must:

- `deallocate` any existing buffer
- `allocate` the buffer using `sizeof(self)` where `self` is declared as
  the concrete type (not the base class)
- fill the buffer with `TRANSFER(self, buf)`

The base class `serialise` exists only as a fallback placeholder.  It must
never be called in production: it serialises by `sizeof` of the base fields
only, which silently truncates all extra fields.

### Message buffer size

The TCP buffer constant `TCP_BUFSZ = 1460` (one Ethernet MTU payload) is
intentionally conservative.  All message structs must fit within that limit.
If a new message type would exceed it, either reduce the payload or change
`TCP_BUFSZ` and update this policy.

---

## Thread-safety obligations

`persistent_worker_server` owns a single mutex (`self%mutex`) that protects:

- the worker registry (`worker_data` array)
- the task queue (`task_queue` array and `job_count`)
- the `l_debug` flag in `worker_data`

Rules:

1. The mutex is initialised in `new()` and destroyed in `kill()` after the
   listener thread has been joined.  No other procedure may init or destroy it.
2. `job_count` increment and `task%job_id` assignment must happen inside the
   mutex in `queue_task`.
3. `set_debug()` acquires the mutex before propagating the flag to all
   `worker_data` entries.
4. The listener thread body is a `bind(c)` subroutine; it must not call any
   Fortran I/O that is not thread-safe in the target runtime.

---

## `nthr_worker` configuration

`nthr_worker` is the number of OpenMP threads (CPU cores) each worker process
is permitted to use.  It is separate from the number of worker processes.

| Symbol | Default | Meaning |
|---|---|---|
| `qsys_ctrl%nthr_worker` | `1` | thread slots requested per job submission |
| `nthr_workers` (local in `qsys_env%new`) | `params%nthr` | resolved at `new()` call time |

Resolution order inside `qsys_env%new()` (highest priority wins):

1. `qsys_nthr` optional argument, if present
2. `params%nthr` (the `nthr=` command-line parameter)

The resolved value is passed to:

- `worker_server%new(nthr_workers)` — advertised thread capacity per worker
- `qscripts%new(…, nthr_worker=nthr_workers)` — controls the `nthr` field
  written into task dispatch messages

The `nthr_worker` field in `qsys_ctrl` must be reset to `1` in `kill()`.

---

## Debug logging flag

`persistent_worker_server` and `persistent_worker_data` each carry an `l_debug` logical
field (default `.false.`).

When `.false.`, only significant events are logged:

- errors
- TERMINATE message receipt
- task dispatch (job id, worker id)
- queue-full rejection
- out-of-range worker id
- changed `worker_uid` for a slot

When `.true.`, high-frequency events are also logged:

- every heartbeat received
- per-task "no tasks available" replies
- heartbeat-skip messages for busy workers

`l_debug` can be set before or after `new()`:

- Before `new()`: set `server%l_debug = .true.` before calling `new()`;
  `new()` copies it into `worker_data%l_debug`.
- After `new()`: call `server%set_debug(.true.)`; the setter acquires the
  mutex before writing to each `worker_data` entry.

Do not add unconditional high-frequency log lines to the listener thread
body.  Performance profiling showed that unbuffered Fortran I/O in the hot
path can stall the accept loop.

---

## Invariants

The following invariants must be preserved by all future changes:

1. The `_worker` suffix is stripped from the qsys name before it reaches
   the factory.  The factory must never see `_worker` as a backend name.
2. `worker_server` is allocated at most once per session.  The guard
   `if( .not. associated(worker_server) )` must not be removed.
3. `serialise` overrides must use `sizeof(self)` where `self` is declared
   at the concrete type, not the base type.
4. The mutex is init in `new()`, destroyed in `kill()` after thread join.
   Never init or destroy it elsewhere.
5. `job_count` increment and `task%job_id` assignment are atomic with
   respect to the mutex in `queue_task`.
6. `self%port = 0` is written only after `ipc_socket%kill()` returns.
7. `nthr_worker` default is `1`; it is reset to `1` in `kill()`.
8. High-frequency log lines in the listener thread body are gated behind
   `if( worker_data%l_debug )`.
