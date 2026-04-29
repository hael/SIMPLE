# Current SIMPLE Worker Backend Policy

## Status

This document describes the worker backend implementation currently present in
the repository. It is intentionally descriptive rather than aspirational.

The current implementation is a `_worker` overlay on the existing qsys system:

- `local_worker`
- `slurm_worker`
- `lsf_worker`
- `pbs_worker`
- `sge_worker`

The overlay starts persistent `simple_worker` processes and routes generated
SIMPLE job scripts to them through a TCP worker server. The existing qsys
scheduler and filesystem completion files still own job-level completion.

---

## What Exists Today

The implementation consists of four concrete pieces.

### 1. Worker overlay selection

Implemented in:

- `src/utils/qsys/simple_qsys_env.f90`
- `src/utils/qsys/simple_qsys_factory.f90`

`simple_qsys_env%new()` reads the configured qsys name, detects `_worker`, and
strips the suffix before constructing the base backend:

```fortran
qsnam = self%qdescr%get('qsys_name')
if( qsnam%has_substr(string('_worker')) ) then
    self%use_workers = .true.
    qsnam = qsnam%substr_remove(string('_worker'))
end if
call self%qsys_fac%new(qsnam, self%myqsys)
```

`simple_qsys_factory` supports a separate internal `worker` qsys type:

```fortran
case('worker')
    allocate(qsys_worker :: self%qsys_base_type)
```

Current behavior:

- The base backend is still created first.
- The base backend is used to launch persistent `simple_worker` processes.
- After workers are launched, `qsys_env` replaces its active qsys object with
  `qsys_worker`, so later script submissions are queued to the worker server.

Current caveat:

- The local variable `qsnam` is stripped, but `qdescr('qsys_name')` is not
  currently normalized to the stripped base name. This matters because script
  generation still tests `qsys_name == 'local'` for local chmod behavior.

---

### 2. Worker server

Implemented in:

- `src/utils/qsys/simple_qsys_worker.f90`
- `src/utils/qsys/worker/simple_qsys_worker_server.f90`
- `src/utils/comm/simple_ipc_tcp_socket.f90`

`simple_qsys_worker.f90` provides a module-level server pointer:

```fortran
type(qsys_worker_server), pointer, public :: worker_server => null()
```

`simple_qsys_env%new()` allocates and starts this server lazily:

```fortran
if( .not. associated(worker_server) ) then
    allocate(worker_server)
    call worker_server%new(nthr_workers)
    port     = worker_server%get_port()
    host_ips = worker_server%get_host_ips()
    n_workers = max(1, params%ncunits)
    if( params%workers > 0 ) n_workers = params%workers
    call self%start_workers(n_workers, nthr_workers, port, host_ips)
end if
```

The server:

- opens a TCP listener
- receives heartbeat messages from workers
- keeps a fixed-size worker registry
- keeps three fixed-size task queues: high, normal, low
- dispatches at most one task in response to a heartbeat

Task dispatch happens in `worker_listener_thread`:

```fortran
if( status%tasks_priority_norm(i)%job_id > 0 .and. &
    .not. status%tasks_priority_norm(i)%submitted ) then
    if( heartbeat_msg%nthr_total - heartbeat_msg%nthr_used < &
        status%tasks_priority_norm(i)%nthr ) cycle
    status%tasks_priority_norm(i)%submitted = .true.
    dispatch_task = status%tasks_priority_norm(i)
    has_task = .true.
end if
```

Current behavior:

- The server sends a task only if the worker reports enough free thread
  capacity.
- The server marks the task as `submitted`.
- The server does not currently receive completion messages from workers.

Current caveat:

- Submitted tasks remain in the server queue because `end_time` is never set by
  a completion acknowledgement. The existing qsys scheduler may still finish
  because it watches filesystem completion files, but the server queue itself
  does not retire dispatched tasks.

---

### 3. Worker process launch

Implemented in:

- `src/utils/qsys/simple_qsys_env.f90`
- `production/simple_worker.f90`

`simple_qsys_env%start_workers()` launches persistent worker processes through
the base backend:

```fortran
executable = 'simple_worker'
call cline%set('nthr',      int2str(nthr))
call cline%set('port',      int2str(port))
call cline%set('server',    host_ips%to_char())
call cline%set('worker_id', int2str(iworker))
call self%exec_simple_prg_in_queue_async(cline, script_name, outfile, executable)
```

The worker executable is built from:

- `production/simple_worker.f90`

`simple_worker` accepts key-value command-line arguments:

- `nthr`
- `port`
- `server`
- `worker_id`

It creates a stable process identifier:

```fortran
worker_uid = trim(hostname) // '_' // trim(int2str(int(pid)))
```

It then sends heartbeat messages to the server:

```fortran
heartbeat_msg%worker_id      = worker_id
heartbeat_msg%worker_uid     = trim(worker_uid)
heartbeat_msg%heartbeat_time = int(c_time(0_c_long))
heartbeat_msg%nthr_total     = nthr
heartbeat_msg%nthr_used      = get_nthr_used()
```

Current behavior:

- A worker is a persistent process.
- It advertises local thread-slot capacity.
- It receives at most one task per heartbeat reply.
- It runs tasks in local pthreads.

---

### 4. Script submission through workers

Implemented in:

- `src/utils/qsys/simple_qsys_ctrl.f90`
- `src/utils/qsys/worker/simple_qsys_worker_message_task.f90`

Existing SIMPLE distributed workflows still generate job scripts in
`qsys_ctrl`. When the active qsys object is `qsys_worker`, the script is not
submitted directly to `nohup`, SLURM, LSF, PBS, or SGE. Instead, the script path
is queued as a worker-server task:

```fortran
class is (qsys_worker)
    call self%submit_task_to_worker_server(script_name, self%nthr_worker)
    cycle
```

For single-script submission:

```fortran
type is (qsys_worker)
    call self%submit_task_to_worker_server(filepath(string(CWD_GLOB),script_name), self%nthr_worker)
    return
```

The task message stores:

```fortran
task_msg%nthr        = nthr
task_msg%script_path = script_path%to_char()
```

The worker receives the task and executes:

```fortran
call execute_command_line('bash ' // trim(thread_args%task_msg%script_path), exitstat=exitstat)
```

Current behavior:

- The generated script remains the unit of scientific work.
- The worker backend changes how scripts are dispatched.
- Completion is still observed through the script's existing filesystem side
  effects, such as `JOB_FINISHED*` and `EXIT_CODE_JOB_*`.

Current caveat:

- `submit_task_to_worker_server` currently logs enqueue failure but does not
  return a success/failure value to the scheduler.
- `submit_scripts` marks jobs submitted before it knows whether worker enqueue
  succeeded.

---

## Message Protocol Currently In Use

Implemented in:

- `src/utils/qsys/worker/simple_qsys_worker_message_types.f90`
- `src/utils/qsys/worker/simple_qsys_worker_message_heartbeat.f90`
- `src/utils/qsys/worker/simple_qsys_worker_message_task.f90`
- `src/utils/qsys/worker/simple_qsys_worker_message_status.f90`
- `src/utils/qsys/worker/simple_qsys_worker_message_terminate.f90`

The protocol currently defines four message types:

| Constant | Direction | Current use |
|---|---|---|
| `WORKER_TERMINATE_MSG` | server to worker | tell worker to shut down |
| `WORKER_HEARTBEAT_MSG` | worker to server | advertise liveness and thread capacity |
| `WORKER_TASK_MSG` | server to worker | hand worker a script path |
| `WORKER_STATUS_MSG` | server to worker | idle or error reply |

Every concrete message type implements its own `serialise` routine using
`sizeof(self)` and `transfer`.

Current caveat:

- There is no completion/status message from worker to server after a script
  finishes.

---

## Asynchronous Behavior

The current worker backend is asynchronous as a dispatch transport.

That means:

- `simple_worker` processes persist
- tasks are queued to the worker server
- workers poll by heartbeat
- workers execute task scripts in local pthreads

It is not currently a replacement for the existing stream scheduler.

The existing scheduling loops still own workflow completion:

- `qsys_ctrl%schedule_jobs`
- `qsys_ctrl%schedule_subproject_jobs`
- `qsys_ctrl%schedule_streaming`

Those loops still watch the same files they watched before:

- `JOB_FINISHED*`
- `EXIT_CODE_JOB_*`
- stream done/fail stacks

So the current model is:

```text
existing SIMPLE scheduler
    -> generated bash scripts
        -> worker-server task queue
            -> persistent simple_worker process
                -> bash script execution
                    -> legacy filesystem completion files
```

The worker server is not currently the source of truth for task completion.

---

## Current Configuration Knobs

Implemented in:

- `src/main/simple_parameters.f90`

Relevant parameters:

```fortran
integer :: workers=0
integer :: worker_nthr=0
```

Current resolution in `simple_qsys_env%new()`:

```fortran
nthr_workers = params%nthr
if( params%worker_nthr > 0 ) nthr_workers = params%worker_nthr
if( present(qsys_nthr)     ) nthr_workers = qsys_nthr

n_workers = max(1, params%ncunits)
if( params%workers > 0 ) n_workers = params%workers
```

Meaning:

- `worker_nthr` controls the thread capacity advertised by each worker process.
- `workers` overrides the number of persistent worker processes launched.
- if `workers` is unset, the worker count defaults to `params%ncunits`.

---

## Concrete Usage Example

Any code path that constructs `qsys_env` and schedules qsys scripts can use the
worker overlay by setting the qsys name to a `_worker` variant.

For example, conceptually:

```text
qsys_name=local_worker
nparts=8
ncunits=4
nthr=8
workers=4
worker_nthr=8
```

Current execution shape:

1. `qsys_env%new()` sees `local_worker`.
2. It strips `_worker` and constructs the base `local` backend.
3. It starts the TCP worker server.
4. It launches `workers` persistent `simple_worker` processes through the
   base `local` backend.
5. It replaces the active qsys object with `qsys_worker`.
6. Later generated partition scripts are queued to the worker server instead
   of being launched directly.
7. Each `simple_worker` heartbeats and receives scripts when it has enough
   free thread slots.
8. The scheduler still waits for legacy completion files.

The same overlay shape is intended for `slurm_worker`, `lsf_worker`,
`pbs_worker`, and `sge_worker`, with the base backend responsible for launching
the persistent worker processes.

---

## What The Current Implementation Does Not Yet Do

This section is deliberately descriptive. These are current implementation
facts, not future requirements.

### No worker-to-server completion acknowledgement

Workers record local script exit status, but they do not send a completion
message back to the worker server.

Consequence:

- qsys scheduling can still complete via filesystem markers
- the worker-server task queue does not know that submitted tasks finished

### Enqueue failure is not propagated to the scheduler

`submit_task_to_worker_server` logs missing-server or queue-full conditions but
does not return failure to `submit_scripts`.

Consequence:

- a failed enqueue can still be marked submitted by `qsys_ctrl`
- the scheduler may wait for a completion file that no worker will create

### Server lifetime is global and not finally owned

`worker_server` is a module-level pointer. It is allocated in `qsys_env%new()`.
`qsys_env%kill()` currently releases qsys objects and scripts but does not call
`worker_server%kill()` or deallocate the global pointer.

Consequence:

- reuse across qsys-env instances is possible
- final teardown ownership is not explicit in the current code

### Production script-path validation is disabled

`production/simple_worker.f90` contains `is_safe_script_path`, but the dispatch
path currently bypasses it:

```fortran
!  safe_path = is_safe_script_path(task_msg%script_path)
safe_path = .true. ! for testing
```

Consequence:

- the current worker accepts any non-empty script path handed to it

### Worker overlay does not change scientific workflow semantics

The worker backend does not alter refine3D, cluster2D, streaming, or any other
scientific algorithm directly.

It only changes how generated qsys scripts are dispatched.

---

## Current Terminology

To describe the current code precisely:

| Term | Current meaning |
|---|---|
| distributed part | one generated qsys script / partition-level work unit |
| worker process | one persistent `simple_worker` executable |
| worker server | the TCP dispatcher in `simple_qsys_worker_server` |
| worker overlay | the `qsys_worker` qsys object used after startup |
| thread slot | local task capacity reported by a worker heartbeat |

The current implementation does not yet make `Nparts`, `Nworkers`, and `Nthr`
a fully independent scheduling model everywhere. It does, however, introduce
the pieces:

- `nparts` controls generated distributed scripts
- `workers` controls persistent worker processes when set
- `worker_nthr` controls worker thread capacity when set

---

## Current Review Checklist

When reviewing changes to the current worker backend, check:

1. Does `_worker` still strip to a valid base backend before factory
   construction?
2. Does worker launch still go through the base backend?
3. Does normal task submission still go through `qsys_worker` only after worker
   startup?
4. Does the scheduler still have a valid filesystem completion path?
5. Can failed enqueue be observed, retried, or made fatal?
6. Can submitted tasks accumulate forever in the worker-server queue?
7. Is the worker server eventually shut down by an explicit owner?
8. Is script-path validation enabled for production use?

---

## Code Review: Make The Persistent-Worker Model Visible

The current implementation works as a qsys `_worker` dispatch layer, but the
code does not yet make the persistent-worker model obvious enough. The most
important cleanup is to make the distinction visible in names and ownership:

```text
base qsys backend
    launches persistent workers

persistent-worker dispatch
    queues generated scripts to those workers

existing SIMPLE scheduler
    still tracks script completion
```

The preferred short term is `persistent_worker`.

Avoid introducing longer names such as `persistent_worker_overlay` unless the
extra precision is genuinely needed. Also avoid relying on the bare word
`worker` where it could mean either a distributed part, a qsys worker object, or
a persistent process.

### Suggested naming changes

Current names:

- `qsys_worker`
- `simple_qsys_worker`
- `worker_server`
- `start_workers`
- `submit_task_to_worker_server`
- `workers`
- `worker_nthr`

Preferred names:

- `persistent_worker`
- `simple_persistent_worker`
- `persistent_worker_server`
- `start_persistent_workers`
- `submit_script_to_persistent_worker`
- `persistent_workers`
- `persistent_worker_nthr`

These do not all need to become command-line names immediately. Backward
compatible aliases are fine. The important point is that internal code should
tell the reader that this is a pool of persistent processes.

### Suggested ownership cleanup

The current `qsys_env` switches `self%myqsys` from the base backend to
`qsys_worker` after launching workers:

```fortran
self%myqsys => self%myqsys_worker
```

That makes it harder to see the two roles:

- the base backend launches persistent workers
- the persistent-worker path dispatches generated scripts

A clearer structure would keep those roles separately named, for example:

```fortran
class(qsys_base), pointer :: base_qsys => null()
class(qsys_base), pointer :: dispatch_qsys => null()
```

or, more explicitly:

```fortran
class(qsys_base), pointer :: launch_backend => null()
class(qsys_base), pointer :: script_dispatch => null()
```

Then `launch_backend` remains `local`, `slurm`, `pbs`, etc., while
`script_dispatch` can be either direct qsys submission or persistent-worker
submission.

### Suggested runtime object

The module-level `worker_server` pointer makes the server shared, but it hides
who owns final teardown and which worker topology it represents. A small runtime
object would make this clearer:

```fortran
type persistent_worker_runtime
    type(qsys_worker_server), pointer :: server => null()
    integer :: nworkers = 0
    integer :: nthr_per_worker = 0
    type(string) :: launch_backend
end type persistent_worker_runtime
```

This object would make it easier to check whether a reused worker server is
compatible with a new qsys-env request.

### Suggested log wording

Worker launch logs should say what is happening in operational terms:

```text
STARTING PERSISTENT WORKERS
launch backend: local
nworkers: 4
nthr_per_worker: 8
```

For a cluster-backed run:

```text
STARTING PERSISTENT WORKERS
launch backend: slurm
nworkers: 8
nthr_per_worker: 16
```

This makes clear that the same persistent-worker model applies whether the
workers are local workstation processes or processes launched through a cluster
backend.

### Required semantic cleanup

The naming cleanup should be paired with one semantic cleanup: task lifetime in
the server queue must be explicit.

Choose one model:

1. Dispatch transfers ownership to the worker, and the server removes the task
   immediately after dispatch. The existing filesystem completion files remain
   the source of truth.
2. Workers send completion acknowledgement back to the server, and the server
   retires tasks on acknowledgement.

The current halfway state, where the server marks tasks `submitted` but never
receives completion and never retires them, obscures the persistent-worker
model and will eventually fill the finite queues.
