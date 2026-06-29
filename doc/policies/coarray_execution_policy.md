# Coarray Execution Policy

This document records durable workflow contracts for SIMPLE's coarray queue
system backend. It is policy, not a line-by-line implementation map.

This document is the single policy authority for the coarray backend.

## 1. Core Model

The coarray backend is an in-process transport for partitioned qsys jobs.

It does not change scientific workflow semantics. The same per-partition
command descriptions that would normally be written into generated qsys scripts
are instead staged as in-memory command-line argument strings and passed to a
multi-image `simple_private_exec --coarray` invocation.

Each coarray image is a private-exec worker. It selects one or more partition
command strings, parses each selected command string in-process, and executes
the ordinary private commander path directly. No secondary `simple_private_exec`
process is spawned for a partition.

Unlike script-based qsys backends, the coarray backend must not write
`distr_simple_script_*` files for each partition. If a workflow needs generated
scripts as the durable unit of work, it should use a script-based backend.

## 2. Ownership

`simple_qsys_env.f90` owns backend selection, validation that `qsys=coarray`
is only used in a `USE_COARRAYS` build, qsys environment creation, and qsys
controller creation.

`simple_qsys_factory.f90` owns construction of the `coarray` qsys backend.

`simple_qsys_coarray.f90` owns the qsys submission object and launcher prefix,
including the `SIMPLE_COARRAY_SUBMIT_CMD` override and the default
`cafrun -np` launcher. It does not own the coarray dispatch protocol.

`simple_qsys_ctrl.f90` owns partition splitting, per-partition command
description augmentation, staged coarray job arguments, coarray submission, and
qsys job-state bookkeeping.

`production/simple_private_exec_driver.f90` owns callable private-exec dispatch
and the `--coarray` execution helper: protocol-argument parsing,
image-to-partition mapping, per-partition output redirection, in-process
private command execution, and per-image synchronization.

`simple_exec_helpers.f90` owns generic restarted, asynchronous, script-based,
and screen-based execution helpers. It must not depend on private commander
dispatch.

`production/simple_private_exec.f90` owns only executable-level dispatch. It
should remain thin: detect `--coarray`, delegate to `run_coarray_direct`, and
otherwise continue through the normal private commander path.

The private exec driver owns normal command-line parsing. CLI invocations use
`cline%parse_private`; coarray in-process invocations use
`cline%parse_private_line` on the staged partition command string. In both
cases, commanders still receive a normal `cmdline`, construct `params` through
`params%new(cline)`, and follow the ordinary workflow path.

## 3. Build And Backend Selection

`qsys=coarray` requires a `USE_COARRAYS` build. Non-coarray builds must fail
early during backend setup rather than reaching a later `--coarray` executable
path.

The coarray backend is selected through the normal qsys name mechanism. It is
not a scientific command option and must not leak into strategy, commander, or
domain-level branching except through normal qsys-controlled distributed
execution.

`SIMPLE_COARRAY_SUBMIT_CMD` may override the launcher prefix. The prefix is
expected to accept the image count immediately after it, matching the default
shape:

```text
cafrun -np <nimages>
```

## 4. Partition Metadata Contract

Partition metadata is propagated through serialized per-partition command-line
strings, not through coarray variables.

For each partition in `fromto_part`, `qsys_ctrl` must stage at least these keys
in the partition command description:

| Key | Meaning |
| --- | --- |
| `fromp` | first particle index handled by this partition |
| `top` | last particle index handled by this partition |
| `part` | global partition number |
| `nparts` | total number of partitions in the distributed job |

Optional partition-local keys such as `outfile` and entries from `part_params`
must be staged with the same semantics used by script-based qsys backends.

The in-process private command's parsed `parameters` object is the
authoritative source for these values. Downstream workflow code should consume
typed fields such as `params%fromp`, `params%top`, `params%part`, and
`params%nparts` after `params%new(cline)`.

Scientific workflow code must not infer `part`, `fromp`, `top`, or `nparts`
from `this_image()` or `num_images()`. Coarray image identity is a transport
detail, not the partition metadata contract.

## 5. Submission Command Protocol

The qsys controller submits the whole active partition range through one
multi-image command:

```text
<submit_cmd> <nimages> <exec_binary> --coarray <from_part> <to_part> <part_job_arg(from)> ... <part_job_arg(to)>
```

Every coarray image receives the same command line.

`from_part` and `to_part` describe the submitted partition range and are used
by the dispatcher for image-to-partition mapping. These protocol values are not
a replacement for the per-partition serialized job arguments.

The dispatcher may infer a local zero-padding width from `to_part` for
diagnostic log names. This inferred width is not the scientific `numlen`
contract. Any workflow-visible `numlen` value belongs in the serialized
per-partition job argument and is parsed by the private command through the
normal `params%new(cline)` path.

`nimages` is bounded by the active partition count and the controller's
available computing-unit count:

```text
nimages = max(1, min(ncomputing_units, to_part - from_part + 1))
```

Every staged job argument must be shell quoted by the submitter before it is
placed on the launcher command line.

## 6. Image Assignment

Each coarray image maps itself to partitions by image number:

```text
do ipart = from_part + this_image() - 1, to_part, num_images()
    part_job_ind = ipart - from_part + 1
    ...
end do
```

For example, with `from_part=1`, `to_part=8`, and `num_images()=3`:

| Image | Partitions |
| --- | --- |
| 1 | 1, 4, 7 |
| 2 | 2, 5, 8 |
| 3 | 3, 6 |

Images that receive no partition simply participate in the final synchronization
after their loop has no iterations.

The `part_job_ind` value selects the serialized command string for `ipart` from
the command-line partition-job argument list. That selected command string is
what carries the partition's `fromp`, `top`, `part`, and `nparts` values.

## 7. In-Process Private Execution

For each assigned partition, the coarray image runs the staged private command
in the current image:

```text
redirect stdout/stderr to simple_private_exec_coarray_part_<part>.out
cline%parse_private_line(part_job_args(part_job_ind))
dispatch private commander
restore stdout/stderr
```

There is no nested shell command and no MPI/PMI environment unsetting in this
path, because there is no child process to protect from the parent coarray/MPI
runtime.

The staged command string has the same key-value shape used by script-based
qsys backends. The coarray protocol parser only selects the string. Scientific
parameters are interpreted by the normal private command-line parser and
commander code.

Per-partition output is appended to:

```text
simple_private_exec_coarray_part_<part>.out
```

## 8. Completion And Error Semantics

The coarray qsys submission is synchronous from the controller's perspective.
`submit_coarray_jobs` marks all active jobs submitted, runs the multi-image
coarray command, and marks all active jobs done only after the command exits
successfully.

If an in-process private command fails, the owning coarray image must raise a
hard failure and the aggregate coarray submission must fail. The per-partition
log path should be preserved in diagnostics whenever possible.

All images must participate in the final `sync all`. A synchronization failure
is a hard failure of the coarray dispatch run.

For `qsys=coarray`, per-partition filesystem completion sentinels such as
`JOB_FINISHED_*` are not the completion authority. Worker-side completion
declaration goes through `qsys_declare_part_finished` / `qsys_job_finished`,
but the coarray branch deliberately avoids touching filesystem sentinels.

Inside the coarray run, each image records completion of the partitions it
executes and all images participate in a final coarray synchronization. That
in-memory declaration plus `sync all` is the readiness barrier for future
in-runtime assembly or reduction work. The parent coarray command's successful
return remains the outer qsys-controller signal that allows
`submit_coarray_jobs` to mark the active partition range done.

Aggregate coarray-submission failures should direct the user to
`SIMPLE_SUBPROC_OUTPUT` and the `simple_private_exec_coarray_part_*.out` logs.

## 9. Architectural Rules

Keep coarray support at the qsys and private-exec driver boundary. Do not place
`this_image()`, `num_images()`, or coarray dispatch rules in commanders,
strategies, numerical kernels, or domain objects.

Keep `simple_private_exec` thin. New `--coarray` protocol code belongs in
`production/simple_private_exec_driver.f90` or the qsys layer, not in the
production executable entry point and not in generic execution helpers.

Keep `simple_qsys_coarray.f90` focused on qsys backend identity and launcher
configuration. It should not grow partition-staging or private-exec dispatch
logic.

Keep `stage_coarray_jobs` semantically aligned with `prep_part_jobs`. If a new
per-partition key is required for script-based execution, decide explicitly
whether the coarray staging path requires the same key.

Do not introduce ad hoc command-line keys for coarray-only behavior unless they
also have typed parameter representation where workflow code consumes them.
After `params%new(cline)`, downstream code should use typed `params` fields,
not raw `cmdline` lookups.

Do not use coarray image identity as a substitute for `part`, `fromp`, `top`,
or `nparts`. The same partition command string should behave the same way
whether it is launched through a script backend or through the coarray backend.

## 10. End-To-End Trace

The intended control flow is:

1. `qsys_env` validates and constructs the `coarray` backend.
2. `qsys_ctrl%prep_part_jobs` detects `qsys_coarray` and calls
   `stage_coarray_jobs`.
3. The shared partition-job helper augments the base job description with
   `fromp`, `top`, `part`, `nparts`, and any partition-local keys.
4. The augmented description is serialized with `chash2str()` and stored in
   `coarray_job_args(ipart)`.
5. `submit_coarray_jobs` launches one multi-image `simple_private_exec
   --coarray` command containing the full active partition range and all staged
   job strings.
6. Every image receives the same command line and selects partitions using
   `this_image()` and `num_images()`.
7. For each selected partition, the image redirects stdout/stderr to that
   partition's coarray log and calls the private exec driver on the staged job
   string.
8. The in-process private command parses the serialized command line into
   `params` and runs the ordinary private commander path.
9. All images synchronize. The controller marks the jobs done only after the
   coarray command exits successfully.

## 11. Current Limits And Follow-Up

Coarray job arguments are staged in memory. They are not durable restartable
script files.

The coarray backend currently has no independent per-partition completion-file
authority. The parent coarray command return status is the completion authority.

The staging path and script-generation path share the same per-partition
job-description augmentation helper. Changes to `fromp`, `top`, `part`,
`nparts`, `outfile`, or `part_params` staging policy must go through that shared
helper to avoid drift between script and coarray backends.

The in-process path intentionally emulates the old shell redirection contract in
the first pass. Future cleanup should decide whether per-partition logs remain
fd-level stdout/stderr redirects or become structured log sinks owned by the
driver.

The current coarray backend is designed for partitioned private-exec workflows,
but the in-process worker model is the intended foundation for future coarray
reductions of partial class averages, reconstructions, and similar distributed
state.
