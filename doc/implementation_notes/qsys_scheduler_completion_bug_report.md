# Bug report: qsys job-completion detection

**Component:** `src/utils/qsys/simple_qsys_ctrl.f90` (with one supporting note in `src/utils/qsys/simple_qsys_funs.f90`)
**Found while:** investigating a distributed `calc_pspec` run that appeared to hang with idle CPUs
**Status:** no code change made — the scheduler is heavily used in production, so this is reported rather than patched
**Baseline:** `master` @ `ec5527df0`

There are two independent issues. Issue 1 is a concrete defect with a confirmed
reproduction. Issue 2 is a missing-capability gap that produces the silent hang.

---

## Issue 1 — recorded exit status is `tee`'s, not the executable's

**Severity:** medium. Silent misclassification; no crash, no error message.

### Where

`generate_script_2`, exit-code block at `simple_qsys_ctrl.f90:588-592`:

```fortran
! exit code
if( present(exit_code_fname) )then
    write(funit,'(a)') ''
    write(funit,'(a)') 'echo $? > '//exit_code_fname%to_char()
endif
```

The output branch immediately above it (`:581-587`) has two forms:

```fortran
if( present(outfile) )then
    write(funit,'(a)') ' > '//outfile%to_char()//' '//STDERR2STDOUT   ! redirect, no pipe
else
    write(funit,'(a)') ' '//STDERR2STDOUT//' | tee -a '//SIMPLE_SUBPROC_OUT   ! pipe
endif
```

In the piped form `$?` is the exit status of the **last** command in the
pipeline, i.e. `tee`, which essentially always succeeds. The executable's status
is discarded.

### Who reaches it

`schedule_streaming` at `simple_qsys_ctrl.f90:1032-1033` calls
`generate_script_2` with `exit_code_fname` but **without** `outfile`, so the
piped branch is taken. This is the production streaming path.

### Reproduction

The generated script has this shape; running it directly:

```bash
#!/bin/bash
sh -c "exit 42" 2>&1 | tee -a SIMPLE_SUBPROC_OUTPUT
echo $? > EXIT_CODE_JOB_1        # writes 0
echo ${PIPESTATUS[0]}            # writes 42
```

Observed: `EXIT_CODE_JOB_1` contains `0`. Verified for statuses 1, 42, 127
and 139 — all recorded as `0`.

### Consequence

The consumer is `update_queue`, streaming branch, `simple_qsys_ctrl.f90:933`:

```fortran
if( (exit_code == 0) .and. (.not. err) .and. file_exists(self%jobs_done_fnames(ipart)) ) then
    ! → done stack
else
    ! → fail stack
end if
```

Because `exit_code` is always `0`, the first term is always true and
classification collapses onto the sentinel-file test alone. Two effects:

1. **A job that fails after writing its `JOB_FINISHED_*` sentinel is recorded as
   successful** and goes to the done stack. The sentinel is written from inside
   the process (`qsys_declare_part_finished`) before the process exits, so a
   failure during teardown/exit falls in this window.
2. **Diagnostic value of `EXIT_CODE_JOB_*` is lost.** Every file contains `0`, so
   the files cannot distinguish `127` (bad binary path), `139` (segfault), `137`
   (OOM kill) or `1` (`THROW_HARD` → `error stop 1`).

Note the common case — a job that dies *before* writing its sentinel — still
lands in the fail stack, because the sentinel test fails. So this is not "all
streaming failures are missed"; it is a narrower window plus a total loss of
failure-cause information.

### Suggested fix

One-line: use `${PIPESTATUS[0]}` in the piped branch, keep `$?` in the `outfile`
branch where there is no pipe. Scripts are `#!/bin/bash`, so `PIPESTATUS` is
available. `generate_scripts_subprojects` avoids the problem a different way
(`set -o pipefail` at `:319` plus an `if`/`else`), which is also fine.

---

## Issue 2 — batch scheduling cannot observe job failure and has no timeout

**Severity:** high for usability. This is what produces the "hangs with idle
CPUs" symptom.

### Where

`schedule_jobs`, `simple_qsys_ctrl.f90:962-971`:

```fortran
subroutine schedule_jobs( self )
    class(qsys_ctrl), intent(inout) :: self
    do
        if( all(self%jobs_done) ) exit
        call self%update_queue
        call self%submit_scripts
        call self%service_persistent_worker_warmup()
        call sleep(SHORTTIME)          ! SHORTTIME = 1
    end do
end subroutine schedule_jobs
```

The only exit condition is `all(jobs_done)`. There is no iteration cap, no
wall-clock budget, no exit-status inspection and no liveness check.

`update_queue`'s batch branch (`:947-956`) polls `file_exists(JOB_FINISHED_<part>)`
and nothing else. `generate_script_1` (`:509-547`) emits no exit-code file at
all — unlike `generate_script_2`, the batch path has no `exit_code_fname`
equivalent, even though `jobs_exit_code_fnames(:)` is allocated for every
partition in `new` (`:173`, `:189`).

### Failure mode

Any partition worker that dies without writing `JOB_FINISHED_<part>` makes the
master loop forever at ~0% CPU with no message:

- `THROW_HARD` anywhere in the worker → `simple_error.f90:26` → `error stop 1`
- segfault, OOM kill
- stale `simple_path` in the project `compenv` segment → bash cannot find
  `simple_private_exec` → exit 127 before SIMPLE runs at all

There is no retry either: `submit_scripts` sets `jobs_submitted(ipart) = .true.`
at `:797` before launching and never clears it within a scheduling round, so a
dead job is never resubmitted.

The worker's error text *is* captured — it goes to `SIMPLE_SUBPROC_OUTPUT` via
`tee` — but nothing surfaces it, and on the local backend there is no per-part
log directory to look in (`qsys_local` creates no `stderrout/`).

### Already-written, unused remedy

`qsys_watcher_diag` in `simple_qsys_funs.f90:124-164` implements exactly the
missing behaviour — periodic reporting of which sentinel files are still
missing, plus a hard-stop timeout. It is referenced **nowhere** in the tree
(0 call sites outside its own definition).

### Suggested direction

Roughly symmetric with what the streaming path already does:

1. Have `generate_script_1` write `EXIT_CODE_JOB_<part>` (using
   `${PIPESTATUS[0]}`, per Issue 1), and delete any stale file when it
   regenerates the script.
2. In `update_queue`'s batch branch, when a submitted partition has no sentinel
   but has a readable **non-zero** exit code, flag it as failed. Deliberately do
   *not* treat "zero exit code, no sentinel yet" as failure — the sentinel is
   written before the shell records the status, so it can legitimately lag; and
   `read_exit_code` reports `err` while the file is mid-write, which should mean
   "retry next poll", not "failed".
3. Abort the loop on any flagged failure, listing the partition, its exit status
   and its script path, and pointing at `SIMPLE_SUBPROC_OUTPUT`.
4. Report outstanding partitions periodically, and add a wall-clock backstop.
   The backstop should be generous (a dead job is caught deterministically by
   step 2, so the timeout only has to bound jobs that hang without exiting) and
   overridable — an env var such as `SIMPLE_QSYS_TIMEOUT` with `0` to disable
   would fit how the rest of qsys reads configuration.

Points 3 and 4 apply equally to `schedule_subproject_jobs` (`:974`) and
`schedule_array_jobs` (`:986`), which have the same unbounded loop.

### Risk note

Changing exit-status handling alters when jobs are declared failed, which can
turn currently-silent partial failures into hard stops. That is the intent, but
it means the change wants a deliberate rollout and testing across backends
(`local`, `slurm`, persistent worker) rather than being folded into unrelated
work — which is why it was not patched here.

---

## Triage checklist for a live hang

1. Read `SIMPLE_SUBPROC_OUTPUT` in the run directory. Worker `ERROR!` lines are
   tagged with the partition number.
2. Compare `ls JOB_FINISHED_*` against `ls distr_simple_script_*` to see which
   partitions never reported.
3. `ps` for `simple_private_exec` — none alive means the master is waiting on
   dead workers rather than slow ones.
4. Check the `simple_path` and `qsys_name` values in the project's `compenv`
   segment; `qsys_env` reads the backend from the project, not the command line
   (`simple_qsys_env.f90:158-164`), so a stray `_worker` suffix or a stale
   install path silently changes behaviour.
