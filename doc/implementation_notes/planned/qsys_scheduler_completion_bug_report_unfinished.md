# Scheduler completion: outstanding diagnostics and timeout

Status, 2026-11-02: reconstructed from the [completed bug report](../completed/qsys_scheduler_completion_bug_report_completed.md).
The earlier split linked here without supplying this file. Exit-code and lost-job handling are
implemented; this record retains the unfinished fourth step. No new scheduler run is claimed.

## Open work

- Report outstanding partition IDs periodically while a scheduler waits.
- Define a generous, overridable wall-clock backstop for a part that remains alive but hangs.
  Its default, configuration surface and disabled behavior require a maintainer decision;
  `SIMPLE_QSYS_TIMEOUT` in the old report was a suggestion, not an accepted key.
- Reconcile the unused `qsys_watcher_diag` with the controller's waiting loops rather than
  introducing another watcher. Current `schedule_jobs` and `schedule_array_jobs` still have no
  periodic diagnostic or timeout; the helper remains defined in `simple_qsys_funs`.

## Validation before closeout

Show bounded failure and useful diagnostics for a hung part, successful completion of a slow
part within the chosen budget, and the explicitly disabled timeout. Preserve nonzero-exit and
lost-part handling, including stale-status cleanup. Cover the local, batch and persistent-worker
paths that use the changed loops. Record commands, source revision, platform and logs here.
These checks remain unrecorded by this documentation change.
