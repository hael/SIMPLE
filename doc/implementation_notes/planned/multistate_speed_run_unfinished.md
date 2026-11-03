# Multi-state speed run: historical closeout

Status, 2026-11-02: reconstructed from the [completed ms_speed record](../completed/multistate_speed_run_completed.md).
Its 2026-10-20 status left the Phase 2 acceptance stop and Phase 4 closeout open. This file
preserves those obligations; it does not assert that the old assignment policy is still current.

## Open decisions and evidence

- Dispose of the two noisy multi-threaded Phase 2 checks, using the recorded baselines and
  measurements. The implementation and one-thread identity evidence were already delivered.
- Reconcile the Phase 4 final comparisons for cases A, MS and M against Phase 0, the simulated
  workflow/FLEX/pyramid guards, and the requested table of wall time, per-state resolution and
  populations. The case definitions and original acceptance rules remain in the completed record.
- Reconcile the policy/week-plan updates and final report with the later speed and regression
  work before creating another report or rerunning an old experiment.

The current `ref_assign` source explicitly replaces the balanced frontier with posterior draws.
Therefore the old frontier's byte-identity contract cannot simply be imposed on today's search.
Mark each remaining obligation as supported by a dated result, superseded by a named later
decision, or still open under the [current testing plan](testing_and_benchmarking_plan.md).
Do not restore the retired policy to satisfy this historical checklist.
