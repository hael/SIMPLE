# Offline stream emulation: remaining validation

Status, 2026-11-02: reconstructed from section 9 of the [completed implementation record](../completed/emulate_solve3D_stream_completed.md).
The supported follow-up is `solve3D_chunks`; `emulate_solve3D_stream` has been removed, as the
completed record already states. No new run is claimed.

## Checks still to reconcile

- Repeat the ribosome case with whole-stack cuts, retaining `SOLVE3D_BENCH.txt` from each run.
- Exercise a base-only run, a split near the 30,000-particle boundary, and a cut inside a stack.
- Check the published project's solution and manifest registration.

The fourth historical check was an `emulate_solve3D_stream` rollback case. It is retired with
that program; it does not require restoring the old command. The completed record also notes
that row growth is not exercised by this fixture; that limitation must remain explicit.

Use current supported entrypoints, record source revision and commands, and reconcile existing
results before deciding whether another run is needed. Close the applicable checks on observed
evidence or an explicit maintainer disposition.
