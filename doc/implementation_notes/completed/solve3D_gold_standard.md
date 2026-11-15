# Gold-standard stage policy for solve2D and solve3D

## Contract

`gold_standard=yes|no` is removed. `gold_std_stage=0` selects the workflow's
default: stage 2 in solve2D, stage 6 in single-state particle solve3D. A
positive integer overrides the start; a value beyond the executed schedule
defers the switch past this run. Negative values are rejected. Activation
depends on the stage number, including skipped stages and seeded restarts,
independently of ML regularization and the objective function.

In solve2D the controller carries the resolved start and current stage to
shared-memory and distributed refine2D workers through existing parameter
bindings. Before the switch they use merged references; from the switch they
use the particle's half-specific reference and the corresponding FRC filtering
convention. Scheduled low-pass limits and the existing low-population merged
fallback are retained. Seed and terminal passes follow their owning stage.
Parameter derivation sets `l_gold_standard_resolved` when the stage and resolved
threshold are present. Both 2D matcher consumers use that bit and
`l_gold_standard`; direct refine2D callers without a stage, including stream
and pool code, retain their existing `l_lpset` behavior.

In solve3D the stage configuration owns activation. NU modes become plain
`nonuniform`, scheduled or explicit `lp` is deleted, and only explicit user
`lpstop` and crop Nyquist cap matching. Non-NU gold-standard stages use the
existing independent-half FSC bandwidth policy. An explicit `lp` applies
before the switch. PCG backend and order-prior start stages remain at 6,
independently of the gold-standard override.
Both PCG transitions use the single `PCG_START_STAGE` constant.

Preserve solve3D FSC=0.5 boundary advancement before the switch, with no new cap at
the planned stage limit. Preserve stage skipping, fractional updates,
trailing accumulators, early stopping and the final reconstruction. Add no
full alignment at the transition. Multi-state and class-average routes stay
on their existing matching topology. All solve3D stages and the ending use
spherical support (`automsk=no`, `envfsc=no`). No combine_eo or solvent-prior
mechanism is added.

## Implementation plan and review

- Replace the boolean field, parser binding and public UI with the integer
  override; validate its domain and remove the former CLI key entirely.
- Add stage-policy routines in both controllers. Reuse the existing `stage`
  worker key in solve2D; introduce no internal CLI switch. Derive the worker's
  reference choice once from its current stage and resolved threshold.
- Keep the 2D matching band independent of reference topology. A deferred
  switch also keeps merged references during an FRC-driven final stage.
- Keep solve3D's existing reference loading, NU algorithm and reconstruction
  backends. Registration explicitly clears activation; multi-state and
  class-average routes keep their existing topology. Explicit positive 3D
  overrides require one state.
- Replace the manifest's recorded/replayed key. Bump its version to 3:
  default-on activation changes the interpretation of an absent old boolean,
  so replaying a version 2 manifest would silently change frozen policy.
- Migrate existing tester references to the removed field/key so they still
  compile. Add and run no tests, as requested by the maintainer.

The switch does not erase correlations in frequencies already aligned with
shared references. Earlier FSCs retain that interpretation. No reference
bootstrap or trailing-chain reset is introduced.

## Validation criteria and status

Review the CLI chain, default resolution, explicit earlier/later stages,
seeded and skipped-stage entry, stale command-line cleanup, and frozen-policy
replay. Build the existing production executables and inspect their compiler
and UI-metadata output. Do not run the fast gate or runtime tests.

Source review and `git diff --check` passed. Both `simple_exec` and
`simple_private_exec` built successfully with GNU Fortran 16 in the existing
Release build, with no compiler warnings or errors. The generated accepted-key
list contains `gold_std_stage` and omits the retired boolean and derived flag;
the generated display default is 0. Repository navigation was regenerated.
Latest build log: `/private/tmp/simple-gold-policy-resolved-build.log`.

No tests were added or run, as requested. Runtime, Linux/BOX and map-quality
validation remain unperformed. Changes remain uncommitted.
