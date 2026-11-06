# solve3D gold-standard matching

## Contract

`gold_standard=yes` is opt-in (default `no`) for single-state particle
solve3D. Every executed stage numbered at least 6 uses half-specific
references: NU modes become plain `nonuniform`, scheduled `lp` is deleted,
and only explicit user `lpstop` and crop Nyquist cap matching. Explicit `lp`
is refused with this option. Non-NU gold-standard stages use the existing
independent-half FSC bandwidth policy.

Preserve FSC=0.5 boundary advancement before the switch, with no new cap at
the planned stage limit. Preserve stage skipping, fractional updates,
trailing accumulators, early stopping and the final reconstruction. Add no
full alignment at the transition. Multi-state and class-average routes stay
on their existing matching topology. All solve3D stages and the ending use
spherical support (`automsk=no`, `envfsc=no`). No combine_eo or solvent-prior
mechanism is added.

## Implementation and review

- Add and bind the public yes/no parameter and its solve3D/chunks UI input.
- The stage configuration owns the guard; registration explicitly clears it.
- Existing reference loading and both reconstruction backends consume the
  emitted policy. No matcher or NU numerical algorithm changes are needed.
- Remove envelope controls and their unused staging state from this workflow.
- Record/replay the option through existing manifest input records so add-ons
  preserve it; there is no new manifest record format.
- Regression assertions exercise stage 5, a direct 5-to-7 transition, stage 8,
  both backends, explicit ceilings, opt-out, state/class-average guards,
  fractional/trailing preservation, spherical support, parsing and replay.

The switch does not erase correlations in frequencies already aligned with
shared references. Earlier FSCs retain that interpretation. No reference
bootstrap or trailing-chain reset is introduced.

## Validation

Source review and `git diff --check` passed. At the maintainer's request,
no build or test is run by the agent. Before scientific use, the
maintainer should build and run the fast gate, `unit_core suite=parameters`,
`unit_reconstruction suite=solve3D_cavgs_stage_rules`, and
`unit_project suite=solve3D_manifest`, followed by representative single-state
shared/distributed and gridding/PCG workflows with the option enabled.
No Linux/BOX or map-quality result is claimed.
