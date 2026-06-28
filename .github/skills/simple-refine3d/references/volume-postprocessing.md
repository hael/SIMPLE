# Volume Postprocessing And Volassemble

## Scope

Use this note when the task involves `volassemble`, automasking, FSC-mask
consumption, or nonuniform filtering.

## Current Policy Boundary

Explicit assembly commanders are the execution point for assembled-reference
work used by refine3D execution modes.

That includes:

- assembly of Cartesian partial reconstructions or polar partial sums
- even/odd handling
- gridding correction
- merged-volume creation
- automask generation
- nonuniform filtering
- projection or reduction of polar reference handoff files

This is deliberate. The expensive shared-memory work is centralized so
shared-memory and distributed refine3D stay behaviorally aligned.
`simple_refine3D_strategy.f90` dispatches assembly but should not absorb the
low-level postprocessing or polar-reference reduction logic.

## Main Code Paths

- `src/main/commanders/simple/simple_commanders_rec_distr.f90`
  Executes Cartesian and polar assembly.
- `src/main/volume/simple_vol_pproc_policy.f90`
  Computes the postprocessing plan for one state and iteration.
- `src/main/volume/simple_reconstructor_eo.f90`
  Consumes compatible state masks when FSC is computed.
- `src/main/strategies/search/simple_matcher_pftc_prep.f90`
  Checks and loads `POLAR_REFS*` for matcher/probability-table consumers.
- `src/main/strategies/search/simple_matcher_refvol_utils.f90`
  Materializes polar references from complete volume sources when a file handoff
  is required; callers settle the 3D k-range with `set_bp_range3D` before
  entering this path.

## Important Current Behavior

- `plan_state_postprocess` decides whether automask runs and whether the state
  mask should feed nonuniform filtering.
- State masks are named `automask3D_stateNN.mrc`.
- Compatibility is box plus sampling equality, not just dimensions.
- Nonuniform filtering prefers the compatible state mask; otherwise it falls
  back to a spherical mask from `mskdiam`.
- `automsk=tight` is a real policy value and should not be silently collapsed away.
- Cartesian assembly always refreshes polar central sections for subsequent 3D
  matching, even when the scientific output is Cartesian volumes.
- Polar assembly writes state-major `POLAR_REFS.bin`, `POLAR_REFS_even.bin`,
  and `POLAR_REFS_odd.bin` for `polar=yes|obsfield`.
- `abinitio3D` stage `_lp` snapshots are FSC-resolution diagnostics. Do not
  replace their cutoff with the planned stage LP except as the existing fallback
  when no valid FSC resolution is available.

## Architecture Guardrails

- Keep Cartesian and polar assembly commanders as the executors of assembled-reference work.
- Push branching policy into helpers or policy objects when possible.
- Keep FSC-mask consumption logic consistent with mask-production compatibility checks.
- Preserve state-specific artifacts and naming conventions.
- Preserve the source-vs-handoff `POLAR_REFS*` producer/consumer contract.

## Common Mistakes

- Reintroducing postprocessing logic back into refine3D strategy code.
- Treating assembly commanders as distributed-only helpers or plain I/O combiners.
- Using a mask compatibility rule in one place and a different one somewhere else.
- Forgetting that nonuniform filtering and FSC may consume the same state mask for different purposes.
- Having probability-table code reproject Cartesian volumes instead of consuming materialized polar reference files.
