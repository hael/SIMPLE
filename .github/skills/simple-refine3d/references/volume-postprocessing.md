# Volume Postprocessing And Volume Assembly

## Scope

Use this note when the task involves volume assembly (`volassemble_exec`),
automasking, FSC-mask consumption, or nonuniform filtering.

## Current Policy Boundary

Volume assembly is the execution point for assembled-volume work in every
refine3D execution mode.

That includes:

- reduction of partition-local Cartesian partial reconstructions
- trailing reconstruction in the accumulator domain
- even/odd handling and restoration (gridding correction, sampling-density
  correction, ML regularization)
- FSC, cFAR and resolution fields on the shipped halves
- merged-volume creation
- density automask generation (`envfsc=yes`)
- nonuniform filtering (`simple_nu_state_filter`) and the NU matching
  low-pass handoff

This is deliberate. The expensive shared-memory work is centralized so
shared-memory and distributed refine3D stay behaviorally aligned.
`simple_refine3D_strategy.f90` dispatches assembly and later postprocessing but
does not absorb restoration or filtering logic. Volume assembly does not make
matcher references: the strategy projects the next reprojection model from the
volumes it writes.

## Main Code Paths

- `src/main/strategies/parallelization/simple_volassemble_strategy.f90`
  `volassemble_exec` on a caller's `(params, build)` (the distributed master's
  own, with `write=.false.`), `run_volassemble` for command-line callers
  (`commander_volassemble`), `restore_state_from_parts` for one state, and
  `filter_pcg_nonuniform_maps` for the PCG master's NU handoff.
- `src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90`
  The PCG backend's master (`execute_rec3D_pcg_distributed_master`), called
  through `assemble_refine3D_pcg` in the refine3D strategy.
- `src/main/volume/simple_nu_state_filter.f90`
  The shared NU competition of both backends.
- `src/main/volume/simple_vol_pproc_policy.f90`
  `state_mask_is_compatible`: box and sampling equality of a state mask.
- `src/main/volume/simple_halfmap_diagnostics.f90`
  Consumes compatible state masks when FSC is computed.
- `src/main/strategies/search/simple_matcher_refvol_utils.f90`
  Reference preparation for the next iteration's reprojection model.

## Important Current Behavior

- State masks are named `automask3D_stateNN.mrc`.
- Compatibility is box plus sampling equality, not just dimensions.
- Nonuniform filtering always uses spherical support from `mskdiam`, constructed
  internally by `setup_nu_dmats`.
- Density automasks remain independent inputs to `envfsc`; they do not
  restrict the NU objective domain and are never multiplied into a matching
  reference.
- `automsk=tight` is a real policy value (Otsu tightness) and should not be
  silently collapsed away; NU refinement rejects it.
- The NU source halves are the restored halves kept in memory as written.
- Every gridding assembly writes the `_unfil` half pair beside its map.
- A state without partial reconstructions in a multi-state run is dropped and
  its previous volume, halves and FSC are carried forward unchanged.
- Assembly records `res`, `res05`, `cfar`, `lp` (NU handoff) and `lp_est` in
  the particle field; on the distributed master they stay in memory for the
  matching low-pass refresh, convergence and sigma commit.
- The distributed master runs per-state `postprocess` after assembly and
  registers state volumes and FSCs in `os_out`.
- `solve3D` stage `_lp` snapshots are FSC-resolution diagnostics. Do not
  replace their cutoff with the planned stage LP except as the existing fallback
  when no valid FSC resolution is available.

## Architecture Guardrails

- Keep `volassemble_exec` (and the PCG master) as the executors of
  assembled-volume work; keep `commander_volassemble` a thin wrapper.
- Route changes to ML ordering, sampling-density correction, FSC handling,
  trailing reconstruction or LP-set behavior through `restore_state_from_parts`
  rather than duplicated branches.
- Keep FSC-mask consumption logic consistent with mask-production compatibility checks.
- Preserve state-specific artifacts and naming conventions.

## Common Mistakes

- Reintroducing postprocessing logic into refine3D strategy code.
- Making volume assembly write or reread the project on the distributed
  master's path.
- Using a mask compatibility rule in one place and a different one somewhere else.
- Reusing a density automask or future NU-evidence envelope as NU support.
- Sending a correlation-derived NU-evidence envelope into FSC correction.
- Having assembly or table code produce matcher references instead of the
  strategy's reprojection model.
