---
name: simple-frac-update-trailing
description: Use when studying, explaining, or modifying SIMPLE's fractional-update and trailing-reconstruction workflow, including update_frac sampling, sampled/updatecnt bookkeeping, matcher reconstruction I/O, partial even/odd reconstruction handoff, rec3D/volassemble ownership, downsampling compatibility, obsfield mirrors, and stage-to-stage previous-reference contracts.
---

# SIMPLE Fractional Update and Trailing Reconstruction

Use this skill for `abinitio3D`, `refine3D`, `update_frac`, `trail_rec`,
downsampling changes, previous-reference compatibility, and obsfield mirrors.

## Core Contract

- Treat the current fractional-update path as the reference behavior. Particle
  selection, partial reconstruction production, trailing reconstruction, FSC
  estimation, and next-stage handoff are separate responsibilities.
- Preserve the online reconstruction I/O contract. The matcher should read each
  active particle batch once, keep raw images needed for reconstruction, and
  produce partial even/odd reconstruction updates from those batch images after
  assignment.
- Do not solve producer/consumer mismatches at the reader by accepting
  incompatible previous artifacts. The producing stage must write the artifact
  expected by the consuming stage.
- Do not use volumes as a hidden source of truth when mirroring this for
  volume-free paths. Mirror the contract: previous artifact, current fractional
  update, update fraction, and compatibility handling.
- After command-line parsing, consume typed `parameters` fields rather than
  ad hoc `cmdline` lookups.

## Reading Order

Read `references/frac-update-contract.md` first.

Then inspect the current code in this order:

1. `src/main/simple_abinitio_controller.f90`
2. `src/main/strategies/search/simple_matcher_smpl_and_lplims.f90`
3. `src/main/ori/simple_oris_sampling.f90`
4. `src/main/ori/simple_oris_getters.f90`
5. `src/main/strategies/search/simple_strategy3D_matcher.f90`
6. `src/main/strategies/parallelization/simple_refine3D_strategy.f90`
7. `src/main/commanders/simple/simple_commanders_rec_distr.f90`
8. `src/main/volume/simple_reconstructor_eo.f90`

## Working Rules

- `volassemble` consumes the realized fractional update; it does not choose the subset.
- Do not move current partial reconstruction into a separate full particle pass
  that re-reads stacks to reduce peak memory; that is a workflow-policy change.
- The selected subset is represented by `sampled`; persistent coverage is
  represented by `updatecnt`.
- For `trail_rec=yes`, the previous even/odd artifact is mandatory. The current
  update is blended with that previous artifact using the realized update
  fraction, or `ufrac_trec` when explicitly provided.
- Downsampling compatibility is handled by the previous-artifact reader/producer
  contract. `reconstructor_eo%read_eos_parallel_io` pads previous smaller
  halfmaps/rhos when `l_update_frac` is active and rejects previous larger
  dimensions.
- When a stage boundary changes the next consumer's representation size, the
  prior stage must write the previous artifact in the next consumer's
  representation, not merely in its own search representation.
- Do not confuse abinitio3D's planned stage LP with the saved `_stageNN_lp.mrc`
  snapshot cutoff. Stage snapshots are filtered to the measured state FSC
  resolution when available; the planned stage LP is only the fallback.
- For obsfield work, use obsfield as the current-update accumulation or
  intermediate representation. Extract an unregularized reprojection model
  before FSC calculation, trailing blend, or ML prior application.
