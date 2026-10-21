---
name: simple-frac-update-trailing
description: Use when studying, explaining, or modifying SIMPLE's fractional-update and trailing-reconstruction workflow, including update_frac sampling, sampled/updatecnt bookkeeping, the particle buffer and cropped records feeding partial reconstruction, the per-state trailing accumulator chain and its population blend in volassemble_exec and the PCG master, downsampling compatibility, and stage-to-stage previous-artifact contracts.
---

# SIMPLE Fractional Update and Trailing Reconstruction

Use this skill for `solve3D`, `refine3D`, `update_frac`, `trail_rec`,
downsampling changes, and previous-artifact compatibility.

## Core Contract

- Treat the current fractional-update path as the reference behavior. Particle
  selection, partial reconstruction production, trailing reconstruction, FSC
  estimation, and next-stage handoff are separate responsibilities.
- Preserve the reconstruction read contract. A 3D pass reads its share once
  into the pass buffer, and the reconstruction after the assignment barrier
  copies its images from it; with `box_crop < box` the search batches keep
  cropped observations as records and both backends take them
  (`simple_matcher_ptcl_io.f90`, `doc/policies/particle_read_path_policy.md`).
- Do not solve producer/consumer mismatches at the reader by accepting
  incompatible previous artifacts. The producing stage must write the artifact
  expected by the consuming stage.
- After command-line parsing, consume typed `parameters` fields rather than
  ad hoc `cmdline` lookups.

## Reading Order

Read `references/frac-update-contract.md` first, then
`doc/policies/importance_sampling_fractional_update_policy.md` section 7.

Then inspect the current code in this order:

1. `src/main/solve/simple_solve3D_controller.f90`
2. `src/main/strategies/search/simple_matcher_smpl_and_lplims.f90`
3. `src/main/ori/simple_oris_sampling.f90`
4. `src/main/ori/simple_oris_getters.f90`
5. `src/main/strategies/search/simple_strategy3D_matcher.f90` and
   `simple_matcher_3Drec.f90`
6. `src/main/strategies/parallelization/simple_refine3D_strategy.f90`
7. `src/main/strategies/parallelization/simple_volassemble_strategy.f90`
   (`volassemble_exec`; `restore_state_from_parts` contains the trailing
   blend, the chain validation and the gridding pair accumulator I/O)
8. `src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90`
   (the PCG master's chain)
9. `src/main/volume/simple_halfmap_diagnostics.f90`

## Working Rules

- Volume assembly consumes the realized fractional update; it does not choose
  the subset.
- The sample, `updatecnt`, poses, states and partial reconstructions of an
  assembly belong to one iteration. Shared-memory gridding assembly
  (`run_volassemble`) reads the project file, so the strategy writes the
  particle field before it (`publish_before_file_assembly`). The distributed
  master runs `volassemble_exec` on its own in-memory particle field after the
  in-memory merge (`write=.false.`); do not add a project reread there.
- Do not move partial reconstruction into a separate pass that re-reads stacks
  to reduce peak memory; that is a workflow-policy change.
- The selected subset is represented by `sampled`; persistent coverage is
  represented by `updatecnt`. Under `cohort_sampling=yes` (`refine3D_states`
  frequency blocks) the same subset is rescored at every iteration of a block,
  so trailing reconstruction gives that cohort the cumulative current-map
  coefficient `1 - (1-u)^k` after `k` iterations; no per-unit weight enters the
  state-level blend.
- For `trail_rec=yes`, the previous artifact is the per-state accumulator chain
  (`trailrec_stateNN_{even,odd}` + rho files + `trailrec_stateNN.txt`
  manifest): blended, unregularized e/o Fourier sums and sampling densities at
  the mass `M` of the population they represent. With `N`/`n` the state's
  updated and sampled populations (`get_group_update_counts`), `f = n/N` the
  realized fraction and `u` the applied weight (`ufrac_trec` for a single
  state when given, else `f`), the population rule
  (`population_blend_weights`, `simple_oris`) scales the current partials by
  `u/f` and the chain by `(1-u)*N/M` (`1-u` when the population is unchanged);
  the chain then records `M = N`. A single restoration after the blend
  produces the trailed halves with current-map coefficient `u`; the FSC is
  estimated post-blend.
- The chain is written before restoration, because regularization mutates rho
  in place. Its filenames deliberately avoid the `recvol_state` stem so
  partial-reconstruction globs and cleanup never match it.
- Every blend happens in the accumulator domain; restored (finished) half maps
  are never blended. When a blend is due but no chain exists, the iteration
  does not blend: the chain is seeded with the current partials scaled by `1/f`
  so it carries full-dataset mass (a fractional-mass seed would make the next
  update far more aggressive than requested), and the iteration restores and
  ships the current sample alone (partials scaled back by `f`), FSC included.
  The previous volumes are neither read nor needed, and no extra reconstruction
  runs. A state whose realized fraction is below 0.001 carries a valid chain
  unchanged and is restored from it; without a chain it is carried forward like
  a dropped state when the directory holds its previous volume, half maps and
  FSC; with neither a chain nor a previous volume the run stops with an error
  (raise `update_frac`). Both backends follow these rules. Stage-boundary full
  reconstructions seed the chain at full-dataset weight via the internal
  `trail_seed` cline handshake, gated on the consuming stage's `trail_rec` so
  early boundaries cannot park stale alignments in the chain.
- The five chain files are one validated artifact set: the manifest is deleted
  first and rewritten last with per-component byte sizes, generation, format
  version, represented population `M`, contributor count, state weight set
  identity and provenance (box, sampling, particle population, state layout).
  `validate_trail_chain` discards the complete set on any mismatch and the
  iteration re-seeds; carry-over between directories copies complete sets
  only, manifest last (`carry_over_trail_rec_chains`).
- Downsampling compatibility is handled by the previous-artifact reader/producer
  contract. `read_gridding_pair_accumulators` pads previous smaller
  halfmaps/rhos when `l_update_frac` is active and rejects previous larger
  dimensions. For the trailing chain, a larger previous grid or a physical
  extent mismatch discards the set and re-seeds instead of failing the run.
- The distributed PCG master applies the same population weights to its raw
  accumulators; its chain halves share the `pcg_trail_stateNN.txt` manifest.
- The recurrence, population rule, `ufrac_trec` weighting and chain-start
  contracts are covered by the deterministic `trailing-reconstruction blend`
  sub-suite of the fast gate
  (`simple_test_exec test=unit_image suite=trailing_reconstruction_blend`).
- When a stage boundary changes the next consumer's representation size, the
  prior stage must write the previous artifact in the next consumer's
  representation, not merely in its own search representation.
- Do not confuse solve3D's planned stage LP with the saved `_stageNN_lp.mrc`
  snapshot cutoff. Stage snapshots are filtered to the measured state FSC
  resolution when available; the planned stage LP is only the fallback.
