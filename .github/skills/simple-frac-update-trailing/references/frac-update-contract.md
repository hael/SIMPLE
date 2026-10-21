# Fractional Update Contract

## Stage Policy

`simple_solve3D_controller.f90` (a submodule of `simple_solve3D_utils`) emits
the per-stage `refine3D` command line:

- `update_frac` is the run's target outer particle subset, fixed for the run:
  `nsample * nstates / active`, capped at `UPDATE_FRAC_MAX` (0.9). When
  `nsample / active` exceeds `FULL_SAMPLE_SWITCH_FRAC` (0.9) every stage is a
  full update: no `update_frac`, `nsample`, `fillin` or trailing on the stage
  lines.
- `trail_rec` turns on the trailing consumer path from stage 5 in single-state
  runs (`TRAILREC_STAGE_SINGLE`) and stage 8 in multi-state runs
  (`TRAILREC_STAGE_MULTI`).
- `fillin=yes` is set on the last active stage of a single-state run only.
- `nspace` and stage search policy are search/reference concerns; they must not
  silently desynchronize the produced handoff artifact from what the next stage
  will consume.

## Particle Selection

`simple_matcher_smpl_and_lplims.f90` dispatches the outer subset:

- `sample_ptcls4iter3D` for a search iteration: `update_missing` runs
  `sample4update_missing`; with `fillin=yes` every fifth iteration
  (`mod(which_iter,5) == 0`) runs `sample4update_fillin`; otherwise
  `sample_ptcls4update3D`.
- `sample_ptcls4update3D`: `sample4update_all` when fractional update is
  inactive; `sample4update_rescore` for the cohort schedule
  (`cohort_sampling=yes`, every iteration of a block after its first);
  `sample4update_class` for the nested equal quota over the sampling units of
  `clssmp.bin` (`balance=class|cavg`); `sample4update_cnt` for
  low-update-count biased sampling (`balance=none`).
- `master_draws_sample3D`: a distributed non-probabilistic fractional pass is
  drawn once by the master over the whole project and written; workers
  reproduce it with `sample4update_reprod`.

`simple_prob_align_strategy.f90` owns the probabilistic pre-step:

- `prob_align_exec` / `prob_align_neigh_exec` sample once (clearing `sampled`
  and `updatecnt` only at `startit == 1` on its first iteration) and write the
  subset to the project.
- `prob_tab`, `prob_tab_neigh` (`simple_commanders_prob.f90`) and
  `refine3D_exec` reproduce that exact subset with `sample4update_reprod`, as
  does the `pose_cont=yes` polish pass.

`simple_oris_sampling.f90` and `simple_oris_getters.f90` define the handshake:

- `sampled` marks the current sampling round.
- `updatecnt` records persistent update history.
- `get_group_update_counts('state', ...)` returns `N(s)` (active rows with
  `updatecnt > 0`) and `n(s)` (those in the latest round);
  `get_state_update_fracs` is `n/N` per state; `get_update_frac` is the global
  summary. Under `m_estimator=flex` the same quantities are applied weight
  masses (`get_update_rows` plus the state weight set).

## Current Partial Artifact

`simple_strategy3D_matcher.f90` writes partition-local current updates after
the assignment barrier (`write_partial_recs`): per state and half, partial
even/odd sums and rho files (`simple_matcher_3Drec.f90`, one standalone
half-map reconstructor at a time), or the PCG worker's raw accumulators.

This is not a final reference yet. It is the current fractional update,
partitioned for assembly.

The reconstruction takes its images from the pass buffer the search filled
(`fill_pass_buffer`), or, with `box_crop < box`, from the cropped records the
search batches made (`take_records`, `get_rec_observation`). Over budget the
buffer stays empty and the reconstruction reads the stacks once. Do not add a
second full particle pass to lower peak memory unless the performance policy
explicitly changes.

## Assembly and Trailing Reconstruction

`simple_volassemble_strategy.f90::volassemble_exec` owns gridding volume
assembly (`commander_volassemble` is a thin wrapper; the distributed master
calls the routine on its own builder):

- It reduces partial reconstructions (`restore_state_from_parts`).
- It always computes realized per-state fractions `f`
  (`determine_trailing_update_fraction`); the applied map-update weight `u`
  equals `f` unless a single-state `ufrac_trec` override is parsed
  (`l_ufrac_trec_defined`).
- If `trail_rec=yes`, trailing happens in the accumulator domain
  (`blend_trailing_accumulators`), mirroring 2D class-average restoration. The
  previous artifact is the per-state chain `trailrec_stateNN_{even,odd}.mrc` +
  `rho_*` files + `trailrec_stateNN.txt` manifest: blended, unregularized e/o
  Fourier sums and sampling densities at the mass `M` the manifest records.
  With `N`, `n` and `f = n/N` from the merged project, the population rule
  (`population_blend_weights`, `simple_oris`) gives:

```text
s = u / f                 (current partials; 1 by default)
w = (1 - u) * N / M       (chain; 1 - u when the population is unchanged)
trailed_sums_eo = s * current_partial_sums_eo + w * chain_sums_eo
trailed_rho_eo  = s * current_partial_rho_eo  + w * chain_rho_eo
M_new = s*n + w*M = N
```

  The restored current-map coefficient is exactly `u` (the historical
  `ufrac_trec` meaning) and the chain keeps the mass of the represented
  population: first-time particles add mass, deactivated rows take their
  share with them.
- The blended accumulators are written back as the new chain BEFORE
  restoration, because ML regularization adds 1/tau^2 to rho in place.
- A single sampling-density correction after the blend restores the trailed
  halves, so each Fourier component is weighted by its accumulated sampling
  density. The FSC is estimated post-blend from the restored trailed halves and
  describes the artifact written to disk; the merged volume and nonuniform-filter
  inputs are derived from the same blended statistics.
- Chain start (chain absent, blend due): no blend. The chain is seeded with the
  current partials scaled by `1/f`, so the stored chain carries full-dataset
  mass (`M = N`), and the current partials are scaled back by `f` and restored,
  so this iteration's halves, FSC, merged volume and nonuniform-filter inputs
  are those of the current sample alone. Finished half maps are never blended
  and the previous volumes are not read. Without the `1/f` normalization a
  chain seeded at fractional mass `f` would give the next iteration an
  effective update weight of `f / (f + (1-f)*f)` — for `f=0.1`, a requested
  10 percent update would act like ~53 percent. From the next iteration every
  blend uses the chain.
- No sample (realized fraction below 0.001 with at least one sampled particle):
  a valid chain carries unchanged and the iteration is restored from it;
  without a chain the state is carried forward (previous volume kept, nothing
  written) when the directory holds its previous volume, half maps and FSC
  file; with neither, the run stops with an error (a map of the sample alone
  would be degenerate). In a multi-state run a state with no partials at all is
  dropped before any trailing code runs (`determine_dropped_states`).
- Full-weight seeding: a stage-boundary full reconstruction writes the chain at
  full-dataset weight when the internal `trail_seed=yes` cline handshake is set
  (`simple_solve3D_utils.f90::calc_rec`, gated on the consuming stage's
  `trail_rec`); the seed records the rows `sample4rec` reconstructs.
- Artifact-set integrity: the manifest is deleted before and rewritten after
  the four data files, recording per-component byte sizes, a generation
  counter, the format version, `M`, the contributor count, the state weight
  set identity, and provenance (box, sampling, particle population, state
  layout). `validate_trail_chain` accepts a chain only when the manifest
  parses, the provenance and weight-set identity match, every component size
  matches, and the grid is not larger than the current one with the same
  physical extent; any failure discards the complete set and re-seeds.
- Continued runs from another directory carry chains over as complete sets
  only, destination cleaned first and manifest copied last
  (`simple_refine3D_strategy.f90::carry_over_trail_rec_chains`).
- The PCG backend's master (`simple_rec3D_pcg_strategy.f90`) follows the same
  rules on its raw accumulators; its halves share a `pcg_trail_stateNN.txt`
  manifest.
- The recurrence, population rule, override weighting, chain-start
  normalization and the chain-start iteration's current-sample-only map are
  covered by the deterministic `trailing-reconstruction blend` sub-suite of the
  fast gate (`simple_test_exec test=unit_image suite=trailing_reconstruction_blend`).

The chain files become the previous artifact for the next iteration or stage.
Their names deliberately avoid the `recvol_state` stem so partial-reconstruction
globs and cleanup never match them.

## Downsampling Compatibility

`read_gridding_pair_accumulators`, contained in `restore_state_from_parts` in
`simple_volassemble_strategy.f90`, is the previous-artifact compatibility
mechanism, used for both partial reconstructions and the trailing
accumulator chain:

- It requires previous MRC and rho files for both even and odd halves.
- When `l_update_frac` is active and previous dimensions are smaller than current
  dimensions, it pads with zeros (the solve3D autoscale ramp).
- It rejects previous dimensions larger than the current dimensions. For the
  trailing chain, `validate_trail_chain` checks the grid and physical extent
  before the read and discards + re-seeds a larger chain instead of failing the
  run.

This is a producer/reader contract for previous artifacts, not a license to
invent a new fallback representation.
