# Refine3D Code Map

## Scope

Use this note to find the right SIMPLE module before making changes.

## Top-Level Flow

For normal 3D refinement, the high-level call chain is:

1. `src/main/commanders/simple/simple_commanders_refine3D.f90`
2. `src/main/strategies/parallelization/simple_refine3D_strategy.f90`
3. `src/main/strategies/search/simple_strategy3D_matcher.f90` (under `refine=prob*` with
   the tables of `src/main/strategies/search/probabilistic/`)
4. `src/main/strategies/parallelization/simple_volassemble_strategy.f90`
5. `src/main/strategies/parallelization/simple_strategy_sigma2.f90`

Read them in that order unless the task is tightly scoped.

For public workflow policy, read `doc/policies/3D/refine3D_policy.md` first.

## Ownership By Layer

### Commanders

- `simple_commanders_refine3D.f90`
  Entry point, defaults, automatic workflow, and execution-strategy selection
  through `create_refine3D_strategy` (`nparts` without `part`: distributed
  master; otherwise shared memory).
- `simple_commanders_rec_distr.f90`
  `commander_volassemble`: a thin wrapper that calls `run_volassemble`.
- `simple_commanders_euclid.f90`
  `commander_commit_sigma2`: a thin wrapper that calls `run_commit_sigma2`.

### Strategy Layer

- `simple_refine3D_strategy.f90`
  Shared-memory versus distributed execution, iteration loop, scheduler
  interaction, per-iteration reprojection model
  (`materialize_reprojection_model`, `materialize_reprojection_model_master`),
  the iteration's sample draw (`sample_ptcls4prob3D` under `refine=prob*`,
  `master_draws_sample3D` otherwise), the run seed (`adopt_run_seed`), the
  merge of the parts' readout records (`merge_prob_readouts`), matcher dispatch,
  the `pose_cont=yes` polish pass, assembly dispatch (gridding or
  `assemble_refine3D_pcg`), the matching low-pass refresh, convergence and the
  sigma2 commit.
- `simple_volassemble_strategy.f90`
  `volassemble_exec` (gridding assembly on a caller's `(params, build)`),
  `run_volassemble` (its command-line form) and `filter_pcg_nonuniform_maps`
  (the PCG master's NU low-pass handoff). `restore_state_from_parts` owns the
  per-state restoration: partial reduction, trailing chain, ML ordering,
  sampling-density correction, FSC, merged map, NU inputs.
- `simple_strategy_sigma2.f90`
  `prepare_canonical_sigma_update` and `commit_sigma2_exec`.
- `simple_rec3D_pcg_strategy.f90`
  PCG workers and the distributed PCG master
  (`execute_rec3D_pcg_distributed_master`).

### Matcher/Search Layer

- `simple_strategy3D_matcher.f90`
  `refine3D_exec`: reads the reprojection model (`read_reprojection_model`),
  reproduces or draws the particle sample (`sample_particles_for_update`),
  under `refine=prob*` opens the pipeline (`begin_pipeline`), fills the pass
  buffer or opens the record buffer, per batch fills and draws
  (`pipeline_batch`: `fill_tab_range`, `assign_batch`) and runs the search
  under a dynamic OpenMP schedule, closes the assignment after the batches
  (`close_assignment`; a part writes its readout record), writes orientations and sigma2 ranges,
  releases the model at the assignment barrier, then writes partial
  reconstructions (`write_partial_recs`: `calc_3Drec` or the PCG worker).
- `simple_matcher_refvol_utils.f90`
  Reference preparation and the reprojection model files
  (`materialize_reprojection_model_from_volumes`, `read_reprojection_model`,
  `adopt_reprojection_model_range`, `read_cart_refvols`).
- `simple_matcher_ptcl_io.f90`, `simple_matcher_ptcl_batch.f90`
  Read sessions, the pass buffer and the cropped-record buffer;
  `build_batch_particles3D`.
- `simple_matcher_smpl_and_lplims.f90`
  `set_bp_range3D`, `sample_ptcls4iter3D`, `sample_ptcls4update3D`,
  `sample_ptcls4prob3D`, `master_draws_sample3D`, `adopt_run_seed`.
- `search/probabilistic/simple_strategy3D_prob.f90`
  Consumer of the batch's assignment map (`assgn_map`, drawn by `assign_batch`
  just before the search) and of the extra poses (`topk > 1`).
- `simple_strategy3D_alloc.f90`, `simple_strategy3D_utils.f90`
  Search-space allocation and shared helpers.
- `simple_strategy3D_*.f90`
  Concrete searches: `greedy`, `greedy_inpl`, `greedy_sub` (`refine=neigh`),
  `shc` (greedy for unsearched particles and at `GREEDY_FREQ`), `eval`,
  `cont` (Cartesian); `refine=sigma` runs no search object.

### Probability Tables (`src/main/strategies/search/probabilistic/`)

- `simple_eul_prob_tab.f90`
  Dense table on the matcher's batches: `new_pipeline`, `set_model_from_field`,
  `fill_tab_range` (keeps the batch's rows), `open_assignment`, `assign_batch`
  (`assign_rows`), `close_assignment`; the weighted insertion's table
  (`insert_tab`).
- `simple_eul_prob_tab_neigh.f90`
  Neighbourhood table (`prob_neigh`): per-thread candidate buffers gathered
  into the batch's candidate store (`assign_batch` -> `batch_store_from_buffers`),
  then drawn (`assign_rows_neigh`).
- `simple_eul_prob_tab_utils.f90`
  The posterior primitives (`posterior_temperature`, `draw_posterior_row`,
  `posterior_row_stats`, `prob_beta_at`), the SHC strategies' samplers
  (`angle_sampling`, `sample_bounded_dists`, `sample_power_dists`), the row and
  candidate entry types (`prob_dense_entry`, `prob_refined_entry`,
  `prob_candidate`), the readout record (`prob_readout`,
  `merge_prob_readouts`), the keyed scope (`prob_rnd_scope`).

### Reconstruction/Postprocessing

- `src/main/simple_final_rec.f90`
  The shared ending of `solve3D`, `refine3D_auto`, `refine3D_states`, and
  `classify3D_refs`. `calc_final_rec` runs the final all-particle
  reconstruction at native sampling from the refinement command line it is
  given, reuses committed canonical sigmas or rebuilds them through
  `bootstrap_rec3D`, registers state maps, FSCs, and envfsc masks in the
  project, and writes the final products and orthogonal reprojections. No
  workflow carries its own copy of this sequence.
- `simple_matcher_3Drec.f90`
  Partial reconstruction written by the matcher: groups the sample by
  `(state, half)` (or by weights under `m_estimator=flex`), one standalone
  half-map reconstructor at a time.
- `src/main/volume/simple_vol_pproc_policy.f90`
  `state_mask_is_compatible` (box and sampling equality of a state mask).
- `src/main/volume/simple_nu_state_filter.f90`
  The NU competition shared by gridding assembly and the PCG master.
- `src/main/volume/simple_halfmap_diagnostics.f90`
  Backend-neutral half-map FSC/cFAR/resolution evaluator.

## Current Artifact Contracts

- Probabilistic passes write no table, assignment or insertion file: a
  distributed part writes `prob_readout_part<NN>.txt` (counts and n_eff sums;
  zero for an empty part), which the master merges in part order, reports and
  deletes (`merge_prob_readouts`).
- The reprojection model (`reprojection_model_{even,odd}.bin`) is written by
  the strategy before the matcher, which reads it; its
  header carries the band (`kfromto`). Cartesian passes read
  `cart_refvols_{even,odd}.bin`. `remove_ref_section_files` deletes both sets.
- Multi-state references are state-major: state `s` occupies
  `(s-1)*nspace+1 : s*nspace`; a dense row holds the active states' references only.
- Cartesian matcher work writes partition-local partial reconstructions;
  stale partials are removed before each matcher iteration. Volume assembly
  reduces them into state/even/odd volumes and does not touch matcher
  references.
- The pass buffer lives for one pass in one process: filled before the search
  batches, copied into the reconstruction's batches after the assignment
  barrier, released after reconstruction. Cropped passes keep records instead.
- Multi-state `solve3D` point-group activation is state-local at the
  symmetry-search stage: each current state map supplies its own axis search,
  and the resulting orientation transform applies only to that state.
- Every gridding assembly leaves a `recvol_stateNN_even_unfil.mrc` /
  `_odd_unfil.mrc` pair that belongs to the map it accompanies: the
  unregularized halves under ML regularization, copies of the halves
  otherwise. `postprocess_nu` consumes that pair, and postprocessing ignores a
  pair whose box differs from the map.

## Questions To Ask Before Editing

- Is this change about candidate generation, assignment, or final orientation stamping?
- Is this change iteration control or search logic?
- Does this affect both shared-memory and distributed refine3D?
- Does it run on the distributed master's `(params, build)`, and does it add a project read or write there?
- Does this change touch particle-domain work or assembled-volume work?
- Does a policy doc in `doc/policies/` need to change with the code?
- Does the reprojection model stay the one source of matching references for the matcher and its pipeline?
- Does reconstruction still take its images from the pass buffer or the records instead of a second stack read?
- If multi-state `solve3D` activates a point group, does the symmetry-axis search remain state-local?
