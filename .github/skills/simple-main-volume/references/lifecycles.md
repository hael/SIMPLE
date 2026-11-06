# Volume lifecycles and artifacts

Reference for `src/main/volume`; the rules are in the guide (`../SKILL.md`). Verify a call order
against the callers named here before relying on it.

## Gridding (`simple_reconstructor.f90`)

The accumulator is an `image` on `box_crop`/`smpd_crop` with a density `rho` (FFTW-aligned, the
packed Fourier shape) and, with `expand=.true.`, the expanded pair `cmat_exp`/`rho_exp` indexed by
native Fourier cells (`h` from `lims(1,1)-wdim`, `k`, `m` symmetric). The KB window is the constant
`kbinterpol(KBWINSZ, KBALPHA)`; the deapodisation `invenv1d` is the discrete transform of that
stencil (`kb_stencil_inv_envelope_1d`).

Worker, per state and half (`calc_3Drec` in `strategies/search/simple_matcher_3Drec.f90`):

1. `new_accumulator(params, spproj, expand=.true.)`.
2. Per batch: `prep_imgs4rec` makes the planes in parallel (`gen_fplane4rec` on the
   `OSMPL_PAD_FAC`-padded `box_croppd` heap, whitened by the particle's half when `objfun=euclid`;
   with `box_crop < box` the cropped observation comes from the pass's records, `take_records`).
3. `insert_plane_oversamp(pgrpsyms, ori, plane[, w])` per plane, serially; under `topk > 1` the
   extra poses re-shift the committed plane (`fplane_shift`) and insert it at their pose and weight.
   Batches are state and half homogeneous (`update_state_half_rec` stops otherwise).
4. `compress_exp`, `write_raw_accum(recvol_stateNN_partPP_{even,odd}.mrc, rho_...)`, `kill`. Both
   halves of a populated state are written even when one is empty.

Assembly, per state (`restore_state_from_parts` in `simple_volassemble_strategy.f90`):

1. `reset`, then each part read into scratch accumulators and `sum_reduce`d in part order.
2. Trailing blend on the raw sums (`apply_weight_sums`, chain read and added, chain set written)
   or the full-mass seed; frozen set added (`add_gridding_set`); carried set written.
3. Without `ml_reg` the merged sum is formed now (`sum_rec`); fractional weights floor the
   densities (`floor_rho_shellwise`).
4. `restore_gridding_pair`: `restore_base` -> `finalize_from_base` (deapodised copy) -> soft sphere
   at `msk_crop` -> write -> FSC on the written halves (`evaluate_halfmap_pair`; the deapodised masked
   halves, not the base, are the intended FSC input). With `ml_reg` the
   base is restored with `preserve_numerator=.true.`, written as `_unfil`, the FSC taken, then
   `add_invtausq2rho(fsc)` and `restore_final(preserve_numerator=.true.)` give the shipped halves,
   and the merged sum is formed after the ridge.
5. `sum_rec%restore_final`, soft sphere, write the merged state volume; NU filtering follows.

`add_invtausq2rho`: per shell SSNR from the FSC (clamped to 0.001-0.999), noise from the mean
density, ridge `1/(tau*tau2)` (`wiener_ridge`; 1000 rho for a shell with no power estimate), none
below `max(6, find(hp))`.

Expanded-lattice users outside assembly: `exp_samples` (a sample set's KB geometry, gathered from
any volume on the lattice), `project_fplane` (central sections from the current `cmat_exp`: callers
`expand_exp` once, one sample set per thread) and `insert_planes_multi` (K targets, shared,
diagonal or packed density) serve `flex/`.

## PCG (`simple_reconstructor_pcg.f90`, `simple_pcg_lattice.f90`)

The unknown lives on the native box (`box_crop`, `smpd_crop`); every Fourier operation runs on the
`OSMPL_PAD_FAC`-padded lattice (centre-pad in, centre-crop out). A plane is the full both-sign-`h`
disk of radius `box/2`. CTF at `smpd_crop`, shifts in `box_crop` pixels (`rec_obs_oris`). The
shipped map is `window*u`, the soft support applied once.

Worker (`execute_rec3D_pcg_worker`, `accumulate_pcg_state_half`), per state and half:
`new(box_crop, smpd_crop, PCG_LAMBDA)` -> `set_sym` -> `prep_particles(selection, use_ctf,
inv_sig2 of that half, weights)` -> `begin_accum` -> `accumulate_batch` per batch ->
`write_raw_accum(state, eo, part, nparts, nptcls, provenance, mass, wset_id)` -> `kill`. Workers
never call `end_accum`.

Master (`prepare_distributed_half_job`): `new` -> support (`set_window_sphere(msk_crop)`, or
`set_window_volume` with the lag-one density envelope or an explicit `pcg_mskfile`; never the
NU-evidence envelope) -> `begin_reduction` -> `add_raw_accum` for parts 1..nparts (or the carried
chain through `add_raw_accum_weighted`) -> trailing: `scale_raw_accum`, chain added, chain half
written -> `end_accum(.true.[, d_raw])` -> `set_op_mode(PCG_OP_KERNEL)` -> solve from zero. The
regularised pair: `set_ml_prior(fsc of the base pair)` -> `regularize_finalized(d_raw)` -> the
closed-form shrink (`shrink_by_ml_prior`). Solvent and atomic priors after `end_accum`. The two
halves are prepared serially and solved in `parallel sections` with nested teams.

Raw artifact (version 3): magic, state, half, part, nparts, nptcls, applied mass, weight-set id,
box, boxpd, padf, lattice, smpd, provenance (`pcg_raw_provenance`: pgrp, objfun, iteration, boxes,
samplings, msk, ctf); payload the `h >= 0` half. The chain's provenance leaves out the iteration
and the crop; a smaller previous chain is embedded at constant field of view.

## Artifacts written here

Names come from `src/defs/simple_refine3D_fnames.f90`; assembly decides when.

| artifact | writer | reader |
| --- | --- | --- |
| `recvol_stateNN_partPP_{even,odd}.mrc` + `rho_*` | `reconstructor%write_raw_accum` | `reduce_partials` |
| `pcg_raw_stateNN_{even,odd}_partPP.bin` | `reconstructor_pcg%write_raw_accum` | master `add_raw_accum` |
| `trailrec_stateNN_*` + `.txt`; `pcg_trail_stateNN_*` + `.txt` | assembly, `trail_chain_manifest` | assembly, refine3D_local staging |
| `frozen_stateNN_boxBBBB_*`, `carried_*` + manifests | `frozen_accum` | `add_gridding_set`, `import_gridding_set` |
| state volume, halves, `_unfil`, `fsc_stateNN.bin`, `RESOLUTION_STATENN` | assembly through `simple_halfmap_diagnostics` | matchers, postprocess, final rec |
| `<vol>_pcg_support.txt` | `write_support_provenance` | postprocess (mask, FSC weighting) |
| `<vol>.provenance` (`.pending` while writing) | `begin_/publish_reference_set` | refine3D_local |
| `_nu_filt`, `_nu_locres`, `automask3D_stateNN`, `nu_envmask3D_stateNN` | `nonuniform_filter_state`, diagnostics | matchers, refine3D_local |

None of these is registered in the project here: `add_vol2os_out`/`add_fsc2os_out` are called by
`simple_final_rec.f90`, `solve/simple_solve3D_utils.f90`, `solve/simple_state_assignment.f90` and
`simple_refine3D_workflow_controller.f90`.

## Other modules

- Wave model: `simple_wave_projector` (projection, first Born, Rytov, multislice on `pcg_lattice`,
  native sampling; its own parallel regions, stops inside another) and `simple_wave_multislice`
  (one thread's sweep). Contract: `doc/policies/3D/wave_model_policy.md`.
- PCG priors: `simple_pcg_solvent_sidecar` (`pcg_solvent=yes`: per-half weight, re-solve, the
  apply pair of NU) and `simple_pcg_atomic_prior` (`atm_prior=yes`, shared by both halves, owned by
  the caller).
- `simple_halfmap_enhancer` (local Wiener sharpening, `halfmap_enhance`), `simple_volpair_metrics`
  (two maps of one frame: FSC and band correlation), `simple_vol_pproc_policy` (mask artifact box and
  sampling check). The dose-series manifest of `refine_motion_model` is in `motion/`.
- Analysis with module state, one use at a time: `simple_volpft_symsrch`, `simple_volanalyzer`,
  `simple_symanalyzer`; also `simple_dock_vols`, `simple_volcluster`, `simple_volpft_corrcalc`,
  `simple_volinterp` (module `simple_simple_volinterp`).
- `simple_reconstructor_openmpoffload`: GPU gridding, compiled only with `USE_OPENMP_OFFLOAD`.
