---
name: simple-main-volume
description: Use when working in SIMPLE's src/main/volume subsystem, including the gridding and PCG reconstructors (accumulators, raw part artifacts, restoration, the PCG lattice and priors), half-map FSC/cFAR diagnostics and the support-provenance sidecar, trailing-chain manifests, the frozen accumulators of solve3D_addon, refine3D reference-set provenance, assembly-owned nonuniform filtering of a state, the wave-model operators, volume analysis, docking, clustering, symmetry search, half-map sharpening and volume pair metrics.
---

# `src/main/volume`: reconstructors, half-map diagnostics and volume analysis

Owns the reconstruction objects (gridding `reconstructor`, PCG `reconstructor_pcg` on `pcg_lattice`,
wave operators) and their raw artifacts, the half-map FSC and support sidecar, the chain manifest,
the solve3D_addon frozen sets, the reference-set record, one state's NU filtering
(`nonuniform_filter_state`, the only NU path of both backends: change it, not a caller) and analysis
tools. Which particles enter, blend weights, restore order and file names are assembly's
(`simple_volassemble_strategy.f90`, `simple_rec3D_pcg_strategy.f90`); planes are made in
`simple_matcher_3Drec.f90`. Nothing here writes the project. The NU algorithm is `nu_filt/`.

## Rules

Silent ones first: nothing reports these.

- Raw before restored: blend, scale, sum and publish (partials, chains, frozen sets) raw accumulators
  only. `restore_base`/`restore_final` divide cmat by rho and inverse-transform the object unless
  `preserve_numerator=.true.`; `add_invtausq2rho` adds the ML ridge to rho in place. What is written
  or summed after them is a map or carries the ridge. PCG `end_accum` consumes its accumulators.
- Gridding: `new_accumulator(expand=.true.)` -> `insert_plane_oversamp` -> `compress_exp` ->
  `write_raw_accum` | `sum_reduce` -> restore -> `kill`. `compress_exp` overwrites the native
  cmat/rho (a sum read in before is lost). The object points to `params` (`tau`, `hp`, `box_crop`
  read at `add_invtausq2rho`): keep that alive.
- `insert_plane_oversamp`, `insert_planes_multi` and PCG `accumulate_batch` open their own parallel
  regions: call them from serial code, one object per caller, or they race. No module state in
  `reconstructor_pcg`: the master solves both halves at once (nested sections).
- One accumulator per (state, half); it does not know its half. Planes are whitened by their half's
  sigma2; PCG `prep_particles` takes one spectrum: one half per call, or the others are whitened
  wrongly. CTF at `smpd_crop`, shifts in `box_crop` pixels, planes padded by `OSMPL_PAD_FAC`.
- A gridding density file (`rho_*`, partials, chains, frozen sets) starts with a header (magic, box,
  density shape, sampling; `RHO_HDR_BYTES`): `read_rho` and the readers of open units
  (`check_rho_header`, serially, before `read_raw_rho`) stop on another box or sampling, and on a file
  without the header; the `update_frac` padding of a smaller chain checks the header's box against its
  numerator and the physical extent. The numerator (`<fbody>_even/odd.mrc`) is still read in the reader's
  shape without its MRC header. PCG raw artifacts stop on any header mismatch. Frozen-set contexts of the
  schema before the header (version 3) are refused (`FROZEN_SCHEMA_VERSION` 4); trailing chains of the
  manifest before it are re-seeded (version 4).
- FSCs are per Fourier shell of their box (`box_crop`) and carry no box: `add_invtausq2rho(fsc)` and
  PCG `set_ml_prior(fsc)` take the overlapping shells of whatever they get.
- PCG setters are read once: `set_sym` by the batches accumulated after it, the support and
  `set_ml_prior` by `end_accum` (the support also by the solve). Late, the result is partial,
  inconsistent (operator masked, RHS not) or unregularised (until `regularize_finalized`). Solvent
  and atomic priors need `end_accum` (stop otherwise).
- A volume and its `_pcg_support.txt` sidecar are one artifact (`copy/rename/remove_support_provenance`):
  without it postprocess masks and FSC-weights again; a stale one skips both.
- `begin_reference_set` before a state's volumes, `publish_reference_set` after all its products. Its
  list is exactly what refine3D_local stages (`stage_reference_derivatives`,
  `simple_refine3D_strategy.f90`), no more: change both together, or refine3D_local stops carrying
  the product.
- Chains and frozen sets: manifest deleted first, written last; frozen sets are gridding only, added
  with weight one before any restoration (skill `simple-frac-update-trailing`, addon policy §5-8).
- Couplings: the expanded-lattice API serves `flex/`; flex PCG and the wave operators extend
  `pcg_lattice`; deapodisation is the transform of the operator's KB stencil (`interp/`).

## Where to start

Call orders (callers `calc_3Drec`, `restore_state_from_parts`, `accumulate_pcg_state_half`,
`prepare_distributed_half_job`), artifacts, other modules: [references/lifecycles.md](references/lifecycles.md).

## Tests

- `build/bin/simple_test_exec test=unit_reconstruction` (about 7 s; `suite=list`, e.g.
  `weighted_insertion`, `frozen_accumulator`); `test=unit_image suite=trailing_reconstruction_blend`.
- `test=pcg_recon` (14 gated stages): a change to the PCG kernel, scale, masks, symmetry, weights,
  preconditioner or box conversion adds a kernel-versus-matrix-free gate there.
- Nightly: `test=lib_reconstruction` (`pcg_half_set`, `pcg_sparse_sampling`,
  `fractional_reconstruction`, ...), `test=pcg_frac_update`; `flex/`: `test=unit_heterogeneity`.
- Behaviour: `simulated_workflow_6vxx`/`_1jxy` (a robustness profile: 1-2 misses in 20 runs is normal),
  `solve3D_addon`, `multistate_workflows suite=states`, `cont_refine3D_1jxy suite=nu_handoff`;
  `test=rec3D_backends` compares the backends on a real project.

## See also

`doc/policies/3D/reconstruct3D_pcg_policy.md`, `wave_model_policy.md`, `solve3D_addon_policy.md`;
skills `simple-main-nu-filt`, `simple-frac-update-trailing`, `simple-nonuniform-regularization`.
