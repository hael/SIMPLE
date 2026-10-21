---
name: simple-main-volume
description: Use when working in SIMPLE's src/main/volume subsystem, including the gridding and PCG reconstructors, half-map diagnostics, trailing-reconstruction chains, the frozen accumulators of solve3D_addon, assembly-owned nonuniform filtering of a state, the PCG priors, volume analysis, docking, clustering, interpolation, polar-Fourier volume correlation, and symmetry search.
---

# SIMPLE `src/main/volume`

This folder owns reconstruction objects and volume analysis.

## Read First

- `simple_reconstructor.f90`: gridding (KB convolution) reconstruction
- `simple_reconstructor_pcg.f90`, `simple_pcg_lattice.f90`: the PCG backend (`rec_backend=pcg`); its solver engine is `opt/simple_pcg_solver.f90`
- `simple_halfmap_diagnostics.f90`: half-map FSC, cFAR and resolution, shared by both backends
- `simple_nu_state_filter.f90`: nonuniform filtering of one state's half-map pair, called by assembly
- `simple_trail_chain_manifest.f90`, `simple_frozen_accum.f90`: trailing-chain identity; `solve3D_addon` frozen sets
- `simple_vol_pproc_policy.f90` (per-state mask compatibility), `simple_volanalyzer.f90`, `simple_dock_vols.f90`
- `doc/policies/3D/reconstruct3D_pcg_policy.md`, `doc/policies/NU/nonuniform_filtering_policy.md`, `doc/policies/3D/automasking_policy.md`, `doc/policies/3D/refine3D_policy.md`

## Ownership

- Volume assembly (even, odd and merged state volumes from partial reconstructions, the trailing blend, the NU filter) is `src/main/strategies/parallelization/simple_volassemble_strategy.f90`: `volassemble_exec` on the caller's `(params, build)` (the distributed refine3D master), `run_volassemble` on its own (shared-memory refine3D, the rec3D service, `simple_commanders_rec_distr.f90`), `filter_pcg_nonuniform_maps` for PCG maps
- The PCG routes (shared-memory solve, worker accumulation, the master's reduction and solve) are `simple_rec3D_pcg_strategy.f90`; reconstruction requests go through `simple_rec3D_service.f90`
- Reconstruction takes the backend-neutral cropped observation (`prep_rec_observation`) or the same observation from the pass's records (`take_records`, `get_rec_observation` in `strategies/search/simple_matcher_ptcl_io.f90`), bit for bit; keep the two identical

## Connections

- Very closely tied to `refine3D`, `rec3D`, `interp/`, `image/`, `nu_filt/`, and `project/`
- Tests: `simple_test_exec test=unit_reconstruction`; `lib_reconstruction` holds the PCG half-set and sparse-sampling suites

## Working Rule

Preserve the distinction between reconstruction mechanics (here) and assembly/workflow policy (strategies).
