---
name: simple-nonuniform-regularization
description: Use when working on SIMPLE's nonuniform filtering and volume-domain regularization path, including filt_mode=nonuniform or nonuniform_lpset, nu_refine, simple_nu_filter, nonuniform mask selection, _nu_filt/_nu_locres outputs, ML-regularized auxiliary replacement coupling, or volassemble postprocessing.
---

# SIMPLE Nonuniform Regularization

Use this skill for SIMPLE's nonuniform volume regularization/filtering path.
Treat it as a volume-domain postprocessing feature: particle search and matchers
consume filtered references, but Cartesian `volassemble` owns generation of the
nonuniform even/odd and merged reference volumes.

## Quick Start

1. Start with the public policy and code map in
   [references/nonuniform-regularization-map.md](./references/nonuniform-regularization-map.md).
2. Identify which layer the task touches:
   - User option and staging policy: `filt_mode`, `nu_refine`,
     `src/main/params`, `simple_abinitio_controller`.
   - Assembly orchestration: `commander_cartesian_volassemble` in
     `simple_commanders_rec_distr.f90`.
   - Mask precedence and automask cadence: `simple_vol_pproc_policy.f90`.
   - Numerical filter implementation: `src/main/nu_filt/simple_nu_filter.f90`.
   - Reference consumption by matchers: `simple_matcher_refvol_utils.f90`.
3. Preserve the execution point: nonuniform filtering runs after even/odd
   assembly, sampling-density correction, low-resolution reinsertion, gridding
   correction, merged-volume write, automask planning, and trailing blend when active.
4. Keep the filtered outputs as derived products. The base state volumes remain
   the primary reconstruction artifacts; `_nu_filt` files are filtered
   references and `_nu_locres` files are diagnostics.

## Working Rules

- Do not move nonuniform filtering into matcher/search code. Matchers may prefer
  NU references when available, but `volassemble` generates them.
- Keep mask policy centralized in `plan_state_postprocess`. Nonuniform mask
  precedence is fresh compatible automask, existing compatible state automask,
  then spherical mask from `mskdiam`.
- Preserve shared-memory and distributed parity by changing the Cartesian
  assembly path, not only a single execution mode.
- Treat `simple_nu_filter` as stateful during a call sequence. The required
  lifecycle is `setup_nu_dmats`, `optimize_nu_cutoff_finds`, optional
  high-resolution shell extension, `nu_filter_vols`, and `cleanup_nu_filter`.
- Always clean up cache/state on normal and error-prone paths. The module uses
  disk cache files named `nu_filter_cache_even_k_*.mrc` and
  `nu_filter_cache_odd_k_*.mrc` plus module-level allocatables.
- When `ml_reg=yes` and `nu_refine=no`, use the `_unfil` even/odd pair as the
  base nonuniform input and let the ML-regularized pair replace the finest
  static-bank member when it is effectively finer.
- When `nu_refine=yes`, do not pass the ML-regularized auxiliary replacement
  pair. High-resolution shell extension owns the resolution-extension experiment.
- Keep `filt_mode=nonuniform` distinct from `filt_mode=uniform` and
  `filt_mode=fsc`. Plain nonuniform mode does not set `l_lpset`; it generates
  local voxelwise-filtered references in assembly.
- Treat `filt_mode=nonuniform_lpset` as NU filtering plus LP-set matching
  topology and selected-LP handoff; do not collapse it into plain `nonuniform`.

## Common Tasks

### Explain the approach

Describe it as local selection among candidate filtered even/odd pairs:

- build a low-pass candidate bank from the unfiltered even/odd pair
- optionally replace the finest bank member with a finer auxiliary pair
- compute mask-packed voxelwise objective costs
- smooth each candidate objective inside the logical support mask
- choose the best candidate per voxel
- apply ordered-label Potts smoothing
- synthesize filtered even/odd outputs, the merged `_nu_filt` volume, and the
  `_nu_locres` diagnostic map

### Modify mask behavior

Read:

- `doc/policies/nonuniform_filtering_policy.md`
- `doc/policies/automasking_policy.md`
- `src/main/volume/simple_vol_pproc_policy.f90`
- `src/main/commanders/simple/simple_commanders_rec_distr.f90`

Keep compatibility checks on dimensions and sampling distance. Do not bypass
state-specific automask naming: `automask3D_stateNN.mrc`.

### Modify filter optimization or performance

Read:

- `src/main/nu_filt/simple_nu_filter.f90`
- `src/main/nu_filt/simple_nu_filter_*.f90`
- `production/tests/simple_test_nu_filter.f90`
- the implementation notes in [references/nonuniform-regularization-map.md](./references/nonuniform-regularization-map.md)

Watch for mask-packed arrays, temporary full-volume buffers, disk-backed cache
traffic, and OpenMP loop shape. If replacing cache files with in-memory buffers,
preserve behavior for auxiliary replacement, stats, local-resolution maps, and cleanup.

### Modify matcher consumption of filtered references

Inspect `src/main/strategies/search/simple_matcher_refvol_utils.f90`.

The intended behavior is:

- in plain `nonuniform`, try `_nu_filt` even/odd files first
- in `nonuniform_lpset` with LP-set matching, prefer the merged `_nu_filt`
  registration reference
- fall back to regular references before the first `volassemble` has generated
  filtered products
- avoid applying another normal low-pass filter when `l_nonuniform_mode` is active

## Language To Prefer

- "volume-domain nonuniform filtering"
- "local low-pass candidate selection"
- "support-mask precedence"
- "derived `_nu_filt` reference products"
- "derived `_nu_locres` diagnostic products"
- "assembly-owned postprocessing"
- "auxiliary replacement pair for ML-regularized volumes"

Avoid:

- describing nonuniform filtering as a particle-domain regularizer
- conflating `filt_mode=nonuniform` with FSC or uniform low-pass filtering
- treating automask and nonuniform-mask policy as separate ad hoc decisions
- implying `_nu_filt` volumes replace the base reconstruction outputs
