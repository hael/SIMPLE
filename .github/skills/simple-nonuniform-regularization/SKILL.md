---
name: simple-nonuniform-regularization
description: Use when working on SIMPLE's nonuniform filtering and volume-domain regularization path, including filt_mode=nonuniform or nonuniform_lpset, simple_nu_filter, spherical NU support, _nu_filt/_nu_locres outputs, ML-regularized auxiliary member coupling, or volassemble postprocessing.
---

# SIMPLE Nonuniform Regularization

Use this skill for SIMPLE's nonuniform volume regularization/filtering path.
Treat it as a volume-domain postprocessing feature: particle search and matchers
consume filtered references, but volume assembly (gridding `volassemble` and the
PCG master, through `simple_nu_state_filter`) owns generation of the
nonuniform even/odd and merged reference volumes.

## Quick Start

1. Start with the public policy and code map in
   [references/nonuniform-regularization-map.md](./references/nonuniform-regularization-map.md).
2. Identify which layer the task touches:
   - User option and staging policy: `filt_mode`,
     `src/main/params`, `simple_abinitio_controller`.
   - Assembly orchestration: `commander_volassemble` in
     `simple_commanders_rec_distr.f90`, the PCG master in
     `simple_rec3D_pcg_strategy.f90`, and the shared per-state
     `simple_nu_state_filter.f90`.
   - Mask artifact compatibility: `simple_vol_pproc_policy.f90`.
   - Spherical NU support: `src/main/nu_filt/simple_nu_filter_bank.f90`.
   - Numerical filter implementation: `src/main/nu_filt/simple_nu_filter.f90`.
   - Reference consumption by matchers: `simple_matcher_refvol_utils.f90`.
3. Preserve the execution point: nonuniform filtering runs after the state's
   half maps are restored (trailing blend applied when active), the FSC is
   computed, and the base volumes are written. Low-resolution even/odd
   insertion is matcher-side reference preparation and never feeds NU filtering.
4. Keep the filtered outputs as derived products. The base state volumes remain
   the primary reconstruction artifacts; `_nu_filt` files are filtered
   references and `_nu_locres` files are diagnostics.

## Working Rules

- Do not move nonuniform filtering into matcher/search code. Matchers may prefer
  NU references when available, but volume assembly generates them.
- Keep spherical support as an invariant of `setup_nu_dmats`: it accepts
  `mskdiam` in Angstrom and constructs the logical sphere internally. Do not
  reintroduce caller-supplied automasks or arbitrary logical support.
- Keep the density automask independent of NU support: assembly generates it
  for `envfsc` and, under active `automsk`, for the filter-field background and
  reference mask; `state_mask_is_compatible` checks reused artifacts.
- Preserve shared-memory and distributed parity by changing the Cartesian
  assembly path, not only a single execution mode.
- Treat `simple_nu_filter` as stateful during a call sequence. The required
  lifecycle is `setup_nu_dmats`, `optimize_nu_cutoff_finds`,
  `nu_filter_vols`, and `cleanup_nu_filter`.
- Always clean up cache/state on normal and error-prone paths. The module uses
  disk cache files named `nu_filter_cache_even_k_*.mrc` and
  `nu_filter_cache_odd_k_*.mrc` plus module-level allocatables.
- When `ml_reg=yes`, use the `_unfil` even/odd pair as the base nonuniform
  input and offer the ML-regularized pair as the auxiliary member. It is
  appended as the last label, sharing the finest retained rung's Potts
  coordinate, once its FSC=0.143 resolution is at or beyond that rung; it never
  replaces a rung. The last label is the finest member and the matching handoff.
- Keep `filt_mode=nonuniform` distinct from `filt_mode=uniform` and
  `filt_mode=fsc`. Plain nonuniform mode does not set `l_lpset`; it generates
  local voxelwise-filtered references in assembly.
- Treat `filt_mode=nonuniform_lpset` as NU filtering plus LP-set matching
  topology and selected-LP handoff; do not collapse it into plain `nonuniform`.

## Common Tasks

### Explain the approach

Describe it as local selection among candidate filtered even/odd pairs:

- build a low-pass candidate bank from the unfiltered even/odd pair
- construct spherical objective support from `mskdiam`
- optionally append the auxiliary pair as the last (finest) bank member
- compute mask-packed voxelwise objective costs
- smooth each candidate objective inside spherical support
- choose the best candidate per voxel
- apply ordered-label Potts smoothing
- synthesize filtered even/odd outputs, the merged `_nu_filt` volume, and the
  `_nu_locres` diagnostic map

### Modify support or automask behavior

Read:

- `doc/policies/NU/nonuniform_filtering_policy.md`
- `doc/policies/3D/automasking_policy.md`
- `src/main/volume/simple_vol_pproc_policy.f90`
- `src/main/commanders/simple/simple_commanders_rec_distr.f90`

For NU support changes, preserve the spherical-only `setup_nu_dmats` API unless
an explicit policy change includes a statistically valid null estimator and a
memory plan. For automask consumers, keep compatibility checks on dimensions
and sampling distance and preserve state-specific naming:
`automask3D_stateNN.mrc`.

### Modify filter optimization or performance

Read:

- `src/main/nu_filt/simple_nu_filter.f90`
- `src/main/nu_filt/simple_nu_filter_*.f90`
- `simple_exec prg=nu_filt3D` on a refined map (the standalone `nu_filter` test
  program was deleted: it asserted nothing; see `doc/policies/test_environment_policy.md`)
- the implementation notes in [references/nonuniform-regularization-map.md](./references/nonuniform-regularization-map.md)

Watch for mask-packed arrays, temporary full-volume buffers, disk-backed cache
traffic, and OpenMP loop shape. If replacing cache files with in-memory buffers,
preserve behavior for the auxiliary member, stats, local-resolution maps, and cleanup.

### Modify matcher consumption of filtered references

Inspect `src/main/strategies/search/simple_matcher_refvol_utils.f90`.

The intended behavior is:

- in plain `nonuniform`, try `_nu_filt` even/odd files first
- in `nonuniform_lpset` with LP-set matching, prefer the merged `_nu_filt`
  registration reference
- fall back to regular references before the first `volassemble` has generated
  filtered products
- avoid applying another normal low-pass filter when `l_nonuniform_mode` is active
  and the `_nu_filt` references exist; the raw half-map fallback gets the
  FSC-optimal filter

## Language To Prefer

- "volume-domain nonuniform filtering"
- "local low-pass candidate selection"
- "spherical NU support"
- "independent automask lifecycle"
- "derived `_nu_filt` reference products"
- "derived `_nu_locres` diagnostic products"
- "assembly-owned postprocessing"
- "auxiliary member (the ML-regularized pair), appended as the last label"

Avoid:

- describing nonuniform filtering as a particle-domain regularizer
- conflating `filt_mode=nonuniform` with FSC or uniform low-pass filtering
- implying density automasks define the NU objective domain
- using the NU-evidence envelope as NU objective support or as the PCG solve support
- implying `_nu_filt` volumes replace the base reconstruction outputs
