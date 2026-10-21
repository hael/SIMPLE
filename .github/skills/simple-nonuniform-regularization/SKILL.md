---
name: simple-nonuniform-regularization
description: Use when working on SIMPLE's nonuniform filtering and volume-domain regularization path, including filt_mode=nonuniform or nonuniform_lpset, simple_nu_filter and simple_nu_state_filter, spherical NU support, _nu_filt/_nu_locres outputs, the NU matching low-pass handoff, ML-regularized auxiliary member coupling, or NU postprocessing in volume assembly (volassemble_exec) and the PCG master.
---

# SIMPLE Nonuniform Regularization

Use this skill for SIMPLE's nonuniform volume regularization/filtering path.
Treat it as a volume-domain postprocessing feature: particle search and matchers
consume filtered references, but volume assembly (gridding `volassemble_exec`
and the PCG master, through `nonuniform_filter_state` in
`simple_nu_state_filter`) owns generation of the nonuniform even/odd and merged
reference volumes.

## Quick Start

1. Start with the public policy and code map in
   [references/nonuniform-regularization-map.md](./references/nonuniform-regularization-map.md).
2. Identify which layer the task touches:
   - User option and staging policy: `filt_mode`,
     `src/main/params`, `simple_solve3D_controller`.
   - Assembly orchestration: `volassemble_exec` in
     `simple_volassemble_strategy.f90` (run by the distributed `refine3D`
     master on its own builder; `commander_volassemble` is a thin wrapper for
     command-line callers), the PCG master in `simple_rec3D_pcg_strategy.f90`
     (project handoff `filter_pcg_nonuniform_maps`), and the shared per-state
     `nonuniform_filter_state` in `src/main/volume/simple_nu_state_filter.f90`.
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
- Preserve parity across execution modes and both reconstruction backends by
  changing the shared `nonuniform_filter_state` (or `volassemble_exec`), not
  one caller.
- Treat `simple_nu_filter` as stateful during a call sequence. The required
  lifecycle is `setup_nu_dmats`, `optimize_nu_cutoff_finds`,
  `nu_filter_vols(vol_even_nu, vol_odd_nu, vol_src_even, vol_src_odd)`, and
  `cleanup_nu_filter`, which releases the module-level allocatables. The
  filter writes no scratch files.
- Keep the source pair alive and unchanged from `setup_nu_dmats` to
  `nu_filter_vols`: the composition filters it again per label. With
  `pcg_solvent=yes` the label field of the base pair is applied to the
  solvent-prior pair instead.
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
- filter each candidate pair and score it at once with mask-packed voxelwise
  objective costs: the squared cross-half prediction error at one global
  noise level (no local noise estimate; the level never changes a voxel's
  winner)
- smooth each candidate objective inside spherical support
- choose the label per voxel coarse to fine, like for like (both costs
  smoothed at the finer candidate's radius)
- apply ordered-label Potts smoothing, which only coarsens
- synthesize filtered even/odd outputs by filtering the source pair again per
  label, then the merged `_nu_filt` volume and the `_nu_locres` diagnostic map

### Modify support or automask behavior

Read:

- `doc/policies/NU/nonuniform_filtering_policy.md`
- `doc/policies/3D/automasking_policy.md`
- `src/main/volume/simple_vol_pproc_policy.f90` (`state_mask_is_compatible`)
- `src/main/volume/simple_nu_state_filter.f90`
- `src/main/strategies/parallelization/simple_volassemble_strategy.f90`

For NU support changes, preserve the spherical-only `setup_nu_dmats` API unless
an explicit policy change includes a statistically valid null estimator and a
memory plan. For automask consumers, keep compatibility checks on dimensions
and sampling distance and preserve state-specific naming:
`automask3D_stateNN.mrc`.

### Modify filter optimization or performance

Read:

- `src/main/nu_filt/simple_nu_filter.f90`
- `src/main/nu_filt/simple_nu_filter_*.f90`
- `src/main/nu_filt/simple_nu_filter_tester.f90` (sub-suite `nonuniform
  filtering` of `unit_reconstruction`) and `simple_exec prg=nu_filt3D` on a
  refined map
- the implementation notes in [references/nonuniform-regularization-map.md](./references/nonuniform-regularization-map.md)

Watch for mask-packed arrays, temporary full-volume buffers, and OpenMP loop
shape. Candidate pairs are never stored: the setup scores each right after
filtering it and the composition filters the source pair again. The level loop
reuses the setup's `dmats_mask` and `norms_mask`, the Potts sweeps walk packed
per-colour voxel lists, and the tent smoothing skips zero lines and is held bit
for bit to a reference by its unit test. Read the filter's split in the
`VOLASSEMBLE` bench report (filtering, objective, optimisation with level
smoothings and Potts sweeps, composition) before and after.

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
