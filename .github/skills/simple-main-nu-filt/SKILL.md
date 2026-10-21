---
name: simple-main-nu-filt
description: Use when working in SIMPLE's src/main/nu_filt subsystem, including simple_nu_filter and its submodules, candidate-bank setup, the source pair of nu_filter_vols, like-for-like label selection, mask-packed objective costs, ordered-label Potts smoothing, tent smoothing, NU sub-timers, NU local-resolution output, the NU-evidence envelope, and cleanup of module-level filter state.
---

# SIMPLE `src/main/nu_filt`

This folder owns the volume-domain nonuniform filtering algorithm. Workflow
ownership belongs to the callers: `nonuniform_filter_state` in
`src/main/volume/simple_nu_state_filter.f90` (volume assembly,
`volassemble_exec` in `simple_volassemble_strategy.f90`, and the PCG master in
`simple_rec3D_pcg_strategy.f90`), `exec_nu_filt3D` in
`simple_commanders_resolest.f90`, and `postprocess_nu`. This subsystem
implements the filter state machine they call.

## Read First

- `simple_nu_filter.f90` (interfaces and constants)
- `simple_nu_filter_vars.f90` (the module-level working state)
- `simple_nu_filter_bank.f90` (`setup_nu_dmats`, `optimize_nu_cutoff_finds`)
- `simple_nu_filter_state.f90` (filtering, smoothing, packing, cleanup)
- `simple_nu_filter_apply.f90` (`nu_filter_vols`)
- `simple_nu_filter_potts.f90`
- `simple_nu_filter_stats.f90`
- `simple_nu_filter_envmask.f90`, `simple_nu_filter_evidence.f90`,
  `simple_nu_filter_sharpen.f90`
- `src/utils/filter/simple_tent_smooth.f90`
- `doc/policies/NU/nonuniform_filtering_policy.md` (sections 8-9)

## Core Lifecycle

Normal callers follow:

```fortran
call setup_nu_dmats(vol_even, vol_odd, mskdiam, aux_resolutions, aux_even, aux_odd, fsc_res=res0143)
call optimize_nu_cutoff_finds()
call nu_filter_vols(vol_even_nu, vol_odd_nu, vol_src_even, vol_src_odd)
call cleanup_nu_filter()
```

The bank is the static ladder `lowpass_limits`. Given `fsc_res` it is cut at
`fsc_res/NU_BANK_FSC_HEADROOM` (at least two rungs); without it (`nu_filt3D`)
it is uncapped. `postprocess_nu` stops after `optimize_nu_cutoff_finds` and
`build_nu_evidence_state`; it never composes.

## Working Rules

- Treat `simple_nu_filter` as stateful for one call sequence. Preserve
  `cleanup_nu_filter` of every module allocatable. The filter writes no
  scratch files: `setup_nu_dmats` filters each candidate pair from the
  transforms of the edge-tapered halves and scores it at once, and
  `nu_filter_vols` filters its source pair again per label with the same
  taper, transform and Butterworth filter, so the composed voxels are the
  candidates'.
- The source pair of `nu_filter_vols` must stay alive and unchanged from
  setup to composition (`nonuniform_filter_state` kills the base pair only
  after it). An apply pair other than the setup pair (the `pcg_solvent=yes`
  solvent-prior pair) receives the label field estimated on the base pair.
- `setup_nu_dmats` owns the NU support mask. It accepts `mskdiam` in Angstrom
  and constructs spherical support internally; callers must not supply density
  automasks or arbitrary logical envelopes.
- Keep objective construction mask-packed; values outside the NU support mask
  must not influence smoothing or label selection.
- Label selection is like-for-like, coarse to fine: at each level the
  incumbent and the candidate are smoothed at the candidate's radius from
  `raw_dmats_mask`, and the candidate wins only when strictly cheaper. The
  candidate's own column is the setup's `dmats_mask`, and the smoothing norm
  of that radius is the setup's `norms_mask` (`load_nu_smooth_norm`); do not
  smooth either again.
- The Potts sweeps walk per-colour voxel lists packed once and may only move
  a voxel to a coarser label than the one it entered with. Their beta comes
  from a real-valued reduction, so the label field is reproducible at a fixed
  thread count, not across thread counts.
- `tent_smooth_3d` writes all-zero rows, planes and slabs as zeros and filters
  Y and Z a plane at a time, each column keeping the one-column operation
  sequence. `test_tent_smooth_exact` holds it bit for bit to a column-by-column
  reference; keep that identity when changing it.
- The auxiliary (ML-regularized) pair is appended as the last label, never in
  a rung's place, and only when its Fourier index is at or beyond the finest
  retained rung. It shares that rung's Potts coordinate and is stashed as a
  copy (`stash_aux_volumes`). The last label is the finest bank member and the
  matching low-pass handoff (`get_nu_filter_bank_finest_lp`).
- Ordered-label Potts smoothing is part of the current algorithm, not an optional
  user-facing mode.
- When the NU-evidence envelope arms the background clamp, derive it after
  `setup_nu_dmats` and before `optimize_nu_cutoff_finds`; the clamp is applied
  before and after the Potts sweeps. Run `build_nu_evidence_state` (or the
  envelope) before `nu_filter_vols`, which releases the unary bank.
- Timings: under `L_BENCH_GLOB` the sequence accumulates `nu_rt_*` (candidate
  filtering, objective, optimisation split into level smoothings and Potts
  sweeps, composition); `get_nu_filter_timings` feeds the `VOLASSEMBLE`
  report through `nu_state_filter_timings`. Time new phases the same way.
- Keep the standalone `nu_filt3D` envelope interface to two shape controls:
  `nu_msk_sig` for evidence threshold and `amsklp` for physical evidence scale.
  Absolute evidence, zero density weighting, MRF beta 1, and the 0.1 component
  fraction are fixed production policy, even though the internal API retains
  diagnostic variants. Preserve their semantics beside the constants: beta
  controls boundary smoothness; density weighting can retain strong but poorly
  ordered density; scale-free evidence can protect weak ordered density from a
  high-contrast core; and the component fraction removes small components
  relative to the largest.
- Express standalone envelope morphology in Angstrom: 1 A binary growth and a
  6 A cosine edge. Convert those lengths to the nearest voxel count at the
  input-map sampling, with a minimum of one voxel.

## Mask Ownership

- Spherical `mskdiam` support owns the NU objective domain.
- The density envelope (`automask3D_stateNN.mrc`) is independent of NU support.
  It feeds `envfsc`, and under `automsk=yes` (or as the `automsk=nu` fallback)
  it fixes the filter-field background and masks the `_nu_filt` references.
- The NU-evidence envelope has its own artifact, `nu_envmask3D_stateNN.mrc`,
  written by `write_nu_evidence_envmask` under active automasking. It is a
  diagnostic under `automsk=yes`. Under `automsk=nu`, when valid, it is the
  background and reference mask, and its lag-one artifact is the gridding
  `envfsc` FSC mask (density fallback). It never replaces spherical NU support
  and is never the PCG solve support.

## Tests

`simple_test_exec test=unit_reconstruction`, sub-suite `nonuniform filtering`
(`simple_nu_filter_tester.f90`): the objective's closed form and noise level,
the label field on a two-resolution phantom, the per-label composition, the
evidence envelope, the tent smoothing identity and the label field repeated at
every processor. End to end: `simple_exec prg=nu_filt3D` on a refined half-map
pair, before and after the change.
