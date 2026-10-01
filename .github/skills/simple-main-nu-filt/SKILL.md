---
name: simple-main-nu-filt
description: Use when working in SIMPLE's src/main/nu_filt subsystem, including simple_nu_filter, candidate-bank setup, mask-packed objective costs, ordered-label Potts smoothing, NU local-resolution output, and cleanup of module-level filter state.
---

# SIMPLE `src/main/nu_filt`

This folder owns the volume-domain nonuniform filtering algorithm. Workflow
ownership still belongs to the callers: `simple_nu_state_filter` (gridding
`volassemble` and the PCG master), `postprocess_nu`, `nu_filt3D` and flex_pca.
This subsystem implements the filter state machine they call.

## Read First

- `simple_nu_filter.f90`
- `simple_nu_filter_bank.f90`
- `simple_nu_filter_potts.f90`
- `simple_nu_filter_apply.f90`
- `simple_nu_filter_stats.f90`
- `simple_nu_filter_state.f90`
- `simple_nu_filter_envmask.f90`
- `simple_nu_filter_evidence.f90`
- `simple_nu_filter_sharpen.f90`

## Core Lifecycle

Normal callers follow:

```fortran
call setup_nu_dmats(vol_even, vol_odd, mskdiam, aux_resolutions, aux_even, aux_odd, fsc_res=res0143)
call optimize_nu_cutoff_finds()
call nu_filter_vols(vol_even_nu, vol_odd_nu)
call cleanup_nu_filter()
```

The bank is the static ladder `lowpass_limits`. Given `fsc_res` it is cut at
`fsc_res/NU_BANK_FSC_HEADROOM` (at least two rungs); without it (`nu_filt3D`,
flex_pca) it is uncapped.

## Working Rules

- Treat `simple_nu_filter` as stateful for one call sequence. Always preserve
  cleanup of module allocatables and scratch cache files.
- `setup_nu_dmats` owns the NU support mask. It accepts `mskdiam` in Angstrom
  and constructs spherical support internally; callers must not supply density
  automasks or arbitrary logical envelopes.
- Keep objective construction mask-packed; values outside the NU support mask
  must not influence smoothing or label selection.
- The auxiliary (ML-regularized) pair is appended as the last label, never in
  a rung's place, and only when its Fourier index is at or beyond the finest
  retained rung. It shares that rung's Potts coordinate and its filtered pair
  is never cached. The last label is the finest bank member and the matching
  low-pass handoff.
- Ordered-label Potts smoothing is part of the current algorithm, not an optional
  user-facing mode.
- When the NU-evidence envelope arms the background clamp, derive it before
  `optimize_nu_cutoff_finds`. Run `build_nu_evidence_state` before
  `nu_filter_vols`, which releases the unary bank.
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
