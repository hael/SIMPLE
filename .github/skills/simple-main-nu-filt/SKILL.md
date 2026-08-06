---
name: simple-main-nu-filt
description: Use when working in SIMPLE's src/main/nu_filt subsystem, including simple_nu_filter, candidate-bank setup, mask-packed objective costs, ordered-label Potts smoothing, high-resolution shell extension, NU local-resolution output, and cleanup of module-level filter state.
---

# SIMPLE `src/main/nu_filt`

This folder owns the volume-domain nonuniform filtering algorithm. Workflow
ownership still belongs to `volassemble`; this subsystem implements the filter
state machine it calls.

## Read First

- `simple_nu_filter.f90`
- `simple_nu_filter_bank.f90`
- `simple_nu_filter_potts.f90`
- `simple_nu_filter_extend.f90`
- `simple_nu_filter_apply.f90`
- `simple_nu_filter_stats.f90`
- `simple_nu_filter_state.f90`
- `simple_nu_filter_envmask.f90`

## Core Lifecycle

Normal callers follow:

```fortran
call setup_nu_dmats(vol_even, vol_odd, mskdiam, aux_resolutions, aux_even, aux_odd)
call optimize_nu_cutoff_finds()
call extend_nu_filter_highres_shell_next(...)
call nu_filter_vols(vol_even_nu, vol_odd_nu)
call cleanup_nu_filter()
```

The extension step is optional and controlled by workflow policy such as
`nu_refine=yes`.

## Working Rules

- Treat `simple_nu_filter` as stateful for one call sequence. Always preserve
  cleanup of module allocatables and scratch cache files.
- `setup_nu_dmats` owns the NU support mask. It accepts `mskdiam` in Angstrom
  and constructs spherical support internally; callers must not supply density
  automasks or arbitrary logical envelopes.
- Keep objective construction mask-packed; values outside the NU support mask
  must not influence smoothing or label selection.
- Auxiliary even/odd pairs replace the finest retained bank member only when
  their effective resolution is finer; they are not sidecar labels.
- Ordered-label Potts smoothing is part of the current algorithm, not an optional
  user-facing mode.
- High-resolution shell extension tests one frontier shell at a time, persists
  accepted depth by state through `volassemble`, and should stop conservatively
  when support is missing or acceptance is too weak.
- NU-evidence envelope derivation must run before `nu_filter_vols`, which
  releases unary storage. The current evidence baseline covers the static bank;
  accepted `nu_refine` extension shells are not yet incorporated.

## Mask Ownership

- Spherical `mskdiam` support owns the NU objective domain.
- Density-derived `automask3D_stateNN.mrc` remains independent and may be
  consumed by `envfsc` or current `envref` behavior.
- A future NU-evidence envelope must use its own artifact and must never feed
  FSC correction or replace spherical NU support.
