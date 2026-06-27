---
name: simple-main-volume
description: Use when working in SIMPLE's src/main/volume subsystem, including reconstruction objects, even/odd reconstruction handling, volume analysis, docking, interpolation, polar-Fourier volume correlation, and postprocessing policy such as automasking and nonuniform filtering.
---

# SIMPLE `src/main/volume`

This folder owns reconstruction and postprocessing objects.

## Read First

- `simple_reconstructor.f90`
- `simple_reconstructor_eo.f90`
- `simple_volume_postprocess_policy.f90`
- `simple_volanalyzer.f90`
- `simple_dock_vols.f90`

## Connections

- Very closely tied to `refine3D`, `rec3D`, `interp/`, `image/`, and `project/`
- Policy docs in `doc/policies/` are especially relevant here

## Working Rule

Preserve the distinction between reconstruction mechanics and workflow/postprocess policy.
