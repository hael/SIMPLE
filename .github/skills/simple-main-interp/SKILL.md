---
name: simple-main-interp
description: Use when working in SIMPLE's src/main/interp subsystem, including gridding, Kaiser-Bessel interpolation, edge/window functions, and the neutral Cartesian Fourier lattice layer (embedding, extraction, packed KB gathers) used by reconstruction, class averaging, polarization and Cartesian alignment.
---

# SIMPLE `src/main/interp`

This folder owns interpolation and gridding primitives.

## Read First

- `doc/policies/KB_Interpolation_Policy.md`: convolution interpolation on the oversampled Fourier grid
- `simple_kbinterpol.f90`
- `simple_gridding.f90`
- `simple_cartesian_fourier.f90`: neutral Cartesian Fourier lattice embedding, extraction and packed KB gathers
- `simple_winfuns.f90`
- `simple_edges_sqwins.f90`

## Connections

- `image/` (polarization in `simple_image_polar`, `simple_projector`), which feeds `pftc/`
- `volume/` (gridding and PCG reconstructors, `simple_pcg_lattice`), `class/` (class averager), `cftc/` (`simple_cartft_calc`), and the `flex` fits
- Tests: `unit_numerics` ("Kaiser-Bessel kernel"), `unit_cart_align3D` ("Cartesian Fourier")

## Working Rule

Small numerical changes here can have wide scientific effects, so verify downstream usage before refactoring.
