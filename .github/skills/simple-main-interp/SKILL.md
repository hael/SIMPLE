---
name: simple-main-interp
description: Use when working in SIMPLE's src/main/interp subsystem, including gridding, Kaiser-Bessel interpolation, edge/window functions, and interpolation support used by reconstruction and Fourier-domain workflows.
---

# SIMPLE `src/main/interp`

This folder owns interpolation and gridding primitives.

## Read First

- `simple_kbinterpol.f90`
- `simple_gridding.f90`
- `simple_winfuns.f90`
- `simple_edges_sqwins.f90`

## Connections

- Heavily used by `image/`, `volume/`, `pftc/`, and reconstruction paths

## Working Rule

Small numerical changes here can have wide scientific effects, so verify downstream usage before refactoring.
