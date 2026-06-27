---
name: simple-main-image
description: Use when working in SIMPLE's src/main/image subsystem, including the core image type, FFT and Fourier-space operations, geometry, filtering, masks, I/O, projector classes, and related image-processing infrastructure.
---

# SIMPLE `src/main/image`

This is one of the central scientific subsystems.

## Read First

- `simple_image.f90`
- `simple_image_core.f90`
- `simple_image_fft.f90`
- `simple_image_access.f90`
- `simple_image_filt.f90`
- `simple_projector.f90`

## Structure

- `simple_image.f90` declares the big public type and procedure surface
- Many behaviors live in submodules such as `*_core`, `*_fft`, `*_filt`, `*_io`, `*_polar`, `*_geom`
- Projector code lives alongside the image primitives because it is tightly coupled

## Connections

- Core dependency for `motion`, `ctf`, `cluster2D`, `refine3D`, `volume`, and `pftc`

## Working Rule

Before adding a new standalone helper, check whether a matching type-bound method or submodule already exists.
