---
name: simple-main-image
description: Use when working in SIMPLE's src/main/image subsystem, including the core image type and its submodules (FFT and Fourier-space operations, geometry, filtering, masks, I/O, polarization), the binary and mask extensions, expanded Fourier transforms for motion correction, projector classes, and related image-processing infrastructure.
---

# SIMPLE `src/main/image`

This is one of the central scientific subsystems.

## Read First

- `simple_image.f90`
- `simple_image_core.f90`: `new`/`kill`, the FFTW buffer and plans
- `simple_image_fft.f90`
- `simple_image_access.f90`
- `simple_image_filt.f90`
- `simple_projector.f90`, `simple_projector_pft_batch.f90`

## Structure

- `simple_image.f90` declares the `image` type and its procedure surface; behaviour lives in submodules (`*_core`, `*_access`, `*_arith`, `*_fft`, `*_filt`, `*_geom`, `*_io`, `*_polar`, `*_norm`, `*_ops`, ...)
- `image_bin` (`simple_image_bin.f90`) and `image_msk` (`simple_image_msk.f90`) extend `image`
- Projector code lives alongside the image primitives because it is tightly coupled
- `simple_ft_expanded.f90` and `simple_ftexp_shsrch.f90` serve motion correction; `simple_fplane_grid.f90` the reconstruction's plane insertion

## Guardrails

- One FFTW buffer is seen as `cmat` (Fourier) and `rmat` (real, with padding rows in the first dimension). In real space the padding holds zeros. Outside the image class use `get_rmat_ptr` (the box only); `get_rmat_ptr_padded` is for the few callers that need the FFT layout
- Plans are made with `FFTW_ESTIMATE` inside an `omp critical` (the planner is not thread-safe); keep plan creation and destruction there
- `wthreads` defaults to on in `new`; an image gets threaded plans only with `wthreads`, more than one plan thread, and a dimension of at least 200 (or `ldim(3)` of at least 100). Per-thread work images are made with `wthreads=.false.`
- The Fourier crop round trip (crop, inverse FFT, forward FFT) is exact in mathematics, not in single precision: never promise bitwise identity across it

## Connections

- Core dependency for `motion`, `ctf`, `class`, `refine2D`, `refine3D`, `volume`, and `pftc`
- The matchers' particle buffer stores decoded images in an `image_slab` (`src/fileio/simple_image_slab.f90`), read back bit for bit; the 3D record buffer of cropped observations lives in `strategies/search/simple_matcher_ptcl_io.f90`
- Tests: `simple_test_exec test=unit_image`

## Working Rule

Before adding a new standalone helper, check whether a matching type-bound method or submodule already exists.
