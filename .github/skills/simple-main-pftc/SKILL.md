---
name: simple-main-pftc
description: Use when working in SIMPLE's src/main/pftc subsystem, including polar Fourier transform calculation, correlation, geometry, memoization, CTF-aware paths, restore/state operations, and search support for alignment workflows.
---

# SIMPLE `src/main/pftc`

This folder is central to polar-Fourier search workflows.

## Read First

- `simple_polarft_calc.f90`
- `simple_polarft_core.f90`
- `simple_polarft_corr.f90`
- `simple_polarft_geom.f90`
- `simple_polarft_memo.f90`
- `simple_polarft_ctf.f90`

## Structure

- The public type is declared centrally and much behavior lives in submodules
- Memoization and restore/state code are first-class concerns here

## Connections

- Tight coupling with `image/`, `interp/`, `ctf/`, `ori/`, and `strategies/search`

## Working Rule

Changes here often affect alignment performance and correctness simultaneously, so treat caching and geometry logic as part of the design, not incidental detail.
