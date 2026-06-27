---
name: simple-main-ori
description: Use when working in SIMPLE's src/main/ori subsystem, including orientation records, orientation collections, sampling, transforms, I/O, statistics, neighborhood relations, and weighting logic used by 2D and 3D search workflows.
---

# SIMPLE `src/main/ori`

This folder owns orientation types and their operations.

## Read First

- `simple_ori.f90`
- `simple_oris.f90`
- `simple_oris_getters.f90`
- `simple_oris_setters.f90`
- `simple_oris_transform.f90`
- `simple_oris_sampling.f90`

## Structure

- `simple_ori.f90` is the elemental record type
- `simple_oris.f90` and submodules handle collections and bulk operations

## Connections

- Central to `project/`, `refine3D`, `cluster2D`, `pick`, and probability/search code

## Working Rule

Be careful about stateful fields and persistence semantics because orientation records are core workflow currency.
