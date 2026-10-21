---
name: simple-main-ori
description: Use when working in SIMPLE's src/main/ori subsystem, including orientation records, orientation collections, point-group symmetry, sampling and update counts, transforms, I/O, statistics, neighborhood relations, and weighting logic used by 2D and 3D search workflows.
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
- `simple_sym.f90`: point-group symmetries
- `src/defs/simple_defs_ori.f90`: the named per-particle slots

## Structure

- `simple_ori.f90` is the record type: named per-particle values in fixed slots (`pparms`), other values in hash tables (`htab`, `chtab` for strings)
- `simple_oris.f90` declares the collection; its submodules (`simple_oris_life`, `_getters`, `_setters`, `_io`, `_sampling`, `_stats`, `_neig`, `_dists`, `_reshape`, `_transform`) hold the bulk operations
- `simple_ori_utils.f90` (Euler/rotation-matrix helpers) and `simple_oris_utils.f90` (stateless utilities over a collection, such as a state's orientation distribution) are plain modules

## Guardrails

- A named slot is an enumerator in `simple_defs_ori.f90` with its flag in `get_oriparam_ind`/`get_oriparam_flag`; append after `I_LAST_NAMED_ORIPARAM` and keep `N_PTCL_ORIPARAMS` within `N_PTCL_RECORD_REALS`, so the on-disk record width does not change and older projects still read (`simple_binoris`)
- Read and write values through the getters and setters, never by parsing strings

## Connections

- Central to `project/`, `refine3D`, `refine2D`, `pick`, and probability/search code
- Tests: `simple_test_exec test=unit_ori`

## Working Rule

Be careful about stateful fields and persistence semantics because orientation records are core workflow currency.
