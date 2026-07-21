---
name: simple-main-nano
description: Use when working in SIMPLE's src/main/nano subsystem, including nanoparticle and molecule data structures, graphene and time-series helpers, and the SINGLE-specific domain logic that supports liquid-cell and nanoparticle workflows.
---

# SIMPLE `src/main/nano`

This is the nanoparticle/SINGLE domain folder.

## Read First

- `simple_nanoparticle.f90`
- `simple_nanoparticle_utils.f90`
- `simple_atoms.f90`
- `simple_molecule_data.f90`
- `single_tseries_tracker.f90`
- `single_tseries_extractor.f90`

## Connections

- Works closely with `single/` commanders and `ui/single`
- Touches `image/`, `volume/`, and time-series workflows

## Working Rule

Expect nanoparticle-specific assumptions here that do not generalize to the main cryo-EM path.
