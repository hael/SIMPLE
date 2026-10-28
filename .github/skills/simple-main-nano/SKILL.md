---
name: simple-main-nano
description: Use when working in SIMPLE's src/main/nano subsystem, including the nanoparticle type and atomic model building, species discovery numerics, the atoms (PDB/CIF, density simulation) and C-alpha finder modules, graphene and time-series helpers, and the SINGLE-specific domain logic that supports liquid-cell and nanoparticle workflows.
---

# SIMPLE `src/main/nano`

This is the nanoparticle/SINGLE domain folder.

## Read First

- `simple_nanoparticle.f90`: the nanoparticle type, automated atomic model building
- `simple_nanoparticle_utils.f90`
- `simple_nano_species.f90`: species discovery numerics (intensity mixture, class choice, threshold calibration, Gaussian width fit)
- `simple_atoms.f90`: atomic models, PDB/CIF I/O and density simulation
- `simple_calpha_finder.f90`: C-alpha candidates in density maps
- `simple_molecule_data.f90`
- `single_tseries_tracker.f90`, `single_tseries_extractor.f90`, `simple_tseries_graphene_subtr.f90`

## Connections

- Works closely with `commanders/single` and `ui/single`; `simple_commanders_atoms.f90` drives atoms, the C-alpha finder and nanoparticle analysis
- `simple_atoms` is not nanoparticle-specific: simulation, masking, the PCG atomic prior and tests use it too
- Touches `image/`, `volume/`, `motion/` (`simple_motion_align_nano`), and time-series workflows
- Tests: `simple_test_exec test=unit_single`; `lib_single` (pdb2mrc)
- Class-average selection for nanoparticle time series is the `nano` context of `model_cavgs_rejection`
  (`src/main/cavg_quality/`, lattice bank in `simple_cavg_quality_lattice.f90`); see
  `doc/policies/sieving_and_rejection/model_cavgs_rejection.md`. The classes of `analysis2D_nano` are not time chunks:
  a fresh 2D start draws a random previous class, so the `tseries` search never reads the chunk initialization

## Working Rule

Expect nanoparticle-specific assumptions here that do not generalize to the main cryo-EM path.
