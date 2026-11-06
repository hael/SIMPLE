---
name: simple-main-nano
description: Use when working in SIMPLE's src/main/nano subsystem (SINGLE), including the nanoparticle type and species-free atom detection (detect_atoms, atoms_stats), the atomic density model with its measured transfer, species numerics and the detect_atomic_species call, the atoms type (PDB/CIF I/O, pseudo-atoms, pdb2mrc density simulation) shared with SIMPLE simulation and the PCG atomic prior, the C-alpha finder, the test molecules, and the time-series tracker and graphene subtraction of liquid-cell workflows.
---

# `src/main/nano`: nanoparticle atomic models and the SINGLE domain objects

Owns `nanoparticle` (detection, the density-model fit, the `atoms_stats` analyses), the density model
(`simple_nano_density`), the species numerics and call, the `atoms` type, the C-alpha finder, the test
molecules and two time-series singletons. Nothing here touches an `sp_project`. Not here: the programs
(`commanders/simple/simple_commanders_atoms.f90`, `commanders/single/`), the PCG atomic prior
(`strategies/parallelization/simple_rec3D_pcg_strategy.f90`), the nano movie aligner (`motion/`), the
`tseries` 2D search, time-series class-average rejection (`cavg_quality`), `trajectory_denoise` (`flex/`),
the element and lattice tables (`src/defs/simple_defs_atoms.f90`).

## Rules

Silent ones first: nothing reports these.

- Units: `atoms` holds Å; `nanoparticle` and `nano_density` hold 1-based voxels. Å = (v - 1) smpd and
  v = Å / smpd + 1 (`set_atomic_coords`, `write_centers`, `call_atomic_species`); `convolve` puts voxel j
  at (j - 1) smpd. `set_atomic_coords` from an array takes voxels as given. Widths B in Å². The box
  centre is ldim/2 smpd from the corner (`center_inbox`, `center_pdbcoord`).
- The PDB B-factor column is overloaded. `_ATMS.pdb` of `detect_atoms` holds the fitted width B, occupancy
  I/S and chain A (strong) or B (weak atoms); other `write_centers` outputs put `valid_corr`, `cn_gen`, ... there.
  `_ATMS.pdb` as written, weak atoms included, is a valid `atoms_stats` input (they get no component).
- The products of `detect_atoms`, named after the map (`_ATMS.pdb`, `_atoms.csv`, `_model.txt`,
  `_transfer.tsv`, `_SIM`/`_FIT`/`_RES`/`_BIN`/`_MSK`/`_CC.mrc`), are a contract: `detect_atomic_species`
  pairs PDB row i with csv row i (count checked, order not); `atoms_stats` pairs it with component i of
  `vol2` (`_CC.mrc`); `autorefine3D_nano` copies them by name; the PCG atomic prior reads `BFAC`
  (`_atoms.csv`) or `POSTERIOR_INTENSITY`/`OCCUPANCY` (`_species.csv`) by column name and, without `BFAC`,
  falls back to a position-only prior at one width. Rename in all of them.
- Detection is species-free: `detect_atoms` refuses `element`, `nspecies`, `vol_even`, `vol_odd` and writes
  pseudo-atoms `X1` (`X1`-`X3` are Z 901-903; `convolve` renders them as occupancy x K(beta), no `lp` blur).
  The element serves `atoms_stats` (`fillin_atominfo` stops without one); `get_lattice_params` stops on an
  element its table lacks.
- Thread-count identity: `convolve` gives each z-plane one writer and adds atoms in input order (`pdb2mrc`
  is byte-identical at 1 and 3 threads); the density fit colours its cells so the result depends on the
  team only through the halo FFTs (capped at `FFT_NTHR_MAX`). Keep both.
- `init_graphene_subtr(params)` runs on the caller's parameters (the cc objective; the calculator points
  at them until `kill_graphene_subtr`): an initialiser leaves the process's threads and RNG alone.
- `single_tseries_tracker` and `simple_tseries_graphene_subtr` keep module state: `init_*` -> use ->
  `kill_*`, one at a time. Legacy singletons: do not copy the pattern; convert to a class when
  touched. `init_tracker` keeps pointers to `params` and the frame-name arrays (they must outlive
  `kill_tracker`) and needs one frame per file (THROW_HARD).
- `nanoparticle`: `new(params, map[, msk])` (`msk` a radius in voxels; without it the mask is estimated
  and a map with no particle stops) -> `identify_atomic_pos` (`detect_atoms`), or `set_atomic_coords`,
  `set_img(cc, 'img_cc')`, `update_ncc`, `fillin_atominfo`, `write_csv_files` (`atoms_stats`) -> `kill`.

## Where to start

- `simple_nanoparticle.f90`: `identify_atomic_pos` (detection), `fit_density_model` (the model and its
  products), `fillin_atominfo` (statistics). `simple_nano_density.f90`: its header gives the model and
  the fit order (`doc/implementation_notes/completed/atomic_density_model.md`).
- `simple_atomic_species.f90` (`call_atomic_species`) on the image-free `simple_nano_species.f90`.
- `simple_atoms.f90` serves all of SIMPLE (simulation, masks, volume ops, `element_exists` in the
  parameter parser, the PCG prior). `convolve`/`pdb2mrc` make the truth maps of the state fixture and
  the simulated workflows; `simple_molecule_data` reads `testdata/pdb` under `SIMPLE_SOURCE_PATH`.

## Tests

- `build/bin/simple_test_exec test=unit_single` (about 1 s; `suite=` `atoms`, `cif2mrc`, `c_alpha_finder`,
  `species`, `density`); `test=lib_single` (nightly, 20 s; `suite=` `pdb2mrc`, `density_particles`).
- High-level, by hand with their command line in `commanders/test/simple_test_registry.f90`:
  `single_atoms_stats` and `_cdsew`/`_cdsez`/`_cdser` (10-46 s), `species_discovery` (1-2 min),
  `single_workflow_fcc`/`_wurtzite` (11-14 min: `trajectory_solve3D`, `autorefine3D_nano` with the
  atomic prior). A `convolve` change also moves `state_fixture` and the simulated workflows.

## See also

`doc/implementation_notes/completed/species_discovery_completed.md`, `.../completed/trajectory_solve3D.md`;
skills `simple-main-motion`, `simple-cluster-cavgs-quality` (`nano` context), `simple-main-flex`.
