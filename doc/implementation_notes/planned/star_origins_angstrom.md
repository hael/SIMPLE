# STAR origins in Angstrom, and poses in export_relion

Status: plan, 2026-10-24. Nothing implemented. Found by the `wave_operators` run
(2026-10-23/24), which used specter as an independent reader of RELION STAR
files; see `ghostbuster_operators.md`, section 5.2.

## What is wrong

SIMPLE stores a particle's shift (`x`, `y` in the `ptcl2D` and `ptcl3D` fields)
in pixels of the particle's stack. RELION's `rlnOriginXAngst` and
`rlnOriginYAngst` are in Angstrom. Every STAR path copies one into the other
without the pixel size:

- import through the flag table of `simple_starproject.f90`
  (`initialise`: the 2D and 3D particle flags of `rlnOriginXAngst/YAngst`;
  `import_stardata`), reached by `import_starproject` and `import_particles`;
- export through the same flag table (`export_ptcls2D`, `export_ptcls3D`,
  `export_stardata`);
- `simple_starfile.f90`, `write_ptcl2D_table`, in its serial and its parallel
  writer;
- `simple_starproject_stream.f90`, `starfile_set_particles2D_table` (the
  stream's export).

Measured through specter at 1.3 A per pixel (two blobs, ten random poses):
imported as is, the blobs land 1.50 pixels off; with the shift divided by the
pixel size, 0.019 pixels. The sign is right, and SIMPLE's Euler angles map to
RELION's exactly (`e1` = `rlnAngleRot`, `e2` = `rlnAngleTilt`, `e3` =
`rlnAnglePsi`).

Two further findings of the same check:

- `export_relion` (`simple_relion.f90`, `write_particles2D_star`) writes no
  poses and no origins: coordinates and CTF only.
- A STAR file without a `data_optics` block (the RELION 3.0 layout, which
  specter writes) stops the import with "oritab empty", a message that does not
  say what is missing.

## Fix

1. **One conversion.** Angstrom = pixels times the pixel size of the particle's
   stack (`os_stk` `smpd`, through the particle's `stkind`). Import divides,
   export multiplies, on every path above and nowhere else. A particle whose
   stack has no positive pixel size stops the import or export with a message
   naming the particle and the stack.
2. **Import.** The flag table says which fields carry Angstrom in STAR and
   pixels in SIMPLE (an attribute of `star_flag`, not a check on the flag's
   name). The conversion happens once per particle, at the point where the
   stack's pixel size is final, on every import path (`import_ptcls2D`,
   `import_ptcls3D`, `import_ptcls`). Where that point is (after
   `check_stk_params` fills the stack's pixel size from the particles or the
   optics block) is settled by reading the code; the contract is one
   conversion per particle per import.
3. **Export.** The flag-table export, both writers of `write_ptcl2D_table` and
   the stream's writer multiply by the stack's pixel size.
4. **export_relion writes poses.** `rlnAngleRot`, `rlnAngleTilt`,
   `rlnAnglePsi`, `rlnOriginXAngst` and `rlnOriginYAngst` from the `ptcl3D`
   field when it holds 3D poses, otherwise `rlnAnglePsi` and the origins from
   `ptcl2D`. Which field holds poses is decided with the project's existing
   accessors, not by testing values for zero.
5. **RELION 3.0 layout.** Refused with a message that names the missing
   `data_optics` block and says that RELION 3.1 or later STAR files are
   required. No reader for the 3.0 layout (release 4 keeps no compatibility
   layers).

The stream's export lives in the stream subsystem; the change there is the
unit conversion only.

## Tests

- In the STAR testers (`unit_project`: "STAR file", "STAR project", "STAR
  stream export"): at a pixel size other than 1 (1.3), pixels to STAR to pixels
  is the identity to float precision on every path; a STAR file with known
  Angstrom origins imports to the known pixel shifts; `export_relion` writes
  the poses and origins of a project with 3D poses and of one with 2D poses
  only; a STAR file without an optics block is refused with the new message.
- External check on the Dell: the specter bridge of the `wave_operators` run
  (`~/agent_runs/wave_operators/scratch/keep/oracle/bridge_positions.py`, the
  specter environment beside it; read only, copy what is needed). Blob
  positions SIMPLE to STAR to specter and specter to STAR to SIMPLE agree within
  0.05 pixel with no conversion by hand, on the paths that write and read
  particle origins.
- The fast gate and `check_test_registry.py`.

## Execution

0. **Fix and tests.** Items 1 to 5 with the tests above. Exit: the testers and
   the external check pass; the fast gate passes; every path in "What is wrong"
   is accounted for in the progress section with the test that covers it.
1. **Documents and report.** Where `doc/user_guide/user_guide.md` describes
   `import_particles`, `import_starproject`, `export_starproject` and
   `export_relion`, it states the units and what is exported; this plan moved
   to `completed/` with the report beside it.

### File table

| File | Purpose | Phases |
| --- | --- | --- |
| `doc/implementation_notes/planned/star_origins_angstrom.md` | this plan: progress | each |
| `src/main/star/simple_starproject_utils.f90` | the unit attribute of `star_flag` | 0 |
| `src/main/star/simple_starproject.f90` | import and export through the flag table | 0 |
| `src/main/star/simple_starfile.f90` | `write_ptcl2D_table`, both writers | 0 |
| `src/main/star/simple_starproject_stream.f90` | the stream's particle export | 0 |
| `src/main/star/simple_relion.f90` | `export_relion` poses and origins | 0 |
| `src/main/commanders/simple/simple_commanders_starproject.f90` | only if the conversion point needs it | 0 |
| `src/main/commanders/simple/simple_commanders_relion.f90` | only if `export_relion` needs it | 0 |
| `doc/user_guide/user_guide.md` | units and exported fields of the four programs | 1 |
| `doc/implementation_notes/completed/star_origins_angstrom.md` | this plan, moved when done | 1 |
| `doc/implementation_notes/completed/star_origins_angstrom_report.md` | the run's report | 1 |

## Progress

(none yet)
