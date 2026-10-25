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
| `src/main/star/simple_starfile_tester.f90` | "STAR file" sub-suite: origins of both writers | 0 |
| `src/main/star/simple_starproject_tester.f90` | "STAR project" sub-suite: import paths, `export_relion`, the 3.0 refusal | 0 |
| `src/main/star/simple_starproject_stream_tester.f90` | "STAR stream export" sub-suite: origins | 0 |
| `doc/user_guide/user_guide.md` | units and exported fields of the four programs | 1 |
| `doc/implementation_notes/completed/star_origins_angstrom.md` | this plan, moved when done | 1 |
| `doc/implementation_notes/completed/star_origins_angstrom_report.md` | the run's report | 1 |

## Progress

### Phase 0 (run star_wave, Phase 0), 2026-10-24: fix and tests done

**In plain words.** SIMPLE now writes a particle's shift to RELION STAR files in
Angstrom and reads it back into pixels, using the pixel size of the stack the
particle's image lives in, exactly once, on every path that carries origins.
`export_relion` writes the particle poses. A STAR file without an optics block
(the RELION 3.0 layout) is refused with a message saying what is missing. An
independent program, specter, agrees with SIMPLE's particle positions to within
0.019 pixel in both directions without any conversion by hand; before the fix
the same check put them 1.50 pixels apart.

**Files changed.** `src/main/star/simple_starproject_utils.f90`,
`simple_starproject.f90`, `simple_starfile.f90`, `simple_starproject_stream.f90`,
`simple_relion.f90`; `src/main/commanders/simple/simple_commanders_starproject.f90`
(the renamed finalising call only); the testers `simple_starfile_tester.f90`,
`simple_starproject_tester.f90`, `simple_starproject_stream_tester.f90` (added
to the file table); this plan and the run plan.

**What was done.**

- The flag table marks `rlnOriginXAngst`/`rlnOriginYAngst` with a new `angst`
  attribute of `star_flag` (Angstrom in STAR, pixels of the particle's stack in
  SIMPLE). Reading such a column marks the imported origins as still in
  Angstrom.
- Import converts in `finalise_ptcl_import` (formerly `check_stk_params`), the
  step every particle import already ends with, after it has completed the
  stacks' pixel sizes and after the imported field has been copied into the
  other particle field: both fields are divided by the pixel size of each
  particle's stack, and the mark is cleared, so a second call cannot convert
  again. A particle with a nonzero shift whose stack has no positive pixel
  size stops the import with a message naming the particle and the stack.
- Export multiplies in both writers of `write_ptcl2D_table` (the serial one,
  used by `export_starproject` and the stream's project snapshots, and the
  parallel one, which checks all pixel sizes before its parallel region), in
  the stream's `starfile_set_particles2D_table`, and in `export_relion`.
- `export_relion` writes `rlnAngleRot`, `rlnAngleTilt`, `rlnAnglePsi` and the
  origins from the 3D field when `has_valid_ptcl3D_alignment` says it holds
  poses, otherwise `rlnAnglePsi` and the origins from the 2D field.
- The RELION 3.0 layout: the reader no longer treats an unnamed `data_` block
  as particles. `read_layout` reads a file's block layout and returns the
  refusal ("<file> has no data_optics block: RELION 3.1 or later STAR files are
  required") without stopping; the micrograph and particle imports stop with
  it. `import_cls2D` keeps the plain reader, because a RELION model file has no
  optics block.

**Decisions taken during the phase.**

- The flag-table export of particles (`export_ptcls2D`, `export_ptcls3D`) had
  no caller anywhere, not even a test. It was deleted as dead code rather than
  fixed, with what only it used (`center_boxes`, `propagate_optics_box`,
  `propagate_optics2D`/`3D`, the image-name branch of `export_stardata`), and
  so were two other uncalled exports of the same module (`export_iter3D`,
  `export_stream2D`, with the 3D-class flag table) and unused parts of the
  utilities (`get_rlnflagindex`, `stk_map`, two unused `star_flag` fields).
  No flag-table path writes particle origins any more, so item 3's flag-table
  export is met by there being none.
- A particle record always reports `x` and `y` as present (zero is a valid
  shift), so the pixel-size check stops a run only when a nonzero shift has to
  be converted; a zero shift is zero in either unit.
- Unit tests may not start child processes and cannot observe a hard stop, so
  the refusal is tested through `read_layout`, the function whose message the
  imports stop with.
- `initialise` never set `initialised` (so a second call appended the flag
  table again) and `kill` left the block layout of the last file behind; both
  are fixed, since `read_layout` relies on them.
- `export_relion` scanned up to 1000 particles past the first for nonzero
  coordinates without stopping at the last particle; on a project without
  coordinates it read past the array (the base build crashes with a
  segmentation fault on the external check's project). The scan now stops at
  the last particle.
- The parallel writer, used only by its tests, was kept and fixed because this
  plan names it.

**Every path that reads or writes particle origins, and its test.**

| Path | What happens now | Covered by |
| --- | --- | --- |
| `import_starproject` (all particle branches, `import_ptcls2D`/`import_ptcls3D` then `finalise_ptcl_import`) | Angstrom divided by the stack's pixel size, both fields, once | "STAR project": `test_import_origins_angstrom` (ptcl2D, ptcl3D), `test_roundtrip_origins`; external check, import_type ptcl3D and ptcl2D |
| `import_particles` (`import_ptcls`) | the same | "STAR project": `test_import_origins_angstrom`; external check |
| `export_starproject` and the stream's snapshots (`write_ptcl2D_table`, serial) | pixels times the stack's pixel size | "STAR file": `test_write_ptcl2D_origins_angstrom`; "STAR project": `test_roundtrip_origins`; external check, export_starproject |
| `write_ptcl2D_table_parallel` | the same, equal to the serial writer | "STAR file": `test_write_ptcl2D_origins_angstrom` |
| stream export (`starfile_set_particles2D_table`) | the same | "STAR stream export": `test_particles2D_origins` |
| `export_relion` | poses and origins written, 3D or 2D field | "STAR project": `test_relion_export_poses`; external check, export_relion (3D poses) |
| flag-table particle export (`export_ptcls2D`/`3D`) | deleted, no caller | none needed |
| `scripts/relion2simplegui.pl` | not changed: reads RELION 3.0 `_rlnOriginX` (pixels) into `x`, consistent in units but a reader of the 3.0 layout; outside this phase's files, left to the maintainer | none |

The origin tests use two stacks at different pixel sizes (1.3 and 0.9 Angstrom
per pixel), so taking the pixel size of the wrong stack fails them. Two
deliberate breakages were run and reverted: without the conversion, 9 checks
of the new tests fail; with stack 1's pixel size used for every particle, 5
fail (logs `scratch/t_unit/unit_project_mutant.log` and
`unit_project_mutant2.log`).

**Evidence** (paths under `/home/elmlundho/agent_runs/star_wave`).

- `unit_project`: all 2914 checks pass on the debug build,
  `scratch/t_unit/unit_project.log`.
- External check through specter, `scratch/keep/bridge/bridge_positions.py`
  (adapted from the previous run's bridge; run with that run's specter
  environment, read only), log `scratch/bridge_run4.log`, report
  `scratch/keep/bridge/positions/report.json`. Two blobs, ten random poses,
  1.3 Angstrom per pixel, box 64. Largest blob position error:
  specter to SIMPLE through `import_particles`, `import_starproject`
  import_type=ptcl3D and import_type=ptcl2D: 0.019 pixel each;
  SIMPLE to specter through `export_relion` (3D poses): 0.012 pixel, and
  through `export_starproject` (in-plane poses): 0.012 pixel. All below 0.05.
  The same script on the base build (`scratch/bridge_base.log`): 1.50 pixels
  on all three imports, and `export_relion` crashes. specter's STAR writer uses
  the 3.0 layout, so its files pass through the previous run's `relion31_star`
  conversion first. specter's reader needs a defocus column and reads no angle
  unless rot, tilt and psi are all present; the check adds those absent
  columns to a copy with RELION's default of zero, leaving every value SIMPLE
  wrote unchanged. `export_starproject` writes nothing for a project without
  micrographs, so the check's project imports one blank micrograph and assigns
  its optics group.
- Base against run on the same project (`export_starproject`,
  `scratch/cmp_star`): the micrograph file is identical; in the particle file
  only the two origin columns differ, by a factor of 1.3000 (1.29997 at the
  file's printed precision).
- `scripts/check_test_registry.py` passes (`scratch/check_registry.log`).
- The fast gate (`scratch/fast_gate2.log`): 14 of 15 suites pass;
  `unit_pftc_align2D3D` fails one check ("a class tie draws one of the tied
  classes", 2D probability table input and output). It fails identically on
  the unmodified base, both the debug build (`scratch/build_debug_initial.log`)
  and the Release build (`scratch/t_base/base_pftc.log`); it is outside the
  STAR code and was found the same way by the previous run. Not fixed here.
- No file mode changes; no `THROW_HARD` argument continued across lines.

**Exit items.** The testers pass; the external check passes on every path that
writes or reads particle origins and can be reached from the command line (the
stream export is covered by its unit test); the fast gate passes except the
failure that the base shows identically; every path of "What is wrong" is in
the table above with its test.

**Left for the maintainer.** `scripts/relion2simplegui.pl` (above). In
`simple_commanders_starproject.f90`, the automatic import tests for
`run_data.star` but then reads `run_it<NNN>run_data.star`, a file name that
cannot exist; this was already so and is outside the units fix. The
`unit_pftc_align2D3D` failure of the base.

### Phase 1 (run star_wave, Phase 3), 2026-10-24: documents and report done

`doc/user_guide/user_guide.md` has a new section, "Exchanging particles with
RELION (STAR files)": the units on each side (shifts in pixels of the
particle's stack against Angstrom origins, defocus, phase shift), the refusal
of the RELION 3.0 layout, and what `import_particles`, `import_starproject`,
`export_starproject`, `export_relion` and the stream's export read or write.
The guide had not described these programs before, so the section is new
rather than amended. This plan moved to `doc/implementation_notes/completed/`
with the run's plan; the report of the whole run, covering this plan, is
`doc/implementation_notes/completed/star_origins_and_wave_simulate_run_report.md`
(the run plan named one report for both plans instead of this plan's own
`star_origins_angstrom_report.md`).
