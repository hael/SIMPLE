---
name: simple-main-star
description: Use when working in SIMPLE's src/main/star subsystem, including STAR file reading and writing, the Fortran wrappers of the C++ starfile library, Relion interoperability (export_relion, RELION 3.1 layout), starproject import/export utilities (import_starproject, export_starproject, the flag table, optics groups from tilt info), unit conversions between SIMPLE and STAR (origins in Angstrom, defocus, phase shift, coordinates), the parallel STAR writers used by sp_project, and the stream's STAR export.
---

# `src/main/star`: STAR files and RELION interoperability

Owns the translation between `sp_project` segments and STAR files: import through the flag table
of `starproject` (`import_mics`, `import_ptcls2D`, `import_ptcls3D`, `import_ptcls`, `import_cls2D`,
`finalise_ptcl_import`), the writers `starfile` (optics, micrograph, particle tables, used by
`sp_project%write_mics_star`/`write_ptcl2D_star`), the stream's export `starproject_stream` (p01,
p02, p04, the pool), `relion_project` (`export_relion`), the batch optics grouping from tilt info
(`starproject%assign_optics`) and the bindings `simple_starfile_wrappers` to
`src/extlibs/starfile` (also used by `motion/` and `ctf/`). Not owned: the commanders
(`simple_commanders_starproject`, `simple_commanders_relion`), the project's fields (`project/`),
the stream's optics groups (`stream/shared/simple_optics_groups`). Keep STAR and RELION
conventions inside this folder.

## Rules

Silent ones first: nothing reports these.

- SIMPLE's `x`, `y` are pixels of the particle's stack; STAR's `rlnOriginX/YAngst` are Angstrom.
  `import_ptcls2D`/`import_ptcls3D` leave origins in Angstrom until `finalise_ptcl_import`, which
  divides both `ptcl2D` and `ptcl3D` by each stack's `smpd` once and clears the mark. Order:
  import one field, copy it into the other, finalise, once per import: a second finalise, a new
  particle import before the finalise, and `export_mics`/`export_cls2D`/`assign_optics` of an
  unfinalised import stop (`ptcl_import_open`). A project written without the finalise is not
  caught. `import_ptcls` finalises itself. The finalise moves the optics constants (`kv`, `cs`,
  `fraca`, `smpd`, `detpix`, `mag`) to the stacks and strips them from the particle rows, which have no
  slot for them.
- Every particle writer multiplies by the stack's `smpd`: `write_ptcl2D_table` (serial and
  parallel), `starproject_stream`'s particle table, `relion_project%write_particles2D_star`. A
  change of a particle column or unit is made in all of them and in the flag table (`initialise`).
- Other units: `dfx`/`dfy` microns in SIMPLE, Angstrom in STAR; `phshift` radians, STAR degrees;
  `e1`, `e2`, `e3` are `rlnAngleRot`, `rlnAngleTilt`, `rlnAnglePsi`. Writers add box/2 to
  `xpos`/`ypos` (SIMPLE's, EMAN's, box corner; STAR's particle centre), and `finalise_ptcl_import`
  takes the stack's box/2 off once (`coords_to_corners`, the same mark-and-clear as the origins).
- Optics rows are found by `ogid`, never by row (stream ids keep gaps); snapshot writers offset
  `ogid` by `optics_offset`.
- Writers take each particle's `stkind` and `indstk` as stored (`write_ptcl2D_star` maps them
  first). Micrograph names are written only when stacks and micrographs are equal in number, and
  then row i of each is taken to be the same micrograph.
- Every writer writes `<name>.tmp` and renames it when complete (`simple_rename` replaces the target
  atomically, so the file is never absent) (`starfile`, `starproject_stream`,
  `export_mics`, `export_cls2D`, `export_opticsgroups`, the `relion_project` writers);
  `export_stardata` stops on a file name without `.tmp`. A new writer does the same. An export does
  not change its project: `export_cls2D` adds its `stk` and `ncls` to a copy of the `cls2D` rows
  (`export_mics` still numbers the micrograph rows it is given: `nmics`, `micid`).
- Batch `assign_optics` groups hierarchically, the stream's p02 by single linkage: ruled to move
  batch to the stream's grouping (after comparing both on one data set with beam tilt). Do not
  extend the hierarchical one.
- Parallel writers: one table per thread (private), rows split by `omp_get_max_threads()`, every
  string read or set inside `!$omp critical`, tables flushed in `!$omp ordered`. A new string column
  goes inside the critical block.
- A wrapped table lives from `starfile_table__new` to `starfile_table__delete`; files open with
  `open_ofile` and close with `close_ofile`.
- A STAR file without `data_optics` (RELION 3.0) is refused (`read_layout`); callers `THROW_HARD`.
  A particle with a shift and no stack pixel size stops import or export (`THROW_HARD`).

## Where to start

- `simple_starproject.f90`: the flag table (`initialise`), imports, `finalise_ptcl_import`,
  `assign_optics`; `simple_starproject_utils.f90`: `star_flag`, line splitting, path resolution.
- `simple_starfile.f90`: the writers; `simple_starproject_stream.f90`: the stream's export;
  `simple_relion.f90`: `export_relion` (poses from `ptcl3D` when valid, else `ptcl2D`).
- Callers that show the order: `commanders/simple/simple_commanders_starproject.f90`.

## Tests

- `build/bin/simple_test_exec test=unit_project` (2-3 s; `suite=star_file`, `star_project`,
  `star_stream_export`): writers, both parallel and serial, origins in Angstrom, round trips, the
  RELION 3.0 refusal, `export_relion` poses. No high-level test imports or exports STAR files; a
  unit change needs a round trip at a pixel size other than 1 in all three sub-suites.

## See also

`doc/implementation_notes/completed/star_origins_angstrom.md` (the conversion contract),
`doc/policies/phase_shift_ctf_policy.md`; skills `simple-main-project`, `simple-main-stream`.
