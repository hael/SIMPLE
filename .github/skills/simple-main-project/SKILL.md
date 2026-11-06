---
name: simple-main-project
description: Use when working in SIMPLE's src/main/project subsystem, including the sp_project type and its mic/stk/ptcl/cls/out/optics segments, binary project I/O, segment writes and projinfo, merging distributed alignment documents, the out-segment artifact registry, the particle-layout digest of side files, reading projects of earlier releases, and project-level helpers (even/odd partition, state weight sets, particle and micrograph subsets, solve3D_addon supersets).
---

# `src/main/project`: the single-particle project (`sp_project`)

`sp_project` is one process's copy of a `.simple` file (segments `os_mic` .. `os_optics` 1-8,
`projinfo`, `jobproc`, `compenv` 11-13; bodies in submodules by segment). It owns what records mean,
their I/O, the part-document merge and the `out` registry; companions own one-maker state
(`simple_eo_partition`, `simple_state_weight_set`, `simple_ptcl_layout`, `simple_project_superset`).
When the file is written is the caller's decision: one owner per file, workers write part files
(`doc/policies/layering_and_project_ownership_policy.md` section 5). Format: `src/fileio/simple_binoris.f90`.

## Rules

Silent ones first: nothing reports these.

- A particle row persists only its named slots: `ptcl2D`/`ptcl3D` records are fixed reals (the
  `get_oriparam_ind` keys of `src/defs/simple_defs_ori.f90`). Any other key, and every string, set
  on a particle row lives in memory only and is gone after a write and read (skill `simple-main-ori`).
- `write`/`write_segment_inside` without a file name go to `projinfo`'s `projfile`. `read` sets it,
  `read_segment`/`read_non_data_segments` do not: after reading a copied project's non-data
  segments call `update_projinfo(projfile)`, or writes go to the source.
- `write_segment_inside` of an empty field leaves the file's old records: empty a segment with
  `clear_segment_inside`. On a missing file it writes the whole project.
- `write` deletes and recreates the file (`tempfile=.true.`: a `.tmp`, then a rename). No locking:
  another process can read it half-written. Full writes end a command-line step or a run; inside an
  iteration write each changed segment once.
- `read_segment(..., fromto)` fills rows `fromto(1):fromto(2)` at their index in the segment; the field keeps
  the segment's size and its other rows are fresh (defaults, 0 on `get`). Ranges are for `ptcl2D`, `ptcl3D` and
  `cls3D` only, in `.simple` and text files; another segment with a range stops.
- Index domains: `os_stk` `fromp`/`top`/`nptcls` count project rows; `indstk` is the image's index
  in its stack file (`1..nptcls_stk`). Active particles are `state > 0`, not `nptcls`. `ptcl2D` and
  `ptcl3D` hold the same rows (several getters stop on a count mismatch). `ptcl2D` `state` is the
  selection authority (`count_state_gt_zero` counts it alone), `ptcl3D` `state` the 3D label.
- `merge_algndocs(nptcls, ndocs, oritype, fbody)` reads `<fbody><i>.simple` (`i` padded to
  `numlen`), each a worker's range written by `binwrite_oritab`, into the in-memory field, and stops
  unless the ranges cover `1..nptcls` without overlap. By default it also writes the segment to the
  projinfo file: a master that owns its project passes `write=.false.` and writes once. A
  partition worker (`part` given, not a stream job) writes part files only: `write`,
  `write_segment_inside`, `clear_segment_inside` and `write_non_data_segments` stop on its input
  project (`worker_projfile_glob`, recorded by `params%new`). A stream job owns its job project.
- `ensure_field_halves(spproj, oritype, projfile)` is the halves' one maker. The master labels
  unlabelled rows by micrograph group (never relabels), the same in `ptcl2D` and `ptcl3D`, writes the
  changed segments and copies the labels into the loaded fields (stops if their rows differ from the
  file's); a worker (`l_distr_worker_glob`) only checks and stops when they are absent. Call it
  before the workers start. `oris%partition_eo` labels by row parity, without groups or provenance
  (only time series use it on particles).
- The `out` registry is keyed by `imgkind` (and `state` for per-state maps and FSCs): `add_*2os_out`
  replaces the entry with that key and stores an absolute path (a volume must exist); `get_vol`,
  `get_fsc` stop on none or several (`vol_msk`: empty). Find artifacts through getters, never by file
  name; drop an emptied state with `remove_state_artifacts_from_osout`.
- Per-particle side files (state weight set, sigma2 state, solve3D manifest, refine3D provenance)
  carry `ptcl_layout_digest` (projname, stack path and `indstk` per row) and are refused once the
  rows or the name change (sigma2 is then rebuilt). A project with new rows drops inherited
  registrations (`cavgs_work_project` deletes `sigma2_state`).
- Files of earlier releases must still read; a contract gap is repaired once by `fix_projfile`, not
  by a run-time fallback.

## Where to start

`simple_sp_project.f90` (the type, every interface), `_io`, `_core` (`merge_algndocs`,
`update_projinfo`), `_out`; a master owning its project: `simple_refine3D_strategy.f90`.

## Tests

- `build/bin/simple_test_exec test=unit_project` (about 2 s): `suite=project_records` (the part
  merge), `binoris`, `project_merge`, `eo_partition`, `state_weight_set` (both with earlier-release
  cases), `project_superset`, `particle_subset`, `micrograph_subset`.
- A format or merge change: `test=simulated_workflow suite=1jxy` (a distributed solve3D merges part
  documents); superset or side files: `test=solve3D_addon`; a refactor: `test=identity` (project
  files are compared record by record).

## See also

`doc/policies/project_stack_indexing_policy.md`, `doc/policies/3D/solve3D_addon_policy.md`
(superset); skills `simple-main-ori` (records), `simple-main-root` (builder), `simple-main-sigma2`.
