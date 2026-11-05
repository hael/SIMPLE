# Stream stages, their folders and the files they share

Reference for `simple-main-stream`: what each process owns and hands on. The contracts are in
`doc/policies/stream/` (index and per-stage diagrams in its README); this page is the map.

## Processes

The master (p00, `simple_stream prg=master`, commander `simple_commanders_stream_p00_master`)
forks one process per stage under the run folder, in the order of `simple_stream_master_stage_ids`.
Each stage is a thin commander in `src/main/commanders/stream` (normalise the command line in
`set_<stage>_cline`, then `new` -> `iterate` until `finished` -> `finalize` -> `kill`) driving a
stage type in `src/main/stream/stages`. A stage runs in its own folder (`params%new` with
`mkdir=yes` moves there) and can run on its own (`simple_stream prg=` `preproc`, `assign_optics`,
`gen_pickrefs`, `pick_extract`, `sieve_cavgs`, `pool2D`, `solve3D_stream`; `dir_exec` for p01, p04,
p06, p07).

| Stage | Type | Folder | Takes from (`dir_target`) | Hands on | GUI updates read |
|---|---|---|---|---|---|
| p01 preprocessing | `stream_stage_preprocess` | `preprocessing/` | movies (`dir_movies`) | job sets renamed into `spprojs_completed/`; `STREAM_IDLE`, `STREAM_FINISHED` | CTF resolution, astigmatism, ice (kept in `gui_thresholds.txt`) |
| p02 optics | `stream_stage_optics` | `optics_assignment/` | p01 completed | `optics_map_<id>.simple`, then `optics_map_<id>.txt` | none |
| p03 initial analysis | `stream_stage_initial_analysis` | `initial_analysis/` | p01 completed | `selected_references.mrcs` (once, `publish_pickrefs`) | picking-reference selection and cycle |
| p04 reference picking | `stream_stage_refpick` | `reference_based_picking/` | p01 completed, p03's references | job sets into `spprojs_completed/`; `moldiam.txt` after the templates; markers | none |
| p05 particle sieving | `stream_stage_sieve` (drives `ptcl_sieve`) | `particle_sieving/` | p04 completed, `moldiam.txt`, p04's markers | chunk copies and `sieve_final_c<n>_f<m>.simple` in `spprojs_completed/` | none |
| p06 pool 2D | `stream_stage_pool2D` (holds `stream_pool2D`) | `classification_2D/` | p05 completed | publications `<id>.simple` with `<id>_cavgs*.mrcs`, `<id>_frcs.bin` (newest 2 kept); `snapshots/`; final project | 2D mask diameter, 2D snapshots |
| p07 multistate 3D | `stream_stage_solve3D` | `solve3D_multistate/` | p06 completed | its project (temporary and rename), `state_refs/`, `assign/it_<n>/`, `snapshots/` | 3D snapshots |

p04, p05 and p06 also read the newest optics map from `optics_assignment/`
(`simple_optics_maps`: readers find a map by its `.txt`, apply it by import index; a micrograph the
map does not list gets group 1).

## Who writes what

- Each stage writes only in its own folder; downstream reads its completed folder. Nothing locks:
  sequencing (temporary name, rename, settle time in `stream_watcher`) is the protocol.
- Exceptions (`doc/policies/layering_and_project_ownership_policy.md`, section 8): the pool and
  its `refine2D_distr` master take turns on `refine2D.simple` and `sigma2_state.bin` (the pool
  writes only while `available()`; the master touches `REFINE2D_FINISHED` when done); the pool's
  final write puts `TERM_STREAM` in the folder of the iteration it cancelled (`async=yes` masters
  exit on it); p01/p04 rewrite the `mic` segment of a finished job's set before moving it; the
  master touches `STREAM_FINISHED` in a linked earlier preprocessing folder (`dir_preprocess`).
- Queued jobs write their record (`<exit file>.job`) and exit status beside them
  (`simple_qsys_job_record`, `qsys_async_job`); a stage's `finalize` cancels its jobs, and a
  restart sets aside a folder whose job recorded itself and wrote no status (`fresh_job_dir`,
  `<folder>_unfinished<k>`).
- What survives a restart (markers, records, numbering): `doc/policies/stream/restart_policy.md`,
  section 4; each stage's module header has a RESTART paragraph.

## The master's parts (`src/main/stream/master`)

- `simple_stream_master_stage`: one stage's fork, command line and master-side pipes; `start`,
  `skip`, `request_stop`, `force_stop`, `discard_pipes`, `send_update` (only when changed).
- `simple_stream_master_stage_ids`: ids, job names, GUI keys, pipe arrays of `simple_stream_state`.
- `simple_stream_master_meta_store`: the latest message per type; `clear_stage` on a restart.
- `simple_stream_master_gui_commands`: parses the GUI's answer; drops oversized selections.
- `simple_stream_master_resources`: threads and parts of every stage and its jobs, with
  `SIMPLE_STREAM_<SET>_*` overrides (README of `doc/policies/stream/`).

## Shared modules (`src/main/stream/shared`) and support

- `simple_stream_pipe` (framing), `simple_stream_state` (pipe descriptors), `simple_stream_sigterm`.
- `simple_stream_job_sets`: p01's and p04's numbered sets (write, submit, collect, complete, restore).
- `simple_optics_groups` (single-linkage grouping, ids kept across passes), `simple_optics_maps`
  (map files, `copy_project_with_optics_map`, `latest_optics_map_table`).
- `simple_mic_import` (append finished upstream `mic` segments), `simple_mic_selection` (the one
  threshold rule; rejection is one-way).
- `simple_stream_gui_senders`, `simple_stream_meta_plots`: GUI messages built from domain data.
- `src/main/stream`: `simple_stream_watcher` (movie and project watcher, history, settle time),
  `simple_stream_utils` (`upstream_done`, `init_stream_qenv`, `import_new_projects`, `log_rss`),
  `simple_mini_stream_utils` (used by the batch `preprocess`, `validate` and the `mini_stream` test,
  not by the stages).

## The pool (`src/main/stream/pool2D`)

`stream_pool2D` keeps its state private; p06 hands it sets (`append_sets`), reads `stats()`, and
calls `iterate`, `update_status`/`update`, `publish`, `write_snapshot`, `finalise`. Each iteration
is one queued script: `calc_pspec`, then `refine2D_distr` through `simple_private_exec` (the batch
`commander_refine2D`, workers `prg=refine2D`) with `stream2d=yes`, `maxits=1` from `startit`, the
working `box_crop`/`smpd_crop`/`msk_crop`. `simple_stream_refine2D_utils` holds the stateless
helpers (clean-up, iteration files, publications, class draws). Schedule, sampling, dimensions,
snapshots: `doc/policies/stream/pool2D_policy.md`.

## p07

Reference volumes are made once from the first publication (solve2D into 150 classes, the pool
model's selection, `cls_expansion` to about 500 sub-classes, `solve3D_cavgs`), checked, retried,
kept in `state_refs/`; each later publication is merged into rows that only grow, split by
`cls_expansion`, matched against the references (`assign_cavgs_states`, `simple_cavgs_state_match`
in `src/main/class`), rejected by fit, mapped to its particles and refined locally with states fixed
(`refine3D refine=cont`) in `assign/it_<n>/`. The method (rules, check, mapping, seeds, command
lines) is `src/main/solve/simple_state_assignment.f90`, shared with `solve3D_pyramid`.
`doc/policies/stream/stream_3D_ingestion_policy.md` sections 4, 9 and 10.
