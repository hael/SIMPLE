---
name: simple-main-stream
description: Use when working on SIMPLE's streaming pipeline, including the p00-p07 commanders in src/main/commanders/stream, the stage types in src/main/stream/stages, the master's parts in src/main/stream/master, the shared stage modules in src/main/stream/shared (pipes, sigterm, job sets, GUI senders, micrograph and optics helpers), the 2D pool in src/main/stream/pool2D, the src/main/stream support layers (watchers, particle-sieve integration, streaming variants of preprocessing and refine2D workflows), stream restart, IPC and process/file ownership.
---

# `src/main/stream`, `src/main/commanders/stream`: the streaming pipeline

The master (p00) forks one process per stage and relays the GUI; each stage (p01-p07) is a type in
`stages/` driven by a thin commander, and owns its folder; `pool2D/` is p06's pool. Stages hand
data on as files, GUI messages go over pipes. Not owned here: the particle sieve (skill
`simple-main-sieve`), the batch programs the stages queue (`refine2D_distr`, `solve2D`,
`solve3D_cavgs`, `refine3D`, `make_pickrefs`, `cls_expansion`, `assign_cavgs_states`), p07's method
(`src/main/solve/simple_state_assignment.f90`), GUI metadata types (`src/utils/gui`), job records
(`src/utils/qsys`), STAR writers (`src/main/star`). Stages, folders, files:
`references/stages_and_files.md`.

## Rules

Silent ones first: nothing reports these.

- Hand files on by temporary name and rename in the destination folder (an optics map: `.simple`
  first, `.txt` last); readers take a file once untouched for the watcher's settle time. Keep both.
- One writer per file: a stage writes only in its own folder (exceptions in the reference); the
  pool writes `refine2D.simple` only while `available()`.
- Stages run in their folder (`params%new`, `mkdir=yes`); relative paths resolve there.
  `projinfo`'s `projfile` may be a bare name: pass a path to `write`/`write_segment_inside`.
- Restart state is on disk: numbering (sets, optics maps, publications) continues after the highest
  on disk, restored inputs go into the watcher's history. A new state file goes in
  `restart_policy.md` section 4 and must survive clean-ups: p06's `cleanup_root_folder` deletes
  every `.mrc .mrcs .txt .star .eps .jpg .jpeg .dat .bin` and chunk folder in its folder.
- Look optics rows up by `ogid`, never by row. p07 matches particles by stack name and `indstk`,
  never by row; its rows only grow.
- The pool may run downscaled: rescale what goes downstream or label its `smpd`. Change its mask
  only through `set_mskdiam` (resets `msk_crop`, which the workers take over `mskdiam`).
- p06's pass is refresh, publish, import, dispatch: a publication is one completed iteration. Only
  the pool sets `stream2d=yes`; batch refine2D changes reach the pool.
- One publisher of the picking references (`publish_pickrefs`, refuses once published); p03
  applies GUI updates before its cycle steps.
- Master: read or discard a stage pipe on the main thread only under `meta_mutex`; fork only before
  the listener starts or under that lock (a child copies no threads and hangs on a held lock).
- A message is a byte copy (`transfer`): nothing allocatable sent (override `serialise`), arrays
  cut to their `MAX_*` capacity (setters `THROW_HARD`). A new GUI answer key: parse it in
  `stream_master_gui_commands`, add it to `ipc_policy.md` section 6.
- Signal handlers only set a flag; poll `sigterm_received()` between steps. `finalize` cancels the
  stage's jobs before `kill`. `forked_process` never restarts a child, and a forked stage exits
  there, never running the tail of `production/simple_stream.f90`.
- GUI domain work stays in stream code, not `src/utils/gui/metadata`. Threads and parts come from
  `simple_stream_master_resources`.

## Order of calls and couplings

- Stage: `new(cline)` (restart read before `params%new`; `init_gui` takes -1 for no pipe) ->
  `iterate` until `finished` -> `finalize` -> `kill`. A new parameter's default goes in the
  commander's `set_<stage>_cline` (chain tests use it), and in p00's `make_stage_clines` if the GUI
  run needs it (forwarded only when given). Restore changes: RESTART header, restart policy §3.
- `simple_state_assignment` is shared with `solve3D_pyramid`, `automask2D_mskdiam` with
  `make_pickrefs`; p04's `moldiam.txt` is p05's mask (p05 stops without it).

## Tests

- `build/bin/simple_test_exec test=unit_stream` (about 4 s; `suite=` preprocessing,
  optics_assignment, initial_analysis, reference_picking, particle_sieving, pool_2d, pool_2d_object,
  solve_3d, stream_master, job_sets, meta_plots, stream_watcher); `unit_ipc suite=stream_pipe`;
  `unit_project suite=micrograph_import|micrograph_selection|optics_groups|optics_maps`;
  `unit_ui` (GUI metadata, assembler); `forked_process` (stream heartbeat);
  `unit_reconstruction suite=state_assignment_rules`.
- `test=lib_stream`: in-process stages (optics assignment waits a minute).
- `test=stream_chains` (`suite=sieve_to_3d|movies_to_3d`): stages chained on local worker jobs,
  after `make -C build install`; tens of minutes each; gated on a planted two-state truth. A chain
  has hung once after "MOVIES TO PROCESS: 0".
- `test=stream_preproc` (82 s), `test=mini_stream suite=6vxx|1jxy` (1 min; batch path via
  `simple_mini_stream_utils`). No test restarts a chained stage or runs the master with stages.

## See also

`doc/policies/stream/` (README, `ipc_`, `restart_`, `pool2D_`, `reference_generation_`,
`stream_3D_ingestion_policy.md`, `pNN_overview.md`), `doc/policies/sieving_and_rejection/ptcl_sieve_policy.md`
(the particle sieve, p03 and p05), `doc/refactoring_notes/completed/stream_refactor.md` (each module, and how
each stage differs from the one it replaced); layering policy section 8; skills
`simple-main-sieve`, `simple-solve2d`, `simple-main-star`.
