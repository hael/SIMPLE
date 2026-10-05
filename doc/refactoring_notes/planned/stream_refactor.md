# Stream refactor: the stages as clients of the library

**Status:** planned follow-up. The stage switch-over described below is
complete, but the note's remaining cleanup, workflow runs, and test-registration
work are still open; those items keep this record in `planned`.

Stream stages rewritten as clients of the library: p01 (preprocessing), p02 (optics
assignment), p03 (initial analysis and picking references), p04 (reference-based picking
and extraction), p05 (particle sieving), p06 (pool 2D), p07 (multistate 3D), and the master
(p00) that runs them for the GUI. This follows recommendations
E5 to E8 of
`doc/refactoring_notes/completed/stream_area_review_2026-09-30.md`. The stream runs on this code:
`production/simple_stream.f90` uses these commanders, and the old stage modules of
`src/main/stream` are deleted. The code was written in `src/main/stream_refactor` and then
moved into the library layout:

- `src/main/commanders/stream`: the commanders `simple_commanders_stream_p00_master` to
  `simple_commanders_stream_p07_solve3D_multistate` (types `commander_stream_pNN_*`);
- `src/main/stream/stages`: the stage types `simple_stream_stage_<x>` (`stream_stage_<x>`)
  and their testers;
- `src/main/stream/master`: the master's parts `simple_stream_master_*` (stage ids and
  pipes, the stage record and its fork, the GUI commands, the metadata store) and its tester;
  the fork runs the commander the master gives it, so these depend on commanders only as
  `commander_base`;
- `src/main/stream/shared`: the modules the stages share (pipe and pipe descriptors
  `simple_stream_state`, sigterm, job sets, GUI senders, meta plots, and the micrograph and
  optics helpers) and their testers;
- `src/main/stream/pool2D`: the 2D pool and chunk layer that stayed (`simple_stream2D_state`,
  `simple_stream_pool2D_utils`, `simple_stream_refine2D_utils`, `simple_stream_chunk`,
  `simple_stream_chunk2D_utils`);
- `src/main/stream`: the watcher, `simple_stream_utils`, `simple_mini_stream_utils`,
  `simple_stream_solve2D_chunks` and the library tests of the stages (`simple_stream_tester`).

**Build note.** `src/CMakeLists.txt` collects `main/*.f90` recursively, so these
modules are compiled into the library on the next configure. The `*_tester` modules
are compiled only with `BUILD_TESTS=ON`.

## Contents

### Shared modules

| file | what it is | replaces |
|---|---|---|
| `simple_stream_pipe.f90` | `stream_pipe`: length-prefixed framed messages on one pipe pair, plus `send_meta` for GUI metadata. A corrupt frame length is reported and the buffered bytes dropped (resync); `discard` empties a pipe before a restart | the `send_to_*_in_pipe` / `receive_from_*_out_pipe` copies in each stage, and the master's `pipe_rx_state`/`pipe_tx_state` framing |
| `simple_stream_master_stage_ids.f90` | the stage ids (1..7, in heartbeat order), their job names, GUI keys and log labels, and their pipes: open, close, close the others' ends in a forked stage, and the master's and the stage's ends, on the named arrays of `simple_stream_state` | the master's hand-written lists of the 14 pipes |
| `simple_stream_master_stage.f90` | `stream_master_stage`: one stage as the master runs it (fork, command line, a reader of its messages and a writer of updates); `stream_master_stage_fork`, the one fork type, whose `execute` runs the stage's commander | the master's 7 fork types, 7 `x*` routines and the per-stage start/stop/restart blocks |
| `simple_stream_master_gui_commands.f90` | `stream_master_gui_commands`: the GUI's answer to a heartbeat (stop all, stop or restart a stage, the updates of this answer only) parsed from its JSON | the JSON handling in the master's loop |
| `simple_stream_master_meta_store.f90` | `stream_master_meta_store`: the latest metadata each stage sent, by type; list items placed by `i`/`i_max`; `assemble` for the heartbeat | the master's metadata variables and the 11 copies of the resize-and-place code in its listener |
| `simple_stream_meta_plots.f90` | fill GUI histogram, windowed time-plot and rate metadata from an `oris` segment; recent shifts per optics group | about 125 lines in p01's main loop; p02's per-group shift scan |
| `simple_mic_selection.f90` | `reject_mics_by_thresholds`: the micrograph ctfres/icefrac/astig rule; `reject_mics_without_particles`: the rule after picking and extraction | 5 inline copies of the first (p01 ×2, p04 ×2, `selection` commander), which disagreed; 2 of the second in p03 |
| `simple_mic_import.f90` | `append_mics_from_projects`: append the `mic` segments of upstream project files, reading each once | the per-stage import loops (p01 restart, p02, p03's cycle-1 list, which read a file once per micrograph; p04 still has its own) |
| `simple_optics_groups.f90` | `assign_optics_groups`: optics groups from beam-image shifts within tilt groups | the assignment hidden in `starproject_stream%stream_export_optics` |
| `simple_optics_maps.f90` | `publish_optics_map`, `import_latest_optics_map`, `latest_optics_map_id`: the optics-map file protocol | p02's pruning loop; for readers, the three "latest map" scans |
| `simple_stream_sigterm.f90` | `install_sigterm_handler`, `sigterm_received`, `restore_sigterm_handler`: a SIGTERM handler that only sets a flag, polled by the commanders' loops and by long stage steps | a handler and flag in each stage commander |
| `simple_stream_gui_senders.f90` | `send_cavgs`: class averages or picking references as sprite-sheet tiles; `send_recent_micrographs`: the latest micrographs' thumbnails, CTF values and picked positions | the copies in p03 (two) and p04, whose tile positions divided by zero for one to three images |
| `simple_stream_job_sets.f90` | `stream_job_sets`: the numbered sets a stage runs as streaming jobs. `write_set`, `submit`, `schedule`, `collect`, `complete`, and `restore` (numbering continued after the highest completed set, unfinished sets dropped, each set's recorded origin returned) | the set naming, submission, collection, move to the completed folder and restart counter repeated in p01 and p04 |

### Stages

| file | what it is | replaces |
|---|---|---|
| `simple_stream_stage_preprocess.f90` | `stream_stage_preprocess`: p01's state and steps as a type (`new`, `iterate`, `finished`, `finalize`, `kill`). `new` is a sequence of named steps (`init_params`, `init_movie_watcher`, `resume_previous_run`, `init_job_dirs`, `init_queue`, `resolve_gain`, `build_worker_cline`, `init_gui`). Gain flip detection and generation only supply movie batches to the motion subsystem. Components and steps are public so the tester can assemble a stage step by step; the waits (`settle_s`, `wait_s`, `sniff_wait_s`) are components that tests set to 0, and `init_gui` takes the pipe ends. | the 520-line `exec_stream_p01_preprocess` body and its internal procedures |
| `simple_commanders_stream_p01_preprocess.f90` | commander `commander_stream_p01_preprocess`: normalises the command line, then loops the stage | `simple_stream_p01_preprocess_new` |
| `simple_stream_stage_optics.f90` | `stream_stage_optics`: p02 as import → assign → STAR → publish → report. Components and steps are public so the tester can run one step at a time; the upstream settle time (`settle_s`) and the idle pause (`wait_s`) are components that tests set. | the `exec_stream_p02_assign_optics` body |
| `simple_commanders_stream_p02_assign_optics.f90` | commander `commander_stream_p02_assign_optics` | `simple_stream_p02_assign_optics_new` |
| `simple_stream_stage_initial_analysis.f90` | `stream_stage_initial_analysis`: p03's two-cycle plan with named steps; jobs through `qsys_async_job`, picking through `segdiam_bin_picker`, selection through `simple_cavg_quality_selection`. `new` is `init_params`, `init_queue`, `init_gui(fd_read, fd_write)`. Components and steps are public so the tester can run one step at a time without a queue; `balance_classes` and `find_final_solve3D_cavgs_dir`, which use no stage state, are bound with `nopass` for the same reason; the settle time and the end-of-pass pause are components. | the 440-line `exec_stream_p03_initial_analysis` body and about 25 internal procedures |
| `simple_commanders_stream_p03_initial_analysis.f90` | commander `commander_stream_p03_initial_analysis` (launched as `prg=gen_pickrefs`, an alias of `initial_analysis`) | the old `simple_commanders_stream_p03_initial_analysis` |
| `simple_stream_stage_refpick.f90` | `stream_stage_refpick`: p04 as watch → one job set per upstream project → pick_extract → import → project write. Public components and steps, waits as components, as in the other stages. | the `exec_stream_pick_extract` body and its internal procedures |
| `simple_commanders_stream_p04_refpick_extract.f90` | commander `commander_stream_p04_refpick_extract` | `simple_stream_p04_refpick_extract_new` |
| `simple_stream_stage_sieve.f90` | `stream_stage_sieve`: p05 as watch → import → particle sieve cycle → report. The sieve hands its finished chunks to the completed folder with the newest optics map's groups. Public components and steps, waits as components, as in the other stages. | the `exec_stream_p05_sieve_cavgs` body and its internal procedures |
| `simple_commanders_stream_p05_sieve_cavgs.f90` | commander `commander_stream_p05_sieve_cavgs` | `simple_stream_p05_sieve_cavgs_new` |
| `simple_stream_stage_pool2D.f90` | `stream_stage_pool2D`: p06 as watch → refresh → publication of the pool's classified state for 3D (once an iteration has come back, before the next dispatch) → import into the pool → pause rule → pool iteration → GUI updates and snapshots. The pool itself is still the module state of `simple_stream_pool2D_utils` and `simple_stream2D_state` (one stage per process); the stage's command line is a pointer component, which the pool reads as `master_cline`. The import is split into `transfer_sets` (into any project, for the tester) and `start_pool`; the pause rule, the particle targets and the default mask diameter are `nopass` functions. Public components and steps, waits as components, as in the other stages. | the `exec_stream_p06_pool2D` body and its internal procedures |
| `simple_commanders_stream_p06_pool2D.f90` | commander `commander_stream_p06_pool2D` (launched as `prg=pool2D`) | `simple_stream_p06_pool2D_new` |
| `simple_stream_stage_solve3D.f90` | `stream_stage_solve3D`: p07 as watch → (no job running) select class averages and import → start or collect a 3D job → report. The jobs are `qsys_async_job`s with an explicit phase (importing, solve3D, idle, addon); which job starts next is the `nopass` function `next_job`, and the mask diameter the class averages get is `fit_mskdiam`. Public components and steps, waits as components, as in the other stages. | the `exec_stream_p07_solve3D_multistate` body and its internal procedures |
| `simple_commanders_stream_p00_master.f90` | commander `commander_stream_p00_master`: parameters, queue, the stage table, a listener thread filling the store, and the heartbeat loop that applies the GUI's commands; the stages' settings are a table of named constants | the old `simple_commanders_stream_p00_master` |
| `simple_commanders_stream_p07_solve3D_multistate.f90` | commander `commander_stream_p07_solve3D_multistate` (launched as `prg=solve3D_stream`), with the 3D job settings as named constants | the old `simple_commanders_stream_p07_solve3D_multistate` |

### Testers (the stage and job-set testers are registered in `unit_stream`; the others are not yet, see "Switching over")

| file | suite it belongs in |
|---|---|
| `simple_stream_pipe_tester.f90`: framing through a real pipe, including the queued-frame regression, the resync after a corrupt length, and `discard` | `unit_ipc`, sub-suite "stream pipe" (**registered**) |
| `simple_stream_master_tester.f90`: the stage names and keys, the GUI's commands (a stop, a restart and updates; each answer's updates alone; text that is not JSON), the store (a status replaced; a list made, an item out of range dropped, the list remade for a new `i_max`; a volume sent with tiles arriving without), a stage's pipes from both sides and the discard before a restart. Forking, the heartbeat and HTTP are left to the high-level tests. | `unit_stream`, sub-suite "stream master" (**registered**) |
| `simple_mic_selection_tester.f90`: the threshold rule and the post-extraction rule | `unit_project` |
| `simple_mic_import_tester.f90`: append order, the accepted-only filter, growth across calls | `unit_project` |
| `simple_optics_groups_tester.f90`: the stream optics test's two-cluster case with and without beam tilt, the group offset, the CTF constants row | `unit_project` |
| `simple_optics_maps_tester.f90`: pruning to the newest five, importing the newest map, an empty directory, a project copied with the map's groups or exactly without a map (both renamed into place) | `unit_project` |
| `simple_stream_stage_preprocess_tester.f90`: p01's steps one at a time, with no queue and no waits: `split_mode` after `init_params`, importing finished jobs (accepted, failed, moved), thresholds and the STAR file after an import, the restart set counter, a movie-set project, the worker command line, a static gain flip, queued GUI updates, one status message, `finished` | `unit_stream`, sub-suite "preprocessing" (**registered**) |
| `simple_stream_stage_optics_tester.f90`: p02 through its own `new`, with no waits: the segment starts empty on restart, a restart removes a leftover `TERM_STREAM`, `beamtilt=yes` from the command line groups by tilt group, waiting for the upstream folder, importing each completed project once (rejected micrographs included), grouping the two shift clusters with the STAR files, the project and the first map, map ids continuing after a restart, the public loop (waiting, importing, idle, finished at `nmics`, `finalize`), one status message, one shift message per group, `finished` | `unit_stream`, sub-suite "optics assignment" (**registered**) |
| `simple_stream_stage_initial_analysis_tester.f90`: p03's steps without a queue: the stage project with its computing environment, the paths, waiting for the upstream folder, one record per accepted imported micrograph, the cycle 1 micrographs, cycle 1 set-up then waiting without a job, the public loop up to that point, a GUI selection of references (written, sent back, the stage finished; no cycle or out-of-range indices ignored), the status messages, class balancing to 501 rows with the even/odd stacks, the solve3D restart directory, `finished`. Picking, extraction, the 2D/3D jobs, the sieve and the reprojected references are left to the high-level tests. | `unit_stream`, sub-suite "initial analysis" (**registered**) |
| `simple_stream_stage_refpick_tester.f90`: p04's steps without a queue: `init_params`, waiting for the upstream folder and putting a restart's upstream projects in the watcher history, waiting for the picking references, a set of an upstream project's accepted micrographs (pixel size, no picking-preprocessing fields; none for a project with nothing accepted), importing finished sets (moved and imported unchanged, the empty set kept), the project write (stacks and particle ranges over two sets), the optics groups of a published map (micrographs, stacks, particles, STAR file, project) and one group without a map, restart (completed set re-imported and its upstream project in the history, unfinished set dropped and its upstream project picked again), restart with `clear=yes`, the status message, the public loop up to the picking references, `finished`. Submission, make_pickrefs and the pick_extract jobs are left to the high-level tests. | `unit_stream`, sub-suite "reference picking" (**registered**) |
| `simple_stream_job_sets_tester.f90`: the folders, set naming and writing (project info, origin, worker command line), completion, restore (numbering, origins, unfinished sets dropped) | `unit_stream`, sub-suite "job sets" (**registered**) |
| `simple_stream_stage_sieve_tester.f90`: p05's steps without the chunk workers: `init_params`, a restart's chunked projects in the watcher history once attached, one record per imported micrograph, the mask diameter from `moldiam.txt`, the sieve made with it on the first import (warm-up cycles, no chunk below the threshold), the status message, the public loop up to attaching, `finished`. Chunk 2D and rejection are covered by the ptcl_sieve tests. | `unit_stream`, sub-suite "particle sieving" (**registered**) |
| `simple_stream_stage_pool2D_tester.f90`: p06's steps without starting the pool: `init_params` (the stage's own command line with `mkdir=no`, the caller's unchanged; stepwise; the optics id offset), the restart clean-up (pool files removed, the project and the snapshots kept), the publication numbering continued after the highest publication, a publication holding only the classified stacks (renumbered, image indices, nothing before any iteration), waiting for the sieve's folder and one record per set with the sieve's mask diameter, the transfer into a project (stacks renumbered, particles as new ones keeping their shifts, the counts, each set once, a later set appended), stepwise import, the sieve's final set, the pause rule and particle targets, the default mask diameter in Angstroms, a GUI mask diameter (taken, resumes, replaces the sieve's; the same again ignored; a snapshot request before the pool ignored), the status message, a snapshot's messages, the public loop up to attaching, `finished`. Pool iterations and snapshots of a running pool are left to the high-level tests. | `unit_stream`, sub-suite "pool 2D" (**registered**) |
| `simple_stream_stage_solve3D_tester.f90`: p07's steps without a queue: `init_params`, a restart removing a leftover `TERM_STREAM`, waiting for pool 2D's folder, publications recorded in order (each once) with the first one's mask diameter, merging publications into rows that only grow (a new stack appended with both particle segments; a stack held matched by name in any order, taking class and selection in place and keeping its 3D parameters; a missing stack deselected; the newest publication's classes), `next_job` and `fit_mskdiam`, the status message, the volume messages of a fixture project (one per state with a volume, three reprojection tiles, the resolution kept), the public loop up to attaching, `finished`. Class-average selection and the 3D jobs are left to the high-level tests. | `unit_stream`, sub-suite "solve 3D" (**registered**) |
| `../cavg_quality/simple_cavg_quality_selection_tester.f90`: the selected/rejected stacks, their order, replacement on rewrite | `unit_numerics`, next to the cavg quality relations |

## Library changes beyond these modules

These are live in the library now.

### Already used by the running pipeline

- **`simple_cavg_quality_selection`** (new, `src/main/cavg_quality`):
  - `score_project_cavgs` reads a project's class averages and scores them, with a quality
    model and the relation parameters built from the box and mask diameter, or with the hard
    gates only.
  - `write_cavg_selection_stacks` / `write_cavg_stack` write the selected and rejected stacks.

  Two live callers have switched to it. Each behaves as before:
  - `ptcl_sieve%reject_cavgs` scores through it. Two small differences: the relation
    parameters are no longer built when model rejection is off, which had no effect; and a
    chunk without class averages now stops with a clear error instead of indexing past an
    empty array.
  - The `model_cavgs_rejection` commander writes its two stacks through it, with the same
    files and the same log lines.
- **`simple_segdiam_bin_picker`** (new, `src/main/pick`): one picker type in place of the two
  300-line routines `segdiampick_mics_multi` and `segdiampick_mics_multi_fixed_bins`. p03
  drives the picker directly; the two wrappers that replaced them were removed on 2026-10-05.
  `MOLDIAMS_PICK` moved with the picker. The only output difference: the two "no accepted
  cluster" warnings became one message.
- **`simple_qsys_env%exec_simple_prg_in_queue_async`** has a new optional `exit_code_fname`.
  Existing callers are unaffected.
- **`forked_process%start`** resets SIGTERM and SIGINT to their defaults in the child,
  before `execute`. Stages the master forks at start-up already had the defaults; a stage
  restarted from the GUI inherited the master's handlers. In such a stage, Ctrl-C ran the
  master's SIGINT handler, which joins a thread the child does not have, and a SIGTERM
  arriving before the stage installed its own handler was swallowed, so the master waited
  for it forever. Now every forked stage starts the same way.
- **The memory monitor** (`simple_posix.c`) is fork-safe: `pthread_atfork` handlers make a
  child forget the parent's monitor, where it used to join the parent's sampler thread and
  hang. With `memreport` on, every GUI restart hung this way, and every stage when
  `memreport=yes` was passed to the master. The master also forks only before its listener
  thread starts or while holding the listener's lock.
- **Stage restarts are the GUI's alone.** A GUI restart used to fork with `forked_process`'s
  auto-restart on, so the stage then re-forked itself on any failure, including the SIGKILL
  that ends a stage which ignores a stop. `stream_master_stage%start` takes no argument, and
  `forked_process` has no auto-restart any more (removed on 2026-10-05).
- **`ptcl_sieve%new`** takes an optional `optics_dir`. With it, a finished chunk handed to
  `completedir` gets the groups of the newest optics map, applied by import index
  (`hand_off`, through `copy_project_with_optics_map`); the chunk's own project is not
  changed, because fine chunks merge coarse ones and the merge offsets group ids per
  source. Callers that do not pass it (the running p05, p03's sieve) get an exact copy as
  before. The sieve policy and tester are updated (`test_hand_off_applies_optics_map`; the
  hard-gates test's fixture is now `make_completed_coarse_chunk`, shared by both).
- **The sieve's hand-off is atomic.** The copy is written as `<stem>.tmp` in `completedir`
  and renamed to `<stem>.simple`, with or without an optics map, so pool 2D never reads a
  partial project (it watched with a 60 s settle time for that reason). The `.tmp` name is
  off the watcher's `\.simple$` filter. Policy and testers check that no `.tmp` is left.
- **`sp_project%write(tempfile=.true.)`** works. It swapped the suffixes the wrong way
  round (`swap_suffix` takes the new suffix first), so the temporary name was `.simple` in
  the current directory; nothing used the option. `copy_project_with_optics_map` now does.
- **Pool 2D fixes** (`simple_stream_pool2D_utils`, `simple_stream_refine2D_utils`, used by
  the running p06 and the new one):
  - A new mask diameter (`update_mskdiam`) also sets the cropped mask radius on the pool's
    command line. The workers take `msk_crop` over the one `parameters` derives from
    `mskdiam`, so GUI and automatic mask updates had no effect on the pool. The uncalled
    `update_user_params2D` gets the same fix for chunks and pool.
  - New particles drawn into a populated class (`iterate_pool`) no longer loop forever when
    no class is populated; they keep the class they have.
  - `get_pool_cavgs_jpeg` and `get_pool_cavgs_mrc` return absolute paths. The JPEG was
    relative after `update_pool` and absolute after `generate_pool_jpeg`, and the running
    p06 prepended the working directory to both; it no longer does.
  - After the first iteration the pool has no classes yet, and the JPEG tile layout divided
    0 by 0 (trapped in Debug builds); at least one column now.
  - A snapshot of an iteration outside the pool's history, or of one pruned from it (after
    five minutes), is reported as failed instead of indexing past the history or copying an
    empty project.
  - A snapshot's class averages are registered with their own pixel size (the pool's
    cropped one); they were given the original pixel size.
- **Pool 2D's exports for 3D** (`write_project_stream2D(export=.true.)`) are written as
  `<id>.tmp` and renamed, and their class averages record the pool's mask diameter. The 3D
  stage reads it (it had none: quality selection got −7.8 Å, solve3D 5000 Å, and both fell
  back to the box default with a warning) and settles for 2 s instead of 60. Changed since
  (2 October, review D3/R3): the deltas are replaced by publications of the pool's classified
  state (`publish_pool_state`); see `doc/policies/stream/stream_3D_ingestion_policy.md`.
- **`nparts3D` and `nthr3D`** (new parameters, 0 by default): the parts and threads of the
  stream's 3D jobs, apart from the 3D stage's own `nparts` and `nthr`.
- **`gui_metadata_vol3D%set_oridist_from_oris`** bins a state's particle orientations into
  the orientation-distribution histogram. `gui_metadata_project` uses it in place of its copy
  of the loop (which was the 3D stage's), and the new 3D stage uses it too.
- **`gui_metadata_vol3D%serialise`** sends a copy without `reprojtiles`. The base
  serialisation copies the object's bytes, which for an allocated component would be a
  descriptor of the sender's memory that the master's `transfer` then follows and frees. No
  sender set the tiles, so nothing changes now; they travel as their own messages, and
  `gui_metadata_project` jsonises its tiles in process.
- **`score_project_cavgs`** takes an optional `smpd`, the class averages' pixel size for the
  relation parameters. They turn the mask diameter into the mask radius in pixels of the
  relational feature's signal statistics (the pairwise correlations take the images' own
  pixel size). Without it the parameters default of 1.3 Å was used, so the radius was off by
  the ratio of the pixel sizes (capped at the box), unlike in the `model_cavgs_rejection`
  commander, whose tables train the models. The 3D stage and p03 pass it. The sieve does not
  need to: it passes a mask diameter of 0, which falls back to the box-based radius whatever
  the pixel size, so its scores do not depend on it.
- **`beamtilt` is a parsed parameter** (`simple_parameters_parse`), and the stream master
  forwards it to the optics-assignment stage when the user sets it. Before, the GUI's "Use
  beamtilts in optics group assignment" option never reached `params%beamtilt`, so the
  stream always grouped as if it were `no`. With `beamtilt=yes`, micrographs are now
  grouped by EPU tilt group before the shift clustering, in the running stage
  (`simple_stream_p02_assign_optics_new`, through `stream_export_optics`) and in the new
  one. The batch `assign_optics_groups` reads its command line directly and is unchanged.
  The master forwards `tilt_thres` the same way, so a threshold changed in the GUI now
  reaches the optics stage; both defaults are 0.05, so runs that leave it alone are
  unchanged.
- **`simple_stream_p02_assign_optics_new`** deletes a leftover `TERM_STREAM` on restart
  after `params%new`, in the stage directory its loop reads, instead of before it, in the
  master's directory.

### Not yet used by the running pipeline

- **`simple_qsys_async_job`** (new, `src/utils/qsys`): one program run started with
  `start(qenv, cline, dir, label)` and polled with `status()` (idle, running, done or
  failed, from the exit status the script writes).
- `simple_motion_gain_analysis`: `gain_flip_analyzer%get_flip_mode()` returns `'no'`,
  `'x'`, `'y'` or `'xy'`. Callers no longer decode `best_idx`.
- `simple_motion_gain_helpers`:
  - `add_movies_to_gain_sum`: adds one batch of movies to a running frame sum.
  - `write_gain_from_sum`: writes the normalised inverse average as a gain reference.
- `simple_motion_gain_tester` (`unit_project`, "motion gain") has tests for all three.
- `simple_starproject_stream`: a new `stream_write_optics` writes `optics.star` from
  the groups already in the project and assigns nothing. (`stream_export_optics`, which
  called it after assigning, was removed on 2026-10-05.)
- `simple_test_utils`: `enter_fixture` / `leave_fixture` (a per-test directory, kept only
  when a check failed) moved here from `simple_stream_tester`, which now uses these.

## Behaviour compared with the old `simple_stream_p01_preprocess_new`

The new p01 does what the current one does, apart from these deliberate changes:

1. **Queued GUI updates.** A complete message already buffered is delivered even when
   the next read would block (EAGAIN). Before, it stayed stuck until more bytes arrived.
2. **Restart set counter.** The movie-set counter is restored from every completed
   project, even when all of them were rejected. Before, it restarted at 0 in that case,
   so new set projects could overwrite completed ones of the same name.
3. **Restart with nothing accepted.** The old code tried to delete those projects but
   built the path twice, so nothing was ever deleted. The dead delete is gone and the
   files are kept, as they are today. Decide whether they should be deleted.
4. **Sliding-window drift check.** Removed. It computed statistics whose only consumers
   (GUI alerts) were commented out, so it had no effect.
5. **Parameters.**
   - `split_mode = 'stream'` is still assigned right after `params%new`, as in the old
     stage. It cannot go on the command line: `derive_parallel_settings` resets it to
     `'even'`, and the queue would then get one partition but `ncunits` computing units.
     (A first version set it on the command line and failed with an out-of-bounds
     `jobs_done` in `schedule_streaming`.)
   - The `ncunits = nparts` assignment is dropped. It restated the default, except that it
     overrode a user-supplied `ncunits`, which is now respected.
   - `params%updated` is no longer set (only p02 read it).
   - Apart from `split_mode`, the only fields changed after `new` are run-time values: the
     resolved gain and GUI threshold updates. Each is changed in one place, together with
     the command line.
6. **Gain options.** `no`, `x`, `y`, `xy` and `yx` are accepted as well as the GUI names
   (`none`, `flip_*`). Before, the parameter default `no` stopped the stage with
   "Unknown gain processing option". Workers receive `flipgain=no`, because the master
   has already resolved and flipped `gainref`. Workers do not flip in the stream path
   today, so this only removes a trap.
7. **Signal handling.** The SIGTERM handler only sets a flag (`simple_stream_sigterm`), and
   the log line is written from the loop. Before, the handler wrote to the log itself. The
   wait for the first movie in `new` also stops on SIGTERM; before, the stage could not be
   stopped until a movie arrived. The default action is restored when the commander
   returns.
8. **Lifecycle.** `qenv`, the watcher, the pipe, the worker command line and the GUI
   metadata objects are killed at the end. Job projects are killed even when none of
   their micrographs is accepted.
9. **Micrograph rule.** Unchanged in effect for p01: strictly greater than the threshold,
   the row's own key, rejected rows left alone. p04 used `threshold - 0.001` and will
   change slightly when it is switched to the shared rule.
10. **User-parameter file.** The final `update_user_params` call is gone. It read
    `stream_user_params.txt`, which nothing in the repository writes, NICE included.
11. **Job sets.** The movie sets go through `stream_job_sets`. A set is written by
    absolute path rather than after changing into the job folder. On restart the
    unfinished sets in `spprojs/` are always dropped; before, only when an accepted
    micrograph came back from the completed sets.

## Behaviour compared with the old `simple_stream_p02_assign_optics_new`

1. **Optics groups** come from `assign_optics_groups`: the same tilt grouping, `h_clust`
   threshold clustering, centroids, ids and names as the STAR exporter's version. Two
   fixes:
   - there is no fixed limit of 10,000 groups;
   - when no micrograph is accepted, the CTF constants for the optics rows are taken
     from the first micrograph. Before, the code read one row past the end.
2. **No hidden writes.** `optics.star` is written by `stream_write_optics`, which only
   writes. The project write the exporter did as a side effect is kept, but done
   explicitly in the stage. Drop it if nothing reads the stage's project during a run.
3. **Import.** Every micrograph of each completed project is imported, rejected ones
   included as before, and each file is read once instead of twice. The old cap of
   `STREAM_NMOVS_SET` per project made no difference, because p01's projects hold
   exactly that many.
4. **Optics maps.** Maps are published in the stage directory, and their ids continue
   from the newest map there. Before, the restart lookup searched `params%outdir`, a
   relative path that does not exist inside the stage directory, so the ids most likely
   restarted at 1. Readers would then keep picking the previous run's newest map. Each
   publication now deletes only map `id-5`, not every older id.
5. **Group shift plot.** Groups are matched on the `ogid` stored in `os_optics`, in one
   pass from the newest micrograph back. The output is the same with the stream's group
   offset of 0.
6. **User-parameter file.** `update_user_params` and the `params%updated` re-export are
   gone, because nothing writes `stream_user_params.txt`. `tilt_thres` and `beamtilt`
   can no longer be changed live, but no GUI offered that. Adding it would mean a NICE
   control, two fields in `gui_metadata_stream_update`, and re-enabling the master's
   send to this stage.
7. **Waiting for preprocessing** happens inside `iterate`, so SIGTERM and `TERM_STREAM`
   are honoured while waiting. The old code gave up after 24 hours and then watched a
   folder that did not exist.
8. **Settle time.** A completed preprocessing project is taken once untouched for
   `SHORTWAIT` (2 s) instead of `LONGTIME` (60 s). Preprocessing moves finished projects
   into `spprojs_completed/` with a rename, so they are complete when they appear, and
   optics groups now follow preprocessing within one pass instead of a minute later.
9. **Smaller fixes.**
   - The restart check only runs for a non-empty `outdir`.
   - On restart, `TERM_STREAM` is deleted after `params%new`, in the stage directory where
     `finished` looks for it. Before, it was deleted before `params%new` moved there, in
     the master's directory, so one left by the previous run stopped the restarted stage
     at once. The running stage (`simple_stream_p02_assign_optics_new`) is fixed the same
     way.
   - The unused `last_micrograph_imported` is gone.
   - SIGTERM is logged.

### Known limitation carried over

Every import reclusters all micrographs. `h_clust` builds an N×N distance matrix for
each tilt group. Without beam tilt there is one group holding every micrograph, so at
50,000 micrographs a single import allocates about 10 GB. Group ids can also change
between imports. Fixing this means clustering incrementally, which is a method change
outside this refactor.

## Behaviour compared with the old `simple_commanders_stream_p03_initial_analysis`

The plan, the method and the files written are the same. The integer plan codes became
named steps. The differences:

1. **Jobs that fail.** Extraction, solve2D and solve3D_cavgs run as
   `qsys_async_job`s, whose exit status tells done from failed. Before, the stage polled
   for a marker each program writes on success, so a crashed job left it waiting forever.
   Now:
   - a failed cycle-1 extraction, either solve2D, or solve3D_cavgs stops the stage
     with the job's log path;
   - a failed extraction of the "all" set is logged and left out, and is counted so the
     sieve's final ingestion is not held up.
2. **Sieve lifecycle.** The sieve is no longer cycled after it is killed (review finding 12).
3. **SIGTERM** sets a flag (`simple_stream_sigterm`). The old handler called `exit(0)`,
   which could truncate a project file mid-write. The flag is also checked between the
   projects of the in-process picking loop and before each cycle step, so a stop does not
   wait for a whole picking pass, and no job is started once it is requested. The default
   action is restored when the commander returns. Jobs already running are not stopped.
4. **Waiting for preprocessing** happens inside `iterate`, as in p02. There is no 24-hour
   give-up. As in p02, a completed preprocessing project is taken once untouched for
   `SHORTWAIT` (2 s) instead of `LONGTIME` (60 s); it arrives by rename, so it is complete.
5. **Cycle-1 micrographs** are built by reading each upstream project once. Before, the
   file was read once per micrograph.
6. **Micrograph rejection.** The post-extraction rule is `reject_mics_without_particles`.
   After class selection the same rule replaces a narrower check (state 0 or no
   particles). The two agree in practice, because extraction already required the box
   files.
7. **Parameters.** `workers` and `worker_nthr` are set on the command line, not on
   `params` after `new`. `worker_nthr` is set only when `nthr` is given, which the master
   always does (32). The `nmics` default is gone: the stage keeps its own target, and
   nothing read `params%nmics`.
8. **Termination message.** The GUI's last 2D status now reports the extraction box. It
   used to send 0 from a variable that was never set.
9. **Lifecycle.** The class-compatibility model is killed after each selection, and the
   reprojection command line after use.
   **Sprite positions.** One to three class averages or references form a single column
   (`mrc2jpeg_tiled` uses `floor(sqrt(n))` columns). The tile position was computed with
   `merge(0.0, x * (100.0 / (xtiles - 1)), xtiles == 1)`, and `merge` evaluates both values, so
   that case divided by zero, which Debug builds trap. The step is now computed only for two or
   more tiles. The live p03 (two places) and p04 still have the old expression.
10. **Not ported:** dead code.
    - `micimporter`, the first `run_cavg_quality_selection` and `run_cavg_size_selection`;
    - the `imgfiles` and `smpd_stk` caches, which nothing read;
    - `NSTAGES3D`, and about 12 unused imports.
11. **Kept as they were:**
    - The loop does not check `TERM_STREAM`; the master stops this stage with SIGTERM.
      The restart branch no longer deletes `TERM_STREAM`: it ran before `params%new`, in
      the master's directory, where nothing writes that file.
    - A GUI selection of references (for cycle 1 or 2) still ends the stage at once,
      wherever the plan is, as the running stage's `exit main_loop` does: the process stops
      after `finalize` and `kill`, the rest of the plan is skipped, and jobs already
      submitted (the "all" extractions, the sieve's chunks, solve2D/3D) keep running
      unattended. Changed since (stream area review, 2 October, D2): the references, from a
      selection or the 3D route, are published once as `OPENING2D_PICKREFS`, the file the
      master points reference picking at, by a rename. A selection is read before the cycle
      steps of each pass, so it pre-empts the 3D route. A selection that publishes nothing no
      longer ends the stage. Published references are final, and a restarted stage that finds
      them is finished at once.
    - Class balancing and references from a 3D volume are unchanged (B1/B2).
    - The sieve still takes a copy of `params` with nine fields edited (E8 proposes a
      `sieve_config`).
12. **GUI senders** are the shared ones (`simple_stream_gui_senders`). Micrograph
    thumbnails are sent when the segment has `thumb`, the key they are read from; the old
    code tested `thumb_den` (with a TODO). A micrograph whose box file is missing is sent
    without positions; before, it was not sent.
13. **Quality selection at the class averages' pixel size.** The relational feature's mask
    radius is computed with the class averages' pixel size, as when the chunk model was
    trained, instead of 1.3 Å (see `score_project_cavgs` above). With class averages coarser
    than 1.3 Å the radius was inflated and usually capped at the box, so the signal
    statistics covered the whole box; finer, it cut into the particle. Borderline classes,
    and so the picking references, can change; worth comparing on one dataset.
14. **The mask diameter of cycle 2 and 3D is estimated** (3 October) from cycle 1's selected
    class averages, by `make_pickrefs`' measure and rule, and capped at the box's default
    (`estimate_mskdiam`; `reference_generation_policy.md` section 3.1). Before, every step
    used the box's default: cycle 1 and cycle 2's `solve2D` through `mskdiam=999`, which fell
    back to it with a warning, the sieve's chunks through `mskdiam=0`, which did the same,
    and 3D and the reprojection through the picker's diameter, which is that default. Cycle 1
    and the sieve now pass the default explicitly. Cycle 2's starting low-pass limit, which
    `solve2D` derives from the mask diameter (`min(max(mskdiam/12, 15), 20)` Å), can drop, to
    no less than 15 Å.
15. **The state of the references is chosen on shape first** (3 October; `choose_state`): the
    fewest connected components of the binarised volume, then the most distinct projection
    directions of the state's classes, then the largest population. Before, the directions
    alone decided (the shape descriptors were only logged), and the choice could fall on a
    state without a volume. The descriptors' mask radius is now in voxels; it was the mask
    diameter in Å, so they covered the whole box.

## Behaviour compared with the old `simple_stream_p04_refpick_extract_new`

1. **No thresholds.** Preprocessing rejects micrographs in the projects it hands on, and
   follows the GUI's threshold updates, which p04 does not receive. Both copies of the
   rule are gone, with the `reject_mics` and threshold defaults (nothing else reads them:
   not pick_extract, not its strategy). Before, a restart re-applied p04's default
   thresholds even with `reject_mics=no`, dropping micrographs, with their particles,
   that preprocessing had accepted after the user loosened a threshold.
2. **Restart history.** Each set records the upstream project it came from, and a
   restart puts those projects in the watcher history. Before, the history got p04's own
   completed set names; as the watcher compares file names, every gap in preprocessing's
   numbering made an upstream project be picked twice or skipped. Sets written before
   this change have no recorded origin: a restart warns, and their upstream projects are
   picked again once.
3. **Restart counter** comes from every completed set (job sets), also when none has a
   micrograph left.
4. **Restart folders.** The job folder is cleared and made again. Before, the restart
   removed `spprojs/` (and with `clear=yes` the completed, picker and extract folders)
   after they had been made, so the next set had nowhere to be written. With `clear=yes`
   all four are now emptied and kept.
5. **Waiting** for the upstream folder and for the picking references happens inside
   `iterate`, so SIGTERM and `TERM_STREAM` are honoured. There is no 24-hour give-up, and
   the wait no longer ends the process with `exit(0)`.
6. **Project write.** The stacks and particles are assembled from the imported sets in
   import order (one stack per micrograph, particle ranges renumbered), as before but per
   set instead of per micrograph record; the record list is gone.
   `sp_project%append_project` was not used: it wipes `os_cls2D`, `os_cls3D` and
   `os_out`, and logs three lines per appended set.
7. **Picking references to the GUI** go through `send_cavgs`: the stack path is absolute
   (was relative), resolution and population are not sent (were sent as 0), and the tile
   positions no longer divide by zero for one to three references.
8. **Pixel size** for make_pickrefs comes from the first set with an accepted micrograph.
   Before, it came from the first upstream project even with nothing selected, which
   would have run make_pickrefs with the default pixel size.
9. **Projects left.** `l_projects_left` now means the last watch was full, so more
   upstream projects may be waiting and the idle project write is held back. Before, it
   was true whenever a project had nothing selected.
10. **Parameters.** `moldiam` is removed before `params%new`. The worker and
    make_pickrefs command lines are built explicitly; the stage's command line is no
    longer changed after `params%new`. `split_mode` is set after `params%new` and the
    `ncunits` override is dropped, as in p01.
11. **Settle time** `SHORTWAIT` (2 s), as in p02 and p03.
12. **Signal handling** through `simple_stream_sigterm`.
13. **Not ported:** `update_user_params` (nothing writes its file), the unused
    `validate_ptcl2D_star_inputs`, the debug timers, the commented STAR writes, and
    `nmics_rejected_glob`, which was never reported.
14. **Optics groups.** Before each STAR export and project write, the newest map that
    optics assignment has published in `optics_dir` (passed by the master) is applied by
    import index (`import_latest_optics_map`): the micrographs, stacks and particles get its
    groups, the project gets its optics segment, and micrographs.star is exported with
    them. Before, the import was commented out and every micrograph was put in one optics
    group, which is still what happens until the first map exists. Micrographs newer than
    the map stay in group 1 until a later map covers them.
15. **No `pickrefs` output entry.** Finished sets are no longer given an `os_out` entry
    for the picking references: its only reader, p05, used one field (`mskdiam`), and now
    reads it from `moldiam.txt`, which make_pickrefs writes in this stage's directory. The
    running p05 still needs the entry, so p04 and p05 must be switched over together.

## Behaviour compared with the old `simple_stream_p05_sieve_cavgs_new`

1. **Mask diameter** from `moldiam.txt` in the reference-picking directory (`dir_target`),
   which make_pickrefs writes before any set is submitted. Before, it came from the
   `pickrefs` entry of the first imported set's output segment, which the running p04
   writes to a stray file (review finding 3), so the running p05 stopped at its first
   import.
2. **Optics groups.** The sieve is made with `optics_dir` (passed by the master, unused
   before), so the chunks it hands to pool 2D carry the newest optics map's groups. The
   chunks themselves stay without groups (see `ptcl_sieve` above). Before, the handed-off
   chunks had no optics groups.
3. **Waiting** for reference picking happens inside `iterate`, so SIGTERM and `TERM_STREAM`
   are honoured; the two 24-hour `wait_for_folder2` loops are gone.
4. **Settle time** `SHORTWAIT` (2 s): reference picking moves finished sets in with a
   rename (`stream_job_sets%complete`).
5. **Parameters.** The `lpstart` default, `workers = nchunks` and `worker_nthr = nthr` are
   set on the command line before `params%new`; before, they were set on `params` after
   it (`lpstart` by reading the command line again). The values are the same: the only
   derivation of `lpstart`, in `derive_image_settings`, is `max(lpstart, 2*smpd)` with no
   pixel size given, and `qsys_env` uses `nthr` when `worker_nthr` is 0. The derived
   `lplims2D` now follows `lpstart`; the sieve does not read it.
6. **GUI.** The pipe framing is `stream_pipe`, and the class averages go through the shared
   `send_cavgs` (now with optional resolution and population arrays). The latest class
   averages are only asked for once the sieve exists, and the final-ingestion flag is only
   set on a sieve that exists; an import still clears it.
7. **Signal handling** through `simple_stream_sigterm`.
8. **Restart** is as before: the sieve restores its chunks from its folders, and the
   projects it has chunked (its `imported_projects.txt`) go into the watcher history once
   the upstream folder is attached; projects imported but not yet chunked are imported
   again. The particle counters still start from 0 after a restart.
9. **Kept as it was:** the stage's own queue environment, which submits nothing but starts
   the persistent workers on the chunk partition that the sieve's queue environment then
   reuses.
10. **Not ported:** the unused read end `ipc_pipe_sieve_cavgs_out`, the commented
    `params%refs = params%pickrefs`, and the pipe and signal imports. `pickrefs` from the
    master stays unused.

## Behaviour compared with the old `simple_stream_p06_pool2D_new`

1. **No match-class selection.** Particles are no longer deselected at import by a GUI
   selection of sieve references (`class_match`); `update_match_class_states` was never
   called and was removed on 2026-10-05; the GUI no longer sends that selection.
2. **Default mask diameter** (no `mskdiam` given) is in Angstroms: the pixel count times
   the pixel size. Before, the pixel count was used as Angstroms.
3. **A GUI mask diameter** also drops the sieve's pending one. Before, the sieve's mask
   diameter replaced the user's at iteration 10.
4. **Waiting** for the sieve's folder happens inside `iterate`, so SIGTERM and `TERM_STREAM`
   are honoured, and GUI updates are drained meanwhile (a snapshot request waits for the
   pool). The two 24-hour `wait_for_folder2` loops are gone.
5. **Settle time** `SHORTWAIT` (2 s) instead of `LONGTIME` (60 s): the sieve now hands sets
   off with a rename.
6. **Restart.** The previous pool's files are removed in the stage's folder after
   `params%new`. With `dir_exec`, the old code removed it from the command line before
   `params%new`, so the stage did not restart in that folder, and cleaned the folder it
   was launched from instead. Snapshots are kept, as before; the pool starts again from
   every handed-off set.
7. **The pool is started once.** Before, a first import whose sets had no selected particle
   left the count at 0, and the next import initialised the pool again.
8. **Import** reads only the sets not yet in the pool, and only their `mic`, `stk`,
   `ptcl2D` and `out` segments. Before, every import allocated a project per set ever seen
   and read the non-data segments too.
9. **GUI.** The pipe framing is `stream_pipe`; the pool's and the snapshots' class averages
   go through the shared `send_cavgs`, with absolute paths (see the pool fixes above).
   `finalize` sends a last status with user input off.
10. **Signal handling** through `simple_stream_sigterm`. The final project is written by
    `terminate_stream2D` when the pool has started; before, it was called in any case,
    and without a pool it wrote nothing.
11. **Kept as they were:** the pause rule (now the `nopass` functions `pause_rate_factor`,
    `target_nptcls` and `runs_to_final`), `stepwise` read from the command line (it is not
    a parameter), the export sequence for 3D, and the project name `stream_solve2D` when
    no `projfile` is given (the GUI looks for it). The initial particle threshold is logged
    when it changes, not on every pass.
12. **Not ported:** the per-iteration snapshots (`L_ITERATION_SNAPSHOTS`, off),
    `extra_pause_iters` and `time_last_import` (set but read only by commented-out code),
    the commented pause rule, and the separate project creation before
    `create_stream_project`, which wrote the same file.

13. **Exports for 3D** continue their numbering after a restart (`restore_export_id`).
    Before, a restarted pool wrote `00001.simple` again; the 3D stage, which takes each
    export once by name, ignored it and every later export up to the old count, and lost
    the particles new in them. The restarted pool exports its stacks again; the new 3D stage
    skips the stacks it has. The exports are written atomically (see above). Changed since
    (2 October): the pool publishes its classified state after each completed iteration, and
    the 3D stage matches a publication's stacks to its rows by name.

### Known limitation carried over

- One stage per process: the pool's state is module state (step B: move it into a type).

## Behaviour compared with the old `simple_commanders_stream_p07_solve3D_multistate`

1. **A failed 3D job stops the stage** (`THROW_HARD` naming the job's log). Before, a crash
   never wrote `TASK_FINISHED`, and the stage waited forever with ingestion paused. The jobs
   run through `qsys_async_job`, which reads the exit status, in the same folders, scripts
   and logs as before.
2. **Mask diameter** from pool 2D's exports (see above); quality selection fits it to the
   class averages' box as the commander's parameters did.
3. **Quality selection** in place through `score_project_cavgs` with the `pool` model and
   `map_cavgs_selection`, instead of running the `model_cavgs_rejection` commander after a
   chdir into `quality_selection/<set>/` (which built its own parameters in process). Same
   selection; the selected and rejected stacks and their JPEGs are written as before. The
   commander's diagnostics (feature table, hard-gate and ranked stacks, quality annotation of
   the classes) and the copy of the project in that folder are not written; nothing read them.
4. **Duplicate stacks skipped.** A stack already in the pool is not imported again.
5. **The GUI gets every run's volumes.** Before, only solve3D's were sent; the addon runs'
   never reached the GUI. The status's state resolutions come from the latest run, read once
   per run instead of on every pass.
6. **The stage's project** holds the latest run's result after each run (written as a
   temporary file and renamed), and the particles imported since at `finalize`. Before, it
   stayed empty.
7. **Job settings** are named constants in the commander (3 states, 5 stages, low-pass 50 to
   10 Å, 8 parts, 8 threads, the values before), overridable with `nstates`, `nstages`,
   `lpstart`, `lpstop`, `nparts3D` and `nthr3D`.
8. **Restart** removes a leftover `TERM_STREAM`, which ended a restarted stage at once. As
   before, a restart imports every export again and runs solve3D again. Changed since
   (2 October): a restart takes the newest publication.
9. **Waiting** for pool 2D's folder happens inside `iterate`; the two 24-hour
   `wait_for_folder2` loops are gone. **Settle time** `SHORTWAIT` (2 s).
10. **GUI.** The pipe framing is `stream_pipe`: a frame is never abandoned once partly sent
    (the copy gave up after 3000 retries, mid-frame, and desynchronised the master's reader).
    The reprojection tiles go through `send_reproj_tiles`, the orientation histogram through
    `set_oridist_from_oris`. `finalize` sends a last status with user input off.
11. **Signal handling** through `simple_stream_sigterm`.
12. **Kept as they were:** the classes taken from the first export only (changed since, 2 October:
    the classes are the newest publication's), the preprocessing
    queue partition (`SIMPLE_STREAM_PREPROC_PARTITION`), the phase texts the GUI shows
    ("running refine3D" for the addon runs), a running job left running at termination.

## Behaviour compared with the old `simple_commanders_stream_p00_master`

1. **One stage table.** The same stages, folders, logs, GUI keys and command lines; the
   settings are named constants with the same values. Dropped: `pickrefs` for particle
   sieving and `nparts=1` for multistate 3D (neither stage reads them; `nparts` defaults to
   1), and the `ref_selection` update (no stage reads it since pool 2D's match selection
   went).
2. **Ctrl-C stops the stream in order**, as SIGTERM does. The SIGINT handler joined the
   listener thread and called `exit(1)` inside the handler, leaving the stages running.
   Both handlers only set a flag now (`simple_stream_sigterm`, `also_sigint`).
3. **A stage that does not stop is killed** (SIGKILL) `STOP_TIMEOUT_S` (10 minutes) after
   the stop; before, the master waited for it forever. The optics assignment still stops
   first, waited for up to a minute.
4. **While stopping, a heartbeat every 5 s**; before, back to back.
5. **Each GUI answer's updates are sent once.** Before, the update object kept every field
   it had received, so each later update sent an old snapshot request, mask diameter or
   reference selection again (the stages ignored repeats).
6. **Pipes** through `stream_pipe` on both sides. A corrupt frame length is reported and
   dropped (as the master did; the stages now do the same instead of stopping). A restart
   empties both of the stage's pipes under the listener's lock; before, the listener could
   be reading them meanwhile.
7. **The listener thread** is a `bind(c)` module procedure given the shared state's address;
   before, an internal procedure was passed to `pthread_create`, which standard Fortran
   does not allow.
8. **Metadata lists.** An item whose index is out of range is dropped with a warning
   (before, written out of bounds); the opening-2D final list is made with its own type
   (it was made with the non-final one); volume messages never carry reprojection tiles
   (`gui_metadata_vol3D%serialise`).
9. **One lock around the whole heartbeat assembly**, instead of seven lock/unlock pairs.
10. **Kept as they were:** the stop order, the persistent worker server and its warm-up, the
    memory log every 12 heartbeats, skipping preprocessing or the initial analysis when
    their outputs are given, the 10 ms pause of the listener, and `exit(EXIT_SUCCESS)` at
    the end.

## Switch-over

Done: production and the tests use the commanders here, the old stage modules of
`src/main/stream` are deleted (with their entries in the unused-function warning list of
`src/CMakeLists.txt`), and the three that carried `_new` have their plain names.

Remaining:

1. **Run:** the `stream_preproc` high-level entry, `lib_stream`, and one GUI session through
   to multistate 3D, including a stop, a stage restart and Ctrl-C.
2. **Register the testers:**
   - the mic selection, mic import, optics groups and optics maps testers in `unit_project`
     (the stage, job-set and master testers are in `unit_stream`, the pipe tester in
     `unit_ipc`);
   - `run_all_cavg_quality_selection_tests` in `unit_numerics`.

   Add each to its suite table in `simple_commanders_test_class` and to the suite's
   `suite=` help in the test UI. `scripts/check_test_registry.py` checks that the table
   and the help agree.
3. **Clean up** (the helpers only the old stages called, the `segdiampick_mics_multi` wrappers
   and `stream_export_optics` with its `assign_optics` / `h_clust` copy were deleted on
   2026-10-05, WS9 of `stream_fix_plan_2026-10-05.md`):
   - Switch `simple_stream_refine2D_utils` to `import_latest_optics_map`, and delete its
     `get_latest_optics_map` copies.
   - Move `simple_mic_selection`, `simple_mic_import`, `simple_optics_groups`,
     `simple_optics_maps`, `simple_stream_gui_senders` and `simple_stream_job_sets` to
     their homes, named in each module header.
   - Decide whether the batch `assign_optics_groups` should use `assign_optics_groups`
     (it needs its XML tilts on `os_mic` first).

None of this has been compiled or run yet.
