# Stream area review: state of the streaming pipeline, October 2026

Baseline: the working tree on 2 October 2026, which is `master` at e443948cd plus the
uncommitted stream refactor (`src/main/commanders/stream`, `src/main/stream/stages`,
`src/main/stream/master`). Line numbers refer to that tree.

Names: this is a dated document and keeps the names of 2 October. On 3 October `abinitio3D`
became `solve3D`, `abinitio2D` became `solve2D` and `cluster2D` became `refine2D` (the stream's
`abinitio2D_stream` program is now `pool2D`). In the code and the living docs, the multistate 3D
stage is `simple_stream_stage_solve3D` (p07 commander `simple_commanders_stream_p07_solve3D_multistate`),
`simple_stream_cluster2D_utils` is `simple_stream_refine2D_utils`, and the `abinitio3D_addon` route is
`solve3D_addon`. Line numbers predate the rename. Also on 3 October, the modules several stages
share (pipe, stream state, sigterm, job sets, GUI senders, meta plots, micrograph and optics
helpers) moved to `src/main/stream/shared`, and the 2D pool and chunk layer (POOL, C2D, ST2D and
the chunk modules) to `src/main/stream/pool2D`; the paths below are those of 2 October.

Scope: the stream commanders and stages, the master and its IPC, the 2D pool and chunk layers in
`src/main/stream`, the particle sieve and the class-average quality code the stream uses, the
stream STAR export, and the GUI metadata that crosses the pipes. Out of scope: the
`abinitio3D_addon` route of p07, which is under construction.

This review follows `stream_area_review_2026-09-30.md`; its findings are cited as S-A1, S-B1,
and so on. Findings marked *traced* were followed through the code but not reproduced; findings
marked *suspected* need a run to confirm. Everything else was read on the code path cited.

| Key | File |
|---|---|
| p00 | `src/main/commanders/stream/simple_commanders_stream_p00_master.f90` |
| p01 to p07 stage | `src/main/stream/stages/simple_stream_stage_{preprocess,optics,initial_analysis,refpick,sieve,pool2D,abinitio3D}.f90` |
| POOL | `src/main/stream/simple_stream_pool2D_utils.f90` |
| C2D | `src/main/stream/simple_stream_cluster2D_utils.f90` |
| ST2D | `src/main/stream/simple_stream2D_state.f90` |
| STAR | `src/main/star/simple_starproject_stream.f90` |
| SIEVE | `src/main/sieve/simple_ptcl_sieve.f90` |
| FEATS, MODEL, LEARN, SEL | `src/main/cavg_quality/simple_cavg_quality_{feats,model,learn,selection}.f90` |
| COMPAT | `src/main/class/simple_class_compatibility.f90` |
| ABI | `src/main/commanders/simple/simple_commanders_abinitio.f90` |

## Verdict

The structure has caught up with the library: every stage is now a type with a lifecycle and a
unit tester, the commanders are thin and live in `commanders/stream`, the seven copies of the
pipe protocol are one module, signal handling is safe everywhere, and five of the six September
runtime defects are fixed.

The science underneath has not changed, and on closer reading parts of it are worse than the
September review said. The picking references still come from a three-state ab initio
reconstruction of duplicated class averages. The pool sends particles to 3D that have never been
through a 2D iteration, and they skip the 3D quality gate. The rejection models were selected on
their own training data, are applied under conditions they were not trained in, and have no
recorded provenance. None of these choices has a policy note or a validation run behind it.

The 2D pool is still a set of module globals with roughly forty dead procedures, files and
histories that grow for the length of the run, and several crash paths. Three defects are serious
enough to fix before the next GUI session: a forked stage can hang for good (D1, since fixed in
the working tree and not yet run), a GUI selection of references never reaches reference picking
(D2, likewise fixed and not yet run), and the 3D stage ingests unclassified particles (D3, fixed
with a new ingestion contract, not yet run). The refactor itself introduced or carried over a few defects; they
are named where they occur.

## 1. Status of the September findings

| Finding | Status | Evidence |
|---|---|---|
| S-A1 `bestvol` uninitialised | Fixed | Falls back to the most populated state with a volume and stops when there is none (p03 stage:760-805) |
| S-A2 GUI selection path disconnected | Fixed in the working tree, 2 October; not yet run | See D2 |
| S-A3 p07 never learns the mask diameter | Fixed | p06's export records `mskdiam` (C2D:716-724); p07 reads it and fits it to the class-average box (p07 stage:296-312, 757-769) |
| S-A4 `increase_nmics` dead | Fixed in Fortran | NICE still sends it (`nice/nice_lite/data_structures/streamjob.py`); the master ignores it |
| S-A5 Unsafe signal handlers | Fixed | Every stage and the master only set a flag (`simple_stream_sigterm`, with SIGINT for the master) |
| S-A6 IPC is the raw memory image | Mitigated | `gui_metadata_vol3D%serialise` drops its allocatable tiles; `gui_metadata_project` has allocatables and no guard (not piped today); no round-trip tests; the rule is not written down |
| S-B1 Class-average duplication | **Open**, unchanged | See M1 |
| S-B2 Reference generation | **Open**, unchanged | See M2 |
| S-B3 3D ingestion contract | Replaced in the working tree, 2 October; not yet run | See D3 and M3 |
| S-B4 Rejection methodology | **Open** | See M5 |
| S-B5 Experimental behaviour in production | Partly | Gain-flip detection and generation are opt-in; `L_ITERATION_SNAPSHOTS` is gone; optics ids are still offset by the GUI display id (p06 stage:214); `stepwise` is forced only for standalone runs (`production/simple_stream.f90:75`) |
| S-C Stage shape | Fixed | Stage types with named steps and testers in `unit_stream` |
| S-C IPC copy-paste | Fixed | `simple_stream_pipe`; the master's store routes list items through four routines |
| S-C Resources | Partly | The master's settings are one table of named constants (p00:56-75); p03, the chunks and the pool keep literals; nothing is derived from the host |
| S-C `parameters` mutated after `new` | Partly | Commanders normalise before `new`; stages and the 2D layers still assign (G4) |
| S-C Imports | Open in the 2D layers | G3 |
| S-C Layering | Open | G2, G3 |
| S-C Dead code and naming | Partly | The old stage modules are deleted and the `_new` suffixes gone; about forty dead procedures remain (G7) |
| S-C Docs and skills | Partly | The skill and `streaming_pipeline.md` point at the new code. The stream policies now exist in `doc/policies/stream/`, stating the current contract (section 5); the review's proposals remain open |

## 2. Runtime defects, ranked

**D1. A forked stage can hang for good.** *Fixed in the working tree, 2 October; not yet
compiled or run.*

The defect, traced:
- A child of `fork` copies the parent's memory but not its threads. The master runs the memory
  monitor's sampler thread (p00:188; `memreport` defaults to `yes`, p00:131).
- In the child, `forked_process%start` calls `mem_monitor_init`
  (`src/utils/simple_forked_process.f90:118`). It sees the inherited `monitor_enabled` and calls
  `mem_monitor_finish` (`src/utils/simple_memory_monitor.f90:35`).
- That locks a mutex copied from the parent in whatever state it was in. It then `pthread_join`s
  the parent's sampler thread, which does not exist in the child
  (`src/fileio/simple_posix.c`, `simple_memory_monitor_stop_c`).
- The child never reaches its commander.
- With the default settings this hits every GUI restart. With `memreport=yes` passed to the master,
  `production/simple_stream.f90:60` starts a monitor before the master runs, and every stage
  hangs on its first start.
- Related: a restart forked while the listener thread could be writing to the log, so a child
  could also inherit the log unit locked.

The fix:
- `simple_posix.c` registers `pthread_atfork` handlers when the monitor first starts. They hold
  the monitor's mutex across the fork, and the child forgets the parent's monitor: it closes its
  copy of the file, clears `running` and reinitialises the condition.
- The master starts its listener thread after the first forks (p00:184) and forks a restart while
  holding the listener's lock (p00:264-267).
- The new test is `test_fork_with_running_monitor` in `simple_forked_process_tester`.

**D1a. A stage restarted from the GUI restarted itself on any failure.** *Fixed in the working
tree, 2 October; not yet run.*
- The master passed `restart=.true.` to `forked_process%start` for a GUI restart, meaning "it
  ran before". `forked_process` reads that flag as "restart on failure".
- From then on, any `is_running` poll that found the stage failed re-forked it, up to eleven
  times, outside the listener's lock. That included the SIGKILL of `force_stop`, so the stream
  could not stop such a stage.
- `stream_master_stage%start` now takes no argument and always forks with `restart=.false.`.
  Restarts are the GUI's call.

**D2. A GUI selection of picking references never reaches reference picking.** *Fixed in the
working tree, 2 October; not yet compiled or run.*

The defect, as found:
- **Wrong file name.** p03's `save_pickrefs_selection` wrote `STREAM_SELECTED_REFS//MRC_EXT`,
  that is `selected_references.mrc`. The master points reference picking at
  `selected_references.mrcs`, the name only the automatic route wrote. p03 then finished, and p04
  logged one line and waited indefinitely.
- **Stale comments.** The p03 header and the routine's comment said the file was
  `STREAM_DESELECTED_REFS`.
- **No precedence rule** between a user selection and a running 3D route. A selection arriving
  in the pass that collected the 3D result was never read.
- **A failed selection ended the stage anyway.** That happened when the cycle had no class
  averages, or no index was in range.
- **A restarted p03 ran its plan again,** and its 3D route would replace references p04 might
  already be picking with.
- **Partial files.** Both routes wrote the file in place, so p04 could read a partial stack.

The fix:
- **One name.** `OPENING2D_PICKREFS` (`simple_defs_stream`) names the file, for p03 and for the
  master (p00:384).
- **Published once.** p03 publishes once per run through `publish_pickrefs`: a complete stack
  renamed into place. Once published, the references are final.
- **The user wins.** Each pass reads the GUI's updates before the cycle steps, so a user selection
  pre-empts the 3D route. A selection that publishes nothing is logged and ignored.
- **Restarts keep them.** A restarted p03 that finds published references sends them to the GUI
  again and is finished at once (`restore_pickrefs`).
- **`nmics.txt` dropped.** The dead copy of it in p04 is gone, along with its constant.
- **Limitation.** A selection made after the 3D route has published cannot replace those
  references, because p03 has stopped. The GUI should stop offering a selection once p03 reports
  `terminating`.

**D3. The pool sends 3D particles that have never been classified.** *Fixed in the working tree,
2 October, by the new ingestion contract (M3); not yet run.* As found: each p06 pass imported new
sets, dispatches the next iteration, then exports (p06 stage:274-280). `iterate_pool` increments
the iteration at dispatch (POOL:258), so the export condition holds straight after every
dispatch (p06 stage:684-697) and the delta (`exported == 0`, C2D:665-700) includes the sets just
imported, with class 0 and no 2D parameters. In p07, `map_cavgs_selection` leaves class-0
particles untouched (`src/main/project/simple_sp_project_ptcl.f90:150`), so they keep their sieve
state and bypass the class-average quality gate. The `exported` flag is sticky, so nothing the
pool learns about them later reaches 3D. S-B3 called the exported labels provisional; for these
particles there are no labels at all. Now the pool publishes its classified state (stacks whose
particles have been through an iteration) right after an iteration comes back and before the next
dispatch; a fresh import is never published (`test_publication_holds_classified_stacks`).

**D4. A snapshot of a recent iteration can kill the pool.** *Fixed in the working tree, 3 October;
not yet run.* The history is a ring of the last `POOL_NHISTORY` (5) iterations. Tidying removes the
files of the iteration that leaves it (including `frcs_iterNNN.bin`), so history and files always
hold the same iterations. A request for an iteration no longer kept, or whose files are missing, is
logged and skipped: the GUI gets a snapshot report with 0 particles and no file. As found: the pool
kept a full project copy
per iteration for five minutes (POOL:215-257), but deletes each iteration's class averages after
five iterations (`tidy_2Dstream_iter`, C2D:377-388). A snapshot of an iteration in between copies
a deleted stack, and `simple_copy_file` stops the process (C2D:622-625).

**D5. A repeated snapshot request turns the next 3D export into a snapshot.** *Fixed in the
working tree. 2 October: publications no longer go through `write_project_stream2D`. 3 October:
the snapshot request reaches the new `write_pool_snapshot` as arguments; the module globals
`snapshot_iteration`, `snapshot_selection` and `snapshot_last_nptcls` are gone, with the dead
`update_user_params2D` that also set them.* As found:
`write_snapshot`
writes the global `snapshot_iteration` before checking that the request is new (p06 stage:658-659).
A repeated id leaves it set, and the next export, which also passes `snapshot_projfile`, takes
the snapshot branch (C2D:581): it writes a full snapshot non-atomically into the folder 3D
watches and sets no `exported` flags. Carried over from the old p06 by this session's port. Less
likely now that the master no longer re-sends earlier GUI fields, but still reachable when the
GUI repeats an id.

**D6. `particles2D.star` writes undefined image and micrograph names.** *Fixed in the working
tree, 2 October; not yet run.* `starfile_set_particles2D_table` tested and used `str_stk` and
`ind_in_stk`, which were never assigned in it; the lookup was removed in 4c004be6e. The final STAR
file of the 2D pool carried whatever was on the stack. The lookup
(`get_stkname_and_ind`) is back. The per-batch copy of the table, which could never run, is
deleted with its batch type. `simple_starproject_stream_tester` (`unit_project`, "STAR stream
export") checks each row's image and micrograph.

**D7. p07 sends stale FSC curves and minimum/maximum values.** *Fixed in the working tree,
2 October; not yet run.* `send_volumes` reused one `gui_metadata_vol3D` across states. `new` and
`kill` reset only the flags of the base type (`src/utils/gui/metadata/simple_gui_metadata_base.f90`),
so a state without an FSC file was sent with the previous state's curve. Introduced by this
session's refactor. Each state now starts from a default-initialised object, and
`test_send_volumes` sends a second state without an FSC curve. The other metadata types still
reset only their flags in `kill`.

**D8. A large GUI selection kills the master.** *Fixed in the working tree, 2 October; not yet
run.*
- **The defect.** The master's parser passed GUI arrays straight to setters that stop on
  overflow: `set_pickrefs_selection` above 500 entries, and `set_snapshot2D_update` above 1000.
  `set_coordinate` stops above 5000 picks. Both its callers, `send_recent_micrographs` and
  `gui_metadata_project`, could reach that limit, which would end a stage over a thumbnail.
- **Named capacities.** The metadata modules now export their capacities
  (`MAX_PICKREFS_SELECTION`, `MAX_SNAPSHOT2D_SELECTION`, `MAX_MIC_COORDINATES`). The wire format is
  unchanged.
- **Selections are dropped whole.** The master drops an oversized selection with a warning,
  because acting on part of it would act on classes the user did not choose.
- **Picks are truncated.** The two coordinate loops stop at the capacity.
- **Setters unchanged.** The setters keep their guards against programming errors.
- **Test.** `test_gui_commands_oversized` covers the parser.

**D9. Optics maps can be read half-written.** *Fixed in the working tree, 2 October; not yet run.*
- **The defect.** `write_optics_map` wrote `optics_map_<id>.txt` and then `optics_map_<id>.simple`.
  Readers choose the newest id by listing `*.txt` and then read both files. p03, p04, p05's
  hand-off and the pool all read maps while p02 publishes them.
- **Writer order.** The `.simple` is now written first and the `.txt` last.
- **Temporary names.** Each file goes to a name the `*.txt` listing cannot match
  (`<prefix>.tmp`, `<prefix>.txt.tmp`) and is then renamed.
- **One reader.** The pool's two private copies of the reader are replaced by
  `import_latest_optics_map` (`simple_optics_maps`), so the map files have one reader.
- **Test.** The optics stage tester checks that both files are published and no temporary file
  is left.

**D10. A restart sends the sieve's final-ingestion chunk to 2D classification.** *Fixed in the
working tree, 2 October; not yet run.*
- **By design.** Final ingestion puts the particles left below the coarse threshold into one
  coarse chunk flagged as already classified and rejected, so they skip the coarse pass. In
  two-tier mode they get the fine pass. In coarse-only mode (p03) the chunk is handed on as it
  is, and p03's cycle 2 classifies everything. The developers confirmed this design, which this
  review first reported as a defect.
- **The defect: flags lost on restart.** The flags lived only in memory.
  `import_existing_chunks_coarse` rebuilds each chunk's state from its marker files, and the
  staged chunk had none, so a restarted sieve submitted it for a coarse 2D run.
- **The defect: complete chunks resubmitted.** `submit` did not check `complete`, so even a chunk
  already merged into a fine chunk was run again.
- **The fix: a marker file.** The staged chunk now gets a `FINAL_INGESTION` marker
  (`FINAL_INGESTION_MARKER`). The import reads it as classified and rejection-complete, and
  cleanup keeps it.
- **The fix: a guard.** `submit` skips complete chunks.
- **Test.** `test_import_restores_final_ingestion_chunk` in the sieve tester.
- **Docs.** `ptcl_sieve_policy.md` now describes final ingestion (section 6.3).

**D11. The master never stops the persistent workers.** *Fixed in the working tree, 2 October;
not yet run.*
- **The defect.** p00 ended with `call exit(EXIT_SUCCESS)`, so the tail of
  `production/simple_stream.f90` never ran for the master. That tail does three things:
  - the persistent worker server's `kill`, which tells every worker to stop;
  - the log close;
  - the timer.
- **The fix.** p00 now returns from `execute` like every other commander. The forked stages still
  exit directly from `forked_process`, as they must: the server is the master's.
- **Test.** No unit test covers this; it is left to a stream run.

**D12. Smaller ones.** *Fixed in the working tree, 2 October; not yet run.*
- **Partition variable.** The sieve read the partition from `SIMPLE_CHUNK_PARTITION`, while p05's
  queue reads `SIMPLE_STREAM_CHUNK_PARTITION`. The sieve now reads `SIMPLE_STREAM_CHUNK_PARTITION`
  only; the old name is no longer read.
- **Heartbeat hashes.** A heartbeat answered with a non-200 status kept the section hashes, so
  what changed was not re-sent until it changed again. Now any heartbeat the GUI did not accept
  clears them.
- **Export FRCs.** Confirmed by reading: the pool's project has an `frc2D` entry only after a
  dimension change. So without downscaling, the export's `get_frcs` stopped the pool at its first
  export. The export now copies the pool's FRCs file, as the scaling branch does.
- **The final project's `ptcl3D`.** The final pool project was written with `ptcl3D` a raw copy
  of `ptcl2D`. It is now prepared as the STAR files have it before the write: 2D clustering
  removed, shifts kept, `updatecnt` and `sampled` removed.

## 3. Methodology

### M1. Class averages are duplicated to make `abinitio3D_cavgs` work (S-B1)

`balance_classes` is unchanged (p03 stage:1153-1269): the selected class averages and their even
and odd stacks are replicated in proportion to population up to `TARGET_NCLS = 501` rows
(p03 stage:132) and fed to `abinitio3D_cavgs`. Three effects the September review did not note:

- Every duplicate row carries its class's full `pop` (p03 stage:1224), so the populations are
  multiplied by the replication counts, not divided among the copies.
- The cycle-2 project on disk is rewritten with the balanced `os_cls2D` (p03 stage:1264): its
  particles' class indices point into a 501-row table they were never assigned to.
- `abinitio3D_cavgs` builds one row per class average with only `class`, `state`, `eo` and
  `stkind` (ABI:236-252). The copies are aligned and assigned to states independently, so one
  class average can vote in several of the three states, and the state chosen for the references
  ("most distinct projection directions", p03 stage:772-797) is counted over copies.

Duplication is integer weighting done by copying images. The weight belongs in the commander and
the reconstruction (R2).

### M2. The picking references are reprojections of a three-state model (S-B2)

Unchanged. `abinitio3D_cavgs` runs with `nstates=3`, `nstages=3`, `nrestarts_collapse=3`,
`lpstop=8`, `lpstart_ini3D=100`, `lpstop_ini3D=20` and `nthr=16` as literals
(p03 stage:1069-1091); the stage then hunts the highest numbered `<n>_abinitio3D_cavgs` directory
for the result (p03 stage:1313-1335) and reprojects the chosen state at `nspace=50`
(p03 stage:814). The 2D runs before it use `mskdiam=999` A whatever the particle size
(p03 stage:1057), with the picker's diameter estimate available. The precedence is now in the
code (D2: a user selection wins until references are published, and published references are
final), but there is still no policy note for the automatic and manual routes.

### M3. The 3D ingestion contract (S-B3)

*Replaced in the working tree, 2 October; not yet compiled or run.* The contract is now in
`doc/policies/stream/stream_3D_ingestion_policy.md`.

As found:
- The contract was deltas of not-yet-exported particles with a sticky flag, and D3 shows the
  deltas included particles with no 2D result.
- p07 ran the pool quality model per delta against the pool's class averages of that moment, and
  took its classes from the first delta only.
- The `exported` bookkeeping bit was a slot of the core particle parameter table (`I_EXPORTED =
  26`), which reused a retired slot. Because the retired slot was excluded from
  `oriparam_isthere` and the new one was not, every particle record of every workflow reported an
  `exported` key.

Now:
- **Publications.** After each completed iteration past 25, and before the next dispatch, p06
  publishes the pool's classified state: the stacks whose particles have been through an
  iteration, with current 2D parameters, `cls2D`, and class averages and FRCs at the native
  sampling. The 2 newest publications are kept.
- **Newest only.** p07 takes the newest publication when no 3D job runs, and makes one quality
  decision per publication.
- **Merge into append-only rows.** p07 merges the publication into rows that only grow, as
  `abinitio3D_addon` requires. Stacks are matched by name, and class, shifts and selection are
  updated in place with the 3D parameters kept. New stacks are appended; a stack missing from a
  publication is deselected. The classes are the publication's.
- **Start trigger.** `abinitio3D` starts on the first publication, which comes after the pool's
  low-pass ramp.
- **`I_EXPORTED`** is gone; slot 26 is `I_RETIRED_W` again, excluded from `oriparam_isthere`.

### M4. The 2D pool's schedule is undocumented

None of the following has a policy note, and most are literals:

- **Pause rule** (p06 stage:737-770): iterations 2 to 20 pause after more than one iteration
  without imports, later after one; a pause waits for 20 particles per class or 50 (500 after
  iteration 20) micrographs' worth. The final sieve set runs the pool uninterrupted to iteration
  25; exports start after iteration 25.
- **Mask diameter:** a box-based default until iteration 10, then the sieve's mask
  (p06 stage:80, 640-645). The classification changes its mask mid-run.
- **Sampling** (POOL:421-441): stacks are shuffled each iteration and taken until 500,000 selected
  particles, with no memory of earlier iterations and no coverage guarantee.
- **Fractional update** (POOL:368-382): *suspected* double sub-sampling, since cluster2D samples
  `update_frac` again inside the sampled stacks
  (`src/main/strategies/search/simple_matcher_smpl_and_lplims.f90:294`), giving about frac²; the
  out-of-sample populations (`prev_pop_even/odd`) are written and never read; `center` is
  switched off for good once frac falls below 1 (POOL:378).
- **New particles** get a random populated class before their first alignment, drawn inside an
  OpenMP loop (POOL:345-365), so runs are not reproducible; the low-pass and Gaussian ramp follows
  the global iteration (POOL:779-800), so particles arriving after iteration 20 never see the
  coarse limits.
- **Dimensions** (`update_pool_dims`, POOL:811-889): the box grows when the resolution sits at
  Nyquist, by Fourier padding the class averages and FRCs, which adds no information; the guard
  `ncls_glob < ncls_max` can never be true.
- **No class rejection in the pool:** `ncls_rejected_glob` stays 0 but is reported to the GUI;
  the only rejection is per export in 3D.
- **Optics group ids** of the STAR exports are offset by the GUI display id
  (`(nicedispid - 1) * 500`, p06 stage:214).

### M5. Rejection and class-average quality (S-B4)

- **The sieve scores with a mask diameter of 0** (SIEVE:1344, 1346) while the chunk's 2D runs
  with the real one from `moldiam.txt`. The features fall back to a box-sized disc (FEATS:211-212),
  the relation parameters to another one with a warning per chunk, and the mask-geometry gate
  ignores the particle size. The models' training conditions are not recorded, so whether this
  matches training cannot be checked.
- **Model selection uses the training data.** The learner fits every feature policy and penalty,
  then picks policy, penalty and threshold (560 candidates) by balanced accuracy on the same
  datasets (LEARN:104-136). There is no hold-out.
- **The feature bank holds exact negations.** `I_NEG_LOCVAR_FG` is `-I_LOCVAR_FG`
  (FEATS:255-262), robust z-scoring keeps the negation, and the sieve preset multiplies the two
  (interaction `[4, 10]`, MODEL:586-587), a hidden squared term.
- **Provenance:** the chunk and pool presets cite a developer's local path (MODEL:74-76); the
  sieve preset cites nothing.
- **Decisions are relative to the batch:** features are robust-z-scored per batch, then a fixed
  surface and threshold are applied, so a clean batch still loses its tail.
- **The compatibility filter fits and filters the same batch** with quantile bounds
  (COMPAT:38-42, 551-552), so tail classes are rejected by construction; it freezes after two
  stable fits (COMPAT:617-623) and works in pixels. p03 does the same per selection
  (p03 stage:1136-1141).
- **Tuning that is ignored:** the sieve uses its own 5,000/10,000 particle thresholds and walltime
  (SIEVE:93-96) and never reads `nptcls_coarse`, `nptcls_fine` or `walltime`, which its commander
  sets and its UI offers. All hard-gate thresholds are literals (FEATS:53-68, 290-343).
- Positive: class selection reaches the particles monotonically; a particle rejected in a coarse
  chunk cannot come back in a fine one.

## 4. Architecture

What is now in line with the library: stage types with named steps and testers; thin commanders;
the master split into a stage table, a metadata store and a GUI command parser; one framing
module; jobs that report failure (`qsys_async_job`); hand-offs that are written and renamed (the
sieve, the pool's exports, the job sets). The rest:

**G1. The pool is module state.** `simple_stream2D_state` holds about thirty-five public globals
(the pool project, its history, both cluster2D command lines, the chunk arrays, the iteration,
the snapshot request), written from POOL, C2D and the chunk utils and re-exported through
`simple_stream_api`. The p06 stage drives it through `get_pool_ptr` and a pointer component that
the pool reads as `master_cline` (p06 stage:87, 302). Nothing resets it, so a process can run one
pool, and nothing about it can be unit-tested (the p06 tester transfers sets into a local project
and skips the pool).

**G2. Library code runs commanders.** `cluster2D_utils` runs `commander_rank_cavgs` (C2D:9); the
p03 stage runs `commander_reproject` and the p04 stage `commander_make_pickrefs` in process;
`ptcl_sieve` still imports `commander_cluster_cavgs`, unused (SIEVE:63); and
`simple_stream_abinitio2D_chunks` is a commander living in `src/main/stream`. Each in-process
commander builds its own `parameters` in the middle of a stage.

**G3. Layering and imports.** Under `src/utils`, the GUI metadata, communicator and NICE modules,
the queue modules, `forked_process` and the memory monitor use `parameters`, `sp_project`,
`cmdline` or `image`; this session added `oris` to `gui_metadata_vol3D`. `sp_project` uses
`simple_gui_utils` with no `only`. The vendored json-fortran carries three local patches
(ab5d1d9c7, 6e53c3180, ce730c851). Most of the 2D layers `use` without `only:`, and
`simple_stream_utils` and `simple_mini_stream_utils` have no `private`, so they re-export what
they import.

**G4. `parameters` discipline.**
- After `new`, the stages still assign:
  - p01 the gain fields;
  - p03 a copy-and-edit `params_sieve` with literal resources (p03 stage:521-529);
  - p06 `smpd`, `box` and `mskdiam`;
  - p01 and p04 `split_mode`.
- The pool reassigns `lpstop` and then re-reads the command line for the value it overwrote
  (POOL:852).
- `ptcl_sieve` and `starproject_stream` each keep a full copy of `parameters`.
- Native `box` and `smpd` go on the pool and chunk command lines (POOL:573, chunk utils:282),
  against the project-metadata rule.
- `stepwise` is read from the command line and is not a parameter (p06 stage:202).

**G5. Encapsulation.** The stage types make every component and step public so their testers can
assemble them step by step. That is a deliberate exception to the rule that a stateful module
keeps its components private, and it needs a decision (keep it and say so in the skill, or test
through narrower interfaces). `ptcl_sieve` makes all 34 methods public, and its `kill` leaves the
tuning defaults, the mask diameter and the latest previews in place.

**G6. Interfaces.**
- `write_project_stream2D` was a routine of fifteen optional arguments and three modes, chosen by
  which optionals were present and by a global. It is now split (2 and 3 October): the final
  project (`write_project_stream2D`, three optionals), the snapshot (`write_pool_snapshot`, its
  request as arguments) and the publication for 3D (`publish_pool_state`).
- `starproject_stream` ignores its `outdir` argument.
- Several names drift:
  - the aliases `initial_analysis`, `opening2D` and `pool2D` in `production/simple_stream.f90` are
    unreachable, since the command-line parser rejects programs missing from the UI;
  - the GUI assembler's heartbeat section is `preprocessing` while the stage's GUI key is
    `preprocess`;
  - `ipc_pipe_abinitio3D_multstate` is misspelt.

**G7. Dead code.** Callers checked with `rg` over `src` and `production`:
- **Chunk utils:** the whole streaming-chunk path in `simple_stream_chunk2D_utils` (`analyze2D_new_chunks`,
  `memoize_chunks`, `update_chunks` and their getters).
- **C2D:** `update_user_params2D` (removed 3 October), `test_repick`, `write_repick_refs`, `set_dimensions`,
  `set_resolution_limits`.
- **POOL:** `update_match_class_states` (so the match-class path can never fire), `biased_stack_sampling`,
  `set_lpthres_type`, and the reference-generation branch.
- **`simple_stream_utils`:** `update_user_params`, `wait_for_folder`, `wait_for_folder2`,
  `process_selected_refs`, `process_selected_refs_2` (the last route into the legacy
  `stream_http_communicator`).
- **`simple_mini_stream_utils`:** the `segdiampick_mics_multi*` wrappers.
- **STAR:** `stream_export_optics` with its `assign_optics`/`h_clust` copy.
- **MODEL and LEARN:** the retired classify cache. `build_classify_cache` has no caller, so the
  cache routines in MODEL and the learner's cached paths can never run.
- **Never read:** the GUI loop for sieve reference selection (`ref_selection`,
  `initial_ref_selection`, the reference class averages the assembler expects). The copy of
  `nmics.txt`, which nothing wrote, was removed with D2.

**G8. Growth for the length of the run.** *Fixed in the working tree, 2 and 3 October; not yet run.*

As found:
- `pool_proj_history` gained a slot per iteration and deep-copied all live entries twice each
  iteration.
- `frcs_iterNNN.bin` and the exports' class-average stacks were never deleted.
- The watcher's history check was O(history) per file, inside an OpenMP loop with a shared flag,
  and the history grew one element at a time.
- p07's stack check was a linear search per stack.

Now:
- The history is a ring of 5 iterations, with one copy per iteration.
- `frcs_iterNNN.bin` is tidied with the other files of an iteration leaving the history.
- Only the 2 newest publications for 3D are kept.
- The watcher keeps a lexically sorted index of its history, searched by binary search, and grows
  it by doubling, with no OpenMP loop. A file added twice is kept once. Its tester is
  `simple_stream_watcher_tester` (`unit_stream`, "stream watcher").
- p07's stack match searches from the last match.

**G9. The UI and the reads disagree.**
- **Read but not exposed:**
  - `abinitio3D_stream`: `nstates`, `nstages`, `lpstart`, `lpstop`, `nparts3D`, `nthr3D`;
  - `abinitio2D_stream`: `stepwise`, `dynreslim`, `projfile_optics`, `optics_dir`;
  - `assign_optics`: `beamtilt`, `tilt_thres`;
  - the master: `beamtilt`, `tilt_thres`;
  - `preproc`: `eer_upsampling`;
  - `gen_pickrefs`: `pcontrast`.
- **Exposed but not read:**
  - `gen_pickrefs`: `nmics`, `optics_dir`;
  - `sieve_cavgs`: `ncls`, `nptcls_per_cls`, `nchunksperset`, `nptcls_fine`;
  - `preproc` offers `dir_prev` while the stage reads `dir_exec`.

**G10. Wire format.** Metadata still crosses the pipes as byte copies, which is safe only because
every stage is forked from the same binary and only while no type that crosses has an allocatable
component. `gui_metadata_project` has five and no guard. Fixed capacities stop the process on
overflow instead of truncating (D8).

**G11. Restart semantics differ per stage.** They are now documented together in
`doc/policies/stream/restart_policy.md`, as well as in each header.

| Stage | On restart |
|---|---|
| p01 | restores its job sets |
| p02 | continues the map ids |
| p03 | starts its plan again |
| p04 | restores its sets |
| p05 | restores its chunks |
| p06 | starts the pool again but numbers its exports on |
| p07 | starts again |

On restart, `cleanup_root_folder` (C2D:90-118) deletes every mrc, mrcs, txt, star, eps, jpeg, jpg,
bin and dat file in the pool's folder, along with the chunk folders.

## 5. Tests and documentation

**Tests.**
- `unit_stream` covers the steps of all seven stages, the job sets and the master's parts.
- `unit_ipc` covers the pipe framing.
- `unit_project` covers the sieve's lifecycle, its hand-off and one hard-gate rejection.
- `lib_stream` tests three helpers.
- The only high-level stream test is preprocessing.

Missing:
- an end-to-end test beyond preprocessing;
- any test of pool iterations, exports, snapshots or dimension changes;
- any rejection decision on realistic or simulated class averages;
- a serialise round-trip for every metadata type.

The fork under a running memory monitor (D1) now has a test, `test_fork_with_running_monitor`.

Five testers compile and are registered in no suite: micrograph selection, micrograph import,
optics groups, optics maps and class-average selection. One sieve test
(`test_new_accepts_tuning_overrides`) checks nothing about the overrides it sets.

**Documentation.**
- The stream policies are written (2 October), in `doc/policies/stream/`: reference generation,
  3D ingestion, the pool schedule, IPC, and restart. They state the current contract and list the
  known gaps; the review's proposals (R2 to R6) are named there as open.
- `model_cavgs_rejection.md` describes a microchunk engine that was deleted and a
  leave-one-dataset-out procedure that has no code.
- `ptcl_sieve_policy.md` says in one place that the model runs in the fine tier only and in
  another that it runs in both (the code: both), and does not mention the zero mask diameter.
  Final ingestion is now described (section 6.3, with D10).
- Links of the form `../../src` in the sieving policies resolve to `doc/src`.

## 6. Recommendations

Each item says what to change, where, and how to know it is done.

### R1. Fixes (days, independent)

1. **D1 and D1a, fork hang and self-restarts.** *Done in the working tree, see D1 and D1a.* It is
   done when `test_fork_with_running_monitor` passes and a GUI restart of a stage with
   `memreport=yes` reaches the stage's commander. What remains is structural: the master still
   forks while threads run (the listener, and the persistent-worker server's socket thread). The
   lasting fix is to start a stage by fork and `exec` of `simple_stream prg=<stage>`, passing the
   pipe descriptors on its command line, so the child inherits nothing from the master's threads.
2. **D2, reference selection.** *Done in the working tree, see D2.* It is done when the p03
   tester passes (`test_gui_selection_ends_stage`, `test_published_pickrefs_are_final`) and a GUI
   selection in a stream session starts reference picking.
3. **D3, unclassified particles.** *Done in the working tree with R3, see D3 and M3.*
4. **D4 and D5, snapshots.** *Done in the working tree, see D4 and D5.*
5. **D6 to D9.** *Done in the working tree, see D6 to D9.* What remains for D7: the other
   metadata types still reset only their flags in `kill`. Making `kill` reset every field in each
   type would make reuse safe everywhere.
6. **D10 and D11.** *Done in the working tree, see D10 and D11.*
7. **D12.** *Done in the working tree, see D12.*

### R2. Weights instead of duplication (M1)

- In `abinitio3D_cavgs`, where the class-average rows are built (ABI:236-252), set a per-row
  weight from `os_cls2D` `pop`, `w = pop / mean(pop of the selected classes)`, identical on the
  even and odd rows.
- Have the reconstruction honour the weight at insertion. If the cavgs reconstruction cannot read
  a per-row weight yet, adding it belongs in `src/main/volume` and the reconstruction strategy,
  not in stream code.
- Delete `balance_classes`, `duplicate_balanced_stack`, `update_os_out_stk` and `TARGET_NCLS`
  from p03, and the balancing tests from its tester.
- If the sampling needs a minimum number of rows, the commander enforces it with a clear error.
- Validation: one dataset, weighted against duplicated, comparing the state agreement, the
  even/odd FSC and the sieve's acceptance rate with the resulting references.

### R3. A 3D ingestion contract (M3), policy first

*Done in the working tree, 2 October; not yet run.* Decided as below, except that p07 keeps its own
append-only rows, matched by stack name, instead of mirroring the pool's index space: a pool
restart reorders the pool's rows, and `abinitio3D_addon` refuses renumbered rows. The start trigger
is the first publication. See M3 and the policy. Proposed:

- The pool publishes its current state at each export point, not a delta; p07 mirrors the pool's
  particle index space, updating `class`, `state` and shifts in place and appending new rows.
- One class-average quality decision per published pool revision, propagated to every row.
- `abinitio3D` starts on a stability trigger (the end of the low-pass ramp or the final sieve
  set), not on the first particles.
- Remove `I_EXPORTED` from `simple_defs_ori.f90`, which also removes the stray `exported` key
  from every particle record.

### R4. Reference generation (M2)

- Write `doc/policies/stream/reference_generation_policy.md`: which classes, the automatic and
  manual routes, their precedence, the output contract, and every parameter with its default.
- Make `nstates`, `nstages`, the low-pass limits and `nspace` inputs of `gen_pickrefs` (advanced
  visibility) instead of literals.
- Run p03's 2D with the mask diameter the picker estimated, not 999 A.
- Have `abinitio3D_cavgs` with restarts publish its final result at a fixed path, so callers stop
  scanning for the highest numbered directory.
- Log coverage and population for every state, and choose the state on view coverage together
  with the state agreement the cavgs route computes.

### R5. Rejection and quality (M5)

- Record each preset's provenance in its spec: datasets, commit and training conditions (pixel
  size, mask diameter, box).
- Score under the training conditions. In the sieve, either pass the real mask diameter or train
  the sieve preset with the box disc it uses, and say which.
- Select models by leave-one-dataset-out, as `model_cavgs_rejection.md` already claims.
- Drop exact-negation features, or reject a feature policy that contains a feature and its
  negation.
- Train the compatibility filter on earlier batches or the references and apply it to the next
  batch, rather than fitting and filtering the same classes.
- Make the hard-gate thresholds fields of the model spec, so they are versioned with the model.
- Either have the sieve read `nptcls_coarse`, `nptcls_fine` and `walltime`, or remove them from
  the commander and the UI.
- Add a regression test that pins each preset's decisions on simulated class averages with known
  good and junk classes.
- Rewrite `model_cavgs_rejection.md` and fix `ptcl_sieve_policy.md`.

### R6. The pool as a type (M4, G1, G6)

- A `stream_pool2D` type owns what `simple_stream2D_state` and the module variables of POOL and
  C2D hold today: the pool project, the history (bounded by iterations, not minutes), the
  dimensions, the queue environment and the command lines.
- The p06 stage holds one; `master_cline` and the pointer component go.
- Split `write_project_stream2D` into `write_pool_project`, `write_pool_snapshot` and
  `publish_pool_for_3D`, with an options type and the iteration as an argument.
- Write `doc/policies/stream/pool2D_policy.md` for the schedule in M4, and make its constants
  parameters.
- Decide and fix the fractional update (measure the effective fraction first), restore `center`
  after it or document why not, and seed new particles' classes deterministically.
- Delete the dead pool code (G7) and bound the per-iteration files (G8).
- Test exports, snapshots and dimension changes on a small pool.

### R7. Layering (G2, G3)

- No library module runs a commander: rank the class averages through a library routine, make
  the reprojection and `make_pickrefs` steps library calls or queued jobs, and drop the sieve's
  stale import.
- Move `simple_abinitio2D_chunks` to `commanders`.
- Move the GUI communicator, NICE and project metadata out of `src/utils/gui` into `src/main`,
  keeping only the metadata types and JSON helpers in utils; move `set_oridist_from_oris` out of
  the vol3D type into a `src/main` helper.
- Keep json-fortran pristine, with the fast printer in a wrapper module.

### R8. Parameters and UI (G4, G9)

- Register `stepwise`, or remove it together with the override in `production/simple_stream.f90`.
- Expose every parameter each stream program reads, and remove what it does not read (G9).
- Stop passing native `box` and `smpd` on the pool and chunk command lines.
- `ptcl_sieve` takes a small configuration type instead of a copy of `parameters`; p03 stops
  copy-and-editing.
- `starproject_stream` drops its copy of `parameters` and honours `outdir`.
- Values only known at run time (`smpd`, `box`, `mskdiam` at the pool's first import) are set
  in one helper that updates the command line and `parameters` together.

### R9. Resources

Derive every stage's `nthr`, `nparts` and `nchunks` from `params%nthr` and the existing
environment overrides in one place, log the budget at start-up, and pass it on. The master's
table of named constants is the first half of this; p03, the chunks and the pool still use
literals.

### R10. Tests and documentation

- Register the five unregistered testers.
- Add a high-level, master-less test that chains p01 to p07 on simulated data with the
  simulation framework the self-contained tests already use.
- Add a serialise round-trip test for every metadata type and a test that feeds the master's parser
  oversized arrays.
- The policies are written as the current contract (`doc/policies/stream/`). When R2 to R6 are
  decided, the policy of each is updated with the code.

### R11. Process

As in September (S-E12): stream work lands from branches through review with `unit_stream`,
`unit_ipc` and the sieve suite green; experimental methodology ships behind a parameter whose
default is the validated behaviour, with a policy note; commit messages state what changed and
why.

### Suggested order

1. R1 is done in the working tree (D1 to D12, D1a), with G8; what remains is compiling and running
   it.
2. R2 and R4: the methodological decisions the rest depends on (R3 is done in the working tree).
3. R5, with its regression test before any retraining.
4. R6: the pool as a type (the project writer is already split, and the history bounded).
5. R7 to R9 as each area is touched; R10 alongside; R11 from now.
