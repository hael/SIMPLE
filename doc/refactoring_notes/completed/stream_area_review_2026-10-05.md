# Stream area review: state of the streaming pipeline, 5 October 2026

Baseline: `master` at 0cb34fedf plus the uncommitted p03 changes of 3 October (the mask diameter
estimated from cycle 1's class averages, and the 3D state choice on connected components) and
the helpers they use (`automask2D_mskdiam`, `vol_shape_descr`). Line numbers refer to that tree.

Since the 2 October review the stream refactor has been committed (61fd9b9ee, "Tests run, stream
runs, testing underway"). Two later commits changed the stream code for compile time: stage
`parameters` and the sieve are allocatable state, and the copies of `parameters` in `ptcl_sieve`
and `starproject_stream` are gone (93e855ada, 8aedefec5). Dead pool code was removed (3a4d57819),
and the notes moved into `doc/refactoring_notes/{completed,planned}/` (b818cd7cc).

Scope: the master and its IPC, the stage commanders and types p01 to p07, the shared modules
(`src/main/stream/shared`), the 2D pool and chunk layer (`src/main/stream/pool2D`), the watcher
and stream utilities, the particle sieve and the class-average quality code it calls, the
stream STAR export, the GUI metadata that crosses the pipes, and the UI entries of the stream
programs. New this time: p07's use of `solve3D_addon` from the stream side (how p07 starts,
feeds and collects addon runs). Out of scope: the internals of the batch `solve3D_addon`
program, and NICE, which was read only where the master's contract depends on it.

Method: six reviews, one per area, read the code; every finding ranked High was then checked
again on its code path. Findings marked *traced* were followed through the code but not
reproduced; *suspected* ones need a run. Everything else was read on the path cited. Nothing
was compiled or run.

Numbering continues the 2 October review: new defects start at D13, methodology at M6,
architecture at G12, recommendations at R12. Items of that review keep their numbers.

| Key | File |
|---|---|
| p00 | `src/main/commanders/stream/simple_commanders_stream_p00_master.f90` |
| p01 to p07 stage | `src/main/stream/stages/simple_stream_stage_{preprocess,optics,initial_analysis,refpick,sieve,pool2D,solve3D}.f90` |
| p01 to p07 cmd | `src/main/commanders/stream/simple_commanders_stream_p0N_*.f90` |
| MSTAGE, GUICMD, STORE | `src/main/stream/master/simple_stream_master_{stage,gui_commands,meta_store}.f90` |
| PIPE | `src/main/stream/shared/simple_stream_pipe.f90` |
| OG, JS | `src/main/stream/shared/simple_{optics_groups,stream_job_sets}.f90` |
| WATCH, SU, MSU | `src/main/stream/simple_{stream_watcher,stream_utils,mini_stream_utils}.f90` |
| POOL, R2D, ST2D | `src/main/stream/pool2D/simple_stream_{pool2D_utils,refine2D_utils}.f90`, `simple_stream2D_state.f90` |
| STAR | `src/main/star/simple_starproject_stream.f90` |
| SIEVE, PU | `src/main/sieve/simple_ptcl_sieve.f90`, `src/fileio/simple_projfile_utils.f90` |
| FEATS, MODEL, LEARN | `src/main/cavg_quality/simple_cavg_quality_{feats,model,learn}.f90` |
| PK, PS | `src/main/pick/simple_segdiam_bin_picker.f90`, `src/main/strategies/parallelization/simple_pick_strategy.f90` |
| FORK, QENV, AJOB | `src/utils/simple_forked_process.f90`, `src/utils/qsys/simple_qsys_{env,async_job}.f90` |
| UI | `src/main/ui/simple/simple_ui_stream.f90` |

## Verdict

The defects the 2 October review ranked first are fixed in the code: the fork hang, GUI-only
restarts, the reference selection, the 3D ingestion contract, the snapshot defects, the STAR
names, stale 3D GUI data, oversized selections, optics maps read half-written, the final
ingestion chunk and the persistent workers (D1 to D12), and the growth of files and histories
(G8). The stages own their parameters, the sieve no longer snapshots `parameters`, and the pool's
history is a bounded ring.

The new findings cluster around four problems that cross stage boundaries, which the
stage-by-stage refactor did not reach:

1. **Jobs outlive the stage that started them.** No stage cancels or adopts its queued jobs when
   it stops or restarts, and most restart into the same folders. p07 starts a second 3D job
   beside the first (D21); the pool can take a stale finished marker for its first iteration
   (D26); p01 and p05 reuse set numbers and chunk folders (D26).
2. **Queued jobs have no failure path,** except where `qsys_async_job` is used (p03, p07). A sieve
   chunk killed by the scheduler holds its slot for good and blocks the final hand-off (D17); a
   failed pool iteration stalls p06 for good (D27).
3. **The master and NICE disagree about updates.** NICE re-sends its fields and restart flags on
   every heartbeat; the master forwards every answer to every stage, so a partly written frame
   can block the master for good (D13), and a restart flag that stays set restarts a stage over
   and over and can stop a stop from ending (D22).
4. **Stages are forked without exec,** so they inherit the master's persistent-worker server
   object (D14), its unflushed log buffers and its process group.

Restarts do not keep the semantics `restart_policy.md` sets out in p01, p03, p05, p06 and p07
(D18, D19, D25, D26, D29). Some defects come from the work of the last days: the GUI route D2
opened reads a cycle 2 selection against the balanced stack and publishes the wrong images
(D15), and the 3 October state choice can prefer a compact junk volume (M6). p02's regrouping
from scratch on every pass, carried over from the old stage, caps the length of a session (D16).

The science is where it was: duplicated class averages (M1), a three-state model with literal
settings (M2), rejection models applied outside their training conditions with no provenance
(M5). The addon route adds open questions of its own: when to run it, when to rebase, and what
to do with its verdict (M8).

Fix before the next GUI session: D13, D14, D15, D17 to D21, D22, D33 and D35 (R12).

## 1. Status of the 2 October findings

| Item | Status | Evidence and what remains |
|---|---|---|
| D1 fork hang | Fixed | `pthread_atfork` for the memory monitor (`simple_posix.c`:1116-1160); listener after the first forks (p00:184); restart forks under the lock (p00:268-271). The persistent-worker server thread still predates the forks (D14), so the comment at p00:158-160 is wrong. |
| D1a self-restarts | Fixed | MSTAGE:105-106 forks with `restart=.false.`. The same loop now comes back through the GUI's restart flag (D22). |
| D2 reference selection | Fixed | One name both sides (p00:388, p03 stage:933); publish once by rename (p03 stage:927-942). New hole: D15. |
| D3 unclassified particles | Partly | Publications precede dispatch, once per iteration past 25 (p06 stage:282, 691-707). Stacks are published whole when any particle was classified, which lets never-aligned particles through under fractional update (D35). |
| D4, D5 snapshots | Fixed | Ring of 5 (POOL:214-224) with matching tidying (POOL:384, R2D:289-300); request as arguments (p06 stage:674-675). One fatal path left (D41). |
| D6 STAR names | Fixed | `get_stkname_and_ind` (STAR:259); tester registered. |
| D7 stale vol3D | Fixed | Fresh object per state (p07 stage:747-750); `test_send_volumes`. |
| D8 large selection | Fixed | Dropped whole (GUICMD:78-85, 97-103); tested. |
| D9 optics maps half-written | Fixed | Temporary and rename (`simple_sp_project_io.f90`:1273-1296); one reader (`simple_optics_maps.f90`:59-65). |
| D10 final-ingestion chunk | Fixed | Marker written and restored (SIEVE:808, 344-349); crash window between project and marker (D41). |
| D11 persistent workers | Fixed | p00 returns; entry point kills the server (`production/simple_stream.f90`:84-86). Untested. |
| D12 smaller ones | Fixed | Heartbeat hashes (p00:206-216), the export FRCs (R2D:616-619), the final `ptcl3D` (R2D:351-354). The sieve reads the partition variable, but only into the chunk projects; its own queue ignores it (D41). |
| M1 duplication | Open | `balance_classes` unchanged (p03 stage:1244-1360); each copy keeps its class's full population. |
| M2 reference model | Partly | Real mask diameters (p03 stage:625, 667, 1071); state choice replaced (M6); settings still literals (p03 stage:140, 827, 1090-1100); still scans for the highest-numbered result folder. |
| M3 3D ingestion | Partly | Implemented as decided, with deviations: the merge overwrites 3D state labels (D33), rows matched by offset not `indstk` (p07 stage:475-481), addon trigger on rows not cohort (D34). |
| M4 pool schedule | Partly | Documented and constants named; double sub-sampling, `prev_pop_*` never read, random classes drawn inside OpenMP (POOL:318-349) remain. |
| M5 rejection | Open | Sieve scores with `mskdiam=0` (SIEVE:1334, 1336); in-sample model selection (LEARN:107-127); feature-negation interactions (MODEL:495, 515); no provenance; fit and filter on the same batch (SIEVE:1350-1383). |
| G1 pool as module state | Open | `simple_stream2D_state` globals re-exported (`simple_stream_api`:7); `kill` resets none of it (p06 stage:306-359). |
| G2 library runs commanders | Open | R2D:7 (`commander_rank_cavgs`); SIEVE:29 imports `commander_cluster_cavgs` and never uses it. |
| G3 layering | Unchanged | `src/utils` GUI metadata and qsys modules still depend on `src/main` (`parameters`, `sp_project`). Not re-reviewed. |
| G4 parameters discipline | Open | p03 copies `params` whole for the sieve (p03 stage:540-551) and fills a local `parameters` by hand (p03 stage:1180); p04 sets `split_mode`, `smpd`, `box` after `new` (p04 stage:190, 454, 482); p01 sets gain fields, thresholds and `split_mode` (p01 stage:202, 382-387, 814-821); p05 sets `mskdiam` (p05 stage:338); p06 sets `smpd`, `box`, `mskdiam` (p06 stage:573-577) and re-reads `lpstop` from the command line (POOL:769-770). |
| G5 public components | Unchanged | A decision is still owed (keep and document, or narrow the testers' access). |
| G6 interfaces | Partly | Writer split done; `starfile_init` ignores `outdir` (STAR:57-76); unreachable aliases in `production/simple_stream.f90`:68; misspelt pipe name fixed. |
| G7 dead code | Partly | See G7 below. |
| G8 growth | Fixed | Bounded ring, tidying, two publications kept, watcher index. p07 now grows without bound (D39). |
| G9 UI vs reads | Open, wider | See G9 below. |
| G10 wire format | Partly | `vol3D` serialises its own way; `gui_metadata_project` still has allocatable components; `kill` resets only flags; the store copies frames without checking their length (D41). |
| G11 restart semantics | Partly | Documented in `restart_policy.md`; violated in five stages (§2). |
| R1 fixes | Done but one | All but fork+exec (G13). |
| R2 weights | Open | |
| R3 ingestion contract | Done as decided | Stability trigger still open (ingestion policy:120-121). |
| R4 reference generation | Partly | Policy written; mask estimate; states logged. Inputs, fixed result path, state agreement not done. |
| R5 rejection | Open | |
| R6 pool as a type | Partly | Writer split, policy, bounded history. No type; iterations untested. |
| R7 layering | Open | `simple_stream_solve2D_chunks` is still a commander in `src/main/stream`. |
| R8 parameters and UI | Partly | Copies dropped in the sieve and STAR export; `stepwise` still unregistered (p06 stage:207-208). |
| R9 resources | Open | Literal threads and parts in p03, the sieve and the pool. |
| R10 tests and docs | Partly | Stage testers registered; four shared testers still not (§5); no chained master-less test; no serialise round trip. |

## 2. Runtime defects, ranked

### High

**D13. A partly written update frame blocks the master for good.** *Traced.*
- NICE returns its whole `master_update` dictionary on every heartbeat
  (`nice/nice_lite/data_structures/job.py`:111-120); the threshold and mask fields are never
  removed (`streamjob.py`:151-157, 460).
- The master turns every answer into an update (GUICMD:52-53) and sends it, about 5.6 KB, to
  every running p01, p03 and p06 every 5 s, with no deduplication (p00:273-278).
- `stream_pipe%send` drops a frame only while nothing of it has been written (PIPE:113). Once part
  of a frame is in the pipe it retries every 10 ms with no limit (PIPE:102-124). The master never
  closes its copy of a stage's read end, so a stage that has exited never turns the write into
  `EPIPE`.
- Scenario: at stop, p06 needs more than about a minute for its last pass and reads no updates.
  About twelve frames fill its 64 KiB pipe and the thirteenth is written in part. The stage then
  exits, and the master's main thread spins for good: no heartbeat, no SIGTERM handling, no exit.
  Short of that, a full pipe stalls the main loop for up to 2 s per stage per pass.
- Fix: send a stage only what differs from the last update it received, send nothing to a
  stopping stage, and give up on a partial frame after a time limit (then discard the channel
  before the next restart).

**D14. Forked stages inherit the master's persistent-worker server.** *Traced.* High when the
environment selects a `*_worker` queue system.
- The master starts the server before the forks (p00:148). The server is a module singleton, and
  `is_running` is just `port > 0`, so every child sees a live server.
- A stage's `qsys_env%new` then takes the reuse branch (QENV:266-271):
  - p03 asks for 32 threads per worker against the inherited 16 and stops at QENV:269 ("cannot
    reuse existing worker server with lower nthr_per_worker"), with its queue built at p03
    stage:264-275;
  - the other stages queue through the master's own client connection, shared by every child, so
    replies can be crossed or tasks lost;
  - `set_warmup_cooldown_enabled` locks a mutex copied in whatever state the master's thread held
    it at the fork.
- The `worker_server` address the master passes to p05, p06 and p07 (p00:404, 420, 432) is never
  used.
- Fix: in the child, before `execute`, drop the inherited server object without killing it and
  connect as a client through `worker_server`. The lasting fix is fork+exec (G13).

**D15. A GUI selection of cycle 2 class averages publishes the wrong images.** *Traced.*
- Cycle 2's selection and its class balancing run in the same pass (p03 stage:679-689).
  Balancing points the project's class averages at the 501-row duplicated stack (p03
  stage:1326-1330).
- The GUI shows cycle 2's classes with their original indices (`simple_stream_gui_senders.f90`:54-57).
  A selection arrives in a later pass, and `save_pickrefs_selection` reads the indices against
  whatever stack the cycle 2 project names (p03 stage:895-916), which by then is the balanced one.
- Scenario: the user picks classes 5, 12 and 30; the references are copies of the first selected
  classes, published as final, and p04 picks the whole session with them.
- Fix: balance into a separate project (for example `all_balanced.simple`) for `solve3D_cavgs`,
  or keep each cycle's selected stack and read selections against it. Add a test that selects
  after balancing.

**D16. p02 regroups every micrograph on every pass with a full distance matrix.** *Traced.*
Carried over from the old stage.
- Each pass that imports anything calls `assign_optics_groups` over all micrographs (p02
  stage:146-148, 233-243), which clusters each tilt group with `h_clust`. That allocates an N×N
  matrix and merges with a `where` over all labels for every close pair of centroids
  (`simple_starproject_utils.f90`:216-273).
- Without `dir_meta`, p01 sets every shift to 0 (p01 stage:591-593), so all centroids coincide and
  the merge step does work of order N³.
- At 20,000 micrographs the matrix alone is 1.6 GB; at 50,000 it is 10 GB, on the node that runs
  the master and every stage. p02 checks no stop request while it clusters.
- Fix: assign new micrographs to the nearest existing group (or a new one) and keep the groups;
  skip clustering when all shifts are equal.

**D17. Sieve chunk jobs have no failure path.** *Read and traced.*
- Nothing writes `REJECTION_FAILED`: it is defined (SIEVE:77) and read on restore (SIEVE:347,
  391), and only the tester creates it.
- A chunk leaves "running" only when `SOLVE2D_FINISHED` appears (SIEVE:1185-1191, 1258-1264); the
  jobs are submitted without an exit-code file (SIEVE:1139, 1157).
- A coarse `solve2D` killed by walltime, node loss or a crash keeps its slot for good (SIEVE:441-451)
  and blocks the final flush, which needs every coarse chunk complete (SIEVE:909). The final set
  never reaches the pool.
- An error inside `reject_cavgs` stops p05; after a restart the chunk has `SOLVE2D_FINISHED` but
  no rejection marker, so rejection runs and stops p05 again.
- Fix: pass an exit-code file and poll it as `qsys_async_job` does; catch rejection failures and
  mark the chunk failed and complete.

**D18. A p05 restart drops the unchunked micrographs of partly chunked sets.** *Traced.*
- Records are per micrograph (SU:87-97) and chunk boundaries fall between micrographs (SIEVE:732-735,
  752), but `imported_projects.txt` lists every set with at least one chunked micrograph
  (SIEVE:769-773).
- On restart those sets go into the watcher history (p05 stage:287-291), so their remaining
  micrographs are never imported. This contradicts `restart_policy.md` §3 and the stage header
  (p05 stage:28-30).
- Fix: record chunked micrographs (`projname`, `micind`), or rebuild the included set from the
  chunk projects, and re-import partial sets with their chunked micrographs marked.

**D19. A p05 restart with nothing new to import never resumes the sieve.** *Traced.*
- The sieve is made only when `project_list` is not empty (p05 stage:201-207), and that list is
  filled only by new imports.
- After a restart in which every record was already chunked (always the case once final ingestion
  has staged the leftovers), or after p04 has finished, `ptcl_sieve%new` never runs. Restored fine
  chunks are never handed off and the final chunk never reaches the pool.
- Fix: on restart, make the sieve as soon as the upstream folder is attached.

**D20. A final flush of only the staged chunk writes an uninitialised FRC object.** *Traced;
crash suspected.*
- The staged chunk has no class averages, so `merge_chunk_projfiles_without_sigma2` never builds
  `frcs` (PU:455-461), yet writes it unconditionally (PU:507). `merge_chunk_projfiles` guards the
  same write (PU:214-217).
- The path is common: the fine loop usually consumes every eligible coarse chunk, so at final
  ingestion the flush often merges only the staged chunk (SIEVE:838-855, 907-939).
- If the write crashes, the folder is left without a project and the flush repeats on restart.
- Fix: guard the write and `add_frcs2os_out` with `frcs_initialised`; add a staged-only flush test.

**D21. A p07 restart runs a second 3D job in the folder of the first.** *Traced.*
- `finalize` leaves a running job to finish (p07 stage:251-258) and `kill` does not stop it; a
  local job is a `nohup … &` process (`simple_qsys_ctrl.f90`:764).
- The restarted p07 takes the newest publication, writes `solve3D/solve3D.simple` over the old
  job's project, and submits a second `solve3D` in the same folder (p07 stage:594-622).
- When the old job ends it writes `EXIT_CODE_solve3D=0`; the new p07 then takes its own running job
  as done and adopts a half-written project, which crashes `finish_run` or makes the next addon
  refuse its frozen project.
- Fix: stop the job in `finalize`, or record it and have a restart wait for or adopt it; give each
  run its own folder and never reuse one without an exit code.

### Medium

**D22. Restart flags that NICE keeps set restart a stage over and over, and can stop a stop
from ending.** *Traced.*
- NICE removes `restart_<key>` only when a heartbeat reports the stage running
  (`streamjob.py`:246-286); p00 restarts any stage that is not running whenever the flag is there
  (p00:265-272), on every pass and during a stop.
- A stage that fails within a heartbeat of starting is forked again on every pass (the D1a loop,
  through the GUI); each start repeats its restart steps.
- A `terminate` together with a pending restart forks the stage and signals it at once
  (p00:295-300); it dies by the default handler before installing its own (FORK:105-106) and is
  reported failed, NICE keeps the flag, and the master never reaches its last loop.
- Fix: apply no restart once stopping; act on a restart request once until the key disappears;
  ask NICE to clear the key once delivered.

**D23. A GUI restart of a skipped stage runs it.** *Traced.*
- `is_running` is false for a skipped stage (MSTAGE:115-118) and `forked_process%start` clears
  `skipped` (FORK:94).
- With `dir_preprocess`, a restart of p01 runs it in the symlinked folder of the user's earlier
  run and rewrites it; with `pickrefs`, a restart runs p03's whole plan for references nobody
  reads.
- Fix: keep `skipped` sticky and ignore restarts of skipped stages.

**D24. The `dir_preprocess` link name has no C terminator, and a master restart fails on the
existing link.** *Read; effect of the first suspected.*
- `symlink` gets `PREPROC_JOB_NAME` without `achar(0)` (p00:162), unlike its first argument.
- A second master start in the same folder stops at p00:163 (`EEXIST`), against
  `restart_policy.md`:13-15.
- Fix: terminate the name; accept an existing link to the same target.

**D25. A p01 restart reprocesses rejected movies and reuses import indices.** *Read.*
- On restore, the import counter becomes the number of accepted micrographs and only accepted
  movies enter the watcher history (p01 stage:356-372), while every movie consumes an index
  (p01 stage:589-590).
- Rejected movies come back from the watcher and are processed again; new micrographs get indices
  already handed downstream. p02 then holds duplicate rows, and the optics maps carry repeated
  `importind` keys, of which the last wins (`simple_sp_project_optics.f90`:38-47): old accepted
  micrographs take a new micrograph's optics group in p04, the sieve and the pool.
- The tester asserts the counter as it is (`simple_stream_stage_preprocess_tester.f90`:164).
- Fix: restore from every micrograph of the completed sets, whatever its state; set the counter to
  the highest `importind`.

**D26. Stopped or crashed stages leave their jobs running, and restarts reuse their folders.**
*Traced or suspected per backend.*
- **p01**: `finalize` only cleans up files (p01 stage:289-295); a restart numbers sets from the
  completed folder and removes the job folder (JS:176-198). A set still queued (say 13, with 12
  the highest completed) is written again under the same number while the old job may still
  read it.
- **p05**: a restart regenerates and resubmits any chunk without `SOLVE2D_FINISHED` into the same
  folder (SIEVE:345-395), although its job may still run on a persistent worker, SLURM or in the
  background.
- **p06**: the pool iteration runs as a background job and outlives a crashed p06; it stops early
  only on `TERM_STREAM`. The restart cleanup keeps `REFINE2D_FINISHED` (p06 stage:227-234;
  R2D:88-116), so the new pool takes iteration 1 as complete at once, then stops copying a
  `frcs.bin` the cleanup deleted; if the old job finishes later, its project is mapped through the
  new pool's stack mask.
- Fix: one job lifecycle for every queued job (R13).

**D27. A queued job killed by the scheduler is never seen to end.** *Traced.*
- The pool is free again only when `REFINE2D_FINISHED` appears (POOL:368, 372, 827-835); a failed
  iteration leaves p06 "finding and classifying particles" for good.
- `qsys_async_job` (p03, p07) reads an exit code the script writes after the program
  (AJOB:86); a walltime kill or lost node writes none, so p07 waits for good.
- `solve2D_chunks` waits for `SOLVE2D_FINISHED` with no failure path
  (`simple_stream_solve2D_chunks.f90`:93-99; `simple_stream_chunk.f90`:390-403).
- Fix: exit codes plus a timeout on the log's modification time, reported in the status (R13).

**D28. p03's cycle 2 can start before the last particles reach its sieve.** *Traced.*
- In one pass p03 cycles the sieve, then sets final ingestion, then moves to `ALL_SIEVE` (p03
  stage:484-491); `run_cycle2` in the same pass checks `get_finished` and combines and kills the
  sieve (p03 stage:650-656). The leftover chunk is only made in the next `cycle` (SIEVE:785-821),
  and `get_finished` looks only at existing coarse chunks (SIEVE:584-590).
- Up to 5,000 leftover particles are then never used, against the policy's "once the sieve has
  every particle".
- Fix: cycle the sieve once more after setting final ingestion, or require no unchunked record
  before combining.

**D29. A p03 restart before publication mixes old and new sieve state.** *Traced.*
- Local copies are numbered from `00001` again over the previous run's files (p03 stage:446,
  1000-1005), while the new sieve takes up the old run's chunks (p03 stage:553-554; SIEVE:237-358).
- If the previous run was still collecting, every micrograph is extracted again under the same
  stack names and fed again: each particle reaches cycle 2 twice, and restored chunk jobs read
  stacks being rewritten.
- If it had already combined, combining does nothing because `all.simple` exists (PU:1057), and
  cycle 2 reruns on the old, already selected and balanced project.
- `reference_generation_policy.md` §3.1 item 6 ("keeps nothing") and `restart_policy.md` (the
  sieve takes up its chunks) disagree.
- Fix: on a restart without published references, clear p03's working folders, then align the two
  policies.

**D30. The automatic reference route never finishes below 500 accepted micrographs.** *Read.*
- Final ingestion of p03's sieve needs 500 picked micrographs (p03 stage:486-487), with no
  timeout, no upstream-finished trigger and no log. The master's `nmics` goes only to
  preprocessing (p00:351); the GUI's `increase_nmics` has no reader.
- A session with 320 accepted micrographs leaves p04 to p07 idle unless the user selects references
  by hand. The policy does not mention it.
- Fix: also set final ingestion when upstream has finished, gone quiet or reached `nmics`; log the
  wait; document the trigger.

**D31. The sieve's 10-minute idle trigger fires mid-session, and its final signal cannot be
taken back.** *Traced.*
- p05 sets final ingestion after 10 minutes without a set (p05 stage:66, 203): a grid exchange is
  enough. Up to 5,000 leftover particles then skip the coarse pass (SIEVE:785-821), a short fine
  chunk is flushed with `sieve_final=yes` (SIEVE:941-942), and the pool latches `l_sieve_final`
  for good (p06 stage:485-488, 590).
- `unset_final_ingestion` (p05 stage:306) stops further staging but cannot recall what was sent.
- Fix: trigger on upstream termination, or let a later non-final set clear the pool's flag.

**D32. The sieve's final signal can be missed.** *Traced.*
- The flush needs pending coarse particles (SIEVE:908); the main loop flags `sieve_final` only when
  every coarse chunk is complete after consuming (SIEVE:883). If nothing is pending at final
  ingestion, or the last coarse chunk finishes with nothing selected after the last fine chunk
  was made (SIEVE:1246-1249), no chunk ever carries the flag.
- Fix: when final ingestion is set and nothing is pending, publish a final marker.

**D33. p07's merge overwrites the multistate 3D labels.** *Read.*
- Every merge sets the `ptcl3D` state of every merged row to the publication's 0 or 1 (p07
  stage:478-481), including rows that already carry a 3D label from a run, against the ingestion
  policy (lines 61-65, 94-95) and `solve3D_addon_policy.md` §1.
- After each run, the next merge erases the labels: the GUI status shows every particle in state 1
  for most of the session (the vol3D messages are right), and `finalize` writes a stage project
  without its state assignments. The addon itself is unaffected; it takes labels from the frozen
  project.
- The tester asserts state 1 after a merge (tester:181).
- Fix: deselect when the publication deselects; otherwise keep an existing label and set 1 only on
  rows without one.

**D34. p07 starts an addon on row growth, and a refused addon stops the stage.** *Read and traced.*
- p07 starts an addon whenever rows were appended (p07 stage:584-590, 814-823), selected or not.
  The addon refuses an empty cohort, fewer than 5 particles per state, or a state without frozen
  particles (`simple_project_superset.f90`:94-106), and any job failure is fatal in p07 (p07
  stage:576-578).
- A publication of rejected stacks, or a small last hand-off, stops p07; a GUI restart then runs
  `solve3D` from scratch and loses all addon progress. Selection-only changes never trigger a
  run.
- Fix: trigger on the cohort as the addon defines it, with a minimum of at least `5*nstates`; treat
  a refusal as "wait for more".

**D35. Never-aligned particles in published stacks reach 3D with a random class.** *Traced.*
- `build_pool_publication` publishes a stack whole when any of its particles has `updatecnt > 0`
  (R2D:545-555).
- Above 500,000 selected particles the pool updates only `update_frac` of each partition; new
  particles first get a random populated class (POOL:314-331) and the unsampled ones come back
  with it and `updatecnt = 0` (POOL:655-658). p07 selects or rejects them on that class.
- This breaks the ingestion policy's "no unclassified particle reaches 3D". The publication test
  builds this case (tester:180) without checking those particles' state.
- Fix: publish selected particles with `updatecnt == 0` as state 0 (rows stay as p07 requires).

**D36. GUI threshold changes made in p01 are lost when it restarts.** *Traced.*
- The thresholds live only in memory (p01 stage:791-826) and are not in `restart_policy.md` §4.
  With D25, reprocessed movies are judged against the command-line defaults and accepted.
- Fix: persist them in p01's folder, or have the master replay the last values to a restarted
  stage.

**D37. p01's gain step blocks start-up with no stop check and no GUI status.** *Read.*
- `resolve_gain` runs before `init_gui` (p01 stage:174-176); its loops check no SIGTERM (p01
  stage:434-450, 476-490); `generate` stops p01 after 40 minutes without a new batch (p01
  stage:480).
- A long grid exchange early in a `flipgain=generate` session kills p01; a restart reads the movies
  again because an existing generated gain is not reused. A master stop waits its full timeout.
- Fix: poll the stop flag, open the pipe first, wait instead of stopping, reuse a generated gain.

**D38. p01 does not recognise a restart started with `dir_exec`.** *Traced.*
- A restart is recognised only by `outdir` existing (p01 stage:188-191), while `simple_stream`
  takes the execution folder from `dir_exec` (`simple_parameters_phases.f90`:128-131).
- Set numbering restarts at 1 and replaces `spprojs_completed/00001.simple` and following in place
  (JS:73-79, 158-164); downstream stages already have those names in their history and never
  import the replacements.
- Fix: recognise a restart by the execution folder, as `restart_policy.md` §2 states.

**D39. p07's disk use grows without bound.** *Read.*
- Every publication taken writes its class averages and JPEGs into `quality_selection/<id>/`
  (p07 stage:410-416), about 40 MB at box 256, never removed; nor is any
  `solve3D_addon/it_<n>/` (p07 stage:632-635).
- Fix: keep the last few quality folders; delete `it_<k>` once `it_<k+1>` is the frozen base.

**D40. Pool mask diameters are not bounded by the box.** *Suspected.*
- The pool passes `msk_crop` explicitly (POOL:526-536; R2D:213), bypassing the cap `parameters`
  applies when it derives it, and the box is 0 at p06's `params%new`.
- A GUI mask diameter above box × smpd gives workers a mask radius beyond half the box; the effect
  is unverified, and a failure would then stall the pool (D27).
- Fix: clamp to `(box_crop - COSMSKHALFWIDTH)/2` pixels and log it.

### Low

**D41. Smaller defects, by area.** *Read unless marked.*
- **Master and IPC**
  - Ctrl-C kills every stage at once: stages share the master's process group and reset SIGINT to
    the default (FORK:106), against p00:17-22 and `ipc_policy.md`:107.
  - Nothing is flushed before `fork`, so the master's start-up lines can appear once per stage in
    its log (p00:132, 151; FORK:107-122). *Suspected.*
  - The store copies any frame of the right size range into the type its tag names, with no
    length check (STORE:119-180, 276-338); after a resync (PIPE:218-224) a mis-sized frame could
    leave an allocatable descriptor undefined. *Traced.*
  - A start-up failure after the forks leaves the stages and the worker server running (p00:181,
    185); `failtime` and `stoptime` are rewritten on every poll (FORK:238; p00:225); no heartbeat
    during the 60 s optics wait (p00:309-319); `force_stop` kills only the stage process.
- **p01 and p02**
  - p06's raw-project path reads the never-updated root stub `optics_assignment.simple` and
    overwrites the optics groups just applied (p01 stage:838-840; STAR:738-741; R2D:270-273);
    p01's own call does nothing. *Traced.*
  - The movie rate spikes after a restart (WATCH:151-154, 210-212).
  - A failed job counts as one movie, not five (p01 stage:634, 661, 777-778); p02 reports accepted
    micrographs as imported (p02 stage:265-269); the last one to four movies of a session are
    never submitted (p01 stage:553).
  - p02 writes its project to the bare basename in `projinfo` (p02 stage:171, 239), against the
    skill's guardrail.
- **p03 and p04**
  - No accepted diameter bin leaves the box at 0 and p03 runs extraction with it (PK:222-231; p03
    stage:592-600). *Suspected.*
  - p04 raises a forced `box_extract` to 128 silently (p04 stage:244-246), while the master's UI
    promises to match an existing dataset (UI:229).
  - The template low-pass clamp is inverted: `min(max(30, x), 15)` is always 15 (PS:199, 333).
  - `make_pickrefs` deletes and rewrites `moldiam.txt` and the templates in place (PS:348-349,
    373-378), which p05 reads once with no retry (p05 stage:335). *Suspected.*
  - `params_sieve%nmics = 100` is never read (p03 stage:547); the 500-micrograph cap is checked
    once per pass; cycle 2's GUI status shows cycle 1's counts (p03 stage:977-982).
- **Sieve and p05**
  - Crash windows: coarse chunks are marked complete before the fine project is written
    (SIEVE:877-893); a coarse chunk is written before the file table (SIEVE:762-773).
  - New chunk ids come from the array size, not the highest folder (SIEVE:687, 739, 790, 863,
    928), so a skipped folder makes the next chunk overwrite an existing one.
  - The sieve's own queue ignores the chunk partition and turns off worker autoscaling when it
    reuses a server (SIEVE:242; QENV:271). *Traced.*
  - Counters reset on restart (the policy calls them cumulative); chunks with nothing selected are
    never counted; `kill` leaves tuning defaults and previews in place (SIEVE:259-292); rejection
    reasons with blanks are truncated on write (SIEVE:1343).
- **p06**
  - `stepwise=yes` imports one set per iteration once the pool is past its first threshold (p06
    stage:463-498).
  - No final project if p06 stops before its first complete iteration (p06 stage:303; R2D:266-278),
    against `pool2D_policy.md` §10.
  - A snapshot of the current iteration registers `frcs.bin` before checking it exists
    (R2D:421-438), which stops p06 when it is missing.
  - Publications carry no optics table and a raw copy of `ptcl2D` as `ptcl3D` (R2D:561, 586).
    *Suspected effect.*
  - Snapshot requests before the pool starts go unanswered (p06 stage:274, 660); a command-line
    `mskdiam` is replaced at iteration 10 (p06 stage:637-642); resolution limits stay on the first
    mask diameter (POOL:74, 711); stopping mid-iteration can pair iteration N-1's class averages
    with a live `frcs.bin` (R2D:236-250, 341-349). *Last one suspected.*
- **p07**
  - `nstates=1` or above 20 stops the stage (p07 cmd:81-82; `simple_commanders_solve3D.f90`:1030,
    1301-1306; `simple_gui_metadata_stream_solve3D_multistate.f90`:81).
  - An addon is reported as "running refine3D"; `particles_imported` counts deselected rows; old
    vol3D entries stay in the master's store after a restart (p07 stage:694-701).
  - `finish_run` adopts the job's `projinfo` wholesale and never carries `os_optics` (p07
    stage:664-665).
  - A malformed publication stops p07, and a restart takes the same one again; the mask diameter
    is read once, from the oldest export (p07 stage:327-348, 405-474).

## 3. Methodology

**M1, M2, M5.** Status in §1. Unchanged in substance.

**M6. The 3 October state choice can prefer a compact junk volume.** *Read.* The first key is
the number of connected components after a double Otsu threshold at 20 Å over the whole box (p03
stage:789, 1222-1240; `simple_image_bin.f90`:579-595). A strict second threshold tends to split an
elongated or multi-domain particle, and any speck counts, while a featureless "sink" state that
collects bad class averages forms one component and then wins whatever its view coverage or
population. The tester asserts that a population of 1 beats 9 (tester:583-584), which is the rule
as decided. Count components inside the mask and above a fraction of the largest, use
connectivity as a veto rather than the first key, add a population floor, and validate on known
datasets.

**M7. The fine tier uses the sieve preset and gates.** *Read.* `reject_cavgs` loads
`CAVG_QUALITY_MODEL_SIEVE_DEFAULT` for both tiers and the `sieve` gate context (SIEVE:1330, 1336),
although `model_cavgs_rejection.md`:62-63 defines the `chunk` context as exactly what a fine chunk
is (cleaned particles, at least 10,000). The sieve gates also drop the shared `res > 40 Å` gate
(FEATS:300-319), and flush chunks of a few hundred particles are scored the same way. Choose the
preset per tier, or retrain the sieve preset on fine-tier tables and record it.

**M8. The addon route has no cadence, rebase or verdict rule.** *Read.*
- Each addon accumulates every frozen particle at every stage box, so its cost grows with the
  stream while the cohort can be tiny; triggered on any growth, total cost grows roughly
  quadratically (p07 stage:814-823).
- Particles are aligned once, against the maps of their time, and never revisited; the first
  `solve3D`'s ladder is inherited by every later run; there is no rule for a fresh `solve3D` or a
  `refine3D` once the data have grown.
- `finish_run` never reads `solve3D_addon_report.txt` (p07 stage:654-670): a REGRESSED state
  silently becomes the frozen base of every later addon, although the addon policy leaves that call
  to the user (lines 300-301).
- Add a cadence rule (minimum cohort fraction), a rebase rule and the verdict to the ingestion
  policy; log and send the verdict.

## 4. Architecture and code hygiene

**G7 (updated). Dead code.** Callers checked with `rg` over `src` and `production`.
- SU: `update_user_params`, `wait_for_folder`, `wait_for_folder2`, `process_selected_refs`,
  `process_selected_refs_2`, `stream_datestr`. They keep `json_module`, the stream communicator
  and `simple_gui_utils` as module imports and are re-exported by `simple_stream_api`:25-31.
- MSU: `segdiampick_mics_multi`, `segdiampick_mics_multi_fixed_bins` ("kept for the current stream
  p03", which no longer calls them).
- WATCH: `write_checkpoint`, `clear_history` (tester only), `lastreporttime`, `ellapsedtime`,
  `chrono`, `fail_cnt`.
- POOL: `update_match_class_states` and its in-pool variant, `set_lpthres_type`,
  `get_pool_cavgs_jpeg_ntiles`, the `reference_generation` and `l_no_chunks` branches,
  `current_jpeg_scale`, `last_iteration_time`; in practice `generate_pool_jpeg` and the image block
  of `generate_pool_stats`, since refine2D writes the sprite sheet.
- C2D and chunk utils: the repick routines, `get_box`, `get_boxa`, `analyze2D_new_chunks`,
  `memoize_chunks`, `update_chunks`, `all_chunks_available`, `get_chunk_rejected_jpeg*`,
  `get_nchunks`, `converged_chunks`, `stream_chunk%remove_folder`.
- STAR: `stream_export_optics` with its own `h_clust` path, `stream_export_pick_diameters`,
  `stream_export_picking_references`.
- SIEVE: unused imports `commander_cluster_cavgs`, `COSMSKHALFWIDTH`, `I_BP40_100_CENTER_EDGE_VAR`;
  constants `SIEVE_BP40_100_CENTER_EDGE_VAR_MIN_LOG` (so a FEATS gate gates nothing) and
  `OVERFIT_CLUSTER_REJECT_FRAC`; reason code 101; `last_import`; `get_n_fine_rejected_ptcls`.
- FORK: the auto-restart machinery (only its tester uses it); the heartbeat always reports 0
  restarts.
- GUI update type: `set_sieverefs_selection` and three accessors; its header says p06 reads
  `sieverefs`, which it does not.
- Outputs with no reader in `src` or `nice`: `streamdata.simple`, `.poolstats`,
  `cls2D_thumbnail.jpeg`.

**G9 (updated). The UI and the reads disagree.**

| Program | Read but not offered | Offered but not read |
|---|---|---|
| `master` | `beamtilt`, `tilt_thres` (p00:361-362) | none found |
| `preproc` | `eer_upsampling`, `icefracthreshold`, `astigthreshold`, `nmics` | `dir_prev`, `tilt_thres`, `beamtilt`; `flipgain` omits `flip_auto` and `generate`; `beamtilt` help says `{yes}` |
| `assign_optics` | `tilt_thres`, `beamtilt`, `nmics` | none |
| `gen_pickrefs` | `pcontrast` | `nmics`, `optics_dir`, and `amsklp`, `ngrow`, `winsz`, `edge`, `pick_roi`, which the commander defaults but p03 does not read (the 3 October estimate uses the `AUTOMASK2D_*` constants) |
| `pick_extract` | `optics_dir` | `moldiam` (deleted by the commander), `nmoldiams`, `moldiam_max`, `pgrp`, the CTF, ice and astigmatism thresholds |
| `sieve_cavgs` | `refs` | `ncls`, `nptcls_per_cls` (required), `nchunksperset`, `nptcls_coarse`, `nptcls_fine`, `nmics`, `maxnchunks`, `dir_exec`; `mskdiam` is always replaced from `moldiam.txt` |
| `pool2D` | `stepwise` (not a parameter), `dynreslim`, `nsample_max`, `update_frac`, `lpstop`, `cenlp`, `center_type` | `projfile_optics` in practice |
| `solve3D_stream` | `nstates`, `nstages`, `lpstart`, `lpstop`, `nparts3D`, `nthr3D` | `nparts`; `dir_target` is described as the pick_extract folder |

**G12. Compile-time policy.** `doc/policies/compile_time_policy.md` asks for heavy components to
be allocatable and forbids snapshots of `parameters`.
- Inline `cmdline`: p01 stage:95, p04 stage:87-88, MSTAGE:45 and FORK:37 (copied into each other
  on every start), `chunk2D_state` in SIEVE:86.
- Inline `sp_project` in p01 to p05 (three of them in p03, stage:155-157), and inline `qsys_env`
  in p01, p03, p04, p05 and `ptcl_sieve` (SIEVE:127).
- Snapshots: p03 copies `params` for the sieve (p03 stage:543); `estimate_mskdiam` declares a local
  `parameters` because `automask2D` takes one (p03 stage:1180).
- Testers with plain stage locals: the master tester (line 203) and the p03 tester's
  `test_estimate_mskdiam` and `test_choose_state` (tester:580, 607), both added on 3 October.
- p00 imports the seven stage commanders at module level for one function (p00:53-59).

**G13. Fork without exec.** R1 item 1 is still open, and three defects come from it: D14 (the
inherited worker server), D41 (unflushed buffers; Ctrl-C to the whole group) and the mutex copied
mid-hold. Short of fork+exec, a child should reset what it inherits before `execute`: drop the
server object, flush and close inherited units, ignore SIGINT or take its own process group.

**G14. The update protocol with NICE is undefined.** `ipc_policy.md`:73 says a field is sent once,
but NICE re-sends fields and flags on every heartbeat, and the master forwards every answer (D13,
D22, D36). The contract should say whether answers are deltas or state, and who clears a restart
request; the master should forward only changes per stage.

**G15. Queued jobs have no common lifecycle.** p03 and p07 use `qsys_async_job` (exit code, no
timeout); the pool, the sieve chunks, p01's sets and `solve2D_chunks` poll marker files with no
failure path; no stage cancels its jobs on stop or adopts them on restart (D17, D21, D26, D27). One
contract for all of them is the main piece of work this review recommends (R13).

## 5. Tests and documentation

**Testers never run.** `simple_mic_import_tester`, `simple_mic_selection_tester`,
`simple_optics_groups_tester` and `simple_optics_maps_tester` are compiled with `BUILD_TESTS=ON`
but no suite calls them; `stream_refactor.md` names `unit_project` as their home.

**Tests that assert the defect.** p01's import counter (preprocess tester:164, D25); state 1 after
a merge in p07 (tester:181, D33); the pool's publication test without the state of never-updated
particles (tester:180, D35).

**Missing coverage.**
- The master loop (stop sequence, restarts during a stop, skipped stages, duplicate updates);
  `stream_pipe%send` against a reader that stops; a serialise round trip per metadata type; store
  frames of the wrong length.
- p01 restore with rejected movies; a job-set restore with an unfinished set above the highest
  completed one; `watch()` itself; optics grouping at scale; the gain loops' stop.
- p03: a cycle 2 selection after balancing; final ingestion against combining; a restart with
  existing chunks or `all.simple`; `vol_shape_descr` on two blobs or a blob and a speck; the picker
  with no accepted bin. p04: `prepare_pickrefs`, `moldiam.txt`, the box raise.
- Sieve: chunk generation from records, the fine merge, final staging and flush (staged-only
  above all), `submit`, the fine tier, partial-set and duplicate restart windows. p05: a restart
  with chunk folders and no new imports.
- p06: iterations, the ring and tidying, snapshots, `update_pool_dims`, `update_mskdiam`, the
  publication trigger and retention, a restart with a stale marker or a running job, a failed
  iteration, the final project.
- p07: `advance_jobs`, `start_solve3D`, `start_addon`, `finish_run` (including that the addon
  command line passes the addon's UI check), the newest-only import, a restart with existing
  run folders.
- STAR: the micrographs and optics tables, `outdir`, origin units.
- `solve2D_chunks`: nothing.
- No chained, master-less test of p01 to p07 on simulated data (R10).

**Documentation drift.**
- `ipc_policy.md`: line 73 (sent once), §3.3 (the resync), §7 (Ctrl-C stops in order).
- `restart_policy.md`: §1.1 (restarts during a stop), §1.2 (`dir_preprocess`), and the p01, p03,
  p05, p06 and p07 rows (D18, D19, D25, D26, D29, D36, D37, D38).
- `reference_generation_policy.md`: §3.1 item 6 against `restart_policy.md` (D29); the
  500-micrograph trigger (D30); the dead `nmics=100`; the claim that p03 and `make_pickrefs` use
  "the same measure", which holds in pixels only when cycle 1's class averages are downscaled.
- `stream_3D_ingestion_policy.md`: still calls the addon route out of scope (lines 15-16) though
  §4.5 drives it; claims matching by image index and untouched 3D parameters (D33); no cadence,
  failure or verdict rules (M8); the publication invariant D35 breaks.
- `solve3D_addon_policy.md`:307-311 describes a p07 import that no longer exists.
- `pool2D_policy.md`: §10 (the raw project), §6.5 (resolution limits), no word on `stepwise`.
- `ptcl_sieve_policy.md`: §2 omits `unset_final_ingestion`; §4 treats `REJECTION_FAILED` as live;
  §6.3, §8, §9 and §11 differ from the code (flush condition, both tiers use the model, coarse
  chunks never terminal in two-tier mode, counters not cumulative).
- `motion_gain_analysis_policy.md` does not cover gain generation or the stream's batch, wait and
  fallback constants.
- `stream_refactor.md`:166-171 (five-minute history, `export=.true.`); the stream skill lists
  "snapshot state" under the pool modules.
- Module headers: p00:158-160 (forks before any thread); the GUI update type (`sieverefs`);
  `simple_optics_maps.f90`:18; MSU ("kept for p03"); p01 has no RESTART section and omits
  `split_mode`; p03's METHOD section cites the 30 September review; p07 cmd:14-16 (the master sets
  no `nparts`).
- UI help: `lpstart {15}` and `nsample_fine {2000}` for `sieve_cavgs`; `beamtilt {yes}`; the
  `solve3D_stream` `dir_target` description.

## 6. Recommendations

### R12. Fixes before the next GUI session

- D13: deduplicate updates per stage, send nothing to a stopping stage, time-limit partial frames.
- D14: drop the inherited server object in the child and connect through `worker_server`.
- D15: balance into a separate project; read selections against the cycle's own stack.
- D17: exit-code files for chunk jobs; rejection failures mark the chunk failed and complete.
- D18, D19: record chunked micrographs; make the sieve on restart.
- D20: guard the FRC write.
- D21: stop or adopt p07's job in `finalize` and on restart; one folder per run.
- D22: no restart while stopping; act on a restart request once.
- D33: keep 3D labels in the merge.
- D35: publish never-updated particles as state 0.

### R13. One lifecycle for queued jobs (G15)

Every queued job gets an exit-code file, a liveness timeout, a cancel on stop and an
adopt-or-cancel on restart, and runs in a folder no restarted job reuses. `qsys_async_job` is the
natural base; it needs the timeout and a cancel. Apply it to p01's sets, the sieve chunks, the
pool iteration, p07's jobs and `solve2D_chunks`. This closes D17, D21, D26 and D27 in one design.

### R14. The update protocol (G14)

Decide whether NICE answers are state or deltas; have the master forward only changes to each
stage, acknowledge restart requests and have NICE clear them; replay the last values to a
restarted stage (D36). Write it into `ipc_policy.md`.

### R15. Fork+exec, or a clean child (G13)

Start stages by fork and exec of `simple_stream prg=<stage>` with the pipe descriptors on the
command line. Until then, reset inherited state in the child before `execute`.

### R16. Restart conformity (G11)

Make each stage keep `restart_policy.md`, then correct the policy where it was wrong:
- p01 restores from every micrograph of the completed sets, sets the counter to the highest
  `importind`, persists the thresholds, recognises `dir_exec` (D25, D36, D38);
- p03 clears its working folders on a restart without published references (D29);
- p05 records micrographs, not sets, and resumes the sieve (D18, D19);
- p06 stops the previous pool job and removes its marker (D26);
- the master accepts its own link on restart and never starts a skipped stage (D23, D24).

### R17. Final-ingestion triggers

p03 (500 micrographs, D30) and p05 (10 idle minutes, D31) both guess when upstream is done. Give
the stages an upstream-finished signal (the upstream stage's `TERM_STREAM`, or a quiet period
measured from upstream's own state), make the pool's final flag retractable until upstream has
terminated, and make sure the final signal is always sent (D28, D32).

### R18. Scaling

Incremental optics grouping in p02 (D16); retention of p07's quality and addon folders (D39); an
addon cadence (M8).

### R19. Methodology

- The state choice as a veto, inside the mask and above a component size, with a population floor
  (M6).
- A fine-tier preset or a retrained sieve preset with recorded provenance (M7, M5, R5).
- An addon cadence, rebase and verdict rule in the ingestion policy (M8).
- Weights instead of duplication (M1, R2); the 3D reference settings as inputs (M2, R4).
- One validation dataset per route, run before methodology changes land (R11).

### R20. Hygiene

- The compile-time policy in the stream code (G12): allocatable `cmdline`, `sp_project` and
  `qsys_env` components; no `parameters` snapshots (give `automask2D` explicit settings); stage
  locals in testers as polymorphic allocatables.
- Remove the dead code in G7 and the re-exports in `simple_stream_api`.
- Bring every UI entry in line with what its stage reads (G9); have p03's mask estimate read the
  `gen_pickrefs` automask inputs it advertises.
- Register the four shared testers; replace the tests that assert defects; add the coverage in §5.
- Correct the documentation drift in §5.

### Suggested order

1. R12, the fixes, each with a test that fails before the fix.
2. R13 and R14: the job lifecycle and the update protocol, policy first.
3. R16 and R17: restart conformity and the final-ingestion triggers.
4. R15: fork+exec.
5. R18 and R19, behind parameters whose defaults are today's behaviour until validated.
6. R20 alongside, file by file as each is touched.
