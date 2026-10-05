# Stream fix plan, 5 October 2026

Follow-up to `doc/refactoring_notes/completed/stream_area_review_2026-10-05.md` (cited as the
review; its defects D13 to D41). Items D3, M2 to M5, G2, G4, G7, G10, G11 and R4 to R10 are
those of `completed/stream_area_review_2026-10-02.md`, with their status as the review gives it.

Scope: D3, M2, M3, M4, M5, G2, G4, G7, G10, G11, R4, R8, R9 and R10 of the 2 October review, and
the review's defects D13 to D15 and D17 to D40, with the fixes it suggests, adapted to the
decisions below. Where a defect is covered by a workstream, the workstream names it; the table
at the end maps every defect to its place. Methodology changes ship
behind a parameter whose default is today's behaviour until a validation run has accepted them
(R11), unless a decision below says otherwise.

## Status, 5 October 2026

Implemented in the order below, except WS4, which was dropped from this round. Nothing has been
compiled or run; the build and the test suites are still to be run. Notes per workstream:

- **WS0, WS9, WS10, WS1, WS8, WS12:** as planned. Beyond the plan: p06 now publishes for 3D from
  iteration `EXPORT_START_ITER` (25) on, not after it, because the final run stops at
  `FINAL_ITER` (25) and a session whose final set came early never reached 3D; the sieve hands on
  an empty final set (`hand_off_final_set`) when no chunk carries the flag; a master that skips
  preprocessing marks the existing folder finished.
- **WS5:** `rank_cavgs_stk` in `simple_imgarr_utils`; p03's reprojection and p04's `make_pickrefs`
  are local queued jobs. The sieve's `commander_cluster_cavgs` import was already gone.
- **WS6:** done, with these limits: the pool's native box, pixel size, mask diameter and hard
  low-pass limit are pool state (`simple_stream2D_state`), not stage components, since the pool's
  utilities are module code; `setup_downscaling` keeps its `parameters` writes for
  `solve2D_chunks` (decision 10), as do the chunk command lines. The UI follows the G9 table,
  except pool2D's `projfile_optics` (never offered) and a `center_type` the pool now forwards.
- **WS7:** `simple_stream_master_resources`; the variable sets are in
  `doc/policies/stream/README.md`. p03's extractions keep a named constant of 4 threads.
- **WS2:** the 3D route's inputs, the fixed result path (`SOLVE3D_CAVGS_FINAL_DIR`) and the
  state veto (10% and 10%, to be set by a validation run). `nstates_pickrefs` must be at least 2,
  because `solve3D_cavgs`' conditional restarts need two states.
- **WS3:** as planned; the dead growth guard is removed.
- **WS4:** not done, except its decision 8 bullet: the fine tier keeps the sieve preset and
  gates, confirmed and documented in `ptcl_sieve_policy.md` §8 and `model_cavgs_rejection.md`.
- **WS11:** the four testers are registered in `unit_project`; the chained tests are
  `lib_stream` "sieve to 3D" and "movies to 3D" (`simple_stream_chain_tester`), never run.

## Decisions

| # | Question | Decision |
|---|---|---|
| 1 | When p07 starts an addon run | When the cohort (selected, not frozen) reaches max(5·nstates, 10% of the frozen particles) |
| 2 | A REGRESSED addon verdict | Logged and sent to the GUI; the result stays the frozen base |
| 3 | Rebasing the 3D solution | None for now; recorded as a known gap |
| 4 | p03's choice of state | Connectivity as a veto (components inside the mask, above a size fraction of the largest) with a population floor, then distinct projection directions |
| 5 | Pool sampling past the cap | The sampled stacks are the update set; `update_frac` is not applied again inside the sample |
| 6 | Pool dimension growth | Kept; its guard is fixed so it can fire |
| 7 | The sieve's scoring mask | The given `mskdiam` when it is positive, otherwise 0 (the box-sized disc) |
| 8 | Fine-tier model | The sieve preset for both tiers, documented as a decision |
| 9 | Commanders run by stages (p03 reproject, p04 `make_pickrefs`) | Queued jobs, always on the local machine |
| 10 | `solve2D_chunks` | Kept as it is for now, with the chunk layer it uses |
| 11 | `stepwise` | Registered as a parameter |
| 12 | Sieve configuration | A small `ptcl_sieve` settings type |
| 13 | UI entries | Inputs no stage reads are removed; inputs read but hidden are exposed at developer visibility |
| 14 | Resources | Per-stage named defaults with per-stage `SIMPLE_STREAM_*` overrides, applied in the master and passed down; no literals in stages; budget logged at start |
| 15 | Jobs at a stop or restart | Cancelled on stop; a restart never reuses a folder that held a job without an exit code |
| 16 | p03 restart before publication | Clear its working folders and run the plan from cycle 1 |
| 17 | p01's GUI thresholds across a restart | p01 persists them in its folder |
| 18 | Wire format | Byte copies, with guards |
| 19 | Chained tests | Both: p05 to p07 from simulated particles in the library suite, and p01 to p07 from simulated movies nightly |
| 20 | p03's 3D settings from the GUI | Passed through the master at developer visibility; names p07 also reads, or that are ambiguous, take a `_pickrefs` suffix end to end (`nstates_pickrefs`, `nstages_pickrefs`, `lpstop_pickrefs`, `nspace_pickrefs`, `nthr3D_pickrefs`); `solve3D_cavgs`-only names (`nrestarts_collapse`, `lpstart_ini3D`, `lpstop_ini3D`) keep theirs |
| 21 | D13: where the update flood stops | In the master: it forwards to each stage only what differs from the last update it sent, and nothing to a stopping stage; NICE unchanged |
| 22 | D14: the inherited worker server | Clean the child now (drop the server object, close its sockets, clear warm-up registrants); fork+exec later |
| 23 | D15: cycle 2 selections | Balance into a project of its own (`all_balanced.simple`) |
| 24 | D17: a failed sieve chunk | Retry once in a fresh folder; on a second failure mark it failed and complete and hand none of its particles on, with a warning and a count in the GUI |
| 25 | D22: restart flags NICE keeps set | Fixed in the master only: a restart request acts once, again only after the key disappears, and never while stopping |
| 26 | D27: jobs killed by the scheduler | Exit codes only; a job killed without writing one is not detected, and that is documented |
| 27 | D30, D31: what ends final intake | An upstream marker, with a quiet fallback triggered when no new movie has been detected for 15 minutes; p03 keeps 500 micrographs as an early trigger |
| 28 | D28: p03's cycle 2 | Cycle the sieve once more right after setting final ingestion |
| 29 | D31: data after a final signal | Retracted: markers removed when movies arrive again, p05 unsets final ingestion, the pool clears its final flag on a non-final set |
| 30 | D37: p01's gain generation | Waits without a limit, polling SIGTERM and reporting to the GUI; an existing generated gain is reused on restart |
| 31 | D39: p07's folders | The last 3 quality folders and those of publications that started a run; an addon iteration folder is deleted once a later run has completed and is the frozen base |
| 32 | D40: a mask diameter beyond the pool's box | Clamped to `(box_crop - COSMSKHALFWIDTH)/2` pixels and logged |

## Workstreams

### WS0. Fixes before the next GUI session (review R12). Medium

Each fix lands with a test that fails before it. D21, D33 and D35, also on the review's list, are
in WS8 and WS1 and should land with this workstream.

- **D13, a partly written update frame blocks the master (decision 21).**
  - The master keeps the last update it sent each stage and sends a new one only when it differs
    (p00:273-278; MSTAGE:151-157). NICE is unchanged.
  - No update goes to a stage that has been asked to stop.
  - `stream_pipe%send` gives up on a partly written frame after a time limit (a named constant of
    about 2 s) and marks the channel broken (PIPE:102-124); the master discards a broken channel
    before the stage's next start.
  - `ipc_policy.md`: NICE's answers are state, re-sent every heartbeat; the master forwards changes
    per stage (G14).
  - Tests: the pipe tester with a reader that stops mid-frame; the master stage sending an
    unchanged update once.
- **D14, forked stages inherit the master's persistent-worker server (decision 22).**
  - In the child, before the stage's commander runs (`stream_master_stage_fork%execute` or
    `forked_process`), drop the inherited `persistent_worker%server` without killing it, close its
    inherited socket descriptors and clear the warm-up registrants.
  - Stages connect to the master's server as clients through the `worker_server` address p00
    already passes (p00:404, 420, 432), and p00 passes it to p03 too.
  - Fork+exec (review R15) stays the lasting fix and is not in this plan.
  - Test: a forked child sees no server object and builds its queue as a client.
- **D15, a cycle 2 GUI selection publishes the wrong images (decision 23).**
  - `balance_classes` writes the balanced project to its own file (`all_balanced.simple` in
    `balance_classes/all`), and `solve3D_cavgs` runs on that file (p03 stage:679-689, 1326-1330).
  - The cycle 2 project keeps its own class averages, so `save_pickrefs_selection` reads a
    selection against the stack the GUI was shown (p03 stage:895-916).
  - `finish_solve3D` reads the balanced project where it reads the 3D result.
  - Test: a cycle 2 selection after balancing publishes exactly the selected classes.
- **D17, sieve chunk jobs have no failure path (decisions 24 and 26).**
  - Chunk jobs are submitted with an exit-code file (SIEVE:1139, 1157) and polled as
    `qsys_async_job` polls its jobs; a non-zero code is a failure.
  - Before rejection runs, the sieve checks what `reject_cavgs` needs (class averages present,
    counts matching); a chunk that lacks it is a failure, instead of rejection stopping the stage.
  - On a first failure the chunk is resubmitted once, in a fresh folder. On a second it is written
    as `REJECTION_FAILED` and `COMPLETE`, frees its slot, and counts as complete for the final
    flush (SIEVE:441-451, 909); none of its particles is handed on. The sieve logs a warning and
    reports the number of failed chunks and their particles in its GUI status.
  - With decision 26, a job the scheduler kills before the script writes its exit code (walltime,
    lost node) is still not detected and keeps its slot; `ptcl_sieve_policy.md` says so.
  - `ptcl_sieve_policy.md` §4 and §6.3.
  - Tests: a chunk with a non-zero exit code, retried and then completing; a chunk failing twice; a
    chunk with no class averages; a final flush with a failed coarse chunk.
- **D20, a staged-only final flush writes an uninitialised FRC object.**
  - Guard `frcs%write` and `add_frcs2os_out` in `merge_chunk_projfiles_without_sigma2` with
    `frcs_initialised`, as `merge_chunk_projfiles` does (PU:507 against PU:214-217).
  - Test: a final flush of only the staged chunk.
- **D22, restart flags that NICE keeps set (decision 25).**
  - The master applies no restart once it is stopping (p00:265-272, 295-300).
  - It acts on a restart request once, and again only after the key has disappeared from an
    answer. NICE is unchanged.
  - `ipc_policy.md` and `restart_policy.md` §1.1.
  - Test: the master's restart handling with a key present on consecutive answers and during a
    stop.

### WS9. Dead code (G7). Small

- Remove what the review's G7 lists, except the chunk modules and routines `solve2D_chunks`
  uses (decision 10): the stream utilities, `simple_mini_stream_utils`' unused pickers, the
  watcher's leftovers, the pool's unused routines and branches, the repick routines, the STAR
  export's unused writers, the sieve's unused imports, constants and reason code, the forked
  process's auto-restart machinery, the unused GUI update accessors, and the outputs nobody reads
  (`streamdata.simple`, `.poolstats`, `cls2D_thumbnail.jpeg`).
- Remove the matching re-exports in `simple_stream_api`, and the module imports that only the
  dead code needed.
- Done when: `rg` finds no caller of anything removed, the forked-process tester no longer tests
  auto-restart, and the build has no new warnings.

### WS10. Wire format (G10). Small

- The master's store checks each frame's length against its tag's type and drops a mismatch with
  a warning (STORE:119-180, 276-338).
- Metadata `kill` resets every field, not only the flags (`simple_gui_metadata_base.f90`:43-48 and
  the extensions).
- `ipc_policy.md` states the rule: a type sent over a pipe has no allocatable or pointer component
  unless it overrides `serialise` (as `gui_metadata_vol3D` does); `gui_metadata_project` does not
  cross the pipes and is outside the rule.
- Tests: a serialise round trip for every type sent by a stage; the store dropping a mis-sized
  frame.

### WS1. The 3D hand-off (D3, M3). Medium

- **p06:** in `build_pool_publication`, a selected particle with `updatecnt == 0` inside a
  published stack is published as state 0; rows stay as p07 requires (review D35).
- **p07 merge:** keep an existing 3D label; deselect when the publication deselects; set state 1
  only on rows without a label (D33). Match rows by `indstk`, with a check that it agrees with the
  offset.
- **p07 addon:** trigger on the cohort per decision 1; a refusal (empty or small cohort, a state
  without frozen particles) waits for more particles instead of stopping the stage (D34).
- **p07 verdict:** read `solve3D_addon_report.txt` after each addon run, log the per-state
  verdict and send it with the status (decision 2).
- **p07 retention (D39, decision 31):** keep the last 3 `quality_selection/<id>` folders and those
  of publications that started a run (p07 stage:410-416); delete an addon iteration folder, and a
  superseded run folder under WS8's fresh-folder rule, once a later run has completed and is the
  frozen base (p07 stage:632-635). `restart_policy.md` and the ingestion policy say what is kept.
- **Policy:** `stream_3D_ingestion_policy.md` takes the addon route into scope with the cadence,
  refusal and verdict rules, states that particles are aligned once and never revisited (decision
  3), and corrects the claims on matching and on 3D parameters. Remove the stale p07 sentence in
  `solve3D_addon_policy.md`:307-311.
- Tests: the pool tester checks the state of never-updated particles in a published stack; the
  p07 tester keeps a label of 2 through a merge, starts no addon below the cohort threshold,
  survives a refused addon, and prunes quality and iteration folders as decided.

### WS8. Restart semantics (G11). Large

The restart and job defects of the review: D18, D19, D21, D23, D24, D25, D26, D27, D29, D36, D37
and D38. D22 is in WS0.

- **Job cancel (decision 15).** Give `qsys_async_job` and the marker-file jobs a cancel: a local
  job by its process group (recorded at submit), a scheduler job by its id (captured at submit),
  a persistent-worker task by a cancel message if the server supports one. Each stage cancels its
  jobs in `finalize`. A restart starts any job in a new run folder whenever the old folder holds a
  job without an exit code. This covers p01's sets, the sieve chunks, the pool iteration and p07's
  jobs; `solve2D_chunks` is left as it is (decision 10).
- **Exit codes for every queued job (D27, decision 26).** The pool iteration is submitted with an
  exit-code file (or through `qsys_async_job`) and p06 polls it: a non-zero code reports the
  iteration failed, frees the pool and retries the iteration once before stopping the stage with
  a clear message, instead of waiting for `REFINE2D_FINISHED` for good (POOL:368, 372, 827-835).
  p01's sets get the same exit-code check. No liveness timeout: a job the scheduler kills before
  its script writes an exit code is not detected, and `restart_policy.md` and the stage policies
  say so.
- **p01:** restore from every micrograph of the completed sets, whatever its state; the import
  counter becomes the highest `importind`; recognise a restart by the execution folder, from
  `outdir` or `dir_exec` (D25, D38). Persist the GUI thresholds in p01's folder, written by
  temporary and rename on each update and read back on restart (decision 17, D36).
- **p01's gain step (D37, decision 30):** open the pipe before `resolve_gain` and report "waiting
  for movies for the gain" to the GUI (p01 stage:174-176); both gain loops poll the SIGTERM flag
  (p01 stage:434-450, 476-490); `generate` waits without a limit instead of stopping after 40
  minutes (p01 stage:480); an existing `gainref_generated.mrc` is reused on restart. The commander
  header and `motion_gain_analysis_policy.md` say so.
- **p03:** on a restart without published references, clear the working folders (stream, extract,
  picker, sieve chunks, `spprojs_sieved`, solve2D, solve3D, quality selection, balancing) and run
  the plan from cycle 1 (decision 16, D29). Published references keep today's behaviour.
- **p05:** record chunked micrographs (`projname`, `micind`) rather than sets, re-import partial
  sets with their chunked micrographs marked, and make the sieve on restart even with nothing new
  to import (D18, D19).
- **p06:** cancel the previous pool job and remove `REFINE2D_FINISHED` before the restart cleanup
  (D26).
- **p07:** cancel the running job in `finalize`; each `solve3D` and addon run gets a folder of its
  own (D21).
- **Master:** accept an existing `dir_preprocess` link to the same target and terminate the link
  name for C (D24); never start a skipped stage, keeping `skipped` sticky and ignoring restart
  requests for it (D23).
- **Policy:** `restart_policy.md` per stage, including the job cancel and the fresh-folder rule;
  `reference_generation_policy.md` §3.1 item 6 aligned with it.
- Tests: per stage, a restart with leftover folders and an unfinished job; p01 restore with
  rejected movies; p05 restart with chunk folders and no new imports; the master's restart of a
  skipped stage; a pool iteration with a non-zero exit code; the gain loops stopping on SIGTERM.

### WS12. Final-ingestion triggers (D28, D30, D31, D32; review R17). Medium

- **Markers (decision 27).** Two markers, written by temporary and rename in a stage's own folder:
  - *finished*: a stage writes it in `finalize`, after handing on its last project;
  - *idle*: p01 writes it once no new movie has been detected for 15 minutes (a named constant,
    `MOVIES_IDLE_TIME_S = 900`) and its own sets are complete. Every later stage writes it in turn
    once its upstream folder holds the idle or finished marker and it has handed on every project
    it took.
- **Retraction (decision 29).** p01 removes its idle marker as soon as a movie arrives again; each
  stage removes its own when new upstream work arrives. p05 then unsets final ingestion (p05
  stage:306), and the pool clears `l_sieve_final` when a non-final set arrives (p06
  stage:485-488, 590).
- **p03 (D30).** Its sieve's final ingestion is set at 500 picked micrographs, as now, or when
  p01's folder holds the idle or finished marker (p03 stage:486-487). The wait is logged.
- **p03 (D28, decision 28).** Right after setting final ingestion, p03 cycles the sieve once more,
  so the leftover chunk exists before `run_cycle2` checks `get_finished` (p03 stage:484-491,
  650-656).
- **p05 (D31).** Final ingestion is set when p04's folder holds the idle or finished marker,
  replacing the 10 idle minutes (p05 stage:66, 203).
- **The final signal is always sent (D32).** When final ingestion is set and nothing is pending,
  the sieve still hands on a final marker set (no particles, `sieve_final=yes`) that the pool reads
  as the end of intake (SIEVE:883, 908, 1246-1249).
- **Restart:** a restarted stage evaluates the markers afresh; a stale idle marker from before the
  restart is removed if its upstream has new work.
- **Policies:** `ptcl_sieve_policy.md` §2 and §6 (triggers, retraction, the empty final set),
  `reference_generation_policy.md` §3 (the triggers of cycle 2), `pool2D_policy.md` (the final flag
  and its retraction), `restart_policy.md` (the markers).
- Tests: the idle marker written after the quiet period and removed on a new movie (a short test
  time); p05 setting and unsetting final ingestion with the marker; the pool clearing its final
  flag; the empty final set; p03's extra cycle before combining.

### WS5. Layering (G2). Small to medium

- Drop the sieve's unused `commander_cluster_cavgs` import.
- Move the core of `commander_rank_cavgs` into a library routine that the pool utilities call
  (R2D:139-150).
- Run p03's reprojection and p04's `make_pickrefs` as queued jobs on the local machine (decision
  9): a `qsys_async_job` with a local queue environment, the exit code polled, the outputs read
  from the job folder. `make_pickrefs`' `moldiam.txt` and templates are then renamed into place
  (review D41).
- `solve2D_chunks` stays where it is (decision 10).

### WS6. Parameters and UI (G4, R8). Medium

- **No assignments to `params` after `new`:** values known only at run time become stage
  components: p01's gain fields, thresholds and `split_mode`; p04's `split_mode`, `smpd` and
  `box`; p05's `mskdiam`; p06's `smpd`, `box`, `mskdiam` and `lpstop` (no re-read from the command
  line), set in one helper when the pool's first import makes them known.
- **Pool and chunk command lines** no longer carry the native box and pixel size.
- **Sieve settings type (decision 12):** `ptcl_sieve%new` takes a small settings type with the
  fields it reads (mask diameter, starting low-pass, resources, tier overrides, final-ingestion
  rules); p03 and p05 fill it, and p03 stops copying `params`.
- **Sieve scoring mask (decision 7):** score with the configured mask diameter when positive,
  otherwise 0. In p05 that is the diameter from `moldiam.txt`, in p03 the box default. Record it in
  the sieve preset's spec and `ptcl_sieve_policy.md`; p05's selections change, so this needs a
  validation run before release.
- **`stepwise` (decision 11):** register it in `simple_parameters.f90` and its parser, move its
  default from `production/simple_stream.f90` into the p06 commander's normalisation, count only
  the particles of the current import against the threshold (review D41), add it to the `pool2D`
  UI at developer visibility, and document it in `pool2D_policy.md`.
- **Automask settings:** `automask2D` takes explicit settings instead of a `parameters` object;
  p03's mask estimate uses `gen_pickrefs`' `amsklp`, `ngrow`, `winsz` and `edge`.
- **UI (decision 13):** per the review's G9 table, remove inputs no stage reads and add the inputs
  read but hidden at developer visibility; correct the help texts it lists.
- **STAR export** honours `outdir`.
- Tests: the commanders' normalisation and range checks; the sieve settings type; `stepwise`'s
  threshold.

### WS7. Resources (R9). Small to medium

- One routine in the master derives, for each stage, `nthr` and `nparts` and the job-level
  resources (2D job threads and parts, 3D job threads and parts, chunk threads, chunk count,
  reprojection threads), from named per-stage defaults overridden by per-stage environment
  variables (decision 14), passes them on the stage command lines, and logs the table at start.
- Stages use no literals: p03's 2D and 3D job resources and reprojection threads, the sieve's
  threads and chunk count, the pool's.
- Environment variables, one set per stage: `SIMPLE_STREAM_<STAGE>_{NTHR,NPARTS,PARTITION}`.
  p03 reads the `REFGEN` set rather than the preprocessing partition; p07 gets a set of its own;
  the sieve reads the `CHUNK` set; the pool's legacy reads of `REFGEN` go with WS9.
- Docs: the variables and the defaults in `doc/policies/stream/README.md` (or a resources section
  of it).

### WS2. Reference generation (M2, R4). Medium

- **3D settings as inputs (decision 20):** `nstates_pickrefs`, `nstages_pickrefs`,
  `nrestarts_collapse`, `lpstart_ini3D`, `lpstop_ini3D`, `lpstop_pickrefs`, `nspace_pickrefs` and
  `nthr3D_pickrefs`, registered in `simple_parameters.f90` and its parser, with defaults in p03's
  commander equal to today's values (3, 3, 3, 100, 20, 8, 50, 16) and range checks
  (`nstates_pickrefs ≥ 1`, `nstages_pickrefs ≥ 1`, `lpstart_ini3D > lpstop_ini3D ≥
  lpstop_pickrefs > 0`, `nspace_pickrefs ≥ 1`). p03 maps them onto the `solve3D_cavgs` and
  `reproject` command lines; `NSTATES3D` becomes `nstates_pickrefs`, and the state arrays become
  allocatable. The master and `gen_pickrefs` expose them at developer visibility; the master
  forwards each to p03 only when set and deletes them from preprocessing's command line.
  `pgrp=c1` and `prune=no` stay fixed.
- **Fixed result path:** the `solve3D_cavgs` restart driver records its final folder in a fixed
  file; p03 reads it instead of scanning for the highest number.
- **State choice (decision 4):** count connected components inside the mask and above a fraction
  of the largest; a state passes the veto when it is one object and holds at least a population
  floor; among those, the most distinct projection directions, then population. When no state
  passes, fall back to the 3 October order. The size fraction and floor are named constants,
  proposed at 10% each, to be set by a validation run.
- **Policy:** `reference_generation_policy.md` §3 (settings, state choice) and §7.
- Tests: `vol_shape_descr` on one blob, two blobs and a blob with a speck; `choose_state` with the
  veto and the fallback; the commander's defaults and checks.

### WS3. Pool schedule (M4). Medium

- **Sampling (decision 5):** when the pool samples stacks past the cap, the sample is the update
  set: refine2D runs with `update_frac=1` on it; out-of-sample particles keep their parameters.
  The user's `update_frac`, when given, still applies.
- Drop `prev_pop_even` and `prev_pop_odd`; reset `center` when the pool returns to a full update.
- Draw new particles' classes from a generator seeded per iteration, outside the OpenMP loop.
- Recompute `lpstart`, `lpcen` and the ramp when the mask diameter changes at iteration 10.
- **Dimension growth (decision 6):** fix the guard (`ncls_glob < ncls_max` can never be true) so
  growth fires as intended, and document what growth does to class averages and FRCs.
- **Mask bound (D40, decision 32):** `update_mskdiam` and `setup_downscaling` clamp the mask
  radius to `(box_crop - COSMSKHALFWIDTH)/2` pixels and log the clamp (POOL:526-536; R2D:213), so a
  GUI or command-line diameter beyond the box never reaches the workers.
- **Policy:** `pool2D_policy.md` §6 and §12.
- Tests: the update set under sampling; reproducible class draws for a fixed seed; the guard; a
  mask diameter beyond the box clamped on the pool's command line.

### WS4. Rejection (M5, R5). Large, needs datasets

- Record each preset's provenance in its spec: datasets, commit, mask diameter, box, pixel size,
  context.
- The learner selects models by leave-one-dataset-out; feature policies that pair a feature with
  its exact negation are rejected.
- The fine tier keeps the sieve preset (decision 8); `ptcl_sieve_policy.md` and
  `model_cavgs_rejection.md` say so and why.
- The compatibility filter's fit-and-filter on one batch is documented as a known gap.
- A regression test pins each preset's decisions on fixed feature tables.

### WS11. Tests and documentation (R10). Medium, alongside

- Register `simple_mic_import_tester`, `simple_mic_selection_tester`,
  `simple_optics_groups_tester` and `simple_optics_maps_tester` in `unit_project`.
- Replace the tests that assert defects (p01 tester:164, p07 tester:181, pool tester:180) as each
  workstream lands.
- Chained tests (decision 19): p05 to p07 from simulated particle stacks in the library suite;
  p01 to p07 from simulated movies in the nightly suite, both master-less, driving the stage types
  in-process.
- Each workstream updates the policy it touches; the review's documentation-drift list is closed
  as part of the workstream that owns each item.

## Order

1. WS0, with D21, D33 and D35 from WS8 and WS1: the fixes before the next GUI session.
2. WS9 and WS10.
3. WS1, WS8 and WS12, with WS8's job cancel first.
4. WS5, WS6 and WS7.
5. WS2, WS3 and WS4, each with a validation run before its default changes.
6. WS11 throughout; the chained tests once WS1, WS8 and WS12 have landed.

## Where each defect of the review is fixed

| Defect | Workstream | Defect | Workstream |
|---|---|---|---|
| D13 update frame blocks the master | WS0 | D27 killed jobs never seen to end | WS8 |
| D14 inherited worker server | WS0 | D28 cycle 2 before the last particles | WS12 |
| D15 cycle 2 selection | WS0 | D29 p03 restart mixes sieve state | WS8 |
| D17 no failure path for chunks | WS0 | D30 no end below 500 micrographs | WS12 |
| D18 p05 restart drops micrographs | WS8 | D31 idle trigger mid-session | WS12 |
| D19 p05 restart never resumes | WS8 | D32 final signal missed | WS12 |
| D20 staged-only flush FRCs | WS0 | D33 merge overwrites 3D labels | WS1 |
| D21 second 3D job after restart | WS8 | D34 addon trigger and refusal | WS1 |
| D22 restart flags kept set | WS0 | D35 never-aligned particles published | WS1 |
| D23 restart of a skipped stage | WS8 | D36 p01 thresholds lost | WS8 |
| D24 `dir_preprocess` link | WS8 | D37 p01 gain step blocks | WS8 |
| D25 p01 restore | WS8 | D38 p01 `dir_exec` restart | WS8 |
| D26 jobs left running | WS8 | D39 p07 disk growth | WS1 |
| | | D40 pool mask beyond the box | WS3 |

## Not in this plan

From the review: D16 (p02's regrouping from scratch; now workstream H of `stream_followup_plan_2026-10-05.md`), most of D41 (parts land in WS5, WS6 and
WS10), M6 beyond the state choice (now workstream G of `stream_followup_plan_2026-10-05.md`), M7 (decision 8), M8 beyond cadence and verdict (now workstream I of `stream_followup_plan_2026-10-05.md`), G12 and G13
(the compile-time policy, fork+exec), and the parts of G14 and G15 that WS0 and WS8 do not cover
(a written update protocol with NICE; liveness timeouts, decision 26). G12, G14 and G15 are now
workstreams J, K and L of `stream_followup_plan_2026-10-05.md`, which also replaces decision 26;
G13 stays deferred there.
