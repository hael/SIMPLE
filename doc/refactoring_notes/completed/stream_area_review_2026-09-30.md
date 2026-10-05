# Stream area review: state of the streaming pipeline, September 2026

Line numbers refer to `master` at commit 4a2f58483.

Out of scope: the `abinitio3D_addon` route in p07 (`start_abinitio3D_addon`,
`finish_abinitio3D_addon` and the add-on loop) is under construction and is not reviewed here.

## Verdict

The streaming layer has become a second application with its own conventions rather than a
client of the library's. Stages are single 300-530-line `execute` routines with host-associated
internal procedures; the IPC protocol code is copy-pasted seven times; resources are literals;
`parameters` is mutated after construction in about 45 places; and workflow-level methodology
(how picking references are generated, how the 3D stage ingests the 2D pool) is decided inside
stage code with no policy note. Underneath that are three runtime defects in the two newest
stages (a variable used before it is set, a GUI selection path that no consumer reads, and an
unset mask diameter), a dead GUI control, and an IPC design that is one setter call away from
undefined behaviour.

The sieve object is properly designed, the IPC framing is mostly careful, and the fixes are
well-defined. Section E gives them.

## A. Runtime defects, ranked

1. **`bestvol` is used uninitialised.** Declared at p03:736
   (`simple_stream_p03_initial_analysis.f90`), assigned only in the success branch (p03:796).
   Both fallback branches print it (p03:800-806) and p03:809 opens
   `recvol_state<garbage>.mrc`. Any project without `state`/`proj` in `os_cls3D` takes that path.

2. **The GUI reference-selection path is disconnected.** `save_pickrefs_selection`
   (p03:1372-1412) writes the *selected* references to `deselected_references.mrc`
   (`STREAM_DESELECTED_REFS//MRC_EXT`) and the main loop exits (p03:503). The consumers, refpick
   and sieve, wait on `../opening2D/selected_references.mrcs`
   (`simple_stream_p00_master.f90:1313` and 1345), which only `finish_abinitio3D` writes
   (p03:851-852). A user who selects references before the 3D route finishes gets a stage that
   terminates and a refpick that waits 24 hours and then throws
   (`simple_stream_p04_refpick_extract_new.f90:152-166`). The p03 header (lines 17-30) still
   documents the user-selection design.

3. **p07 never learns the mask diameter.** The quality gate receives `params%mskdiam`
   (`simple_stream_p07_abinitio3D_multistate.f90:245`); nothing on p07's command line sets it
   (the master builds that line from scratch, p00:1396-1413), so it is the default `0.`
   (`simple_parameters.f90:616`): a zero mask radius in `model_cavgs_rejection`. `abinitio3D`
   then gets the literal `5000` A (p07:153). The pool already knows the value
   (`get_mskdiam('cavg', ...)`, as p06 does at `simple_stream_p06_pool2D_new.f90:209`).

4. **`increase_nmics` is dead.** The master forwards it (p00:514-515); no stage reads it; p03's
   receive loop handles only `pickrefs_selection` (p03:497-506).

5. **Signal handling is inconsistent and partly unsafe.** p03's handler calls `exit(0)` inside
   the handler (p03:1673-1676): a SIGTERM during `spproj%write` truncates a project file. The
   master's SIGINT handler does `pthread_join` plus `exit` in-handler (p00:741-751). p05, p06 and
   p07 use the correct volatile-flag pattern.

6. **The IPC wire format is the raw memory image, and one type now has an allocatable
   component.** `serialise` is `transfer(self, buffer)` with `sizeof(self)`
   (`src/utils/gui/metadata/simple_gui_metadata_base.f90:102-109`); the master reconstitutes
   with `transfer` (p00:769-1054). `gui_metadata_vol3D%reprojtiles` is `allocatable`
   (`simple_gui_metadata_vol3D.f90:83`, setter at 178). Today no sender calls `set_reprojtiles`
   (p07 ships tiles as separate messages), so the descriptor crosses the fork unallocated and the
   copy is a no-op; the first caller of that setter hands the master a pointer into a child's
   heap. Related: p07's pipe writer diverges from the other six copies (`MAX_RETRIES = 3000`,
   and it abandons mid-frame, p07:589-601), and the master grew a "drop and resync" branch to
   survive it (p00:1157-1166). The protocol was patched at the reader instead of fixed at the
   writer.

## B. Methodological choices to push back on

### B1. Class-average duplication (`balance_classes`, p03:1172-1331)

Selected class averages are replicated to `TARGET_NCLS = 501` in proportion to population and
written to `cavgs_balanced.mrc`; `os_cls2D` is rebuilt with `n_balanced` rows (p03:1253-1263),
`os_cls3D` overwritten from it (1264-1269), `os_out` re-pointed, the even/odd stacks duplicated
(1278-1289), and the sigma2 companion copied without expansion (1291-1301).

Consequences: population becomes an integer replication count (rounding, and a `pop`/`state`
mismatch between rows and particles); the `ptcl2D%class -> os_cls2D` index invariant is broken
for the cycle project; and `abinitio3D_cavgs` treats every row as an independent observation,
so its even/odd statistics see duplicated "particles" and the resolution estimate is
overconfident. The sigma2 file is dead rather than harmful: `abinitio3D_cavgs` deletes inherited
sigma2 state and builds its own (`simple_commanders_abinitio.f90:222-230`).

The honest fix is a per-row weight where the commander builds the class-average rows
(`simple_commanders_abinitio.f90:235-251` sets only `class`, `state`, `eo`, `stkind`) and the
cavgs reconstruction honouring it. If the real driver was a minimum row count for sampling, that
limit belongs in the commander, not in a stream workaround. Introduced 11 June (dd1450261) and
7 August (6f2f34dbb). See E2.

### B2. Picking references are reprojections of a 3-state ab initio volume

The state with the most distinct `proj` indices is chosen (p03:774-807), reprojected at
`nspace=50`, rescaled and written as `selected_references.mrcs` (p03:816-863). The idea (a
template bank from a 3D model) is defensible; the implementation is not validated and silently
replaced the documented human-in-the-loop design (`doc/algorithms/streaming_pipeline.md`,
"Bootstrap"; commits 0b63eab5d "auto reference generation enabled" and bc63531a0). Hardcoded:
`nstates=3`, `nstages=3`, `nrestarts_collapse=3`, `lpstop=8`, `lpstart_ini3D=100`, `nthr=16`.
`mkdir` was left at its default, so the stage has to hunt the highest-numbered
`<n>_abinitio3D_cavgs` restart directory afterwards (`find_final_abinitio3D_cavgs_dir`,
p03:655-677). And it consumes the duplicated stack from B1. See E3.

### B3. The 3D ingestion contract treats provisional 2D labels as final

p06 exports a delta of not-yet-exported particles each pool iteration after iteration 5
(p06:367-376; the `export` branch of `write_project_stream2D` in
`simple_stream_cluster2D_utils.f90`), with a sticky `exported` flag. Everything the pool does to
those particles afterwards, class reassignment, class rejection, the user's sieve-reference
deselection (`update_match_class_states`, p06:320), never reaches 3D.

p07 then re-runs `model_cavgs_rejection` (pool preset) per delta against the full pool class
averages (p07:238-264), so a class accepted in delta 1 can be rejected in delta 3 and the 3D pool
carries inconsistent particle states; `os_cls2D` is seeded from the first delta and never
refreshed (p07:280-282). `abinitio3D` runs once, on the first delta, with ingestion paused, and
`NSTATES3D = 3` is fixed. How the model is revised as the pool grows is left to the add-on
route, which is under construction and outside this review.

The `exported` flag was added to the core `ori` parameter enum (slot 26,
`src/defs/simple_defs_ori.f90`): a stream lifecycle bit in every particle record.

This contract is the thing to decide on paper before more code goes in. See E4.

### B4. Rejection methodology was iterated on master

Between 6 and 22 July: hard gates, then consensus/overfit clustering, then the Feret-based
`class_compatibility` model, then the "trained neighbourhood model" preset, then "removed
initialisation of compatibility model from references". Along the way
`simple_microchunked2D_fast.f90` (2,449 lines, created 10 June, deleted 9 July) and
`simple_cavg_compatibility_analysis.f90` (2,215 lines) were created and deleted. The commit
messages say it: "under development so results are not guaranteed", "sync pre testing", "lots
of redundant code needs removing".

In `run_cavg_quality_selection_2` the compatibility model is trained on and then applied to the
same batch (p03:1068-1076). Acceptable as an outlier filter, but it should be named as one.

### B5. Experimental behaviour lives in the production stages

Gain-flip auto-detection and gain generation from movies
(`simple_stream_p01_preprocess_new.f90:210-224`); the `L_ITERATION_SNAPSHOTS` compile-time flag
(p06:66); "secret" options (`dir_preprocess`, `thres`); and optics-group numbering derived from
a GUI display id (`optics_id_offset = (nicedispid-1)*500`, p06:167), a presentation identifier
deciding a scientific grouping key.

## C. Architectural drift

- **Stage shape.** One `execute` per stage (p01: 526 lines before the first internal
  procedure, p03: 530, p04: 320, p06: 310) with 10-27 host-associated internal procedures sharing
  dozens of locals. Nothing is unit-testable. The handover note
  (`stream_area_tests_handover.md`) records that the stream tests could only check counts and
  files, and `simple_stream_tester` tests helpers, not stages.
- **Copy-paste.** Seven `send_to_*_in_pipe` copies (about 65 lines each, p01-p07) plus the
  master's writer; three `receive_from_*_out_pipe`; twelve near-identical allocate-or-resize
  `case` blocks in the listener (p00:800-1055); `send_micrograph_meta` and
  `send_micrograph_meta_part` (p03:1465 and 1509); `send_available_cavgs2D` and
  `send_selected_pickrefs`. `write_project_stream2D` has 17 optional arguments and the comment
  "NEEDS A TYPE DEFINITION FOR ARGUMENTS !!!!!!!!!!!!!" (`simple_stream_cluster2D_utils.f90:523`).
- **Resources.** 53 literal `nthr`/`nparts`/`nchunks`/`ncls` settings in the stream directory.
  The master alone hands stages 4/1/32/8/16/8/8 threads plus their children (abinitio2D 16,
  abinitio3D 4x16, sieve 4x16, pool 6 parts), with `params%ncunits = 8 ! set to 8 for now`
  (p00:302). No relation to the host or to `params%nthr`.
- **`parameters` mutated after `new`** in about 45 sites: p01 gain and threshold fields; p03
  `params_sieve = params` followed by edits (p03:286-294); p05; p06 `smpd`/`box`/`mskdiam`; the
  master's `qsys_name`/`ncunits`. This bypasses validation, and the sieve keeps a third full
  copy (`simple_ptcl_sieve.f90:151`). CLAUDE.md asks for command-line normalisation before
  `params%new`.
- **Imports.** Twelve stream modules `use` without `only:`; `simple_stream_api` re-exports the
  global-state module `simple_stream2D_state` into every stage (a pre-existing pattern, but p03
  and p07 now inherit it needlessly); p03 carries about 12 unused imports (p03:41-66: six
  commanders, `abinitio_rec_fbody`, `merge_selected_project_files`, `MSK_EXP_FAC`,
  `BOX_EXP_FAC`, `COSMSKHALFWIDTH`, `PREPROCESS_MORPH_SIZE`).
- **Layering.** `src/utils/gui/simple_gui_communicator.f90` uses `simple_parameters` and
  `simple_sp_project` (utils depending on main; it extends an inversion `simple_nice.f90`
  already had) and is wired into nine batch commanders (`gui_comm%new`/`add_metadata`/`kill` in
  refine3D, abinitio, abinitio2D, preprocess, motion, pick, project_*). The vendored
  json-fortran was patched (+413 lines, ab5d1d9c7), an upgrade hazard.
- **Dead code and naming.** `micimporter` (p03:877), `run_cavg_quality_selection` v1 (p03:945),
  `run_cavg_size_selection` (p03:1098), `process_selected_refs`/`process_selected_refs_2` and
  `wait_for_folder` v1 in `simple_stream_utils.f90`, `nptcls_glob` in p07. p03 answers to
  `gen_pickrefs`/`initial_analysis`/`opening2D`/`OPENING2D_JOB_NAME`/"SIMPLE_GEN_PICKREFS"; the
  `_new` suffixes now name the only versions; `ipc_pipe_abinitio3D_multstate`.
- **Docs and skills are behind the code.** `doc/algorithms/streaming_pipeline.md` describes six
  stages and class averages as references; `.claude/skills/simple-main-stream` lists files that
  no longer exist; the `simple-microchunk-rejection` skill and the CLAUDE.md routing rule point
  at code and docs deleted on 9 July. Joe's own `ptcl_sieve_policy.md` and
  `class_compatibility_policy.md` are good and should be the template for the missing
  3D-ingestion and reference-generation policies.

## D. What is good

The message-queue-to-pipes migration (Mac support); the master's memory-leak and
allocation-churn work; EINTR/EAGAIN-aware framing with a length prefix; volatile-flag SIGTERM
handling in most stages; sentinel-driven restart recovery; and above all `ptcl_sieve`, a real
type with `new`/`kill`, a policy doc and a tester. That is the pattern the stages should follow.

## E. Recommendations

Each item names the change, where it goes, and how to know it is done. E1 and E9 are small and
independent; E2-E4 are the methodological decisions; E5-E8 are the structural work; E10-E12
keep it from recurring. A suggested order closes the section.

### E1. Immediate fixes (days)

1. **`bestvol` (A1).** Initialise it, and make the fallback explicit: when `os_cls3D` cannot be
   ranked, choose the state with the largest `os_cls3D` population among those whose
   `recvol_state<nn>.mrc` exists, and `THROW_HARD` when none does. Move the ranking out of
   `finish_abinitio3D` into a pure function (`rank_states_by_view_coverage(os_cls3D)`) in a stream
   utils module so the `lib_stream` suite can test it on a synthetic `oris`.
2. **Selection path (A2).** One routine writes the references file both routes produce, under
   the name the consumers wait for (`STREAM_SELECTED_REFS//STK_EXT`; replace the literal
   `'selected_references.mrcs'` in p03 and p00 with the constant). Decide precedence and encode
   it: a user selection wins over a running 3D route (terminate or ignore the `abinitio3D` job
   and do not overwrite the user's file), and p03 exits its loop only after the consumers' file
   exists. Add a `lib_stream` test that drives `save_pickrefs_selection` on a fixture and asserts
   the consumer path exists with the selected count.
3. **Mask diameter (A3).** p06 stamps `mskdiam` into the exported project's `os_out` cavg entry;
   p07 reads it from the first set with `get_mskdiam('cavg', ...)` and passes it to both the
   quality gate and `abinitio3D`. Delete the literal `5000`. Fail fast in p07 when `mskdiam <= 0`.
4. **`increase_nmics` (A4).** Remove it from the master, the update metadata type and NICE
   until the p03 plan (`NMICS_PLAN`) is a parameter; or implement step 5a of the p03 header.
   Do not leave a control the GUI offers and nothing consumes.
5. **Signal handlers (A5).** p03 and the master set the volatile flag like p05-p07; the master
   joins the listener in its normal shutdown path, never inside the handler.
6. **IPC memory image (A6).** Remove `reprojtiles` from `gui_metadata_vol3D` (p07 already sends
   tiles as separate `cavg2D` messages) or make it a fixed-size array with a count. Add a
   round-trip test to `simple_gui_metadata_tester` for every transferable type (`serialise`,
   `transfer` back, compare), and state the rule in `src/utils/gui/metadata/metadata.inf`: no
   allocatable or pointer components in a type that crosses a pipe. Restore p07's writer to the
   shared semantics (drop only whole frames, never abandon mid-frame) as part of E5, and then
   remove the master's resync branch, which only exists to cover it.

### E2. Population weights instead of duplication (B1)

- In `abinitio3D_cavgs`, where the class-average rows are built
  (`simple_commanders_abinitio.f90:235-251`), set a per-row weight from `os_cls2D` `pop`:
  `w = pop / mean(pop over selected classes)`, identical on the even and odd rows, so the mean
  weight is 1 and the effective number of observations is unchanged.
- Have the cavgs reconstruction honour that weight at insertion. If the reconstruction path does
  not yet read a per-row weight, adding it is part of this change and belongs in
  `src/main/volume` and the reconstruction strategy, not in stream code.
- Delete `balance_classes`, `duplicate_balanced_stack` and `update_os_out_stk` from p03; the
  quality-selected cycle project goes to `abinitio3D_cavgs` unchanged, and the
  `ptcl2D%class -> os_cls2D` invariant holds again.
- If the sampling scheme needs a minimum number of rows, enforce it in the commander with a
  clear error or a documented fallback, and say so in the program's UI description.
- Validation: one dataset, references produced with weights versus with duplication; the
  weighted run's even/odd resolution should be lower (it is the honest one) and the picking
  references at least as good by the sieve's acceptance rate.

### E3. Reference-generation policy (B2)

- Write `doc/policies/stream/reference_generation_policy.md` on the `ptcl_sieve_policy.md`
  template: inputs (which cycle's class averages), the automatic route (`abinitio3D_cavgs`,
  state selection, reprojection, rescaling), the manual route (GUI selection), precedence between
  them, the output contract (`selected_references.mrcs`, `STREAM_NMICS`), and every parameter with
  its default and where it is set.
- Make `nstates`, `nstages`, `lpstop`, `nspace` inputs of `gen_pickrefs` in
  `src/main/ui/simple/simple_ui_stream.f90` (advanced visibility) and read them from `params` in
  p03 instead of literals.
- Make the state choice auditable: log coverage and population for every state, and prefer a
  criterion that combines view coverage with the even/odd state-agreement score the cavgs route
  already computes (`simple_commanders_abinitio.f90:506`) over "most distinct `proj`" alone.
- Have the conditional-restarts commander publish its final result at a fixed path (a sentinel
  or a copy of the final project beside the numbered directories) so callers never scan for the
  highest `<n>_abinitio3D_cavgs`; delete `find_final_abinitio3D_cavgs_dir`.
- Keep the manual route working (E1.2) and make the automatic route the default only once the
  validation in E2 has been run.

### E4. 3D ingestion contract (B3)

Write `doc/policies/stream/stream_3D_ingestion_policy.md` before p07 grows. Proposed contract:

- **The pool publishes state, not deltas.** At each export point p06 publishes its current
  project (mic, stk, ptcl2D, cls2D, out) under `spprojs_completed/<id>.simple`; this is a copy of
  the project `write_project_stream2D` already writes, so the cost is a file copy. p07 mirrors
  the pool's particle index space: on each import it updates `class`, `state` and shifts in
  place for rows it already holds and appends new rows. The `exported` flag becomes
  unnecessary; remove it from `simple_defs_ori.f90`.
- If publishing whole projects proves too heavy at the 10M-particle scale, the alternative is
  change records: re-export rows whose `state` or `class` changed since the last export with a
  revision counter, and have p07 apply them. Only worth its bookkeeping if the snapshot route is
  measured to be a problem.
- **One quality decision per pool revision.** The class-average quality gate runs once per
  published pool (on the pool's `os_cls2D`), and its decision propagates as particle state to
  every row, current and past. Not once per delta.
- **Model lifecycle.** `abinitio3D` starts when the pool is stable (the low-pass schedule has
  ended, iteration 20 or `l_sieve_final`), not on the first export; the policy names the trigger
  for re-deriving the model (a growth factor or a user request) and, when the add-on route is
  finished, how the two relate.
- `nstates` is an input of `abinitio3D_stream` (default 3), not `NSTATES3D`.
- The mask diameter, box and sampling come from the published project, never from literals.

### E5. One IPC module (C, copy-paste)

- Create `simple_stream_ipc` (beside `src/utils/comm/simple_stream_communicator.f90`, or in
  `src/main/stream`) with a `stream_pipe` type holding the descriptors and the rx/tx state, and
  two procedures: `send_framed(pipe, buffer)` with the agreed semantics (length prefix; retry
  budget with `EINTR`/`EAGAIN` handling; a frame is dropped whole or delivered whole), and
  `recv_framed(pipe, buffer) result(got)`.
- Replace the twelve listener `case` blocks (p00:800-1055) with one generic receive into a slot
  array keyed by `get_i`/`get_i_max`.
- `simple_stream_state` holds `stream_pipe` objects instead of raw descriptor pairs.
- Add `simple_stream_ipc_tester` (a `lib_stream` sub-suite): round-trip through a real pipe,
  including a full-pipe partial-write case and a multi-frame read.
- Delete the seven `send_to_*_in_pipe` and three `receive_from_*_out_pipe` copies and the
  master's `send_framed_to_pipe`.

### E6. Stage structure and testability (C, stage shape)

- Each stage becomes a type with explicit state (`params`, `qenv`, the global project, the set
  list, stage flags, metadata objects) and small type-bound procedures (`init`, `watch`,
  `import`, `step`, `broadcast`, `terminate`); `execute` is a short loop over them. Start with
  p07 (smallest, newest), then p05 (already backed by `ptcl_sieve`), then p06 (the pool utils
  exist), then p03.
- p03 runs a whole pipeline (segmentation picking, extraction, 2D, sieve, 2D again, 3D,
  reprojection). Split its `init` and `all` cycles into named steps, and move
  "references from a volume" (rank states, reproject, rescale) into a batch program under
  `simple_exec`, following `ui -> exec -> commander`. It is a legitimate operation outside
  streaming and becomes testable that way.
- With the state explicit, `lib_stream` tests call `stage%step` on fixtures without forking; the
  chunk-path test the handover note asks for becomes possible (factor the chunk command-line
  construction out of `simple_stream_chunk2D_utils` as it suggests).

### E7. Resources (C, resources)

- One derivation, `stream_resources` in `simple_stream_utils` or the master, from `params%nthr`,
  `params%nparts` and the existing environment overrides (`SIMPLE_STREAM_PREPROC_NTHR` and
  friends), hands every stage its `nthr`/`nparts`/`nchunks`, and stages pass those through to
  the command lines they build. Log the total budget once at start-up.
- Remove the 53 literals and `params%ncunits = 8`.

### E8. `parameters` discipline (C, mutation)

- Normalise the command line before `params%new`, as CLAUDE.md requires; delete the post-hoc
  assignments that only restate defaults.
- Where a value is only known at run time (`smpd`, `box`, `mskdiam` from the first import), set
  it in one helper that updates `cline` and `params` together, so the derived fields stay
  consistent, and call it from one place per stage.
- `ptcl_sieve` takes a small `sieve_config` type (the dozen fields it uses) instead of a copy of
  `parameters`; p03 then stops building `params_sieve` by copy-and-edit.

### E9. Dead code, imports, naming (C)

- Delete: `micimporter`, `run_cavg_quality_selection` v1, `run_cavg_size_selection`,
  `process_selected_refs`/`process_selected_refs_2`, `wait_for_folder` v1, the unused imports in
  p03, `nptcls_glob` in p07, `L_ITERATION_SNAPSHOTS` (make it a parameter or drop it),
  `send_micrograph_meta_part` (pass the project as an argument), and the
  `send_available_cavgs2D`/`send_selected_pickrefs` duplicate.
- Rename the `_new` modules and files to their plain names; one name per stage (module, program,
  job name, log banner, metadata type); fix `multstate`.
- Split `write_project_stream2D` into `write_pool_project`, `write_pool_snapshot` and
  `publish_pool_for_3D`, with an options type instead of 17 optionals.
- `use ..., only:` everywhere; stop re-exporting `simple_stream2D_state` through
  `simple_stream_api` (stages that need pool state import it explicitly).

### E10. Layering (C, layering)

- `src/utils/gui` must not depend on `src/main`. Move `simple_gui_communicator` and
  `simple_nice` to a `src/main/gui` (or `src/main/comm`) directory; keep only the metadata types
  and JSON helpers in utils.
- The batch commanders keep exactly three calls each (`gui_comm%new`, `add_metadata`, `kill`),
  and `gui_communicator` gates everything on `niceserver` being set, so GUI concerns cannot grow
  inside commanders.
- Move the fast JSON-to-string routine out of `src/extlibs/json` into a wrapper module
  (`src/utils/json_utils`) so the vendored library stays pristine, or upstream it.
- Remove `exported` from `simple_defs_ori.f90` when E4 lands; stream lifecycle bits do not belong
  in the core parameter table.

### E11. Docs and skills

- `doc/algorithms/streaming_pipeline.md`: seven stages, the reference-generation route, the 3D
  ingestion contract.
- `.claude/skills/simple-main-stream/SKILL.md`: the real file list and the read-first order;
  delete `simple-microchunk-rejection` and its CLAUDE.md routing bullet.
- New policies under `doc/policies/stream/`: reference generation (E3), 3D ingestion (E4), and
  IPC (framing rules, the no-allocatable rule, retry semantics), each with its tester named.

### E12. Process

- Stream work lands from branches through review, with the 24 September test environment as the
  gate: `lib_stream` and the `unit_project` sieve suite green, plus a stage-level test for
  whichever stage was touched.
- Experimental methodology (rejection models, gain generation, automatic references) ships
  behind an explicit parameter whose default is the validated behaviour, with a policy note,
  before it reaches master.
- Commit messages state what changed and why; no "sync" commits on master.

### Suggested order

1. E1 and E9 (fixes and deletions): small, independent, immediately reduce risk and noise.
2. E2 (weights) with E3's policy: the methodological correction the rest depends on.
3. E5 (IPC module): the first shared module, and it removes A6 for good.
4. E4 (3D contract) before any further p07 work.
5. E6 stage by stage (p07, p05, p06, p03), E7 and E8 alongside.
6. E10 and E11 as each area is touched; E12 from now.
