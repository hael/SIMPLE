# GUI metadata and assembler refactoring plan, 6 October 2026

This plan implements `gui_area_review_2026-10-06.md` (the review). Finding IDs (D, C, T, H) are
the review's. It has five phases, and each one builds and tests on its own, so the work can stop
after any of them. Phases 1 to 4 stay inside the existing modules; phase 5 splits the project
type.

The goal, besides the defects: **the GUI metadata types and the assembler import nothing that
reaches `src/main`**, so that stages and the master can use them cheaply (review C4). The folder
stays where it is: no rule ties layers to folders, and moving the folder gains no compile time.

## Decisions

1. **The project type (C3): a record plus a builder, as the last phase.**
   `gui_metadata_project` keeps its fields, its stage slots and its JSON, and takes plain values
   through setters. A new stateless `simple_gui_project_builder`, beside the communicator (its
   only caller), reads the `sp_project`, makes the previews and fills the record.
2. **The orientation binning: a new small module in `src/main/ori`.**
   - `simple_oris_utils` holds `oridist_from_oris(os, state, hist)`. The bins follow the shape of
     `hist`, and `gui_metadata_vol3D` makes its grid sizes public. `src/main/ori` therefore
     imports nothing from the GUI.
   - It is not in `simple_ori_utils`, which cannot import `simple_oris`: `simple_oris` imports
     `simple_ori_api`, which imports `simple_ori_utils`.
3. **`user_input`: removed where SIMPLE never sets it true.** That is the initial-analysis, sieve
   and multistate types. NICE tests `jobstats.user_input is True`, so a missing key reads as false
   and its panels behave as now. Only the pool's `user_input`, which p06 sets, stays.
4. **The layering rule is transitive.** The stream heartbeat takes plain stage-status records in
   place of `class(forked_process)`. `simple_forked_process` imports `simple_cmdline`, and through
   it `simple_ui`, so the assembler, and the master's meta store with it, waits today for the
   whole UI layer.

Defaults taken without a question (each can be changed when its step comes):
- **D2.** A snapshot name must be a bare file name ending in `.simple`, with no `/`, and at most
  128 characters long. NICE sends `snapshot_<id>.simple` (`nice_lite/data_structures/streamjob.py:469`).
- **D3.** The windowed time plots widen their window to fit, so they keep the whole run at a
  coarser step. The rate plot keeps its newest points, labelled with their hour.
- **D4.** The particle montage keeps the process id in its final names.
- **Empty `reprojtiles` arrays are left out.** NICE tests `{% if state.reprojtiles %}`, which
  treats an empty list and a missing key alike.
- **`initial_ref_selection` is removed.** Nothing in `nice/` reads it.

## Phase 1. Defects from outside input

1. **D1, the missing class-average stack** (PROJ:382-383). When `get_cavgs_stk` returns no stack,
   or a count different from `os_cls2D`'s, warn and leave out the 2D section instead of stopping.
   Phase 5 moves this code into the builder.
2. **D2, the snapshot name** (GUICMD:90-103, S-UPD).
   - Add a public `MAX_SNAPSHOT2D_FNAME_LEN = 128` to S-UPD, beside `MAX_SNAPSHOT2D_SELECTION`,
     and make the name field that long.
   - In `parse`, drop a snapshot request whose name breaks the rule (decision defaults), with a
     warning beside `warn_dropped`. The rest of the answer is still applied.
   - Test: in the master tester's GUI answers, a name with `/`, an overlong name and a name
     without `.simple` are each dropped, and the rest of the answer is applied.
   - IPC policy §6: add the rule to the `snapshot2D` row and to §6.3.
3. **D3, the time plots** (TPLOT, `simple_stream_meta_plots`).
   - TPLOT gets a public `MAX_TIMEPLOT_POINTS = 512`.
   - `set_timeplot_from_oris` uses `max(window, ceiling(real(n) / MAX_TIMEPLOT_POINTS))`
     micrographs per point. Its labels stay the point index.
   - `set_rate_timeplot` sends the newest `MAX_TIMEPLOT_POINTS` rates, labelled with their hour.
   - Test: a new `simple_stream_meta_plots_tester` covers a run longer than the capacity and a
     rate history longer than the capacity. Register it as the `unit_stream` sub-suite "meta
     plots", in the test commander and in the test UI's list of sub-suites.
4. **D5, the JSON root.** In ASM:39, write `type(json_value), pointer :: json_root => null()`.

## Phase 2. Imports and the layering rule

1. **The orientation binning (C4).**
   - New file `src/main/ori/simple_oris_utils.f90`, a stateless module.
     `subroutine oridist_from_oris(os, state, hist)` takes `class(oris), intent(in)` and the state,
     and returns `integer, intent(out) :: hist(:,:)`. The azimuth runs over -180..180 in
     `size(hist,1)` bins and the elevation over -90..90 in `size(hist,2)` bins. Its imports are
     `simple_oris` and `simple_linalg` (`rad2deg`).
   - VOL loses `set_oridist_from_oris` and its imports of `simple_oris` and `simple_linalg`. It
     makes `ORIDIST_NBINS_X` and `ORIDIST_NBINS_Y` public.
   - The callers declare `hist(ORIDIST_NBINS_X, ORIDIST_NBINS_Y)`, call `oridist_from_oris`, then
     `set_oridist`. There are two: p07 (`simple_stream_stage_solve3D.f90:1510`, and its header at
     :29) and PROJ:577, which phase 5 moves into the builder.
   - Test: in the oris tester, orientations with known normals land in the expected bins, and
     particles of other states are ignored.
   - Replace the comment at `simple_stream_meta_plots.f90:10-11`: its helpers live on the stream
     side because the metadata types import nothing from `src/main`.
2. **The heartbeat records (decision 4).**
   - In ASM, add a public `type gui_stage_status` with `pid`, `queuetime`, `starttime`,
     `failtime`, `stoptime` and `status`. Give the status public constants
     `GUI_STAGE_STATUS_{RUNNING,FAILED,FINISHED,SKIPPED,UNKNOWN}`.
   - `assemble_stream_heartbeat` takes seven `type(gui_stage_status), intent(in)` arguments, one
     per stage, plus `n_active_persistent_workers`. It keeps the JSON keys, including
     `initial_picking` and `opening2D`, which both come from p03. ASM no longer imports
     `simple_forked_process`.
   - `simple_stream_master_stage` gets a public `fork_gui_status(fork) result(st)` for a
     `class(forked_process)`. It maps `FORK_STATUS_*` onto the record. p00's `send_heartbeat`
     passes `fork_gui_status(shared%stages(id)%fork)`.
   - Tests, split by what they need:
     - ASMT checks the running and finished payloads from plain records, with
       `heartbeat_mismatch`, in `run_all_gui_assembler_tests`. It needs no forks.
     - The live-fork test moves to the master tester and checks `fork_gui_status` on running and
       stopped children. It stays registered as "stream heartbeat" in the forked_process entry
       (`simple_commanders_test_class.f90:49,745`), now imported from the master tester.
3. **The assembler's imports (C1).**
   - ASM imports from the defining modules with `only:`: json, `CK`, string, unix, the base type,
     the display types and the stream types in its signatures.
   - `assemble_batch_metadata` takes `class(gui_metadata_base)`, so ASM imports neither the
     project type nor the umbrella.
   - UTILS and the master's meta store import from the defining modules too. The testers may keep
     the umbrella, which stays as an existing umbrella.
4. **Other imports (C2).** Remove the unused imports the review lists. Remove `use json_kinds`
   from the six modules that use nothing from it.
5. **Profile.** Run a clean `scripts/profile_build.sh` before phase 2 and again after it, and
   record the result in this plan. Claim a compile-time gain only from that pair.

## Phase 3. The pipe contract and dead state

1. **T1, the vol3D tiles.** A test sets tiles with `set_reprojtiles`, serialises the volume,
   receives it by `transfer`, and checks that the copy has no tiles and the same other fields.
2. **T3, the multistate status.** Add set/get and `jsonise` tests, and put its tests in order in
   `run_all_gui_metadata_tests`.
3. **`initial_ref_selection`.** Remove the field, `set_initial_ref_selection`, the JSON key and
   their tests from S-POOL. In the header, drop the sentence that names it.
4. **`user_input` (decision 3).** Remove the field, `set_user_input`, the JSON key and the `get`
   argument from S-IA, S-SIEVE and S-3D. That reaches:
   - the callers p05 (`simple_stream_stage_sieve.f90:300`) and p07
     (`simple_stream_stage_solve3D.f90:332`);
   - the tests MDT:1248-1260 and 1330-1356, `simple_stream_stage_sieve_tester.f90:299-320` and
     `simple_stream_stage_solve3D_tester.f90:490-509`.

   S-POOL keeps its `user_input`.
5. **C6.** Give PROJ a `serialise` override that stops with "gui_metadata_project holds
   allocatable components and is jsonised in-process only".
6. **C7, the multistate guard.** S-3D's `set` stops when `nstates` exceeds
   `MAX_STATES_SOLVE3D_MULTISTATE`. p07 already checks the bound on its command line.
7. **VOL's `l_cfar`.** Drop it: `jsonise` and `get` both use the documented rule, `cfar` when
   positive.
8. **IPC policy §9.** Remove the project type gap (closed by step 5) and the
   `initial_ref_selection` gap (closed by step 3).

## Phase 4. Assembler and metadata hygiene

1. **T2 first.** Before the assembler changes, give ASMT checks of fixed fields through
   `json%get`: the multistate section with tiles nested per state, the project section's counts,
   and an unchanged section being suppressed.
2. **H3, the assembler's duplication.**
   - Two private helpers replace the repeated blocks. `add_assigned_array(self, parent, key,
     items)` takes `class(gui_metadata_base), intent(inout) :: items(:)`, and
     `commit_section(self, json_ptr, hash)` does the hashing.
   - The per-state tile filter stays as its own loop, and an empty `reprojtiles` array is left
     out.
3. **C5, named capacities.** Make public `MAX_HISTOGRAM_BINS` (HIST), `MAX_SIEVE_SELECTION`
   (S-SIEVE), `MAX_FSC_VOL3D` (VOL) and `MAX_OPTICS_SHIFTS` (OG). The last replaces
   `get_max_points`, whose only caller is `simple_stream_stage_optics.f90:333`.
4. **C7, intents and alignment.**
   - `jsonise` becomes `intent(in)` in BASE and in every override, and the assembler's metadata
     arguments follow. So does MDT's `json_text`.
   - The getters of S-IA, S-PP and S-OA become `intent(in)`.
   - Realign the columns the rename left out of line, and rename `meta_initial_picking` in UTILS.
5. **H2, the comments.** Fix the stale or wrong comments the review lists: S-UPD:4 and 31, VOL:4
   and 41, PROJ:160 and 200, ASM:44 and 468.
6. **H4, the reprojection file name.** Add `refine3D_reprojs_fname(state)` to
   `simple_refine3D_fnames`, beside `refine3D_oris_heatmap_fname`. Use it in the writer
   (`simple_solve3D_utils.f90:966`) and the two readers (`simple_stream_stage_solve3D.f90:1474`
   and PROJ:493, which phase 5 moves).

## Phase 5. The project type as a record plus a builder (C3)

1. **The record.** PROJ keeps its fields, the two stage-slot types and `jsonise`, and drops
   `set(spproj, ...)`. One setter per section takes plain values:
   - `set_project(projname, projfile)`, which stamps `created` the first time;
   - `set_movies(movies)`;
   - `set_micrographs(micrographs, nmics, nmics_selected, xdim, ydim, smpd, pspec_size)`;
   - `set_particles(ptcls_jpg, ptcls, nstks, nptcls)`;
   - `set_cavgs2D(stage, cavgs, ncls2D, ncls2D_selected, nptcls_selected, dim_cavgs, mskdiam,
     mskscale)`, where stage 0 is the final slot;
   - `set_vols3D(stage, vols, nstates3D)`.

   The final-slot reuse and the growth of the stage arrays stay in the record. PROJ then imports
   only json, defs, string, `simple_string_utils`, error, unix and the metadata modules.
2. **The builder.**
   - New stateless `src/utils/gui/simple_gui_project_builder.f90`, with
     `build_project_metadata(meta, spproj, oritype, stage, selection)`. It does what `set` does
     today and fills the record through its setters. It may import `src/main`, being neither a
     metadata type nor the assembler.
   - Moved with it are:
     - the fix of D1;
     - D6: tiles carry `mrcpath=volpath_out`;
     - D4: the montage's final names carry the process id;
     - the binning through `simple_oris_utils`;
     - the reprojection name through H4's function;
     - the FSC capacity through `MAX_FSC_VOL3D`.
   - Dropped along the way:
     - the `movthumb` guard (PROJ:170), which never fires and would empty the section if it did;
     - the commented-out lines (PROJ:156-159) and the unread `n_valid_cavgs` counter;
     - the name `micrograph_indices`, which becomes `ptcl_indices`.
3. **The communicator.** `add_metadata_1` and `add_metadata_2` call `build_project_metadata`
   under the metadata mutex, as they call `set` today.
4. **Tests.** ASMT's `test_project` builds its project through the builder, and checks the
   record's fields as T2 does. MDT gains set/get and `jsonise` tests of the record's setters.
5. **The rule, written down.** After phase 5 the metadata types and the assembler import nothing
   that reaches `src/main`. With approval, add this to the change rules of the IPC policy (§8).
   Then propose updates to the `simple-main-stream` skill if its GUI lines need them.

## Checks

- Each phase builds without warnings in the files it touches.
- Tests:
  - `unit_ui`: "GUI metadata" and "GUI assembler";
  - `unit_stream`: "stream master", the sieve, multistate and optics stage tests, and the new
    "meta plots";
  - the forked_process entry: "stream heartbeat";
  - the oris tests;
  - `scripts/check_test_registry.py`.
- The build profiles of phase 2.
- After phases 2 and 5:
  - one stream session from NICE: heartbeat, the stages' panels, a 2D snapshot, the 3D views;
  - one batch job of each kind from NICE: micrographs, particles, 2D classes, 3D, and a
    `selection oritype=cls2D`.

## Behaviour changes

- **D1:** a `selection oritype=cls2D` from NICE on a project without a class-average stack no
  longer fails; the 2D section is left out, with a warning.
- **D2:** a snapshot request with an unsafe name is dropped with a warning. NICE's names are not
  affected.
- **D3:** past their capacity, the windowed time plots coarsen instead of stopping p01. The rate
  plot shows its newest 512 hours.
- **Removed JSON keys:** `user_input` from the initial-analysis, sieve and multistate sections, and
  `initial_ref_selection` from the pool's. NICE reads the first as false, as now, and never reads
  the second.
- **Empty `reprojtiles` arrays are left out**, which NICE treats as before.
- **D4:** the montage files are named `ptcls_sample_<pid>.jpg` and `ptcls_sample_<pid>_lp.jpg`.
- **D6:** the tiles of non-final 3D stages carry an empty `mrcpath`.
- **Unchanged:** the heartbeat JSON, and every other key.

## Out of scope

- `simple_nice` and `simple_guistats`, which are A1 and C13 of the release 4 inventory.
- NICE. Its prompts for the initial-analysis, sieve and multistate panels stay inert, as they are
  today; it could drop them.
- Making the previews before the communicator takes the metadata mutex. Today the communication
  thread waits while a movie is summed. The builder makes this possible, as a later step.
  *Done as a follow-up (6 October 2026):* the builder became a class, `gui_project_builder`. Its
  `build` reads the project and writes the previews into its own state, without the lock; its
  `apply` copies the result into the record, and only `apply` runs under the mutex.
- `simple_forked_process`'s own import of `simple_cmdline`.

## Status

Phases 1 to 5 are implemented (6 October 2026). The debug build (gfortran 15.2) passes `unit_ui`,
`unit_ori`, `unit_stream` and the forked_process entry. The first runs found one fault: the
assembler passed the pointer result of a polymorphic `jsonise()` call straight to `json%add`,
which corrupted the JSON tree; `add_item` now takes it into a local pointer first. Still open:
- the build profiles of phase 2, step 5;
- the checks above.

A script over every module's imports confirms the goal: no metadata module (the project record
included) and not the assembler reaches `src/main` through any chain of imports.

Where the code differs from the steps above:
- **D1, D4, D6 and H4** are not applied to the project type first and then moved: they went
  straight into `simple_gui_project_builder`. Phase 1 therefore does not touch the project type.
- **The assembler's section hashes** are an array indexed by section (`SECTION_*`).
  `commit_section` takes the index rather than a hash component of `self`, which avoids passing
  `self` and one of its components to the same procedure.
- **`add_assigned`** takes `l_object`, since the histograms and time plots go out as objects
  rather than arrays.
- **The record's setters** are `set_project`, `start_update`, `set_summary` (the counts and sizes,
  all optional), `set_movies`, `set_micrographs`, `set_particles`, `set_cavgs2D` and `set_vols3D`.
  `start_update` keeps the old behaviour: each update drops the movie, micrograph and particle
  lists, and stage 1 drops the 2D and 3D stages.
- **`fork_gui_status`** lives in `simple_stream_master_stage`. The stream master's tester checks it
  on one live child.
- **Comments and literal capacities fixed along the way:** in p07, a comment saying that `kill`
  resets only the flags, and the FSC capacity written as the literal 1000.
- **`pool2D_policy.md`** is untouched; IPC policy §4, §6, §8 and §9, and the GUI onboarding page,
  are updated. The rule of phase 5, step 5 is in IPC policy §8 (approved 6 October 2026).
