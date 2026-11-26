# select2D in its own executable, simple_viewer, saving into a new selection folder

Status, 2026-11-25: approved by the owner; Phases 1-3 written, not yet built (Phase 0 and every exit run on the
owner's Mac). Base: `2750c34ac` (clean tree). Supersedes:
`select2D.md` decisions 2 and 4, the `mkdir=yes` part of decision 3, and the executable/router/category part of
decision 1 (its banner says so in the same commit). Owner: Dibyajyoti Maity. Depends on:
- `select2D` as committed in `2750c34ac` (built and run by the owner on the Mac, 2026-11-25; validated: yes);
- `prg=selection oritype=cls2D` (`exec_selection`, `simple_commanders_project_core.f90`) as the reference
  semantics of a class selection and of its output folder (validated: yes, in production).

## 1. Why

`select2D` today is `simple_exec prg=select2D`: `simple_exec` links GTK 3 for one program, and with `mkdir=yes`
`params%new` makes `<n>_select2D/` and copies the project there the moment the program starts, so a session that
is only looked at, or closed without saving, leaves a folder and a project copy behind (a 435 MB `1_select2D/`
was left in the checkout on 2026-11-25). Save then writes four segments into that copy with
`write_segment_inside` and a `selection.txt` beside it. The owner's ruling (2026-11-25): select2D becomes a
separate executable, `simple_viewer`, which for now launches the select2D GTK application and nothing else;
nothing is written until **Save selection**; Save writes a new project into a selection folder named and placed
as `prg=selection` names and places its run folder, through the routines `prg=selection` uses; no
`selection.txt`.

## 2. Contract

- `simple_viewer prg=select2D projfile=<p>.simple` opens the select2D window (`prg=` first and required since
  decision 15); there is no `append=` (decision 14: `simple_viewer` stops on one).
- `simple_exec` no longer offers or routes `select2D` and no longer links GTK (`otool -L` / `ldd` of
  `simple_exec` lists no `gtk`). `simple_exec prg=select2D` stops in `cmdline%parse` with "select2D is not
  executed by simple_exec but by simple_viewer" (no alias; release 4 keeps no compatibility).
- Opening a project and closing the window without a Save creates no directory and writes no file. The input
  project is never written (byte-identical before and after a session, saved or not).
- The first Save produces `<n>_selection/<input basename>.simple` in the input's workspace root (decision 13);
  later Saves of the session rewrite that file. Its segments equal those
  `prg=selection oritype=cls2D infile=<the same 0/1 states> append=no` writes from the same input, except
  the declared exclusions of decision 12 (the folder in the projinfo paths, and the jobproc row that records the
  job: `prg=select2D` here, `prg=selection` from `simple_exec`'s front door there).
- `prg=selection` is unchanged (no edit to `exec_selection`).
- Off-state: no program other than `select2D` changes behaviour; only `simple_viewer` and `simple_test_exec`
  link GTK.

## 3. Decisions and rulings

1. (2026-11-25, owner) select2D runs in a new executable `simple_viewer` only. Its UI entry and router case leave
   `simple_exec`; `production/CMakeLists.txt` links `SIMPLE_GTK_LIBRARIES` to `simple_viewer` (and
   `simple_test_exec`, as now) instead of `simple_exec`. Replaces `select2D.md` decision 2.
2. (2026-11-25, owner: "simple_viewer is the executable that launches the select2D gtk application for now")
   Command line: `simple_viewer projfile=<p>.simple`, no `prg=` (the `[append=yes]` first written here is removed
   by decision 14). `cmdline%parse` gets one case for
   `simple_viewer`: its program is the one program the registry holds for that executable (new
   `get_executable_prg` in `simple_ui.f90`, which stops unless there is exactly one), so the name `select2D`
   lives only in the UI descriptor; every argument is a `key=value` (no program slot in the count of required
   arguments); `prg` is set on the `cmdline` after parsing so `params%new` sees `select2D`; a `prg=` on the line
   is accepted only if it names that program (so a launcher that always writes `prg=` still works), any other
   stops; no argument at all prints the usage and stops, as `simple_exec prg=<x>` alone does; `describe=yes`
   and `usage=yes` work as elsewhere. A second view later is the owner's decision then (this rule then stops
   with "more than one program").
3. (proposed) Front door `production/simple_viewer.f90` (about 35 lines, the shape of `single_exec.f90`):
   `make_ui`, `cline%parse`, the router, git version and timer. It does not call
   `update_job_descriptions_in_project` (layering policy section 8 lists front-door project writes as a deviation
   not to extend; decision 8 records the job instead), `script_exec` (a window is not queued) or the memory
   monitor. One router, `src/main/exec/simple_viewer_views.f90` (`exec_viewer_commander`). No API umbrella
   module (compile-time policy): the front door imports with `only:`.
4. (proposed) UI: the `select2D` descriptor moves from `simple_ui_cavgproc.f90` to a new category module
   `src/main/ui/simple/simple_ui_viewer.f90`, category `('viewer', 'Viewer', 10)`, constructed from
   `simple_ui_simple_group.f90`, executable `'simple_viewer'`. Inputs: `projfile`, required as now (`append`
   removed by decision 14); visibility Standard (`select2D.md` decision 19 stands); the help says that nothing is written
   before Save selection and where the selection goes.
5. (2026-11-25, owner) Nothing on disk at start: the commander sets `mkdir=no` unconditionally (over any
   `mkdir=` given), so `params%new` neither makes a directory nor copies the project. Replaces the `mkdir=yes`
   of `select2D.md` decision 3. The display check before `params%new` (`select2D.md` decision 15) stays.
6. (proposed) Start-up reads only the `cls2D` and `out` segments of the input (`read_segment`), all the window
   uses (`select2D.md` decision 8: rows, `pop`, `res`, `state`, the cavg stack of `out`); the particle segments
   are read only at Save. Today the whole project is read at start.
7. (2026-11-25, owner) The first Save makes `<n>_selection/` in the directory the viewer was started from
   (`params%cwd`), with `n = find_next_int_dir_prefix(<that directory>)`: the name, number and place
   `setup_execution_context` gives `prg=selection`'s run folder under `mkdir=yes` (`int2str(idir)//'_'//prg`).
   Later Saves of the same session rewrite the project in that same folder (one folder per session). The place
   (the launch directory) is replaced by decision 13.
8. (2026-11-25, owner) What a Save writes, all in `cls2D_selection%save`: read the input project in full from
   disk (every Save starts from the input, as `prg=selection` starts from its copy, so repeated Saves do not
   build on each other and the display project stays untouched); stop with `THROW_HARD`, writing
   nothing, if its `cls2D` row count is no longer the session's; `report_state2cls2D(states, append=.false.)`
   (decision 14; the routine `exec_selection`'s cls2D branch calls: `cls2D` `state`, then `map2ptcls_state`,
   which sets the `ptcl2D` and `ptcl3D` `state` of every classified particle to its class's state, whatever its
   state in the input, sets particles without a class to 0, and copies
   the class states into `cls3D` when it has as many rows as `cls2D`; `mic`, `stk`, `out` and `optics` are not
   touched, as in `prg=selection`); `update_projinfo(<dest>)` (projinfo names the new file, as
   `prg=selection`'s read of its copy does); append one job description to `jobproc` (`prg=select2D`, the
   viewer's keys, `projfile=<dest>`); one full `sp_project%write(<dest>)`. No `write_segment_inside`, no
   `selection.txt` (owner), no write to the input. Replaces `select2D.md` decision 4.
9. (2026-11-25, owner) Project file only: the JPEG montages `cls2D_selected.jpg` / `cls2D_deselected.jpg` that
   `exec_selection` writes for NICE are not written, and `exec_selection` is not edited.
10. (proposed) The selection session `cls2D_selection%new(spproj, projfile, append, rootdir, job_descr)`: the
    display project, the input project file, `append`, the directory for the selection folder and the job
    description; new getter `get_saved_projfile()` (empty before the first Save). The commander passes
    `params%projfile`, `params%cwd` and `cline%gen_job_descr(job_descr, string('select2D'))`. (`rootdir` is
    removed by decision 13, `append` by decision 14: now `new(spproj, projfile, job_descr)`.)
11. (proposed) The window: Save selection and the close dialog's Save call the new save; after a Save the
    toolbar shows `<input>  →  <dest>`; the title drops its `*` as now; Discard before any Save logs that nothing
    was written.
12. (proposed) Test: sub-suite `select2D save` of `unit_project`, tester
    `src/main/commanders/simple/simple_select2D_tester.f90`, headless (no GTK call). Fixture in `enter_fixture`:
    6 classes (one with `pop` 0), 40 particles assigned to them (one with `state` 0 in a class that stays
    selected), one stack record, a 16-pixel cavg stack registered in `out`. Checks: (a) after `new` no
    `*_selection` directory exists; (b) the first Save makes `1_selection/<p>.simple`; the input file is
    byte-identical before and after; (c) identity: `commander_selection` run in process from the fixture root
    with `infile=` the session's states, `oritype=cls2D`, `append=no` (it makes `2_selection/`); segments
    `mic`, `stk`, `ptcl2D`, `ptcl3D`, `cls2D`, `cls3D`, `out`, `optics`, `compenv` equal entry by entry;
    `projinfo` equal but for the folder in `projfile`/`cwd`; `jobproc` equal plus one `prg=select2D` row with
    `projfile=<dest>`; (d) a second Save with another selection rewrites `1_selection/<p>.simple` with it and
    makes no `3_selection/`; (e) the state-0 particle in a kept class is at state 1 in every saved project
    (decision 14; first written as "(c) for `append=no` and `yes`"). `make_ui` first, as
    `test_execution_project_copy` does, so `prg=selection` copies its project. No new CTest entry (budget 47
    stays). The `cmdline%parse` case of decision 2 is checked by hand in Phase 2 (no unit suite covers the front
    doors).
13. (2026-11-25, owner) The selection folder goes to the input's workspace root, not the launch directory:
    the folder that holds the input project, or that folder's parent when its name is a numbered run folder
    (`<digits>_<name>`, the form `setup_execution_context` and NICE give run folders; new `is_run_dir_name`
    beside `find_next_int_dir_prefix` in `simple_syslib`). `/ws/workspace.simple` saves to `/ws/<n>_selection/`,
    `/ws/7_solve2D/workspace.simple` to `/ws/8_selection/` when 7 is the highest number in `/ws`. Replaces the
    place of decision 7 (its name, numbering and one folder per session stand); `cls2D_selection%new` derives the
    root from the absolute input path and loses its `rootdir` argument. Decision 12 follows: one run puts the
    input in `7_solve2D/` (first save `8_selection/`, `prg=selection` from the fixture root `9_selection/`), the
    other in a plain `ws/` (`1_selection/` and `2_selection/` inside it); first written as the `append=no` and
    `append=yes` runs, both without append since decision 14. Noted:
    NICE numbers jobs from its database, so it may reuse the number of a folder made here; it refuses a job only
    whose folder name already exists (`8_selection` for a NICE `selection` job numbered 8).
14. (2026-11-25, owner: "We do not want any append route here. Makes no sense") select2D has no `append`: the
    UI input goes, `cls2D_selection` loses its `append` component and argument, and every Save calls
    `report_state2cls2D(states, append=.false.)`, so the class selection alone sets the particle states (every
    classified particle takes its class's state; particles without a class go to 0), as
    `prg=selection oritype=cls2D append=no` does. Because `cmdline%parse` accepts every key SIMPLE knows, the
    commander stops on an `append=` on the command line rather than ignore it (agent's default, not ruled).
    `prg=selection`, `report_state2cls2D` and `map2ptcls_state` keep their `append`. Replaces the `append` of
    decisions 2, 4, 8, 10 and 12, and the append split of the tester runs in decision 13.
15. (2026-11-25, owner: "keep the select2D commander but make it a commander for simple_viewer prg=select2D")
    `simple_viewer` takes `prg=` like `simple_exec`: `simple_viewer prg=select2D projfile=<p>.simple`, with
    `prg=` first. `commander_select2D` and the router `simple_viewer_views` stay. `cmdline%parse` puts
    `simple_viewer` in the `prg=`-required case and takes the program from the first argument, as for the other
    executables; the single-program path of decision 2 (`l_single_prg`, `nprg_args`, the `prg=` check after the
    key loop) and `get_executable_prg` are removed, so `parse` and `simple_ui.f90` are back to their `2750c34ac`
    form but for `simple_viewer` in its executable cases. `simple_viewer` alone stops asking for `prg=` first
    (since decision 16 in `cmdline_err`, as `simple_exec` alone does); `simple_viewer prg=select2D` alone prints
    the usage; a program of another executable stops in the executable check of `parse`. Replaces the command
    line of decision 2.
16. (2026-11-25, owner) `simple_viewer prg=list` lists the viewer programs, and an unknown name stops with
    `Program "<x>" not recognized. Use prg=list to see available programs.`, both in the front door before
    `cline%parse`, as `simple_exec` and `single_exec` do: the front door reads the first argument
    (`cmdline_err`), builds the UI, lists on `prg=list` with the new `list_viewer_prgs_in_ui` of `simple_ui.f90`
    (`print_registered_programs(prgtab, prgnames, 'simple_viewer')`) and looks the name up with `get_prg_ptr`.
    The `prg=list` block inside `parse` is never reached for any executable (the program lookup stops on `list`
    first); its `simple_viewer` case stays a no-op with a comment.
17. (2026-11-25, owner reported "the right click does not work sometimes"; ruled: fix in the grid module, a
    left click on the grid while the popover is up only closes it) The popover is modal (GTK 3's default), so
    while it is up GTK's grab sends every grid click to the popover and none reaches the grid's `on_press`; the
    grid selects no release events, so the popover's own close-on-release likely never fires for grid clicks
    either, and a right-click elsewhere did nothing until Esc or a click outside the grid. `gtkui_image_grid`
    gets `attach_popover(popover)`, which connects the popover's `button-press-event` and
    `button-release-event`. For events whose widget is the grid's drawing area: a context-menu press (right
    click; Ctrl-click on macOS) calls the delegate's `context`, which re-points the open popover at the new
    position in one click; any other single press closes it without touching the selection; the extra events of
    a double or triple click are dropped; the matching release is swallowed so the popover's own handler cannot
    close a popover that has just moved. Presses inside the popover, Esc, Shift+F10/Menu and the popover's
    buttons are unchanged. select2D calls it once in `build_popover`; the popover stays modal.

Defaults taken without a question: the file and module names in decisions 3, 4 and 12; the category id,
display name and order; `get_executable_prg`; the sub-suite name; the folder suffix `'selection'` as a private
constant of `simple_cls2D_selection` (replacing `SELECTION_FNAME`); the toolbar wording; the jobproc row's keys.

## 4. Not in this plan

NICE (whether it can launch `simple_viewer`; it reads the descriptor's executable), other views, any change to
selection, sorting or display behaviour of the window, `exec_selection`, the living policy of `select2D.md`
Phase 6, skill and `CLAUDE.md` edits (proposed in Progress only), GTK 4. Files outside the table are out of
bounds.

## 5. Change types and gates

| Phase | Kind | Gate |
| --- | --- | --- |
| 0 | measurement | incremental build of the base, fast gate, `test=unit_project`, `otool -L simple_exec` lists gtk |
| 1 | new behaviour (Save) | `test=unit_project` with `select2D save`; `test=gtk`; fast gate |
| 2 | new executable, routing | build (UI JSON validation included); fast gate; the checklist of Phase 2 |

## 6. Acceptance

- Folder creation, names, counts, states and the identity of 12(c): exact (*integers and strings*); exclusions
  declared in decision 12 and frozen here.
- Input project unchanged: byte comparison (*exact*).
- `simple_exec` without GTK: no `gtk` line in `otool -L` (Mac) / `ldd` (Linux) (*exact*).
- Reported, not gated: the owner's session on the betagal solve2D project of 2026-11-25.

## 7. Common rules and stop conditions

No new parameter key. No assertion relaxed. Fix as you go, with a file-table row before the edit. Tests at one
OpenMP thread. Nothing is compiled in the agent's sandbox: every build and run is the owner's (Mac, gfortran 16,
with a display); the agent reads the logs.

Stop and report when: a file outside the table needs an edit; 12(c) fails beyond its declared exclusions;
`simple_exec` still needs GTK to link; `prg=selection`'s output would have to change; a key outside `projfile`,
`append`, `mkdir`, `oritype` would be needed. Gates no ruling may waive: the fast gate; `unit_project` with
`select2D save`; the input project byte-unchanged.

## 8. Phases

- **Phase 0, baseline.** On `2750c34ac`: `make -C build -j`, `scripts/run_fast_gate.sh build`,
  `build/bin/simple_test_exec test=unit_project`, `otool -L build/bin/simple_exec | grep -i gtk`. Exit: the
  three results recorded in Progress; base failures named and left.
- **Phase 1, Save into a new selection folder.** Decisions 5-12. Exit:
  `simple_test_exec test=unit_project suite=select2D_save` passes, the whole `unit_project` and `test=gtk`
  pass, fast gate passes.
- **Phase 2, the executable.** Decisions 1-4: front door, CMake, router, UI module, the `cmdline%parse` case and
  `get_executable_prg`; `select2D` out of `simple_ui_cavgproc` and `simple_exec_cavgproc`. Exit: the build
  (with `validate_complete_ui_json`) and fast gate pass; `otool -L build/bin/simple_exec` has no gtk;
  `simple_exec prg=select2D projfile=x.simple` stops naming `simple_viewer`; `simple_viewer` alone prints the
  usage; `simple_viewer prg=other projfile=x.simple` stops; on the owner's betagal project, from a scratch
  directory: `simple_viewer projfile=<p>.simple` opens the window and `ls` shows no new folder; closing
  unchanged leaves none; toggling and Save makes `1_selection/<p>.simple` whose
  `simple_exec prg=print_project_field oritype=cls2D` shows the saved states; a second Save rewrites it and makes
  no second folder; the input's checksum is unchanged; without a display it stops with `select2D.md` decision 15's message and
  makes no folder.
- **Phase 3, records.** Progress entries; `select2D.md` banner naming the superseded decisions. Exit: both
  written.

### File table

| File | What | Phases |
| --- | --- | --- |
| `src/main/project/simple_cls2D_selection.f90` | `new` arguments, `save` to `<n>_selection/`, getter; no `selection.txt` | 1 |
| `src/fileio/simple_syslib.f90`, `simple_syslib_tester.f90` | `is_run_dir_name` and its checks (decision 13) | 1 |
| `src/main/commanders/simple/simple_commanders_select2D.f90` | `mkdir=no`, segment reads, job description, root dir | 1 |
| `src/main/gtk/simple_select2D_window.f90` | toolbar path after Save, Discard message; `attach_popover` call (decision 17) | 1 |
| `src/main/gtk/simple_gtk_image_grid.f90` | `attach_popover` and the popover's press/release handlers (decision 17) | 1 |
| `src/main/commanders/simple/simple_select2D_tester.f90` | new: `select2D save` | 1 |
| `src/main/commanders/test/simple_commanders_test_class.f90` | suite row | 1 |
| `src/main/ui/simple_test/simple_test_ui_class.f90` | `suite=` list of `unit_project` | 1 |
| `production/simple_viewer.f90` | new front door | 2 |
| `production/CMakeLists.txt` | the executable; GTK link moves from `simple_exec` to it | 2 |
| `src/main/exec/simple_viewer_views.f90` | new router | 2 |
| `src/main/exec/simple_exec_cavgproc.f90` | `select2D` case removed | 2 |
| `src/main/ui/simple/simple_ui_viewer.f90` | new category module with `select2D` | 2 |
| `src/main/ui/simple/simple_ui_cavgproc.f90` | `select2D` removed | 2 |
| `src/main/ui/simple_ui_simple_group.f90` | constructs the new category | 2 |
| `src/main/ui/simple_ui.f90` | `get_executable_prg` (removed by decision 15); `list_viewer_prgs_in_ui` (decision 16) | 2 |
| `src/main/simple_cmdline.f90` | the `simple_viewer` case of `parse` (`prg=` required since decision 15) | 2 |
| `doc/implementation_notes/planned/select2D_viewer.md` | this plan, Progress | all |
| `doc/implementation_notes/planned/select2D.md` | banner: superseded decisions | 3 |

## Progress

### Phase 0, 2026-11-25 (pending, owner)

The base `2750c34ac` was built and its `simple_exec prg=select2D` run by the owner on the Mac the same day; the
fast-gate, `unit_project` and `otool -L` numbers of the exit were not recorded before the edits of Phases 1-2,
which are now in the working tree. Compare the Phase 1-2 runs against the last green fast gate of the base, or
record the base by building `2750c34ac` in a separate checkout (not in `build/`).

### Phases 1-3, 2026-11-25, written in the agent's sandbox (no compiler there), not built

Files changed, all in the table: `simple_cls2D_selection.f90` (`new(spproj, projfile, append, rootdir,
job_descr)`, `save` to `<n>_selection/`, `get_saved_projfile`; `selection.txt`, the `spproj` pointer and
`SELECTION_FNAME` gone), `simple_commanders_select2D.f90` (`mkdir=no` forced, `cls2D` and `out` read with
`read_segment`, the job description from `cline%gen_job_descr`), `simple_select2D_window.f90` (toolbar shows
`<input>  →  <saved>` after a Save; Discard message), new `simple_select2D_tester.f90` with sub-suite
`select2D save` (decision 12, `append=no` and `yes`) registered in `suites_project` and the `suite=` list
(`select2d_save`); new `production/simple_viewer.f90`, `src/main/exec/simple_viewer_views.f90`,
`src/main/ui/simple/simple_ui_viewer.f90`; `get_executable_prg` in `simple_ui.f90`; the `simple_viewer` case of
`cmdline%parse` (`nprg_args` replaces the literal 1 of the program slot in the usage and required-argument
counts, so the other executables are unchanged); `simple_ui_cavgproc.f90` and `simple_exec_cavgproc.f90` are
byte-identical to their state before `2750c34ac`; `production/CMakeLists.txt` builds and installs
`simple_viewer` and links GTK to it instead of `simple_exec`.

Static checks run in the sandbox, all passing: `check_test_registry.py`, `check_descr.py`, `check_flex_dag.py`,
`check_guide_imports.py`, `check_cli_chain.py` (598 accepted keys, 1832 UI inputs), `check_submodule_linkage.py`;
no multi-line `THROW_` macro. Exits still to run on the Mac: the build, `scripts/run_fast_gate.sh build`,
`simple_test_exec test=unit_project` (and `suite=select2d_save`), `test=gtk`, `otool -L build/bin/simple_exec`,
and the checklist of Phase 2.

Proposed, not applied (owner approval needed): `simple-main-exec` and `cli_change_chain.md` name the executables
and the rule that a router module goes into `src/main/apis/<exe>_api.f90`; `simple_viewer` has no API module
(compile-time policy) and runs one program without `prg=`. `doc/code_overview/code_base_map.md` and the
generated Fortran indexes do not list `simple_viewer` yet.

What it means, in short: select2D now lives in `simple_viewer`, opens a project without writing anything, and
each Save writes the input with the selection applied, by `prg=selection`'s routine, to one `<n>_selection/` per
session; nothing of it has been compiled yet.

### First build, 2026-11-25 13:26-13:30, owner's Mac (gfortran 16), read by the agent from `build/`

The incremental build compiled every changed file and linked all seven executables, `simple_viewer` included;
`validate_complete_ui_json` passed (`select2D`: executable `simple_viewer`, category `viewer`). The load commands
of the new binaries name `libgtk-3.0.dylib` in `simple_viewer` and `simple_test_exec` and not in `simple_exec`
(the `otool -L` exit, read from the binaries). Fast gate: 14 of 15 suites passed; `unit_project` failed in its
last sub-suite, `select2D save`, after its 23 other sub-suites passed: the first Save worked (`1_selection/`
written), then `prg=selection`, run in process, stopped in `params%new` with "box needs to be larger than 26",
because the fixture's stack record had box 16. Fixed in the tester: the fixture box is 32 (decision 12 said a
16-pixel stack; the floor of 26 is parameters'). The gate's failure stopped the command chain before
`make install`, so `build/bin` still held the 11:27 binaries and the later commands of that run used them.
Next: rerun the build, gate and install.

### Second build, 2026-11-25 13:39-13:40, owner's Mac, read by the agent from `build/` and the checkout

Build, fast gate and install passed: 15 of 15 fast suites (`build/test_runs/ctest_fast.log`), `select2D save`
94 checks in 0.09 s for `append=no` and `yes` (first save to `1_selection/`, input byte-identical, segments equal
to `prg=selection`'s, second save into the same folder; `build/Testing/Temporary/LastTest.log`). `build/bin` holds
`simple_viewer` (13:40:07). `simple_test_exec test=gtk`: PASS, cairo tile 9 and window 5 assertions
(`SIMPLE_TEST_gtk_20261125/simple_test_gtk_report.txt`, run from the checkout root). Phase 1's exit is met. Of
Phase 2's exit, met: build, UI JSON, fast gate, no GTK in `simple_exec`; not yet recorded: the refusal of
`simple_exec prg=select2D` (printed to the terminal only), `simple_viewer` alone printing the usage, and the
owner's window checklist on the betagal project.

Owner's terminal output, same day: `test=gtk` reports GTK 3.24.52; `simple_exec prg=select2D projfile=x.simple`
stops with "select2D is not executed by simple_exec but by simple_viewer" (`simple_cmdline.f90`, the executable
check of `parse`; the backtrace after it is how every `THROW_HARD` ends). Left of Phase 2's exit: `simple_viewer`
alone printing the usage, and the window checklist.

### Decision 13, 2026-11-25, written in the agent's sandbox, not built

The owner's first window session, run from the checkout root, saved `1_selection/workspace.simple` (2.9 MB) into
`~/SIMPLE` (the launch directory of decision 7); the owner then ruled decision 13. Changed: `is_run_dir_name` in
`simple_syslib.f90` with a check in `simple_syslib_tester.f90` (`unit_core`, sub-suite syslib);
`cls2D_selection%new` without `rootdir`, the root derived from the absolute input path; the commander no longer
passes `params%cwd`; the UI help states the rule; the tester runs `append=no` with the input in `7_solve2D/`
(expects `8_selection/`, `prg=selection` from the root makes `9_selection/`) and `append=yes` with the input in
`ws/` (`1_selection/`, `2_selection/` inside it). Static checks pass. Exits to run on the Mac: build, fast gate
(`unit_core`, `unit_project`), install, and the window checklist with the folder now beside the input's run
folder. The `1_selection/` in the checkout root is the owner's test output and untracked; it is not this
change's to delete.

### Decision 14, 2026-11-25, written in the agent's sandbox, not built

The owner ruled out `append` for select2D. Changed: `simple_ui_viewer.f90` (the `append` input gone; the help
says each particle takes its class's state); `simple_commanders_select2D.f90` (no `append` default; stops on an
`append=`; `sel%new` without it); `simple_cls2D_selection.f90` (no `append` component or argument; `save` calls
`report_state2cls2D(states, append=.false.)`); `simple_select2D_tester.f90` (`test_save(run_layout)`, both
layouts compared with `prg=selection append=no`; the fixture's state-0 particle in kept class 3 must be at
state 1 after each save; fixture folders `select2D_save_run_folder` and `select2D_save_plain_folder`).
`prg=selection` and `report_state2cls2D` are unchanged. Exits to run on the Mac: build, fast gate
(`unit_project`, sub-suite `select2d_save`), install; `simple_viewer projfile=<p>.simple append=yes` stops with
decision 14's message.

### Decision 15, 2026-11-25, written in the agent's sandbox on `595499064`, not built

The owner first asked to remove `commander_select2D`, then ruled to keep it and give `simple_viewer` a `prg=`.
Changed: `simple_cmdline.f90` (`parse` differs from `2750c34ac` only by `simple_viewer` in the `prg=`-required
and `get_prg_ptr` cases and a no-op case in the `prg=list` block); `simple_ui.f90` (`get_executable_prg` and its
`public` gone: byte-identical to `2750c34ac`); `simple_ui_viewer.f90` (help: `simple_viewer prg=select2D
projfile=<project>.simple`); `production/simple_viewer.f90` (`!@descr` and the comment on `parse`). Unchanged:
the commander, the router, CMake, the tester. The build is the proof that nothing else called
`get_executable_prg` (skill `simple-dead-code`). Exits to run on the Mac: build, fast gate, install;
`simple_viewer` alone stops asking for `prg=`; `simple_viewer prg=select2D` prints the usage;
`simple_viewer prg=select2D projfile=<p>.simple` opens the window; `simple_viewer prg=selection projfile=x.simple`
stops naming `simple_exec`.

### Decision 16, 2026-11-25, written in the agent's sandbox on top of decision 15, not built

Changed: `simple_ui.f90` (`list_viewer_prgs_in_ui`, public); `production/simple_viewer.f90` (reads the first
argument with `cmdline_err`, lists on `prg=list`, looks the name up with `get_prg_ptr` before `parse`; the
router's "not recognized" message also points to `prg=list`); `simple_cmdline.f90` (comment of the
`simple_viewer` case in the `prg=list` block only). Exits to run on the Mac: build, fast gate, install;
`simple_viewer prg=list` prints the `Viewer:` heading and `select2D`; `simple_viewer prg=foo` stops with the
"not recognized. Use prg=list" message; `simple_viewer` alone stops in `cmdline_err`; `simple_viewer prg=select2D`
still prints the usage.

### Decision 17, 2026-11-25, written in the agent's sandbox on `0b43d033f`, not built

Changed: `simple_gtk_image_grid.f90` (`attach_popover`; private `on_popover_press` and `on_popover_release` on the
grid's user data; comments of the module, `context` and `on_press`); `simple_select2D_window.f90` (`build_popover`
calls `w%grid%attach_popover(w%popover)`; comment of `show_popover`). gtk-fortran interfaces checked against the
3.24.49 source (`gtk_get_event_widget`, `gtk_popover_popdown`). Not reproduced here (no display, no compiler): the
cause is read from the code and GTK 3's documented popover grab. Exits to run on the Mac: build, fast gate,
`test=gtk`, install; then in the window on betagal: right-click opens the popover; a right-click elsewhere on the
grid moves it there (one click); Ctrl-click does the same; a fast double right-click leaves it open; a left
click on a class while it is up closes it and the footer count does not change; the next left click toggles;
Esc, Sort, Select All and Deselect All still close it; brightness, contrast and the histogram still work inside.
