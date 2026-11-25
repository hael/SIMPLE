# select2D: interactive class-average selection with a GTK 3 window, written to the project

Partly superseded, 2026-11-25, by `select2D_viewer.md`: decisions 2 and 4, the `mkdir=yes` of decision 3 and the
executable, router and category of decision 1 are replaced (select2D runs in its own executable `simple_viewer`,
writes nothing until Save selection, and saves a new project into `<n>_selection/` as `prg=selection` does, without
`selection.txt`). The rest stands.

Status, 2026-11-25: implemented in one pass by the owner's ruling to skip the phased run (Phase 0 ran on the Mac;
no tester or Progress entries for Phases 1-6; the decisions below are what the code follows), built and tried by the
owner. Restructured the same day, before the first build, into shared GTK modules (simple_defs_gtk, simple_gtk_utils,
simple_gtk_draw, simple_gtk_image_grid, simple_gtk_display_panel; the display model renamed simple_stack_display) with
select2D as an app on top; the file names of decision 13 are superseded by them. Decisions 2, 6 and 19 ruled by the
owner on 2026-11-24. Base:
`b374f51b4` (clean tree). Supersedes: none. Owner: Dibyajyoti Maity. Depends on:
- GTK 3 and gtk-fortran 3.24.49 in every build (`cmake/Dependencies.cmake`, static library only; validated: yes,
  builds on the Mac and on Oracle Linux 8 per the owner, 2026-11-24);
- the platform test `simple_test_exec test=gtk` (`src/main/gtk/simple_gtk_tester.f90`; validated: yes on the
  Mac with a display);
- `prg=selection oritype=cls2D` as the reference semantics of a class selection (validated: yes, in production).

## 1. Why

Class averages are selected today with `select2d.py`, a 2 075-line Tkinter script outside SIMPLE: it reads the
highest `cavgs_iter<N>.mrc` and `classdoc_ranked.txt` of a solve2D folder, shows the averages in a zoomable grid
with population and resolution captions and a display popup (brightness, contrast, a min/max histogram), and
writes a 0/1 `selection.txt` that the user then feeds to `simple_exec prg=selection oritype=cls2D infile=`. It
needs a Python with Tk 8.6, numpy, mrcfile, Pillow and tkinterdnd2, it re-implements the class document and MRC
readers, and the selection reaches the project only through a second program. `select2D` moves it into SIMPLE as a
commander: it makes an execution directory, copies the project there as every project program does, shows the
project's class averages from its `cls2D` and `out` segments, and writes the selection into that local project
copy with the same routine `prg=selection` uses.

## 2. Contract

- `select2D` is a `simple_exec` program (decision 2): `simple_exec prg=select2D projfile=<p>.simple`.
  With the default `mkdir=yes`, `params%new` creates `<n>_select2D/`, enters it and copies the project there
  (`setup_execution_context`, `simple_parameters_phases.f90`); every read and write of the run is on that copy.
  The input project is never written.
- A saved selection has exactly the effect of `prg=selection oritype=cls2D infile=<the same 0/1 file>` with the
  same `append`: `os_cls2D` `state`, then `map2ptcls_state(append)` into `ptcl2D`, `ptcl3D` and `cls3D`, through
  one shared routine (decision 5). The `selection.txt` written beside the project reproduces the selection with
  `prg=selection`.
- `prg=selection` is byte-identical after the refactor of Phase 1 (its project file, at one thread, on the Phase 0
  fixture).
- Off-state: no other program changes behaviour. Without `prg=select2D` no GTK code runs in any executable; the
  only GTK call before the window is the display check of decision 15.
- The display maths of decision 9 follow `select2d.py` number for number; deviations are listed in decision 12.

## 3. Decisions and rulings

1. (2026-11-24, proposed) The program is `select2D` (the 2D naming of `solve2D`/`refine2D`), category Class
   Average Processing (`simple_ui_cavgproc.f90`), display name "Select Class Averages", router
   `simple_exec_cavgproc.f90`. Inputs: `projfile` and `append` (existing key, yes|no, default no). No new
   parameter keys; `oritype` is set to `cls2D` by the commander.
2. (2026-11-24, owner) Executable: `simple_exec` links `SIMPLE_GTK_LIBRARIES`, so every node that runs
   `simple_exec` needs the GTK 3 runtime libraries (not a display); GTK is already a build requirement of every
   build. The alternative is a separate front door (`simple_gui` with its own router, a case in
   `cmdline%parse`, its own `prg=list`): more machinery, GTK confined to one binary. Phase 0 records whether
   `ldd` of a GTK-linked binary resolves on the compute nodes; if not, this decision is reopened (stop condition).
3. (proposed) The commander: set `mkdir=yes` and `append=no` when absent, set `oritype=cls2D`, check the display
   (decision 15), then `build%init_params_and_build_spproj(cline, params)` (one `params%new`, one full project
   read of the copy), build the selection session and the display state, run the window until it closes,
   `simple_end`. About 40 lines; the window module is imported inside the procedure (compile-time policy).
4. (proposed) Writing: nothing is written until the user presses Save selection. Save applies the states
   (decision 5) and writes the four changed segments `cls2D`, `cls3D`, `ptcl2D`, `ptcl3D` with
   `write_segment_inside` to the local copy, and writes `selection.txt` (one 0/1 per `cls2D` row, row order).
   Closing with unsaved changes asks Save / Discard / Cancel. The front door's job-description write after the
   commander is the existing deviation of the layering policy (section 8) and is not extended.
5. (proposed) One routine: `sp_project%report_state2cls2D(states, append)` in `simple_sp_project_cls.f90`, beside
   `report_state2mic`/`report_state2stk`, holding the two statements of `exec_selection`'s cls2D branch
   (`os_cls2D%set_all('state', real(states))`, `map2ptcls_state(append=)`). `exec_selection` calls it;
   `select2D` calls it. Its NICE extras (thumbnails, JPEG montages, `gui_comm`) stay in `exec_selection`.
6. (2026-11-24, owner) Initial selection: the project's current states (a visible class with `state > 0` starts
   selected), so reopening shows the saved selection and a fresh solve2D project opens all selected. The
   script's behaviour (nothing selected) was not chosen.
7. (proposed) Visible classes follow the script: a class is shown unless its `pop` is present and 0. Hidden
   classes count in the totals and are always written 0.
8. (proposed) Data: the `cls2D` rows (row index = class) with `pop` and, when present, `res`; the stack of the
   `out` entry `imgkind=cavg`, read with `read_cavgs_into_imgarr` (`simple_imgarr_utils`). The row count must
   equal the stack's image count, else `THROW_HARD` before the window opens. The script's discovery of
   `cavgs_iter<N>.mrc` and its `classdoc_ranked.txt` reader are not ported.
9. (proposed) Display numbers, from `select2d.py`: brightness -100..100 %, integer, step 1; contrast 0.1..5,
   two decimals, step 0.01; min/max inside the stack's finite range, step (range)/1000, min < max, a slider
   clamps to one ulp short of the other limit; stack range over finite pixels of all classes, a constant stack
   gets upper = min + max(1, 1e-6·|min|); automatic window = medians, taken separately, of each visible class's
   1st and 99.9th percentiles (numpy's linear interpolation) of its finite pixels, classes whose two
   percentiles coincide skipped, fallback the stack range; histogram of 256 equal bins over the stack range of
   all classes (hidden ones included), last bin closed; mapping v = mid + c·(x − mid) + (b/100)·span with
   mid, span of the current min/max, then (v − min)/(max − min) clipped to [0, 1], non-finite pixels 0; the
   real image is resampled to the tile (Lanczos-3 when shrinking, bilinear when enlarging, as Pillow), clipped,
   and only then converted to 8 bits with rint(255·v); the display-range box spans the raw values
   mid + (level − 0.5 − b/100)·span/c for level 0 and 1, clipped to the histogram range. Opening resets
   min/max to the automatic window; Reset restores brightness 0, contrast 1 and the automatic window.
10. (proposed) Grid, from the script: margins and spacing 12 px; tile 240 px at start; Ctrl- or Cmd-scroll zooms
    by 1.15 per wheel step between 50 px and min(1000, width − 24) px; columns =
    (width − 24 + 12) div (tile + 12), at least 1; captions below the square image (rows layout under 180 px
    tiles, two columns from 180 px), population with thousands separators, resolution in three significant
    digits, the font shrinking until "999,999" fits, then numeric abbreviation, then an ellipsis; a selected
    tile gets an accent outline, a tinted caption and a check mark. The grid is one drawing area beside its own
    scrollbar and draws only the visible rows, so the work of a frame is bounded by the window, not by the
    number of classes, and the layout and hit test are the script's canvas arithmetic.
11. (proposed) Interaction, from the script: a click toggles one class; Shift-click adds the inclusive range in
    display order between the previous and the current click (every click becomes the next anchor); sort by
    resolution ascending (the start), population descending or class ascending, stable on the current order;
    right-click, Ctrl-click on macOS, Shift+F10 or the Menu key open a popover at the pointer with Sort by
    (radio), Select All, Deselect All, Brightness and Contrast (slider and entry; Enter or leaving the entry
    commits, an invalid entry reverts), Reset and the histogram (drag the Min/Max lines or handles; Up/Down
    choose, Left/Right move by one step); the footer shows "Selected: n/N classes, p/P particles".
12. (proposed) Not ported: Open folder, the path entry, drag-and-drop and the folder chooser (a session is one
    project, given by `projfile`); the Tk theme colour lookup (fixed colours, light background); the Tk/Quartz
    popup workarounds (a GTK popover dismisses itself and places itself on the right monitor); the 33 ms
    coalescing timer (GTK's frame clock draws at most once per frame; only visible tiles are rendered, cached
    per class for the current display revision and tile size).
13. (proposed) Homes (all GTK code under `src/main/gtk/`, everything else GTK-free and testable headless):
    - `percentile` (numpy linear) added to `simple_stat`, the owner; the two private copies in
      `simple_class_compatibility` are left alone (not in this plan);
    - `src/main/image/simple_cavg_display.f90`: class `cavg_display` (decision 9: range, automatic window,
      histogram, settings with validation, mapping box, `render(icls, tile, gray8)`), holding the raw images
      as single-precision arrays and computing in double precision;
    - `src/main/project/simple_cls2D_selection.f90`: class `cls2D_selection` (decisions 6-8, 11: records,
      visibility, display order, toggle/range/all/none, counts, `states()`, `save(spproj, projfile, append)`
      doing decision 4's writes);
    - `src/main/gtk/simple_select2D_window.f90`: class `select2D_window` (`new/run/kill`), the GTK 3 view;
      its `bind(c)` callbacks are private procedures that recover the object from the user data; layout
      arithmetic (columns, tile rectangles, hit test, content height, zoom clamp) as procedures that need no
      GTK, so they test headless;
    - `src/main/commanders/simple/simple_commanders_select2D.f90`: `commander_select2D` (decision 3).
14. (proposed) Compile time: only files under `src/main/gtk/` import gtk-fortran modules, always with `only:`;
    the commander imports the window at procedure level. If Phase 0 measures the GTK tester's compile above 30 s
    (guessed), the window is split into a parent and submodules by concern (grid, popover, histogram) in
    Phase 4, each importing only what it calls.
15. (proposed) Without a display (`gtk_init_check` false) the commander stops with "select2D needs a desktop
    session: run it where a display is available (locally or over ssh -X)" before `params%new`, so no
    directory is made.
16. (proposed) GTK API floor: GTK 3.22 (Oracle Linux 8). Variadic GTK functions are absent from gtk-fortran: the
    close dialog is `gtk_dialog_new` with `gtk_dialog_add_button`; event data comes through the
    `gdk_event_get_*` accessors, not through `gdk_events` struct layouts.
17. (proposed) Tests: `unit_numerics` sub-suite statistics gains the percentile checks; new sub-suites
    `cavg display` (`unit_image`) and `cls2D selection` (`unit_project`); the platform suite `gtk` gains
    `select2D layout` (headless) and `select2D window` (a fixture project, a window driven by the test's calls,
    saved and closed by itself, states read back from disk; skipped without a display). No new CTest entry.
18. (proposed) Living description: `doc/policies/2D/select2D_policy.md`, written in Phase 6 (what is read and
    written, decisions 4-11 in present tense); the UI help points to it.
19. (2026-11-24, owner) UI visibility: Standard. A NICE job of `select2D` on a server without a display ends in
    decision 15's message; if NICE should not offer it at all, that is a NICE-side filter, not in this plan.

Defaults taken without a question: the file names of decision 13; window title "Select 2D - <project name>";
`selection.txt` as the file name of decision 4; the GApplication id `org.simple.select2D` with
`G_APPLICATION_NON_UNIQUE`; the fixed colours: text #20242b, histogram bars #748ba7, Min #3378bd, Max #c06532 and
range box #8b4a9e from the script, background #f6f6f6 and accent #3584e4 chosen here (the script took both from
the Tk theme).

## 4. Not in this plan

NICE (its GUI, its JSON client and the `gui_comm` metadata of `prg=selection`); the stream; other GUI programs;
GTK 4; any change to `prg=selection` beyond Phase 1's extraction; `oritype` other than `cls2D`; particle-level
selection; the private percentile copies in `simple_class_compatibility`; `select2d.py` itself (it stays where it
is until the owner retires it after Phase 7). Files outside the table are out of bounds.

## 5. Change types and gates

| Phase | Kind | Gate |
| --- | --- | --- |
| 0 | measurement | base build, fast gate, `test=gtk`, two `prg=selection` runs byte-compared, compile time recorded |
| 1 | identical | `prg=selection` project byte-identical to Phase 0; `unit_project` |
| 2 | new code, tests-only gate | `unit_numerics`, `unit_image` |
| 3 | new code, tests-only gate | `unit_project` |
| 4 | new code | `test=gtk` (all sub-suites with a display on the Mac; headless ones on Oracle Linux 8) |
| 5 | new code | `test=gtk`; manual checklist |
| 6 | new program | fast gate (registry, CLI chain, descr, guide imports, submodule linkage), UI JSON validation, the commander run of Phase 6's exit |
| 7 | reported, not gated | the owner's real-data read |

## 6. Acceptance

- Percentile and display mapping: equal to numpy (`numpy.percentile(..., method='linear')`, the formulas of
  decision 9 in float64) within 1e-12 relative (*oracle*), the expected values pasted in the tester with the
  numpy line that made them.
- Resampling: a constant image stays constant within 1e-6 (*closed form*); bilinear enlargement of a linear ramp
  is linear within 1e-6 (*closed form*); the output size is the tile size (*exact*). Agreement with Pillow's
  Lanczos is reported in Phase 7, not gated.
- Selection, counts, written states and `selection.txt`: exact (*integers*). `report_state2cls2D` on the
  fixture equals the two statements it replaces, compared field by field (*exact*).
- `prg=selection` identity: the project file byte-identical to the Phase 0 runs at one thread (*measured*:
  Phase 0 runs it twice and records whether the file is deterministic; if not, the exclusions are declared
  there and frozen).
- Window: GApplication status 0, at least one draw, the tile area drawn at least at its requested size, the
  scripted selection on disk after the scripted Save (*exact*).
- Reported, not gated (Phase 7, by the owner, on a real solve2D project whose path Phase 0 records): the same
  visible classes, sort orders and footer counts as `select2d.py` on the same run; `selection.txt` of an
  identical selection byte-equal to the script's; screenshots of the grid and the popover in
  `doc/implementation_notes/planned/select2D_screens/`.

## 7. Common rules and stop conditions

No new parameter key. No floor, tolerance or assertion relaxed. Fix as you go, with a file table row before the
edit. Tests at one OpenMP thread; the display render has no OpenMP region in this plan. Every phase is built and
run before it is marked done: the GUI phases need the Mac (gfortran 16, a display) and the Oracle Linux 8 box
(gcc-toolset-15, headless, and once over `ssh -X`).

Stop and report when: a needed GTK function is missing from GTK 3.22 or from gtk-fortran 3.24.49; a key outside
`projfile`, `append`, `mkdir`, `oritype` would be needed; `prg=selection` is not identical after Phase 1; a GTK-linked
binary does not start on a compute node (decision 2 reopens); a GTK file's compile exceeds 60 s after the split of
decision 14 (guessed); a file outside the table needs an edit. Gates no ruling may waive: the fast gate, the
`prg=selection` identity, `test=gtk` with a display on the Mac.

## 8. Phases

- **Phase 0, baseline.** Build Release in `build/` (`./compile_clean.sh`) on the Mac and on Oracle Linux 8; fast
  gate; `simple_test_exec test=gtk` on both (window sub-suite only with a display); record
  `pkg-config --modversion gtk+-3.0 glib-2.0` on both; `ldd build/bin/simple_test_exec` on an Oracle compute node
  (decision 2); `scripts/profile_build.sh touch src/main/gtk/simple_gtk_tester.f90` for the cost of a file that
  uses `gtk` (decision 14); make the Phase 1 fixture (a small project: 12 classes with two of `pop` 0, 300
  particles, a cavg stack of box 64) and run `prg=selection oritype=cls2D infile=` on it twice, keep both project
  files and compare them. Record the real solve2D project for Phase 7 (owner). Exit: all recorded in Progress,
  base failures named and left.
- **Phase 1, the shared routine.** `report_state2cls2D` (decision 5); `exec_selection` calls it. Exit: the
  `prg=selection` project byte-identical to Phase 0; `test=unit_project` passes.
- **Phase 2, display maths.** `percentile` in `simple_stat`; `simple_cavg_display` and its tester. Exit:
  `test=unit_numerics suite=statistics` and `test=unit_image suite=cavg_display` pass with the acceptance of
  section 6.
- **Phase 3, the selection session.** `simple_cls2D_selection` and its tester (fixture project in
  `enter_fixture`; visibility, three sort orders and their stability, toggle, range across a re-sort, all/none,
  counts, save with `append=no` and `yes` read back, `selection.txt`). Exit: `test=unit_project
  suite=cls2D_selection` passes; `SIMPLE_UNIT_ORDER=reverse` too.
- **Phase 4, the window.** `simple_select2D_window`: grid, captions, check marks, click and Shift-click, zoom,
  scroll, resize, footer, Save, the close dialog; tester sub-suites `select2D layout` and `select2D window`. Exit:
  `test=gtk` passes on the Mac with a display and its headless sub-suites on Oracle Linux 8.
- **Phase 5, the popover.** Sort radio, Select/Deselect All, brightness and contrast, Reset, histogram with
  markers, keys and the range box; the window sub-suite drives a sort, a brightness change and Reset. Exit:
  `test=gtk` passes; the checklist of Progress walked by hand on the Mac.
- **Phase 6, the program.** `commander_select2D`, UI entry, router case, `simple_exec` link (per decision 2),
  `doc/policies/2D/select2D_policy.md`. Exit: fast gate passes; UI JSON validation passes; on the Phase 0
  fixture `simple_exec prg=select2D projfile=<fixture>` makes `<n>_select2D/` with a copy of the project, the
  window opens, a selection saved there and `simple_exec prg=print_project_field oritype=cls2D` in that folder
  shows it, and the input project is byte-unchanged; without a display it stops with decision 15's message and
  makes no folder.
- **Phase 7, real-data read and report.** The owner's read of section 6 on the recorded project;
  `select2D_report.md`. Exit: the report is written; reported items listed for the owner.

### File table

| File | What | Phases |
| --- | --- | --- |
| `src/main/project/simple_sp_project.f90` | interface of `report_state2cls2D` | 1 |
| `src/main/project/simple_sp_project_cls.f90` | body of `report_state2cls2D` | 1 |
| `src/main/commanders/simple/simple_commanders_project_core.f90` | `exec_selection` calls it | 1 |
| `src/utils/math/simple_stat.f90`, `simple_stat_tester.f90` | `percentile` and its checks | 2 |
| `src/main/image/simple_cavg_display.f90`, `simple_cavg_display_tester.f90` | display state | 2 |
| `src/main/project/simple_cls2D_selection.f90`, `simple_cls2D_selection_tester.f90` | selection session | 3 |
| `src/main/gtk/simple_select2D_window*.f90` | the window (submodules only if decision 14 asks) | 4, 5 |
| `src/main/gtk/simple_gtk_tester.f90` | `select2D layout`, `select2D window` | 4, 5 |
| `src/main/commanders/test/simple_commanders_test_class.f90` | suite rows | 2, 3, 4 |
| `src/main/ui/simple_test/simple_test_ui_class.f90` | `suite=` lists | 2, 3, 4 |
| `src/main/commanders/simple/simple_commanders_select2D.f90` | `commander_select2D` | 6 |
| `src/main/ui/simple/simple_ui_cavgproc.f90` | the `select2D` program | 6 |
| `src/main/exec/simple_exec_cavgproc.f90` | router case | 6 |
| `production/CMakeLists.txt` | GTK link of `simple_exec` (decision 2) | 6 |
| `doc/policies/2D/select2D_policy.md` | living description | 6 |
| `doc/implementation_notes/planned/select2D_report.md`, `select2D_screens/` | report, screenshots | 7 |

## Progress

(one entry per phase, written by the run)
