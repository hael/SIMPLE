# Interactive tracking-position picker

`simple_viewer prg=pick_positions stk=<MRC>` replaces the removed Python
manual picker. Box size is chosen entirely in the window: the slider and
numeric field start at 64 pixels (or the smaller image dimension), resize
all picks, preserve centres where possible, and move boxes inward at edges.
`box=` is not a program option and is explicitly rejected with slider guidance.

The middle frame (`nframes/2+1`, matching Python's zero-based middle) is shown
with `fromf=0`; a positive `fromf` selects a one-based frame. Left click adds
or selects, left drag moves, right drag pans, the wheel zooms about the
pointer, and Delete/Backspace removes the selected box. Contrast starts at
the finite-pixel 1st/99th percentiles.

Save and normal window close write `outfile=positions_all.box` as four integer
columns: zero-based lower-left x/y, width, height. Pick order and Python's
ties-to-even rounding are preserved. Save errors keep the picks and window.

The GTK-free manual_picker owns selection, resizing, clamping, viewport
transforms and box output. The GTK window reuses stack_display, the display
panel and Cairo helpers. The commander loads one frame and uses a GTK-owned
initial-size constant, independent of parameters%box. Shared parameter phases
and their tests match HEAD; there are no viewer exceptions to grid validation.
UI/default generation covers viewer commanders. No new CLI key, executable
or CTest entry was added.

Validation of the implementation:

- Incremental Linux Release build/install passed without compiler warnings.
- Previous full fast gate: 15/15 passed in 23.17 s with local sockets enabled.
- Manual-picking suite: 28 assertions; reverse unit_image: 16 suites and
  1031 assertions; offscreen GTK Cairo: 9 assertions passed.
- Synthetic 160x120x3 MRC opened at frame 2. Lower-origin orientation and
  controls were visually inspected. User tried the window; normal close
  saved three in-bounds boxes all resized to 44 pixels. Automated mouse
  testing stopped when the user confirmed they were trying the window;
  no full automated GUI gesture pass is claimed.
- Model assertions cover geometry, resize, selection/deletion, zoom/pan,
  output rounding, empty output and failed saves. Viewer listing/help and
  invalid-frame/no-desktop error paths were checked.
- Slider-only box sizing: build/install passed, all 15 fast suites passed
  in 23.04 s, help has no box option, and explicit box= is rejected before GTK
  initialization. Logs: `/tmp/simple-pick-positions-slider-only-*.log`.

Logs and synthetic artifacts: `/tmp/simple-pick-positions-*`. The Python
source and its matching installed copy were removed. Changes are uncommitted.

```sh
build/bin/simple_viewer prg=pick_positions stk=frames.mrc
```

Use `fromf=10` for frame 10 and `outfile=seeds.box` for a different destination.
Feed the saved file to `tseries_track_particles boxfile=seeds.box` with the
appropriate tracking start frame. The viewer does not launch tracking.
