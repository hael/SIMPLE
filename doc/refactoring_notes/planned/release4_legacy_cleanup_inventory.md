# Release 4: Legacy cleanup, outstanding items

**Status:** outstanding items only, 2026-10-03. Everything implemented is
recorded in
[`../completed/release4_legacy_cleanup_report.md`](../completed/release4_legacy_cleanup_report.md),
including what was deliberately kept. The IDs below are labels for this list
(the letter keeps the original grouping: A dead code, B old-data support, C
design decisions, E checked and kept, V verification). Line numbers are
approximate, from 2026-10-03.

Statuses:

- `open`: ruled, not yet implemented;
- `deferred`: blocked by the scope rulings below;
- `confirm`: implemented, but needs the developer's confirmation or run.

## Rulings that still apply

- **Clean break.** "No backwards compatibility. This is the current standard."
  The one exception is reading `.simple` project files from earlier releases
  (narrower particle records).
- **GUI and stream code is out of scope** while Joe works in that area:
  `src/main/stream/`, `production/simple_stream.f90`, `src/utils/gui/` and
  `nice/`.
- **GPU offloading is work in progress and stays as it is.**

## 1. Before committing this change set

| ID | What | Status |
| --- | --- | --- |
| V1 | Final build, fast gate, and the project, sigma, Cartesian and UI/NICE tests on the developer's machine. | open |
| V2 | **Confirm the split of the particle record** (raised by the independent review: the 64-slot record also cost memory). In memory a particle holds the 53 named slots (`N_PTCL_ORIPARAMS = 53`, slot 42 spare); on disk the record is `N_PTCL_RECORD_REALS = 64` reals with zero padding (`src/defs/simple_defs_ori.f90`, `src/fileio/simple_binoris.f90`). This revises the original decision of 64 slots in memory as well, and saves about 420 MiB per particle segment at ten million particles. Reverting to 64 in memory is a two-line change. | confirm |
| V3 | Build and link on the Dell (Oracle Linux) without librt (no longer required by the build). | confirm (after the push) |
| V4 | When v4.0.0 is tagged, update `README.md` (version badge, release link, tarball name) and `doc/installation.md` (tarball name); they still name v3.0.0 (the CMake project version is already 4.0.0). | open (at tagging) |

## 2. Planned

None. The last planned item, C6 (remove the finished-halfmap trailing blend),
was implemented and validated on 2026-10-04; see
[`../completed/trailing_reconstruction_without_halfmap_blend_report.md`](../completed/trailing_reconstruction_without_halfmap_blend_report.md).

## 3. Deferred: GUI and stream (Joe's area)

| ID | What | Where | Size | Notes | Status |
|---|---|---|---|---|---|
| A1 | Old NICE message channel `simple_nice`. It sends `{"jobid":..,"job":..}` without `version`. `nice_lite/api.py:189` answers every such message with HTTP 400, so it does nothing. | `src/utils/gui/simple_nice.f90`. Still instantiated in the `simple_commanders_` modules `refine2D`, `project_core`, `project_cls`, `project_ptcl`, `starproject` and `imgops`, in `single_commanders_nano2D`, and in the `commanders`/`stream` APIs | ~1,750 + 9 call sites | `simple_gui_communicator` is the live channel. Before deleting, check each commander for reporting that only goes through `simple_nice` (that reporting fails today anyway). | deferred |
| A2 | Stream subroutines with no callers. | `progress_estimate_preprocess_stream` (re-exported by `simple_stream_api`). The stream utilities, chunk utilities and repick routines listed here before were removed on 2026-10-05 (WS9 of `stream_fix_plan_2026-10-05.md`); the p03 routines went with the stage refactor. | small | Several of these modules can then come off the `-Wno-unused-function` list (`src/CMakeLists.txt`). | deferred |
| A8 | The stream aliases `initial_analysis`/`opening2D` in `simple_stream.f90` ("to maintain GUI support"; NICE never launches them). The non-stream aliases are done. | `production/simple_stream.f90:71` | ~3 | | deferred |
| A13 | Unused constants. | `simple_defs_stream.f90`: `CHUNK_CLS_REJECTED, REJECTED_CLS_STACK, SIEVING_REFS_FNAME, FLUSH_TIMELIMIT, NMICS_DELTA, POOL_FREQ_REJECTION, SIEVING_MATCH_CAVGS_MAX, SIEVING_REF_CAVGS_MAX, STREAM_NMOVS_SET_TIFF, STREAM_SRCHLIM, FRAC_SKIP_REJECTION`. `simple_defs_fname.f90`: `PAUSE_STREAM` (`STREAM_DESELECTED_REFS` was removed on 2026-10-05) | ~15 | `POLAR_REFS_FBODY` is already removed. | deferred |
| A14 | Commented-out blocks and the dead routine behind one of them. | `p03_initial_analysis.f90:484-495`, which also prints the stale `SIMPLE_GEN_PICKREFS NORMAL STOP` at 497 and 856. `p04_refpick_extract_new.f90:347-349, 448-456`, which makes `validate_ptcl2D_star_inputs` (465-505) dead. `p06_pool2D_new.f90:224-232`, `simple_stream_utils.f90:566-593`, `pool2D_utils.f90:1004,1007`, `refine2D_utils.f90:92`. | ~120 | Keep the intent note in `p00_master.f90:1443-1447` if you like it. The dead `validate_ptcl2D_star_inputs` also holds the last copy of the stack-index tolerance removed elsewhere in release 4 (it accepts a particle without `indstk`), so removing it completes the removal of that tolerance. | deferred |
| A19 | NICE leftovers. `WorkspaceModel.nstr` ("backward-compatibility payload"; never read or written). `Project.new` has a two-signature compatibility form that only tests use. `index.html:138-146` migrates an old session-storage key. `batch_views.py:708` falls back to `latest_micrographs`, a stream key. Unused images: `simple_stream.png`, `simple_stream_min.png`, `single_logo.png`. | `nice/nice_lite/` | ~40 + a migration | `nstr` is best removed together with C12. | deferred |
| B13 | **File-based stream parameter updates** (`stream_user_params.txt`). Nothing in SIMPLE or NICE writes it any more; NICE pushes updates over HTTP. | `update_user_params` and `USER_PARAMS`, removed on 2026-10-05 (WS9 of `stream_fix_plan_2026-10-05.md`) | ~95 | Only someone hand-writing that file is affected. | confirm |
| C12 | **NICE migrations 0001-0007**: squash them into a single fresh `0001_initial`. Existing installs then need a fresh database or `migrate --fake-initial`. Also remove the build-time `makemigrations` (`nice/CMakeLists.txt:111`), which can silently create migrations on user installs. | `nice/nice_lite/migrations/` | 7 files to 1 | 0006 and 0007 (program renames) were added during this cleanup. | deferred |
| C13 | **Old GUI stats channel** (`.guistats`/`.poolstats`): nothing in SIMPLE or NICE reads `.poolstats`. Keep `generate_pool_jpeg` (`pool_jpeg_map` is still used). | `src/utils/gui/simple_guistats.f90`. The `.poolstats` writer and `cls2D_thumbnail.jpeg` were removed on 2026-10-05 (WS9 of `stream_fix_plan_2026-10-05.md`); `generate_pool_stats` still uses `guistats%generate_2D_jpeg` as the sprite-sheet fallback | ~600 | | deferred |
| C14 | **Stream parameter shims marked "backwards compatibility"**: `ncls_start`, `nparts_pool`, `nthr2D`, `nparts_chunk`. These are a refactor rather than a deletion; the chunk path still reads them. | `commanders_solve2D.f90:816,863,1120,1121,1197` (the `solve2D_chunks` program), `solve2D_chunk.f90:94,96,266,311`, `stream_pool2D.f90:210,1147` | small, but cascades into params and UI | | deferred |
| C16 | **`*_new` suffixes on stream modules p01/p02/p04/p05/p06.** The old versions are gone, so this is a rename only. | 5 files, ~14 `use` lines, CMake | rename | | deferred |
| C19s | `moldiam_max` "kept for API compatibility". | `simple_mini_stream_utils.f90:257,624` | small | | deferred |
| C20s | The `compile_gui.sh` TMPDIR workaround; `nice/.../updatesubmissiontemplate.py`; `LD_LIBRARY_PATH` in `wsgi.py.in`. | `compile_gui.sh`, `nice/` | small | | deferred |
| E1 | **The `opening2D` naming in GUI keys and NICE stats** is a coordinated SIMPLE plus NICE rename. Optional; leave it out of release 4 unless wanted. | stream and NICE | small | | deferred |

## 4. Deferred: GPU offloading

| ID | What | Status |
| --- | --- | --- |
| C20g | `compile_gpu_debug.sh` is referenced nowhere; revisit with the GPU work. | deferred |
