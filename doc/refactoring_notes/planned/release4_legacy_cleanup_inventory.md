# Release 4: Legacy cleanup, outstanding items

**Status:** outstanding items only, 2026-10-03; every entry re-verified
against the code on 14 October 2026 (below: two entries are done, one is
reduced to a leftover, line numbers refreshed). 15 entries remain: 2 `open`
(V1, V4), 3 `confirm` (V2, V3, B13) and 10 `deferred` (A2, A8, A13, A19, C12,
C14, C19s, C20s, E1, C20g). Everything implemented is recorded in
[`../completed/release4_legacy_cleanup_report.md`](../completed/release4_legacy_cleanup_report.md),
including what was deliberately kept. The IDs below are labels for this list
(the letter keeps the original grouping: A dead code, B old-data support, C
design decisions, E checked and kept, V verification). Line numbers are
approximate, from 2026-10-14.

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
| V1 | Final build, fast gate, and the project, sigma, Cartesian and UI/NICE tests on the developer's machine. (Not verifiable from the code; nothing was built or run for the 2026-10-14 re-verification.) | open |
| V2 | **Confirm the split of the particle record** (raised by the independent review: the 64-slot record also cost memory). In memory a particle holds the 53 named slots (`N_PTCL_ORIPARAMS = 53`, slot 42 spare); on disk the record is `N_PTCL_RECORD_REALS = 64` reals with zero padding (`src/defs/simple_defs_ori.f90`, `src/fileio/simple_binoris.f90`). This revises the original decision of 64 slots in memory as well, and saves about 420 MiB per particle segment at ten million particles. Reverting to 64 in memory is a two-line change. Re-verified 2026-10-14: implemented as described (`N_PTCL_ORIPARAMS = 53`, `N_PTCL_RECORD_REALS = 64` at `simple_defs_ori.f90:67-68`; `simple_binoris_tester.f90` asserts 53, 64 and the 256-byte record); the decision is the developer's. | confirm |
| V3 | Build and link on the Dell (Oracle Linux) without librt (no longer required by the build). No `librt`/`-lrt` is named in any `CMakeLists.txt` or `*.cmake` file in the tree (2026-10-14); the link itself can only be confirmed on the Dell. | confirm (after the push) |
| V4 | When v4.0.0 is tagged, update `README.md` (version badge, release link, tarball name) and `doc/installation.md` (tarball name); they still name v3.0.0 (the CMake project version is already 4.0.0). Re-verified 2026-10-14: `README.md` (badge line 12, release link 50-52, tarball 96-98) and `doc/installation.md` (27-29) say 3.0.0; `CMakeLists.txt:13` says 4.0.0. | open (at tagging) |

## 2. Planned

None. The last planned item, C6 (remove the finished-halfmap trailing blend),
was implemented and validated on 2026-10-04; see
[`../completed/trailing_reconstruction_without_halfmap_blend_report.md`](../completed/trailing_reconstruction_without_halfmap_blend_report.md).

## 3. Deferred: GUI and stream (Joe's area)

| ID | What | Where | Size | Notes | Status |
|---|---|---|---|---|---|
| A2 | Stream subroutines with no callers. | `progress_estimate_preprocess_stream` (`simple_progress.f90:~89`; imported only in `simple_stream_api.f90:~13`; no callers, 2026-10-14). The stream utilities, chunk utilities and repick routines listed here before were removed on 2026-10-05 (WS9 of `stream_fix_plan_2026-10-05.md`); the p03 routines went with the stage refactor. | small | Several of these modules can then come off the `-Wno-unused-function` list (`src/CMakeLists.txt`). | deferred |
| A8 | The stream aliases `initial_analysis`/`opening2D` in `simple_stream.f90` ("to maintain GUI support"; NICE never launches them). The non-stream aliases are done. | `production/simple_stream.f90:~68-69` (the `case( 'gen_pickrefs','initial_analysis', 'opening2D' )`) | ~3 | | deferred |
| A13 | Unused constants. | `simple_defs_stream.f90`: `CHUNK_CLS_REJECTED, REJECTED_CLS_STACK, SIEVING_REFS_FNAME, FLUSH_TIMELIMIT, NMICS_DELTA, POOL_FREQ_REJECTION, SIEVING_MATCH_CAVGS_MAX, SIEVING_REF_CAVGS_MAX, STREAM_NMOVS_SET_TIFF, STREAM_SRCHLIM, FRAC_SKIP_REJECTION`. `simple_defs_fname.f90`: `PAUSE_STREAM` (`STREAM_DESELECTED_REFS` was removed on 2026-10-05) | ~15 | `POLAR_REFS_FBODY` is already removed. Re-verified 2026-10-14: each of the 12 names still occurs once in `src`, `production` and `nice`, its declaration. | deferred |
| A19 | NICE leftovers. `WorkspaceModel.nstr` (`models/workspacemodel.py:20`, "backward-compatibility payload"; never read or written). `Project.new` has a two-signature compatibility form that only tests use. `index.html:~138-146` (`legacyWorkspaceId`) migrates an old session-storage key. `batch_views.py:~704` falls back to `latest_micrographs`, a stream key. Unused images: `simple_stream.png`, `simple_stream_min.png`, `single_logo.png`. | `nice/nice_lite/` | ~40 + a migration | `nstr` is best removed together with C12. | deferred |
| B13 | **File-based stream parameter updates** (`stream_user_params.txt`). Nothing in SIMPLE or NICE writes it any more; NICE pushes updates over HTTP. | `update_user_params` and `USER_PARAMS` are gone (confirmed 2026-10-14). Left over: the neighbouring `USER_PARAMS2D = 'stream2D_user_params.txt'` at `simple_defs_fname.f90:~125`, still deleted in `cleanup_root_folder` (`stream/pool2D/simple_stream_refine2D_utils.f90:~32`) and touched by `simple_stream_stage_pool2D_tester.f90` (~13, 114, 125). Whether to remove it is the maintainer's call. | ~95 | Only someone hand-writing that file is affected. | confirm |
| C12 | **NICE migrations 0001-0008**: squash them into a single fresh `0001_initial`. Existing installs then need a fresh database or `migrate --fake-initial`. Also remove the build-time `makemigrations` (`nice/CMakeLists.txt:~111`, still runs `manage.py makemigrations`), which can silently create migrations on user installs. | `nice/nice_lite/migrations/` | 8 files to 1 (`0001_initial` plus `0002_dispatchmodel_simple_path` to `0008_snapshot_types`) | 0006 and 0007 (program renames) were added during this cleanup; 0008 (snapshot types) since. | deferred |
| C14 | **Stream parameter shims marked "backwards compatibility"**: `ncls_start`, `nparts_pool`, `nthr2D`, `nparts_chunk`. These are a refactor rather than a deletion; the chunk path still reads them. | `commanders_solve2D.f90:~830` (`nthr2D`), `~1087-1088` (`ncls_start`, `nparts_chunk`) and `~1163` (`nthr2D`) (the `solve2D_chunks` program), `solve2D_chunk.f90:~94,96,266,311`, `stream_pool2D.f90:~212` (`nparts_pool`, read at ~1154). `nthr2D` is also read in `simple_stream_stage_initial_analysis.f90:~713,1342` | small, but cascades into params and UI | | deferred |
| C19s | The unused dummy `moldiam_max`. The "kept for API compatibility" comment and the old line numbers are gone (the file is now ~246 lines), but the dummy argument remains and is never read. | `simple_mini_stream_utils.f90:~19-23` (`segdiampick_preprocess`) and `~77-81` (`segdiampick_mics`) | small | | deferred |
| C20s | The `compile_gui.sh` TMPDIR workaround; `nice/.../updatesubmissiontemplate.py`; `LD_LIBRARY_PATH` in `wsgi.py.in`. | `compile_gui.sh`, `nice/` | small | | deferred |
| E1 | **The `opening2D` naming in GUI keys and NICE stats** is a coordinated SIMPLE plus NICE rename. Optional; leave it out of release 4 unless wanted. | stream and NICE | small | | deferred |

## 4. Deferred: GPU offloading

| ID | What | Status |
| --- | --- | --- |
| C20g | `compile_gpu_debug.sh` is referenced nowhere; revisit with the GPU work. | deferred |

## 5. Done since 2026-10-03 (re-verified 14 October 2026)

Removed from the lists above; nothing further to do unless a residue below
is wanted gone.

- **A14**, the commented-out blocks and the dead routine behind one of them.
  `validate_ptcl2D_star_inputs` and the stale `SIMPLE_GEN_PICKREFS NORMAL STOP`
  text no longer exist anywhere in `src` or `production`. The p03, p04 and p06
  modules the entry pointed at were replaced by the stream refactor, and
  `simple_stream_utils.f90` is now ~146 lines. Residue: one commented-out line,
  `! call simple_rmdir(DIR_SNAPSHOT) ! need snapshots keeping`, at
  `src/main/stream/pool2D/simple_stream_refine2D_utils.f90:~31`. The copy of
  the stack-index tolerance that the entry named went with the routine; note
  for the maintainer that `simple_stream_stage_solve3D.f90` (~660 and ~826)
  takes a row's position in its set as the image index when the row has no
  `indstk`, which is new stage code and not the removed copy.
- **C16**, the `*_new` suffixes on stream modules. No `_new` module or file
  name remains in `src/main/commanders/stream`, `production/` or
  `src/CMakeLists.txt` (the commanders are `simple_commanders_stream_p0N_*.f90`
  without the suffix). Residue outside this inventory's scope: a template
  comment at `nice/nice_lite/templates/jobbuilder.html:128` still names
  `simple_stream_p04_refpick_extract_new.f90`, and
  `doc/refactoring_notes/planned/stream_refactor.md` (~80-86) lists the old
  `_new` module names as its record.
