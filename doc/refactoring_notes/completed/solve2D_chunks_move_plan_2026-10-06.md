# solve2D_chunks move plan, 6 October 2026

`solve2D_chunks` splits a project into stack-bound, particle-balanced subsets and runs one
`solve2D` per subset. It is not a stream program: the refine2D exec router runs it
(`prg=solve2D_chunks`), yet it lives in `src/main/stream` as `stream_solve2D_chunks`, reaches
everything through the `simple_stream_api` umbrella, and uses a chunk layer in
`src/main/stream/pool2D` that nothing in the stream uses any more. This plan moves the commander
beside `commander_solve2D` and the chunk layer beside the solve2D controller. It resolves C15 of
`release4_legacy_cleanup_inventory.md` as "kept, moved".

## Status

Steps 1–6 are implemented (6 October 2026), not yet compiled or run (§3).
- C15 has left the inventory for §12 of `release4_legacy_cleanup_report.md`; C14 names the new
  files.
- `pool2D_policy.md` named `setup_downscaling` for the pool's downscaling, which is
  `stream_pool2D%set_dimensions` since the pool became a type; it now names that.
- `simple_solve2D_chunk` keeps `simple_core_module_api`, and its other imports are explicit.
  `copy` allocates the project and command line before their defined assignments.

## Decisions

1. **The commander's home.** `commander_solve2D_chunks` (`execute => exec_solve2D_chunks`), in
   `simple_commanders_solve2D` beside `commander_solve2D`. `simple_exec_refine2D` keeps routing
   `prg=solve2D_chunks` and imports it from there. Not in the exec router itself: the router only
   routes (ui -> exec -> commander). The chunk layer's imports sit inside the procedures that use
   them, so the module's importers (the exec router, `simple_commanders_validate`, the high-level
   tests) gain no build-order edge.
2. **The chunk layer.** `simple_stream_chunk` moves to `src/main/solve` as `simple_solve2D_chunk`,
   holding the type `solve2D_chunk` (named after its module). `init_chunk_clustering` and the
   downscaling it calls become private helpers of the commander module, and
   `simple_stream_chunk2D_utils` goes. Nothing stream-named stays outside the stream.
3. **Imports.** Explicit `use ..., only:` lists (compile-time policy); `simple_core_module_api`
   stays, as for the pool modules. The moved code no longer uses `simple_stream_api`.

## 1. Where things stand

| File | What | Used by |
|---|---|---|
| `src/main/stream/simple_stream_solve2D_chunks.f90` | `stream_solve2D_chunks` (`execute => exec_stream_solve2D_chunks`) with its internal steps (subset projects, chunk start-up, the submission loop) | `simple_exec_refine2D` |
| `src/main/stream/pool2D/simple_stream_chunk.f90` | `stream_chunk`: one subset's project, queue environment, command line and job | the program, `simple_stream_api` (re-export) |
| `src/main/stream/pool2D/simple_stream_chunk2D_utils.f90` | `init_chunk_clustering` (the chunks' `solve2D` command line, dimensions, low-pass limits) | the program, `simple_stream_api` (re-export) |
| `src/main/stream/pool2D/simple_stream_refine2D_utils.f90` | `setup_downscaling`, used only by `init_chunk_clustering` | the chunk layer, `simple_stream_api` (re-export) |

- The program also uses stream constants for its clean-up (`POOL_DIR`, `POOL_PROJFILE`,
  `REFINE2D_FINISHED`, `STDERROUT_DIR`, `DIR_SNAPSHOT`). They stay in the definition modules, which
  `simple_core_module_api` provides.
- `stream_chunk` holds `sp_project`, `qsys_env` and `cmdline` as plain components, which the
  compile-time policy forbids.
- No test covers the program or the chunk layer.

## 2. Steps

1. **The chunk type.**
   - `git mv` `src/main/stream/pool2D/simple_stream_chunk.f90` to
     `src/main/solve/simple_solve2D_chunk.f90`, so the history follows.
   - Rename the module `simple_solve2D_chunk` and the type `solve2D_chunk`; update the header.
   - Make its heavy components (project, queue environment, command line) allocatable, allocated
     in `init_chunk` and released in `kill`.
   - Replace any umbrella import with `only:` lists.
2. **The commander,** in `simple_commanders_solve2D`:
   - `type, extends(commander_base) :: commander_solve2D_chunks` with
     `procedure :: execute => exec_solve2D_chunks`, public beside `commander_solve2D`;
   - `exec_solve2D_chunks` is today's `exec_stream_solve2D_chunks` with its internal procedures,
     its imports (`solve2D_chunk`, `rec_list`, `project_rec`, `qsys_cleanup`, ...) inside the
     procedure;
   - `init_chunk_clustering` and `setup_chunk_downscaling` (today's `setup_downscaling`) become
     private module procedures, each importing what it uses.
3. **The exec router:** `use simple_commanders_solve2D, only: commander_solve2D,
   commander_solve2D_chunks`; `type(commander_solve2D_chunks) :: xsolve2D_chunks`.
4. **Removals:**
   - `simple_stream_solve2D_chunks.f90` and `simple_stream_chunk2D_utils.f90`;
   - `setup_downscaling` from `simple_stream_refine2D_utils`;
   - `stream_chunk`, `init_chunk_clustering` and `setup_downscaling` from `simple_stream_api`
     (its other users, `simple_ptcl_sieve_utils` and the high-level tests, use none of them);
   - the comment in `simple_stream_pool2D` that names `setup_downscaling`.
5. **Unchanged:** the program name and its UI (`simple_ui_refine2D`), NICE, CMake (it globs the
   folders), the stream, and the C14 shims the chunk path still reads (`ncls_start`,
   `nparts_chunk`, `nthr2D`).
6. **Docs:**
   - `stream_refactor.md` (the folder lists);
   - `pool2D.inf` (no chunk layer; and it still says "module state");
   - `solve.inf` (the chunks);
   - `release4_legacy_cleanup_inventory.md` (C15 kept and moved; C14's file references);
   - the skills, with approval: `simple-main-stream` (lines on the chunk layer and on
     `simple_stream_solve2D_chunks`) and `simple-solve2d` (the program and its chunk type).
   - The code overview indexes are regenerated by the build.

## 3. Checks

- A build without warnings in the touched files.
- `scripts/check_test_registry.py`, which does not change.
- One run of `simple_exec prg=solve2D_chunks` on a small project. No test covers it, so the run
  is the check that the subsets, the chunk jobs and their outputs (`chunk_<n>/chunk.simple`) come
  out as before.

## 4. Behaviour

None changes: the same subsets, command lines, jobs and outputs. The allocatable components of
`solve2D_chunk` and the procedure-level imports are structural only.
