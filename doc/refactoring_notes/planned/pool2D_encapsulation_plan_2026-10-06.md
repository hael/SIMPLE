# Pool 2D encapsulation plan, 6 October 2026

The 2D pool of stream p06 is module state: about thirty-five public variables in
`simple_stream2D_state` and twenty more inside `simple_stream_pool2D_utils`, read and written by
three modules. One process can hold one pool, `kill` resets none of it, and the pool cannot be
unit-tested. This plan makes the pool a type, `stream_pool2D`, following the encapsulated-class rule
of the `simple-modern-fortran` skill. It closes `pool2D_policy.md` §12 "Module state" (review G1,
proposal R6) and step B of `stream_refactor.md`.

## Status

Step 1 is implemented (6 October 2026), not yet compiled:
- `simple_stream_pool2D` holds `stream_pool2D`; `simple_stream2D_state` is deleted;
  `solve2D_chunks` owns its chunk state; p06 holds the pool.
- `simple_stream_pool2D_tester` runs within the "pool 2D" sub-suite (called from p06's tester),
  so the test UI's sub-suite list is unchanged; step 2 can register it on its own.
- Where the code differs from §3:
  - p06's `apply_final_mskdiam` takes the pool's iteration as an argument, so its test sets the
    iteration instead of the pool's state;
  - the snapshot test asks an empty pool for its current iteration (0);
  - p06's `cline` stays, since it starts the pool;
  - `get_mask` is a query added for the mask test.

Step 2 is implemented (6 October 2026), not yet compiled:
- `stats()` returns `stream_pool2D_stats`, which p06's `send_status` and `send_pool_cavgs` read;
  the single getters are gone.
- The append is the pool's; `project_ptr` and p06's `target` are gone.
- The public procedures are those p06 calls, and `init_state` for the tester.
- "pool 2D object" is a registered sub-suite of `unit_stream`.
- Where the code differs from §4:
  - the append takes the sets of one import together (`append_sets(sets, nmics, nsel)`), so the
    pool's rows are reallocated once per import as before, and it returns the pool's micrographs
    and selected particles;
  - the rows are written by a stateless helper (`append_project_sets` in
    `simple_stream_refine2D_utils`), which the sub-suite checks on a plain project;
  - p06 takes the data's box and pixel size from the first set with stacks (they were read from
    the pool's project, whose first stack is that set's);
  - `draw_new_classes` moved to `simple_stream_refine2D_utils` (stateless), and `stats` carries
    the mask, which replaces step 1's `get_mask`;
  - the history ring and `stats` after an iteration need a queue, so they have no unit test.

## Decisions

1. **Shape.** One module, `simple_stream_pool2D`, holding the type `stream_pool2D` with private
   components. Every pool procedure that reads pool state becomes type-bound: those of
   `simple_stream_pool2D_utils`, and the snapshot, publication and final-project procedures of
   `simple_stream_refine2D_utils`. Stateless helpers stay free in `simple_stream_refine2D_utils`;
   no submodules (the pool modules are leaves, so the compile-time policy gives no reason for them).
2. **Chunk state.** The chunk-layer variables, used only by the `solve2D_chunks` program, become
   that program's own. `simple_stream2D_state` is deleted, and `simple_stream_api` stops
   re-exporting it.
3. **The stage's command line.** The pool copies what it reads from p06's command line when it is
   made; the `master_cline` pointer goes.
4. **Two steps.** Step 1 moves the state into the type with behaviour unchanged; step 2 tidies
   the public surface and adds the pool's own tests. Each step builds and passes the tests on its
   own.
5. **Public surface.** The GUI getters fold into one `pool%stats()` returning a small derived
   type; the rest of the public surface is the calls p06 makes.
6. **Tests.** A pool sub-suite (`simple_stream_pool2D_tester`, `unit_stream`); the p06 tester
   stops touching pool internals.
7. **The pool owns its rows.** In step 2 the append of imported sets moves into the pool
   (`pool%append_set`) and p06 no longer holds a pointer to the pool's project.

## 1. Where things stand

| Module | State | Used by |
|---|---|---|
| `simple_stream2D_state` | pool: `cline_refine2D_pool`, `pool_proj` (target), `pool_proj_history`(5) and `pool_history_iter`, `starproj`, `pool_dims`, `pool_native_box`, `pool_native_smpd`, `pool_mskdiam`, `pool_user_lpstop`, `pool_lpstop`, `l_pool_available`, `pool_iter`, `last_complete_iter`, `ncls_glob`, `lpstart`, `lpstop`, `lpcen`, `pool_jpeg_map`/`_pop`/`_res`, `refs_glob`, `orig_projfile`; chunks: `chunks`, `cline_refine2D_chunk`, `chunk_dims`, `glob_chunk_id`, `nptcls_per_chunk`; both: `master_cline`, `l_scaling`, `l_stream2D_active`, `numlen` | `pool2D_utils` (through the `simple_stream_api` umbrella), `refine2D_utils`, `chunk2D_utils`, p06 and its tester, `solve2D_chunks` |
| `simple_stream_pool2D_utils` | its own: `pool_qenv`, `pool_job`, `pool_nattempts`, `l_pool_failed`, `pool_center`, `pool_stacks_mask`, `lim_ufrac_nptcls`, the counts (`nptcls_glob`, `nptcls_rejected_glob`, `ncls_rejected_glob`), `current_resolution`, `resolutions`, `current_jpeg*`, `conv_*` | p06 and its tester |
| `simple_stream_refine2D_utils` | `starproj_stream` | p06, `chunk2D_utils` (`setup_downscaling`) |

- **p06's calls.** About twenty-five kinds of call into the two utils modules. Beyond those, p06
  reads `last_complete_iter` and the JPEG maps straight from the state module, takes a raw
  pointer to `pool_proj` (`get_pool_ptr`) to append imported sets, and nullifies
  `master_cline` in `kill`.
- **p06's tester** writes `pool_iter`, `pool_dims`, `pool_mskdiam` and `cline_refine2D_pool`
  directly, and calls `draw_new_classes`.
- **What the pool reads from p06's command line:** whether `nsample_max`, `update_frac` and
  `cenlp` were given. The chunk layer reads whether `lp`, `lpstart`, `lpstop` and `cenlp` were
  given from the `solve2D_chunks` command line.
- **The chunk layer's use of `pool_proj`:** a placeholder for the computing environment its
  chunks copy. The pool's only use of chunk state is copying `pool_dims` into `chunk_dims`, a
  leftover.
- **Procedures in `refine2D_utils`:**
  - read pool state: `write_project_stream2D`, `terminate_stream2D`, `write_pool_snapshot`,
    `publish_pool_state`, `rescale_cavgs`, `rank_cavgs`, `apply_snapshot_selection`;
  - stateless: `cleanup_root_folder`, `tidy_2Dstream_iter`, `build_pool_publication`,
    `delete_pool_publication`, `pool_publication_names`, `snapshot_cavgs_meta`, `log_rss`,
    `debug_print`;
  - shared with the chunk layer: `setup_downscaling`, which only sets `l_scaling`.

## 2. The type

- **`stream_pool2D`** in `simple_stream_pool2D`:
  - all components private;
  - the heavy ones allocatable, allocated in `start` and released in `kill` (compile-time
    policy): the refine2D command line, the pool project, the history ring
    (`type(sp_project), allocatable :: history(:)`), the queue environment, the STAR writers;
  - `POOL_NHISTORY` and the history slot function become the module's private parameter and
    function.
- **Lifecycle:**
  - p06 allocates an empty pool when it starts. `pool%start(...)` is what
    `init_pool_clustering` does today, at the first import; `kill` releases and resets
    everything.
  - Queries made before the first import (iteration, available, failed) answer 0 or `.false.`,
    as today.
  - `start` is split into the state part (dimensions, mask, low-pass limits, the refine2D
    command line) and the queue part (queue environment, job), so tests can build a pool
    without a queue.
- **The stage's command line:** `start` takes the three "was it given" flags (`nsample_max`,
  `update_frac`, `cenlp`) as arguments, and the pool keeps them as logical components.
- **Procedure map** (step 1 keeps every meaning):

| Today | Type-bound |
|---|---|
| `init_pool_clustering` | `start` (state part) and its queue part |
| `iterate_pool`, `update_pool`, `update_pool_aln_params`, `update_pool_status` | `iterate`, `update`, `update_aln_params`, `update_status` |
| `cancel_pool_job`, `is_pool_failed`, `is_pool_available`, `get_pool_iter` | `cancel`, `failed`, `available`, `iteration` |
| `update_mskdiam`, `set_pool_resolution_limits` | `set_mskdiam`, private `set_resolution_limits` |
| `generate_pool_stats` | `write_stats` |
| GUI getters (`get_pool_assigned`, `_rejected`, `_resolution`, `_cavgs_jpeg`, `_ntilesx`/`_ntilesy`, `_cavgs_mrc`), `last_complete_iter`, the JPEG maps | step 1: one getter each; step 2: `stats()` |
| `get_pool_ptr` | step 1: `project_ptr`; step 2: gone, replaced by `append_set` (§4) |
| `draw_new_classes` | `procedure, nopass` (stateless; the tester calls it) |
| `write_project_stream2D`, `terminate_stream2D` | `write_project`, `finalise` |
| `write_pool_snapshot`, `publish_pool_state` | `write_snapshot`, `publish` |
| `rescale_cavgs`, `rank_cavgs`, `apply_snapshot_selection` | private type-bound |
| `setup_downscaling` | stays free, returns the scaling flag instead of setting module state |

## 3. Step 1: the state moves, behaviour unchanged

1. **The pool module:**
   - `git mv` `simple_stream_pool2D_utils.f90` to `simple_stream_pool2D.f90`, so the history
     follows;
   - declare the type and turn the procedures into type-bound ones (`self%` for module state);
   - move in the stateful `refine2D_utils` procedures and `starproj_stream`;
   - replace the `simple_stream_api` umbrella with `use ..., only:` imports.
2. **`solve2D_chunks`:**
   - its commander owns `chunks`, its refine2D command line, `chunk_dims`, the chunk id counter,
     `nptcls_per_chunk`, `numlen`, its scaling flag, its low-pass limits, and a template project
     for the computing environment (the placeholder `pool_proj` was);
   - `init_chunk_clustering` takes them as arguments and reads the "was it given" flags from
     the commander's own command line;
   - the "called once" guard (`l_stream2D_active`) goes, since the state is local.
3. **Deletions:** remove `simple_stream2D_state`; drop it from `simple_stream_api`; drop the
   pool's copy into `chunk_dims`.
4. **p06:**
   - `type(stream_pool2D), allocatable, target :: pool` (target, because of `project_ptr` until
     step 2), allocated in `init_params` and killed in `kill`;
   - every call becomes `self%pool%...`;
   - the `master_cline` handling goes, and so does the stage's `cline` pointer if nothing else
     reads it.
5. **Tests:** the p06 tests that set pool internals (the mask clamp, the iteration-10 case) move
   to a first `simple_stream_pool2D_tester` that builds a pool with the state part of `start`.
   Everything else in the p06 tester keeps passing unchanged.
6. **Docs:**
   - `pool2D_policy.md` §1 (files) and §12 (the gap closed);
   - `stream_refactor.md` (step B done);
   - the `simple-main-stream` skill's line on module state (a skill edit, to be approved).

**Checks:** a build without warnings in the touched files; `unit_stream` (pool 2D, and the new
pool sub-suite); the `lib_stream` sieve-to-3D chain. No test covers `solve2D_chunks`, so step 1
needs one run of it on a small project (or a test added for it).

## 4. Step 2: the surface and the tests

- **`stats()`** returns a public `stream_pool2D_stats`:
  - fields: the last complete iteration, the assigned and rejected counts, the resolution, the
    class-average JPEG and MRC paths with their tile counts, and the JPEG maps (index,
    population, resolution);
  - p06's `send_status` and `send_pool_cavgs` read it;
  - the single getters go.
- **The project pointer goes; the append moves into the pool** (decision 7):
  - `pool%append_set(set)` appends one sieve set to the pool's project: its micrographs, its
    stacks renumbered after the pool's, and its particles as new rows with no 2D parameters but
    their shifts and the stack index rewritten. It returns the counts p06 keeps: micrographs,
    particles, and selected particles.
  - p06's `transfer_sets` keeps what is p06's business: which sets are taken (hand-off order,
    stepwise import, the sieve's final set and its retraction), the counts and the log. It
    calls `append_set` for each set it takes.
  - `project_ptr` and p06's `target` attribute go; nothing outside the pool writes its project.
  - The p06 tests of the transfer (stacks renumbered, particles new with their shifts, each set
    once, a later set appended) keep passing; the row-level checks move to the pool sub-suite
    as tests of `append_set`.
- **Visibility:** every type-bound procedure p06 does not call becomes private.
- **The pool sub-suite** (`unit_stream`, registered):
  - a pool built without a queue;
  - `set_mskdiam` clamping to the box and moving the low-pass ramp;
  - the resolution limits;
  - the history ring (slots, overwriting, the iterations kept);
  - `stats()` before and after an iteration's results are taken in;
  - a publication built from a synthetic pool project;
  - `kill` resetting everything, and a second `start` after it.
- **The p06 tester** goes through the pool's public procedures only.
- **Docs:** `pool2D_policy.md` names the type's surface and the tests (§11 change rules).

## 5. Notes and risks

- **Still one pool per folder:** the pool's files have fixed names (the pool folder, its exit
  status, its iteration files). The type lifts the one-pool-per-process limit, not this one.
- **Order:** queries before `start` must keep answering as today (iteration 0, not available),
  since p06 asks before the first import.
- **Restart:** behaviour is unchanged: a restarted p06 starts a fresh pool, now a fresh object.
- **Compile time:** the pool module imports only what it uses, and p06 no longer pulls the state
  module. Any compile-time claim needs a `scripts/profile_build.sh` comparison; none is planned.
