# Layering and Project Ownership Policy

Status 2026-10-21, drafted for the maintainer's review. Rules tagged **(ruled)**
come from existing policies, skills or maintainer rulings; rules tagged
**(proposed)** are new and take effect when the maintainer accepts them.

Read this before you add or change a commander, a strategy routine, a worker
step, or any code that reads or writes an `sp_project`. It exists because
step bodies and project I/O once accumulated in the commanders: each
probability-table, assembly and commit commander built its own `params` and
builder, re-read the project, and wrote it again, so one iteration paid for four
copies and several full writes. The restructuring of 2026-10-20 removed that
(`doc/refactoring_notes/completed/structure2D_and_project_owner.md`); this policy
keeps it removed. The per-workflow contracts stay in the workflow policies
(`3D/refine3D_policy.md` section 4, `2D/solve2D_policy.md`); this one states the
rules they share.

## 1. Before you add or change a step

1. **Search first.** The step you need is usually there. Run `rg` for the verb
   and for the data it touches (`sample4update_reprod`, `ensure_field_halves`,
   `merge_algndocs`, `read_reprojection_model`, `build_batch_particles3D`,
   `qsys_declare_part_finished`) and read the callers. If a routine does most of
   what you need, give it an argument or a hook; do not copy it.
2. **Name the step's kind** (section 3) and put the body where that kind lives.
3. **Write the project contract in one line** at the routine's header: what it
   reads, what it writes, how often, or "none". If you cannot say it in one
   line, the step is doing too much.
4. **State who owns the project** in the process you are writing (section 5).
   If the answer is "the caller", the routine takes `(params, build)` and does
   not read or write the file.
5. **Test the routine on its own `(params, build)`**, not only through the
   command line, and for a change that should not move results run the
   identity check (`simple_test_exec test=identity`,
   `test_environment_policy.md`).
6. **Run `scripts/layering_audit.py`** once it exists (section 9) and do not
   raise any count.

## 2. Layers

The flow is `ui -> exec -> commander -> strategy/domain object`.

- **exec routers** dispatch a program name to a commander. Nothing else.
- **Commanders** are command objects: `commander_base` defers `execute(cline)`
  and nothing else **(ruled)**. A commander normalizes and validates `cline`,
  calls `params%new(cline)` once, builds the builder, calls a routine, and
  ends. (In the preprocessing family the strategy's `initialize` does this
  instead; the commander only creates the strategy.) It does not own numerics, project choreography or a step that another
  workflow runs.
- **Strategy routines** (`src/main/strategies`, `src/main/solve`) own "how a
  workflow runs": shared-memory versus distributed policy, iteration, the
  master's steps.
- **Domain objects** (`image`, `project`, `ori`, `pftc`, `volume`, `class`,
  `sigma2`, ...) own what the data means and know nothing of commanders or
  strategies.

Imports point downward. A strategy or domain module does not import a concrete
commander module (`simple_commanders_<topic>`) **(proposed)**; the controllers
(`src/main/solve`, the multi-state controller, the final reconstruction, the
external-reference pose initialization) may, to run whole programs (section 4).
The strategy modules of the six master steps never import the commander that
wraps them: a body that needs an in-process commander takes it as a
`commander_base` dummy **(ruled)**. Shared-memory `refine3D` and `refine2D` run
those steps through the commanders' `execute(cline)` **(ruled 2026-10-20)**; that
is the one accepted exception, and it ends if the steps get a shared-memory
entry on the caller's `(params, build)`.

`simple_commanders_api` is an import convenience, not a layer: it re-exports
`builder`, `parameters`, `cmdline`, `sp_project` and the queue helpers. New
strategy code imports those directly (`compile_time_policy.md` forbids
widening the umbrella).

## 3. Kinds of step, and where the body lives

| Kind | Examples | Home | Commander |
|---|---|---|---|
| Program | `refine3D`, `reconstruct3D`, `preprocess` | a controller or strategy (`create_<x>_strategy`) | `execute(cline)` |
| Master step: a part of a stage that the master or another workflow runs in process | volume assembly, `prob_align`(`_neigh`), `prob_align2D`, class-average assembly, sigma2 commit | strategy routine `<x>_exec(params, build, explicit args)` in `strategies/parallelization` with a `run_<x>(cline)` entry **(ruled)** | thin wrapper that calls `run_<x>` |
| Worker step: the per-partition program a scheduler starts | the probability-table workers, the matchers' `refine3D_exec`/`refine2D_exec`, `calc_pspec` parts | the matcher or strategy module beside the master routine for the same step **(proposed)** | thin shell: `params`, builder, body, `qsys_declare_part_finished`, `simple_end` |
| Domain operation | restore class averages, assemble volumes' reduction, NU filter | domain module | none |

A step that is *part of a stage* is never a commander called from a strategy.
Typed entry points beside `execute(cline)` (`execute_from_project`,
`execute_hand_on`) were deleted and do not come back **(ruled)**.

**Variants share one skeleton (proposed).** The probability-table workers for
3D, 3D-neighbourhood and 2D, and any future worker, repeat the same frame:
`params`/builder set-up, reproduction of the master's sample
(`sample4update_reprod`), the empty-partition contribution, batching with
`split_nobjs_even`, the benchmark timers, and the completion call. Write that
frame once; a variant supplies only what differs (which table, which
references, which fill). A new variant that copies the frame is rejected in
review.

**A worker reproduces, never redraws (ruled).** The master draws the sample and
writes it; a worker reproduces it (`sample4update_reprod`) and stops with an
error if the field was not sampled.

**Whole programs inside a program.** A controller (`solve3D`, `refine3D_auto`,
the multi-state controller) may run whole programs as stages by calling their
commanders. Section 4 states what that obliges.

## 4. A commander calling a commander

Allowed only to run a **whole program** as a stage, from a controller. The
caller then owes:

1. **Hand the project over through the file.** Write the segments the callee
   reads before the call (`write_segment_inside(oritype)`), and read back the
   segments it wrote after (`read_segment`). Do not assume the callee shares
   your in-memory project; it builds its own.
2. **Build the child command line sparsely** and set `mkdir=no` once
   `params%new` has created the execution directory. Do not pass native `box`
   or `smpd` for a project run.
3. **Expect `params%new` side effects**: it reseeds the random generator
   (`seed_rnd`) and calls `omp_set_num_threads(nthr)`. The caller saves and
   restores its own thread count. Moved master steps do not call
   `params%new`; they set the threads they need, and the caller restores its own
   **(ruled 2026-10-20)**.
4. **Do not redo the callee's work** in the caller (reference volumes,
   sigma2 set-up, even/odd labelling): if the caller needs the product, ask the
   callee for it or call the routine they share.
5. **A step that is a part of a stage is not called this way** (section 3).

## 5. Project ownership and I/O

The project file is durable workflow state; `build%spproj` is one process's
copy of it. `sp_project` has no file locking: a file written while another
process reads it can be seen half-written.

1. **One owner per project file at a time (ruled for refine3D/refine2D,
   proposed elsewhere).** The distributed master owns its project for the run.
   A stream stage and the 2D pool each own the files they write. Two processes
   never write one file.
2. **Workers read the file and write part files** (orientation documents,
   probability tables, partial sums); a worker does not write the project
   **(proposed)**. The master reads only what the workers write, and merges
   their documents in memory with `merge_algndocs(..., write=.false.)`
   **(ruled)**.
   **Part files have one shape per step family (proposed).** A step that only
   changes micrograph rows (`motion_correct`, `ctf_estimate`, `pick`,
   `gen_pspecs_and_thumbs`, `preprocess`) has each worker write a mic-row part
   document that the master merges with `merge_algndocs`, which checks
   coverage and overlap. A step that creates stacks and particles (`extract`,
   `reextract`, `refine_motion_model`) has each worker write a part project and
   the master assemble them with one shared routine. The stream's job project,
   which the worker updates in place, is the documented exception.
3. **Write only the segments you changed**, with
   `write_segment_inside(oritype)`. A full `write` rewrites the whole file
   (it deletes and recreates it) and belongs to the end of a command-line step
   or of a run **(ruled)**. Note that `write_segment_inside` falls back to a full
   `write` when the file does not exist yet.
4. **One write of a segment per iteration (ruled).** A second write needs a
   reader that must see the file (the `out` segment before postprocessing; the
   particle field before a polish pass); the owning workflow policy names it.
5. **Do not read what you hold.** A routine that has `build` does not read the
   project again; pass `build`, not a file name **(proposed)**. Read a segment
   from the file only when another process has changed it, and read that
   segment, not the project (`read_segment`, `read_non_data_segments`).
6. **Keep `projinfo` honest.** `read` names the file it read. A builder made on
   a project copied or relocated into another directory calls `update_projinfo`
   so that writes go to the run's own project, not the original **(ruled)**; the
   regenerated fixture of the 2026-10-20 run is the incident.
7. **Shared state that must agree across processes has one maker and
   verifiers (proposed).** The master makes it once, in one transaction;
   workers check it and stop if it is absent or inconsistent. The model is
   `ensure_field_halves` (`simple_eo_partition.f90`): the master labels the
   even/odd halves, a worker only verifies. The sample, the sigma2 state path,
   the reprojection model and the class assignment of distributed `make_cavgs`
   (class remap, expansion, time-series split: made by its master and written
   before the workers start) follow the same rule. Never re-derive such state
   in a worker "to be safe".
8. **Register run artifacts in the `out` segment** through its getters and
   adders (`add_vol2os_out`, `add_fsc2os_out`, ...), not by file-name
   convention **(ruled)**.

## 6. Parameters and builders

- One `params%new(cline)` per commander execution path; normalize and validate
  `cline` first **(ruled)**.
- After `params%new`, `params%...` is authoritative; do not re-read parsed
  values from `cline` **(ruled)**.
- Do not build a second builder for the same project in one process. A
  different toolbox on the same project is a builder method (`new_vols`, the
  search toolbox, the band, the model), not a second builder **(proposed)**.
- A routine on the caller's `(params, build)` runs at the thread count it is
  given, saves and restores nothing but its own, and draws from the random
  generator without a reseed **(ruled)**.

## 7. Red flags in review

Each is a signal to ask "what already does this?":

- `params%new` or `init_params_and_build_*` inside a strategy routine that
  already receives `params` and `build`.
- `spproj%read(` of a whole project in code that holds `build%spproj`, or two
  `write_segment_inside` of one segment on one path.
- A worker that calls `write`/`write_segment_inside` on the project file.
- `sample4update_reprod`, `split_nobjs_even` and `qsys_declare_part_finished`
  together in a new commander: a third copy of the worker frame.
- A message or `qsys_declare_part_finished` label that names another module
  (evidence of a copy; `exec_prob_tab` labels itself `simple_commanders_refine3D`).
- `use simple_commanders_<topic>` of a concrete commander module in a strategy
  or domain module, or a routine added to a commander module for a strategy to
  call (`postprocess_volume_from_files`).
- A commander routine of more than about 150 lines that contains numerics or
  project choreography rather than sequencing.
- A new `optional` argument that only one caller omits, added to avoid reading
  the caller's code.

## 8. Known deviations (2026-10-21)

Measured on the working tree at `fd5b5bc8e`. They are the ledger of what the
rules above do not yet describe; each one's fix is a planned item, not a
precedent.

- **Worker frames in commanders.** `exec_prob_tab`, `exec_prob_tab_neigh` and
  `exec_prob_tab2D` (`commanders/simple/simple_commanders_prob.f90`) each carry
  their own frame (100, 113 and 137 lines), and `write_prob_tab_bench` and
  `write_empty_prob_tab` are local to that module. In shared memory,
  `prob_align_exec` writes the sample to the project file and runs the table
  through `xprob_tab%execute(cline)`, which reads the file again.
- **Concrete commander modules imported by strategies.**
  `simple_refine3D_strategy` (`commit_sigma2`, `prob_align`(`_neigh`) by the
  shared-memory ruling; also `rec3D` and `postprocess`),
  `simple_refine2D_strategy` (`commit_sigma2`, `prob_align2D` by the same ruling;
  also `make_cavgs`, `make_cavgs_distr`, `scale`), `simple_calc_pspec_strategy`
  (`calc_pspec_assemble`), `simple_pick_strategy` (`cluster_cavgs`) and
  `simple_rec3D_service` (`postprocess_volume_from_files`, a routine that lives
  in `simple_commanders_volops`). The controllers import them by design
  (`simple_solve3D_utils`, `simple_refine3D_states_controller`,
  `simple_external_reference_pose_initialization`, `simple_sigma2_bootstrap`,
  `simple_flex_pca_project_gateway`).
- **The PCG master writes the project.** `execute_rec3D_pcg_distributed_master`
  writes the particle field once per call, with no `write` flag like the
  gridding path has (`simple_rec3D_pcg_strategy.f90`, after the `res`, `res05`,
  `cfar` update), so a distributed `rec_backend=pcg` iteration writes the field
  once more than `3D/refine3D_policy.md` section 4 counts.
- **Shared memory hands over through the file.** A shared-memory `refine3D`
  iteration writes the particle field at up to five points (after the model,
  before `prob_align`, before a file-based assembly, around a polish pass, at
  the iteration's end) and re-reads it after the `prob_align` commander, because
  each commander it calls builds its own `params` and builder from the file.
- **The 2D master writes `ptcl2D` more than once.** The distributed `refine2D`
  master writes it in `merge_algndocs` (called without `write=.false.`), after
  assembly (`run_distributed_cavg_assembly`) and in `distr_finalize_iteration`,
  and once more after sampling in `prob_align2D_exec`: three writes per
  iteration, four in probabilistic mode, against the one of section 5.4.
  Shared-memory `refine2D` in probabilistic mode writes it twice and re-reads it
  once.
- **A domain module imports a strategy module.** `ptcl_read_session` lives in
  `strategies/search/simple_matcher_ptcl_io.f90` but `class/`, `volume/` and
  `flex/` use it.
- **The preprocessing steps use four patterns and no policy states them.** Mic-row
  part documents merged by `merge_algndocs` (`motion_correct`, `ctf_estimate`,
  `pick`, `gen_pspecs_and_thumbs`, `preprocess`); whole part projects assembled
  by hand in three near-identical copies (`extract`, `reextract`,
  `refine_motion_model`), with count checks only and `extract` copying particle
  fields one by one where the others use `transfer_ori`; a job project updated
  in place by the worker (stream); and stage-owned projects. In the first group
  the worker reads the whole project to update a few rows and the master reads
  twice (a segment to count, then the project for the merge); `extract`'s
  master writes the full project before the workers and again after.
  `segdiam_bin_picker%pick` writes the caller's project twice in one call, and
  the pick master's `make_pickrefs_impl` calls `params%new` a second time and
  runs a commander. In this family the strategy's `initialize`, not the
  commander, calls `params%new`.
- **The stream shares files by sequencing, not by lock.** The master owns no
  project; each stage owns its folder and project and hands on by writing a
  `.tmp` file and renaming it. The exceptions: the pool and its `refine2D` master
  take turns on `refine2D.simple`, and the pool's `finalise` writes the job
  protocol's marker files into the folder of the job it has just cancelled;
  `p01`/`p04` rewrite the `mic` segment of a finished job's project before moving
  it; the master touches `STREAM_FINISHED` in a symlinked earlier folder. Nothing
  but `is_running` and the exit-status files stops two processes from writing one
  stage folder. None of this is written in `stream/`.
- **solve3D add-on.** The final reconstruction rewrites the working project in
  full three times (`switch_addon_membership`, `calc_final_rec_addon`).
- **Long commander routines.** 33 `exec_*` routines exceed 150 lines and 14
  exceed 300 (`exec_refine3D_auto` 707, `exec_solve3D` 597,
  `exec_bootstrap_cavgs` 553, `exec_stream_p00_master` 490,
  `exec_refine3D_states_cavgs` 459, `exec_classify3D_refs` 392, ...). Most
  are controllers that sequence programs; the numerics and the I/O inside them
  are the part to move out when they are next changed.
- **Counts are static.** No test measures project reads and writes per
  iteration; the refine3D contract (section 4 of its policy) was checked by the
  2026-10-20 run, and nothing keeps it checked.

## 9. Checks

Existing: the fast gate (`scripts/run_fast_gate.sh`) runs
`check_test_registry.py`, `check_descr.py` and `check_flex_dag.py`; the
identity and truth suites (`test_environment_policy.md`) catch a moved step that
changes results.

Proposed, not built:

- `scripts/layering_audit.py`, report-only with a committed baseline, so a
  count may fall but not rise: strategy/domain modules that import commander
  modules; commander `exec_*` routines over 150 lines; `params%new` and
  builder construction in strategy modules; `sample4update_reprod` outside the
  matchers and the shared worker frame; `write`/`write_segment_inside` in code
  that runs as a worker. Add it to the fast gate once the baseline is clean of
  false positives.
- Project I/O counters in `sp_project` (reads and writes by segment, printed
  under `L_BENCH_GLOB`), and a high-level test that runs one distributed
  `refine3D` and one `refine2D` iteration and asserts the budgets of the two
  workflow policies.
