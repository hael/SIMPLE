# Stream 3D first run plan, 7 October 2026

p07 (stream multistate 3D) stops re-classifying its first publication with `solve2D`. The first
publication is taken like every later one. The first `solve3D` runs on at most a set number of
selected particles, 100,000 by default. Every other selected particle waits for the first
`solve3D_addon` run.

## Decisions

1. **The cap applies to the first `solve3D` only.** Every selected particle it leaves out joins
   the first addon run's cohort, however many there are. Later addon runs are unchanged.
2. **The first run's particles are drawn class-balanced.** Each selected 2D class gets an equal
   share, capped at its population, with its best 2D scores first.
3. **p06 keeps publishing iteration 10.** p07 scores that publication's class averages with the
   pool model, as it scores every later one. They were made with the previous mask diameter and
   before the low-pass ramp ended, so the first selection is less settled than later ones.
   Later publications correct it, because they now govern every row (item 4 of section 2.1).
4. **The cap is a p07 parameter only.** It is a developer setting on p07's UI entry, like
   `nstates` and `nparts3D`. The master does not pass it, so a stream run uses the default.

Defaults taken without a question:
- **The name** is `nptcls3D_max`, after the stage's 3D job settings (`nparts3D`, `nthr3D`) and
  the `_max` limits (`nsample_max`, `nboxes_max`). Change it before implementing if another name
  reads better.
- **A selection within the cap** gives the first run every selected particle, as today, and
  leaves nothing queued.
- **The first run still starts** once `MIN_PTCLS_PER_STATE` (5) particles per state are
  selected. It takes the first publication's selection alone, because p07 takes no publication
  while a job runs.

## 1. Where things stand

- **The first set** (`take_first_set`, `merge_first_set`; follow-up plan, decisions 31 to 38):
  - The first publication into a stage without rows is merged as it selects its particles.
    The pool model's selection is kept only as a fallback (`first_fallback`).
  - Those rows are the first set. `solve2D` classifies them again in `solve2D/`
    (`start_solve2D`, `PHASE_SOLVE2D`). The chunk model and the class-compatibility filter select
    from its class averages (`finish_solve2D`, `quality_selection/first_set/`), and
    `map2Dshifts23D` makes its shifts the 3D start. A failed `solve2D` falls back to the pool
    model's selection (`fallback_first_set`).
  - Later publications change only the first set's 2D parameters, not its selection
    (`in_first_set` in `merge_publication`).
  - Imports pause from the first set until `solve3D` starts (`l_solve2D_due` in `iterate`).
- **The first `solve3D`** (`start_solve3D`) runs on every selected row and records
  `nptcls_at_full`. The final refine3D is skipped when the selection still has that many
  particles.
- **Addon runs:** a result's active rows are the frozen particles. The cohort is the selected
  rows that are not frozen (`count_cohort`). An addon run starts once the cohort reaches
  max(5 × `nstates`, 10 % of the frozen particles) (`next_job`).
- **Class-balanced draws already exist:**
  - `oris%get_class_sample_stats` lists each class's selected particles, best 2D score
    (`corr`) first.
  - `oris%sample4update_class`, run greedily (`l_greedy`), takes each class's equal quota from
    the top. It also increments `updatecnt` and sets `sampled` on what it draws.

## 2. SIMPLE

### 2.1 Retire the first set (`simple_stream_stage_solve3D.f90`)

1. **Remove:**
   - procedures: `take_first_set`, `merge_first_set`, `in_first_set`, `start_solve2D`,
     `finish_solve2D`, `fallback_first_set`;
   - components: `first_set`, `first_fallback`, `l_solve2D_due`;
   - constants: `PHASE_SOLVE2D`, `SOLVE2D_DIR`, `FIRST_SET_DIR`, `NPTCLS_PER_CLS2D`,
     `NCLS2D_MIN`, `NCLS2D_MAX`, `NSAMPLE2D`, `LPSTOP2D`;
   - the `l_compat` path of `select_with_model`, with the imports only it uses
     (`class_compatibility`, `CAVG_QUALITY_MODEL_CHUNK_DEFAULT`).
2. **`import_sets`** takes every publication the same way: `select_cavgs`, `merge_publication`,
   then `take_cavgs`. The pool model's selection applies to the first publication too.
3. **`iterate`** no longer pauses imports for a due `solve2D`. **`advance_jobs`** has no
   `solve2D` branch. The `solve2D` cases go from **`send_status`** and **`prune_run_dirs`**.
4. **`merge_publication`** has no first-set exception: every row follows each publication's
   selection.
5. **Restart** is unchanged, except that the newest publication is taken as a normal one.

### 2.2 The cap

1. **The parameter** `nptcls3D_max`, following the parameter lifecycle (`simple-main-params`):
   - declared in `simple_parameters.f90` as `0` ("the stage's default"), registered with
     `add_int` in `simple_parameters_parse.f90`. The generated argument list follows from
     `simple_parameters.f90`;
   - given its default `NPTCLS3D_MAX = 100000` in the p07 commander's `set_solve3D_cline`, beside
     `NPARTS3D` and the others, with the header comment updated;
   - checked in `init_params`: a cap under `MIN_PTCLS_PER_STATE` × `nstates` is refused;
   - added to `solve3D_stream`'s UI entry (`simple_ui_stream.f90`), developer visibility,
     `preserve_default`, `{100000}`.
2. **`start_solve3D`**, when more particles are selected than the cap:
   1. A new `draw_first_run(cap)` picks the run's rows:
      - It works on a scratch copy of `os_ptcl2D`, so the rows' `updatecnt` and `sampled` stay
        unchanged.
      - It builds `get_class_sample_stats` over the selected classes and calls greedy
        `sample4update_class` with `update_frac = cap / nselected`.
      - The draw may exceed the cap by one particle per class, so the lowest 2D scores are
        trimmed until the run has exactly the cap.
   2. The selected rows the draw leaves out are recorded in a new component,
      `queued(:)` (logical per row).
   3. They are deselected (state 0) in both `ptcl2D` and `ptcl3D` of the run's project only.
      solve3D's class sampling reads `ptcl2D` states (`make_class_samples` without
      `l_drop_inactive3D`). The stage's rows keep their selection.
   4. `nptcls_at_full` is set to -1, because the run does not align every selected particle.
      The final refine3D is then never skipped on the strength of this run.
   5. The log line gives the selected, drawn and queued counts.
3. **`finish_run`**, for the first `solve3D`, follows this order:
   1. It reads the result as now: the queued rows come back deselected.
   2. It records `frozen_active` from the result, so the queued rows are not frozen.
   3. It selects the queued rows again (state 1 in both segments, the label a selected row
      without one gets) and releases `queued`.

   The queued rows are then cohort (`count_cohort`), so they join the first addon run. That
   run starts by the usual rule: right away when there are enough of them, otherwise once
   later publications add more.
4. **3D snapshots** take the result project, so the queued particles, which have no 3D label
   yet, are not in a snapshot taken before the first addon run.

## 3. Docs

- **`stream_3D_ingestion_policy.md`:**
  - §1 scope: no `solve2D`; `draw_first_run`.
  - §3 item 1: iteration 10 is published early so that 3D starts sooner, and the pool model
    scores its classes (decision 3).
  - §4 item 3: the first publication is taken like the others, and is replaced by the cap.
  - §4 item 5: the first-set exception goes, and so does "the first set's `solve2D` replaces
    them" in the `take_cavgs` bullet.
  - §4 item 6: the first run, the cap, the draw and the queue.
  - §4 item 12: the folders kept, without `solve2D/` and `quality_selection/first_set/`.
- **`restart_policy.md`:** the p07 row and "p07 runs `solve2D` and `solve3D` again".
- **`pool2D_policy.md`:** the iteration-10 publication is no longer "3D's first set".
- **`doc/policies/stream/README.md`:** the p07 row gains `nptcls3D_max`.
- **Code comments** that cite decisions 31 to 38, in the stage and the commander.
- **The `simple-main-stream` skill**, with approval: p07's line, should it mention the cap.

## 4. Tests (`unit_stream`, "solve 3D")

- **Remove** `test_first_set` and `test_first_set_fallback`.
- **Add:**
  - the first publication into a stage without rows merges as it selects, without exceptions,
    and a later publication deselects one of its rows;
  - `draw_first_run` on a fixture with uneven classes and scores:
    - exactly the cap;
    - an equal share per class, capped at a small class's population;
    - best scores first;
    - `updatecnt` and `sampled` unchanged;
  - `finish_run` on a fixture result with queued rows: they are not frozen, they are selected
    again, and `count_cohort` counts them;
  - a selection within the cap queues nothing;
  - the cap check in `init_params`.
- **`test_rules`:** unchanged; `next_job` keeps its signature.
- **The stream chain tester** runs far below the cap and needs no change. Check that it doesn't
  look for `solve2D/`.

## 5. Order

1. Section 2.1 and its tests. Behaviour: no `solve2D`, an uncapped first run.
2. Section 2.2 with its tests, the parameter and its UI entry.
3. The docs.

## Status

Implemented, 7 October 2026: sections 2.1 and 2.2, their tests, and the docs of section 3
(the skill line included) except the README row. Nothing has been compiled or run yet. The plan moves to
`completed/` once `unit_stream` passes and a stream has run its first `solve3D` and addon run.

Where the implementation differs from the plan:
- **The cap is checked in the commander** (`set_solve3D_cline`), beside the `nstates` check,
  rather than in `init_params`. The stage exports `MIN_PTCLS_PER_STATE` for it. The stage itself
  treats a cap below 1 as none, so a tester that builds a stage without the commander draws
  nothing.
- **The draw, the queue's states and its release are three steps** (`draw_first_run`,
  `set_queued_states`, `release_queue`), so the tester can run them without a 3D job.
  `start_solve3D` deselects the queued rows, writes the job's project, and selects them again.
- **`doc/policies/stream/README.md` is unchanged:** its p07 row lists the stage's resources, and
  the cap is not one.
