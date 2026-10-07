# Stream 3D Ingestion Policy

How the 2D pool hands particles to multistate 3D, how 3D takes them in, and what a change must
preserve.

## 1. Scope

- The pool's publications (p06): `export_pool_state` in
  `src/main/stream/stages/simple_stream_stage_pool2D.f90`, the pool's `publish` in
  `src/main/stream/pool2D/simple_stream_pool2D.f90`, and `build_pool_publication` and
  `delete_pool_publication` in `src/main/stream/pool2D/simple_stream_refine2D_utils.f90`.
- Their import (p07): `import_sets`, `select_cavgs`, `merge_publication` and `take_cavgs` in
  `src/main/stream/stages/simple_stream_stage_solve3D.f90`.
- p07's runs: when it starts `solve3D` and `solve3D_addon` (`next_job`), the first `solve3D`'s
  cap and queue (`draw_first_run`, `release_queue`), what it does with a failed run and an addon
  verdict, and the folders it keeps. The addon program itself is
  governed by `doc/policies/3D/solve3D_addon_policy.md`; this policy relies on its row contract
  (section 4 there) and its report (section 11 there).
- p07's 3D snapshots (section 8): `apply_gui_updates` and `write_snapshot` in the same file.

## 2. What the pool publishes

The pool publishes its state, not deltas. A publication is the pool's classified set as one
completed iteration left it.

1. **The stacks whose particles have been through an iteration** (`updatecnt > 0`), with their
   micrographs and all their particles. Stacks never classified, such as the sets imported since
   the last dispatch, are never published. Within a published stack, a particle never updated (a
   fractional update did not sample it) keeps its row but is published deselected (state 0): its
   class is the random one it was given on import. A later publication selects it once updated.
2. **2D parameters at publication time:** class, in-plane angle, shifts and state in `ptcl2D`.
   `ptcl3D` is prepared as the pool's final project prepares it: a copy of `ptcl2D` with the 2D
   clustering removed (class, in-plane angle, corr, frac), shifts kept, and no `updatecnt` or
   `sampled`. A publication carries the newest optics map's groups (micrographs, stacks,
   particles) and its optics table (`build_pool_publication`); the pool keeps none of its own.
3. **Self-consistent indexing.** Stacks are renumbered from 1, and their particle ranges and the
   particles' stack indices follow. Each particle carries its image index in its stack
   (`indstk`).
4. **Classes:** the pool's `cls2D` table, all classes.
5. **Class averages and FRCs.** The pool's class averages (and their even and odd halves) and
   FRCs at the native sampling: rescaled when the pool runs downscaled, copied otherwise. They are
   registered in `os_out` with the pool's mask diameter.

## 3. When and how it is published

1. **When:** once per completed pool iteration from `EXPORT_START_ITER` (25) on, and once after
   iteration `FIRST_EXPORT_ITER` (10) in a pool that has published nothing yet (`exports_after`,
   follow-up plan, decisions 31 and 37), in the pass where that iteration's results have come back
   and before the next iteration is dispatched. The publication is therefore one consistent
   iteration. The iteration-10 publication carries the sieve's mask diameter: the pool takes it
   when it dispatches iteration 10 (`MSKDIAM_SWITCH_ITER`), so that iteration's class averages were
   still made with the previous one. 3D scores it with the pool model as it scores every
   publication (section 4, item 4): it starts 3D sooner, on a less settled selection that later
   publications correct (section 4, item 5). A restarted pool with publications on disk skips it and resumes at 25, so 3D never
   takes the early classification of a restarted pool as a later publication. Iterations 11 to 24
   publish nothing.
2. **Nothing classified:** a pool with no classified stack publishes nothing for that iteration.
3. **File names:** a publication is `<id>.simple` (5-digit id) in p06's completed folder
   (`DIR_STREAM_COMPLETED`). Its class averages and FRCs are written beside it first
   (`<id>_cavgs.mrcs` with `_even` and `_odd`, `<id>_frcs.bin`), and the project last, under a
   temporary name renamed into place. A reader that finds the project finds it complete.
4. **Ids** increase by one per publication. A restarted p06 continues after the highest id on
   disk (`restore_export_id`).
5. **Retention:** the newest `NPUBLICATIONS_KEPT` (2) publications are kept. When publication `n`
   is written, publication `n - 2` and its files are removed.

## 4. How 3D takes them in

1. **Watching:** p07 watches p06's completed folder and records each publication once.
2. **Only the newest:** in a pass where no 3D job runs, p07 takes the newest publication it has
   not taken. Older ones not taken are passed over, because the newest holds what they held.
   - **A publication p07 cannot use** (no class averages, FRCs or `cls2D`, a count mismatch, a stack
     that changed size, image indices that disagree with the rows: `publication_problem`) is
     passed over with a warning before anything is merged, and listed in
     `rejected_publications.txt` in p07's folder, so a restart passes it over too. The stage
     waits for the next publication instead of stopping.
   - **The mask diameter** is taken from every publication used (`take_mskdiam`): a change is
     logged and applies from the next run.
3. **The first publication** p07 takes (in a fresh session, the iteration-10 one; after a p07
   restart, the newest) is taken as every later one (items 4 and 5). `solve3D` starts in the pass
   that takes it when it selects enough particles (item 6), and no publication is taken while a
   job runs, so the first run sees the first publication's selection alone and later particles
   come in through the addon runs. When that selection is too small for `solve3D`, publications
   are taken until it is not. Until 7 October 2026 p07 classified the first publication's
   particles again with a `solve2D` of its own (the first set; follow-up plan, decisions 31 to
   38), which `doc/refactoring_notes/planned/stream_3D_first_run_plan_2026-10-07.md` retired.
4. **One quality decision per publication.** The publication's class averages are scored
   once with the pool quality model (`CAVG_QUALITY_MODEL_POOL_DEFAULT`), with a mask diameter
   fitted to their box (`fit_mskdiam`). The selection is mapped to the publication's particles by
   class (`map_cavgs_selection`). The selected and rejected class averages are written with JPEGs
   in `quality_selection/<id>/`.
5. **Merge** into the stage's rows (`merge_publication`):
   - **A stack the stage holds** is matched by its stack name, in whatever order the publication
     lists it, and each particle by its image index in the stack (`indstk`; the stage's row must
     record the same image, or the stage stops). Its particles take the publication's 2D
     parameters (`transfer_2Dparams`) and 2D state in place. The 3D state is a run's multistate label: a particle the publication
     deselects gets 0, a selected particle keeps its label, and a selected particle without one
     (state 0) gets 1. They keep their 3D parameters, CTF and optics group. The first
     publication's rows follow every later publication's selection as all others do.
   - **A new stack** is appended with its micrograph and particles (2D and 3D), after the rows
     the stage holds.
   - **A stack the publication lacks** keeps its rows, deselected (state 0 in both segments). This
     happens after a pool restart, until the restarted pool has classified it again.
   - **The stage's `cls2D`** is the publication's, and so is its optics table when the
     publication carries one (the newest optics map's; group ids are kept across maps).
   - **The class averages and FRCs come with the classes** (`take_cavgs`). The publication's are
     copied into `quality_selection/<id>/` and replace the stage's `cavg` and `frc2D` entries in
     its out segment; the volumes and FSCs stay. A run's class-average balancing (`balance=cavg`,
     solve3D's default) reads the class averages and FRCs of the classes its rows are labelled
     with, and p06 removes a publication two publications later (section 3, item 5), possibly
     while a run reads it.
6. **The first run:** `solve3D` starts once at least `MIN_PTCLS_PER_STATE` (5) particles per
   state are selected (`next_job`), in the pass that takes them; with fewer, it waits for later
   publications to select more.
   - **The cap:** it runs on at most `nptcls3D_max` of the selected particles (`draw_first_run`):
     an equal share per selected 2D class, capped at the class's population, best 2D score
     (`corr`) first, trimmed to the cap by the lowest scores. The draw is a greedy
     `sample4update_class` on a copy of `ptcl2D`, so the rows' `updatecnt` and `sampled` stay.
   - **The queue:** the selected particles left out are deselected in the job's project only, so
     they are not among the result's active rows, and selected again once the run is done
     (`release_queue`). They are then the first addon run's cohort (item 7), which starts by the
     usual rule. A run that queued particles has not aligned every selected particle, so the
     final run is never skipped on its account (item 13).
   - **Settings:** the commander checks 2 <= `nstates` <= 20 (the GUI status holds
     `MAX_STATES_SOLVE3D_MULTISTATE`) and `nptcls3D_max` >= 5 × `nstates` before the stage starts.
     The settings come from the command line, defaulting to `NSTATES3D` (3), `NSTAGES3D` (5),
     `LPSTART3D` (50), `LPSTOP3D` (10), `NPTCLS3D_MAX` (100,000), `NPARTS3D` (8) and `NTHR3D` (8),
     with the pool's mask diameter. `nptcls3D_max` is p07's own: the master does not set it.
   - **Selections:** the pool model scores the first publication on the pool's iteration-10
     classification; the later publications start past the pool's low-pass ramp (iteration 20),
     so from then on it scores settled classifications.
7. **Addon runs:** once a result exists, its active rows are the frozen particles, and the
   cohort is the selected rows that are not (appended since, or selected again). An
   `solve3D_addon` run starts when the cohort reaches max(`MIN_PTCLS_PER_STATE` × nstates,
   `ADDON_COHORT_FRAC` (10%) of the frozen particles). The first bound is the addon's own floor
   (its refusals of an empty or small cohort cannot happen); the second bounds the total cost,
   since every run accumulates all frozen particles.
8. **Failures:** a failed `solve3D` stops the stage (`THROW_HARD`). A failed addon run (refused,
   for example for a state without frozen particles, or crashed) leaves the latest result the
   base, and the next run waits for a cohort larger than the one that failed.
9. **The verdict:** after each addon run its report (`solve3D_addon_report.txt`) is read before
   anything is adopted; the verdict per state (IMPROVED, UNCHANGED, REGRESSED) is logged and shown
   in the GUI status. A run with a REGRESSED state is rolled back (follow-up plan, decision 22): the
   frozen base, the stage's project and the GUI's volumes stay the previous run's, and the next
   attempt waits until the cohort has grown by the cadence step again (`retry_cohort`). The status
   says "rolled back" and how many in a row. A run that keeps regressing never replaces the base;
   the final run (item 13) realigns everything at the end. A run that wrote no report is adopted.
10. **Particles are aligned once during a session.** A particle's pose and state, once frozen, are
   not revisited by the addon runs, and they inherit the first `solve3D`'s ladder. There is no
   rebase during a session (decision 21); the final run (item 13) is the one realignment.
11. **After each run** the job's data segments (micrographs, stacks, particles, classes, outputs)
    become the stage's; the stage keeps its own project information, computing environment and
    optics table (`finish_run`). The GUI gets each state's volume, population, resolution, FSC
    curve, reprojections and orientation distribution, all from the result's rows: after the
    first `solve3D`, before its queued particles are selected again. The stage's project is then
    written (under a temporary name, renamed).
    - **The status** reports the particles selected now (`particles_imported`, shown as "particles
      selected"), the particles the latest run took (`particles_at_last_refine`: those selected
      at its start, less the first `solve3D`'s queue), and each state's population and resolution
      in the latest result, kept from `send_volumes`. The rows' current labels are not used, since
      the particles merged since a run carry state 1 before any map holds them. The per-state
      stats appear only once a run has set them.
12. **What p07 keeps:** the newest `NQUALITY_KEPT` (3) quality folders and those of the
    publications a run started from (listed in `quality_selection/runs.txt`, the first run's
    publication among them), with the class averages and FRCs copied into them, so every run's
    project and every 3D snapshot finds its own; the latest addon
    iteration folder (the frozen base), its predecessors removed once it has completed; no folder
    set aside for a job left unfinished once a later run has completed. A stop cancels the
    running job (`doc/policies/stream/restart_policy.md`).
13. **The final run** (decisions 23, 24): the pool flags its final publication (`pool_final=yes`
    in its out segment: the publication after `FINAL_ITER` while the sieve's final set is in the
    pool). Once p07 has merged it and no job runs, it starts one multistate `refine3D` of every
    selected particle from the current state volumes (`vol<s>`), capped by its `lpstop`, in
    `refine3D_final/`, unless the last run that aligned every particle (a `solve3D` that queued
    none, or an earlier final run) had the same number of selected particles. Its result becomes the stage's project
    and the GUI's volumes, with each state's resolution logged before and after; the addon runs
    keep their own base (`solve3D_addon` needs a `solve3D` base), so when the sieve takes finality
    back they go on from it, and the next final publication starts a new final run. A failed final
    run leaves the latest result. A stop cancels it like any job.

## 5. Invariants a change must keep

- **No unclassified particle reaches 3D:** a publication holds only stacks whose particles have
  been through an iteration, and publishes deselected any of their particles never updated.
- **A publication is one completed iteration**, published before the next dispatch.
- **A publication is complete when its name appears** (temporary name, rename), and its class
  averages and FRCs exist before it.
- **Publication ids only grow.**
- **p07's rows only grow and are never renumbered:** row *i* is the same particle image for the
  life of the stage. This is what `solve3D_addon` requires between the frozen and the current
  project.
- **Particles are matched by stack name and image index, never by row**, so a pool restart that
  reorders its stacks changes nothing in p07.
- **A row's 3D parameters, its multistate label, CTF and optics group are changed only by a 3D
  run**, never by a publication, which only deselects or selects (a deselected row loses its
  label).
- **One quality authority per row:** every row's selection is the newest publication's. The
  first `solve3D`'s queue deselects rows in its job's project only, never in the stage's rows.

## 6. Change rules

- A change of what a publication holds, or when it is written, updates sections 2 and 3, and
  checks p07's merge against it.
- A change of the merge keeps the row invariants of section 5 and the addon policy's section 4.
- Tests:
  - `unit_stream` "pool 2D": `test_publication_holds_classified_stacks` (classified stacks only,
    renumbering, image indices, nothing to publish before any iteration),
    `test_export_numbering`, and `test_pause_rules` (`exports_after`: iteration 10 in a fresh pool
    only, every iteration from 25);
  - `unit_stream` "solve 3D": `test_merge_publications` (appending, matching in another order
    and by image index, selection and class in place with the 3D parameters and multistate label
    kept, deselection of a missing stack, the classes), the watch order and mask diameter, the job
    rule (first run, cohort thresholds, a failed run), the cohort count, the retention, the
    volumes sent, `test_take_cavgs`; `test_first_publication` (taken as every later one, its rows
    following a later publication's selection) and `test_first_run_draw` (nothing queued within
    the cap; the class-balanced draw, a small class taken whole, best scores first, the trim to
    the cap, the rows' update counts untouched; the queue released as the first addon run's
    cohort).
  - The class-average selection of a publication has no unit test, and neither has the first
    `solve3D`'s job project with its queued rows deselected.

## 7. Known gaps

- **Cost of a publication:** each one writes the pool's whole classified set, about the size of
  the pool project, once after iteration 10 and once per iteration from 25.
- **Stack matching** is a search from the last match; it is linear per stack only when the
  publication's order differs from the stage's (after a pool restart).
- **p07 never removes rows:** a stack the pool drops for good stays deselected in p07's rows.
- **Retention of 2:** a reader that lags two publications behind finds the files of the older
  one gone. p07 reads the newest in the pass it sees it.
- **The first publication's selection is the least settled:** the pool model scores the
  iteration-10 classes, whose class averages predate the sieve's mask diameter. Later publications
  correct it, but the first `solve3D` has run on it.
- **A p07 restart** takes the newest publication first, so after a restart past iteration 25 the
  first `solve3D` draws from the whole pool of that publication, at most `nptcls3D_max`.
- **The first addon run's cohort is not capped:** it holds every particle the first `solve3D`
  queued, however many.
- **No rebase during a session** (decided): the addon runs never realign frozen particles; only
  the final run does (section 4, items 10 and 13).
- **The final run's command line** (`refine3D` with `vol<s>` and `lpstop`) relies on refine3D's
  own low-pass schedule from the FSC; it has not been run end to end in the stage.
- **The addon run is not tested end to end in the stage:** its command line, frozen project and
  report reading are covered by the addon's own tests, not by p07's.

## 8. 3D snapshots

A 3D snapshot is a particle set the GUI asks for from multistate 3D's latest result: the particles
of one or more selected states. The request (`snapshot3D`) and its report are in the IPC policy
(section 6).

In NICE the states come from the multistate page's state tiles (`_cls3D_state_selector.html`). A
click on a tile only shows its state; the checkbox in its corner puts it in or out of the
selection. Every state starts in. The states left out are kept in `sessionStorage`, because the
zoom page reloads itself every 10 s without an interaction. The button names the states it sends
and is disabled after one click.

1. **The source** is the latest finished run's project (`result_projfile`: solve3D, an addon pass
   or the final refine3D), not the stage's rows. Rows imported since that run have no 3D
   assignment yet. An addon run that is rolled back leaves the previous result as the source.
2. **What it holds:**
   - the particles of the selected states, merged into state 1 in `ptcl3D` and `ptcl2D` alike,
     with their 3D orientations;
   - every other particle deselected;
   - no state artifacts in `out`: no volumes and no FSCs (`remove_state_artifacts_from_osout`, as
     the batch `selection oritype=ptcl3D states=...` gives).
3. **Where:** `snapshots/<name>/<name>.simple` in p07's folder, written with:
   - `<name>_micrographs.star` and `<name>_particles.star`, with the optics table the result
     holds (the newest publication's) and p06's optics offset (`OPTICS_ID_DELTA` per GUI display);
   - each selected state's volume, copied as `vol_state<NN>.mrc` with its original state number,
     which the project does not register.
4. **When:** in the pass after the request arrives, between steps, so never while a job is
   started or collected. Each id is written once.
5. **Not written:** a request before any result, or one whose states hold no particle, is
   answered with no particles and no file, as p06 answers a 2D snapshot of an iteration it no
   longer keeps.
6. **After the stream has finished,** p07 no longer runs. NICE then starts the batch
   `selection oritype=ptcl3D states=...` on p07's final project
   (`solve3D_multistate/solve3D_multistate.simple`). Its result is a batch job, with no snapshot
   folder or STAR files.
7. **Tests:** `test_snapshot3D` in p07's tester: before a result; from a fixture result with
   three states, two of them selected; the same request again.
