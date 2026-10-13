# Stream 3D Ingestion Policy

How the 2D pool hands particles to multistate 3D, how 3D takes them in, and what a change must
preserve.

## 1. Scope

- The pool's publications (p06): `export_pool_state` in
  `src/main/stream/stages/simple_stream_stage_pool2D.f90`, the pool's `publish` in
  `src/main/stream/pool2D/simple_stream_pool2D.f90`, and `build_pool_publication`,
  `pool_publication_nselected`, `build_sieve_publication`, `combine_sieve_classes` and
  `delete_pool_publication` in `src/main/stream/pool2D/simple_stream_refine2D_utils.f90`; the
  first set from the sieve's 2D, the stage's `publish_first_set`, `first_set_due`,
  `collect_sieve_set` and `publish_sieve_set` (sections 3 and 9).
- Their import (p07): `import_sets`, `take_first_publication`, `select_cavgs`,
  `merge_publication` and `take_cavgs` in `src/main/stream/stages/simple_stream_stage_solve3D.f90`.
- p07's runs: when it starts `solve3D` and `solve3D_addon` (`next_job`), the first set's
  initial volumes from class averages before the first `solve3D` (`ini3D_cavgs`: solve2D, the
  pool model's selection, cls_expansion, solve3D_cavgs; section 9), the first `solve3D`'s cap
  and queue (`cap_first_run`, `release_queue`), what it does with a failed run and an addon
  verdict, and the folders it keeps. The addon program itself is
  governed by `doc/policies/3D/solve3D_addon_policy.md`; this policy relies on its row contract
  (section 4 there) and its report (section 11 there).
- p07's 3D snapshots (section 8): `apply_gui_updates` and `write_snapshot` in the same file.
- The optional first 3D from the sieve's class averages (section 9): p06's sieve publication and
  p07's `start_cavgs3D`, `finish_cavgs3D` and `take_cavgs3D_result`.

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

1. **When:**
   - **The first set**, in a pool that has published nothing yet, is the sieve's own 2D of the
     sets it has handed off (section 9, items 1 and 2), with their classes and class averages,
     so the first `solve3D` can balance over them (section 4, item 6). It is published as soon
     as it is due, in the pass that finds it so, without waiting for a pool iteration
     (`publish_first_set`, `first_set_due`), at the first of:
     - the sets handed off, in the pool yet or not, hold `NPTCLS_FIRST3D` (100,000) accepted
       (state 1) particles;
     - they hold `NMICS_FIRST3D` (500) micrographs' worth of accepted particles at their average
       rate per micrograph, that is, they cover 500 micrographs;
     - the sieve's final set has come, after which no more particles come;
     - pool iteration `FIRST_EXPORT_ITER` (10) has completed (`export_pool_state`).

     Decided 12 October 2026; until then the first publication was the pool's (the sieve's only
     with `sieve_ini3D=yes`), after the first iteration that selected 100,000 particles or
     iteration 10. The first set's size is p07's cap (`nptcls3D_max`, section 4, item 6), not
     this trigger: a set past the cap has its later stacks deselected for a later publication.
   - **When the first set cannot be built** (section 9, item 2), the pool's classified state
     goes instead, once (`exports_after`): after the first iteration whose publication would
     select `NPTCLS_FIRST3D` particles (`npublishable`, counted as `build_pool_publication`
     selects: updated and kept by the pool), or after iteration 10, whichever comes first. At
     iteration 10 it carries the sieve's mask diameter, which the pool takes when it dispatches
     iteration 10 (`MSKDIAM_SWITCH_ITER`); earlier, the pool's starting diameter.
   - **Later publications** are the pool's, once per completed pool iteration from
     `EXPORT_START_ITER` (25) on. Each is published in the pass where that iteration's results
     have come back and before the next iteration is dispatched, so a publication is one
     consistent iteration. Iterations 11 to 24 publish nothing.

   3D takes the first publication as published (section 4, item 3). A restarted pool with
   publications on disk skips the first and resumes at 25, so 3D never takes the early
   classification of a restarted pool as a later publication.
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
3. **The first publication** p07 takes, into a stage without rows (in a fresh session, the first
   set, section 3, item 1; after a p07 restart, the newest), is taken as published, without the
   pool model (`take_first_publication`):
   - its selected particles (the sieve's selection; the pool's updated and kept particles when
     p06 published its own) are merged as selected;
   - at most `nptcls3D_max` stay selected, whole stacks in order (`cap_first_run`, item 6); the
     others are deselected, for a later publication's model selection to take (they then come in
     through an addon run);
   - the rows left selected are the **first set**: their selection stays for the session, and
     later publications change only their 2D parameters (item 5).

   `solve3D` starts in the pass that takes it when it selects enough particles (item 6), and no
   publication is taken while a job runs, so the first run sees the first set alone. When the
   first set is too small for `solve3D`, later publications are taken until it is not. Until 7
   October 2026 p07 classified the first publication's particles again with a `solve2D` of its
   own (follow-up plan, decisions 31 to 38); that was retired, and on 8 October 2026 the first
   publication stopped going through the pool model.
4. **One quality decision per later publication.** The publication's class averages are scored
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
     (state 0) gets 1. They keep their 3D parameters, CTF and optics group. **A first-set row**
     takes the 2D parameters only: its selection and multistate label stay.
   - **A new stack** is appended with its micrograph and particles (2D and 3D), after the rows
     the stage holds.
   - **A stack the publication lacks** keeps its rows, deselected (state 0 in both segments). This
     happens after a pool restart, until the restarted pool has classified it again. The first
     set's rows keep their selection.
   - **The stage's `cls2D`** is the publication's, and so is its optics table when the
     publication carries one (the newest optics map's; group ids are kept across maps).
   - **The class averages and FRCs come with the classes** (`take_cavgs`). The publication's are
     copied into `quality_selection/<id>/` and replace the stage's `cavg` and `frc2D` entries in
     its out segment; the volumes and FSCs stay. p07's runs balance their samples over them
     (solve3D's default `balance=cavg`, item 6), and any reader of a run's project or a 3D snapshot needs the
     class averages and FRCs of the classes its rows are labelled with, and p06 removes a
     publication two publications later (section 3, item 5).
6. **The first run:** `solve3D` starts once at least `MIN_PTCLS_PER_STATE` (5) particles per
   state are selected (`next_job`), in the pass that takes them; with fewer, it waits for later
   publications to select more. With `ini3D_cavgs=yes` (the default) the first set's initial
   volumes come first, and `solve3D` starts from them when they are done (section 9).
   - **The cap:** the first set already holds at most `nptcls3D_max` particles (item 3). The cap
     (`cap_first_run`) takes whole stacks in row order (the pool's stack order, the earliest
     imported first) up to the first whose selected particles no longer fit, so no micrograph is
     split and the run may hold fewer than the cap. When not even the first stack fits (a cap
     below one micrograph's particles), it takes that stack's first selected particles up to the
     cap. No 2D class or score is consulted.
   - **The queue:** when the first set was too small and later publications took the selection
     past the cap before `solve3D` started, `solve3D` applies the cap again. The selected particles left
     out are then deselected in the job's project only, so they are not among the result's
     active rows, and selected again once the run is done (`release_queue`). They are the first
     addon run's cohort (item 7), which starts by the usual rule. A run that queued particles has
     not aligned every selected particle, so the final run is never skipped on its account (item
     13).
   - **Settings:** the commander checks 2 <= `nstates` <= 20 (the GUI status holds
     `MAX_STATES_SOLVE3D_MULTISTATE`) and `nptcls3D_max` >= 5 × `nstates` before the stage starts.
     p07 gives no `balance`, so its runs sample by solve3D's default (`set_balance_default`):
     `balance=cavg`, since every publication brings classes, class averages and FRCs. The
     fractional-update samples are balanced over groups of class averages (`nclust`, 20, by
     default): the first `solve3D` over the first set's, the sieve's classes of the first
     publication. The manifest records only keys given on the base run's command line, so it
     holds no `balance`, and each addon run takes the same default over the classes of the
     newest publication taken. Decided 12 October 2026; from 8 October every run was given
     `balance=none`.
     The settings come from the command line, defaulting to `NSTATES3D` (3), `NSTAGES3D` (5),
     `NPTCLS3D_MAX` (30,000), `NPARTS3D` (8) and `NTHR3D` (8), with the pool's mask diameter.
     p07 passes no low-pass limits: `solve3D` and the final `refine3D` set their own. `nptcls3D_max` is p07's own: the master does not set it, and
     it is independent of p06's first-set trigger (section 3, item 1), which may publish more.
   - **Selections:** the first set is the sieve's selection (its chunks' classification and
     rejection), with no pool model decision; the later publications start past the pool's
     low-pass ramp (iteration 20), so the pool model scores settled classifications.
7. **Addon runs:** once a result exists, its active rows are the frozen particles, and the
   cohort is the selected rows that are not (appended since, or selected again). An
   `solve3D_addon` run starts when the cohort reaches max(`MIN_PTCLS_PER_STATE` × nstates,
   `ADDON_COHORT_FRAC` (10%) of the frozen particles). The first bound is the addon's own floor
   (its refusals of an empty or small cohort cannot happen); the second bounds the total cost,
   since every run accumulates all frozen particles.
   - **The cohort cap:** a run aligns at most `nptcls_addon_max` (`NPTCLS_ADDON_MAX`, 60,000) of
     its cohort (`cap_addon_cohort`, through the same whole-stacks rule as the first run's cap,
     `cap_rows`): the oldest stacks first, up to the first whose cohort rows no longer fit, so no
     micrograph is split and a run may hold fewer. The frozen particles are always in the run.
     The cohort left out is deselected in the job's project only and selected again when the run
     is done (`release_queue`), so it is the next run's cohort, which starts by the usual rule; a
     backlog drains one capped run after another. A failed or rolled-back run drops the queue
     (the stage's rows are unchanged). `nptcls_addon_max` is p07's own, checked to be at least
     `MIN_PTCLS_PER_STATE` × nstates by the commander; the master does not set it.
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
    selected particle from the current state volumes (`vol<s>`), in
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
- **One quality authority per row:** a first-set row's selection is the first publication's,
  as the sieve made it (the pool, when p06 published its own); every other row's is the newest
  publication's, as the pool model made it.
  The first `solve3D`'s queue deselects rows in its job's project only, never in the stage's rows.

## 6. Change rules

- A change of what a publication holds, or when it is written, updates sections 2 and 3, and
  checks p07's merge against it.
- A change of the merge keeps the row invariants of section 5 and the addon policy's section 4.
- Tests:
  - `unit_stream` "pool 2D": `test_publication_holds_classified_stacks` (classified stacks only,
    renumbering, image indices, the selected particles counted as published, nothing to publish
    before any iteration), `test_export_numbering`, `test_pause_rules` (`exports_after`: in a
    fresh pool only, iteration 10 or earlier once `NPTCLS_FIRST3D` are selected; every iteration
    from 25; `first_set_due`: `NPTCLS_FIRST3D` accepted particles handed off, 500 micrographs,
    or the final set) and `test_sieve_final_set` (the accepted particles, their micrographs and
    the final set noted as handed off; a
    first set that cannot be built fails once, publishes nothing and leaves it to the pool);
  - `unit_stream` "solve 3D": `test_merge_publications` (appending, matching in another order
    and by image index, selection and class in place with the 3D parameters and multistate label
    kept, deselection of a missing stack, the classes), the watch order and mask diameter, the job
    rule (first run, cohort thresholds, a failed run), the cohort count, the retention, the
    volumes sent, `test_take_cavgs`; `test_first_publication` (taken as the pool selects it, the
    stack that no longer fits under the cap deselected, the first set kept through later
    publications and a missing stack, the deselected rows following them) and
    `test_first_run_cap` (nothing queued within the cap; whole stacks in order, deselected rows
    skipped, fewer than the cap rather than a split stack, a first stack over the cap cut to it;
    the queue released as the first addon run's cohort) and `test_addon_cohort_cap` (an addon
    run's cohort capped in whole stacks, the oldest first, frozen particles never queued).
  - The class-average selection of a publication has no unit test, and neither has the first
    `solve3D`'s job project with its queued rows deselected.

## 7. Known gaps

- **Cost of a publication:** each one writes the pool's whole classified set, about the size of
  the pool project, once per iteration from 25; the first set holds the sieve's sets instead.
- **Stack matching** is a search from the last match; it is linear per stack only when the
  publication's order differs from the stage's (after a pool restart).
- **p07 never removes rows:** a stack the pool drops for good stays deselected in p07's rows.
- **Retention of 2:** a reader that lags two publications behind finds the files of the older
  one gone. p07 reads the newest in the pass it sees it.
- **The first set's selection is the least settled and stays:** it is the sieve's selection,
  from its chunks' classification of the first sets handed off, with no pool model
  decision, and later publications do not change it.
- **The first run balances over unaligned classes:** the sieve's chunks are classified
  separately, so the same view appears in the class averages of many chunks. `balance=cavg`
  groups them by correlation (`cluster_cavgs`, `nclust` groups), which merges such repeats, but
  its cost grows with the square of the number of class averages. Every addon run groups the
  newest publication's classes again.
- **A p07 restart** takes the newest publication first, as the pool selects it, so after a
  restart past iteration 25 the first set is that publication's earliest stacks (at most
  `nptcls3D_max`), again without the pool model.
- **A large backlog takes several addon runs:** with `nptcls_addon_max` (60,000), the particles
  the first `solve3D` queued, or a large publication, are aligned over successive capped runs, the
  oldest first; the failure and rollback rules count the whole cohort, not the capped one.
- **No rebase during a session** (decided): the addon runs never realign frozen particles; only
  the final run does (section 4, items 10 and 13).
- **The final run's command line** (`refine3D` with `vol<s>`) relies on refine3D's own
  low-pass schedule from the FSC; it has not been run end to end in the stage.
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

## 9. The first set from the sieve's 2D, and its initial volumes from class averages (`ini3D_cavgs`)

The first publication is the sieve's own 2D (items 1 and 2; since 12 October 2026). p07 then makes
the first set's initial volumes from class averages of its own (items 3 to 9), unless its
`ini3D_cavgs=no`. Decided 12 October 2026: this replaced the route of 8 October, on which p06
marked the first publication (`sieve_ini3D=yes`, a master input passed to p06) and p07 ran
`solve3D_cavgs` on the sieve's combined classes
(`doc/refactoring_notes/planned/stream_sieve_ini3D_plan_2026-10-08.md`). The marker, the master's
input and its forwarding are gone.

1. **p06 collects the sieve's 2D.** While it has published nothing, each set the sieve hands off
   (as it is handed off, before the pool takes it) has its class averages (main, `_even`,
   `_odd`) and FRCs copied into `sieve_cavgs/<set stem>/` (`collect_sieve_set`), with the first
   mask diameter and pixel size the sets carry; its accepted particles and micrographs are
   counted (`nptcls_sieved`, `nmics_sieved`) and the sieve's final set is noted
   (`l_sieve_final_seen`) for the trigger
   (section 3, item 1). A restarted pool with publications on disk collects nothing.
2. **The first publication is the sieve's 2D** (`publish_sieve_set`), when the first set is due
   (section 3, item 1):
   - every set the sieve has handed off, in the pool yet or not, re-read whole, combined (`build_sieve_publication`): the
     stacks renumbered with each particle's image index, the sieve's 2D parameters and selection,
     each set's class labels offset by the classes of the sets before it, the class tables
     concatenated with the sieve's states and the populations of the selected particles;
   - the class averages, halves and FRCs concatenated in set order beside the publication under
     the usual names (`combine_sieve_classes`), refused when a file is missing, a stack does not
     hold its set's classes, or the boxes or FRC sizes differ;
   - the newest optics map's groups and the sieve's mask diameter.

   It is tried once: when it cannot be built, the pool's publication goes instead, at the pool's
   trigger (section 3, item 1), with a warning. `sieve_cavgs/` is removed after the attempt.
   Later publications are the pool's.
3. **p07 takes it as any first publication** (section 4, item 3: the first set, at most
   `nptcls3D_max`, whole stacks in order). With `ini3D_cavgs=yes` (p07's input, default yes) the
   first run that is due (section 4, item 6) is the route below instead of `solve3D`. Each step
   runs as a job on the 3D jobs' parts and threads, in the pass the previous one finished, and no
   publication is taken while one runs.
4. **solve2D** (`start_solve2D_first`, `solve2D/`) classifies every selected row into
   `NCLS2D_FIRST` (150) classes with the pool's mask diameter and `cache=yes`, and solve2D's
   defaults otherwise. Every selected row is classified, not only the
   first run's: `solve3D` caps them afterwards (item 9), so the rows it queues carry the class
   averages' poses too. Its classes, class averages, halves and FRCs and the rows' 2D parameters
   become the stage's (`adopt_job_2D`); a result without them, or with other rows, counts as
   failed.
5. **The pool model's selection** (`select_first_set_classes`): the classes are scored with the
   pool quality model (`select_with_model`, `CAVG_QUALITY_MODEL_POOL_DEFAULT`), and the
   selection is mapped to their particles; the selected and rejected stacks go to
   `solve2D/quality_selection/`. **A rejected class's particles are deselected for the
   session:** they stay first-set rows, whose selection later publications do not change. When
   the model cannot score the classes or keeps none, every populated class is kept, with a
   warning.
6. **cls_expansion** (`start_clsexp`, `cls_expansion/`) splits each kept class into the same
   number of subclasses, about `NCLS_EXPANDED` (500) in all: `nint(500 / kept classes)`, at
   least 2 (`subclasses_per_class`), with the mask diameter fitted to the particles' box. Its
   subclasses, their class averages, halves and FRCs, and the particles' subclass labels become
   the stage's (`adopt_job_2D`).
7. **solve3D_cavgs** (`start_cavgs3D`, `solve3D_cavgs/`) makes `nstates` volumes from those class
   averages (`pgrp=c1`, `prune=no`, the mask diameter, its own low-pass limits from the FRCs; one
   run, no collapse restarts).
8. **Its result** (`take_cavgs3D_result`): each class's pose and state is mapped to its selected
   particles (`map2ptcls`); the states the selected particles hold are renumbered from 1
   (`compact_states`) in both particle segments and the classes, and become the stage's states
   (`params%nstates`) for the session, addon runs and the final run included. A first-set row
   the mapping deselects (its class has no state) stays a first-set row, deselected for the
   session.
9. **solve3D starts in the same pass** from those poses and states (`cavg_ini_ext=yes`, after the
   symmetry stage, which needs every state populated), with the cap, the queue and the default
   `balance=cavg`, over the subclasses, as always (section 4, item 6).
10. **A failed step is skipped,** with a warning: a failed or unusable solve2D, or a failed
    `solve3D_cavgs` or one that leaves no class with a state, starts `solve3D` on the first set as
    without the route; a failed or unusable cls_expansion runs `solve3D_cavgs` on the classes the
    model kept.
11. **The GUI** shows "running solve2D on the first set", "running cls_expansion on the first set"
    and "running solve3D_cavgs", each with the solve3D progress of 1. When `solve3D_cavgs` is
    done, its initial volumes (`vol_cavg`) go to the GUI under the stage's renumbered states
    (`send_cavgs3D_volumes`), with their low-pass copies, reprojections and the orientation
    distributions and populations of the particles mapped to them, but no FSC or resolution yet.
    They are shown only: the stage's project does not register them, and the first `solve3D`'s
    result replaces them.
12. **Restart:** a p07 restart takes its first publication anew and runs the route again; a job
    folder a stopped step left is set aside (`fresh_job_dir`).
13. **Tests:** `test_sieve_publication` in p06's tester (the offsets, selection, class tables and
    populations; the combined files; a box mismatch and a missing half refused); in p07's,
    `test_cavgs3D` (the subclasses per class, the compaction, the mapping onto the first set),
    `test_first_set_deselection_stays` (a deselected first-set row is not selected again by a
    later publication), `test_send_cavgs3D_volumes` (the initial volumes sent under the
    renumbered states, without FSCs, not registered) and `test_take_cavgs`. The jobs and the pool model's selection of the
    first set run only in a stream.
