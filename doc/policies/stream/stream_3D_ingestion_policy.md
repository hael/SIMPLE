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
- Their import (p07): `import_sets`, `take_first_publication`, `merge_publication` and
  `take_cavgs` in `src/main/stream/stages/simple_stream_stage_solve3D.f90`.
- p07's jobs (`next_job`): the route to the reference volumes on the first set, their check,
  the further attempts and the fallback (section 9); each later publication's assignment to the
  references, by cls_expansion, matching, rejection by fit, mapping and a local refinement of the
  particles' poses with the states fixed (section 10); what it does with a failed job, and the folders it keeps. The matching itself is
  `simple_cavgs_state_match` (`src/main/class/`), run by the `simple_private_exec` program
  `assign_cavgs_states` (`exec_assign_cavgs_states` in `simple_commanders_cavgs`).
- p07's 3D snapshots (section 8): `apply_gui_updates` and `write_snapshot` in the same file.
- The method's rules, constants, steps on projects and job command lines, without the stage:
  `simple_state_assignment` (`src/main/solve/`, since 15 October 2026), shared with the batch
  program `solve3D_pyramid` (`doc/policies/3D/solve3D_pyramid_policy.md`). The stage keeps its
  phases, jobs, retries, fallbacks, publications and GUI, and calls the shared module for the
  rest: the pool model's scoring (`score_cavgs_with_model`, `select_route_classes`,
  `model_rejected_rows`), the check (`check_state_refs`, `refs_usable`, `docked_pair_cc`,
  `read_refs_record`), the assignment (`map_assignment`, `reject_by_fit`, `flag_ambiguous`,
  `write_matched_cavgs`), the local refinement (`set_refine_seeds`, `write_refine_project`,
  `map_refine_result`) and the `set_*_cline` builders. A change there changes both programs.

Since 13 October 2026 p07 runs no `solve3D`, `solve3D_addon` or final `refine3D`
(`doc/refactoring_notes/completed/stream_p07_state_assignment_plan_2026-10-13.md`); the inputs
`nptcls3D_max`, `nptcls_addon_max`, `nstages` and `ini3D_cavgs` went with them.

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
     from which p07 makes its reference volumes (section 9). It is published as soon
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
     iteration 10. p07 takes the first set whole (section 4, item 3); from 12 to 13 October 2026
     it kept at most `nptcls3D_max` of it.
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
2. **Only the newest:** in a pass where no job runs, p07 takes the newest publication it has
   not taken. Older ones not taken are passed over, because the newest holds what they held.
   - **A publication p07 cannot use** (no class averages, FRCs or `cls2D`, a count mismatch, a stack
     that changed size, image indices that disagree with the rows: `publication_problem`) is
     passed over with a warning before anything is merged, and listed in
     `rejected_publications.txt` in p07's folder, so a restart passes it over too. The stage
     waits for the next publication instead of stopping.
   - **The mask diameter** is taken from every publication used (`take_mskdiam`): a change is
     logged and applies from the next job.
3. **The first publication** p07 takes, into a stage without rows (in a fresh session, the first
   set, section 3, item 1; after a p07 restart, the newest), is taken as published, without the
   pool model (`take_first_publication`): its selected particles (the sieve's selection; the
   pool's updated and kept particles when p06 published its own) are merged as selected, all of
   them, and are the **first set**: their selection stays for the session, and later
   publications change only their 2D parameters (item 5). Without references the route to them
   is due on it (section 9); with references (a restart that found them) its assignment is
   (section 10). No publication is taken while a job runs, so the route sees the first set
   alone. Until 7 October 2026 p07 classified the first publication's particles again with a
   `solve2D` of its own; on 8 October 2026 the first publication stopped going through the pool
   model; from 12 to 13 October 2026 at most `nptcls3D_max` of it stayed selected.
4. **The pool model is logged, not applied, on later publications** (since 13 October 2026). A
   publication's class averages are scored once with the pool quality model
   (`score_cavgs_with_model`, `CAVG_QUALITY_MODEL_POOL_DEFAULT`), with a mask diameter fitted to their
   box (`fit_mskdiam`), and the selected and rejected class averages written with JPEGs in
   `quality_selection/<id>/`. The verdict changes no row: the selected rows a rejected class
   holds are flagged (`flag_model_rows`) and counted beside the fit's rejection once the
   publication is assigned (section 10, item 5), so the two can be compared before the model is
   dropped. `L_MODEL_LOGGED` turns the scoring off. The rows' selection is the pool's, as
   published, less the fit's rejection.
5. **Merge** into the stage's rows (`merge_publication`):
   - **A stack the stage holds** is matched by its stack name, in whatever order the publication
     lists it, and each particle by its image index in the stack (`indstk`; the stage's row must
     record the same image, or the stage stops). Its particles take the publication's 2D
     parameters (`transfer_2Dparams`) and 2D state in place. The 3D state is an assignment's: a
     particle the publication deselects gets 0, a selected particle keeps its state, and a
     selected particle without one (state 0) gets 1 until it is assigned. They keep their 3D
     parameters, CTF and optics group. **A first-set row** takes the 2D parameters only: its
     selection and state stay, except that a row the latest assignment deselected (its
     sub-class rejected by fit, or without a class of the matching; `fit_desel`) is selected
     again, as state 1, and judged again by the next assignment.
   - **A new stack** is appended with its micrograph and particles (2D and 3D), after the rows
     the stage holds.
   - **A stack the publication lacks** keeps its rows, deselected (state 0 in both segments). This
     happens after a pool restart, until the restarted pool has classified it again. The first
     set's rows keep their selection.
   - **The stage's `cls2D`** is the publication's, and so is its optics table when the
     publication carries one (the newest optics map's; group ids are kept across maps).
   - **The class averages and FRCs come with the classes** (`take_cavgs`). The publication's are
     copied, with their halves, into `quality_selection/<id>/` and replace the stage's `cavg` and
     `frc2D` entries in its out segment; the volumes and FSCs stay. The route's solve2D, the
     fallback's `solve3D_cavgs` and the expansion read them, and an assignment whose expansion
     fails matches them; p06 removes a publication two publications later (section 3, item 5).
6. **The jobs** (`next_job`), each in the pass that makes it due, none while another runs:
   - without references, once `MIN_PTCLS_PER_STATE` (5) particles per state (`nstates`) are
     selected: the route when it is due (the first set taken, or an attempt failed with attempts
     left), else the fallback when it is due (section 9);
   - with references, a publication's assignment when one is due and 5 particles per reference
     are selected (section 10).

   **Settings:** the commander checks 2 <= `nstates` <= 20 (the GUI status holds
   `MAX_STATES_SOLVE3D_MULTISTATE`); `nstates` (default `NSTATES3D`, 3) is the reference
   volumes' states at most, and the states that hold class averages are kept. The jobs run on
   `nparts3D` parts (default 8; a job is distributed only with more than one) and `nthr3D`
   threads (default 8), with the pool's mask diameter; the matching runs in shared memory on
   `nthr3D` threads. Every other setting is a named constant of the stage.
7. **Failures** never stop the stage: a failed step of the route is skipped (section 9, item 11);
   a failed or unusable `solve3D_cavgs` is an attempt without volumes (section 9, item 10); a
   publication whose expansion fails has its pool class averages matched; a failed matching
   leaves the states as they were until the next publication; a failed local refinement leaves
   the assignment at the classes' poses, with the previous volumes.
8. **After each local refinement** the volumes and FSCs refine3D registered in its project, and
   the refined poses, become the stage's, in place of every earlier state's
   (`take_refine_result`); the stage's project is written into the refinement's folder as the
   snapshots' source (`result_projfile`) and as the stage's project (under a temporary name,
   renamed); the GUI gets each state's volume, population, resolution, FSC curve and orientation
   distribution (`send_volumes`).
   - **The status** reports the particles selected now (`particles_imported`, shown as "particles
     selected"), the particles the latest assignment gave a state (`particles_at_last_refine`),
     the assignments started (`refine_iteration`), and, once the references exist, each state's
     population and resolution as sent with its volume (none for the references). The rows'
     current labels are not used, since the particles merged since an assignment carry state 1
     until the next.
9. **What p07 keeps:** the newest `NQUALITY_KEPT` (3) quality folders and those of the
   publications the route or the fallback started from (listed in `quality_selection/runs.txt`),
   with the class averages and FRCs copied into them; `state_refs/` (each attempt's volumes in
   `try_<n>/`, the references and their record); the latest assignment's folder
   (`assign/it_<n>/`), its predecessors removed once it is refined; no folder set aside for a job
   left unfinished once an assignment is refined. A stop cancels the running job
   (`doc/policies/stream/restart_policy.md`).
10. **The pool's final publication** (`pool_final=yes` in its out segment) is assigned like any
    other, and nothing follows it. Taken before the references exist, it is assigned once they
    do (`accept_refs`); while the route's attempts fail it also ends the wait for a fallback
    publication (section 9, item 10).
11. **Retired on 13 October 2026:** `solve3D` on the first set (with the first set's cap and
    queue, `cavg_ini_ext` from the route's poses), the `solve3D_addon` runs on each growth (their
    cohort cap, verdict and rollback) and the final multistate `refine3D`. The policy as it was is
    in the git history of this file.

## 5. Invariants a change must keep

- **No unclassified particle reaches 3D:** a publication holds only stacks whose particles have
  been through an iteration, and publishes deselected any of their particles never updated.
- **A publication is one completed iteration**, published before the next dispatch.
- **A publication is complete when its name appears** (temporary name, rename), and its class
  averages and FRCs exist before it.
- **Publication ids only grow.**
- **p07's rows only grow and are never renumbered:** row *i* is the same particle image for the
  life of the stage. Every job's project is a copy of the rows, and the assignment maps its
  result back by row.
- **Particles are matched by stack name and image index, never by row**, so a pool restart that
  reorders its stacks changes nothing in p07.
- **A row's 3D parameters, its state, CTF and optics group are changed only by an assignment**,
  never by a publication, which only deselects or selects (a deselected row loses its state).
- **The references never change once accepted**: every publication is matched
  against the same volumes, and the refinements from the particles never feed back into them.
- **One quality authority per row:** a first-set row's selection is the first publication's, as
  the sieve made it (the pool, when p06 published its own), less the route's model selection,
  which deselects it for the session; every other row's is the newest publication's, as the pool
  made it. The fit's rejection is a per-publication judgement for every row: it deselects until
  the next publication, which brings the row back to be judged again (since 14 October 2026; a
  first-set row was deselected for the session before).

## 6. Change rules

- A change of what a publication holds, or when it is written, updates sections 2 and 3, and
  checks p07's merge against it.
- A change of the merge keeps the row invariants of section 5.
- A change of the check, the rejection or the ambiguity flag keeps its rule a pure function
  (`refs_usable`, `reject_by_fit`, `flag_ambiguous` in `simple_state_assignment`; `best_attempt`
  of the stage), tested on fixed numbers; it is shared with `solve3D_pyramid`.
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
    and by image index, selection and class in place with the 3D parameters and state kept,
    deselection of a missing stack, the classes), the watch order and mask diameter,
    `test_take_cavgs`, `test_first_publication` (taken whole as the pool selects it, the first
    set kept through later publications and a missing stack, the other rows following them),
    `test_rules` (the job rule, the best attempt, the final publication), `test_reference_volumes`,
    `test_failed_attempts` (section 9), `test_assignment` (section 10), the retention, the status,
    the volumes and the references sent;
  - `unit_reconstruction` "state assignment rules" (`simple_state_assignment_tester`): the check,
    the renumbering, the sub-class budget, the mask diameter, the fit rules (the outlier cut, the
    30% cap, the worst first, no non-candidate rejected, none when all fit; the ambiguity flag),
    the first set's cut at whole stacks, the class averages by state and the jobs' command lines;
  - `unit_pftc_align2D3D` "cavgs state match": the matching (section 10, item 3);
  - `lib_stream` "sieve to 3D" and "movies to 3D" (`simple_stream_chain_tester`) run the jobs,
    until the first assignment is refined.
  - The scoring of a publication with the pool model has no unit test.

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
  decision but the route's, and later publications do not change it.
- **The assignment is only as good as the references:** the check refuses halves that disagree,
  nearly empty states and near-duplicate states (at 15 Å), but the references are never
  refreshed.
- **The rejection by fit may follow the sub-classes' sizes:** smaller sub-class averages are
  noisier and correlate less at 15 Å whatever their fit. Each assignment logs the fit's median,
  deviation and threshold, each rejected class's correlation and population, the median
  population of the kept and rejected candidates, and the correlation of the best correlation
  with the log population over the candidates; a population-aware statistic waits for those
  numbers from a validation run.
- **The particles' poses start from the class averages,** composed with each particle's 2D
  parameters, and are refined only locally (`refine=cont`: `athres_cont`, 10 degrees, and `trs`):
  a particle further off than that, or in the wrong state, stays wrong. Each assignment starts
  again from the class poses, so the refinement never accumulates over the session. It is one
  iteration (`maxits=1`) over every selected particle of every publication, still the most
  expensive step of an assignment.
- **The matching is at a fixed 15 Å** (`LP_ASSIGN`): states that differ only at higher
  resolution are not told apart. Its shift limit is the program's default (10 pixels), not the
  pool's.
- **Every publication is expanded,** which reads every particle of it, much like a pool
  iteration; when a publication's three jobs take longer than a pool iteration, p07 assigns only
  the newest waiting, so the assignment lags the pool rather than queueing behind it.
- **The rejection by fit discards what the references lack** (a conformation or views the first
  set did not hold), and blurred junk may fit low-pass projections well enough to be kept; the
  pool model's verdict is logged beside it until the comparison settles which to keep.
- **A p07 restart** takes the newest publication first, as the pool selects it, without the pool
  model; with references on disk it is assigned at once, otherwise the route runs on it.
- **The jobs are not tested end to end in the stage's unit tests:** the command lines of
  cls_expansion, `assign_cavgs_states` and refine3D, and the reading of their results, are
  exercised only by the chained stream tests.

## 8. 3D snapshots

A 3D snapshot is a particle set the GUI asks for from multistate 3D's latest result: the particles
of one or more selected states. The request (`snapshot3D`) and its report are in the IPC policy
(section 6).

In NICE the states come from the multistate page's state tiles (`_cls3D_state_selector.html`). A
click on a tile only shows its state; the checkbox in its corner puts it in or out of the
selection. Every state starts in. The states left out are kept in `sessionStorage`, because the
zoom page reloads itself every 10 s without an interaction. The button names the states it sends
and is disabled after one click.

1. **The source** is the latest local refinement's rows (`result_projfile`,
   `assign/it_<n>/assignment.simple`), not the stage's current rows. Rows imported since that
   assignment have none yet. A failed refinement leaves the previous one as the source.
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

## 9. The first set from the sieve's 2D, and the reference volumes from class averages

The first publication is the sieve's own 2D (items 1 and 2; since 12 October 2026). p07 then makes
the states' reference volumes once, from class averages of its own (items 3 to 13). Decided 12
October 2026 as the initial volumes of the first `solve3D` (`ini3D_cavgs`), and since 13 October
2026 the references every publication is assigned to (section 10); it replaced the route of 8
October, on which p06 marked the first publication (`sieve_ini3D=yes`) and p07 ran
`solve3D_cavgs` on the sieve's combined classes
(`doc/refactoring_notes/rejected/stream_sieve_ini3D_plan_2026-10-08.md`). The marker, the master's
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
3. **p07 takes it as any first publication** (section 4, item 3), and the route below is due on
   it. Each step runs as a job on the 3D jobs' parts and threads, in the pass the previous one
   finished, and no publication is taken while one runs.
4. **solve2D** (`start_solve2D_first`, `solve2D/`) classifies every selected row into
   `NCLS2D_FIRST` (150) classes with the pool's mask diameter and `cache=yes`, and solve2D's
   defaults otherwise. Its classes, class averages, halves and FRCs and the rows' 2D parameters
   become the stage's (`adopt_job_2D`); a result without them, or with other rows, counts as
   failed (`usable_job_2D`).
5. **The pool model's selection** (`select_first_set_classes`): the classes are scored with the
   pool quality model (`select_route_classes` in `simple_state_assignment`, `CAVG_QUALITY_MODEL_POOL_DEFAULT`), and the
   selection is mapped to their particles; the selected and rejected stacks go to
   `solve2D/quality_selection/`. **A rejected class's particles are deselected for the
   session:** they stay first-set rows, whose selection later publications do not change. When
   the model cannot score the classes or keeps none, every populated class is kept, with a
   warning.
6. **cls_expansion** (`start_clsexp`, `cls_expansion/`) splits the classes of the selected
   particles into `NCLS_EXPANDED` (500) subclasses in all, and never fewer than the classes, which
   cls_expansion refuses (`expansion_budget`; `ncls` is cls_expansion's total since 14 October
   2026). cls_expansion spreads the total over groups of similar class averages by population and
   keeps a class too small to split whole. Its sub-class averages are kernel pre-images of members
   (kernel-weighted, CTF-corrected means of the particles, the latent rank and bandwidth chosen by
   cross-validation; cls_expansion's only route since 14 October 2026), for the route and every
   publication alike; its settings are its defaults (`minpop` 25). The mask diameter is fitted to the
   particles' box. Its
   subclasses, their class averages, halves and FRCs, and the particles' subclass labels become
   the stage's (`adopt_job_2D`).
7. **solve3D_cavgs** (`start_cavgs3D`, `solve3D_cavgs/`) makes `nstates` volumes from those class
   averages (`pgrp=c1`, `prune=no`, the mask diameter; one run, no collapse restarts). Its
   low-pass limits are a fixed ramp from `LPSTART_CAVGS3D` (20 Å) to `LPSTOP_CAVGS3D` (12 Å) over
   its stages (`lpstart`, `lpstop`, `force_lp_range=yes`), not one derived from the class FRCs.
   Its angular sampling window is `solve3D_cavgs`' default, 10 degrees, as in particle `solve3D`
   (`solve3D_cavgs_policy.md`, section 2). Its result is not mapped onto the particles: they get
   their states from the assignments (section 10).
8. **The check** (`finish_cavgs3D`, `register_attempt`, `refs_usable`). `solve3D_cavgs` records
   the even/odd angular distance and state agreement of the class averages on its `vol_cavg`
   entries (`solve3D_cavgs_policy.md`, section 7). The angular distance checked is the median over
   the class averages whose halves are in one state (`eo_angdist_median`): the states are in
   frames of their own, so the mean over every class (`eo_angdist`, logged beside it, and checked
   only for a result without the median) mixes in poses compared across frames. The states
   that hold class averages are renumbered from 1 (`compact_states`), and their volumes, with the
   reprojection and orientation JPEGs beside them, copied into `state_refs/try_<n>/` with a
   record (`state_refs.txt`: the count, box and pixel size, the even/odd values and the largest
   correlation between two docked states, the file names). solve3D_cavgs refines its states
   independently, so each volume is in a frame and hand of its own: every pair of states is
   compared after docking the second onto the first, as it is and mirrored, over rotations and
   shifts, from 100 Å to 15 Å, the matching's limit (`docked_pair_cc`, `dock_vols`); the better
   of the two hands counts, and a correlation that is not finite counts as identical. The volumes
   pass when the angular distance is at most `EO_ANGDIST_MAX` (50 degrees since 14 October 2026,
   20 before; the first runs measured about 40 on volumes that looked usable), the agreement at least
   `EO_OVERLAP_MIN` (70%), no two docked states correlate above `REFS_PAIR_CC_MAX` (0.97: near-duplicate
   states would split the particles by noise), and at least `NSTATES_MIN` (2) states each hold
   `STATE_FRAC_MIN` (5%) of the class averages (provisional constants, until a validation run
   sets them). The run of 13 October 2026 that collapsed ended with the halves 90.6 degrees
   apart, as random orientations would be. The class averages each state was made from are
   written beside its volume, a stack and a tiled JPEG per state (`cavgs_state<NN>`,
   `write_state_cavgs`), when solve3D_cavgs' class averages are on disk. The GUI status is
   "checking the initial volumes".
9. **Passed:** the attempt's volumes, JPEGs (with the per-state class averages) and record are
   copied into `state_refs/` and are the
   references for the session (`accept_refs`, `restore_refs`): the stage's states
   (`params%nstates`) are theirs, and the GUI gets them (`send_refs`) under their states, without
   FSCs, populations or orientations, since they come from class averages. The stage's project
   does not register them.
10. **Failed** (`next_attempt`): the route runs again from its solve2D (a new random start), at
    most `ROUTE_TRIES` (2) attempts in all; the first set's rows the earlier attempt deselected
    stay deselected. Then the fallback: once the next publication is merged, `solve3D_cavgs` runs
    once on its pool class averages (`start_fallback`), checked the same way. When the fallback
    fails too, or the pool has published its last, the attempt with the smallest even/odd angular
    distance (`best_attempt`; one not measured last) gives the references, with a warning, rather
    than run with no states. Without any volumes there are none, and nothing is assigned.
11. **A failed step is skipped,** with a warning: a failed or unusable solve2D runs
    `solve3D_cavgs` on the first set's own class averages; a failed or unusable cls_expansion on
    the classes the model kept; a failed `solve3D_cavgs`, or one that leaves no class with a state
    or no volume for a state that holds class averages, is an attempt without volumes.
12. **The GUI** shows "running solve2D on the first set", "running cls_expansion on the first set"
    and "running solve3D_cavgs", each with the 3D progress of 1; with the references the progress
    is 2.
13. **Restart:** a restarted p07 takes the references recorded in `state_refs/` when every volume
    is there, and skips the route and the check (`init_params`); without them it runs the route on
    its first publication. A job folder a stopped step left is set aside (`fresh_job_dir`).
14. **Tests:** `test_sieve_publication` in p06's tester (the offsets, selection, class tables and
    populations; the combined files; a box mismatch and a missing half refused); in p07's,
    `test_rules` (the check, the best attempt, the renumbering, the sub-class budget),
    `test_reference_volumes` (an attempt that fails kept aside and the route run again, one that
    passes renumbered with its JPEGs and taken, a restart that takes it),
    `test_failed_attempts` (the fallback waiting for a publication, the best attempt after it or
    after the last publication, none without volumes), `test_send_refs` (the references sent
    under their states, without FSCs or populations, not registered) and `test_take_cavgs`; in the
    nightly `lib_stream` "solve 3D reference docking", `test_reference_docking` (two different
    structures distinct, a rotated and mirrored copy a near-duplicate). The fast tests compare no
    states (`l_dock_refs` off: the docking takes too long for them). The jobs and the pool
    model's selection of the first set run only in a stream.

## 10. The assignment of a publication

From the first publication after the first set (in p06's schedule, pool iteration 25), every
publication merged once the references exist is assigned to them: its classes are split into
sub-classes, whose averages are matched against the references, each particle takes the
state and pose of its sub-class, and the poses are refined locally with the states fixed.
Decided 13 October 2026
(`doc/refactoring_notes/completed/stream_p07_state_assignment_plan_2026-10-13.md`, decisions 1 to
12); the local refinement replaced the reconstruction from the class poses (reconstruct3D) on
14 October 2026.

1. **Due:** each publication merged while the references exist (`import_sets`), the first set of
   a restart that found them, and the final publication when it came before them
   (`accept_refs`). Only the newest waiting is assigned (section 4, item 2). Each assignment works
   in `assign/it_<n>/` (`fresh_job_dir`; `n` counts the session's assignments).
2. **The expansion** (`start_expand`, `assign/it_<n>/expand/`): cls_expansion of the
   publication's selected classes into `NCLS_EXPANDED` (500) sub-classes in all, never fewer than
   the classes (`expansion_budget`), with kernel sub-class averages, both as
   for the route, with the mask diameter fitted to the
   particles' box, on the 3D
   jobs' resources. Pool classes group particles mainly by view, so particles of different states
   with similar projections often share a class; the expansion separates them, and the matching
   compares sub-class averages with references made from sub-class averages. A usable result
   (`usable_job_2D`) is matched; otherwise the pool's class averages are, with a warning. The
   expansion is not adopted: the rows keep the pool's 2D parameters, which the next merge
   refreshes.
3. **The matching** (`start_assign`, `assign_cavgs_states`, a `simple_private_exec` program run in
   `assign/it_<n>/` in shared memory on the 3D jobs' threads), on a copy of the expansion's (or
   the stage's) project, with the references as `vol1..volN`, `lp=15` (`LP_ASSIGN`),
   `nspace=1000` (`NSPACE_ASSIGN`) and the mask diameter fitted to the class averages' box
   (`simple_cavgs_state_match`):
   - both sets are brought to one pixel size, the coarsest of lp/3 and the two sets' own (with
     boxes of 32 at least), then to one box by clipping, so references made from other class
     averages (the fallback's) match too;
   - each reference is reprojected over an even grid of directions on the whole sphere
     (`oris%spiral`, `reproject`); a mirrored view is the view from the opposite side, so the
     grid holds every mirror and none is searched (`match_imgs` with `l_nomirr`);
   - both sets are normalised within, and softly masked to, the particle, and every class
     average is aligned with every reprojection by a full in-plane search with shifts
     (`match_imgs`, the core of `match_cavgs` and `cluster_cavgs`);
   - per class: its state, the one whose best reprojection correlates best; its correlation; its
     margin, that correlation less the best of any other state; its pose, that reprojection's
     direction with euler(3) = 360 − the in-plane rotation (refine3D's convention), and its shift
     rotated to the reference frame, in the class averages' pixels (`reduce_state_matches`,
     `matched_pose`).

   The project's `cls3D` is written in place: per class its state (0 for a class `cls2D`
   deselects), Euler angles, shift, `corr`, `margin` and `class`.
4. **Rejection and flags** (`take_assignment`): over the candidates, the classes with a state and
   particles, a class whose correlation falls below the candidates' median less `FIT_NMAD` (3)
   median absolute deviations is rejected (`reject_by_fit`: state 0), the worst first, at most
   `FIT_REJECT_MAX` (30%) of them; a kept class with a margin below `MARGIN_AMBIGUOUS` (0.02) is
   flagged (`flag_ambiguous`: `cls2D` key `ambiguous` 1 in the matched project) and assigned as
   any other. The matched project is written back with both, and the matched (sub-)class
   averages are grouped into the assignment's folder (`write_assignment_cavgs`): a stack and a
   tiled JPEG per state (`cavgs_state<NN>`, the ambiguous ones included), one of those the fit
   rejected (`cavgs_rejected`) and one of the ambiguous (`cavgs_ambiguous`). Every JPEG p07
   writes is logged with its absolute path, then once more on a line of its own,
   `>>> JPEG <path>`, as the other stages log theirs. The log carries the fit's median,
   deviation, threshold and cap, each rejected class's correlation, population and state, the
   median population of the kept and rejected candidates, and the correlation of the best
   correlation with the log population over the candidates, so a validation run shows whether
   the rejection follows the classes' sizes (section 7).
5. **Mapping:** each kept class's pose and state are composed with its particles' 2D parameters
   (`map2ptcls` on the matched project), and each row selected before takes its particle's state,
   Euler angles, shift and correlation, in both particle segments: a particle follows its
   (sub-)class, so a later assignment that moves the class moves it. A row of a rejected class,
   or without a class of the matching, is deselected until the next publication, and recorded
   (`fit_desel`) so that a first-set row is selected again by that publication's merge (section
   4, item 5) and judged again. The log counts each state's particles, the
   rejected and ambiguous classes, and the particles rejected by the fit and the pool model, by
   the model alone (not applied) and by the fit alone.
6. **The local refinement** (`start_refine`, `assign/it_<n>/refine3D/`): one iteration
   (`maxits=1`) of `refine3D refine=cont` over every selected row from its mapped pose
   (`pgrp=c1`, `objfun=cc`, the pool's mask diameter, the stage's `nspace`, the 3D jobs'
   resources, distributed even with one part), each pose refined continuously against its own
   state's reference within refine=cont's bounds (`athres_cont`, `trs`). refine3D reconstructs
   its starting volumes from the mapped poses (it clears the update counts at its first
   iteration, so every row with a state is taken), sets its low-pass limits from the FSC, and
   ends the iteration with the states reconstructed from the refined poses and postprocessed.
   The states never change: refine=cont searches none, and a multi-state start seeds no labels
   for it. One iteration bounds the job's cost per publication (`refine=cont` would otherwise
   run until its convergence rule, up to `maxits=100`).
   - **The seeds** (`prepare_refine_seeds`): the half-sets of the project's partition
     (`eo_transaction`: one half per micrograph, recorded on the stage's stack rows, so the same
     for every assignment; refine3D stops on halves that disagree with that record) and the
     projection direction closest to the mapped pose in the job's grid (`set_projs`), which
     refine=cont requires of a seed. A given `nspace` must be
     even (the reference spiral's); the commander refuses an odd one.
   - **The states:** the job's project holds the states the particles hold, renumbered from 1
     (`compact_states`), since refine3D reconstructs no volume for an empty state;
     `take_refine_result` maps them back, so a reference whose sub-classes were all rejected has
     no volume until a later assignment gives it particles.
   - When done, its volumes, FSCs and poses become the stage's (section 4, item 8), the stage's
     rows are written as `assign/it_<n>/assignment.simple`, the snapshots' source, beside
     refine3D's own project, which is kept as the job left it, and older assignment folders go.
7. **The GUI** shows "expanding the classes", "assigning class averages" and "refining the poses
   locally" while the jobs run, and the refinement's states when done.
8. **Tests:** `test_assignment` (sub-classes of one pool class in two states, a badly fitting
   sub-class rejected and its first-set row selected again by the next publication and judged
   again, the margin flagged, a later matching moving particles, a matching of other rows
   refused),
   `test_local_refinement` (the assigned rows valid refine=cont seeds; a result taken under the
   stage's states with an empty state left out, with the refined poses; a result of other rows
   refused); `unit_reconstruction` "state assignment rules": `test_state_cavgs` (a stack and a JPEG
   per state, of the rejected and of the ambiguous; none for a state without class averages) and
   `test_fit_rules`; and `unit_pftc_align2D3D` "cavgs state match" (`reduce_state_matches` and
   `matched_pose` on fixed numbers; projections of two volumes of Gaussian blobs, each matched to
   its own with a positive margin at its own direction, in-plane angle and no shift, and a
   mirrored one at the opposite view).
