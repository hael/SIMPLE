# Stream 3D Ingestion Policy

How the 2D pool hands particles to multistate 3D, how 3D takes them in, and what a change must
preserve.

## 1. Scope

- The pool's publications (p06): `export_pool_state` in
  `src/main/stream/stages/simple_stream_stage_pool2D.f90`, and `build_pool_publication`,
  `publish_pool_state` and `delete_pool_publication` in
  `src/main/stream/pool2D/simple_stream_refine2D_utils.f90`.
- Their import (p07): `import_sets`, `select_cavgs` and `merge_publication` in
  `src/main/stream/stages/simple_stream_stage_solve3D.f90`.
- p07's runs: when it starts `solve3D` and `solve3D_addon` (`next_job`), what it does with a
  failed run and an addon verdict, and the folders it keeps. The addon program itself is
  governed by `doc/policies/3D/solve3D_addon_policy.md`; this policy relies on its row contract
  (section 4 there) and its report (section 11 there).

## 2. What the pool publishes

The pool publishes its state, not deltas. A publication is the pool's classified set as one
completed iteration left it.

1. **The stacks whose particles have been through an iteration** (`updatecnt > 0`), with their
   micrographs and all their particles. Stacks never classified, such as the sets imported since
   the last dispatch, are never published. Within a published stack, a particle never updated (a
   fractional update did not sample it) keeps its row but is published deselected (state 0): its
   class is the random one it was given on import. A later publication selects it once updated.
2. **2D parameters at publication time:** class, in-plane angle, shifts and state. `ptcl3D` is a
   copy of `ptcl2D`.
3. **Self-consistent indexing.** Stacks are renumbered from 1, and their particle ranges and the
   particles' stack indices follow. Each particle carries its image index in its stack
   (`indstk`).
4. **Classes:** the pool's `cls2D` table, all classes.
5. **Class averages and FRCs.** The pool's class averages (and their even and odd halves) and
   FRCs at the native sampling: rescaled when the pool runs downscaled, copied otherwise. They are
   registered in `os_out` with the pool's mask diameter.

## 3. When and how it is published

1. **When:** once per completed pool iteration from `EXPORT_START_ITER` (25) on, in the pass where
   that iteration's results have come back and before the next iteration is dispatched. The
   publication is therefore one consistent iteration.
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
3. **One quality decision per publication.** The publication's class averages are scored once
   with the pool quality model (`CAVG_QUALITY_MODEL_POOL_DEFAULT`), with a mask diameter fitted
   to their box (`fit_mskdiam`). The selection is mapped to the publication's particles by class
   (`map_cavgs_selection`). The selected and rejected class averages are written with JPEGs in
   `quality_selection/<id>/`.
4. **Merge** into the stage's rows (`merge_publication`):
   - **A stack the stage holds** is matched by its stack name, in whatever order the publication
     lists it, and each particle by its image index in the stack (`indstk`; the stage's row must
     record the same image, or the stage stops). Its particles take the publication's 2D
     parameters (`transfer_2Dparams`) and 2D state in place. The 3D state is a run's multistate label: a particle the publication
     deselects gets 0, a selected particle keeps its label, and a selected particle without one
     (state 0) gets 1. They keep their 3D parameters, CTF and optics group.
   - **A new stack** is appended with its micrograph and particles (2D and 3D), after the rows
     the stage holds.
   - **A stack the publication lacks** keeps its rows, deselected (state 0 in both segments). This
     happens after a pool restart, until the restarted pool has classified it again.
   - **The stage's `cls2D`** is the publication's.
5. **The first run:** `solve3D` starts once at least `MIN_PTCLS_PER_STATE` (5) particles per
   state are selected (`next_job`). Its settings come from the command line, defaulting to
   `NSTATES3D` (3), `NSTAGES3D` (5), `LPSTART3D` (50), `LPSTOP3D` (10), `NPARTS3D` (8) and
   `NTHR3D` (8), with the pool's mask diameter. Publications start past the pool's low-pass ramp
   (iteration 20), so the first is already a settled classification.
6. **Addon runs:** once a result exists, its active rows are the frozen particles, and the
   cohort is the selected rows that are not (appended since, or selected again). An
   `solve3D_addon` run starts when the cohort reaches max(`MIN_PTCLS_PER_STATE` × nstates,
   `ADDON_COHORT_FRAC` (10%) of the frozen particles). The first bound is the addon's own floor
   (its refusals of an empty or small cohort cannot happen); the second bounds the total cost,
   since every run accumulates all frozen particles.
7. **Failures:** a failed `solve3D` stops the stage (`THROW_HARD`). A failed addon run (refused,
   for example for a state without frozen particles, or crashed) leaves the latest result the
   base, and the next run waits for a cohort larger than the one that failed.
8. **The verdict:** after each addon run its report (`solve3D_addon_report.txt`) is read; the
   verdict per state (IMPROVED, UNCHANGED, REGRESSED) is logged and shown in the GUI status. A
   REGRESSED state is warned about and the result stays the base: rejecting it is the user's
   call (addon policy, section 11).
9. **Particles are aligned once.** A particle's pose and state, once frozen, are never revisited,
   and later runs inherit the first `solve3D`'s ladder. There is no rebase (a fresh `solve3D` or
   a `refine3D` once the data have grown); that is a known gap (section 7).
10. **After each run** the stage's project is written (under a temporary name, renamed). The GUI
    gets each state's volume, resolution, FSC curve, reprojections and orientation distribution.
11. **What p07 keeps:** the newest `NQUALITY_KEPT` (3) quality folders and those of the
    publications a run started from (listed in `quality_selection/runs.txt`); the latest addon
    iteration folder (the frozen base), its predecessors removed once it has completed; no folder
    set aside for a job left unfinished once a later run has completed. A stop cancels the
    running job (`doc/policies/stream/restart_policy.md`).

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

## 6. Change rules

- A change of what a publication holds, or when it is written, updates sections 2 and 3, and
  checks p07's merge against it.
- A change of the merge keeps the row invariants of section 5 and the addon policy's section 4.
- Tests:
  - `unit_stream` "pool 2D": `test_publication_holds_classified_stacks` (classified stacks only,
    renumbering, image indices, nothing to publish before any iteration) and
    `test_export_numbering`;
  - `unit_stream` "solve 3D": `test_merge_publications` (appending, matching in another order
    and by image index, selection and class in place with the 3D parameters and multistate label
    kept, deselection of a missing stack, the classes), the watch order and mask diameter, the job
    rule (first run, cohort thresholds, a failed run), the cohort count, the retention, the
    volumes sent.
  - The class-average selection of a publication has no unit test.

## 7. Known gaps

- **Cost of a publication:** each one writes the pool's whole classified set, about the size of
  the pool project, once per iteration past 25.
- **Stack matching** is a search from the last match; it is linear per stack only when the
  publication's order differs from the stage's (after a pool restart).
- **p07 never removes rows:** a stack the pool drops for good stays deselected in p07's rows.
- **Retention of 2:** a reader that lags two publications behind finds the files of the older
  one gone. p07 reads the newest in the pass it sees it.
- **Not yet decided:** a stability trigger other than "past iteration 25" for the first
  `solve3D`.
- **No rebase:** frozen particles are never realigned, and the first run's ladder is inherited
  (section 4, item 9).
- **The addon run is not tested end to end in the stage:** its command line, frozen project and
  report reading are covered by the addon's own tests, not by p07's.
