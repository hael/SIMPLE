# Stream multistate 3D by class-average assignment, plan, 13 October 2026

**Status:** implemented on 13 October 2026, not yet built or run; the living description is
`doc/policies/stream/stream_3D_ingestion_policy.md`, sections 4, 9 and 10. Where the code differs
from this plan:
- the matching brings the class averages and the references to one pixel size and box (the
  coarsest of lp/3 and the two sets' own), so references made from pool class averages (the
  fallback) match sub-class averages of another sampling;
- the pool model's verdict is compared with the fit's by particles (both, the model alone, the
  fit alone), since the fit rejects sub-classes and the model pool classes;
- the matching's shift limit is `assign_cavgs_states`' default (10 pixels), not the pool's;
- the pool's final publication, when taken before the references exist, is assigned once they do,
  and it ends the wait for a fallback publication (the best attempt is taken instead);
- on 14 October 2026 the reconstruction from the class poses (decision 3) was replaced by the
  local refinement of section 11, as the execution path: `refine3D refine=cont`, states fixed,
  from the mapped poses until refine=cont converges, on the populated states renumbered in the
  job's project; refine3D now seeds no state labels for refine=cont;
- after the review of 14 October 2026 (`stream_p07_review_2026-10-14.md`): the refinement runs
  one iteration (`maxits=1`); the check also refuses near-duplicate states (two states
  correlating above 0.97 at 15 Å after docking, in both hands, since the states are refined in
  frames of their own); the even/odd angular distance checked is the median over the class
  averages whose halves are in one state (solve3D_cavgs logs it beside its original mean over
  every class), at most 50 degrees (the first runs measured about 40 on usable volumes, on the
  original mean); the fit's deselection lasts until the next publication for
  every row, the first set's included; each assignment logs the fit's statistics with the
  rejected classes' populations; the snapshots' source is `assign/it_<n>/assignment.simple`.

p07 stops refining particles. It makes the states' volumes once, from class averages (the route
`ini3D_cavgs` already runs: solve2D, the pool model's selection, cls_expansion, `solve3D_cavgs`), and
checks them before use. From the pool's publication at iteration 25 on, each publication's classes
are split by cls_expansion into about 500 sub-classes, whose averages are matched against those
volumes: a sub-class that fits no state well is rejected, every other goes to the state it fits best,
and its particles take that state. The states' volumes are then
rebuilt from the particles, without refinement. The `solve3D` and `solve3D_addon` runs, and the final
multistate `refine3D`, are retired.

## Decisions

First round:

1. **Later publications are matched again.** Every publication from iteration 25 on is matched
   against the volumes. The pool's 150 classes keep their numbers, but their contents drift, so a
   map made once would age.
2. **A particle follows its (sub-)class.** At every assignment a particle takes its current
   sub-class's state (decision 12). A particle that moves to another class, or whose sub-class is
   matched to another state, changes state.
3. **The states' volumes are rebuilt from the particles** after every assignment: each class's 3D pose
   and state are mapped onto its particles, and each state is reconstructed from them, without
   refinement. This gives particle-level volumes, FSCs and resolutions for the GUI and the snapshots.
4. **No final multistate refine3D.** The session ends with the last assignment and its reconstruction.
5. **The references are always the checked initial volumes from class averages.** Later
   publications are matched against them, a fixed reference set, so the assignment cannot drift or
   collapse through feedback (assign, rebuild, assign). The particle reconstructions are for display,
   FSCs and resolutions only.
6. **Without usable initial volumes** (the route failed, or its volumes failed the check of decision
   9 on every try), `solve3D_cavgs` runs on the first later publication's pool class averages, and its
   volumes, checked the same way, are the references.
7. **Matching at a fixed 15 Å** (a named constant): every state at the same limit, so none wins on
   sharpness or smoothness, coarse enough for noisy class averages.

Second round:

8. **The matching is a program of `simple_private_exec`,** `assign_cavgs_states`, which p07 runs as a
   job on the 3D jobs' threads. It works in memory with existing routines: each state volume
   reprojected over a grid of directions (`reproject`, `simple_volinterp`), the class averages
   prepared as `cluster_cavgs` prepares them, and every reprojection aligned with every class average
   by a full in-plane search (`match_imgs`, `simple_strategy2D_utils`, the core of `match_cavgs` and
   `cluster_cavgs`).
9. **The initial volumes are checked before use,** and the route runs again when they fail. The
   checks use what `solve3D_cavgs` already measures: the even/odd angular distance of the class
   averages, the even/odd state agreement, and each state's class averages. The run of 13 October
   2026 that collapsed ended with the halves 90.6 degrees apart, as random orientations would be.
10. **Publications are rejected by fit, not by the pool model.** A class whose best correlation
    against every state is an outlier of its publication is rejected. The pool model stays in the
    route (its classes are selected before any volume exists), and at first it is also run on each
    publication and logged beside the fit's decision, without being applied, so the two can be
    compared before the model is dropped there.
11. **The margin flags ambiguous sub-classes:** the difference between a sub-class's best and
    second-best state. An ambiguous one is assigned to its best state as any other, and flagged and
    counted in the log, so a later change can deselect them instead.

Third round:

12. **Every publication's classes are expanded before matching.** cls_expansion splits each class of
    the publication into the same number of sub-classes, about `NCLS_EXPANDED` (500) in all
    (`subclasses_per_class`), as the route does for the first set, and the sub-class averages are
    matched, rejected and flagged instead of the pool's class averages. Pool classes group particles
    mainly by view, so particles of different states with similar projections often share a class;
    the expansion separates them, and the matching then compares sub-class averages with references
    made from sub-class averages. Each publication's assignment is then three jobs: cls_expansion,
    the matching, the reconstruction.

Defaults taken without a question (change any before implementing):
- **The route always runs.** `ini3D_cavgs` goes: without `solve3D` there is no other way to make volumes.
- **The states** are those of solve3D_cavgs that hold class averages, renumbered from 1, so `nstates`
  (3) is an upper bound. The first set's particles are no longer mapped onto the route's subclasses:
  they get their state at the first assignment, like every other particle.
- **The first assignment** is made on the first publication after the first set (iteration 25), or on
  the newest one once the references exist, if they come later.
- **The first set** keeps its session-long selection (a first-set row the route deselected stays
  deselected, and so does one whose class the fit rejects).
- **Provisional constants,** each a named constant until a validation run sets it:
  - the matching: 15 Å, 1,000 directions per state (about 6 degrees apart), the pool's shift limit;
  - the check: even/odd angular distance at most 20 degrees, even/odd state agreement at least 70%,
    at least 2 states each holding at least 5% of the class averages; the route tried at most twice;
  - the rejection: a best correlation below the publication's median minus 3 times the median
    absolute deviation, at most 30% of a publication's sub-classes rejected by fit (the worst first);
  - the margin: ambiguous below 0.02 in correlation;
  - the model logged beside the fit for every publication until the constant that runs it is turned
    off.

## 1. Where things stand

- p07 takes the first publication as the first set and, with `ini3D_cavgs=yes` (the default), runs the
  route to initial volumes; `take_cavgs3D_result` maps the classes' poses and states onto the first
  set, and `solve3D` starts from them (`cavg_ini_ext`).
- Every later publication is scored with the pool model (`select_cavgs`), merged into the rows
  (`merge_publication`), and grows the cohort of the next `solve3D_addon` run; a run's verdict can roll
  it back. The pool's final publication triggers one multistate `refine3D`.
- The GUI gets each run's volumes, FSCs, populations and orientation distributions (`send_volumes`),
  and the initial volumes as soon as `solve3D_cavgs` is done (`send_cavgs3D_volumes`). 3D snapshots
  are cut from the latest run's project (`result_projfile`).
- `solve3D_cavgs` computes the even/odd angular distance (`conv_eo`) and state overlap
  (`conv_eo_states`) of its last stage, but only logs them.
- `simple_private_exec` runs internal programs: a spec in `make_private_ui`
  (`simple_private_prgs`), the commander imported through `simple_private_exec_api`, and a case in
  `simple_private_exec_driver`. p07 can start a job with it (`exec_bin`).

## 2. The new flow in p07

1. **The first publication** is taken as the first set, as now (section 4, item 3 of the ingestion
   policy).
2. **The route** runs as now: solve2D, the pool model's selection, cls_expansion, `solve3D_cavgs`. Its
   result is no longer mapped onto the particles.
3. **The check** (decision 9, section 4) reads the result's even/odd diagnostics and states. Passed:
   the volumes of the states that hold class averages are copied into p07's `state_refs/`, renumbered
   from 1, and registered as the reference set; the GUI gets them as now (without FSCs). Failed: the
   route runs again, at most twice in all, then the fallback (decision 6).
4. **Each later publication,** once the references exist:
   - merged into the rows, as now;
   - the pool model's verdict scored and logged (decision 10), not applied;
   - **expanded** (cls_expansion, in `assign/it_<n>/expand/`): every class of the publication split
     into `subclasses_per_class(ncls)` sub-classes, with the mask diameter fitted to the particles'
     box, on the 3D jobs' resources; the sub-classes, their averages (with halves) and the selected
     particles' sub-class labels and 2D parameters are the job's result;
   - **matched** (`assign_cavgs_states`, in `assign/it_<n>/`): every sub-class average against the
     references; the result is each sub-class's best state, its 3D pose, its best correlation and its
     margin (section 3);
   - **rejected and flagged** (section 5): sub-classes whose best correlation is an outlier are
     rejected, their particles deselected; ambiguous sub-classes flagged;
   - **mapped:** each kept sub-class's pose and state onto its particles (`map2ptcls` on the
     expansion's labels), overwriting their previous states (decision 2); the stage's rows keep the
     pool's 2D parameters, which the next publication's merge refreshes;
   - a failed or unusable expansion is skipped, with a warning: the pool's class averages are
     matched instead, as without decision 12;
   - **reconstructed** (`reconstruct3D`, in the same folder): each state from its particles at the
     mapped poses; the volumes, FSCs, populations and orientation distributions go to the GUI
     (`send_volumes`), and the folder's project becomes the snapshots' source (`result_projfile`).
5. **No publication is taken while a job runs,** as now; only the newest waiting one is taken.
6. **The pool's final publication** is assigned and reconstructed like any other; nothing follows it.

## 3. The matching program, `assign_cavgs_states` (simple_private_exec)

- **Inputs:** a project with the averages to match (`cavg` entry and `cls2D`: the expansion's
  sub-class averages, or the pool's class averages without them), the reference volumes
  `vol1..volN`, `mskdiam`, `lp` (15 Å), `trs`, `nspace` (1,000), `nthr`.
- **Steps:**
  1. each reference volume reprojected over an even grid of directions (`oris%spiral`, `reproject`),
     at the class averages' sampling;
  2. the class averages read, masked and normalised as `cluster_cavgs` prepares them;
  3. every reprojection aligned with every class average (`match_imgs` at `lp`, with `trs`): a
     correlation, in-plane angle, shift and mirror flag per pair;
  4. for each class: its best correlation per state; its state, the best; its margin, the best minus
     the second best (0 with one state); its 3D pose, the best reprojection's direction composed
     with the in-plane angle and shift, the mirror converted to the matching 3D pose.
- **Output:** the project's `cls3D`, written in place: per class the state, Euler angles, shift,
  correlation and margin.
- **Registration:** the commander beside the class-average commanders, imported through
  `simple_private_exec_api`, its keys in `make_private_ui`, its case in
  `simple_private_exec_driver`; no public UI entry.
- **Cost:** averages × states × directions aligned pairs: 500 × 3 × 1,000 = 1.5 million at 15 Å,
  a few `cluster_cavgs` runs' worth. It may need fewer directions, or nparts, once measured.

## 4. The check of the initial volumes

- `solve3D_cavgs` writes the even/odd angular distance and state overlap it already computes into its
  project's out segment, beside the `vol_cavg` entries, so p07 can read them.
- The volumes pass when (provisional constants above):
  - the even/odd angular distance is at most 20 degrees;
  - the even/odd state agreement is at least 70%;
  - at least 2 states each hold at least 5% of the class averages.
- A pure function decides it from those values, so the rule is unit-tested.
- Failed: the route runs again from its solve2D (a new random start), at most twice in all; then the
  fallback of decision 6, checked the same way. When every attempt fails, p07 takes the attempt with
  the smallest even/odd angular distance, with a warning, rather than run with no states.

## 5. Rejection and ambiguity

- **Fit-based rejection** (decision 10), in p07 from the program's output: over the publication's
  sub-classes, one whose best correlation falls below the median minus 3 times the median absolute
  deviation is rejected (state 0, its particles deselected), at most 30% of them, the worst first.
  A pure function decides it, unit-tested.
- **The model beside it:** the pool model scores the publication as `select_cavgs` does now, writes
  its selected and rejected stacks to the publication's quality folder, and is not applied. The log
  counts the classes both reject, the model alone, and the fit alone.
- **Ambiguity** (decision 11): a kept sub-class with a margin below 0.02 is flagged (`cls2D` key
  `ambiguous` in the expansion's project) and counted in the log; it is assigned as any other.

## 6. What goes from p07

- `start_solve3D`, `start_addon`, `start_final_run`, `finish_run`'s run handling, `roll_back_addon`,
  `read_addon_verdict`, the addon verdict in the status;
- the first run's cap and queue (`cap_first_run`, `cap_addon_cohort`, `cap_rows`,
  `set_queued_states`, `release_queue`, `queued`), the cohort rules (`count_cohort`, `count_frozen`,
  `retry_cohort`, `ncohort_refused`, `n_rollbacks`, `ADDON_COHORT_FRAC`), the frozen base
  (`frozen_projfile`, `frozen_active`), `nptcls_at_full`, `l_final_pending`;
- `take_cavgs3D_result`'s particle mapping and `l_cavg_ini_ext` (the route keeps `compact_states` for
  the references' renumbering);
- the pool model's selection of publications as the decision (`select_cavgs` becomes log-only);
- the inputs `nptcls3D_max`, `nptcls_addon_max`, `nstages` (p07's solve3D stages) and `ini3D_cavgs`,
  from the parameters, the parser, the UI and the commander's defaults (`NPTCLS3D_MAX`,
  `NPTCLS_ADDON_MAX`, `NSTAGES3D`);
- the folders `solve3D/`, `solve3D_addon/`, `refine3D_final/` and their clean-up; `next_job` is
  rewritten for the new jobs (the route's steps, the matching, the reconstruction).

Kept: the publication import and merge, `take_cavgs`, the route (its cls_expansion step serves the
publications too), the snapshots (from the latest assignment's project), the GUI messages.

## 7. The GUI

- The stage's status: "running solve2D on the first set", "running cls_expansion on the first set",
  "running solve3D_cavgs", "checking the initial volumes", "expanding the classes", "assigning
  class averages", "reconstructing the states", "idle".
- The states: each state's population (the particles assigned) and resolution (FSC 0.143 of its
  reconstruction); "particles at last refine" becomes the particles of the latest assignment.
- The volumes: the references until the first reconstruction, then each assignment's reconstruction.
- The p07 NICE page needs no new fields; its labels may change ("assigned" for "refined").

## 8. Restart

- The references are kept in `state_refs/` with a record of their states and their check. A restarted
  p07 with them takes them and skips the route and the check; without them it runs the route again on
  the newest publication (taken as the first set).
- The newest publication is then assigned and reconstructed again. Job folders left unfinished are
  set aside (`fresh_job_dir`).

## 9. Docs

- `stream_3D_ingestion_policy.md`: section 4 (runs) rewritten for the assignment; section 9 (the
  route's result checked and kept as references, not mapped); the rejection by fit; the fallback; the
  scope; tests.
- `restart_policy.md`: p07's row.
- `solve3D_cavgs_policy.md`: the even/odd diagnostics written to the project.
- A short policy for `assign_cavgs_states`, or a section of the ingestion policy.
- The `simple-main-stream` skill's p07 line, with approval.

## 10. Tests

- `assign_cavgs_states`: class averages projected from two distinct small volumes, each assigned to
  the volume it was projected from with a positive margin; a mirrored projection given the matching
  3D pose. A library test (it reprojects and aligns), or `unit_reconstruction` if small enough.
- Pure rules, in `unit_stream` "solve 3D": the check (each threshold, the two-state minimum); the
  fit-based rejection (the outlier cut, the 30% cap, the worst first, none rejected when all fit);
  the ambiguity flag; the composition of a matched direction, in-plane angle, shift and mirror into
  a 3D pose.
- p07 (`unit_stream` "solve 3D"):
  - the references registered from a fixture `solve3D_cavgs` result that passes the check,
    renumbered, states without class averages left out; a failing one leading to a second route;
  - an assignment from a fixture `cls3D` mapped onto the rows; a later assignment that moves a class
    to another state moves its particles (decision 2); a rejected class's particles deselected;
  - the job rules: the expansion due once the references exist and a later publication is taken,
    the matching after it, the reconstruction after that, no job before the references; a failed
    expansion leading to matching the pool's class averages;
  - the mapping by sub-class: particles of one pool class split across two sub-classes matched to
    different states take different states;
  - the fallback: no usable references after the route, `solve3D_cavgs` on the pool's class averages;
  - the restart: references on disk taken, the route skipped.

## 11. Later, not in this plan

- The particle-vote purity of each class (its particles' states before the assignment) as a check.
- Ambiguous classes deselected instead of assigned, once the margin's constant is validated.
- The pool model dropped from publications once the logged comparison supports it.
- References refreshed from a new `solve3D_cavgs` on a larger publication.
- A bounded particle polish (one local refinement per assignment, states fixed) if the
  reconstructions from class poses are too blurred. Done on 14 October 2026 (status note above).

## 12. Order

1. `solve3D_cavgs` writing its even/odd diagnostics to its project; `assign_cavgs_states` in
   `simple_private_exec`, with its tests.
2. p07: the check and the references (registration, renumbering, `state_refs/`, the second route,
   the fallback); tests.
3. p07: the expansion job per publication, the matching job, the rejection by fit with the model
   logged beside it, the ambiguity flag, the mapping by sub-class, the reconstruction job, the GUI,
   the restart; tests.
4. p07: the retirement of the solve3D/addon/final paths and their inputs; tests updated.
5. The docs and the skill line.

## 13. Risks, and what lowers them

- **The assignment is only as good as the references.** The check (decision 9) refuses volumes whose
  halves disagree or whose states are nearly empty, and the route runs again; `prob_athres_taper`
  (already on) lowers the chance of a collapse. Near-duplicate states are not yet detected.
- **Particle poses come from the class averages,** composed with each particle's 2D in-plane
  parameters, so the reconstructions are limited by the class averages' alignment and the matching
  grid (about 6 degrees). A finer grid, or a bounded particle polish (section 11), would help.
- **Mixed classes** are split by the expansion (decision 12) before matching; a sub-class that is
  still mixed goes to one state, and the margin flags it (decision 11).
- **The expansion is the heaviest step:** it reads every particle of the publication, much like a
  pool iteration. When a publication's three jobs take longer than a pool iteration, p07 takes only
  the newest publication waiting, so the assignment lags the pool rather than queueing behind it.
- **Sub-class averages hold about a third of a class's particles,** so they are noisier, and the
  match at 15 Å is less certain than for whole classes.
- **Fit-based rejection discards whatever the references lack,** such as a conformation or views the
  first set did not hold, and blurred junk may fit low-pass projections well enough to be kept. The
  model's verdict is logged beside it until the comparison settles which to keep.
- **Matching at 15 Å** may not separate states that differ only at higher resolution.
- **The pool's class numbers persist but their contents change;** re-matching every publication
  (decision 1) follows that, at the cost of a job per publication.
