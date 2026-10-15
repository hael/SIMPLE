# solve3D_pyramid: batch multistate 3D by class-average assignment, plan, 15 October 2026

A `simple_exec` program that does offline, in one pass, what stream p07 does session-long: it
makes the states' reference volumes once from class averages of a first subset of the particles,
assigns every particle to a state by matching sub-class averages against those references, and
refines the poses locally with the states fixed. It has none of p07's stream machinery (watching,
publications, retries, fallbacks, GUI): every failed step is a hard exit.

Status (15 October 2026): implemented, not yet compiled or run. The contract is
`doc/policies/3D/solve3D_pyramid_policy.md`. Where the implementation departs from this plan:
- The result is written over the execution directory's project, which keeps the input's name
  (`mkdir=yes` copies it), as `solve3D_chunks` publishes. There is no separate
  `solve3D_pyramid.simple` (sections 2.12 and 5).
- The report is rewritten after every step and on a failure, with a status line, so a hard exit
  leaves the numbers so far.
- The phantom test takes the whole phantom as the first set, not `nptcls_refs=1200`. The
  phantom is one stack ordered by direction, so the first third would be a cap of 20 of the 60
  directions (section 6).
- `solve3D_pyramid_blobs` is its own CTest entry (highlevel, 7200 s, 8 processors), and
  `SIMPLE_CTEST_BUDGET` was raised from 38 to 39 for it. This was the owner's decision on
  15 October 2026, chosen over folding it into `flex_pca_blobs`, whose wall time it would double.
- The docking test stays in the nightly `lib_stream` "solve 3D reference docking", through p07's
  `register_attempt`, now a thin call to the shared `check_state_refs` (open question 6).
- p07's `register_attempt`, `restore_refs`, `take_assignment`, `write_assignment_cavgs`,
  `prepare_refine_seeds` and `take_refine_result` stay as stage wrappers of the shared
  procedures (`check_state_refs`, `read_refs_record`, `map_assignment`, `write_matched_cavgs`,
  `set_refine_seeds`, `map_refine_result`), so the stage tests keep their names. Its
  `score_with_model`, `select_with_model`, `docked_pair_cc` and the `nopass` rules went.

Decided 15 October 2026 (the four questions asked):
1. **One pass over everything.** The references come from a first subset; then one expansion,
   matching and refinement of all particles. No growing levels.
2. **The input is 2D-classified**: `ptcl2D` holds a refine2D/solve2D result and the class
   averages (with halves and FRCs) are registered in `os_out`, as a pool publication is.
3. **The references are made by p07's full route**: solve2D into 150 classes on the first
   subset, the pool model's selection, cls_expansion to 500, solve3D_cavgs; checked, hard exit on
   failure.
4. **The pass ends with one `refine3D refine=cont` iteration** over every particle.

Decided 15 October 2026 (the second round):
5. **`nptcls_refs` is a count**, default 100,000 (`NPTCLS_REFS`): the first 100,000 selected
   particles, whole stacks in arrival order, as p07's first set (`NPTCLS_FIRST3D`); a project
   with fewer uses them all. A new input (parameter, parser, UI).
6. **The pool model is logged, not applied, on the full set**, as in p07: the fit's rejection is
   the only rejection there, and the model's verdict is logged beside it.
7. **The report is a text table** in the style of `solve3D_chunks_report.txt`.
8. **The name stays `solve3D_pyramid`.**
9. **The matching is the `assign_cavgs_states` program, executed in-process** (its commander on
   a child command line), as every other step: one code path for the matching, shared with p07,
   at the cost of writing and reading a 500-class project.
10. **The pool model's selection of the first subset's classes is optional**: `use_model=yes|no`
    (the existing sieve input, default yes). With `no`, every populated class of the first
    subset's solve2D goes into the expansion and the references; the particles of none are
    deselected. The model's scoring of the full set (decision 6) is unaffected: logged either way.

Defaults taken without a question (change any before implementing):
- **Every rule and constant is p07's**, shared, not copied (section 4): 150 classes, 500
  sub-classes, solve3D_cavgs' 20 → 12 Å ramp with the angular taper, the check (same-state
  median ≤ 50°, agreement ≥ 70%, docked pair correlation ≤ 0.97, ≥ 2 states with ≥ 5%), the
  matching at 15 Å over 1,000 directions, the fit rejection (median − 3 MAD, ≤ 30%), the
  ambiguity margin 0.02, `maxits=1`, `objfun=cc`, `pgrp=c1`.
- **The pool model** selects the first subset's classes for the route, as in p07, unless
  `use_model=no` (decision 10); on the full set it scores the classes and is logged beside the
  fit (decision 6), not applied.
- **The mask diameter** comes from the project's class-average entry (`get_mskdiam('cavg')`), or
  `mskdiam` on the command line, which wins; the route's jobs fit it to the box as p07 does.
- **Steps run in-process**, each in a subdirectory of the execution directory, as
  `solve3D_chunks` runs `solve3D` and `solve3D_addon` (commander `execute` on a child command
  line, `simple_chdir` back), not as queued asynchronous jobs. Distribution comes from `nparts`
  on the children (`nparts` makes `solve2D`, `cls_expansion`, `solve3D_cavgs` and `refine3D`
  distributed).

## 1. Where things stand

- p07 (`simple_stream_stage_solve3D`) holds the whole method, but entangled with the stage:
  the route's steps (`start_solve2D_first`, `select_first_set_classes`, `start_clsexp`,
  `start_cavgs3D`), the check (`register_attempt`, `refs_usable`, `docked_pair_cc`), the
  assignment (`start_expand`, `start_assign`, `take_assignment`, `reject_by_fit`,
  `flag_ambiguous`), the seeds and the refinement (`prepare_refine_seeds`, `start_refine`,
  `take_refine_result`), the per-state class-average stacks (`write_state_cavgs`, `log_jpeg`) and
  the small rules (`compact_states`, `expansion_budget`, `count_active_classes`, `fit_mskdiam`,
  `median_of`). Its stateless rules are `nopass` type-bound procedures; the rest read and write
  the stage.
- `assign_cavgs_states` (`commander_assign_cavgs_states`, `simple_cavgs_state_match`) is already
  a program and can be executed in-process.
- `solve3D_chunks` (`simple_commanders_solve3D_chunks`, `simple_solve3D_chunks_plan`) is the
  model for a batch workflow that runs child commanders in subdirectories, cuts whole stacks in
  arrival order (`partition_selected_balanced`), and reports.

## 2. The program

```text
simple_exec prg=solve3D_pyramid projfile=<2D-classified project> nstates=3 [nptcls_refs=100000]
            [use_model=yes|no] [mskdiam=<A>] [nparts=<n>] nthr=<n>
```

Execution directory `N_solve3D_pyramid/` (`mkdir=yes`, refused otherwise, as solve3D_chunks),
the input project copied into it. Steps, each a subdirectory, each a hard exit on failure:

1. **`refs_subset/`:** the first `nptcls_refs` selected particles, whole stacks in arrival
   order, selected in a copy of the project (the others deselected). Fewer than
   `MIN_PTCLS_PER_STATE × nstates` selected particles: exit.
2. **`refs_subset/solve2D/`:** solve2D into 150 classes, the mask diameter, `cache=yes`, nparts.
   No usable result (no class averages, halves or FRCs): exit.
3. **The pool model's selection** of its classes (`select_with_model`'s logic), with
   `use_model=yes` (the default): the particles of rejected classes deselected for the rest of
   the route; the selected and rejected stacks with JPEGs in `refs_subset/quality_selection/`.
   The model keeping none: every populated class is used, with a warning (the one lenience
   kept, since it is not a failure of the data). With `use_model=no` the step is skipped and
   every populated class goes on.
4. **`refs_subset/cls_expansion/`:** cls_expansion into `expansion_budget` (500) sub-classes in
   all. Unusable: exit.
5. **`refs_subset/solve3D_cavgs/`:** solve3D_cavgs, `nstates`, c1, `prune=no`, the fixed
   20 → 12 Å ramp, `prob_athres_taper=yes`. No project, or no class with a state: exit.
6. **The check** (`register_attempt`'s logic): the states that hold class averages renumbered
   from 1; the same-state even/odd median, the state agreement, every pair docked in both
   hands, the populations; the volumes, their JPEGs and the per-state class averages into
   `state_refs/` with the record. **Fails: exit**, the record and the per-state stacks left for
   inspection (no second attempt, no fallback, no best attempt).
7. **The pool model scored on the full set's classes**, logged only (the selected and rejected
   stacks with JPEGs in `assign/quality_selection/`), so its verdict can be compared with the
   fit's (step 9), as p07 compares them.
8. **`assign/cls_expansion/`:** cls_expansion of the project's classes (every selected particle,
   the project's own 2D) into 500 sub-classes. Unusable: exit (no matching of the pool classes).
9. **`assign/`:** `assign_cavgs_states` on the expansion's project against `state_refs/`
   (`lp=15`, `nspace=1000`, the mask diameter fitted to the class averages' box), in-process.
10. **Rejection, flags, mapping** (`take_assignment`'s logic on the full project): the fit's
    outliers rejected (state 0), ambiguous sub-classes flagged, each particle given its
    sub-class's state and pose by `map2ptcls`, the fit's statistics and the particles the fit
    and the model reject logged; the per-state, rejected and ambiguous stacks and JPEGs in
    `assign/`. No candidate sub-class: exit.
11. **`refine3D/`:** seeds (half-sets by row parity, projection indices on the job's grid),
    the particles' states renumbered to the populated ones, one `refine3D refine=cont maxits=1`
    iteration, `objfun=cc`, c1, nparts. No volume registered: exit.
12. **The result:** the execution directory's project (`solve3D_pyramid.simple`) with every
    particle's state and refined pose, the states' volumes and FSCs registered (`vol`, `fsc`),
    and `state_refs/` beside it; a report `solve3D_pyramid_report.txt`, a text table as
    `solve3D_chunks_report.txt`: the first subset's size, the check's numbers (median, agreement,
    docked correlations, populations), the fit's statistics (median, MAD, threshold, rejected
    and ambiguous counts, the model's and the fit's rejections), and per state its class
    averages, particles and FSC 0.143 resolution.

## 3. Errors

`THROW_HARD` at every step above, after the step's own log and files are on disk, so the
failure can be read. Nothing is retried, skipped or replaced. A later version may add p07's
retries and fallbacks behind an input; this one does not.

## 4. Shared code (the main refactoring)

The method must not be copied out of p07. The stateless parts move to a module both use,
`src/main/solve/simple_state_assignment.f90` (`simple_state_assignment`), free procedures,
`use ..., only:`:

- rules: `refs_usable`, `reject_by_fit` (with its statistics), `flag_ambiguous`,
  `compact_states`, `expansion_budget`, `count_active_classes`, `fit_mskdiam`, `median_of`;
- the check's constants (`EO_ANGDIST_MAX`, `EO_OVERLAP_MIN`, `STATE_FRAC_MIN`, `NSTATES_MIN`,
  `REFS_PAIR_CC_MAX`, `DOCK_HP_REFS`), the matching's (`LP_ASSIGN`, `NSPACE_ASSIGN`), the fit's
  (`FIT_NMAD`, `FIT_REJECT_MAX`, `MARGIN_AMBIGUOUS`), the route's (`NCLS2D_FIRST`,
  `NCLS_EXPANDED`, `LPSTART_CAVGS3D`, `LPSTOP_CAVGS3D`), `MIN_PTCLS_PER_STATE`;
- helpers that read and write files but no stage: `docked_pair_cc` (now free, the test seam
  for skipping the docking becomes an argument), `write_state_cavgs`, `log_jpeg`,
  `cavgs_fbody`, the attempt's registration (`register_attempt` minus the stage's
  `try_angdist`: a function returning the check's values and whether it passed, writing
  `state_refs/try_<n>/`), the fit statistics log (`log_fit`), the assignment's mapping onto a
  project (`take_assignment`'s core: given the matched project and the rows, returns the new
  states and poses, the counts and the deselected rows), the refinement seeds
  (`prepare_refine_seeds`'s body on an `oris`), the refinement result's mapping
  (`take_refine_result`'s core);
- the command lines of the children as builders (`solve2D` for the first subset,
  `cls_expansion`, `solve3D_cavgs`, `assign_cavgs_states`, `refine3D refine=cont`), each a
  function of the few settings that vary (project file, nstates, mask diameter, nparts, nthr),
  so p07's `start_*` and the batch program set the same keys.

p07 keeps: the stage, its phases and job handling, the publications, the retries and
fallbacks, the GUI. Its `nopass` procedures become thin calls or go. p07's tests of the moved
rules move with them into a `simple_state_assignment_tester` (fast; the docking case stays in
the nightly `lib_stream` suite, or moves to `lib_reconstruction`).

This is a refactoring of p07 before the new program is written, and the step most likely to
break something; it is done first and gated by the existing `unit_stream` tests.

## 5. Registration

- Commander `commander_solve3D_pyramid` / `exec_solve3D_pyramid` in a new
  `src/main/commanders/simple/simple_commanders_solve3D_pyramid.f90` (CMake re-run), importing
  the child commanders (`commander_solve2D`, `commander_cls_expansion`, `commander_solve3D_cavgs`,
  `commander_assign_cavgs_states`, `commander_refine3D`) and `simple_state_assignment`.
- `simple_exec_solve3D.f90`: `case('solve3D_pyramid')`.
- `simple_ui_solve3D.f90`: `new_solve3D_pyramid`, UI_VIS_ADVANCED, inputs `nstates`,
  `nptcls_refs`, `use_model` (the existing parameter, as the sieve's UI declares it), `mskdiam`,
  `nparts`, `nthr` (and the constants in the help text).
- `simple_parameters.f90` / `_parse.f90`: `nptcls_refs` (int, default 0 = the program's
  constant).

## 6. Tests

- Unit (fast): the moved rules in `simple_state_assignment_tester` (today's p07 cases:
  `test_rules`, `test_fit_rules`, the budget, the first-subset cut on fixture rows); p07's
  tester keeps its stage tests on the shared module.
- High level (gated, nightly): `solve3D_pyramid_blobs`, the three-state phantom of 3,600
  particles in 60 exact views that `cls_expansion` and `refine3D_states_cavgs` already use
  (`flex_pca_blobs`), run end to end with `nstates=3`, `nptcls_refs=1200` (below the default,
  so the references come from a third of the phantom): the check passes,
  every reference docks onto one planted truth (the hand free), the assignment's purity against
  the planted states is at least a threshold to be set from the first run, and the refined
  maps correlate with the truths. This is the first end-to-end run of the method anywhere; p07's
  jobs have only run in the chained stream tests.

## 7. Docs

- `doc/policies/3D/solve3D_pyramid_policy.md`: the contract above (steps, exits, outputs,
  constants), pointing at the ingestion policy for the method's rationale.
- `stream_3D_ingestion_policy.md`: the shared module named as the owner of the rules.
- `doc/policies/3D/README` or the solve3D policy index, and the `simple-main-commanders` /
  `simple-refine3d` skill lines (with approval).

## 8. Order

1. Extract the shared module from p07 (section 4), p07 on it, `unit_stream` green.
2. The commander, exec case, UI, parameter; the report.
3. The high-level test on the phantom; set its thresholds from the first run.
4. Docs and skill lines.

## 9. Open questions

The design decisions are made (the nine above). What stays open is settled by the
implementation or by the first runs, not by a decision now:

1. **The phantom test's thresholds** (section 6): the assignment's purity against the planted
   states, and the docked correlation of each reference with its truth, are set from the first
   run, as the other high-level gates were.
2. **A first subset that is the whole project** (`nptcls_refs` at or above the selected
   particles): the route then classifies and expands every particle, and the assignment expands
   them again from the project's own 2D. Correct, but the second expansion repeats work; a
   later version may reuse the route's expansion when the subset is the whole set.
3. **The cut at `nptcls_refs`** follows whole stacks (the stack that crosses the count is left
   out, as p07's cap did), so the subset may be smaller than the count; a project whose first
   stack alone exceeds the count takes that stack's first particles up to it.
4. **No mask diameter anywhere** (no `cavg` entry with one in the project and none on the
   command line): an exit, since every step needs it; the message says which input to give.
5. **The matching's parts:** `assign_cavgs_states` runs in shared memory on `nthr` threads
   (its reprojections and alignments are in memory); `nparts` goes to the other children only.
6. **Where the docking test lives** once `docked_pair_cc` is shared: the nightly `lib_stream`
   "solve 3D reference docking" as now, or `lib_reconstruction` beside the addon's docking test;
   decided when the module moves (section 4). Settled: it stays in `lib_stream`, through p07's
   `register_attempt`.
