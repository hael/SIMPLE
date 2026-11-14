# Student prediction in parts

Status, 2026-11-14: delivered, uncommitted (report: `student_predict_parts_report.md`). Base: `38e033cbf` with the class-student diff uncommitted in the working tree (the
2D plan, `completed/class_student_2D.md`, through its ruling 18). Supersedes: none. Depends on: the pose student
(validated: yes, `completed/pose_student_training_speed_report.md`) and the class student (validated through Phase 4
of the 2D plan, uncommitted).

## 1. Why

Route B of solve2D at scale (the 2D plan's Phase 4: teacher, student, `center_ptcls`, prediction, the seeded sparse
stage) costs about 3.4 ms per particle at 24 threads, against route A's 13 ms. Every per-particle step of it runs
distributed except the prediction: `center_ptcls` hands `nparts` to its `refine2D`, the sparse stage is
`refine2D`, but `predict_classes` has no part mode and `predict_poses` refuses parts when it registers shifts.
At 0.35 to 0.6 ms per particle (trpv1, 24 threads) a million particles is 6 to 10 minutes on one node, read through
one node while the rest of the route runs on the cluster. The validation controllers (`validate_class_student`,
`validate_pose_student`) accept no `nparts` at all. Training stays in one process (2D plan, ruling 18's record):
its cost is set by the teacher's sample, not by the data set.

## 2. Contract

- **Off state.** Without `nparts`, both programs run as on the base, byte for byte at one thread and a fixed seed.
- **Identity in parts.** `nparts=2` and `nparts=3` give the shared-memory run's particle records, prediction table
  and (2D) candidate artifact byte for byte at one thread per part. Each particle's prediction depends only on its
  image, the model and, when registering, the references and sigma2: nothing crosses particles.
- **One project owner.** Workers read the project and write part files only (`binwrite_oritab` to `outfile`, their
  table, their candidates). The master makes the halves and sigma2 before scheduling, merges, and writes the
  particle segment and the `out` registration once (`doc/policies/layering_and_project_ownership_policy.md`,
  section 5).
- **No new command-line key.** `nparts`, `part`, `fromp`, `top` and `outfile` exist and `predict_poses` already reads
  them.

## 3. Decisions and rulings

1. (2026-11-14, Hans) Prediction in parts is built for both students. Training is not distributed: neither the
   teacher's rows (3% of the 2D learn, 15% of the 3D one) nor one design per node nor gradient exchange.
2. (2026-11-14) **The worker** is the existing pass with a range. `predict_classes_exec` takes `predict_poses_exec`'s
   part mode: `student_part_pinds` cuts the selection to `fromp:top`, the table is `<fname>_part<i>_pred.txt`, the
   range's records go to `outfile` with `binwrite_oritab`, an empty range writes its records unchanged. Both
   workers end with `qsys_declare_part_finished`, as `make_cavgs`' worker does: `predict_poses`' part mode does not
   today (its parts were only run in process, by the pose workflow test), and a master would wait for it forever.
3. (2026-11-14) **Registration in parts.** The refusal in `predict_poses_exec` ("the shift registration (vol1) runs
   in shared memory") goes. A worker does not make shared state: `ensure_field_halves` already only checks in a
   worker (`l_distr_worker_glob`); for sigma2 a worker calls `sigma2_estimate_available` and stops when it is
   false, instead of `ensure_sigma2_for_iteration`, which would run `calc_pspec` and publish.
4. (2026-11-14) **The master** is the program itself with `nparts > 1` and no `part`, decided as
   `create_make_cavgs_strategy` decides. In order: read the project; when registering, `ensure_field_halves` and
   `ensure_sigma2_for_iteration`; `qenv%new(params, nparts)`, the job description from the command line with
   `mkdir=no`, `gen_scripts_and_schedule_jobs(job_descr, algnfbody=ALGN_FBODY, array=L_USE_SLURM_ARR, ...)`;
   `merge_algndocs(nptcls, nparts, oritype, ALGN_FBODY, write=.false.)`; the parts' tables concatenated in part
   order into `<fname>_pred.txt` under the one header; (2D) the candidates merged and published; the particle
   segment written once; the part files deleted after a successful merge.
5. (2026-11-14) **The candidates in parts.** `class_candidates` gains a part file: the artifact's three arrays over
   the part's rows, with the range in its header, written by the worker. The master's publish reads the parts in
   part order, stops unless they cover `1..nrows` without overlap (the rule of `merge_algndocs`), and writes the one
   artifact under the master's layout digest. The artifact's format (`CLSCAND1`) does not change.
6. (2026-11-14) **The callers.** `validate_class_student` and `validate_pose_student` take `nparts` (optional, in
   their UI) and hand it to every child that runs distributed: the teacher (`solve2D` / `solve3D`),
   `center_ptcls`, the prediction, the dense pass and the averages or reconstruction. The learn programs run
   without it.
7. (2026-11-14) **Tests.** The registered `class_student` and `pose_student` workflow tests each add one master run
   at `nparts=2` (the local queue, as `simulated_workflow suite=1jxy` runs its distributed solve3D), compared byte
   for byte with their shared-memory prediction. The pose test's in-process `parts` stage stays. No new CTest entry.

Defaults taken without a question:
- the master's table keeps the shared-memory name `<fname>_pred.txt`; the parts' candidate files are
  `<fname>_part<i>_cands.bin`;
- the master is a routine in the student's strategy module (`predict_classes_master`, `predict_poses_master`)
  chosen in `run_predict_*`, not a new strategy type: the programs are one pass with no iteration or hooks;
- the 3D shift-check table (`<fname>_shiftcheck.txt`) stays per part and is not merged: it samples the first
  particles of whatever range it is given.

## 4. Not in this plan

- Distributing `learn_class_student` or `learn_pose_student` in any form (decision 1).
- The stream (`planned/class_student_stream.md`, owners Joe and Cyril): its learn and predict jobs may use the
  master later; nothing under `src/main/stream/` or `src/main/commanders/stream/` changes.
- The sparse 2D table's reading of the artifact (`eul_prob_tab2D`), refine2D, refine3D and `center_ptcls`.
- The queue system (`src/utils/qsys/`).

## 5. Change types and gates

| Phase | Kind | Gate |
| --- | --- | --- |
| 0 | measurement | builds, fast gate, the base's runs recorded |
| 1 | identical without `part` | off-state identity; in-process parts equal the whole run |
| 2 | identical without `nparts` | off-state identity; master at `nparts` 2 and 3 equal to shared memory |
| 3 | identical without `nparts` | both workflow tests pass with their new master arm |
| 4 | measurement | reported, not gated |

## 6. Acceptance

- **Identity** (decision-free, from the contract): at one thread per part and a fixed seed, the merged particle
  segment (record by record, `test=identity`'s comparison), `<fname>_pred.txt` and the candidate artifact's bytes
  equal the shared-memory run's. Exclusions, declared now: `BENCH/` files, logs, the 3D shift-check tables, the
  part files. Each identity case runs twice on the base first (Phase 0); an item that cannot be made identical is
  restored and recorded, never rerun until it passes.
- **Off state:** the base's `predict_classes` and `predict_poses` outputs, with and without registration, byte for
  byte at one thread, on the workflow tests' fixtures.
- **Reported, not gated** (Phase 4): the wall time of predicting trpv1's 25,645 unseen particles at `nparts` 1, 2
  and 4 (24, 12 and 6 threads) on the Dell, with the master's scheduling and merge times. One node cannot show a
  gain; it shows the overhead. The run at real scale is Hans's (a data set of millions on the cluster), named in the
  report with its command line.

## 7. Common rules and stop conditions

- No new parameter; no floor or tolerance relaxed; a file outside the table gets a row before the edit and the
  report lists it. Thread and OpenMP rules of `simple-modern-fortran`; no new parallel region is planned.
- Stops: a worker would have to write the project or the sigma2 store; identity fails on an item outside the
  declared exclusions after one diagnosis; a change is needed in a file of section 4; a new command-line key would
  be needed.
- Gates no ruling may waive: the off-state identity, the identity in parts, the fast gate, both workflow tests.

## 8. Phases

- **Phase 0, baseline.** Release and Debug builds, the fast gate. `class_student` and `pose_student` workflow tests
  twice each. `predict_classes` (with and without `refs`) and `predict_poses` (with and without `vol1`) twice at
  one thread on the workflow fixtures, outputs kept as the identity reference; once at 8 threads, compared
  (informative). `predict_classes` on trpv1's unseen particles at 24 threads for the Phase 4 base.
  Exit: the reference outputs and their two-run comparison recorded; the exclusions confirmed as the only
  differences, or added here with the reason before Phase 1.
- **Phase 1, the 2D worker** (decisions 2, 3, 5). `predict_classes_exec`'s part mode; `class_candidates`' part write
  and merging publish with their unit tests (`test=unit_project`, the `class candidates` suite: coverage, overlap, an empty part, a part of
  unselected rows); the commander checks `part` with `nparts`, `fromp`, `top` (as `predict_poses` does); the
  registration's worker check of decision 3 in both students. The class workflow test runs parts in process, as the
  pose test's `parts` stage does. Exit: the unit suite passes; in-process parts give the whole run's table rows and
  records; off-state identity.
- **Phase 2, the masters** (decision 4). Both masters; the 3D refusal removed. Exit: at one thread per part,
  `nparts=2` and `nparts=3` (one part holding only unselected rows) equal to shared memory, both students, with and
  without registration; no part file left after the run; off-state identity.
- **Phase 3, the callers and tests** (decisions 6, 7). Exit: both workflow tests pass with their master arm under
  CTest, the times recorded; `validate_class_student nparts=2` and `validate_pose_student nparts=2` run on bgal.
- **Phase 4, measurement and report.** Exit: section 6's table; the algorithm pages (`doc/algorithms/class_student.md`,
  `pose_student.md`) say how prediction runs in parts; the report beside this plan.

### File table

| File | What | Phases |
| --- | --- | --- |
| `src/main/strategies/parallelization/simple_strategy_class_student.f90` | 2D worker and master | 1, 2 |
| `src/main/strategies/parallelization/simple_strategy_pose_student.f90` | 3D master; registration in parts | 1, 2 |
| `src/main/strategies/parallelization/simple_strategy_student_common.f90` | shared master helpers (table concatenation), if needed | 2 |
| `src/main/project/simple_class_candidates.f90` and its tester | the part file and merging publish | 1 |
| `src/main/commanders/simple/simple_commanders_class_student.f90` | part keys check; `validate_class_student` | 1, 3 |
| `src/main/commanders/simple/simple_commanders_pose_student.f90` | `validate_pose_student` | 3 |
| `src/main/ui/simple/simple_ui_refine2D.f90`, `simple_ui_solve3D.f90` | `nparts` on the validation programs | 3 |
| `src/main/ui/simple_private_prgs.f90` | the private programs' keys, if needed | 1 |
| `src/main/commanders/test/simple_commanders_test_highlevel_class_student.f90`, `_pose.f90` | in-process parts; master arm | 1, 3 |
| `doc/algorithms/class_student.md`, `doc/algorithms/pose_student.md` | prediction in parts | 4 |
| this plan, its report, `doc/implementation_notes/README.md` | record | 0-4 |

## Progress

### Phase 0, baseline (2026-11-14, the Dell, interactive session)

- Base: `38e033cbf` with the class-student diff and ruling 18 (`CLASS_STUDENT_NDESIGNS` = 1) in the working tree,
  Release, GCC 15.2, 24 cores, load under 4. Fast gate passed. The base binaries are kept in `$TMPDIR/spp_base/bin`.
- `class_student` and `pose_student` ran twice each at one thread on the base (`$TMPDIR/spp_p0/`): every read of
  `metrics.tsv` equal between the two runs except the wall times and the pose test's `prediction_ms_per_particle`.
  Those are the exclusions of the off-state comparison.
- What the code check changed in the plan before Phase 1: `predict_poses`' part mode never declared its part
  finished (decision 2, added to the plan); the pose registration's reprojection model is written by every part
  into the shared directory (found in Phase 2, below).
- Not run: the separate one-thread `predict_*` reference outputs of the plan's Phase 0. The workflow tests delete
  their directories when they pass, so the identity is checked inside them instead (Phase 2), and the off state by
  the one-thread reads above.

### Phase 1, the 2D worker (2026-11-14)

- `class_candidates` gained `write_part` (a part's rows with its range in a `CLSCANDP` header) and `publish_parts`
  (the parts read in order, coverage and overlap checked as `merge_algndocs` does, one `CLSCAND1` artifact), with
  `publish` and both sharing `fill_rows` and `write_and_register`. New unit test `test_parts` (`unit_project`,
  `class candidates`): three parts, one with no predicted row, merge into the bytes of the whole publish. Passes.
- `predict_classes_exec` has the part mode of `predict_poses_exec`: the range by `student_part_pinds`, the table
  `<fname>_part<i>_pred.txt`, the range's records to `outfile` (`binwrite_oritab`), the candidates part file; a range
  without a selected row writes its header, records and an empty candidates part (`write_empty_part`).
- Both workers end with `qsys_declare_part_finished`; both take sigma2 through `student_sigma2_for_pass`, which in
  a part only checks that the store is consumable (`sigma2_estimate_available`) and stops otherwise.
- The commander checks a part's keys; `simple_private_prgs.f90` lists them for `predict_classes`.

### Phase 2, the masters (2026-11-14)

- `student_predict_parts` (`simple_strategy_student_common.f90`) is the masters' shared frame: `qsys_env` over the
  project's rows balanced over the selected ones (`l_active`, as refine2D's master), the command line as the job,
  `merge_algndocs(write=.false.)`, the tables joined byte for byte (`concat_part_tables`), the part tables deleted,
  `qsys_cleanup`. `predict_classes_master` and `predict_poses_master` add their setup and write the field once;
  the 2D master publishes the candidates from the part files and deletes them.
- **A finding the plan missed:** with `vol1`, every part ran `materialize_reprojection_model_from_volumes`, which
  removes and rewrites `reprojection_model_{even,odd}.bin` in the shared directory; part 2 found them gone and
  stopped. Now the master builds the model once (after the halves and sigma2, as the shared-memory run orders
  them) and a part only reads it. The 2D preparation (`prep_pftc4align2D`) writes no file.
- Exit: the `pose_student` test's new `master` stage (shm/ against distr/, both registering against the truth map)
  and the `class_student` test's `parts_check` (shared memory against the master, from copies of the predicted
  project and its sigma2) give the same table, candidates and records, byte for byte, at 8 and at 1 thread; the
  masters leave no part table, candidates part, algndoc, finished flag or script. Off state: the one-thread reads of
  both tests equal the base's (Phase 0's exclusions). Fast gate passed.
- Not as planned: the identity is shown at two parts, not at three with one part of unselected rows. The master
  balances the parts over the selected rows, so it makes no such part; the empty part is checked in process
  instead (the class test's `parts_check/empty`, the pose test's `parts` stage), and `test_parts` merges one.

### Phase 3, the callers and tests (2026-11-14)

- `validate_class_student` and `validate_pose_student` hand `nparts` to every step but the training (`set_parts`);
  the teacher's `solve2D` / `solve3D` already took it from the entry's line. `validate_class_student`'s UI offers
  `nparts`; `validate_pose_student`'s already did (`add_solve3D_base_inputs`).
- Workflow tests at 8 threads: `class_student` 18.7 s (22.7 s before the master arm, on a busier machine),
  `pose_student` 31.0 s (31.9 s in the training-speed report). At one thread: 106.2 s and 232.3 s, logs in
  `$TMPDIR/spp_p3/`.
- Phase 3's real-data runs (`/data2/agent_runs/student_predict_parts/runs/`, `tools/run_p34.sh`):
  `validate_class_student` on bgal (seed 1) at 2 parts x 12 threads ran every step in 62.6 s; its reads are in the
  range of the earlier shared-memory seed-1 run (the teacher's solve2D now runs in parts, so they cannot be equal).
  `validate_pose_student` on bgal at 2 x 12 ran every step in 585 s, its prediction in parts with registration in
  5.3 s; its reads are in the report.

### Phase 4, measurement and report (2026-11-14)

- `predict_classes` on trpv1's 25,645 unseen particles (the 2D run's `p3_validate_trpv1_s1` model, averages and
  project copy), registering: 11.5 s at one part x 24 threads, 10.8 s at 2 x 12, 11.3 s at 4 x 6; tables and
  candidate files identical across the three. One node shows no overhead and no gain, as section 6 expected.
- The algorithm pages say how prediction runs in parts (`doc/algorithms/class_student.md`, `pose_student.md`).
- Exit: the report beside this plan. The run at real scale stays with Hans (the report's last section).

