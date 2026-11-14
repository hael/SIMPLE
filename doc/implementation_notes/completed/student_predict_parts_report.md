# Student prediction in parts, report

An interactive session on 2026-11-14, on the Dell (Oracle Linux 8.10, Core Ultra 9 285, 24 cores, GCC 15.2,
Release), carried out `student_predict_parts.md` beside this report. Its Progress section holds each phase's
detail. Nothing is committed: the change is part of the working-tree diff on `38e033cbf` that also holds the
class student.

## What was learned

- **Prediction in parts gives the shared-memory result byte for byte**, as the contract required: the joined table,
  the candidate artifact and the particle records. The two workflow tests check it at two parts, at 8 and at 1
  thread, with the shift registration on. On trpv1 (25,645 unseen particles, 100 classes, registering), one, two
  and four parts at 24, 12 and 6 threads gave identical tables and candidate files.
- **The pose registration could not run in parts as the plan described it.** Every part rebuilt the reprojection
  model in the shared directory, removing the files another part was reading. The master now builds it once, as
  refine3D's master does, and the parts read it.
- **`predict_poses`' existing part mode could never have run under a master**: it did not declare its part
  finished, and nothing but the pose test (in process) had run it.
- **On one node, parts cost nothing and gain nothing**: trpv1's prediction took 11.5, 10.8 and 11.3 s at one, two
  and four parts. The gain is for a cluster, where the prediction no longer reads every particle through one node.

## What changed, phase by phase

- Phase 1 (the 2D worker): `class_candidates` part file and merging publish, with a unit test; `predict_classes`'
  part mode; parts declare themselves finished; parts check sigma2 instead of making it. 5 files.
- Phase 2 (the masters): the shared master frame and byte-exact table joining in
  `simple_strategy_student_common.f90`; `predict_classes_master` and `predict_poses_master`; the reprojection model
  built by the pose master; master arms in both workflow tests. 5 files.
- Phase 3 (the callers): `nparts` handed to every step but the training in both validation programs; the UI input
  on `validate_class_student`; the private programs' part keys for `predict_classes`. 4 files.
- Phase 4: the measurements below; the algorithm pages.

## Decisions taken under delegation

- The master is a routine chosen in `run_predict_*` (the plan's default), sharing one frame,
  `student_predict_parts`, between the students, so the frame exists once.
- The parts are balanced over the selected rows (`qsys_env` with `l_active`, as refine2D's master). A consequence:
  the master makes no part of unselected rows only, so the plan's "three parts, one empty" identity case was
  replaced by in-process empty parts (Progress, Phase 2).
- The pose shift-check table stays per part, as the plan's default said. `EXIT_CODE_JOB_*` and
  `fft_wisdoms_part*.dat` stay after a master run, as after the other masters (`qsys_cleanup` does not remove
  them).

## Deviations from the plan

- Phase 0's separate one-thread reference outputs were not made: the workflow tests delete their directories on a
  pass, so the identity is checked inside the tests and the off state by the one-thread reads of both tests against
  the base build.
- The reprojection model's ownership (above) was not in the plan; `simple_strategy_pose_student.f90` was in the
  file table. No file outside the table was changed.

## Evidence against the baseline

- Fast gate: passed after every phase.
- `unit_project`: 3,066 assertions passed, `test_parts` among them.
- `class_student` at 8 threads: all gates passed (18.7 s); the master checks: table, candidates and records as shared
  memory, no part file left, the empty part's contributions. At 1 thread: the same, and every read the base's
  (exclusions: wall times).
- `pose_student` at 8 threads: all gates passed (31.0 s); master table and records as shared memory, registering;
  no part file left. At 1 thread: the same, every read the base's (exclusions: wall times,
  `prediction_ms_per_particle`).
- Real data, reported, not gated (`/data2/agent_runs/student_predict_parts/runs/`):

| run | parts x threads | wall | result |
| --- | --- | --- | --- |
| `predict_classes`, trpv1 25,645 unseen, registering | 1 x 24 | 11.5 s | reference |
| same | 2 x 12 | 10.8 s | table and candidates identical |
| same | 4 x 6 | 11.3 s | table and candidates identical |
| `validate_class_student` bgal, seed 1 | 2 x 12 | 62.6 s | every step ran; reads in the range of the earlier shared-memory run (top-1 0.33 against 0.30, dense class among candidates 0.74 against 0.70); learn 11.9 s against 45.8 s (width 16 alone, ruling 18 of the 2D plan) |
| `validate_pose_student` bgal, 1,500 training particles | 2 x 12 | 585 s | every step ran (solve3D 494 s, learn 28 s, prediction in parts with registration 5.3 s, six reconstructions 57 s); training particles' registered shift median 0.59 px from the teacher's; unseen map from the student's poses and registered shifts 7.59 Å between its halves, the teacher's control 4.13 Å |

## Left open for the maintainer

- **The run at real scale** (plan section 6): prediction in parts on a cluster over a data set of millions, with
  `nparts` on `predict_classes` (or the validation programs). Command line, for example:
  `simple_private_exec prg=predict_classes projfile=<p> infile=<model> refs=<cavgs> mskdiam=<d> nthr=<n> nparts=<N>`.
- The stream plan (`class_student_stream.md`) can now start its predict job as a master; that is its owners' call.

## Commit message

```
Student prediction in parts

predict_classes and predict_poses run as distributed masters with nparts:
the parts predict their ranges, the master merges their records, joins
their tables and (2D) publishes one candidate artifact, and writes the
project once. The result is the shared-memory run's byte for byte; both
workflow tests check it. Shift registration now runs in parts: the
master makes the halves, sigma2 and (3D) the reprojection model, which
the parts only read. Parts declare themselves finished. The validation
programs hand nparts to every step but the training.
```
