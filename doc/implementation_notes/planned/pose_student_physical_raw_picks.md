# The pose student revisited: the physical teacher, raw picks and alignment over the cone

Status, 2026-11-16: planned, not started. Base: `1ea0f8340` (pinned at launch), the commit that makes the 3D in-plane insertion the
cone with no key (`inpl_posterior` removed). Depends on: the pose student (`doc/algorithms/pose_student.md`), what the
class student established (`completed/student_inplane_posterior.md` with its report: the physical teacher P;
`completed/class_student_raw_picks.md` with its report: no centring, alignment after prediction over the candidate
list), and the current solve3D, whose solutions are the best yet (Hans, 2026-11-16).

## 1. Why

The pose student has not changed since its training-speed and working-grid runs except its prediction seed
(`euler_inpl_closest`, 2026-11-15). Meanwhile the class student showed four things that carry over:

- **The teacher should be what the run ends at, the physical posterior.** The class student's tempered teacher
  (`prob_beta` x T) made it over-cautious (97.7% where its threshold was calibrated for 95%); the physical teacher
  held the target and shortlisted 5 to 10 more particles in 100. Under the cone (refine3D since 2026-11-16) the table's
  view weights stand and the in-plane angle only spreads locally, so the pose student's matching teacher is today's
  one-rotation rows read at the physical temperature T, not the in-plane marginal (R, rejected for 3D on 2026-11-15).
- **Judge against the physical posterior's best pose**, not a tempered draw or the teacher's stored pose: the class
  student's old reference made safe shortlists look unsafe (72% against 97.7%).
- **No centring; predict, then align over the candidate set.** On raw picks the class student's mode matched the
  reference for 23 to 35% of particles while its list held it for 90 to 96%; aligning over the list beat the centred
  route on every set. The pose student's cone of directions (with their mirrors) plays the list's part, each
  direction's reprojection searched at every in-plane rotation over the shift disk (Hans's neighbourhood idea,
  `planned/polar_student_future_directions.md` section 6). Today's registration (a 5 x 5 grid over 25 Å at the read-out
  pose) cannot cover picking errors of 16 to 22 Å at the 90th percentile.
- **Costs on small sets are measured and reported, never gated** (Hans, 2026-11-15).

What is known about the pose student on uncentred particles: it seeds about 1% (`doc/algorithms/pose_student.md`), and
a 2 px shift moves its posterior by a total variation near 0.7. Whether its cone still holds the truth on uncentred
particles was never measured; that is Phase 0's first question (decision 9).

## 2. Contract

- Off state: refine2D and refine3D byte-identical at one thread with a fixed seed (`simulated_workflow`, both suites,
  `test=identity`); the class student byte-identical to the base (model header and weights).
- The training set stays the teacher's (solve3D's) particles at the teacher's poses and shifts; only the unseen
  particles change.
- No new command-line key or parameter; no second code path kept for a rejected variant.
- The networks keep their symmetries (rotation equivariance on the grid, mirror invariance, the conjugation symmetry)
  and training stays bit-identical at a fixed thread count.
- The refine3D route that consumes the student's cones (`prob_neigh geom`) is unchanged; whether the alignment over
  the cone can replace it is reported, not decided.

## 3. Decisions

1. (2026-11-16, Hans) The physical teacher **P3** for the pose student: the dense 3D table's rows read at the physical
   temperature T, no `prob_beta`, one rotation per direction as today; the threshold calibrated on it as now. The
   teacher audits itself on real rows as the class student's does (the row's masses against exp(-T (L - L_min))
   normalised over the directions at T; exact up to single-precision rounding).
2. (2026-11-16, Hans) **The validation's reference is the physical posterior's best pose**: from the dense 3D table at T
   against the teacher's map on the same unseen particles, with its shift disk, the direction and in-plane angle of
   the largest mass (mirror pairs folded as the student folds them).
3. (2026-11-16, Hans) **No pre-centring of unseen particles**: they enter at their picked positions (shifts zero).
4. (2026-11-16, Hans) **Alignment after prediction** against the teacher map's reprojections, every in-plane rotation of
   the polar grid and shifts over the disk of the reference pass, at the student's band, replacing today's
   registration. Two variants, both measured: **S1**, the mode's direction and its mirror; **SK**, every direction of
   the student's cone and their mirrors (the cone as the class student's candidate list), the best alignment kept,
   its pose committed. The search reuses the polar shift map (`shiftmap_srch`, as the class student's alignment does)
   unchanged; a change in `src/main/pftc` is a stop.
5. **Teachers** (decided by Hans on the proposal, 2026-11-16):
   - the simulated medium-assembly truth sets of the earlier runs (`T200`, `T1000`, `T5000`, `T` and the test half `S`,
     `/data2/agent_runs/student_speed/scratch/keep/sets/`, read in place, checksums checked), the truth map fixed;
     eight training seeds. The uncentred test half: `S`'s particles displaced by shifts drawn from the real sets'
     picking-error distribution (trpv1's: median 8 Å, 90th percentile 16 Å, 99th 40 Å; `completed/class_student_2D.md`),
     built with SIMPLE's programs and frozen with checksums;
   - three real sets, each with **three independent solve3D teachers** (teacher seeds 1 to 3, the current solve3D with
     the cone) and **three training seeds** per teacher, nine paired observations per set:
     bgal `/data2/bgal_addon/3_selection/bgal_addon.simple` (`pgrp=d2 pgrp_start=c1 mskdiam=180 nsample=1500`, as the
     earlier validations), trpv1 `/data2/trpv1_c2/3_selection/trpv1_c2.simple` (`mskdiam=175`, its point group and
     start group read from the project's earlier solve3D job; `nsample=5000`), streptavidin
     `/data/10335-streptav/tst/3_selection/tst.simple` (`pgrp=d2 mskdiam=85`, `nsample=5000`). `validate_pose_student`
     draws the training set from the selection's classes as now; its unseen particles enter at zero shift (decision 3)
     in the raw-pick arms and at their 2D shifts in the yardstick C0.
6. **Arms.** C0: today's pose student (teacher at `prob_beta` x T, today's registration) on unseen particles at their 2D
   shifts, the yardstick. R0: the same student on raw picks (what dropping the 2D shifts alone does). Then on raw picks
   with the alignment after prediction: A (today's teacher) and P3 (the physical teacher), each with S1 and SK.
7. **Reads**, on the unseen particles, against the reference (decision 2) and on the truth sets also against the truth:
   - seeded (route 1) share; the cone holding the reference's direction (all particles); the mode's agreement;
   - after the alignment: the committed direction within 6° of the reference (and of the truth), the in-plane angle's
     error, the shift's distance to the reference's shift (median and 90th percentile, Å);
   - the unseen map: `reconstruct3D` from the committed poses, its FSC against the teacher's map (cross-validated: the
     two share no particle), FSC 0.143 and 0.5;
   - cost per particle: prediction, alignment, their sum against the reference's dense search.
8. **Verdicts, written now** (paired by seed: the truth sets' eight training seeds, the real sets' nine teacher and
   training pairs; a read is worse when the paired mean difference lies in the worse direction beyond twice its standard
   error):
   - (a) **S1 or SK**: the one with the higher committed-direction agreement; within 2 SE, the cheaper (S1);
   - (b) **The teacher**: P3 is adopted unless it is worse than A on a read of decision 7 with the chosen alignment
     (the class student's case showed the physical teacher's conditional reads can dip while the route gains: if P3 is
     worse only on a conditional read and better on the committed-direction agreement and the unseen map, stop with the
     table for the maintainer rather than reject);
   - (c) **The raw-pick route**: the best raw-pick arm against C0 on every real set: committed-direction agreement and the
     unseen map's FSC not worse than C0's. Then the route needs no prior 2D shifts; otherwise stop with the table.
   Costs are reported, not gated.
9. **Early stop**: after Phase 0, if R0's cone holds the reference's direction for fewer than half of the particles on
   every set, no alignment over the cone recovers it: skip SK, continue with S1, record it.

Defaults taken without a question: S1 and SK measured on one prediction in one validation run; the FRC and FSC reads of
the earlier runs; `KMAX_STUDENT` as today's cone cap.

## 4. Not in this plan

- The class student (byte-identical, a guard); the amplitude channel (tested and rejected 2026-11-16, not to be
  revisited); the continuous refinement after prediction (`polar_student_future_directions.md` section 5).
- Replacing refine3D's `prob_neigh geom` pass for route 1 by the alignment: reported only.
- solve3D, refine3D and the stream; any file outside the table.

## 5. Change types and gates

| Phase | Kind | Gate |
| --- | --- | --- |
| 0 | measurement | the uncentred test half and the real-set teachers frozen; C0 and R0 tables; decision 9 |
| 1 | result-changing for the pose student's prediction and validation | testers of the alignment (a planted shift recovered over the disk; the mirror of a cone direction searched; SK keeps the best); the teacher's audit and its plants; fast gate, `unit_numerics`, `unit_pftc_align2D3D`; `pose_student` and `class_student` workflow tests; off-state identity; class student identity |
| 2 | experiment | decision 8 (a), decision 9 if it applies |
| 3 | experiment | decision 8 (b) |
| 4 | adoption or record | decision 8 (c); Phase 1's gates on the adopted tree; `validate_pose_student` on the three real sets at seed 1 |

## 6. Acceptance

- Paired standard errors come from the runs (measured).
- Off-state identity: `simulated_workflow` both suites on the base and the final tree at one thread, seed fixed,
  `test=identity`, exclusions only from the base-against-base comparison.
- Reported, not gated: the unseen maps of C0 and the adopted arm for the maintainer (three real sets, teacher 1, seed 1).

## 7. Common rules and stop conditions

- No floor relaxed; a file outside the table gets a row before the edit; the thread policy.
- Stops: the off-state identity breaks; the class student changes; the alignment needs a change in `src/main/pftc`,
  `src/main/volume` or `src/main/class`; a symmetry the networks would give up; decision 8 (b) or (c) ends at the
  maintainer.
- Gates no ruling may waive: the off-state identity, the class student's identity, thread identity of training.

## 8. Phases

- **Phase 0, the baseline.** Pin the base; Release and Debug builds, fast gate. Build and freeze the uncentred test half
  and the nine real-set teachers (record each solve3D's wall time and final FSC). C0 and R0 on every set at their
  seeds; the reads of decision 7 that exist without the alignment, including the cone's hold on the reference
  (from the prediction tables); decision 9. Exit: the tables frozen in Progress.
- **Phase 1, the physical teacher, the reference and the alignment.** In the tree: P3 behind the learn program's
  teacher (selected per arm only in a variant build for Phase 3's comparison; the tree carries one teacher after
  Phase 3), the validation's reference (decision 2), `predict_poses` aligning every predicted particle (S1 and SK) in
  place of the registration, `validate_pose_student` without the 2D shifts for the unseen particles in the raw-pick
  arms, decision 7's reads, a report that says what the cone and the alignment are. Testers and plants.
  Exit: section 5's Phase 1 gates.
- **Phase 2, the alignment's read.** A with S1 and SK on raw picks on every set. Exit: decision 8 (a); the losing
  alignment removed.
- **Phase 3, the teacher's read.** P3 against A with the chosen alignment on every set. Exit: decision 8 (b); the losing
  teacher removed.
- **Phase 4, adoption.** Decision 8 (c): the route kept, `doc/algorithms/pose_student.md` and the test policy updated;
  or the stop with the table. Exit: section 5's gates, the report.

### File table

| File | What | Phases |
| --- | --- | --- |
| `src/main/strategies/parallelization/simple_strategy_pose_student.f90` (and tester) | the teacher's temperature, prediction with alignment | 1, 2, 3, 4 |
| `src/main/strategies/parallelization/simple_pose_student_align.f90` (and tester) | the alignment over the cone (new) | 1, 2, 3, 4 |
| `src/main/commanders/simple/simple_commanders_pose_student.f90` | the validation: reference, raw picks, reads, report | 1, 2, 3, 4 |
| `src/main/learn/simple_pose_student.f90`, `simple_polar_student.f90` (and testers) | the teacher's rows and calibration at T | 1, 3, 4 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90` (and tester) | the reference's mode at T (read-only use of the table where possible) | 1, 4 |
| `src/defs/simple_defs_hhcnn.f90` | named constants | 1, 2, 3, 4 |
| `src/main/commanders/test/simple_commanders_test_highlevel_pose.f90`, `src/main/commanders/test/simple_commanders_test_class.f90`, `src/main/ui/simple_test/simple_test_ui_class.f90` | the workflow test, the unit suite of the alignment | 1, 2, 3, 4 |
| `doc/algorithms/pose_student.md`, `doc/policies/test_environment_policy.md` | the route, the test | 1, 4 |
| `doc/implementation_notes/README.md`, `doc/implementation_notes/planned/polar_student_future_directions.md` | record | 0, 1, 2, 3, 4 |

## Progress

(one entry per phase, written by the run)
