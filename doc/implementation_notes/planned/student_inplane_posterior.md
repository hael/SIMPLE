# The students learn the in-plane posterior

Status, 2026-11-15: planned, not started. Base: pinned in Phase 0, at or after `86babf86b` (the in-plane posterior
insertion by default in refine3D/solve3D; in solve2D since `d6cf4ae25`). Supersedes: none. Depends on: the in-plane
posterior weighting in 2D and 3D (landed; validated by Hans on solve2D and solve3D), the pose student
(`completed/pose_student_training_speed_report.md`) and the class student (`completed/class_student_2D.md`), which
share one core (`simple_polar_student`). Runs before `planned/pose_student_amplitude_channel.md`, which shares this
plan's Phase 0.

## 1. Why

Both students learn from a teacher row that gives each direction (3D) or class (2D) the loss at one in-plane
rotation: `score_ref` returns one distance and one rotation per reference, and `get_posterior_row` turns those into
the row's posterior. The students' readout matches that: for each direction or class it takes the maximum over the
rotation grid, so the training signal flows through one rotation, and the rotation the students read out is never
trained. The refinement no longer works that way: it inserts each particle under its in-plane posterior
(`inpl_quadrature`, whose `lnz` is the log of the sum over rotations, and `inpl_view_weights`, the in-plane
marginals). Hans (2026-11-15): the in-plane posterior should improve the network.

Three expected gains, each read separately:
- **a cleaner teacher**: the marginal over rotations is the posterior of a direction or class; the loss at one
  rotation out of a couple of hundred adds selection noise at low SNR. Better-calibrated targets should raise the
  share seeded at the calibrated target, the students' weakest read (pose student about 0.42 on bgal);
- **a matching readout**: a log-sum-exp over the rotations at a learned temperature in place of the maximum, so the
  student computes the teacher's quantity and every rotation carries gradient;
- **rotation supervision**: the teacher's in-plane posterior at the top directions or classes as a target, the
  first supervision of the rotation the students read out.

## 2. Contract

- Off state: refine2D and refine3D are byte-identical at one thread with a fixed seed. The teacher's marginal and its
  in-plane rows are computed only when a learn program keeps its rows (`keep_batch_rows`).
- The experiment lives outside the tree until the verdict: every variant, the teacher's extension included, is a
  patch on a `git archive` export of the base. Only what is adopted enters the tree, in Phase 4.
- The networks keep their symmetries: rotation equivariance on the rotation grid, mirror invariance (3D), the
  conjugation symmetry, bit-identical training at any thread count.
- No new command-line key or parameter. If adopted, the row file and both model files raise their versions; other
  versions are refused, as now.
- Neither student is closed by cost on small data sets (Hans, 2026-11-15): costs are reported, not judged, here.

## 3. Decisions

1. (2026-11-15, Hans) The in-plane posterior for both students, the 2D class student and the pose student, through
   their shared core.
2. **The teacher's marginal.** For every reference of a kept row, `score_ref` (and the 2D table's scorer) also returns
   the log-sum over the rotation grid of exp(-T (d_r - d_min)) at the particle's temperature T, the quantity
   `inpl_quadrature` calls `lnz`, at the shift the fill scored at. The row's posterior then reads
   p(ref) ∝ exp(lnz(ref) - T d_min(ref)), the in-plane marginal, as `inpl_view_weights` forms it. The loss at the
   one rotation stays in the row for the base route.
3. **The rotation targets.** For the `KMAX_INPL_ROW` = 8 references of largest marginal mass, the in-plane
   posterior over the teacher's rotation grid, binned onto the student's rotation grid by angle. In 3D a mirror
   pair's targets follow the readout's rule (the unmirrored direction's rotation, or the reflected one mapped through
   the pair's frame angle); Phase 1's tester pins the mapping at 0, 75 and 180°, as the HHCNN tester pins the readout.
4. **Variants.** A, the base. M, marginal targets (decision 2) and a log-sum-exp readout at a learned temperature,
   initialised to the base's scale. R, M plus a rotation loss: the cross-entropy between the student's in-plane
   softmax and decision 3's targets at the same 8 references, weight 1 (guessed). M and R are read apart so that each
   gain is attributed.
5. **Teachers**, rebuilt on the base in Phase 0 and frozen with checksums, three training seeds each:
   - 3D: the medium truth set's `T200`, `T1000`, `T5000`, `T` (24,000) and bgal, as in the amplitude plan;
   - 2D: bgal (2,756 training particles, 50 classes) and trpv1 (10,000, 100 classes), the solve2D teachers of the
     class student's validation, rebuilt with the current solve2D.
6. **Reads.**
   - 3D, held-out: within 6°, the truth in the cone, the share seeded, held-out recovery; the readout's in-plane
     error against the teacher's mode rotation (median, agreeing direction);
   - 2D, held-out: classes recovered, the share seeded and the threshold, the in-plane error; and
     `validate_class_student`'s report: particles given a shortlist, shortlists holding the full search's class
     (centred and perfectly centred), shortlist length, class averages against the teacher's;
   - both: training and prediction time, reported, not gated;
   - beside the 2D reads, reference matching at the working box (`completed/class_student_vs_reference_matching.md`),
     reported.
7. **Verdict, written now.** M (or R) is adopted for a student only if:
   - (a) no held-out read's mean over the three seeds falls below A's range over its seeds, on any of that student's
     sets;
   - (b) the share seeded at the calibrated target, or within 6° (3D) / classes recovered (2D), rises above A's
     range by more than twice that range, on at least three of the five 3D sets for the pose student, on both sets
     for the class student (margin guessed, re-derived from Phase 0's spread before Phase 2);
   - (c) for R over M: the in-plane error falls by more than twice M's range with (a) holding against M.
   Otherwise the variant is rejected; the record goes to `rejected/` (Phase 4).
8. **Early stop.** Phase 2 reads M and R on `T1000`, bgal (3D) and bgal (2D). If neither moves any read of (b) beyond
   A's range, the experiment stops as a negative result.
9. **One schedule.** All variants train with the base's schedule, designs and shift augmentation; the class
   student at width 16 (`CLASS_STUDENT_NDESIGNS`), the pose student with its design choice.

Defaults taken without a question: `KMAX_INPL_ROW` = 8; the rotation loss's weight 1; the log-sum-exp temperature
learned per student, one scalar.

## 4. Not in this plan

- The amplitude channel (its own plan, after this one; its base is this plan's outcome).
- The insertion, the tables' draws and the refinement's use of the cones or shortlists.
- The stream (`planned/class_student_stream.md`).
- Any file outside the table.

## 5. Change types and gates

| Phase | Kind | Gate |
| --- | --- | --- |
| 0 | measurement | the base's reads at three seeds, frozen; shared with the amplitude plan |
| 1 | patch outside the tree | the teacher's marginal against a brute-force sum; off-state identity of refine2D and refine3D at one thread; the readout's symmetries and thread identity in the HHCNN tester |
| 2 | experiment | decision 8 |
| 3 | experiment | decision 7 |
| 4 | result-changing for the students (if adopted) | fast gate; `unit_numerics` (HHCNN, polar student, pose student, class student suites), `unit_pftc_align2D3D` (both probability tables); `pose_student` and `class_student` workflow tests; `validate_pose_student` and `validate_class_student` on bgal |

## 6. Acceptance

- Spreads come from Phase 0 (measured), recorded before Phase 2. The margin of 7 (b) and the rotation loss's weight
  are guessed and say so.
- Off-state identity: refine2D and refine3D on the base and on the patched build at one thread, seed fixed, the
  project records and the class averages or maps byte for byte (`test=identity`).
- Reported, not gated: the class averages of the student's assignments (2D) and the unseen map (3D) for A and the
  adopted variant, for Hans to inspect.

## 7. Common rules and stop conditions

- No floor relaxed; a file outside the table gets a row before the edit; the thread policy.
- Stops: the off-state identity breaks; the brute-force check of the marginal fails; a teacher cannot be rebuilt;
  a symmetry the readout would have to give up; the row's memory beyond 2 GB for `T` (8 references × the rotation
  grid per row is planned at under 1 GB).
- Gates no ruling may waive: the off-state identity, the thread identity of the training step, decision 7 (a).

## 8. Phases

- **Phase 0, baseline (shared with the amplitude plan).** Pin the base. Release and Debug builds, fast gate. Rebuild
  the seven teachers; train A at three seeds; every read of decision 6; reference matching beside the 2D reads.
  Exit: the reads and spreads frozen in Progress, the teachers' checksums.
- **Phase 1, the teacher and the readout.** On an export: the marginal and the in-plane rows in both tables and the
  row file; M and R in the core and the HHCNN; the testers (brute-force marginal, the mirror mapping of decision 3,
  the readout's equivariance, invariance and thread identity).
  Exit: the testers pass; the off-state identity holds.
- **Phase 2, first read.** M and R on `T1000`, bgal (3D) and bgal (2D), three seeds.
  Exit: decision 8 applied, continue or stop recorded.
- **Phase 3, the matrix.** M and R on every set.
  Exit: decision 7 applied per student, the verdict recorded.
- **Phase 4, adoption or rejection.** Adopted: the patch into the tree with its tests, the versions raised, the
  algorithm pages (`doc/algorithms/pose_student.md`, `class_student.md`). Rejected: the record beside this plan as `student_inplane_posterior_rejected.md`, which the
  maintainer moves to `rejected/`. Exit: section 5's gates, or the rejected note.

### File table

| File | What | Phases |
| --- | --- | --- |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90`, `simple_eul_prob_tab2D.f90` (and testers) | the marginal and in-plane rows of kept rows | 1 (export), 4 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_utils.f90` | `prob_row_record`, `prob_row_file` (format, version) | 1 (export), 4 |
| `src/main/learn/simple_hhcnn.f90`, `simple_hhcnn_bench.f90` (and tester) | the log-sum-exp readout, the rotation loss | 1 (export), 4 |
| `src/main/learn/simple_polar_student.f90`, `simple_pose_student.f90`, `simple_class_student.f90` (and testers) | targets, training, the model files' versions | 1 (export), 4 |
| `src/main/strategies/parallelization/simple_strategy_pose_student.f90`, `simple_strategy_class_student.f90`, `simple_strategy_student_common.f90` (and testers) | the learn programs ask for the rows' extension | 1 (export), 4 |
| `doc/algorithms/pose_student.md`, `doc/algorithms/class_student.md` | the teacher and the readout | 4 |
| `doc/implementation_notes/planned/student_inplane_posterior_rejected.md` | the rejected record, if both variants are rejected (the maintainer moves it to `rejected/`) | 4 |
| `doc/implementation_notes/README.md`, `doc/implementation_notes/planned/pose_student_amplitude_channel.md` | record (this plan and its report are always allowed) | 0, 1, 2, 3, 4 |

## Progress

(one entry per phase, written by the run)
