# The students learn the in-plane posterior

Status, 2026-11-15: delivered: P adopted for the class student, 2026-11-15; pose student unchanged (ruling R7). P
is the in-plane marginal at the physical temperature as the teacher (no `prob_beta`), with the log-sum-exp readout
and the rotation loss, and `validate_class_student` judges against the physical posterior's best class. M and R
were rejected: by decision 7 in Phase 3, and R again under ruling R4. The report is
`student_inplane_posterior_report.md` beside this plan, and the plan stays in `planned/` until the maintainer moves
it. Bases: `86babf86b` for Phases 0 to 3 and R4's first audit (ruling R1), `4282d58a3` for R4's readings (ruling
R5), `603c44490` for R6's readings and the adoption (rulings R6 and R7). Supersedes: none. Depends on: the in-plane
posterior weighting in 2D and 3D (landed; validated by Hans on solve2D and solve3D), the pose student
(`completed/pose_student_training_speed_report.md`) and the class student (`completed/class_student_2D.md`), which
share one core (`simple_polar_student`).

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
   Made precise in Phase 1 (ruling R2, from the code): T is the temperature the table both draws its in-plane
   rotations at and reads its rows at, the tempering factor times the particle's own (`prob_beta` × T, the base
   targets' too). The fill already computes each entry's whole in-plane row to draw its rotation; the marginal's
   route computes that row once, draws the same rotation with the same keyed variate (the entry unchanged, bit for
   bit) and sums it, so the calculator is untouched (`simple_eul_prob_marginal`). An entry the fill refined in shift
   takes the row at its refined shift. The row's losses become d = d_min - lnz / T, whose tempered weight is the
   marginal, and the row record's mode rotation becomes the mode reference's in-plane mode (its most probable
   rotation, no longer one draw). The drawn losses stay in the table for its own assignment; the rows file carries
   the training losses with lnz and d_min beside them.
   **Revised in Phase 4** (ruling R2, on the audit ruling R4 asked for): the teacher's in-plane weighting must be the
   refinement's, and it was not. The learn programs run at `prob_beta` = 0.3, so the marginal above was summed at
   0.3 T, while the in-plane posterior insertion weighs every row at the particle's physical temperature T (3D
   `inpl_pose_row`, 2D `insert_inplane`) and, in 3D, adds the shift prior's term to each pose's view weight
   (`inpl_prior_penalty`); the sum was also taken in double where `inpl_quadrature` rounds each term to single
   precision first. Now each entry's row is drawn as before (the table's entries unchanged, bit for bit) and weighed by
   `inpl_quadrature` itself (an optional output `post`, the row's posterior, added; no existing result changes) at the
   insertion's temperature and at the shift the insertion would take the entry at (the refined entry's, else the seed
   as `record_seed_shift` records it); in 3D the entry's shift-prior term is `inpl_prior_penalty`'s at that shift and
   centre, in 2D none (the 2D insertion has no view weights). A reference's training loss is d = d_min - (lnz - pen)
   / T, whose weight at T is the insertion's view weight exp(lnz - T d_min - pen); the row is read at `prob_beta` T
   as A's drawn losses are, so the class targets are the insertion's view weights tempered as the table tempers its
   rows. The rotation targets (decision 3) are the insertion's in-plane posteriors. The learn programs audit this on
   their first batch of real rows against the insertion's own route, and the targets' mapping on every row (Progress,
   Phase 4).
3. **The rotation targets.** For the `KMAX_INPL_ROW` = 8 references of largest marginal mass, the in-plane
   posterior over the teacher's rotation grid, binned onto the student's rotation grid by angle. In 3D a mirror
   pair's targets follow the readout's rule (the unmirrored direction's rotation, or the reflected one mapped through
   the pair's frame angle); Phase 1's tester pins the mapping at 0, 75 and 180°, as the HHCNN tester pins the readout.
   Made precise in Phase 1 (ruling R2): a reference's targets are its marginal mass in the row times its in-plane
   posterior (the joint posterior over class and rotation), and a pair's two members go into one class's targets,
   the unmirrored on the map's rotations, the mirror on the reflection's.
4. **Variants.** A, the base. M, marginal targets (decision 2) and a log-sum-exp readout at a learned temperature,
   initialised to the base's scale. R, M plus a rotation loss: the cross-entropy between the student's in-plane
   softmax and decision 3's targets at the same 8 references, weight 1 (guessed). M and R are read apart so that each
   gain is attributed.
   Made precise in Phase 1 (ruling R2): the learned temperature is the readout's existing scale, so a class's logit is
   ln Σ_r exp(scale × s_c(r)) and the class posterior is the marginal of a softmax over (class, rotation); with
   decision 3's joint targets, R's loss is the cross-entropy of that joint softmax at the 8 references, weight 1 per
   unit of the teacher's mass. The early stopping reads the held-out class cross-entropy in M and R alike (the
   base's schedule, decision 9). The variants differ in two constants of `simple_defs_hhcnn` (the marginal with its
   readout; the rotation loss), the patches otherwise identical.
5. **Teachers**, rebuilt on the base in Phase 0 and frozen with checksums, three training seeds each:
   - 3D: the medium truth set's `T200`, `T1000`, `T5000`, `T` (24,000) and bgal;
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

   The references, as Phase 0 measured them (Progress, Phase 0, "What each read is"): in 3D the readout's pose
   (`predict_poses`' Euler angles) against the truth pose on the test half of the truth sets, and against the
   teacher's mode rotation of the kept rows (`student_rows.txt`) on bgal's training rows; in 2D the held-out in-plane
   error is `learn_class_student`'s own (its log), and a prediction read is added on the frozen teacher's training
   particles (`predict_classes` against the teacher's stored class and angle). "The share seeded" is, in 3D, the
   prediction's route-1 share (bgal: the training rows; truth sets: the test half), the read the pose student's 0.42
   in section 1 refers to.
   Revised in Phase 2, before any variant was read (ruling R2): M and R change what the rows' mode rotation is (the
   mode reference's in-plane mode, no longer one drawn rotation), so a read against it would judge each variant
   against its own reference. The in-plane reads that decision 7 compares across variants are therefore taken
   against the stored poses, one reference for all: on bgal in 3D, the training rows' read-out pose against the
   training set's stored pose (the teacher map's refinement), agreeing directions within 6°; in 2D, the held-out
   rows whose student class is the teacher's stored class, |student e3 - stored e3|. A's values, from Phase 0's kept
   predictions (`R/scratch/keep/fixed_inpl.tsv`): bgal 3D 0.866, 0.880, 0.883° (seeds 1 to 3; 1,251 to 1,258 agreeing
   rows), mean 0.876; bgal 2D 0.934, 0.839, 0.928° (188 to 198 rows), mean 0.900. The reads against each run's own rows
   stay in the tables, reported. The truth sets' in-plane read was already against a fixed reference (the truth).
7. **Verdict, written now.** M (or R) is adopted for a student only if:
   - (a) no held-out read's mean over the three seeds falls below A's range over its seeds, on any of that student's
     sets;
   - (b) the share seeded at the calibrated target, or within 6° (3D) / classes recovered (2D), rises above A's
     range by more than twice that range, on at least three of the five 3D sets for the pose student, on both sets
     for the class student (margin guessed, re-derived from Phase 0's spread before Phase 2);
   - (c) for R over M: the in-plane error falls by more than twice M's range with (a) holding against M.
   Otherwise the variant is rejected; the record goes to `rejected/` (Phase 4).
   **Replaced by ruling R4** (Hans, 2026-11-15, after Phase 3's verdict, applied here for the first time): "accept
   unless worse". The variant is R, for both students (M is contained in it); (b) no longer applies. First an audit
   that the teacher's in-plane weighting is the refinement's (a difference is a bug, fixed in R, and R's three seeds
   re-run on every set). Then "worse" is read beyond the noise of three seeds: the reads where R fell below A's range
   are re-measured at five more seeds (4 to 8) for A and R, and R is worse on a read when the mean of the eight paired
   differences (R minus A at one seed) lies in the worse direction by more than twice its standard error. Not worse on
   any read: adoption. Worse: look for an implementation error, else stop with the paired table.
   **Ruling R5** (Hans, 2026-11-15, applied here for the first time): the run's base moves to `4282d58a3` (master
   after the day's fixes to the in-plane insertion's collapsed quadrature, the Euler poles and `predict_poses`' seed)
   before R4 goes on. The corrected variants are carried onto the new base, its off-state identity and the audit are
   redone there, and A and R are re-measured on it: seeds 1 to 3 on every set, then the paired seeds 4 to 8 on the sets
   with a read below A's range (R4's list as the start). Phases 0 to 3 and R4's first audit stay on the old base and
   are reported as such; the frozen 3D sets and 2D teachers stay.
   **Ruling R6** (Hans, 2026-11-15, after R4's stop, applied here for the first time): the class student learns the
   physical posterior, `prob_beta` dropped from its teacher; the pose student keeps A (3D's in-plane insertion is off
   by default on the new base `603c44490`). Variant P is R's teacher and network for the class student with the
   class distribution read at the physical temperature too, and `validate_class_student`'s reference becomes the
   physical posterior's best class. A and P are read on both 2D sets at eight paired seeds by R4 item 3's rule: P is
   adopted for the class student unless it is worse on a read; if it is, an implementation error is looked for, then
   the run stops with the paired table.
   **Ruling R7** (Hans, 2026-11-15, after R6's stop, applied here for the first time): "adopt P for the class
   student". The four reads R6's paired test flagged (shortlists holding the physical reference's class) are
   accepted as the calibration's trade: P holds the 95% target where A is over-cautious, and shortlists more
   particles. P goes into the repository for the class student only; the pose student keeps A.

   **(b)'s margin re-derived from Phase 0's spread** (2026-11-15, before any variant was built or read; the rule as
   written, not changed): a read passes (b) on a set when the variant's three-seed mean exceeds A's maximum over its
   seeds by more than twice A's range, so the bar is A's max + 2 × (A's max − A's min). (a)'s floors are A's ranges in
   Progress, Phase 0 (for the reads where lower is better, the in-plane error, the distance and the cone size, "below
   A's range" means above A's maximum).

   | student | set | read | A's range over seeds 1-3 | (b)'s bar |
   | --- | --- | --- | --- | --- |
   | pose | `T200` | seeded / within 6° | 0.1158-0.3003 / 0.6904-0.7922 | 0.6693 / 0.9958 |
   | pose | `T1000` | seeded / within 6° | 0.2266-0.2611 / 0.7562-0.8228 | 0.3301 / 0.9560 |
   | pose | `T5000` | seeded / within 6° | 0.3678-0.4812 / 0.8363-0.8919 | 0.7080 / 1.0031 (unreachable) |
   | pose | `T` | seeded / within 6° | 0.4671-0.5840 / 0.9125-0.9458 | 0.8178 / 1.0124 (unreachable) |
   | pose | bgal | seeded / within 6° | 0.3873-0.4540 / 0.8720-0.8793 | 0.5874 / 0.8939 |
   | class | bgal | held-out classes recovered / held-out seeded | 0.9673-0.9818 / 0.9673-0.9818 | 1.0108 / 1.0110 (both unreachable) |
   | class | bgal | `validate_class_student`: particles given a shortlist | 71.4-80.5% | 98.7% |
   | class | trpv1 | held-out classes recovered / held-out seeded | 0.8530-0.8780 / 0.8560-0.8790 | 0.9280 / 0.9250 |
   | class | trpv1 | `validate_class_student`: particles given a shortlist | 67.1-72.6% | 83.6% |

   What the bars imply, recorded for Hans and not acted on (the margin moves only by this rule): on bgal the class
   student's held-out reads sit at their ceiling, so (b) can pass there only through the validation's share of
   particles given a shortlist, at 98.7%; the pose student can pass (b) on `T5000` and `T` only through the share
   seeded.
8. **Early stop.** Phase 2 reads M and R on `T1000`, bgal (3D) and bgal (2D). If neither moves any read of (b) beyond
   A's range, the experiment stops as a negative result.
9. **One schedule.** All variants train with the base's schedule, designs and shift augmentation; the class
   student at width 16 (`CLASS_STUDENT_NDESIGNS`), the pose student with its design choice.

Defaults taken without a question: `KMAX_INPL_ROW` = 8; the rotation loss's weight 1; the log-sum-exp temperature
learned per student, one scalar.

## 4. Not in this plan

- The insertion, the tables' draws and the refinement's use of the cones or shortlists.
- The stream (`planned/class_student_stream.md`).
- Any file outside the table.

## 5. Change types and gates

| Phase | Kind | Gate |
| --- | --- | --- |
| 0 | measurement | the base's reads at three seeds, frozen |
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

- **Phase 0, baseline.** Pin the base. Release and Debug builds, fast gate. Rebuild
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
| `src/main/strategies/search/probabilistic/simple_eul_prob_marginal.f90` (new; tested in `simple_eul_prob_tab_tester.f90`) | the marginal's per-entry work both tables call, on the calculator (added under ruling R3) | 1 (export), 4 |
| `src/main/learn/simple_hhcnn.f90`, `simple_hhcnn_bench.f90` (and tester) | the log-sum-exp readout, the rotation loss | 1 (export), 4 |
| `src/main/learn/simple_polar_student.f90`, `simple_pose_student.f90`, `simple_class_student.f90` (and testers) | targets, training, the model files' versions | 1 (export), 4 |
| `src/main/strategies/parallelization/simple_strategy_pose_student.f90`, `simple_strategy_class_student.f90`, `simple_strategy_student_common.f90` (and testers) | the learn programs ask for the rows' extension | 1 (export), 4 |
| `src/defs/simple_defs_hhcnn.f90` | `KMAX_INPL_ROW` and the variant's switches (added under ruling R3) | 1 (export), 4 |
| `doc/code_overview/code_base_map.md` | regenerated by the build from the changed and new `!@descr` lines (added under ruling R3) | 1 (export), 4 |
| `doc/algorithms/pose_student.md`, `doc/algorithms/class_student.md` | the teacher and the readout | 4 |
| `src/main/commanders/simple/simple_commanders_class_student.f90` | the physical teacher's `prob_beta`, `validate_class_student`'s reference (added under rulings R3 and R6) | 4 |
| `src/main/commanders/test/simple_commanders_test_highlevel_class_student.f90` | the `class_student` workflow test reads the reference's renamed timing line (added under rulings R3 and R7) | 4 |
| `doc/implementation_notes/planned/student_inplane_posterior_rejected.md` | the rejected record, if both variants are rejected (the maintainer moves it to `rejected/`); deleted on the adoption (ruling R7) | 4 |
| `doc/implementation_notes/README.md` | record (this plan and its report are always allowed) | 0, 1, 2, 3, 4 |

## Progress

### Phase 0, the baseline (2026-11-15, the Dell, done)

**The base and the machine.** HEAD `86babf86b66c26a6e521ca61309e3f78c3745f0c` with the launch's uncommitted
documentation (this plan, untracked, and five documentation files: `doc/algorithms/class_student.md`,
`doc/implementation_notes/README.md`, `completed/class_student_vs_reference_matching.md`,
`planned/class_student_stream.md` and a planned note since deleted, 29 lines added and 19 removed); no
source file differs from HEAD, so the source the run builds is `git archive 86babf86b` (ruling R1: the run works on
this base and does not follow master). The Dell: Oracle Linux 8.10, Intel Core Ultra 9 285 (24 cores, no
hyperthreading), 62 GB, GCC 15.2 of gcc-toolset-15, CMake 3.26. Run root `/data2/agent_runs/student_inplane`
(`R` below); every path under `R/scratch/keep/` unless absolute.

**Builds.** Both configured as `cmake -S src_base -B <tree> -DCMAKE_BUILD_TYPE=<Release|debug> -DBUILD_TESTS=ON
-DNICE=OFF -DUSE_ARCHOPT=ON -DFORTLS=OFF -DPython3_EXECUTABLE=/usr/bin/python3.12`, `make -j20`, `make install`
(`tools/build_base.sh`): `build_base` (Release, make 83 s) and `build_base_debug` (Debug, 41 s). The fast gate
(`scripts/run_fast_gate.sh`), 15 of 15 on both (11 s and 33 s; `R/scratch/logs/gate_env_*.log`). The login shell
of the Dell points `SIMPLE_PATH` and `PATH` at the maintainer's `~/src/SIMPLE/build`, so a program that starts
workers (here the identity runs' motion correction) ran the maintainer's binaries until every script sourced
`tools/env.sh <tree>`; the first identity start was killed and discarded, and the gates were re-run under it.

**Load.** Every timed job ran alone at 24 threads; the one-minute load during them was 15 to 23, which is the job's
own threads (nothing else heavy ran; four single-thread jobs beside a 24-thread training had made it 30 times slower,
so no timed job overlapped another).

**Off-state identity, base against itself.** `simulated_workflow` (the end-to-end simulated test: refine2D inside
solve2D, then solve3D's refine3D stages, one thread, its own fixed seed), both suites, twice on `build_base`
(`identity/base1`, `identity/base2`; 1jxy 16:32 and 15:37, 6vxx 14:41 and 14:07); `test=identity` of `base1`
against `base2`: 1jxy 223 files the same, 0 differing; 6vxx 186 and 0 (`identity/cmp_base1_base2/<suite>/`). The
exclusion lists are empty for both suites. `base1` stays for Phases 1 and 4.

**Thread identity of the training** (`thread/`): at a fixed thread count the training is reproducible (pose student,
`T1000` seed 1, 24 threads twice: the same weights; class student, bgal seed 1, 24 threads twice: the same weights).
Between one and 24 threads it is not, on the base, for both students: the weights differ (pose: the teacher's rows,
the sigma2 and the reprojection models are byte-identical and the feature cache differs only in its header; epoch 1
agrees, epoch 2's training loss is 1.71403 against 1.71404, the held-out loss 1.60900 against 1.60840, the best epoch
3 in both). This is the finding the class-student run recorded in its Phase 2 (`T200`, left for Hans); the network
tester's team checks pass. The plan's gate is the tester's thread identity (section 5), so Phase 1 holds the
variants to that and to the base's behaviour at the program level (the same at a fixed thread count).

**The teachers.**
- 3D: the fixed sets of `student_speed`, read in place (`/data2/agent_runs/student_speed/scratch/keep/sets/`); their
  `MD5SUMS` checked: every file matches its latest line (`T200` a91c2eb3, `T1000` 70575505, `T5000` 30c7998f, `T`
  62ce365b with `tmap.mrc` 85008090, the centred test half `S` 4f75a6f6, bgal `training_set` b353e24e,
  `all_particles` 6c76a170, `teacher_map` babde2eb); the one stale line is `S.simple`'s first, the copy without
  shifts now kept as `S_zero_shifts.simple`. The teacher's rows are computed inside `learn_pose_student` by the base
  build; the rows' fill takes 0.1, 0.5, 2.6 and 56 s on `T200`, `T1000`, `T5000` and `T`, 0.8 s on bgal.
- 2D (`teachers/{bgal,trpv1}`, `tools/teachers2d.sh`): the import projects, a random subset by
  `selection nran=2756|10000 prune=no` under `SIMPLE_SEED=1` (the other particles stay unseen, state 0), then
  `solve2D ncls=50|100 mskdiam=180|175 nthr=24` with every default (the in-plane posterior insertion and the sigma2
  annealing). bgal: 33.1 s, 30 iterations, resolution 12.02 Å (the frozen class-student teacher: 35.4 s, 11.18 Å);
  trpv1: 136.4 s, 30 iterations, 12.18 Å (139.3 s, 13.03 Å). Frozen with `MD5SUMS`: bgal project 33efdf45,
  `cavgs_iter030.mrc` 31004589 (even 41f5e739, odd 3aac5b54), `sigma2_state.bin` 90b7d7ee; trpv1 project a04dcef4,
  `cavgs_iter030.mrc` 93268ef3 (even 76da7232, odd 474193f6), `sigma2_state.bin` 23e5c9f3.

**What each read is** (the tools in `tools/`; they drive SIMPLE's programs and read their outputs).
- 3D (`harness.sh`, student_crop's harness): `learn_pose_student` on each set (truth sets against the training half's
  map, `lp=15 nspace=1000 shift_sigma=2 mskdiam=300 pgrp=c1`; bgal against its teacher map, `mskdiam=180 pgrp=d2`),
  `predict_poses` on the centred test half `S` (24,000 unseen particles with their truth poses) or, for bgal, on the
  training set and on all particles. Within 6°, truth in the cone, mode = truth and the like are the prediction
  table against the stored pose; "seeded" is the share the prediction seeds (route 1: its confidence at or above
  the threshold calibrated on the held-out rows); held-out recovered and seeded are the model header's.
- The 3D in-plane error (`inpl_read.py`): the readout's pose against a reference pose, on the particles whose view
  direction agrees within 6°, as the angle between the two images' x axes after carrying one view axis onto the
  other; the truth sets against the truth pose of `S` (no teacher rows exist for the test half, so the truth's
  in-plane angle stands in for the teacher's mode rotation), bgal against the teacher's mode rotation of the kept
  rows on every training row and on the held-out rows. Checked before use: the base's earlier `T` prediction against
  the truth, median 0.83°; the teacher's mode rotations of `T1000` against the truth, median 1.125°, one step of the
  teacher's grid (360°/320).
- 2D (`learn2d.sh`): `learn_class_student` on the frozen teacher (`refs=` its final averages, `lp=15`); held-out
  reads from the model header and the log (top-1, seeded, the teacher's class among a seeded row's candidates,
  their product "classes recovered", the threshold, the in-plane error on rows of agreeing class); then
  `predict_classes` on a copy of the teacher's project (its training particles at the teacher's shifts) against the
  teacher's stored class and angle.
- `validate_class_student` (`validate2d.sh`, the command lines of the reference-matching run: `nsample=2756 ncls=50
  mskdiam=180` on bgal, `nsample=10000 ncls=100 mskdiam=175` on trpv1, `nthr=24`), its own solve2D teacher inside,
  seeds 1 to 3; reference matching beside it at seed 1 (`refmatch2d.sh`: one `refine2D refine=prob prob_beta=100`
  pass of the unseen particles at their picked shifts against the validation teacher's averages at the working box
  of 120 px, then `make_cavgs`).

**The pose student, A at seeds 1 to 3** (mean [min, max]; `summary.tsv`, `spread3d_base_p0.txt`, the runs in
`runs/base_p0/`):

| set | within 6° | truth in the cone | seeded | held-out recovered | held-out seeded | in-plane error (°) | agreeing share |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `T200` | 0.743 [0.690, 0.792] | 0.968 [0.962, 0.973] | 0.202 [0.116, 0.300] | 1.000 [1.000, 1.000] | 0.550 [0.500, 0.650] | 0.97 [0.92, 1.01] | 0.743 [0.690, 0.792] |
| `T1000` | 0.797 [0.756, 0.823] | 0.981 [0.974, 0.986] | 0.245 [0.227, 0.261] | 1.000 [1.000, 1.000] | 0.500 [0.470, 0.540] | 0.93 [0.89, 0.96] | 0.796 [0.755, 0.822] |
| `T5000` | 0.870 [0.836, 0.892] | 0.984 [0.978, 0.989] | 0.434 [0.368, 0.481] | 0.997 [0.996, 0.997] | 0.596 [0.554, 0.626] | 0.89 [0.86, 0.91] | 0.870 [0.836, 0.892] |
| `T` | 0.930 [0.912, 0.946] | 0.995 [0.991, 0.996] | 0.529 [0.467, 0.584] | 0.999 [0.998, 1.000] | 0.535 [0.474, 0.587] | 0.82 [0.80, 0.83] | 0.930 [0.912, 0.946] |
| bgal | 0.876 [0.872, 0.879] | 0.959 [0.957, 0.961] | 0.426 [0.387, 0.454] | 1.000 [1.000, 1.000] | 0.349 [0.313, 0.380] | 0.70 [0.70, 0.70] | 0.854 [0.847, 0.857] |

| set | mode = truth | within 12° | median distance (°) | median cone | seeded in the cone | unseen seeded | in-plane q90 (°) | held-out rows' in-plane (°) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `T200` | 0.237 [0.216, 0.261] | 0.839 [0.800, 0.870] | 3.79 [3.51, 4.08] | 17 [9, 29] | 0.982 [0.981, 0.983] | - | 2.45 [2.32, 2.56] | - |
| `T1000` | 0.253 [0.235, 0.263] | 0.894 [0.869, 0.912] | 3.58 [3.46, 3.78] | 12 [11, 13] | 0.987 [0.979, 0.993] | - | 2.36 [2.27, 2.43] | - |
| `T5000` | 0.281 [0.264, 0.291] | 0.952 [0.934, 0.963] | 3.31 [3.21, 3.47] | 6 [6, 7] | 0.985 [0.981, 0.988] | - | 2.22 [2.16, 2.28] | - |
| `T` | 0.307 [0.295, 0.320] | 0.993 [0.991, 0.995] | 3.10 [3.02, 3.18] | 5 [5, 6] | 0.992 [0.987, 0.995] | - | 2.07 [2.00, 2.12] | - |
| bgal | 0.172 [0.166, 0.175] | 0.928 [0.924, 0.933] | 2.48 [2.46, 2.51] | 20 [19, 20] | 0.997 [0.997, 0.998] | 0.062 [0.043, 0.082] | 1.96 [1.90, 2.06] | 0.75 [0.70, 0.79] |

Training at 24 threads 8.2 [7.9, 8.7] s on `T200`, 18.4 [17.9, 18.9] on `T1000`, 86.4 [82.5, 90.3] on `T5000`,
368.8 [358.0, 380.2] on `T`, 23.7 [23.5, 23.9] on bgal; prediction 8 to 16 s for the 24,000 particles of `S`, 2 to
3 s for bgal's 2,756. The reads of `T200`, `T1000` and bgal agree with the student_crop run's base to the third
decimal (`T200` seed 1: 0.690, 0.962, 0.116). bgal's in-plane median is 0.70° at every seed: both angles lie on
rotation grids, so the median takes few values; its spread is zero and (c) of decision 7 reads against M's range.

**The shift check** of `predict_poses` (256 particles moved by 2, 5 and 10 px in four directions; reported, not
judged), mean over the seeds [min, max]; the truth sets' student on `S`, bgal's on all particles
(`shiftcheck_base_p0.txt`, each prediction's `*_shiftcheck.txt` in its run):

| set | px | mode kept | cone holds the mode | cones' Jaccard index | posterior TV (mean) |
| --- | --- | --- | --- | --- | --- |
| `T200` | 2 / 5 / 10 | 0.206 / 0.006 / 0.003 | 0.919 / 0.279 / 0.163 | 0.295 / 0.079 / 0.046 | 0.854 / 0.985 / 0.993 |
| `T1000` | 2 / 5 / 10 | 0.224 / 0.006 / 0.003 | 0.929 / 0.300 / 0.168 | 0.296 / 0.073 / 0.042 | 0.834 / 0.985 / 0.993 |
| `T5000` | 2 / 5 / 10 | 0.268 / 0.006 / 0.002 | 0.926 / 0.243 / 0.145 | 0.284 / 0.040 / 0.024 | 0.847 / 0.994 / 0.997 |
| `T` | 2 / 5 / 10 | 0.371 / 0.003 / 0.004 | 0.968 / 0.208 / 0.146 | 0.259 / 0.023 / 0.015 | 0.800 / 0.996 / 0.997 |
| bgal | 2 / 5 / 10 | 0.720 / 0.216 / 0.007 | 0.998 / 0.955 / 0.266 | 0.781 / 0.409 / 0.089 | 0.259 / 0.793 / 0.986 |

The per-seed ranges are within 0.04 of these means (`shiftcheck_base_p0.txt`).

**The decentring curve** of the `T` students (reported): `S` decentred by exactly r px
(`sets_extra/S_decentred/`, built with SIMPLE's `new_project`, `import_particles` and `selection` from an orientation
table whose truth shifts are moved by r along a seeded direction per particle, the pose_student run's recipe;
checked row by row against `S.simple`: r = 0 identical, every r moved by r ± 0.0001 px with every other key equal;
`MD5SUMS` beside them), mean [min, max] over the seeds (`decentring.tsv`):

| r (px) | 0 | 1 | 2 | 3 | 5 | 8 |
| --- | --- | --- | --- | --- | --- | --- |
| truth in the cone | 0.995 [0.991, 0.996] | 0.991 [0.988, 0.994] | 0.955 [0.951, 0.959] | 0.710 [0.704, 0.716] | 0.216 [0.213, 0.217] | 0.153 [0.149, 0.157] |
| mode within 12° | 0.993 [0.991, 0.995] | 0.981 [0.978, 0.985] | 0.848 [0.834, 0.866] | 0.394 [0.373, 0.415] | 0.048 [0.047, 0.049] | 0.030 [0.029, 0.030] |

**The class student, A at seeds 1 to 3** (`summary2d.tsv`, `spread2d_base_p0.txt`, runs in `runs2d/base_p0/`), on
the frozen teachers. Held-out rows (every tenth selected particle):

| set | top-1 | seeded | recovered among seeded | classes recovered | threshold | in-plane (°) |
| --- | --- | --- | --- | --- | --- | --- |
| bgal | 0.661 [0.651, 0.665] | 0.975 [0.967, 0.982] | 1.000 [1.000, 1.000] | 0.975 [0.967, 0.982] | 0.9900 [0.9900, 0.9901] | 0.74 [0.70, 0.82] |
| trpv1 | 0.486 [0.475, 0.497] | 0.869 [0.856, 0.879] | 0.997 [0.997, 0.999] | 0.866 [0.853, 0.878] | 0.9900 [0.9900, 0.9900] | 0.92 [0.85, 0.97] |

The prediction on the teacher's training particles, against its stored class and angle:

| set | top-1 | seeded | stored class among a seeded particle's candidates | among every particle's | median candidates | in-plane (°) |
| --- | --- | --- | --- | --- | --- | --- |
| bgal | 0.714 [0.709, 0.717] | 0.954 [0.941, 0.962] | 0.991 [0.989, 0.993] | 0.989 [0.988, 0.992] | 3 [3, 3] | 0.89 [0.82, 0.95] |
| trpv1 | 0.647 [0.640, 0.652] | 0.894 [0.880, 0.903] | 0.996 [0.996, 0.997] | 0.993 [0.992, 0.994] | 7 [6, 7] | 1.20 [1.17, 1.21] |

Training (width 16 only) 11.3 [11.2, 11.4] s on bgal, 58.9 [43.5, 71.7] s on trpv1, of it the rows' fill 0.33 and
2.87 s; prediction 2.2 and 5.5 s.

**`validate_class_student`** at seeds 1 to 3 (each seed its own training subset and solve2D teacher;
`summary_validate.tsv`, `spread_validate_base_p0.txt`, reports in `validate/base_p0/<set>_s<seed>/kept/`, the seed-1
class averages of the teacher, the full search and the student beside them for inspection):

| set | shortlisted (%) | shortlist holds the full search's class (%) | same, perfectly centred (%) | shortlist length | class FRC 0.5 student / full search (Å) | top-1 (%) | in-plane (°) | held-out seeded | total (s) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| bgal | 75.6 [71.4, 80.5] | 76.3 [73.4, 79.0] | 87.0 [83.2, 90.1] | 5 [4, 5] | 17.98 [16.32, 20.40] / 14.62 [14.19, 14.84] | 32.8 [32.1, 33.9] | 1.01 [0.94, 1.05] | 0.969 [0.964, 0.975] | 75 [74, 76] |
| trpv1 | 70.1 [67.1, 72.6] | 90.9 [90.5, 91.7] | 96.7 [96.5, 97.1] | 7 [7, 7] | 22.02 [20.75, 23.08] / 18.01 [16.83, 18.88] | 36.9 [36.1, 37.8] | 1.14 [1.07, 1.24] | 0.865 [0.858, 0.876] | 352 [337, 365] |

Its times (s): solve2D 37.4 / 146.2, learning 12.1 / 72.2, `center_ptcls` 14.9 / 53.5, prediction 1.63 / 10.68,
the full search 6.1 / 51.7; cost per particle, student against full search, 3.29 against 2.21 ms (bgal) and 1.92
against 2.01 ms (trpv1).

**Reference matching** beside the seed-1 validations (reported; `refmatch/base_p0/<set>_s1/reads.txt`): the same
class as the full search for 87.0% of bgal's 2,757 unseen particles (the student 33.9%) and 87.4% of trpv1's 25,645
(36.1%), in 1.9 and 12.9 s; class FRC 0.5 against the teacher's, matching / full search / student: bgal 14.84 /
14.84 / 20.40 Å, trpv1 22.23 / 18.31 / 22.23 Å.

**Decision 7 (b)'s bars** re-derived from these spreads by the rule as written, before any variant: in decision 7.
On bgal the class student's held-out reads are at their ceiling (bars 1.011), so (b) can pass there only through the
validation's share of particles given a shortlist (bar 98.7%).

**Exit.** The reads and spreads frozen above (one table per student, mean and range over the seeds per read and
set); decision 7's margins written; the teachers' checksums recorded; the off-state identity of the base holds with
empty exclusion lists. Files changed: this plan (Progress, the base pin, decision 6's references, decision 7's
bars). Logs:
`R/scratch/logs/` (`build_base*.log`, `gate_env_*.log`, `identity_base*.log`, `thread3d.log`, `chain_p0.log`), the
tools and tables under `R/scratch/keep/`.

### Phase 1, the teacher and the readout (2026-11-15, the Dell, done)

**What was built, outside the tree.** The variants are copies of the base export with a patch each
(`R/scratch/keep/variants/M.patch`, `R.patch`, a `diff -ruN` of `src_base` against the variant's source; 17 files,
about 1,500 lines added). They differ in one constant: R is M with the rotation loss switched on. Release builds are
in `variants/{M,R}/build`. The repository's sources are untouched. The patch, file by file (lines added / removed):
- The teacher (the probability tables): the new module `simple_eul_prob_marginal.f90` (+173), `simple_eul_prob_tab.f90`
  (+107 / -14), `simple_eul_prob_tab2D.f90` (+109 / -8), `simple_eul_prob_tab_utils.f90` (+110 / -17). A learn program
  calls `set_inpl_marginal`. The fill then computes each entry's in-plane row once, draws the same rotation from it
  with the same keyed variate (the entry unchanged, bit for bit) and keeps the row's log-partition sum `lnz` and its
  minimum. A refined entry is rescored at its refined shift, and the in-plane posteriors of the 8 heaviest references
  are kept. `get_posterior_row` reads the marginal (decision 2 as made precise). The rows file is at version 2: the
  record, the training losses, `lnz` and the minimum, and the in-plane rows.
- The network, `simple_hhcnn.f90` (+168 / -31): the log-sum-exp readout as a construction flag beside the
  reflection, its backward pass, and the optional rotation targets of R. The weights file is at version 3.
- The students: `simple_polar_student.f90` (+67 / -13; the readout recorded in the model header and checked on load,
  the rotation targets through training, `teacher_rot_bin`), `simple_pose_student.f90` (+87 / -14; the pair's targets,
  `readout_e3`, model version 5), `simple_class_student.f90` (+51 / -9; model version 2).
- The strategies: `simple_strategy_pose_student.f90` (+38 / -6) and `simple_strategy_class_student.f90` (+39 / -6) ask
  for the marginal, write the rows file and, in R, pass the in-plane rows (rows header `rows_version=2`).
- The constants, `src/defs/simple_defs_hhcnn.f90` (+5): `KMAX_INPL_ROW` = 8, `STUDENT_INPL_MARGINAL` (M and R) and
  `STUDENT_INPL_ROTLOSS` (R).
- Testers: `simple_hhcnn_tester.f90` (+193 / -41), `simple_eul_prob_tab_tester.f90` (+278 / -4),
  `simple_pose_student_tester.f90` (+62 / -1), `simple_polar_student_tester.f90` (+3 / -2), and the benchmark's
  constructor calls (+3 / -3).
- `doc/code_overview/code_base_map.md` (+4 / -3), regenerated by the build from the changed `!@descr` lines.
Three files came in under ruling R3, each given its row in the file table before its first edit: the constants module,
the new marginal module and the generated map. Decisions 2, 3 and 4 were made precise under ruling R2, in their own
text: the temperature, the refined entries, the joint targets, the readout's temperature being its scale, and the
early stopping on the class cross-entropy.

**The checks** (each new one written with its violation planted first; 18 plants in `R/scratch/keep/plants/`,
`plants.log` holds every failure message). 17 of the 18 were caught. The one that was not, a thread-dependent change
of 1e-6 in the scale's gradient, was too weak to test anything, since Adam's first steps do not react to it. It was
replaced by a thread-dependent change of the scores' gradient, which the check catches at both teams:
- The teacher's marginal against a brute-force sum (`simple_eul_prob_tab_tester.f90`, `test_inpl_marginal`, on a
  calculator of random references and particles). An entry drawn through the marginal's route is the fill's own draw,
  bit for bit, on the shifted-reference kernel and on the bank. The row's `lnz` and minimum equal a brute-force sum
  over every rotation in double, also at a refined shift. The kept references are those of largest marginal mass,
  with their brute-force in-plane posteriors. `marginal_dist`'s tempered weights are the marginal masses. At a team
  of four the entries and marginals equal those at one thread, bit for bit. The rows file at version 2 round-trips
  exactly (`test_row_file_marginal`).
- The mirror mapping of decision 3 (`simple_pose_student_tester.f90`, `test_rotation_targets`). At frame angles of 0,
  75 and 180°, for the unmirrored and the reflected member, on a grid as fine as the teacher's and on one a third as
  fine, a teacher rotation mapped onto the readout's grid reads out the teacher's angle within half a grid step. A
  pair's two members go into one class's targets, the reflected one on the reflection's rotations. The plant with the
  wrong sign of the frame is caught at 75° (104 mismatches).
- The readout (`simple_hhcnn_tester.f90`):
  - Each class's logit equals a brute-force log-sum over its rotations (and reflections), and each class's rotation
    posterior sums to one.
  - The class posterior is unchanged by a rotation of one grid step, by conjugation and, with the reflection, by the
    calculator's mirror at 0, 75 and 180°.
  - The gradient check passes with the log-sum-exp readout, without and with rotation targets, with and without the
    reflection.
  - The thread identity of the training step holds with the log-sum-exp readout and rotation targets at teams of 4
    and 16, bit for bit.
- A finding about the gradient check itself. Its central differences fail on some inputs, on either readout alike,
  across a modulus near zero inside the step. Over ten fixed seeds (`R/scratch/keep/phase1_evidence/gradcheck_seeds.txt`):
  - the base's max readout fails on 1 of 10 per reflection setting;
  - the log-sum-exp readout fails on 1 of 10;
  - the log-sum-exp readout with rotation targets fails on 2 of 10, with 1 to 3 entries of 359 each time;
  - seeds 3, 6 and 7 fail for both readouts; the other seven pass every configuration.
  The base's cases keep their inputs. The new cases run on the fixed seed 20261101, and no tolerance was changed.
- Suites on both variant builds (Release): `unit_numerics` 1,905 assertions, `unit_pftc_align2D3D` 6,317, all passing;
  the fast gate 15 of 15. M's Debug build (bounds checks, the Debug arithmetic path): both suites pass.
- An independent read of the patch found no error and no effect on the off state. It named latent traps, now fixed:
  - the row is read at the fill's `prob_beta` or refused;
  - a table with the marginal switched on but never filled is refused;
  - the in-plane mode is left as the draw for a particle without a temperature;
  - the heaviest references are ranked by the marginal loss the row is read with;
  - copying weights between readouts is refused;
  - the rotation loss without the marginal is refused;
  - the tester's calculator and parameters are allocatable objects, and the global thread count is restored.

**Off-state identity** (refine2D inside solve2D and refine3D at one thread, the suites' fixed seed;
`simulated_workflow` on each variant's final build against `base1`, `test=identity`): M 1jxy 223 files the same, 0
differing; M 6vxx 186, 0; R 1jxy 223, 0; R 6vxx 186, 0 (`R/scratch/keep/identity/cmp_base1_{M,R}/`). The builds
before the review's fixes gave the same (`cmp_base1_{M,R}_build1`). No exclusion was needed.

**End to end** (`R/scratch/logs/smoke_p1.log`; not reads, no line in the summaries):
- `learn_pose_student` and `predict_poses` ran on `T200`, and `learn_class_student` and `predict_classes` on the bgal
  teacher, on both variants.
- Every row of the version-2 rows files passed `check_rows_v2.py`: the training loss is the marginal's, the kept
  references are the heaviest in order, each in-plane row sums to one, and the record's mode and in-plane mode are
  the row's.
- R trained on the whole set `T` (24,000 rows) at a peak of 2.51 GB for the process, against 1.72 to 1.79 GB for the
  base in Phase 0. The rows' in-plane extension is about 0.75 GB: 0.25 GB of in-plane rows and 0.49 GB of targets.
  That is within the planned 1 GB and below section 7's stop at 2 GB.
- Its teacher's fill took 55 s, against the base's 56 s.

**Exit.** The testers pass on both variants (with the plants recorded) and the off-state identity holds. Files
changed in the repository this phase: this plan only (the precisions of decisions 2 to 4, three file-table rows under
ruling R3, this entry). Logs: `R/scratch/logs/` (`chain_p1_final.log`, `identity_MR.log`, `make_variants.log`,
`smoke_p1.log`), `R/scratch/keep/variants/{M,R}/tests/`, `R/scratch/keep/plants/`.

### Phase 2, the first read (2026-11-15, the Dell, done: decision 8 says continue)

**What ran** (`R/scratch/logs/chain_p2.log`, `R/scratch/keep/tools/chain_p2.sh`): M and R, each at seeds 1, 2 and 3,
with Phase 0's commands for A. In 3D, the harness on `T1000` and on bgal (`learn_pose_student`, then `predict_poses`).
In 2D, `learn_class_student` on the frozen bgal teacher with its prediction read, and `validate_class_student` on
bgal. One job at a time, 24 threads; the one-minute load was 15 to 22, which is the jobs' own threads. The lines
are in `summary.tsv`, `summary2d.tsv` and `summary_validate.tsv` under the builds `M` and `R`, phase 2. The comparison
with A, every read, is in `R/scratch/keep/phase2_compare.txt` (`tools/compare.py`).

**Decision 8** reads M and R on these sets: the experiment stops if neither moves any read of decision 7 (b) beyond
A's range. Mean over the three seeds [range]; A's range from Phase 0:

| set | read of (b) | A | M | R |
| --- | --- | --- | --- | --- |
| `T1000` (3D) | seeded | 0.227-0.261 | 0.273 [0.245, 0.295], above | 0.303 [0.254, 0.376], above |
| `T1000` (3D) | within 6° | 0.756-0.823 | 0.774 [0.732, 0.817], inside | 0.778 [0.742, 0.817], inside |
| bgal (3D) | seeded | 0.387-0.454 | 0.403 [0.353, 0.447], inside | 0.492 [0.428, 0.576], above |
| bgal (3D) | within 6° | 0.872-0.879 | 0.864 [0.842, 0.877], below | 0.888 [0.869, 0.909], above |
| bgal (2D) | held-out classes recovered | 0.967-0.982 | 0.965 [0.949, 0.975], below | 0.979 [0.967, 0.986], inside |
| bgal (2D) | held-out seeded | 0.967-0.982 | 0.965 [0.949, 0.975], below | 0.984 [0.971, 0.993], above |
| bgal (2D) | validation: particles given a shortlist (%) | 71.4-80.5 | 75.1 [74.4, 76.1], inside | 80.2 [74.1, 85.4], inside |

M moves one read of (b) above A's range (`T1000` seeded), and R moves four. **The experiment continues to Phase 3.**
None of these means yet reaches decision 7 (b)'s bars (for example `T1000` seeded 0.330, bgal seeded 0.587, bgal within
6° 0.894; the class student's held-out reads 1.011). The verdict is Phase 3's, over every set.

**The other reads**, for information (the verdict of decision 7 (a) is taken in Phase 3 over every set).
- **M** is below A's range on:
  - `T1000` and bgal: the truth or teacher in the cone, and seeded rows in the cone;
  - bgal 3D: within 6° and within 12°;
  - 2D: held-out classes recovered and held-out seeded;
  - the validation: shortlists holding the full search's class (73.3% against 73.4-79.0%).
- **R**:
  - In 3D it is above A's range on bgal: mode = truth, within 12°, the agreeing in-plane share, and below in median
    distance, cone size and in-plane error (all better).
  - It is below on: the cone holding the truth on `T1000` and bgal, seeded rows in the cone, bgal held-out recovered
    (0.994 against 1.000), the 2D held-out recovered (0.995 against 1.000), the 2D prediction's stored class among the
    candidates (0.978 against 0.988-0.992), and the validation's shortlists holding the full search's class (72.0%
    against 73.4-79.0%).
- **The fixed-reference in-plane error** (decision 6 as revised before these runs; mean [range] in degrees):

  | read | A | M | R |
  | --- | --- | --- | --- |
  | bgal 3D, training rows against the stored pose | 0.876 [0.866, 0.883] | 0.892 [0.880, 0.907], worse | 0.835 [0.818, 0.859], better |
  | bgal 2D, held-out rows against the stored class's angle | 0.900 [0.839, 0.934] | 0.912 [0.887, 0.934], inside | 0.845 [0.809, 0.887], inside |

- **Costs** (reported, not judged): the teacher's fill is 0.1 to 0.2 s longer on these small sets. Training on bgal
  takes 26 to 28 s against 24 s.
- The 2D held-out reads (classes recovered, seeded, recovered) are taken against each variant's own teacher, the
  marginal's mode class in M and R. The prediction's and the validation's reads are against fixed references (the
  stored class, the full search).

**Exit.** Decision 8 applied: continue. Files changed this phase: this plan (decision 6's revision of the in-plane
reads, made before any variant run; this entry). Evidence: `R/scratch/logs/chain_p2.log`, `R/scratch/keep/
phase2_compare.txt`, `R/scratch/keep/fixed_inpl.tsv`, the kept runs in `R/scratch/keep/runs/{M,R}_p2/`,
`runs2d/{M,R}_p2/`, `validate/{M,R}_p2/`.

### Phase 3, the matrix and the verdict (2026-11-15, the Dell, done: both variants rejected for both students)

**What ran** (`R/scratch/logs/chain_p3.log`, `chain_p3b.log`, `chain_p3c.log`; no run failed):
- M and R at seeds 1 to 3 on the remaining sets, with Phase 0's commands and one job at a time at 24 threads (the
  one-minute load 15 to 24, the jobs' own threads). In 3D, the harness on `T200`, `T5000` and `T`, where `T` also
  gives the decentring curve. In 2D, `learn_class_student` on the frozen trpv1 teacher and `validate_class_student`
  on trpv1.
- Reference matching beside the seed-1 validations of both sets. On bgal it needed a fresh seed-1 validation of
  each variant (phase label `3m`), because Phase 2's had been trimmed. The validation's own solve2D teacher at 24
  threads is not the same twice (the seed-1 teachers of A, M and R differ in their checksums), so the matching
  has to sit beside its own validation. The `3m` validations serve the matching only; the verdict keeps Phase 2's.
- `validate_pose_student` on bgal at seed 1 for A and for R, the better variant, for inspection.

**The verdict, decision 7 as written** (`R/scratch/keep/tools/verdict7.py`, output `R/scratch/keep/verdict7.txt`).
- (a): no held-out read's three-seed mean may fall below A's range on any of the student's sets; for a read where
  lower is better, above A's range.
- The held-out reads are decision 6's. 3D: within 6°, the truth (or the teacher on bgal) in the cone, the share
  seeded, held-out recovery, and the in-plane error (the truth sets against the truth; bgal against the stored pose,
  decision 6 as revised). 2D: held-out classes recovered, held-out seeded, the in-plane error against the stored
  angle, and the validation's report (particles given a shortlist, shortlists holding the full search's class,
  centred and perfectly centred, shortlist length, the class averages' FRC 0.5 against the teacher's). FRC is the
  Fourier ring correlation; 0.5 is the threshold its resolution is read at.
- Taken as written but not judged: the confidence threshold. It is a calibration value with no better or worse
  direction; it is reported in the tables.
- (b): decision 7's bars, A's max + 2 × A's range, on at least three of the five 3D sets or on both 2D sets.

| student | variant | (a): reads worse than A's range | (b): sets past the bars | verdict |
| --- | --- | --- | --- | --- |
| pose (3D) | M | 5: truth in the cone on `T1000` (0.971 against 0.974-0.986) and bgal (0.952 against 0.957-0.961); held-out recovered on `T` (0.997 against 0.998-1.000); bgal within 6° (0.864 against 0.872-0.879); bgal in-plane error (0.892° against 0.866-0.883°) | 0 of 5 | rejected |
| pose (3D) | R | 3: truth in the cone on `T1000` (0.968) and bgal (0.952); bgal held-out recovered (0.994 against 1.000) | 0 of 5 | rejected |
| class (2D) | M | 3: bgal held-out classes recovered and held-out seeded (0.965 against 0.967-0.982); bgal shortlists holding the full search's class (73.3% against 73.4-79.0%) | 0 of 2 | rejected |
| class (2D) | R | 4: bgal shortlists holding the full search's class (72.0%); trpv1 held-out seeded (0.856 against 0.856-0.879, by 0.0003); trpv1 in-plane error (1.201° against 1.082-1.198°); trpv1 shortlist length (7.7 against 7) | 0 of 2 | rejected |

Decision 7 (c), R over M, has nothing to decide: neither variant is adopted.

**The reads of decision 7 (b)** against the bars (three-seed means):

| set | read | A | bar | M | R |
| --- | --- | --- | --- | --- | --- |
| `T200` | seeded / within 6° | 0.202 / 0.743 | 0.669 / 0.996 | 0.221 / 0.732 | 0.223 / 0.735 |
| `T1000` | seeded / within 6° | 0.245 / 0.797 | 0.330 / 0.956 | 0.273 / 0.774 | 0.303 / 0.778 |
| `T5000` | seeded / within 6° | 0.434 / 0.870 | 0.708 / 1.003 | 0.463 / 0.869 | 0.477 / 0.875 |
| `T` | seeded / within 6° | 0.529 / 0.930 | 0.818 / 1.012 | 0.576 / 0.928 | 0.561 / 0.932 |
| bgal 3D | seeded / within 6° | 0.426 / 0.876 | 0.587 / 0.894 | 0.403 / 0.864 | 0.492 / 0.888 |
| bgal 2D | classes recovered / held-out seeded / shortlisted | 0.975 / 0.975 / 75.6% | 1.011 / 1.011 / 98.7% | 0.965 / 0.965 / 75.1% | 0.979 / 0.984 / 80.2% |
| trpv1 2D | classes recovered / held-out seeded / shortlisted | 0.866 / 0.869 / 70.1% | 0.928 / 0.925 / 83.6% | 0.869 / 0.870 / 67.9% | 0.855 / 0.856 / 68.0% |

**What the reads say.**
- The marginal as teacher with its readout (M) changes nothing measurable. M's means sit inside A's ranges except a
  few reads at their edges, about half better and half worse; the share seeded rises on every truth set but stays
  inside the range on three of the four.
- The rotation loss (R) helps the pose student on bgal, the real set with the most direction-to-rotation coupling:
  - seeded 0.492 against 0.426, within 6° 0.888 against 0.876, the in-plane error 0.835° against 0.876°;
  - mode equal to the teacher's 0.180 against 0.172, and the median cone 18 against 20 directions.
- R costs on that same set: the cone holds the teacher less often (0.952 against 0.957-0.961), and about one seeded
  held-out row in 180 loses its recovery (0.994 against 1.000).
- On the truth sets R is within the noise; on `T1000` its seeded share is above A's range (0.303) but its cone holds
  the truth less often.
- In 2D, R seeds more held-out rows on bgal, but its shortlists hold the full search's class less often on bgal.
  On trpv1 it is level or slightly worse.

**The decentring curve of the `T` students** (the truth in the cone / the mode within 12°,
mean over the seeds, on the test half decentred by r px):

| r (px) | 0 | 1 | 2 | 3 | 5 | 8 |
| --- | --- | --- | --- | --- | --- | --- |
| A | 0.995 / 0.993 | 0.991 / 0.981 | 0.955 / 0.848 | 0.710 / 0.394 | 0.216 / 0.048 | 0.153 / 0.030 |
| M | 0.993 / 0.993 | 0.988 / 0.981 | 0.951 / 0.847 | 0.715 / 0.395 | 0.217 / 0.048 | 0.152 / 0.030 |
| R | 0.995 / 0.993 | 0.990 / 0.980 | 0.947 / 0.839 | 0.689 / 0.376 | 0.214 / 0.045 | 0.154 / 0.029 |

**Costs** (reported, not judged; training wall time at 24 threads, the mean over the seeds in seconds):

| set | A | M | R |
| --- | --- | --- | --- |
| `T200` | 8.2 | 9.6 | 10.2 |
| `T1000` | 18.4 | 20.8 | 24.2 |
| `T5000` | 86.4 | 95.8 | 106.0 |
| `T` | 368.8 | 800.5 | 411.1 |
| bgal 3D | 23.7 | 26.3 | 27.6 |
| bgal 2D | 11.3 | 11.2 | 14.8 |
| trpv1 2D | 58.9 | 64.9 | 59.2 |

- M's mean on `T` holds one seed that tried all three widths (1,542 s, keeping 32 channels against 16 at the other
  seeds; the rule tries the widest design from 10,000 rows when the narrower ones gain).
- The teacher's fill on `T`: 62 s against 56 s.
- Peak memory on `T`: M 1.8 to 2.4 GB, R 2.5 GB, A 1.76 GB.
- Prediction is unchanged.

**Reference matching** beside the seed-1 validations (reported; `R/scratch/keep/refmatch/*_s1/reads.txt`). Each line
gives the share of unseen particles in the full search's class, by matching and by the student, then the class FRC
0.5 of matching / full search / student against the teacher's:
- bgal:
  - A: 87.0% and 33.9%, 14.84 / 14.84 / 20.40 Å;
  - M: 86.5% and 34.2%, 19.80 / 15.19 / 15.19 Å;
  - R: 85.3% and 34.1%, 20.40 / 20.40 / 18.13 Å.
- trpv1:
  - A: 87.4% and 36.1%, 22.23 / 18.31 / 22.23 Å;
  - M: 87.3% and 35.9%, 22.23 / 19.45 / 22.23 Å;
  - R: 86.6% and 35.6%, 23.94 / 17.29 / 19.45 Å.
- The matching is variant-independent; its spread here is that of the validations' own teachers.

**For Hans to inspect** (not gated).
- The class averages of the student's assignments: seed 1, under `R/scratch/keep/validate/<build>_p<phase>/
  <set>_s1/kept/` (`student_cavgs.mrc`, beside the teacher's and the full search's).
- The unseen map from `validate_pose_student` on bgal, seed 1, under `R/scratch/keep/validate_pose/{base,R}_s1/kept/`
  (the report and the maps). Command:
  `SIMPLE_SEED=1 simple_exec prg=validate_pose_student projfile=<copy of /data2/bgal_addon/3_selection/bgal_addon.simple>
  nsample=1500 mskdiam=180 pgrp=d2 pgrp_start=c1 nthr=24`.
- The unseen particles at the student's poses and registered shifts reach FSC 0.143 at 4.87 Å against the teacher's
  map for A and 5.10 Å for R. At the 2D shifts they reach 5.44 Å and 6.04 Å. Between the halves: 7.59 against 7.77 Å
  registered, 8.82 against 7.96 Å at the 2D shifts.
- R seeds 27.4% of all particles against A's 20.3%, and 45.5% of the training particles against 34.6%, at the same
  calibration (95% of seeds recovered). The readout's in-plane angle against the teacher's: median 1.00° for both.
- The two runs used the same teacher map (the control map reaches 2.55 Å in both).

**Exit.** Decision 7 applied per student and per variant: both variants are rejected for both students. Phase 4
writes the rejected record and changes no source file. Files changed this phase: this plan (this entry).
Evidence: the logs named above, `R/scratch/keep/verdict7.txt`, `summary.tsv`, `summary2d.tsv`,
`summary_validate.tsv`, `fixed_inpl.tsv`, `decentring.tsv`, the kept runs under `R/scratch/keep/runs/{M,R}_p3/`,
`runs2d/{M,R}_p3/`, `validate/{M,R}_p3/`.

### Phase 4, the rejected record and the report (2026-11-15, the Dell, done)

Decision 7 rejected both variants for both students in Phase 3, so this phase took the plan's rejected path:
- no source file changed;
- the record is `student_inplane_posterior_rejected.md` beside this plan (the maintainer moves it to `rejected/`);
- the report is `student_inplane_posterior_report.md`, also beside this plan, which stays in `planned/`;
- this plan's status line carries the verdict and the date.

The patches the record names were checked on a fresh `git archive 86babf86b`: the 16 source files of
`R/scratch/keep/variants/M.patch` and `R.patch` apply cleanly, and R's replay equals the source that was built and
tested. Each patch also holds a hunk for `doc/code_overview/code_base_map.md`. That file is git-ignored and
regenerated by the build (`scripts/CMakeLists.txt`), so the hunk is not needed. The file-table row Phase 1 gave it
was unnecessary, and the report lists it among the deviations.

No source changed, so the fast gate was not re-run, as this phase's rule allows. Exit: the rejected note and the
report written. Files changed this phase:
- this plan (status line and this entry);
- the new `planned/student_inplane_posterior_rejected.md` and `planned/student_inplane_posterior_report.md`.

### Phase 4 reopened under rulings R4 and R5 (2026-11-15, the Dell, stopped for the maintainer)

R4 replaced decision 7 with "accept unless worse" for R. It asked first for an audit of the teacher's in-plane
weighting against the refinement's, then for "worse" read over eight paired seeds. R5 moved the base to `4282d58a3`
part-way through. This entry supersedes the rejected path above; the rejected record stays until the maintainer
rules.

**The audit found a difference, which was fixed (decision 2's revision).**
- The learn programs run at `prob_beta` = 0.3 (every rows header says so). The teacher's in-plane marginal and
  rotation targets were summed at 0.3 T.
- The in-plane posterior insertion weighs every row at the physical temperature T. In 3D it also adds the shift
  prior's term to each view weight (`inpl_prior_penalty`).
- The old teacher also summed in double, where `inpl_quadrature` rounds each term to single precision.
- Planted back into the corrected code, the old temperature moves a reference's lnz by up to 1.94 on bgal, a factor
  of about 7 in its mass.
- Fixed in M and R (`R/scratch/keep/variants*/`; the uncorrected patches are kept in `variants/before_audit/`):
  - the teacher calls `inpl_quadrature` itself (new optional output `post`), at the insertion's temperature
    and shift, with `inpl_prior_penalty`'s term in 3D;
  - the training loss is d_min − (lnz − pen)/T, read at `prob_beta` × T as A's rows are (decision 2);
  - the rows file carries pen and T.

**The audit itself is now part of the learn programs.** On the first batch of real rows, every entry is
recomputed by the insertion's own route and compared with the teacher's:
- the route: the row at the shift the insertion would take the entry at (`record_seed_shift`'s rule), the
  particle shifted once above `SHERRSQ`, the insertion's temperature, `inpl_quadrature`, `inpl_prior_penalty` and
  `inpl_view_weights`;
- the targets' mapping is read on every row, from the teacher's grid angles to the student's rotation and back
  through `readout_e3`, the mirror's through its pair's frame;
- any difference stops the program (`<fname>_teacher_audit.txt`).

On bgal, seed 1, the corrected R gave the same audit on both bases:

| read | bgal 3D (teacher map) | bgal 2D (frozen teacher) |
| --- | --- | --- |
| entries compared (first batch × references) | 750,000 | 45,950 |
| in-plane rows compared | 6,000 | 7,352 |
| largest difference: shift, temperature, lnz, d_min, shift-prior term, in-plane posterior | all exactly 0 | all exactly 0 |
| view weights against the training targets, in units of their single-precision bound | 0.103 | 0.092 |
| mapping: references read (on the mirror's side) | 12,000 (6,215) | 22,048 (0) |
| mapping: largest error against half a student step | 0.5845° < 0.6429° | 0.5844° < 0.6429° |
| the teacher grid against the uniform grid the targets assume | 3.1e-5° | 3.1e-5° |

Every pair on the bgal d2 grid has frame angle 0, so a wrong frame sign cannot show on these rows; the tester pins
frames of 75° and 180°. Nine checks were planted first and every one was caught, on both bases
(`R/scratch/keep/plants/plants_p4*.log`):
- the old temperature: the unit suite, and both learn audits;
- the shift prior left out of the ranking;
- the audit blind to lnz, and blind to the view weights;
- an unnormalised posterior;
- the 3D shift-prior term zeroed: the learn audit, pen off by 9.53;
- the 2D refined classes left at the seed: the learn audit, shift off by 7.07 px;
- a target one step off: both learn audits, 1.87°;
- the mirror's reflection ignored: the 3D learn audit, 180°.

**The new base (R5).**
- `4282d58a3` exported to `R/scratch/keep/src_base_4282d58a3`, built Release and Debug (`build_base_n`,
  `build_base_n_debug`); the fast gate passed 15 of 15 on both.
- The corrected patches were carried onto the new base (`variants_n/{M,R}.patch`, 17 files). One hunk was rejected:
  the tester's import list, where the base gained `report_prob_readout`. It was resolved by hand, dropping a
  duplicate import the old patch carried.
- Both variants passed the fast gate (15 of 15), `unit_numerics` (1,905 assertions) and `unit_pftc_align2D3D`
  (6,366).
- Off-state identity, `simulated_workflow` at `SIMPLE_SEED=7`:

  | comparison | 1jxy files | 6vxx files | exclusions |
  | --- | --- | --- | --- |
  | `base2n` against `base1n` | 223 of 223 | 186 of 186 | none |
  | M against `base1n` | 223 of 223 | 186 of 186 | none |
  | R against `base1n` | 223 of 223 | 186 of 186 | none |

- The base's 1jxy shared-memory solve3D fails its pose gate at this seed (19.75 against a floor of 15), as R5 said.

**Seeds 1 to 3 on every set, new base.**
- Both arms ran at 16 threads (A as `An`, R as `Rn`, phase label 5). The maintainer's seed profile held about six
  cores beside the run, and a 24-thread team stalled beside it (T200 took 13 times as long).
- A's reads equal the old base's, seed for seed.
- R falls below A's range on two reads only (`R/scratch/keep/r5_screen_s123.txt`):
  - the truth in the cone on `T1000`: 0.9718 against 0.9738 to 0.9864;
  - bgal 2D, shortlists holding the full search's class when perfectly centred: 85.7% against 86.0 to 90.1%.
- Both lie on R4's sets, so the paired list stayed R4's: `T1000`, bgal 3D, and bgal and trpv1 in 2D, learn and
  validate.

**The paired test, seeds 1 to 8** (`R/scratch/keep/r5_paired.txt`; R − A at the same seed, its standard error;
worse = beyond twice the standard error in the worse direction; † a read R4 named):

| set | read | A | R | R − A | SE | verdict |
| --- | --- | --- | --- | --- | --- | --- |
| `T1000` | within 6° | 0.7871 | 0.7789 | −0.0082 | 0.0037 | **worse** |
| `T1000` | truth in the cone † | 0.9768 | 0.9723 | −0.0045 | 0.0032 | within 2 SE |
| `T1000` | seeded | 0.2555 | 0.2787 | +0.0232 | 0.0227 | within 2 SE |
| `T1000` | held-out recovered | 1.0000 | 1.0000 | 0 | 0 | level |
| `T1000` | in-plane error (°) | 0.9366 | 0.9319 | −0.0048 | 0.0032 | within 2 SE |
| bgal 3D | within 6° | 0.8641 | 0.8716 | +0.0075 | 0.0084 | within 2 SE |
| bgal 3D | teacher in the cone † | 0.9481 | 0.9541 | +0.0060 | 0.0035 | within 2 SE |
| bgal 3D | seeded | 0.3611 | 0.4114 | +0.0503 | 0.0404 | within 2 SE |
| bgal 3D | held-out recovered † | 1.0000 | 1.0000 | 0 | 0 | level |
| bgal 3D | in-plane error against the stored poses (°) | 0.8969 | 0.8464 | −0.0505 | 0.0076 | better |
| bgal 2D | classes recovered | 0.9737 | 0.9846 | +0.0109 | 0.0033 | better |
| bgal 2D | held-out seeded | 0.9736 | 0.9877 | +0.0141 | 0.0035 | better |
| bgal 2D | in-plane error, held-out against the stored poses (°) | 0.9006 | 0.8410 | −0.0596 | 0.0225 | better |
| bgal 2D | particles given a shortlist (%) | 74.66 | 77.98 | +3.31 | 2.42 | within 2 SE |
| bgal 2D | shortlist holds the full search's class (%) † | 74.84 | 71.64 | −3.20 | 1.09 | **worse** |
| bgal 2D | the same, perfectly centred (%) | 87.56 | 86.85 | −0.71 | 0.98 | within 2 SE |
| bgal 2D | shortlist length (median) | 4.875 | 4.250 | −0.625 | 0.324 | within 2 SE |
| bgal 2D | class FRC 0.5 (Å) | 17.01 | 17.42 | +0.42 | 0.98 | within 2 SE |
| trpv1 | classes recovered | 0.8600 | 0.8812 | +0.0213 | 0.0058 | better |
| trpv1 | held-out seeded † | 0.8615 | 0.8845 | +0.0230 | 0.0057 | better |
| trpv1 | in-plane error, held-out against the stored poses (°) † | 1.1532 | 1.1615 | +0.0083 | 0.0250 | within 2 SE |
| trpv1 | particles given a shortlist (%) | 69.39 | 71.90 | +2.51 | 0.65 | better |
| trpv1 | shortlist holds the full search's class (%) | 91.08 | 90.01 | −1.06 | 0.36 | **worse** |
| trpv1 | the same, perfectly centred (%) | 96.90 | 96.61 | −0.29 | 0.23 | within 2 SE |
| trpv1 | shortlist length (median) † | 7.125 | 7.125 | 0 | 0 | level |
| trpv1 | class FRC 0.5 (Å) | 20.82 | 21.10 | +0.28 | 0.85 | within 2 SE |

**R is worse on three reads.** One of them is a read R4 named (bgal 2D, shortlists holding the full search's class).

**The look for an implementation error behind them found none:**
- The teacher equals the insertion's weighting exactly on real rows, on both bases. The plants show the audit
  catches each way it could differ.
- The adopted form of R (the switches removed, `R/scratch/keep/adopt/`) trains the same models as the tested R:
  learn on bgal 3D and 2D at seed 1 gives identical headers, apart from the project path.
- The shortlist deficit is not a worse ranking (`R/scratch/keep/topk2d.tsv`, `tools/topk2d.sh`). On all unseen
  particles, R ranks the full search's class as A does:
  - bgal: top 1 +0.004 (SE 0.009), top 6 −0.008 (SE 0.009);
  - trpv1: top 1 to 4 level, top 5 and 6 −0.006 (SE 0.003).
- R shortlists more particles (bgal +91, trpv1 +643 per run), the extra ones harder, and recovers its own teacher's
  class more often (classes recovered +0.011 and +0.021). Its teacher is now the insertion's marginal. The full
  search's class is one draw from the table's single-rotation losses, a different distribution that R no longer
  learns.
- The 3D mirror mapping is symmetric: R lowers the bgal in-plane error alike on both members of a pair (median
  0.649° unmirrored, 0.617° mirrored, against A's 0.786° and 0.747°; `inpl_split.tsv`).
- `T1000`'s within 6° deficit is small: R's median distance to the truth is −0.06° to +0.17° from A's, seed by
  seed. R's cones are narrower on six of the eight seeds, and it seeds more on seven.

So the phase stops with the paired table for the maintainer, as R4 item 4 says. The adoption copy is ready if the
maintainer accepts R anyway: `R/scratch/keep/adopt/src`, built and tested.

**A finding on the new base, not the variant's.** One validation of R (bgal, seed 5) crashed after 4.5 s with
SIGSEGV. The crash was in `sp_project`'s segment writer inside the solve2D teacher's first refine2D iteration;
no student code runs there.
- Run again with the same command, it completed; it was re-run through the tool for the table.
- The log is `R/scratch/keep/crash_Rn_p5_bgal_s5_validate.log`.

**Deviations.**
- Both arms ran at 16 threads, beside the maintainer's profile.
- An old-base run of R at 24 threads (`R4_p4`, seed 1 and part of seed 2) was stopped by R5 and is kept as
  reference only.

No repository source changed, so no repository build or gate was needed. Files changed this phase:
- this plan: decision 2's revision, R4 and R5 in decision 7, the status line and this entry;
- the report, rewritten.

### Phase 4 under ruling R6: the class student's physical teacher (2026-11-15, the Dell, stopped for the maintainer)

R6 moved the base to `603c44490` and narrowed R4 to the class student:
- variant **P**: R's teacher and network, with the class distribution read at the physical temperature as well, so
  the learn program uses no `prob_beta`;
- `validate_class_student` judged against the physical posterior's best class;
- the pose student keeps A, and the 3D marginal is not adopted.

R stays as recorded above; it was not re-run.

**The sources** (`R/scratch/keep/variants_r6/`, built Release on `git archive 603c44490`):
- **Carrying R's patch:** two import lines conflicted with the new base, and were resolved by hand.
- **The 3D files restored to the base:** `simple_eul_prob_tab.f90`, `simple_strategy_pose_student.f90`, and
  `simple_pose_student.f90` with its tester.
- **The core's readout** is chosen by the owner in code: an optional last argument of `polar_student%new`,
  recorded as `readout=lse` only for the log-sum-exp. The weights file stays at version 2, its format unchanged. The
  pose student's code and model files are therefore the base's. A new test pins that the maximum writes no readout
  key.
- **One constant selects the class student's path,** `CLASS_STUDENT_PHYSICAL`:
  - false gives the A arm, the base's class student with the new validation;
  - true gives P, where `learn_class_student` sets `prob_beta` to 1 (`CLASS_STUDENT_PROB_BETA`).
- **The reference** is `run_class_posterior_mode`, the learn routine's table in a reference mode on the centred
  unseen particles at the full box:
  - each particle's class of largest in-plane marginal mass at T is written as refine2D writes a pose: the stored
    shift plus the refined increment, and e3 = 360 − rot of the class's in-plane mode;
  - the class averages of that assignment come from `make_cavgs`;
  - it replaces the one refine2D pass at `prob_beta` 100, whose constant `CLASS_DENSE_PROB_BETA` was dropped;
  - its first batch is audited, and every written class is checked against its row's marginals.

**The audit, extended to the class distribution (R6 item 2).** P's learn on bgal 2D at seed 1:
- the in-plane terms equal the insertion's exactly (45,950 entries, 7,352 rows);
- all 919 rows of the first batch are read at the physical temperature;
- the class distribution the student trains on equals exp(lnz − T d_min) normalised over the classes to 9.1e-6,
  0.047 of the bound its single-precision losses allow;
- the targets' mapping lies within half a rotation step.

**Planted violations (`R/scratch/keep/plants/plants_r6.log`), all caught:**
- the learn program's `prob_beta` left out: the class-mass audit, "not at the physical temperature";
- the class distribution read at 1.1 T, and the rows left at the drawn losses: the learn audits;
- the reference's class moved off the mode: the reference's own check;
- the readout key written for every student: the new tester check. It was not caught at first, and the test was
  added for it.

**The reference's check stopped one validation** (A arm, trpv1, seed 6). It had compared the class written, the
argmin of the row's single-precision losses, with the exact double-precision argmax of lnz − T d_min. Either a
near-tie reversed by rounding, or a particle without a temperature (read at its argmin), can separate the two.
- The check now applies the row's own rules: a gap within the losses' rounding counts as a near-tie and is
  reported; a particle without a temperature is checked by the argmin.
- Re-run, the validation passed. Its solve2D teacher differs from the first run's (it is not reproducible at 24
  threads), so which cause it was is not known.

**Gates and identity on `603c44490`:**
- fast gate 15 of 15 on the base (Release and Debug) and both arms;
- `unit_numerics` (1,899 assertions) and `unit_pftc_align2D3D` (6,384) pass on both arms;
- off-state identity, `simulated_workflow` at `SIMPLE_SEED=7`: `base2r6`, A and P each equal `base1r6` (1jxy 223 of
  223 files, 6vxx 186 of 186, no exclusions).

**Guards** (`R/scratch/keep/guard_r6/`, learn on bgal at seed 1, 24 threads):
- the pose student of both arms is the base's, header and weights byte for byte;
- the A arm's class student is the base's except its model version (2 against 1); its weights are identical.

**The paired test, seeds 1 to 8, 24 threads** (`R/scratch/keep/r6_paired.txt`; P − A at the same seed and its
standard error; the validations against the physical posterior's best class):

| set | read | A | P | P − A | SE | verdict |
| --- | --- | --- | --- | --- | --- | --- |
| bgal | held-out classes recovered | 0.9737 | 0.9814 | +0.0077 | 0.0039 | within 2 SE |
| bgal | held-out seeded | 0.9736 | 0.9836 | +0.0100 | 0.0041 | better |
| bgal | in-plane error, held-out against the stored poses (°) | 0.9006 | 0.8423 | −0.0584 | 0.0205 | better |
| bgal | particles given a shortlist (%) | 70.64 | 80.69 | +10.05 | 3.17 | better |
| bgal | shortlist holds the reference's class (%) | 97.71 | 94.59 | −3.13 | 0.45 | **worse** |
| bgal | the same, perfectly centred (%) | 96.38 | 94.88 | −1.50 | 0.58 | **worse** |
| bgal | shortlist length (median) | 5.25 | 4.38 | −0.88 | 0.30 | better |
| bgal | class FRC 0.5 (Å) | 17.58 | 15.90 | −1.68 | 0.90 | within 2 SE |
| trpv1 | held-out classes recovered | 0.8600 | 0.9001 | +0.0401 | 0.0027 | better |
| trpv1 | held-out seeded | 0.8615 | 0.9055 | +0.0440 | 0.0028 | better |
| trpv1 | in-plane error, held-out against the stored poses (°) | 1.1532 | 1.1732 | +0.0200 | 0.0253 | within 2 SE |
| trpv1 | particles given a shortlist (%) | 70.03 | 75.59 | +5.56 | 1.34 | better |
| trpv1 | shortlist holds the reference's class (%) | 99.09 | 98.14 | −0.95 | 0.11 | **worse** |
| trpv1 | the same, perfectly centred (%) | 99.68 | 99.45 | −0.23 | 0.07 | **worse** |
| trpv1 | shortlist length (median) | 7.125 | 6.50 | −0.63 | 0.26 | better |
| trpv1 | class FRC 0.5 (Å) | 20.87 | 21.73 | +0.86 | 0.47 | within 2 SE |

**P is worse on four reads,** the shortlists holding the reference's class, centred and perfectly centred, on both
sets. It is better on nine. Against the new reference both arms' shortlists hold the class far more often than
against the old one (bgal A 97.7% against 74.8%).

**The look for an implementation error found none:**
- **The teacher:** the audits hold exactly on real rows, and the class distribution is the physical posterior's to
  its rounding.
- **The prediction** reads the readout it was trained with (the header check enforces it). P recovers its teacher's
  class on held-out rows more often than A (trpv1 +0.044).
- **On all unseen particles** (`R/scratch/keep/topk2d_r6.tsv`), P ranks the reference's class as A does on bgal:
  top 1 +0.004, top 6 −0.003, within noise. On trpv1 it ranks it 0.6 to 1.0 points lower at top 3 to 6, beyond 2 SE.
- **The shortlists:** P shortlists many more particles (bgal +277 per run, trpv1 +1,423) with shorter lists. The
  added particles, centred by `center_ptcls` rather than by the teacher, hold the reference's class less often.
- **Reading:** P's sharper posterior is calibrated on the teacher's centred rows and is more confident on unseen
  ones than their centring supports. This is the trade R showed, measured now against the posterior P learns.

So the phase stops with the paired table, as R6 item 4 says.
- **The adoption copy is ready** if the maintainer accepts P anyway: `R/scratch/keep/adopt_r6/src`, P with its switch
  removed.
- **The class page's text** for it is drafted in `R/scratch/keep/tools/docs_adopt_r6.py`.
- **Costs:** the validations take as long in both arms (bgal 71 s, trpv1 320 against 310 s); the reference pass is
  cheaper than the refine2D pass it replaced.

No repository source changed. Files changed:
- this plan (R6 in decision 7, a file-table row, the status line, this entry);
- the report.

### Phase 4 under ruling R7: P adopted for the class student (2026-11-15, the Dell, done)

R7 adopted P for the class student. It accepted R6's four flagged reads as the calibration's trade: P holds the 95%
target where A is over-cautious, and shortlists more particles. The adoption is the class student's path only, on
`603c44490`, from the prepared copy (`R/scratch/keep/adopt_r6/src`, P with its switch removed).

**What went into the repository:**
- 13 source files changed and one new, `simple_eul_prob_marginal.f90` (the report lists them);
- `doc/algorithms/class_student.md`: the teacher, the readout, the rotation loss, the audit, the validation's
  reference, the reads.

The class model file is at version 2, its readout recorded. The pose student, its strategy and the 3D table are
unchanged: `git diff` is empty for `simple_pose_student.f90` and its tester, `simple_strategy_pose_student.f90` and
`simple_eul_prob_tab.f90`.

**One test file was added to the table and changed.** The `class_student` workflow test failed its check
`the_report_holds_every_read` on the first Release run: the test reads the report's timing line by its label, and
R7's rewording had renamed it. The test's label became "Reference: the physical posterior's best class"
(`simple_commanders_test_highlevel_class_student.f90`, a row added under R3 and R7). The chain was stopped by PID and
re-run on the fixed sources, so both trees carry the same code.

**The gates** (`R/scratch/keep/r7/`, the chain `tools/chain_r7.sh`), on the Release tree `repo/build` and the Debug
tree `repo/build_debug`:

| gate | Release | Debug |
| --- | --- | --- |
| fast gate | 15 of 15 | 15 of 15 |
| `unit_numerics` | 1,900 assertions pass | 1,900 pass |
| `unit_pftc_align2D3D` | 6,384 pass | 6,384 pass |
| `pose_student` workflow test (CTest, 8 threads) | passed, 37 s | passed, 121 s |
| `class_student` workflow test (CTest, 8 threads) | passed, 23 s | passed, 91 s |

- **`validate_class_student` on bgal at seed 1** (Release, 24 threads, 63 s):
  - shortlists for 79.4% of the unseen particles;
  - they hold the physical posterior's best class for 95.9%, and 96.1% when perfectly centred, both meeting the 95%
    goal;
  - the teacher's audits hold, and the reference checks all 2,757 classes it writes, with no near-tie.
- **Off-state identity, the Release tree against `base1r6`** (`simulated_workflow`, `SIMPLE_SEED=7`): 1jxy 223 of 223
  files, 6vxx 186 of 186, no exclusions.
- **The pose student is the base's.** `learn_pose_student` on bgal at seed 1 from the repository tree gives the
  base's model header and weights byte for byte (md5 `8f0d176e…`; `R/scratch/keep/r7/pose_check/`).

The rejected record (`student_inplane_posterior_rejected.md`) was deleted, as R7 directs. The report was rewritten
as the adoption, with the paired tables of R4 and R6 and the base each came from. This plan's status line says
what was delivered. Exit: the gates of section 5 for the class student pass on both trees,
and the `pose_student` guard is unchanged.
