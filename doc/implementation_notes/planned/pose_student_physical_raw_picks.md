# The pose student revisited: the physical teacher, raw picks and alignment over the cone

Status, 2026-11-18: **done; the centring policy adopted** (ruling R7; report: `pose_student_physical_raw_picks_report.md`).
Phase 2 kept the alignment over the cone (SK) and removed the mode's pair (S1); Phase 3 adopted the physical teacher P3
on the maps (ruling R5); Phase 4's raw-pick route was not adopted (its map coarser than the 2D-centred route's on
streptavidin), and under ruling R7 the policy (training and predicted particles centred by the prior solve2D, shift
refinement only after the prediction, over the cone) was measured as C3SK and adopted. Phase 4 ran on master's
`68638be28` (ruling R6), where the gates, the identity and the
validations ran. Measurements of Phases 0 to 3 are on `fd9af0add` (ruling R4,
2026-11-17), the commit that removes the 3D in-plane posterior insertion and `simulated_workflow_1jxy`; the run
launched on `1ea0f8340` (the 3D in-plane insertion the cone, `inpl_posterior` removed), whose Phases 0 and 1 are
recorded below as superseded. Depends on: the pose student (`doc/algorithms/pose_student.md`), what the
class student established (`planned/student_inplane_posterior.md` with its report: the physical teacher P;
`planned/class_student_raw_picks.md` with its report: no centring, alignment after prediction over the candidate
list), and the current solve3D, whose solutions are the best yet (Hans, 2026-11-16).

## 1. Why

The pose student has not changed since its training-speed and working-grid runs except its prediction seed
(`euler_inpl_closest`, 2026-11-15). Meanwhile the class student showed four things that carry over:

- **The teacher should be what the run ends at, the physical posterior.** The class student's tempered teacher
  (`prob_beta` x T) made it over-cautious (97.7% where its threshold was calibrated for 95%); the physical teacher
  held the target and shortlisted 5 to 10 more particles in 100. 3D refinement has no in-plane insertion (removed in
  `fd9af0add`, ruling R4: each pose is inserted at its own angle), so the pose student's matching teacher is the
  one-rotation table rows read at the physical temperature T, not the in-plane marginal (R, rejected for 3D on
  2026-11-15).
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

- Off state: refine2D and refine3D byte-identical at one thread with a fixed seed (`simulated_workflow` suite 6vxx,
  `test=identity`; 1jxy no longer exists on the base, ruling R4); the class student byte-identical to the base (model
  header and weights).
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
   How the run computes it (Phase 0, for the maintainer's review): learn_pose_student's rows (the dense table at the
   student's grid and band, prob_beta 1, the class reference's disk of three 10 Å widths, every unseen particle not
   yet centred) against the teacher's map, started from the most centred position the particle has (the 2D solution's
   shift on the real sets, the stored shift of the centred test half on the truth sets); the row's mode with its
   in-plane angle and refined shift. One reference per particle and teacher, frozen, serves every arm: the arms see
   the same images. Started from the raw picks instead, the same table put its mode within 6° of the truth for 72.5%
   of the truth sets' test half against 99.1% from the centred start (it agreed with the centred reference for 64%):
   the dense search's screen misses off-centre modes, so a raw-pick start would judge every raw-pick arm against a
   broken reference. The raw-pick dense pass is kept as a read of its own (D0, decision 7).
3. (2026-11-16, Hans) **No pre-centring of unseen particles**: they enter at their picked positions (shifts zero).
   **Withdrawn by ruling R7 (Hans, 2026-11-18) and replaced by the centring policy**, in Hans's words: "(1) centering
   comes from the solve2D run prior to solve3D on both particles used for training and particles being predicted,
   (2) shift refinement only enters after the prediction on the reference obtained from the training set, thus
   mirroring the class student. This is the policy most consistent with what we want to accomplish in the stream with
   the HHCNN." (In the stream the unseen particles' 2D shifts come from the class student's own alignment after its
   prediction; in this run they are the solve2D shifts, as C0's.) The raw-pick route was not adopted; its tables of
   Phases 3 and 4 stay below as the reason.
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
   - Set from the data (ruling R2, Phase 0): trpv1 `pgrp=c2 pgrp_start=c2`, from the project's `5_solve3D` job;
     streptavidin `pgrp_start=c1`, from its `4_solve3D` job. The uncentred test half: each of S's 24,000 truth shifts
     moved by a displacement whose radius is resampled from trpv1's picking error (the radii of the converged solve2D
     shifts of all 35,645 particles, `class_student/teachers/trpv1_full`: median 7.98, 90th percentile 15.95, 99th
     39.98 Å) and whose direction is uniform (seed 20261116); the stored shift carries the displacement (the particle
     preparation shifts before it masks, so this is a displaced pick up to the box's wrap), built with `new_project`,
     `import_particles` and `selection` as the earlier decentred halves.
6. **Arms.** C0: today's pose student (teacher at `prob_beta` x T, today's registration) on unseen particles at their 2D
   shifts, the yardstick. R0: the same student on raw picks (what dropping the 2D shifts alone does). Then on raw picks
   with the alignment after prediction: A (today's teacher) and P3 (the physical teacher), each with S1 and SK.
   (Phase 0: C0 and R0 are one model predicted twice, so they pair exactly. Both run `predict_poses` with centring off,
   learn's setting: today's `predict_poses` given the map leaves `center` at its default and, for a cyclic point
   group, re-centres the map and moves every particle's stored shift by the offset, which ruins its prediction on the
   c1 truth sets: 2.4% of modes within 6° of the truth instead of 75.6%. Phase 1 fixes the commander.)
7. **Reads**, on the unseen particles, against the reference (decision 2) and on the truth sets also against the truth:
   - seeded (route 1) share; the cone holding the reference's direction (all particles); the mode's agreement;
   - after the alignment: the committed direction within 6° of the reference (and of the truth), the in-plane angle's
     error, the shift's distance to the reference's shift (median and 90th percentile, Å);
   - the unseen map: `reconstruct3D` from the committed poses, its FSC against the teacher's map (cross-validated: the
     two share no particle), FSC 0.143 and 0.5;
   - cost per particle: prediction, alignment, their sum against the reference's dense search.
   - context (Phase 0): D0, the dense search itself on the raw picks (the reference's pass started from the picked
     positions, no student), read as a committed pose. The reference's direction is a grid direction (about 6° apart),
     so "within 6° of the reference" is close to "the same grid direction"; 12° is reported beside it.
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
   **Revised by ruling R5 (Hans, 2026-11-17 evening): the map judges.** A student that regularizes can make maps as good
   as the full search's while agreeing with its single answer less often, so (b) and (c) are decided on the unseen map:
   on the truth sets its FSC at 0.143 and 0.5 against the truth map, on the real sets its FSC against the teacher's map
   read beside the full search's own map from the same particles (the reference's poses) as the yardstick; paired by
   seed, worse or better beyond twice the paired standard error. (b): P3 is adopted unless its map is worse beyond 2 SE
   on a set. (c): the raw-pick route against C0 on the unseen map per real set and per truth set. The agreement reads of
   decision 7 (the cone holding the reference, the committed direction within 6° or 12°, the mode's agreement, the
   shift and in-plane errors) are reported as diagnostics and decide nothing. (a) stands: SK also gave the better map
   on every set.
   **(c) revised by ruling R7 (2026-11-18)**: the arm judged is C3SK, the physical teacher's students of Phase 3 (no
   retraining) predicted on the unseen particles at their 2D shifts (C0's input) and aligned over the cone against the
   teacher's map, on `fd9af0add`'s build so that it pairs with C0 and P3SK; C3SK is adopted unless its unseen map is
   worse than C0's beyond 2 SE on a set (real sets against the teacher's map, beside the full search's own map; truth
   sets against the truth map). If adopted, the tree implements the policy and nothing else: no raw-pick or
   uncentred-half path left in the code.
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
- Off-state identity: `simulated_workflow` suite 6vxx (ruling R4) on the base and the final tree at one thread, seed fixed,
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
  arms, decision 7's reads, a report that says what the cone and the alignment are. Testers and plants. The
  commander's `predict_poses` sets `center=no` as learn's does (Phase 0's finding, decision 6). The frozen teachers
  (Phase 0) enter `validate_pose_student` through the student's commander (ruling R3).
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
| `doc/algorithms/pose_student.md`, `doc/policies/test_environment_policy.md` | the route, the test (Phases 2 and 3: the losing alignment and the losing teacher removed from the description; ruling R3) | 1, 2, 3, 4 |
| `src/main/ui/simple/simple_ui_solve3D.f90` | `validate_pose_student`'s help and its frozen teacher (`dir_refine`, an existing key; ruling R3, Phase 1) | 1, 4 |
| `doc/implementation_notes/README.md`, `doc/implementation_notes/planned/polar_student_future_directions.md` | record | 0, 1, 2, 3, 4 |

## Progress

Ruling R4 (Hans, 2026-11-17) moved the base to `fd9af0add` at the start of Phase 2: the two entries that follow were
made on the launch base `1ea0f8340` and are **superseded** (a record: the old base's teachers inserted each pose with
the in-plane posterior, the new base's at its own angle). Their artefacts are under `scratch/keep/old_1ea0f8340/`. Phase 0
is redone on the new base (its entry after them), then Phase 1 re-applies its work onto the new base and reruns its gates.

### Phase 0, the baseline (2026-11-16 to 2026-11-17, the Dell, done; superseded by ruling R4)

**The base and the machine.** HEAD `1ea0f8340051911e069f7087da7d3165064f6bec` with the launch's one uncommitted change
(this plan's status line); no source file differs from HEAD, so the source the run builds is `git archive 1ea0f8340`,
exported to `scratch/keep/src_base`. The Dell: Oracle Linux 8.10, Intel Core Ultra 9 285 (24 cores, no
hyperthreading), 62 GB, GCC 15.2.1 of gcc-toolset-15, CMake 3.26.5. Run root `/data2/agent_runs/pose_raw` (`R` below;
paths under `R/scratch/keep/` unless absolute). The maintainer's launch rulings apply: the run builds and runs on the
Dell without following master, never uses the maintainer's checkout or its build, and records the machine's load with
every timing. Every script points `SIMPLE_PATH` and `PATH` at the tree it means (`tools/env.sh <tree>`), since the
login shell points them at the maintainer's build. The one-minute load before the run started was 2.2.

**Builds.** Configured as `cmake -S src_base -B <tree> -DCMAKE_BUILD_TYPE=<Release|debug> -DBUILD_TESTS=ON -DNICE=OFF
-DUSE_ARCHOPT=ON -DFORTLS=OFF -DPython3_EXECUTABLE=/usr/bin/python3.12`, `make -j20`, `make install`
(`tools/build_base.sh`): `build_base` (Release, make 90 s) and `build_base_debug` (Debug, 45 s). The fast gate
(`scripts/run_fast_gate.sh`): 15 of 15 on both (9.8 s and 36.5 s; `R/scratch/logs/build_base{,_debug}.gate.log`).

**Off-state identity, the base against itself.** `simulated_workflow` (the end-to-end simulated test: solve2D's
refine2D stages, then solve3D's refine3D stages in shared memory and over two parts), both suites, at one thread and
`SIMPLE_SEED=7`, twice on `build_base` side by side (`identity/base1`, `identity/base2`; 1jxy 15:59 and 14:08,
6vxx 12:03 and 12:05 minutes, load 4 to 7 from the four jobs); `test=identity` of `base1` against `base2`
(`tools/identity_cmp.sh`): 1jxy 223 files the same and 0 differing, 6vxx 186 and 0
(`identity/cmp_base1_base2/<suite>/identity_report.txt`). **The exclusion lists are empty for both suites.** 1jxy's
shared-memory solve3D misses its pose gate in both runs (pair pose error 19.76° against 15°), the known miss at this
seed. `base1` stays for Phases 1 and 4; `base2`'s outputs were deleted, its logs kept.

**The truth sets.** The simulated medium-assembly sets of the earlier runs, read in place: the training sets T200,
T1000, T5000 and T (24,000) with their truth poses, the truth map `T/tmap.mrc`, the test half S (24,000; box 128,
3.0 Å per pixel, 300 Å mask, c1). Their `MD5SUMS` verify (the one failing line is the file's known stale first line
for `S/S.simple`; its second line holds); the run's own list, stacks included, is `truth_sets.md5`. On the truth sets
the teacher is the truth itself, as in the earlier runs: the student learns the training set's truth poses against
the truth map.

**The uncentred test half** (decision 5; `tools/make_S_unc.sh`, `sets/S_unc/`, `MD5SUMS` beside it). Each of S's
truth shifts moved by a displacement drawn from trpv1's picking error, as decision 5 now records: displacements of
median 7.83, 90th percentile 15.92, 99th 40.12 and largest 84.2 Å; every other key of every row equals S's
(`check_S_unc.log`). Built beside it: `S_check` (S rebuilt by the same route, equal to `S.simple` on every key) and the
reference passes' copies with the disk flag set (`shsig=-1`, a particle not yet centred).

**The teachers** (decision 5): base `validate_pose_student` at `SIMPLE_SEED=<t>`, t = 1 to 3, per real set
(`tools/teacher.sh`): its class-balanced greedy sample, its solve3D on that sample (the teacher), then the base's own
learn, prediction and maps, kept as a check. Each teacher was then frozen (`tools/freeze_teacher.sh`,
`teachers/<set>_t<t>/`, `MD5SUMS`): the solve3D result project, its map and halves, the FSC records, the base
validation's report and model; reconstructions and intermediate volumes deleted. Command lines: bgal
`nsample=1500 mskdiam=180 pgrp=d2 pgrp_start=c1`, trpv1 `nsample=5000 mskdiam=175 pgrp=c2 pgrp_start=c2`,
streptavidin `nsample=5000 mskdiam=85 pgrp=d2 pgrp_start=c1`, all `nthr=24`. Unseen particles: bgal 1,256 of 2,756
selected, streptavidin 1,821 of 6,821, trpv1 30,504 of 35,504.

| set | teacher seed | solve3D (s) | final FSC 0.5 / 0.143 (Å) | validation wall (min:s) | load at start / end |
| --- | --- | --- | --- | --- | --- |
| bgal | 1 | 495 | 7.10 / 4.13 | 9:37.93 | 22.31 / 11.06 |
| bgal | 2 | 507 | 7.59 / 4.18 | 9:47.95 | 1.20 / 12.90 |
| bgal | 3 | 835 | 7.10 / 4.24 | 15:13.20 | 12.90 / 13.25 |
| strep | 1 | 848 | 3.61 / 3.19 | 15:15.09 | 13.25 / 18.84 |
| strep | 2 | 882 | 3.61 / 3.19 | 15:49.81 | 18.84 / 22.75 |
| strep | 3 | 847 | 3.61 / 3.12 | 15:11.70 | 22.75 / 20.47 |
| trpv1 | 1 | 1498 | 8.89 / 7.59 | 28:21.43 | 20.47 / 13.56 |
| trpv1 | 2 | 1480 | 9.43 / 7.41 | 28:03.33 | 13.56 / 16.84 |
| trpv1 | 3 | 1424 | 9.15 / 7.78 | 27:25.30 | 15.57 / 13.55 |

bgal's third teacher ran beside the truth harness's smoke run and a debugging run (its solve3D 835 s against 495 and
507 s); the others ran with nothing else of the run's on the machine, the one-minute loads being the runs' own. The
checksums of every frozen teacher still hold after Phase 0's 27 validations read them.

**How the arms were run.**
- **Real sets: one variant build, V0** (`variants/V0/V0.patch`, 625 lines, on a copy of the export in
  `variants/V0/src`, Release in `variants/V0/build`, fast gate 15 of 15). It changes `validate_pose_student` and adds
  one early exit to learn's rows; the student's training and prediction are the base's. Given a frozen teacher's
  solve3D directory (`POSE_TEACHER_DIR`), the validation skips its sample's solve3D and checks the teacher was trained
  on the same sample; the learn at seed s then reproduced the base validation's own model bit for bit (bgal teacher 1,
  seed 1: the weights files compare equal). It computes the reference (decision 2, as recorded there: the dense table
  at the physical temperature from the 2D solution's shifts, cached beside the teacher and frozen) and D0 (the same
  pass from the raw picks); learn's rows suffice, so the variant stops learn after its rows when the model file is the
  reference's. It predicts the unseen particles twice with one model: C0 at their 2D shifts and R0 at their picked
  positions (shift zero), both with the reference's pose stored so that the prediction table's reference columns read
  against it, both with centring off (below). It writes the reads (`V0READS` lines of the report, the per-particle
  `v0_details.txt`) and four maps: C0's poses at the 2D shifts and at the registered shifts, R0's at its registered
  shifts, and the reference's own. One run per real set, teacher t and training seed s (`tools/v0run.sh <set> <t>
  <s>`): `POSE_TEACHER_DIR=<teacher>/1_validate_pose_student/1_solve3D POSE_REF_CACHE=<teacher>/reference
  SIMPLE_SEED=<s> simple_exec prg=validate_pose_student projfile=<copy of the selection project> <the teacher's
  arguments> nthr=24`, kept in `v0/<set>_t<t>_s<s>/kept/`.
- **Truth sets: the harness of the earlier runs** (`tools/truth.sh <set> <seed>`): `learn_pose_student` on the
  training set (`vol1=tmap.mrc lp=15 nspace=1000 shift_sigma=2 mskdiam=300 pgrp=c1 nthr=24`, `SIMPLE_SEED=<seed>`),
  then `predict_poses ... vol1=tmap.mrc center=no` four times: C0 on S, R0 on the uncentred half, each with the
  reference's pose stored and with the truth's; then the unseen maps of C0 and R0 (`reconstruct3D` at the readout
  poses and registered shifts, `fsc` against the truth map). The references (`tools/truth_ref.sh`,
  `tools/truth_ref_unc.sh`, `truth_ref/`): one pass from S's centred shifts serves both arms (the uncentred half holds
  S's images); the pass from the uncentred half's shifts is D0. Seeds 1 to 8 for each of the four training sets.
- C0 and R0 share the model, so each pair differs by the input shifts alone.

**A defect of today's `predict_poses`, found here.** Given the map (`vol1`, the shift registration), `predict_poses`
leaves `center` at its default `yes` (learn sets `no`); for a cyclic point group the reprojection model's setup then
re-centres the map and moves every particle's stored shift by the map's centre-of-mass offset, which the prediction
then reads and the project keeps. On the c1 truth sets (the truth map lies about 9 pixels off centre) it ruined the
prediction: 2.4% of S's modes within 6° of the truth with the map against 75.6% without it, and 1.1% against 91.7% on
T1000's own training rows; with `center=no` the prediction with the map equals the one without it row for row, and
equals learn's own held-out predictions (100 of 100). bgal and streptavidin (d2) are not centred by that path, and
trpv1's maps are centred (its training rows read 0.59 in the base's own validation). The arms run with `center=no`;
the commander's fix is Phase 1's (section 8, decision 6). The diagnosis: `scratch/tmp/vol1_bug/` (runs A, B, C) and
`scratch/tmp/pred_dbg/`.

**Two further observations for the maintainer.** The dense table's rows are not reproducible bit for bit: two passes
over S (the base build's learn and V0's, the same source compiled twice) gave different modes for 74 of 24,000 rows, so each reference is
computed once and frozen. And on trpv1 the base student is weak: held-out mode agreement 0.13 and a threshold above
every confidence, so it seeds nobody (the same on streptavidin); on bgal it seeds 4% of the unseen particles.

**What each read is** (over the unseen particles; the reference's direction is a grid direction, about 6° from its
neighbours). Seeded: the share the student would send to the neighbourhood search (its cone within the cap of 64
directions and the bound, its confidence at the threshold). Capped: the share whose 99% cone needs more than 64
directions. Cone holds the reference: the share whose folded cone contains the reference's direction (mirror pairs
folded). Mode within 6° or 12°: the student's most likely direction pair against the reference's (mirror and point
group folded). Committed direction within 6° or 12°: the readout's full pose (direction with the readout's mirror,
in-plane angle; no search) against the reference's, point group folded, mirror not; the in-plane error where it is
within 6°; the shift's distance to the reference's (median and 90th percentile, Å), the registered shift (today's 5 by
5 grid over 25 Å then the gradient refinement at the readout pose). The unseen map: reconstruct3D from the committed
poses and registered shifts, its FSC against the teacher's map at 0.143 and 0.5 (cross-validated: no shared
particle). D0 is read like a committed pose. On the truth sets every direction read also against the truth.

**Real sets, C0 → R0** (mean over the nine teacher-and-seed pairs of each set; the paired difference R0 − C0 and its
standard error in brackets; `reads_real.tsv`, `tables_p0.md`):

| read | bgal | strep | trpv1 |
| --- | --- | --- | --- |
| seeded (%) | 4.2 → 0.8 (-3.5, 0.7) | 0.0 → 0.0 (+0.0, 0.0) | 0.0 → 0.0 (+0.0, 0.0) |
| cone capped at 64 (%) | 39.3 → 86.0 (+46.6, 1.5) | 100.0 → 100.0 (+0.0, 0.0) | 93.0 → 98.4 (+5.4, 0.5) |
| cone holds the reference (%) | 91.3 → 48.6 (-42.8, 0.4) | 70.8 → 20.6 (-50.1, 0.5) | 90.3 → 64.8 (-25.5, 0.7) |
| mode within 6° (%) | 70.3 → 20.4 (-50.0, 0.4) | 16.9 → 3.5 (-13.4, 0.3) | 40.8 → 19.3 (-21.5, 0.3) |
| mode within 12° (%) | 86.5 → 37.3 (-49.1, 0.4) | 45.2 → 12.1 (-33.1, 0.6) | 62.2 → 35.9 (-26.3, 0.5) |
| committed within 6° (%) | 68.9 → 19.0 (-49.9, 0.4) | 13.7 → 2.2 (-11.5, 0.3) | 34.7 → 14.1 (-20.6, 0.3) |
| committed within 12° (%) | 85.6 → 34.3 (-51.3, 0.4) | 38.4 → 8.5 (-30.0, 0.5) | 50.1 → 24.6 (-25.6, 0.6) |
| committed, median (°) | 3.6 → 21.2 (+17.6, 0.3) | 15.9 → 43.4 (+27.4, 0.3) | 12.1 → 39.6 (+27.5, 0.6) |
| in-plane error, median (°) | 1.60 → 2.25 (+0.65, 0.06) | 4.17 → 77.66 (+73.49, 7.13) | 2.29 → 3.50 (+1.21, 0.20) |
| shift to the reference, median (Å) | 0.4 → 8.0 (+7.6, 0.2) | 1.6 → 12.6 (+11.0, 0.1) | 1.7 → 5.7 (+4.0, 0.1) |
| shift to the reference, q90 (Å) | 6.7 → 22.4 (+15.8, 0.1) | 6.6 → 37.3 (+30.7, 0.3) | 6.7 → 14.8 (+8.1, 0.2) |
| unseen map, FSC 0.143 (Å) | 4.90 → 8.42 (+3.52, 0.16) | 6.96 → 13.44 (+6.48, 0.50) | 7.92 → 8.60 (+0.68, 0.12) |
| unseen map, FSC 0.5 (Å) | 7.79 → 13.67 (+5.88, 0.10) | 9.66 → 13.72 (+4.06, 0.10) | 9.13 → 14.08 (+4.95, 0.11) |

**D0, the dense search on the raw picks**, against the reference (one pass per teacher, the mean of three):

| set | within 6° (%) | within 12° (%) | median (°) | shift to the reference, median / q90 (Å) | C0's committed within 6° (%) |
| --- | --- | --- | --- | --- | --- |
| bgal | 52.7 | 71.3 | 4.4 | 0.4 / 9.0 | 68.9 |
| strep | 5.7 | 13.5 | 40.3 | 5.2 / 35.6 | 13.7 |
| trpv1 | 54.7 | 61.3 | 4.7 | 0.8 / 9.2 | 34.7 |
| truth sets (S's images from the uncentred half's shifts; within 6° of the truth 72.5, median 3.6°) | 64.1 | 84.4 | 0.0 | 0.4 / 7.3 | - |

**Truth sets, C0 → R0** (eight training seeds each; `reads_truth.tsv`):

| read | T200 | T1000 | T5000 | T |
| --- | --- | --- | --- | --- |
| seeded (%) | 21.8 → 5.7 (-16.1, 1.5) | 25.7 → 6.3 (-19.3, 1.0) | 45.3 → 11.1 (-34.2, 1.2) | 53.9 → 10.5 (-43.4, 1.6) |
| cone capped at 64 (%) | 21.2 → 52.4 (+31.1, 1.2) | 18.6 → 55.3 (+36.7, 1.9) | 6.9 → 45.7 (+38.8, 1.5) | 3.2 → 57.4 (+54.2, 1.1) |
| cone holds the reference (%) | 97.0 → 59.9 (-37.1, 0.5) | 98.0 → 63.2 (-34.9, 0.4) | 98.9 → 65.2 (-33.7, 0.2) | 99.9 → 70.3 (-29.5, 0.1) |
| cone holds the truth (%) | 96.5 → 59.3 (-37.3, 0.5) | 97.7 → 62.6 (-35.1, 0.4) | 98.4 → 64.5 (-33.8, 0.2) | 99.5 → 69.8 (-29.7, 0.1) |
| mode within 6° of the truth (%) | 74.2 → 23.5 (-50.7, 0.7) | 78.7 → 26.3 (-52.4, 0.5) | 87.6 → 31.2 (-56.4, 0.4) | 93.4 → 36.9 (-56.5, 0.3) |
| committed within 6° of the reference (%) | 50.5 → 12.6 (-37.9, 0.8) | 53.3 → 14.2 (-39.0, 0.7) | 61.4 → 17.5 (-43.9, 0.6) | 68.5 → 22.1 (-46.4, 0.6) |
| committed within 12° of the reference (%) | 82.9 → 37.8 (-45.2, 0.4) | 88.1 → 41.5 (-46.6, 0.4) | 95.3 → 47.4 (-47.9, 0.2) | 99.2 → 51.8 (-47.5, 0.3) |
| committed within 6° of the truth (%) | 74.1 → 23.2 (-50.8, 0.7) | 78.7 → 26.0 (-52.6, 0.5) | 87.5 → 31.0 (-56.6, 0.4) | 93.4 → 36.6 (-56.7, 0.3) |
| committed against the truth, median (°) | 3.8 → 25.0 (+21.2, 0.8) | 3.6 → 18.0 (+14.3, 0.6) | 3.3 → 12.4 (+9.1, 0.3) | 3.1 → 10.1 (+7.0, 0.1) |
| shift to the truth, median (Å) | 1.2 → 6.3 (+5.1, 0.1) | 1.1 → 5.7 (+4.5, 0.1) | 1.0 → 4.8 (+3.8, 0.1) | 1.0 → 4.3 (+3.4, 0.0) |
| shift to the reference, q90 (Å) | 3.9 → 20.7 (+16.8, 0.1) | 3.2 → 19.8 (+16.6, 0.1) | 2.6 → 18.4 (+15.7, 0.1) | 2.3 → 17.2 (+14.9, 0.0) |
| unseen map, FSC 0.143 (Å) | 6.06 → 8.54 (+2.48, 0.10) | 6.04 → 7.93 (+1.90, 0.14) | 6.00 → 7.53 (+1.53, 0.15) | 6.00 → 7.15 (+1.15, 0.08) |
| unseen map, FSC 0.5 (Å) | 8.47 → 13.30 (+4.83, 0.42) | 8.16 → 11.73 (+3.57, 0.10) | 7.61 → 11.46 (+3.85, 0.08) | 7.47 → 11.02 (+3.55, 0.13) |

The reference itself against the truth (S): 99.1% within 6°, median 2.70°, its shift 0.93 Å from the truth's at the
median (1.91 at the 90th percentile).

**Costs per particle** (ms, wall time at 24 threads unless marked; the prediction's parallel region from
`predict_poses`' bench, the forward pass and today's registration in thread-milliseconds; the dense passes' table fill
and their setup, particle preparation and fill together; mean over the runs):

| set | prediction, C0 | prediction, R0 | forward pass (thread-ms) | registration (thread-ms) | reference: fill | reference: all | D0: fill | D0: all |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| bgal | 0.432 | 0.454 | 4.280 | 0.391 | 3.457 | 7.196 | 3.920 | 6.524 |
| strep | 0.079 | 0.079 | 0.616 | 0.123 | 2.200 | 2.512 | 2.320 | 2.603 |
| trpv1 | 0.194 | 0.197 | 3.119 | 0.447 | 3.167 | 3.546 | 3.380 | 3.740 |

On the truth sets the whole `predict_poses` run costs 0.46 to 0.61 ms per particle of S (setup included) and the
reference's table fill 3.0 ms.

**What it shows.**
- **Dropping the 2D shifts alone costs every read on every set, by far more than its noise.** The cone's hold on the
  reference falls by 25 to 50 points on the real sets (to 21 to 65%) and by 30 to 37 on the truth sets (to 60 to 70%);
  the committed direction within 6° by 12 to 57 points; the unseen map coarsens by 0.7 to 6.5 Å at FSC 0.143 and by
  3.5 to 5.9 Å at 0.5. Every paired difference lies beyond twice its standard error (the shares' paired errors are 0.1
  to 1.9 points: C0 and R0 share their model, so their difference is the shift's alone).
- **On raw picks the cone still holds the reference for most particles where the student is informative**: 65% on
  trpv1, 60 to 70% on the truth sets, 49% on bgal, against a mode within 6° of the reference for 19% (trpv1) and 20% (bgal) and
  within 6° of the truth for 23 to 37% (the truth sets). That gap is
  what an alignment over the cone (SK) can recover and the mode alone (S1) cannot. Streptavidin is the exception (21%):
  there even C0's student is weak (its cone capped for every particle, its mode within 6° for 17%).
- **The raw picks' cones are wide**: capped at 64 directions for 46 to 57% of the truth sets' particles and 86 to 100%
  of the real sets' (C0: 3 to 39%, except trpv1 and streptavidin, already above 90%). SK will mostly search the cap's
  64 directions and their mirrors, so its cost is set by the cap.
- **Today's registration does not centre a raw pick.** The registered shift lies a median 6 to 13 Å from the
  reference's on the real sets (C0: 0.4 to 1.7 Å) and 4.3 to 6.3 Å from the truth's on the truth sets, with a 90th
  percentile of 15 to 37 Å: the class student's finding again.
- **The dense search itself on raw picks (D0) is no remedy**: its mode agrees with the reference for 53 to 55% on bgal
  and trpv1, 6% on streptavidin and 64% on the truth sets' images (73% within 6° of the truth, against 99% from the
  centred start). On trpv1 it still beats C0's committed direction (55% against 35%): the base student is weak there
  (it seeds nobody on trpv1 and streptavidin, 4% of bgal's unseen particles).
- **Costs**: the prediction takes 0.08 to 0.45 ms per particle (wall, 24 threads), R0 the same as C0; the reference's
  dense table fill takes 2.2 to 3.5 ms. The alignment's cost is Phase 2's measurement.

**Decision 9** (if R0's cone holds the reference's direction for fewer than half of the particles on every set, skip
SK). R0's cone holds the reference for 48.6% on bgal and 20.6% on streptavidin, but for 64.8%
on trpv1 and 59.9, 63.2, 65.2 and 70.3% on T200, T1000, T5000 and T: not below half on every set, so the early stop
does not apply. SK is measured in Phase 2 as planned.

**Exit.** Each item of Phase 0's exit and its evidence:
- The base pinned, Release and Debug builds, the fast gate on both: above (`R/scratch/logs/build_base*.log`).
- The truth sets' checksums: `truth_sets.md5` and the earlier `MD5SUMS`, verified.
- The off-state identity of the base against itself, both suites: empty exclusion lists
  (`identity/cmp_base1_base2/`).
- The uncentred test half built and frozen with checksums: `sets/S_unc/MD5SUMS`, `check_S_unc.log`.
- The nine real-set teachers built and frozen, each solve3D's wall time and final FSC recorded: the teachers' table,
  `teachers/<set>_t<t>/MD5SUMS` (re-verified after Phase 0's validations), the references frozen beside them
  (`teachers/<set>_t<t>/REFERENCE_MD5SUMS`; the truth sets' in `truth_ref/`).
- C0 and R0 on every set at their seeds: 27 real-set validations (three sets, nine teacher-and-seed pairs) and 32
  truth-set runs (four training sets, eight seeds), no failure (`scratch/logs/chainB.log`, `summary.tsv`).
- The reads of decision 7 that exist without the alignment, the cone's hold on the reference included: the tables
  above, frozen in `reads_real.tsv`, `reads_truth.tsv` and `tables_p0.md`.
- Decision 9: not triggered.

Files changed in the repository: this plan only (decisions 2, 5, 6 and 7 annotated with what Phase 0 found and how it
was done, section 8's Phase 1 with the centring fix and the frozen teachers, and this entry). Kept for the later
phases under `R/scratch/keep/`: the base export and both base builds, `identity/base1`, the frozen teachers with their
references, the truth sets' references and projects (`truth_ref/`, `sets/S_unc/`), the V0 variant (its patch is the
draft of Phase 1's frozen-teacher path and reference), the harness (`tools/`), the runs' kept outputs (`v0/`,
`truth/build_base/`) and the tables. Bulky intermediates (reprojection models, feature caches, maps of finished runs,
project copies) were deleted after their numbers were recorded. Not done in Phase 0: the run's own trees
(`repo/build`, `repo/build_debug`), whose source equals the base until Phase 1 changes it.

### Phase 1, the physical teacher, the reference and the alignment (2026-11-17, the Dell, done; superseded by ruling R4)

**What the tree now does.**
- **The alignment after prediction** (`simple_pose_student_align.f90`, new). Given the teacher's map (`vol1`),
  `predict_poses` treats every particle as a raw pick, not yet centred, and aligns it at its picked position against
  the map's reprojections at the student's grid and band: every in-plane rotation of the polar grid and shifts over a
  30 Å disk (three widths of a 10 Å prior) about the stored shift, the classic searches' polar shift map
  (`shiftmap_srch`, unchanged, as the class student's alignment uses it), at the physical temperature. Two lists are
  aligned on the same prediction, each timed on its own: S1, the mode's mirror pair (two references); SK, every
  direction of the cone with its mirror partner (up to 128 references). SK's best alignment is committed to the
  particle's record (direction, in-plane angle, stored shift plus the found increment, and the dynamic centring's state
  for a following refine3D); S1's pose and the readout's are in the prediction table beside it, the walls in
  `<fname>_times.txt`. Phase 2 decides between S1 and SK and removes the loser. Today's shift registration (a 5 by 5
  grid over 25 Å at the readout's pose) is gone, with its two constants.
- **Centring off in `predict_poses`** (Phase 0's finding): the commander sets `center=no` as learn does; with the map it
  also sets the alignment's search keys (`refine=prob`, `trs`, `shift_sigma` 10 Å, `prob_beta` 1).
- **The physical teacher P3** behind one named constant, `POSE_STUDENT_PHYSICAL_TEACHER`, off in the tree (today's
  teacher A, which Phase 2 runs); a variant build flips it for Phase 3 (`variants/P3/P3.patch`). Under P3 learn's
  commander sets `prob_beta` to 1, the class student's pattern. Every row learn reads is audited against the physical
  posterior recomputed from its losses at the table's own temperature (`physical_row_audit`): the training's targets
  must equal `exp[-T (L - L_min)]` normalised, within the rounding bound of the single-precision losses (twice T times
  their spacing). A physical teacher stops on a row read at another temperature or beyond the bound; today's teacher
  logs its departure. The rows header records the audit.
- **The reference** (decision 2, as Phase 0 recorded it): `run_pose_posterior_mode`, learn's rows with no student at
  `prob_beta` 1 over the 30 Å disk; each row's mode with its in-plane angle and refined shift goes to the particle
  field and to `<base>_poses.txt`.
- **`validate_pose_student`** rewritten. It takes a frozen teacher through the existing key `dir_refine` (a solve3D
  directory of an earlier validation; the sample must match), and the reference and the dense search on the raw picks
  (D0) from `pose_reference.txt` and `pose_dense_raw.txt` in that directory when present (Phase 0's frozen ones were
  copied there), else it makes them. The unseen particles are predicted at their picked positions only. The report
  gives, per arm (the readout RO, S1, SK, D0), the reads of decision 7 against the reference with a `READS` line each,
  the per-particle details, the costs per particle (prediction, each alignment, the dense search) and three maps (SK's,
  S1's and the reference's poses). `validate_pose_student`'s help and its `dir_refine` input are in the UI (a file-table
  row added under ruling R3). `doc/algorithms/pose_student.md` describes the alignment, the physical teacher with its
  audit, and the validation.
- **A defect fixed on the way**: `predict_poses`' distributed master bootstrapped sigma2 over its parts, summing the
  spectra in another order than a shared-memory run; the old registration hid it, the alignment did not (the
  `pose_student` workflow test's master-against-shared-memory checks failed on one particle's S1 shift and on the
  records' shifts). The master now bootstraps in its own process; its table and sigma2 then equal the shared-memory
  run's byte for byte (`scratch/tmp/p1/order_dbg`, one and eight threads also identical).

**Testers and plants.** `test=unit_numerics suite=pose_student_alignment` (new, five tests, 230 assertions on a fixture
of four white references in two mirror pairs, particles planted rotated and shifted 5 px beyond a 2 px window): the
planted reference recovers its rotation and shift over the disk; a list naming only the planted reference's mirror
partner finds it; the cone in either order keeps the planted direction and another pair alone scores lower; a
reference named three times is searched once, and a list beyond the object's capacity aligns to nothing (the caller's
status after the parallel region); a team of four aligns as one thread, bit for bit. `suite=pose_student` gains the
audit's test (rows at the physical temperature within the bound; the same rows read at half or twice it beyond; invalid
entries without mass). Nine plants, each a copy of the tree with one check's violation planted, rebuilt and its suite
run (`tools/plants_p1.py`, `scratch/keep/plants/plants.log`, the patches beside it): the mirror not searched, the first
reference kept, the window instead of the disk, the refinement ignored, a reference searched twice, a thread-dependent
score, the audit's targets tempered, the physical masses unnormalised, a bound blind to the temperature. **All nine are
caught.**

**Gates on the final code.**
- Builds: the run's Release tree `repo/build` and Debug tree `repo/build_debug`, configured as the base trees
  (Phase 0's configure line), `make`, `make install`.
- The fast gate (`scripts/run_fast_gate.sh`): 15 of 15 on both trees. `test=unit_numerics` (32 groups, 2,423
  assertions) and `test=unit_pftc_align2D3D` (15 groups, 6,384 assertions): passed on both trees
  (`scratch/keep/gates_p1/`, `tools/gates_p1.sh`).
- CTest `pose_student` and `class_student` on the Release tree: both passed (36.6 s and 14.6 s;
  `scratch/keep/ctest_p1.log`), the former after the sigma2 fix above.
- **The off-state identity**: `simulated_workflow` both suites at one thread and `SIMPLE_SEED=7` on the final tree
  (`identity/tree`; 1jxy 15:28, 6vxx 14:04 minutes) against Phase 0's `base1` (`test=identity`,
  `identity/cmp_base1_tree/`): 1jxy 223 files the same and 0 differing, 6vxx 186 and 0. 1jxy's known pose-gate miss
  as on the base.
- **The class student's identity** (`tools/guard_p1.sh`, `scratch/keep/guard_p1/`): `learn_class_student` on bgal (the
  frozen 2D teacher of the student_inplane run) at seed 1, 24 threads, on the base and on the tree: the model headers
  equal and the weights' checksums equal. The same for `learn_pose_student` on bgal under today's teacher (the header's
  new audit keys aside): the training is the base's, bit for bit.

**Smoke runs** (one seed each, read only as checks of the machinery; the runs ran beside the identity jobs, so their
times are inflated):
- Truth set T1000, seed 1 (Phase 0's model; `scratch/tmp/p1/align_smoke`): on the uncentred half, SK's committed
  direction lies within 6° of the reference for 62.8% of S's 24,000 particles and within 6° of the truth for 72.1%
  (S1 12.8% and 24.0%, the readout 12.7% and 23.8%); SK's shift lies 0.39 Å from the reference's at the median.
- bgal, teacher 1, seed 1 (`tools/vrun.sh`, `scratch/keep/val/smoke/bgal_t1_s1`): the frozen teacher and references
  are read; the prediction equals Phase 0's R0 (the cone holds the reference for 0.4984 in both); SK's committed
  direction within 6° of the reference 59.4%, S1 20.5%, D0 52.4% (Phase 0's C0: 69.0%); the unseen map at SK's poses
  4.41 Å at FSC 0.143 against 4.35 Å at the reference's poses (Phase 0's C0 at this seed: 5.26 Å).
- The physical teacher (the P3 variant, learn on bgal at seed 1; `scratch/keep/p3_smoke`): 1,500 rows, none off the
  physical temperature, the largest mass difference 3.7e-8 (0.0002 of the rounding bound); today's teacher's rows
  depart by up to 0.60 in mass.

**For Phase 2.** `tools/vrun.sh <label> <tree> <set> <t> <s>` runs the tree's validation on a frozen teacher;
`tools/truth_tree.sh <label> <tree> <set> <seed>` the truth sets (the uncentred half with the reference's poses, SK's
and S1's maps); `tools/aln_reads_truth.py` their reads. SK costs 3.4 ms per particle on S (24 threads, wall) against
0.06 ms for S1 and 3.0 ms for the dense table's fill: its cost is set by the cone cap, since most raw picks' cones are
capped (Phase 0).

**Exit** (section 5's Phase 1 gates): the testers of the alignment (a planted shift recovered over the disk, the
mirror of a cone direction searched, SK keeps the best) and the teacher's audit with their plants; the fast gate,
`unit_numerics` and `unit_pftc_align2D3D` on both trees; the `pose_student` and `class_student` workflow tests; the
off-state identity; the class student's identity. All met, as listed above. Files changed: `src/defs/simple_defs_hhcnn.f90`;
`src/main/learn/simple_pose_student.f90` and its tester; `src/main/strategies/parallelization/simple_strategy_pose_student.f90`;
`simple_pose_student_align.f90` and its tester (new); `src/main/commanders/simple/simple_commanders_pose_student.f90`;
`src/main/commanders/test/simple_commanders_test_class.f90`; `src/main/ui/simple_test/simple_test_ui_class.f90`;
`src/main/ui/simple/simple_ui_solve3D.f90`; `doc/algorithms/pose_student.md`; this plan;
`doc/implementation_notes/README.md`.

### Phase 0 on the new base, the baseline (2026-11-17, the Dell, done; ruling R4)

**The base and the machine.** HEAD `fd9af0add88cc50d677162e1cd811c2e4c056d3b` (ruling R4: the 3D in-plane posterior
insertion and `simulated_workflow_1jxy` removed), with the run's carried documents (this plan and its row in
`doc/implementation_notes/README.md`); no source file differs from HEAD, so the source built is `git archive fd9af0add`
in `scratch/keep/src_base`. The machine, the run root, the environment script and the rulings on the machine's use are
as in the superseded entry above; the old base's artefacts moved to `scratch/keep/old_1ea0f8340/` (the old build trees
and the old identity runs' outputs deleted, their logs kept).

**Builds.** The same configure line as before: `build_base` (Release, make 87 s) and `build_base_debug` (Debug, 38 s);
the fast gate 15 of 15 on both (10.9 s and 32.7 s; `R/scratch/logs/build_base_r4.log`, `build_base{,_debug}.gate.log`).

**Off-state identity, the base against itself, suite 6vxx only** (ruling R4): `simulated_workflow suite=6vxx` at one
thread and `SIMPLE_SEED=7`, twice on `build_base` side by side (`identity/base1`, `identity/base2`; 10:44 and 12:56
minutes); `test=identity` of `base1` against `base2`: 186 files the same and 0 differing
(`identity/cmp_base1_base2/6vxx/identity_report.txt`). **The exclusion list is empty.** `base1` stays for Phases 1
and 4; `base2`'s outputs were deleted, its logs kept.

**The truth sets and the uncentred test half** are reused (no solve3D in them, ruling R4): their checksums hold
(`truth_sets.md5`, 13 files; `sets/S_unc/MD5SUMS`, 5 files). Their references were recomputed on the new base's
dense table (`tools/truth_ref.sh`, `tools/truth_ref_unc.sh`, the V0 variant's rows, `truth_ref/`): the reference lies
within 6° of the truth for 99.13% of S (median 2.70°), as on the old base; the dense search from the uncentred half's
raw picks (D0) agrees with it for 64.2% (72.5% within 6° of the truth). 72 of S's 24,000 modes differ from the old
base's reference, the same order as the dense table's irreproducibility between two runs found before (74).

**The teachers**, rebuilt on the new base exactly as before (`tools/teacher.sh`, `tools/freeze_teacher.sh`,
`teachers/<set>_t<t>/`, `MD5SUMS`; the references and D0 frozen beside them by the first V0 run of each teacher,
`REFERENCE_MD5SUMS`):

| set | teacher seed | solve3D (s) | final FSC 0.5 / 0.143 (Å) | validation wall (min:s) | load at start / end |
| --- | --- | --- | --- | --- | --- |
| bgal | 1 | 456 | 7.10 / 4.13 | 9:00.06 | 19.50 / 14.97 |
| bgal | 2 | 497 | 7.59 / 4.18 | 9:38.89 | 14.97 / 9.85 |
| bgal | 3 | 492 | 7.10 / 4.13 | 9:33.65 | 9.85 / 10.10 |
| strep | 1 | 783 | 3.61 / 3.19 | 14:10.40 | 10.10 / 19.23 |
| strep | 2 | 787 | 3.61 / 3.12 | 14:14.70 | 19.23 / 23.80 |
| strep | 3 | 791 | 3.61 / 3.12 | 14:17.18 | 23.80 / 19.93 |
| trpv1 | 1 | 1277 | 8.89 / 7.59 | 24:44.01 | 19.93 / 16.37 |
| trpv1 | 2 | 1253 | 9.15 / 7.78 | 24:13.63 | 16.37 / 13.64 |
| trpv1 | 3 | 1226 | 9.15 / 7.59 | 23:51.00 | 14.59 / 16.49 |

The teachers ran one at a time with nothing else of the run's on the machine (the one-minute loads are the runs'
own). Their solve3D took 2 to 15% less wall time than the old base's (bgal's third, contended there, aside); their final
FSC values lie within the spread of the old base's teachers (a 24-thread solve3D does not repeat itself). **The first run of trpv1's third teacher crashed** in the base's solve3D, half a second
into its first refine3D stage, on a corrupted heap (`malloc(): unaligned tcache chunk detected`, SIGABRT; the log is
`teachers/trpv1_t3_crash1/`). The same command at the same seed then ran to the end (the table's row) on the same tree:
the fault is not deterministic, it lies in base code outside this plan, and it is left for the maintainer, like the
non-deterministic crash of the class_raw run. A bounds-checked Debug build of the same case (`build_base_debug`, seed 3, 24 threads; `debug_crash/trpv1_t3/validate.log`) passed the point of the crash and entered solve3D's second stage with no runtime error before it was stopped.

**How the arms were run**: as in the superseded entry. V0 is the same patch re-applied to the new export (it applied
cleanly with the same changes; `variants/V0/V0.patch`, Release build make 116 s, fast gate 15 of 15). C0 and R0 share
one model per run, `predict_poses` with centring off.

**Real sets, C0 → R0** (mean over the nine teacher-and-seed pairs of each set; the paired difference R0 − C0 and its
standard error in brackets; `reads_real.tsv`, `tables_p0.md`):

| read | bgal | strep | trpv1 |
| --- | --- | --- | --- |
| seeded (%) | 4.1 → 0.8 (-3.3, 0.5) | 0.0 → 0.0 (+0.0, 0.0) | 0.0 → 0.0 (+0.0, 0.0) |
| cone capped at 64 (%) | 38.0 → 85.1 (+47.1, 0.6) | 100.0 → 100.0 (+0.0, 0.0) | 92.4 → 98.3 (+5.9, 0.5) |
| cone holds the reference (%) | 91.1 → 49.3 (-41.8, 0.3) | 70.4 → 20.5 (-49.9, 0.5) | 90.4 → 65.1 (-25.3, 0.6) |
| mode within 6° (%) | 69.9 → 20.6 (-49.3, 0.7) | 16.6 → 3.6 (-13.0, 0.2) | 41.2 → 19.5 (-21.7, 0.3) |
| mode within 12° (%) | 86.2 → 37.7 (-48.5, 0.4) | 44.5 → 11.9 (-32.7, 0.4) | 62.7 → 36.2 (-26.5, 0.6) |
| committed within 6° (%) | 68.1 → 19.2 (-48.9, 0.7) | 13.5 → 2.3 (-11.2, 0.3) | 35.2 → 14.2 (-21.0, 0.3) |
| committed within 12° (%) | 85.3 → 34.5 (-50.7, 0.5) | 37.8 → 8.2 (-29.6, 0.5) | 50.6 → 24.7 (-25.9, 0.6) |
| committed, median (°) | 3.6 → 21.0 (+17.4, 0.4) | 16.2 → 43.8 (+27.6, 0.3) | 11.7 → 39.2 (+27.5, 0.5) |
| in-plane error, median (°) | 1.49 → 2.08 (+0.58, 0.06) | 3.78 → 59.15 (+55.37, 11.60) | 2.27 → 3.50 (+1.23, 0.18) |
| shift to the reference, median (Å) | 0.4 → 8.0 (+7.6, 0.2) | 1.6 → 12.5 (+10.8, 0.1) | 1.6 → 5.7 (+4.0, 0.1) |
| shift to the reference, q90 (Å) | 6.5 → 22.3 (+15.8, 0.1) | 6.6 → 37.2 (+30.6, 0.3) | 6.7 → 14.8 (+8.1, 0.2) |
| unseen map, FSC 0.143 (Å) | 5.12 → 8.56 (+3.44, 0.16) | 7.02 → 13.58 (+6.56, 0.22) | 7.74 → 9.36 (+1.62, 0.51) |
| unseen map, FSC 0.5 (Å) | 7.86 → 13.22 (+5.36, 0.25) | 9.51 → 13.72 (+4.21, 0.11) | 9.01 → 14.15 (+5.14, 0.06) |

**D0, the dense search on the raw picks**, against the reference (one pass per teacher, the mean of three):

| set | within 6° (%) | within 12° (%) | median (°) | shift to the reference, median / q90 (Å) | C0's committed within 6° (%) |
| --- | --- | --- | --- | --- | --- |
| bgal | 51.9 | 70.9 | 4.5 | 0.4 / 9.0 | 68.1 |
| strep | 5.6 | 13.5 | 40.0 | 5.2 / 35.0 | 13.5 |
| trpv1 | 54.9 | 61.5 | 4.6 | 0.8 / 9.1 | 35.2 |
| truth sets (S's images from the uncentred half's shifts; within 6° of the truth 72.5, median 3.6°) | 64.2 | 84.4 | 0.0 | 0.4 / 7.3 | - |

**Truth sets, C0 → R0** (eight training seeds each; `reads_truth.tsv`). The truth sets hold no solve3D: their reads
equal the old base's to the reported decimals but for one tenth (the reference moved in 72 of 24,000 modes):

| read | T200 | T1000 | T5000 | T |
| --- | --- | --- | --- | --- |
| seeded (%) | 21.8 → 5.7 (-16.1, 1.5) | 25.7 → 6.3 (-19.3, 1.0) | 45.3 → 11.1 (-34.2, 1.2) | 53.9 → 10.5 (-43.4, 1.6) |
| cone capped at 64 (%) | 21.2 → 52.4 (+31.1, 1.2) | 18.6 → 55.3 (+36.7, 1.9) | 6.9 → 45.7 (+38.8, 1.5) | 3.2 → 57.4 (+54.2, 1.1) |
| cone holds the reference (%) | 97.0 → 59.9 (-37.1, 0.5) | 98.0 → 63.2 (-34.9, 0.4) | 98.9 → 65.2 (-33.7, 0.2) | 99.9 → 70.3 (-29.5, 0.1) |
| cone holds the truth (%) | 96.5 → 59.3 (-37.3, 0.5) | 97.7 → 62.6 (-35.1, 0.4) | 98.4 → 64.5 (-33.8, 0.2) | 99.5 → 69.8 (-29.7, 0.1) |
| mode within 6° of the truth (%) | 74.2 → 23.5 (-50.7, 0.7) | 78.7 → 26.3 (-52.4, 0.5) | 87.6 → 31.2 (-56.4, 0.4) | 93.4 → 36.9 (-56.5, 0.3) |
| committed within 6° of the reference (%) | 50.5 → 12.6 (-37.9, 0.8) | 53.3 → 14.2 (-39.0, 0.7) | 61.4 → 17.5 (-43.9, 0.6) | 68.5 → 22.1 (-46.4, 0.6) |
| committed within 12° of the reference (%) | 82.9 → 37.8 (-45.1, 0.4) | 88.1 → 41.5 (-46.6, 0.4) | 95.3 → 47.4 (-47.9, 0.2) | 99.2 → 51.8 (-47.5, 0.3) |
| committed within 6° of the truth (%) | 74.1 → 23.2 (-50.8, 0.7) | 78.7 → 26.0 (-52.6, 0.5) | 87.5 → 31.0 (-56.6, 0.4) | 93.4 → 36.6 (-56.7, 0.3) |
| committed against the truth, median (°) | 3.8 → 25.0 (+21.2, 0.8) | 3.6 → 18.0 (+14.3, 0.6) | 3.3 → 12.4 (+9.1, 0.3) | 3.1 → 10.1 (+7.0, 0.1) |
| shift to the truth, median (Å) | 1.2 → 6.3 (+5.1, 0.1) | 1.1 → 5.7 (+4.5, 0.1) | 1.0 → 4.8 (+3.8, 0.1) | 1.0 → 4.3 (+3.4, 0.0) |
| shift to the reference, q90 (Å) | 3.9 → 20.7 (+16.8, 0.1) | 3.2 → 19.8 (+16.6, 0.1) | 2.6 → 18.4 (+15.7, 0.1) | 2.3 → 17.2 (+14.9, 0.0) |
| unseen map, FSC 0.143 (Å) | 6.06 → 8.54 (+2.48, 0.10) | 6.04 → 7.93 (+1.90, 0.14) | 6.00 → 7.53 (+1.53, 0.15) | 6.00 → 7.15 (+1.15, 0.08) |
| unseen map, FSC 0.5 (Å) | 8.47 → 13.30 (+4.83, 0.42) | 8.16 → 11.73 (+3.57, 0.10) | 7.61 → 11.46 (+3.85, 0.08) | 7.47 → 11.02 (+3.55, 0.13) |

The reference itself against the truth (S): 99.1% within 6°, median 2.70°, its shift 0.93 Å from the truth's at the
median (1.91 at the 90th percentile).

**Costs per particle** (ms, wall time at 24 threads unless marked; as in the superseded entry):

| set | prediction, C0 | prediction, R0 | forward pass (thread-ms) | registration (thread-ms) | reference: fill | reference: all | D0: fill | D0: all |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| bgal | 0.393 | 0.405 | 3.832 | 0.384 | 3.463 | 6.195 | 3.960 | 6.639 |
| strep | 0.076 | 0.076 | 0.596 | 0.119 | 2.180 | 2.591 | 2.313 | 2.694 |
| trpv1 | 0.193 | 0.196 | 3.115 | 0.441 | 3.173 | 3.558 | 3.380 | 3.748 |

**What it shows.**
- **The new base changes nothing that matters here.** Every real-set read lies within about one point (or 0.1 Å; trpv1's
  raw-pick map 0.8 Å) of the old base's, the truth sets' within a few tenths.
- **Dropping the 2D shifts alone costs every read on every set, far beyond its noise.** The cone's hold on the
  reference falls by 25 to 50 points on the real sets (to 21 to 65%) and by 30 to 37 on the truth sets (to 60 to 70%);
  the committed direction within 6° by 11 to 57 points; the unseen map coarsens by 1.2 to 6.6 Å at FSC 0.143 and by
  3.6 to 5.4 Å at 0.5. The shares' paired errors are 0.1 to 1.9 points: C0 and R0 share their model.
- **On raw picks the cone still holds the reference for most particles where the student is informative**: 65% on
  trpv1, 60 to 70% on the truth sets, 49% on bgal, against a mode within 6° of the reference for 21% (bgal) and 20%
  (trpv1) and within 6° of the truth for 24 to 37% (the truth sets). Streptavidin is the exception (21%): even C0's
  student is weak there (every cone capped, the mode within 6° for 17%).
- **The raw picks' cones are wide**: capped at 64 directions for 46 to 57% of the truth sets' particles and 85 to 100%
  of the real sets' (C0: 3 to 38%, except trpv1 and streptavidin, already above 90%).
- **Today's registration does not centre a raw pick**: its shift lies a median 6 to 13 Å from the reference's on the
  real sets (C0: 0.4 to 1.6 Å), 4.3 to 6.3 Å from the truth's on the truth sets, 15 to 37 Å at the 90th percentile.
- **The dense search on raw picks (D0) is no remedy**: 52 to 55% within 6° of the reference on bgal and trpv1, 6% on
  streptavidin, 64% on the truth sets' images; on trpv1 it beats C0's committed direction (55% against 35%), where the
  base student seeds nobody.
- **Costs**: the prediction takes 0.08 to 0.41 ms per particle, R0 as C0; the dense table's fill 2.2 to 3.5 ms.

**Decision 9.** (if R0's cone holds the reference's direction for fewer than half of the particles on every set, skip SK). R0's
cone holds it for 49.3% on bgal and 20.5% on streptavidin, but for 65.1% on trpv1 and 59.9, 63.2, 65.2 and 70.3% on T200,
T1000, T5000 and T: not below half on every set, so the early stop does not apply; SK is measured in Phase 2.

**Exit.** Each item of Phase 0's exit on the new base and its evidence:
- The base pinned (`fd9af0add`, ruling R4), Release and Debug builds, the fast gate on both: above.
- The truth sets' and the uncentred half's checksums: verified (reused, ruling R4).
- The off-state identity of the base against itself, suite 6vxx (ruling R4): empty exclusion list.
- The nine teachers rebuilt and frozen, each solve3D's wall time and final FSC recorded (one crashed run repeated);
  their references frozen beside them; the old base's teachers kept as a record in `old_1ea0f8340/teachers/`.
- C0 and R0 on every set at their seeds: 27 real-set pairs and 32 truth-set runs (`scratch/logs/chain_r4.log`,
  `chain_r4b.log`).
- The reads of decision 7 that exist without the alignment, the cone's hold on the reference included: the tables
  above, frozen in `reads_real.tsv`, `reads_truth.tsv`, `tables_p0.md`.
- Decision 9: not triggered.

Files changed in the repository: this plan (ruling R4 recorded in the status line, sections 1, 2 and 6; the old base's
entries marked superseded; this entry) and its row in `doc/implementation_notes/README.md`, carried from the old base.

### Phase 1 on the new base, the physical teacher, the reference and the alignment (2026-11-17, the Dell, done; ruling R4)

Phase 1's work on the old base (its entry above) re-applied onto `fd9af0add`: the tracked files' changes with
`patch` (no conflict; the removal of the 3D in-plane code touched none of the pose student's lines), the two new files
copied; every source file equals its final version on the old base. Phase 0's frozen references copied into the nine
new teachers' solve3D directories. The gates, the plants, the variant and the smoke runs were all rerun on the new base.

**What the tree now does.**
- **The alignment after prediction** (`simple_pose_student_align.f90`, new). Given the teacher's map (`vol1`),
  `predict_poses` treats every particle as a raw pick, not yet centred, and aligns it at its picked position against
  the map's reprojections at the student's grid and band: every in-plane rotation of the polar grid and shifts over a
  30 Å disk (three widths of a 10 Å prior) about the stored shift, the classic searches' polar shift map
  (`shiftmap_srch`, unchanged, as the class student's alignment uses it), at the physical temperature. Two lists are
  aligned on the same prediction, each timed on its own: S1, the mode's mirror pair (two references); SK, every
  direction of the cone with its mirror partner (up to 128 references). SK's best alignment is committed to the
  particle's record (direction, in-plane angle, stored shift plus the found increment, and the dynamic centring's state
  for a following refine3D); S1's pose and the readout's are in the prediction table beside it, the walls in
  `<fname>_times.txt`. Phase 2 decides between S1 and SK and removes the loser. Today's shift registration (a 5 by 5
  grid over 25 Å at the readout's pose) is gone, with its two constants.
- **Centring off in `predict_poses`** (Phase 0's finding): the commander sets `center=no` as learn does; with the map it
  also sets the alignment's search keys (`refine=prob`, `trs`, `shift_sigma` 10 Å, `prob_beta` 1).
- **The physical teacher P3** behind one named constant, `POSE_STUDENT_PHYSICAL_TEACHER`, off in the tree (today's
  teacher A, which Phase 2 runs); a variant build flips it for Phase 3 (`variants/P3/P3.patch`). Under P3 learn's
  commander sets `prob_beta` to 1, the class student's pattern. Every row learn reads is audited against the physical
  posterior recomputed from its losses at the table's own temperature (`physical_row_audit`): the training's targets
  must equal `exp[-T (L - L_min)]` normalised, within the rounding bound of the single-precision losses (twice T times
  their spacing). A physical teacher stops on a row read at another temperature or beyond the bound; today's teacher
  logs its departure. The rows header records the audit.
- **The reference** (decision 2, as Phase 0 recorded it): `run_pose_posterior_mode`, learn's rows with no student at
  `prob_beta` 1 over the 30 Å disk; each row's mode with its in-plane angle and refined shift goes to the particle
  field and to `<base>_poses.txt`.
- **`validate_pose_student`** rewritten. It takes a frozen teacher through the existing key `dir_refine` (a solve3D
  directory of an earlier validation; the sample must match), and the reference and the dense search on the raw picks
  (D0) from `pose_reference.txt` and `pose_dense_raw.txt` in that directory when present (Phase 0's frozen ones were
  copied there), else it makes them. The unseen particles are predicted at their picked positions only. The report
  gives, per arm (the readout RO, S1, SK, D0), the reads of decision 7 against the reference with a `READS` line each,
  the per-particle details, the costs per particle (prediction, each alignment, the dense search) and three maps (SK's,
  S1's and the reference's poses). `validate_pose_student`'s help and its `dir_refine` input are in the UI (a file-table
  row added under ruling R3). `doc/algorithms/pose_student.md` describes the alignment, the physical teacher with its
  audit, and the validation.
- **A defect fixed on the way**: `predict_poses`' distributed master bootstrapped sigma2 over its parts, summing the
  spectra in another order than a shared-memory run; the old registration hid it, the alignment did not (the
  `pose_student` workflow test's master-against-shared-memory checks failed on one particle's S1 shift and on the
  records' shifts). The master now bootstraps in its own process; its table and sigma2 then equal the shared-memory
  run's byte for byte (`scratch/tmp/p1/order_dbg`, one and eight threads also identical).

**Testers and plants.** `test=unit_numerics suite=pose_student_alignment` (new, five tests, 230 assertions on a fixture
of four white references in two mirror pairs, particles planted rotated and shifted 5 px beyond a 2 px window): the
planted reference recovers its rotation and shift over the disk; a list naming only the planted reference's mirror
partner finds it; the cone in either order keeps the planted direction and another pair alone scores lower; a
reference named three times is searched once, and a list beyond the object's capacity aligns to nothing (the caller's
status after the parallel region); a team of four aligns as one thread, bit for bit. `suite=pose_student` gains the
audit's test (rows at the physical temperature within the bound; the same rows read at half or twice it beyond; invalid
entries without mass). Nine plants, each a copy of the tree with one check's violation planted, rebuilt and its suite
run (`tools/plants_p1.py`, `scratch/keep/plants/plants.log`, the patches beside it): the mirror not searched, the first
reference kept, the window instead of the disk, the refinement ignored, a reference searched twice, a thread-dependent
score, the audit's targets tempered, the physical masses unnormalised, a bound blind to the temperature. **All nine are
caught.**

**Gates on the final code** (`tools/chain_p1r4.sh`, `scratch/logs/chain_p1r4.log`).
- Builds: the run's Release tree `repo/build` and Debug tree `repo/build_debug`, configured as the base trees
  (Phase 0's configure line), `make`, `make install`; both clean.
- The fast gate (`scripts/run_fast_gate.sh`): 15 of 15 on both trees (11 s and 39 s). `test=unit_numerics` (32 groups,
  2,423 assertions) and `test=unit_pftc_align2D3D` (15 groups, 6,324 assertions): passed on both trees
  (`scratch/keep/gates_p1/`, `tools/gates_p1.sh`).
- CTest `pose_student` and `class_student` on the Release tree: both passed (41.7 s and 21.6 s;
  `scratch/keep/ctest_p1.log`).
- **The plants**, rebuilt on the new base (`variants/plants`, `scratch/keep/plants/plants.log`): the clean tree passes
  both suites; each of the nine plants builds and fails its suite (3 to 141 failed assertions). All nine are caught.
- **The off-state identity**: `simulated_workflow` suite 6vxx (ruling R4) at one thread and `SIMPLE_SEED=7` on the final
  tree (`identity/tree`, 13:40 minutes, beside the gates) against Phase 0's new `base1` (`test=identity`,
  `identity/cmp_base1_tree/`, `scratch/logs/identity_cmp_r4.log`): 186 files the same, 0 differing.
- **The class student's identity** (`tools/guard_p1.sh`, `scratch/keep/guard_p1/`): `learn_class_student` on bgal (the
  frozen 2D teacher of the student_inplane run) at seed 1, 24 threads, on the base and on the tree: the model headers
  equal and the weights' checksums equal. The same for `learn_pose_student` on bgal under today's teacher (the header's
  new audit keys aside): the training is the base's, bit for bit.
- **The P3 variant** rebuilt on the new base (`variants/P3`, its patch the one constant flipped; fast gate 15 of 15).

**Smoke runs** (one seed each, read only as checks of the machinery; each ran alone, the load its own 24 threads,
`scratch/logs/smoke_p1r4.log`). Every read equals the old base's to within a few tenths of a point:
- Truth set T1000, seed 1 (`tools/truth_tree.sh`, `scratch/keep/truth/smoke_r4/T1000_s1`, reads in `aln_reads.txt`): on
  the uncentred half, SK's committed direction lies within 6° of the reference for 62.8% of S's 24,000 particles and
  within 6° of the truth for 72.1% (S1 12.8% and 24.0%, the readout 12.7% and 23.8%); SK's shift lies 0.39 Å from the
  reference's at the median; the map at SK's poses 6.00 Å at FSC 0.143 against the truth map (S1's 6.98 Å).
- bgal, teacher 1, seed 1 (`tools/vrun.sh`, `scratch/keep/val/smoke_r4/bgal_t1_s1`): the frozen teacher and references
  are read; the prediction equals Phase 0's R0 on the new base (seeded 0.40%, capped 85.3%, the cone holds the
  reference for 51.0%, the mode within 6° for 20.1%); SK's committed direction within 6° of the reference 60.2%, S1
  20.1%, D0 52.4% (C0 at this seed: 67.4%); the unseen map at SK's poses 4.41 Å at FSC 0.143 against 4.35 Å at the
  reference's poses (C0 at this seed: 5.26 Å, S1's 6.94 Å).
- The physical teacher (the P3 variant, learn on bgal at seed 1; `scratch/keep/p3_smoke`): 1,500 rows, none off the
  physical temperature, the largest mass difference 3.7e-8 (0.0002 of the rounding bound); today's teacher's rows
  depart by up to 0.60 in mass (the guard's run).

**For Phase 2.** `tools/vrun.sh <label> <tree> <set> <t> <s>` runs the tree's validation on a frozen teacher;
`tools/truth_tree.sh <label> <tree> <set> <seed>` the truth sets (the uncentred half with the reference's poses, SK's
and S1's maps); `tools/aln_reads_truth.py` their reads. SK costs 4.3 ms per particle on S (24 threads, wall) against
0.07 ms for S1 and 3.0 ms for the dense table's fill (Phase 0); on bgal 5.2 ms against 0.08 ms and 3.9 ms. Its cost is
set by the cone cap, since most raw picks' cones are capped (Phase 0).

**Exit** (section 5's Phase 1 gates): the testers of the alignment (a planted shift recovered over the disk, the
mirror of a cone direction searched, SK keeps the best) and the teacher's audit with their plants; the fast gate,
`unit_numerics` and `unit_pftc_align2D3D` on both trees; the `pose_student` and `class_student` workflow tests; the
off-state identity (6vxx, ruling R4); the class student's identity. All met, as listed above. Files changed:
`src/defs/simple_defs_hhcnn.f90`; `src/main/learn/simple_pose_student.f90` and its tester;
`src/main/strategies/parallelization/simple_strategy_pose_student.f90`;
`simple_pose_student_align.f90` and its tester (new); `src/main/commanders/simple/simple_commanders_pose_student.f90`;
`src/main/commanders/test/simple_commanders_test_class.f90`; `src/main/ui/simple_test/simple_test_ui_class.f90`;
`src/main/ui/simple/simple_ui_solve3D.f90`; `doc/algorithms/pose_student.md`; this plan;
`doc/implementation_notes/README.md`.

### Phase 2, the alignment's read (2026-11-17, the Dell, done)

**What ran.** Today's teacher A on the tree of Phase 1, every unseen particle predicted at its picked position and
aligned twice on the same prediction: over the mode's mirror pair (S1) and over the cone's directions with their
mirrors (SK). The real sets: `validate_pose_student` on each of the nine frozen teachers per set at three training
seeds (`tools/vrun.sh`, `tools/chain_p2.sh`, `scratch/keep/val/A/`), 27 runs. The truth sets: the student trained at
eight seeds per set, then `predict_poses` with the truth map on the uncentred half's 24,000 raw picks
(`tools/truth_tree.sh`, `scratch/keep/truth/A/`), 32 runs. The reads (`tools/reads_tree.py`,
`tools/reads_truth_tree.py`, `tools/collect_p2.sh`) are in `scratch/keep/p2/`: `reads_real_A.tsv`, `reads_truth_A.tsv`, the paired comparisons
`paired_real_S1_SK.md` and `paired_truth_S1_SK.md`, and the tables below (`tables_p2.md`). One real run (bgal, teacher 3,
seed 1) died of a segmentation fault 0.27 s after its start, inside the base's `oris%get_pinds` parallel loop at
learn's first step (a heap fault surfacing, as Phase 0's malloc abort in the base's solve3D); twelve restarts of the
same command did not repeat it and the rerun passed (`val/A/failed_bgal_t3_s1_segv` kept). The runs after 17:00 shared
the machine with another session's jobs (load 30 to 55), so their walls are inflated; costs are reported, not gated.

**Real sets** (nine teacher-and-seed pairs each; means, S1 → SK, then the paired difference and its standard error):

| read | bgal | strep | trpv1 |
| --- | --- | --- | --- |
| committed within 6° of the reference (%), S1 → SK | 20.3 → 59.0 (+38.7, 0.33) | 2.7 → 9.1 (+6.4, 0.11) | 18.1 → 65.5 (+47.4, 0.25) |
| committed within 12° (%), S1 → SK | 37.2 → 75.2 (+38.1, 0.26) | 9.1 → 18.7 (+9.6, 0.16) | 30.8 → 70.5 (+39.7, 0.40) |
| committed against the reference, median (°) | 19.0 → 3.3 (-15.73, 0.45) | 43.7 → 37.3 (-6.34, 0.23) | 32.7 → 0.0 (-32.65, 0.41) |
| shift to the reference, median (Å) | 3.51 → 0.31 (-3.20, 0.09) | 8.58 → 4.72 (-3.86, 0.06) | 3.56 → 0.35 (-3.22, 0.03) |
| shift to the reference, q90 (Å) | 19.7 → 9.7 (-10.00, 0.30) | 35.8 → 34.7 (-1.17, 0.20) | 12.5 → 8.2 (-4.30, 0.07) |
| unseen map, FSC 0.5 (Å) | 10.17 → 7.77 (-2.40, 0.03) | 13.72 → 12.61 (-1.11, 0.14) | 13.53 → 9.01 (-4.52, 0.06) |
| unseen map, FSC 0.143 (Å) | 7.07 → 4.67 (-2.40, 0.10) | 12.26 → 9.46 (-2.79, 0.38) | 8.46 → 7.87 (-0.59, 0.11) |
| alignment, ms per particle (wall, 24 threads) | 0.094 → 6.164 (+6.07, 0.30) | 0.048 → 3.402 (+3.35, 0.01) | 0.076 → 5.409 (+5.33, 0.01) |
| D0 within 6° of the reference (%) | 51.9 | 5.6 | 54.9 |
| unseen map at the reference's poses, FSC 0.143 (Å) | 4.31 | 5.06 | 7.78 |
| prediction, ms per particle | 0.653 | 0.116 | 0.191 |
| dense table fill, ms per particle | 3.960 | 2.313 | 3.380 |

**Truth sets** (eight training seeds each; the reference is S's dense table from the centred start):

| read | T200 | T1000 | T5000 | T |
| --- | --- | --- | --- | --- |
| committed within 6° of the reference (%), S1 → SK | 12.6 → 61.2 (+48.6, 0.46) | 14.2 → 64.1 (+49.9, 0.33) | 17.5 → 66.4 (+48.9, 0.22) | 21.9 → 70.1 (+48.1, 0.17) |
| committed within 6° of the truth (%) | 23.4 → 70.9 (+47.6, 0.43) | 26.1 → 73.2 (+47.1, 0.29) | 31.1 → 75.2 (+44.1, 0.27) | 36.6 → 77.4 (+40.9, 0.20) |
| committed against the truth, median (°) | 23.3 → 3.7 (-19.58, 0.91) | 17.0 → 3.6 (-13.41, 0.57) | 12.1 → 3.5 (-8.64, 0.30) | 10.0 → 3.3 (-6.67, 0.13) |
| shift to the truth, median (Å) | 5.84 → 1.35 (-4.49, 0.11) | 5.00 → 1.31 (-3.70, 0.10) | 3.96 → 1.28 (-2.68, 0.08) | 3.35 → 1.24 (-2.11, 0.03) |
| map against the truth map, FSC 0.5 (Å) | 11.21 → 7.53 (-3.68, 0.05) | 10.76 → 7.53 (-3.23, 0.18) | 10.11 → 7.53 (-2.58, 0.00) | 9.82 → 7.53 (-2.29, 0.03) |
| map against the truth map, FSC 0.143 (Å) | 7.17 → 6.00 (-1.17, 0.07) | 6.75 → 6.00 (-0.75, 0.08) | 6.46 → 6.00 (-0.46, 0.05) | 6.14 → 6.00 (-0.14, 0.04) |
| alignment, ms per particle | 0.074 → 4.184 (+4.11, 0.10) | 0.076 → 4.237 (+4.16, 0.10) | 0.080 → 3.916 (+3.84, 0.16) | 0.077 → 4.093 (+4.02, 0.07) |

**What it shows.**
- **SK agrees with the reference far more often than S1 on every set**: by 39 and 47 points on bgal and trpv1, 6 on
  streptavidin, 48 to 50 on the truth sets, each at least 58 of its standard errors. It also agrees with the truth
  (41 to 48 points more on the truth sets), centres the raw picks (the median shift to the reference 0.3 Å on bgal and
  trpv1, 1.2 to 1.4 Å to the truth on the truth sets, against 3.4 to 8.6 Å for S1), and gives the better unseen map on
  every set (by 0.6 to 2.8 Å at FSC 0.143 on the real sets). On the truth sets SK's map reaches the curve's last shell
  (6.00 Å, Nyquist at 3 Å per pixel) at every seed, so its 0.143 read there is a floor, not a measure.
- **S1 is cheap but rarely right**: the mode is within 6° of the reference for a fifth of the raw picks on bgal and
  trpv1 and fewer on streptavidin (Phase 0), and aligning its pair alone cannot reach a direction the mode missed.
- **The cost**: SK takes 3.4 to 6.2 ms per particle on the real sets and 3.9 to 4.2 ms on the truth sets (24 threads,
  wall), fifty to seventy times S1's, and as much as the dense table's fill (2.3 to 4.0 ms). Its cost is set by the cone
  cap (most raw picks' cones hold 64 directions, Phase 0).
- **For decision 8 (c)** (Phase 4, not decided here): against Phase 0's C0 (the prior 2D shifts), SK's committed
  direction is within 6° of the reference less often on bgal (59.0% against 68.1%) and streptavidin (9.1% against
  13.5%) and more often on trpv1 (65.5% against 35.2%); the unseen map is finer on bgal (4.67 Å against 5.12 Å), coarser
  on streptavidin (9.46 Å against 7.02 Å) and about equal on trpv1 (7.87 Å against 7.74 Å).

**Decision 8 (a): SK.** Its committed-direction agreement is higher than S1's on every set, far beyond twice the paired
standard error; the cheaper arm wins only within 2 SE.

**The losing alignment removed.** `predict_poses` aligns over the cone alone (SK, committed as before): the S1 pass, its
pose and shift in the prediction table (`e1_s1` to `y_s1`), its wall in `<fname>_times.txt` (now `predict_s align_s`)
and its bench step are gone. `validate_pose_student` reports the arms RO, SK and D0 (`COSTS predict sk dense`) and two
maps (SK's and the reference's). `doc/algorithms/pose_student.md` says so (its file-table row now lists Phase 2,
ruling R3). Removing S1 leaves SK's results unchanged bit for bit: bgal, teacher 1, seed 1 on the new tree
(`val/p2_sk/bgal_t1_s1`) gives the A run's reads and the same committed poses and shifts in the prediction table.

**Gates** (Phase 1's again, on the final tree; `tools/chain_p2g.sh`, `scratch/logs/chain_p2g.log`):
- Builds: `repo/build` (Release) and `repo/build_debug` (Debug), `make` and `make install`, no warning.
- The fast gate 15 of 15 on both trees; `test=unit_numerics` (32 groups, 2,423 assertions) and
  `test=unit_pftc_align2D3D` (15 groups, 6,324 assertions) passed on both (`scratch/keep/gates_p1/`).
- CTest `pose_student` and `class_student` on the Release tree: both passed (41.5 s and 19.6 s;
  `scratch/keep/ctest_p1.log`).
- The plants, rebuilt from the final tree (`variants/plants`, `scratch/keep/plants/plants.log`): the clean tree passes
  both suites; all nine plants are caught (3 to 141 failed assertions each).
- **The off-state identity**: `simulated_workflow` 6vxx at one thread and `SIMPLE_SEED=7` on the final tree (14:09
  minutes) against Phase 0's `base1` (`identity/cmp_base1_tree/`, `scratch/logs/identity_cmp_p2.log`): 186 files the
  same, 0 differing. (Phase 1's tree run is kept as `identity/tree_p1`.)
- **The class student's identity** (`tools/guard_p1.sh`): `learn_class_student` on bgal at seed 1 on the base and the
  tree: the headers and the weights' checksums equal; `learn_pose_student` under today's teacher the same.
- The P3 variant refreshed from the final tree (its patch the one constant); its fast gate 15 of 15; learn on bgal at
  seed 1: 1,500 rows audited, none off the physical temperature, the largest mass difference 3.7e-8 (0.0002 of the
  rounding bound) (`scratch/keep/p3_smoke`).

**Exit.** Decision 8 (a) taken (SK); decision 9 did not apply (Phase 0); the losing alignment removed; Phase 1's gates
met on the final tree, as listed above. Files changed in this phase:
`src/main/strategies/parallelization/simple_strategy_pose_student.f90`; `src/main/commanders/simple/simple_commanders_pose_student.f90`;
`src/defs/simple_defs_hhcnn.f90` (a comment);
`doc/algorithms/pose_student.md`; this plan (the file table's row of the algorithm's document, this entry, the status
line) and `doc/implementation_notes/README.md`.

### Phase 3, the teacher's read (2026-11-17, the Dell, done; decided under ruling R5)

**What ran.** The physical teacher P3 (the dense table's rows at the physical temperature T, no tempering factor
`prob_beta`, every row audited against the physical posterior) against today's tempered teacher A, both with the
alignment over the cone (SK, decision 8 (a)). The P3 arm is the variant build of Phase 2's final tree with the one
constant flipped (`variants/P3`, `P3.patch`), run on the same frozen teachers, sets and seeds as Phase 2's A arm, which
is its pair: on the real sets `validate_pose_student` on the nine teacher-and-seed pairs per set (`tools/chain_p3.sh`,
`scratch/keep/val/P3/`, 27 runs), on the truth sets the student trained at eight seeds and `predict_poses` with the truth
map on the uncentred half's raw picks (`scratch/keep/truth/P3/`, 32 runs). Removing S1 in Phase 2 left SK's results
unchanged bit for bit, so Phase 2's A runs pair with these. The reads are in `scratch/keep/p3/` (`reads_real.tsv`,
`reads_truth.tsv`, the SK rows `sk_real.tsv` and `sk_truth.tsv`, the paired tables `paired_real_A_P3.md`,
`paired_truth_A_P3.md`, the tables below `tables_p3.md`; `tools/collect_p3.sh`). The audit on real rows passed in every
P3 run (bgal, teacher 1, seed 1: 1,500 rows, none off the physical temperature, the largest mass difference 3.8e-8,
0.0002 of the rounding bound).

**Ruling R5 arrived during the runs** (recorded in decision 8): the unseen map decides (b) and (c); the agreement reads
are diagnostics. A removal of P3's switch drafted on the agreement reads (P3 lower by 0.5 to 1.6 points beyond 2 SE on
bgal and trpv1) was undone before any build, the four truth runs the ruling's delivery stopped were rerun, and the
teacher was judged on the maps.

**Real sets, the decisive reads** (nine teacher-and-seed pairs; A → P3, the paired difference and its SE; * beyond 2 SE)

| read | bgal | strep | trpv1 |
| --- | --- | --- | --- |
| unseen map against the teacher's map, FSC 0.143 (Å) | 4.67 → 4.63 (-0.04, 0.16) | 9.46 → 8.93 (-0.54, 0.24)* | 7.87 → 7.79 (-0.08, 0.14) |
| the same, FSC 0.5 (Å) | 7.77 → 7.77 (+0.00, 0.00) | 12.61 → 12.26 (-0.35, 0.24) | 9.01 → 9.01 (+0.00, 0.04) |

The yardstick, the full search's map from the same particles (the reference's poses), FSC 0.143: bgal 4.31 Å,
streptavidin 5.06 Å, trpv1 7.78 Å.

**Real sets, diagnostics**

| read | bgal | strep | trpv1 |
| --- | --- | --- | --- |
| cone holds the reference's direction (%) | 49.3 → 47.8 (-1.5, 0.5)* | 20.5 → 19.1 (-1.4, 0.5)* | 65.1 → 64.1 (-1.0, 0.3)* |
| committed direction within 6° of the reference (%) | 59.0 → 57.4 (-1.6, 0.4)* | 9.1 → 9.1 (+0.1, 0.1) | 65.5 → 65.0 (-0.5, 0.2)* |
| committed direction within 12° (%) | 75.2 → 74.0 (-1.3, 0.7) | 18.7 → 19.8 (+1.1, 0.2)* | 70.5 → 70.1 (-0.5, 0.2)* |
| mode within 6° (%) | 20.6 → 19.2 (-1.4, 0.3)* | 3.6 → 3.3 (-0.2, 0.2) | 19.5 → 18.9 (-0.6, 0.2)* |
| shift to the reference, median (Å) | 0.31 → 0.34 (+0.02, 0.01)* | 4.72 → 4.37 (-0.34, 0.07)* | 0.35 → 0.37 (+0.02, 0.01)* |
| cone capped at 64 directions (%) | 85.1 → 87.5 (+2.4, 1.2)* | 100.0 → 100.0 (+0.0, 0.0) | 98.3 → 97.6 (-0.7, 0.2)* |
| alignment, ms per particle | 6.16 → 6.08 (-0.08, 0.22) | 3.40 → 3.40 (-0.00, 0.01) | 5.41 → 5.39 (-0.02, 0.01)* |

**Truth sets, the decisive reads** (eight training seeds; the map of the uncentred half's 24,000 raw picks against the
truth map)

| read | T200 | T1000 | T5000 | T |
| --- | --- | --- | --- | --- |
| FSC 0.143 (Å) | 6.00 → 6.00 (+0.00, 0.00) | 6.00 → 6.00 (+0.00, 0.00) | 6.00 → 6.00 (+0.00, 0.00) | 6.00 → 6.00 (+0.00, 0.00) |
| FSC 0.5 (Å) | 7.53 → 7.55 (+0.02, 0.02) | 7.53 → 7.53 (+0.00, 0.00) | 7.53 → 7.53 (+0.00, 0.00) | 7.53 → 7.53 (+0.00, 0.00) |

**Truth sets, diagnostics**

| read | T200 | T1000 | T5000 | T |
| --- | --- | --- | --- | --- |
| mean FSC over the shells finer than 20 Å | 0.6545 → 0.6533 (-0.0012, 0.0006) | 0.6615 → 0.6612 (-0.0003, 0.0008) | 0.6666 → 0.6653 (-0.0013, 0.0005)* | 0.6737 → 0.6737 (+0.0000, 0.0004) |
| committed direction within 6° of the truth (%) | 70.9 → 70.6 (-0.3, 0.1)* | 73.2 → 73.1 (-0.1, 0.2) | 75.2 → 74.8 (-0.4, 0.2)* | 77.4 → 77.4 (-0.0, 0.1) |
| committed direction within 6° of the reference (%) | 61.2 → 60.9 (-0.4, 0.1)* | 64.1 → 63.9 (-0.1, 0.4) | 66.4 → 65.6 (-0.8, 0.2)* | 70.1 → 69.9 (-0.2, 0.1) |
| cone holds the reference's direction (%) | 61.9 → 61.6 (-0.4, 0.1)* | 64.8 → 64.7 (-0.1, 0.4) | 67.2 → 66.4 (-0.8, 0.2)* | 70.9 → 70.7 (-0.2, 0.1) |
| seeded (%) | 3.7 → 4.1 (+0.4, 0.2) | 4.0 → 3.8 (-0.2, 0.4) | 7.5 → 8.2 (+0.8, 0.3)* | 6.4 → 6.1 (-0.4, 0.5) |
| shift to the truth, median (Å) | 1.35 → 1.35 (+0.00, 0.00)* | 1.31 → 1.31 (+0.00, 0.00) | 1.28 → 1.28 (+0.01, 0.00) | 1.24 → 1.24 (+0.00, 0.00) |

The truth sets' FSC reads sit at the curve's last shells: 0.143 at 6.00 Å (Nyquist at 3 Å per pixel) and 0.5 at
7.53 Å in 63 of 64 runs, so they cannot separate the teachers there. As a diagnostic beside them (not a read of the
ruling) the table gives the mean FSC over the shells finer than 20 Å (`tools/fsc_area.py`): P3's is lower by 0.0013 on
T5000 (beyond 2 SE), by 0.0012 on T200 (2.0 SE) and equal on T1000 and T, a difference of 0.2% of the value.

**What it shows.**
- **On the maps P3 is never worse**: on streptavidin its unseen map is finer by 0.54 Å at FSC 0.143 (beyond 2 SE) and
  0.35 Å at 0.5 (within noise); on bgal and trpv1 the two teachers' maps agree within noise (0.04 and 0.08 Å finer under
  P3); on the truth sets the maps reach the same last shells. trpv1's map under either teacher equals the full search's
  own from the same particles (7.8 Å).
- **The agreement reads dip slightly under P3**: the cone holds the reference's direction 1.0 to 1.5 points less often on
  the real sets and up to 0.8 on the truth sets, and the committed direction agrees within 6° 0.5 to 1.6 points less on
  bgal, trpv1, T200 and T5000, all beyond 2 SE; streptavidin's committed direction is as good and within 12° better
  (+1.1). This is the class student's pattern: the physical teacher's agreement with the full search's single answer
  dips while the maps hold or improve.
- **Costs** are unchanged (the alignment 3.4 to 6.2 ms per particle on the real sets, the prediction 0.1 to 0.7 ms).

**Decision 8 (b), on the maps (ruling R5): P3 adopted.** No set's unseen map is worse beyond twice the paired standard
error; streptavidin's is finer.

**The losing teacher removed.** `POSE_STUDENT_PHYSICAL_TEACHER` is gone: `learn_pose_student` always reads the rows at
the physical temperature (its commander sets `prob_beta` to `POSE_STUDENT_PROB_BETA`, 1, as the class student's learn
does), and the audit always enforces it, stopping on a row read at another temperature or beyond the rounding bound
(the rows header's `teacher_physical` key and the tempered teacher's departure log went with it). The adopted tree
reproduces the P3 variant bit for bit (bgal, teacher 1, seed 1, `val/p3_tree/bgal_t1_s1`: the reads, the prediction
table and the model header equal). `doc/algorithms/pose_student.md` describes the physical teacher and its audit (its
file-table row now lists Phase 3, ruling R3).

**Gates** (Phase 1's again, on the final tree; `tools/chain_p3g.sh`, `scratch/logs/chain_p3g.log`):
- Builds: `repo/build` (Release) and `repo/build_debug` (Debug), `make` and `make install`, no warning.
- The fast gate 15 of 15 on both trees; `test=unit_numerics` (32 groups, 2,423 assertions) and
  `test=unit_pftc_align2D3D` (15 groups, 6,324 assertions) passed on both (`scratch/keep/gates_p1/`).
- CTest `pose_student` (now with the physical teacher) and `class_student` on the Release tree: both passed (36.7 s and
  18.5 s; `scratch/keep/ctest_p1.log`).
- The plants, rebuilt from the final tree (`scratch/keep/plants/plants.log`): the clean tree passes both suites; all
  nine are caught (3 to 141 failed assertions each), the audit's three among them.
- **The off-state identity**: `simulated_workflow` 6vxx at one thread and `SIMPLE_SEED=7` on the final tree (13:00
  minutes) against Phase 0's `base1` (`identity/cmp_base1_tree/`, `scratch/logs/identity_cmp_p3.log`): 186 files the
  same, 0 differing. (Phase 2's tree run is kept as `identity/tree_p2`.)
- **The class student's identity** (`tools/guard_p1.sh`, `scratch/keep/guard_p1/`): `learn_class_student` on bgal at
  seed 1 on the base and the tree: the headers and the weights' checksums equal. The pose student's training now differs
  from the base's, as the adopted teacher intends.
- The audit on real rows with the final tree (`learn_pose_student` on bgal at seed 1, `scratch/keep/p3_smoke_tree`):
  1,500 rows, none off the physical temperature, the largest mass difference 3.7e-8 (0.0002 of the rounding bound).

**Exit.** Decision 8 (b) taken on the maps (ruling R5): P3; the losing teacher removed; Phase 1's gates met on the final
tree, as listed above. Files changed in this phase: `src/defs/simple_defs_hhcnn.f90`;
`src/main/commanders/simple/simple_commanders_pose_student.f90`;
`src/main/strategies/parallelization/simple_strategy_pose_student.f90`; `doc/algorithms/pose_student.md`; this plan
(decision 8's revision by ruling R5, the file table's row of the algorithm's document, this entry, the status line) and
`doc/implementation_notes/README.md`.

### Phase 4, the route's read (2026-11-18, the Dell, stopped at decision 8 (c); superseded by ruling R7 and the entry below)

**The base sync (ruling R6).** Ruling R6 (Hans, 2026-11-17: "We can sync code as well") moves the copy to master's HEAD
between Phases 3 and 4. It reached the rulings a minute after Phase 3's last session had read them, so Phase 3 ended
DONE; Phase 4's first session then stopped before any work, awaiting the sync. The maintainer's session moved the copy
to `68638be28` (eight commits since `fd9af0add`, among them the static link of the Fortran runtime that cures the
intermittent heap fault Phases 0 and 2 met), carrying the run's work as a patch (`state/carry_r6/`). The run
re-applied it with `patch` (no conflict: none of the 13 files changed on master) and copied the two new files; all 13
equal their final versions on the old copy. The frozen teachers and Phases 0 to 3's measurements stay on `fd9af0add`;
the gates, the identity and the three validations below ran on `68638be28`.

**Decision 8 (c) on the maps (ruling R5)**, from Phase 0's yardstick C0 (today's teacher on the particles at their 2D
shifts) and Phase 3's adopted raw-pick route (the physical teacher P3, the particles at their picked positions,
aligned over the cone), both on `fd9af0add`, paired by teacher and seed (nine pairs per real set) and by seed (eight
per truth set): the means C0 → raw picks, then the paired difference and its standard error, * beyond twice it
(`scratch/keep/p4/`: `route_real.tsv`, `route_truth.tsv`, `paired_real_C0_P3SK.md`, `paired_truth_C0_P3SK.md`,
`tables_p4.md`; `tools/collect_p4.py`).

| read (C0 → raw picks, paired difference, SE) | bgal | streptavidin | trpv1 |
| --- | --- | --- | --- |
| unseen map, FSC 0.143 (Å) | 5.12 → 4.63 (-0.49, 0.13)* | 7.02 → 8.93 (+1.91, 0.49)* | 7.74 → 7.79 (+0.04, 0.14) |
| unseen map, FSC 0.5 (Å) | 7.86 → 7.77 (-0.09, 0.07) | 9.51 → 12.26 (+2.75, 0.20)* | 9.01 → 9.01 (+0.00, 0.09) |
| committed within 6° of the reference (%), diagnostic | 68.1 → 57.4 (-10.72, 0.70)* | 13.5 → 9.1 (-4.36, 0.25)* | 35.2 → 65.0 (+29.83, 0.54)* |
| cone holds the reference (%), diagnostic | 91.1 → 47.8 (-43.29, 0.45)* | 70.4 → 19.1 (-51.28, 0.22)* | 90.4 → 64.1 (-26.36, 0.84)* |
| the full search's map from the same particles, FSC 0.143 (Å) | 4.31 | 5.06 | 7.78 |

| read (C0 → raw picks) | T200 | T1000 | T5000 | T |
| --- | --- | --- | --- | --- |
| map against the truth map, FSC 0.143 (Å) | 6.06 → 6.00 (-0.06, 0.02)* | 6.04 → 6.00 (-0.04, 0.02)* | 6.00 → 6.00 (+0.00, 0.00) | 6.00 → 6.00 (+0.00, 0.00) |
| the same, FSC 0.5 (Å) | 8.47 → 7.55 (-0.92, 0.06)* | 8.16 → 7.53 (-0.63, 0.11)* | 7.61 → 7.53 (-0.08, 0.08) | 7.47 → 7.53 (+0.06, 0.03)* |
| committed within 6° of the truth (%), diagnostic | 74.1 → 70.6 (-3.47, 1.07)* | 78.7 → 73.1 (-5.56, 0.79)* | 87.5 → 74.8 (-12.75, 0.65)* | 93.4 → 77.4 (-15.97, 0.40)* |

- **bgal**: the raw-pick route's unseen map is finer than C0's by 0.49 Å at FSC 0.143 (beyond 2 SE), equal at 0.5.
- **trpv1**: equal (7.79 against 7.74 Å; both at the full search's own 7.78 Å), though the committed direction agrees
  with the reference far more often (65.0% against 35.2%).
- **Streptavidin**: coarser by 1.91 Å at FSC 0.143 and 2.75 Å at 0.5, both beyond 2 SE. Its student is weak on raw
  picks (every cone capped at 64 directions, the mode within 6° of the reference for 3%; Phase 0), and the alignment
  over a cone that holds the reference for 19% of the particles cannot recover what C0's 2D shifts give.
- **The truth sets**: finer than C0 at 0.5 on T200 and T1000 (by 0.92 and 0.63 Å), equal on T5000, coarser by 0.06 Å
  on T (beyond 2 SE); at 0.143 every map sits at Nyquist (6.00 Å) or one shell from it.
- **Costs** (reported, not gated; per particle, 24 threads, wall): C0's prediction 0.08 to 0.39 ms
  (Phase 0); the raw-pick route's prediction 0.1 to 0.7 ms plus the alignment over the cone 3.4 to 6.1 ms, about the
  dense table's own fill (2.3 to 4.0 ms; Phase 3).

**Decision 8 (c): not met.** The raw-pick route's unseen map is worse than C0's beyond twice the paired standard error
on streptavidin (and by one shell at FSC 0.5 on T). As the plan says, the run makes no further source change and stops
with the table: the route does not yet replace the prior 2D shifts on every set. The tree keeps what Phases 1 to 3
established (the physical teacher, the alignment over the cone given the map, the validation against the physical
posterior's best pose); `predict_poses` without `vol1` predicts as before at the stored shifts. Whether to keep, gate
or drop the raw-pick route is the maintainer's decision.

**Verification on `68638be28`** (`tools/chain_p4b.sh`, `tools/chain_p3g.sh`; `scratch/logs/chain_p4b.log`,
`chain_p4g.log`):
- Builds: the base export `build_base` (Release, 130 s, its fast gate 15 of 15) and the run's `repo/build` (Release) and
  `repo/build_debug` (Debug), configured as Phase 0's; no warning from the run's files (the link's five or six
  executable-stack warnings come with master's static runtime and appear in the base build too).
- The fast gate 15 of 15 on both trees; `test=unit_numerics` (32 groups, 2,423 assertions) and
  `test=unit_pftc_align2D3D` (15 groups, 6,327 assertions) passed on both (`scratch/keep/gates_p1/`).
- CTest `pose_student` and `class_student`: both passed (40.2 s and 19.7 s; `scratch/keep/ctest_p1.log`).
- The nine plants rebuilt from the tree and all caught (`scratch/keep/plants/plants.log`).
- **The off-state identity**: `simulated_workflow` 6vxx at one thread and `SIMPLE_SEED=7` on the new base build and on
  the tree, side by side (11:51 and 12:30 minutes), `test=identity` (`scratch/logs/identity_cmp_p4.log`): 185 files the
  same, 0 differing (master's workflow writes one file fewer than `fd9af0add`'s 186).
- **The class student's identity** (`tools/guard_p1.sh`): the model header and the weights' checksum equal the new
  base's (the base's own weights differ from `fd9af0add`'s: master's changes move training numerics). The teacher's
  audit on real rows: 1,500 rows, none off the physical temperature, 0.0002 of the rounding bound
  (`scratch/keep/p3_smoke_tree`).
- **`validate_pose_student` on the three real sets** at teacher 1, seed 1, with the final tree, each alone
  (`scratch/keep/val/p4/`, `scratch/logs/val_p4.log`), beside the same pair on `fd9af0add` (Phase 3) and C0 there:

| set | unseen map FSC 0.143, Å: here (Phase 3; C0) | committed within 6° of the reference, % | the full search's map, Å | cost: prediction, alignment, dense fill (ms) | wall |
| --- | --- | --- | --- | --- | --- |
| bgal | 4.53 (4.60; 5.26) | 59.8 (57.6; 67.4) | 4.35 | 0.95, 5.09, 3.93 | 45 s |
| streptavidin | 8.58 (7.22; 8.58) | 8.4 (8.5; 13.1) | 4.16 | 0.12, 3.45, 2.31 | 43 s |
| trpv1 | 7.78 (7.78; 7.98) | 65.1 (65.1; 36.3) | 7.78 | 0.17, 5.40, 3.39 | 5:16 |

  Single runs, read as checks: the route works on the new base; each map lies inside the range of its set's nine
  pairs on `fd9af0add` (bgal 4.41 to 5.02, streptavidin 7.22 to 10.56, trpv1 7.41 to 8.19 Å).

**Exit** (section 5's Phase 4 gates): decision 8 (c) taken on the maps (not met; the stop with the table, as the plan
says); Phase 1's gates on the final tree (above, on `68638be28`); `validate_pose_student` on the three real sets at
seed 1 (above); the report beside this plan. Files changed in this phase, against the synced copy: the 13 files of
Phases 1 to 3 re-applied unchanged (ruling R6); this plan (the status line, this entry);
`doc/implementation_notes/README.md`; the report `pose_student_physical_raw_picks_report.md` (new). No new source
change.

### Phase 4 under ruling R7, the centring policy (2026-11-18, the Dell, done)

**Ruling R7** (Hans, answering the stop above) withdraws decision 3 and sets the centring policy (decision 3 records it
in Hans's words): training and predicted particles centred by the prior solve2D, shift refinement only after the
prediction, over the cone, against the training set's reference. The raw-pick route is not adopted. Decision 8 (c)
judges one new arm, C3SK: Phase 3's physical-teacher students (no retraining) predicted on the unseen particles at
their 2D shifts (C0's input: the real sets' solve2D shifts, the truth sets' centred test half S) and aligned over the
cone against the teacher's map, on `fd9af0add` so that it pairs with C0 and P3SK.

**How it ran.** The real sets through a measurement variant of `fd9af0add`'s final tree (`variants/C3`, `C3.patch`,
nine lines: the unseen particles selected at their 2D shifts, the one-line change the policy makes, and `infile` to
predict with a trained model instead of training, Phase 3's own model of each pair; built 80 s, its fast gate 15 of
15; `tools/vrun_c3.sh`, `scratch/keep/val/C3/`, 27 runs). The truth sets with `fd9af0add`'s tree build itself
(`tools/truth_c3.sh`: Phase 3's model of each set and seed on S with the reference's eulers stored, C0's input;
`scratch/keep/truth/C3/`, 32 runs). Every run passed (`scratch/logs/chain_c3.log`); the logs show no training and the
models' thresholds equal Phase 3's. The tables: `scratch/keep/p4/` (`c3_real.tsv`, `c3_truth.tsv`, `route_real.tsv`,
`route_truth.tsv`, `paired_real_C0_C3SK.md`, `paired_truth_C0_C3SK.md`, `paired_*_P3SK_C3SK.md`, `tables_c3.md`;
`tools/collect_c3.sh`, `tools/collect_p4.py`).

**Real sets, C0 → C3SK** (nine teacher-and-seed pairs; the means, the paired difference and its SE; * beyond 2 SE)

| read | bgal | strep | trpv1 |
| --- | --- | --- | --- |
| unseen map, FSC 0.143 (Å) | 5.12 → 4.32 (-0.80, 0.09)* | 7.02 → 4.90 (-2.12, 0.38)* | 7.74 → 7.77 (+0.02, 0.12) |
| unseen map, FSC 0.5 (Å) | 7.86 → 7.59 (-0.27, 0.07)* | 9.51 → 8.58 (-0.93, 0.11)* | 9.01 → 8.72 (-0.28, 0.09)* |
| committed within 6° of the reference (%), diagnostic | 68.1 → 90.8 (+22.64, 0.59)* | 13.5 → 62.1 (+48.55, 0.20)* | 35.2 → 86.6 (+51.40, 0.25)* |
| cone holds the reference (%), diagnostic | 91.1 → 90.8 (-0.31, 0.21) | 70.4 → 72.2 (+1.76, 0.49)* | 90.4 → 89.4 (-0.97, 0.10)* |
| shift to the reference, median (Å), diagnostic | 0.41 → 0.00 (-0.41, 0.01)* | 1.64 → 0.13 (-1.51, 0.02)* | 1.65 → 0.00 (-1.64, 0.05)* |

The yardstick, the full search's map from the same particles, FSC 0.143: bgal 4.31, streptavidin 5.06, trpv1 7.78 Å.

**Truth sets, C0 → C3SK** (eight seeds; the map against the truth map)

| read | T200 | T1000 | T5000 | T |
| --- | --- | --- | --- | --- |
| FSC 0.143 (Å) | 6.06 → 6.00 (-0.06, 0.02)* | 6.04 → 6.00 (-0.04, 0.02)* | 6.00 → 6.00 (+0.00, 0.00) | 6.00 → 6.00 (+0.00, 0.00) |
| FSC 0.5 (Å) | 8.47 → 6.86 (-1.61, 0.06)* | 8.16 → 6.86 (-1.30, 0.11)* | 7.61 → 6.86 (-0.75, 0.08)* | 7.47 → 6.86 (-0.61, 0.03)* |
| committed within 6° of the truth (%), diagnostic | 74.1 → 97.3 (+23.22, 0.97)* | 78.7 → 98.1 (+19.48, 0.85)* | 87.5 → 98.9 (+11.33, 0.62)* | 93.4 → 99.1 (+5.75, 0.40)* |
| shift to the truth, median (Å), diagnostic | 1.19 → 0.95 (-0.25, 0.01)* | 1.14 → 0.94 (-0.20, 0.01)* | 1.05 → 0.93 (-0.11, 0.01)* | 0.99 → 0.93 (-0.06, 0.01)* |

**The input shifts alone, P3SK (raw picks) → C3SK (2D shifts)**, the same models

| real set | bgal | strep | trpv1 |
| --- | --- | --- | --- |
| unseen map, FSC 0.143 (Å) | 4.63 → 4.32 (-0.31, 0.07)* | 8.93 → 4.90 (-4.02, 0.41)* | 7.79 → 7.77 (-0.02, 0.14) |
| committed within 6° of the reference (%) | 57.4 → 90.8 (+33.36, 0.45)* | 9.1 → 62.1 (+52.91, 0.14)* | 65.0 → 86.6 (+21.57, 0.62)* |
| cone holds the reference (%) | 47.8 → 90.8 (+42.98, 0.46)* | 19.1 → 72.2 (+53.04, 0.51)* | 64.1 → 89.4 (+25.39, 0.81)* |
| alignment, ms per particle | 6.08 → 3.64 (-2.44, 0.14)* | 3.40 → 3.20 (-0.20, 0.02)* | 5.39 → 4.99 (-0.40, 0.01)* |

| truth set | T200 | T1000 | T5000 | T |
| --- | --- | --- | --- | --- |
| FSC 0.5 (Å) | 7.55 → 6.86 (-0.69, 0.02)* | 7.53 → 6.86 (-0.67, 0.00)* | 7.53 → 6.86 (-0.67, 0.00)* | 7.53 → 6.86 (-0.67, 0.00)* |
| mean FSC over the shells finer than 20 Å | 0.6533 → 0.7265 (+0.07, 0.00)* | 0.6612 → 0.7281 (+0.07, 0.00)* | 0.6653 → 0.7294 (+0.06, 0.00)* | 0.6737 → 0.7299 (+0.06, 0.00)* |
| committed within 6° of the truth (%) | 70.6 → 97.3 (+26.70, 0.23)* | 73.1 → 98.1 (+25.04, 0.27)* | 74.8 → 98.9 (+24.08, 0.19)* | 77.4 → 99.1 (+21.72, 0.12)* |
| alignment, ms per particle | 4.16 → 2.54 (-1.62, 0.08)* | 4.33 → 2.44 (-1.90, 0.12)* | 3.77 → 1.37 (-2.40, 0.07)* | 4.11 → 1.17 (-2.94, 0.05)* |

**What it shows.**
- **The policy's map is as good as C0's or finer on every set**: finer at FSC 0.5 everywhere (0.27 to 0.93 Å on the real
  sets, 0.61 to 1.61 Å on the truth sets) and at 0.143 on bgal (0.80 Å) and streptavidin (2.12 Å, finer than the full
  search's own map from the same particles, 5.06 Å); trpv1 equal at 0.143 (7.77 against 7.74 Å; the full search 7.78).
- **The alignment over the cone recovers the full search's answer**: the committed direction agrees with the reference
  for 91% (bgal), 62% (streptavidin) and 87% (trpv1) of the unseen particles, against 68, 14 and 35% for C0, and with the
  truth for 97 to 99% on the truth sets; the committed shift equals the reference's at the median on bgal and trpv1.
- **The input shifts alone** (the same models, raw picks against 2D shifts) account for the difference with Phase 3's
  raw-pick route: the cone holds the reference 25 to 53 points more often from the 2D shifts, and streptavidin's map
  is 4.0 Å finer.
- **Costs** (reported, not gated; per particle, 24 threads): the alignment over the cone 3.2 to 5.0 ms on the real sets
  and 1.2 to 2.5 ms on the truth sets (cheaper than from raw picks, the cones being narrower), beside the full search's
  table fill of 2.3 to 4.0 ms and C0's prediction and registration of 0.08 to 0.39 ms.
- The truth sets' FSC reads are near their ceiling: 0.143 at Nyquist (6.00 Å) for every C3SK run and 0.5 at 6.86 Å for
  every C3SK run; the mean FSC over the shells finer than 20 Å, a diagnostic, is 0.06 to 0.07 above P3SK's.

**Decision 8 (c) under ruling R7: C3SK adopted.** No set's unseen map is worse than C0's beyond twice the paired
standard error.

**The policy in the tree, and nothing else.** `validate_pose_student` predicts its unseen half at the 2D shifts (the
selection's raw-pick switch is gone); the dense search on the raw picks (D0) is gone with its frozen file
(`pose_dense_raw.txt`, `POSE_DENSE_RAW_FNAME`), its report arm and its cost line (the costs line now reads
`COSTS predict sk reference`, the reference pass being the full search the route replaces); the readout's shift read
compares the stored (2D) shift. `predict_poses` is unchanged in behaviour (it aligns about the stored shift, which is
now the 2D shift); its comments, the alignment module's description, the validation's help and `dir_refine`'s help say
so. `doc/algorithms/pose_student.md` states the policy in R7's words. The earlier verification of the raw-pick tree on
`68638be28` moved to `scratch/keep/old_68638be28_rawpick/` (record).

**Verification on `68638be28`** (`tools/chain_p3g.sh`; `scratch/logs/p4r7_build*.log`, `chain_p4r7g.log`):
- Builds: `repo/build` (Release) and `repo/build_debug` (Debug), no warning; the base build of `68638be28` from the
  earlier entry.
- The fast gate 15 of 15 on both trees; `test=unit_numerics` (32 groups, 2,423 assertions) and
  `test=unit_pftc_align2D3D` (15 groups, 6,327 assertions) on both (`scratch/keep/gates_p1/`).
- CTest `pose_student` and `class_student`: passed (43.2 s and 18.2 s; `scratch/keep/ctest_p1.log`).
- The nine plants rebuilt from the tree and all caught (`scratch/keep/plants/plants.log`).
- **The off-state identity**, `simulated_workflow` 6vxx at one thread, `SIMPLE_SEED=7`, the tree (11:44 minutes) against
  the base build's run of the earlier entry (`scratch/logs/identity_cmp_p4r7.log`): 185 files the same, 0 differing.
- **The class student's identity**: the header and the weights' checksum equal the base's (`tools/guard_p1.sh`); the
  teacher's audit on real rows: 1,500 rows, none off the physical temperature, 0.0002 of the rounding bound.
- **`validate_pose_student`** on the three real sets at teacher 1, seed 1 with the final tree, each alone
  (`scratch/keep/val/p4/`, `scratch/logs/val_p4r7.log`):

| set | unseen map FSC 0.143, Å: here (C3SK on `fd9af0add`; C0) | committed within 6° of the reference, % (C3SK; C0) | the full search's map, Å | cost: prediction, alignment, the full search's fill (ms) | wall |
| --- | --- | --- | --- | --- | --- |
| bgal | 4.41 (4.29; 5.26) | 90.7 (90.4; 67.4) | 4.35 | 0.93, 3.51, 3.45 | 44 s |
| streptavidin | 4.57 (4.73; 8.58) | 62.3 (61.1; 13.1) | 4.16 | 0.11, 2.90, 2.18 | 39 s |
| trpv1 | 7.78 (7.98; 7.98) | 86.9 (86.9; 36.3) | 7.78 | 0.17, 5.00, 3.18 | 4:58 |

  Single runs, read as checks: the policy works on the new base. Streptavidin's and trpv1's maps lie inside their
  sets' nine-pair ranges on `fd9af0add` (4.16 to 6.53 and 7.24 to 8.19 Å), bgal's one shell above its range (4.24 to
  4.35 Å); the base's training numerics moved between the two bases (the base's own model weights differ), so one run
  on the new base is not a paired read.

**Exit** (section 5's Phase 4 gates): decision 8 (c) taken on the maps under ruling R7 (C3SK adopted); the docs (the
algorithm document, this plan's decisions 3 and 8); Phase 1's gates on the final tree on `68638be28` (above);
`validate_pose_student` on the three real sets at seed 1 (above); the report beside this plan. Files changed in this
phase, against the synced copy: the 13 files of Phases 1 to 3 re-applied (ruling R6), then
`src/main/commanders/simple/simple_commanders_pose_student.f90`, `src/defs/simple_defs_hhcnn.f90`,
`src/main/strategies/parallelization/simple_strategy_pose_student.f90`,
`src/main/strategies/parallelization/simple_pose_student_align.f90`, `src/main/ui/simple/simple_ui_solve3D.f90`,
`doc/algorithms/pose_student.md` (the policy), this plan, `doc/implementation_notes/README.md` and the report (new).
