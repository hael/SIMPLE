# The pose student with an amplitude channel

Status, 2026-11-15: planned, not started. Runs after `planned/student_inplane_posterior.md` and shares its Phase 0;
its base is that plan's outcome (the adopted network, or the base if both variants are rejected), pinned when this
plan starts. Supersedes: none. Depends on: the 3D posterior weighting (landed in `86babf86b`; it changes the teacher's
rows), the pose student (validated: yes, `completed/pose_student_training_speed_report.md`).

## 1. Why

The pose student tolerates about one pixel of shift error at 15 Å: on the decentred test half of the medium truth
set, the truth in the cone fell from 84.6% at r = 0 to 26.8% at 3 px and 10.6% at 5 px for a student trained on
centred particles (`completed/learned_pose_state_assignment.md`, section *The centring*). Unseen particles are
predicted at their 2D shifts, whose error is several pixels on real data. The Fourier amplitude of a particle is
exactly invariant to its shift; NUViT (arXiv 2610.32041; `completed/class_student_vs_reference_matching.md`, section 7)
gives its transformer an amplitude channel beside the phase. The prior evidence is against amplitudes alone: a
student on power spectra in place of the complex values reached 0.006 within 6° on T1000 against 0.80 with phases.
The question here is narrower: does an amplitude path beside the complex one make the student more robust to shift
error without losing what the phases carry?

## 2. Contract

- The experiment lives outside the tree: each variant is a patch on a `git archive` export of the base. Only an
  adopted variant enters the tree, in Phase 4.
- The base student is untouched until then. If adopted, the model file's version rises
  (`POSE_STUDENT_VERSION`); a model of another version is refused, as now.
- No new command-line key and no new parameter; the variant is the network's structure.
- The network keeps its symmetries: rotation equivariance on the rotation grid, mirror invariance, conjugation
  symmetry of the complex path, and bit-identical training at any thread count.

## 3. Decisions

1. (2026-11-14, Hans) The experiment: an amplitude channel beside the complex input of the pose student, judged on the
   held-out reads and the shift check, after the 3D posterior weighting lands.
2. **The amplitude path.** Per ring and angle of the network's half ring, the modulus of the complex input (after the
   CTF multiplication and the shell scaling the complex channels get). It runs through its own three blocks of the
   same shape (3 by 3 circular convolutions along angle and ring, kernels symmetric along the angle, a ReLU in place
   of the modulus gate), then its own readout: the angular spectrum of the full ring, even orders only (the modulus
   repeats every half turn), against a per-class amplitude template. The two paths' scores are added at every
   rotation before the maximum. The amplitude term cannot tell a rotation r from r + 180°; the complex term does.
3. **Variants.** A, the base. B, the base plus the amplitude path. B', B with each particle's amplitude channel taken
   from another particle (a fixed permutation): a gain B' also shows comes from capacity, not from the amplitudes.
   C, the amplitude path alone: a wiring check that must reproduce the earlier negative (within 6° at most 0.05 on
   bgal).
4. **Teachers.** Rebuilt on the base in Phase 0 and then frozen with checksums: the medium truth set's subsets
   `T200`, `T1000`, `T5000`, the whole training half `T` (24,000), and bgal
   (`/data2/agent_runs/pose_student/scratch/keep/medium`, `/data2/agent_runs/validate_ps/run5`). Three training seeds
   each.
5. **Reads.**
   - held-out: within 6°, the truth in the cone, the share seeded, held-out recovery (as the training-speed report);
   - the shift check (`predict_poses`, `POSE_SHIFT_CHECK_*`): at 2, 5 and 10 px, the mode kept, the cone holding the
     unshifted mode, the cones' Jaccard index, the posteriors' total variation;
   - **the decentring curve**: prediction on the test half decentred by exactly r = 0, 1, 2, 3, 5, 8 px
     (`S_decentred/`): the truth in the cone and the mode within 12°, for students of the whole half `T`;
   - cost: training and prediction wall time at 24 threads.
6. **Verdict, written now.** B is adopted only if all four hold:
   - (a) every held-out read's mean over the three seeds lies within A's range over its seeds, on every set;
   - (b) on the decentring curve, the truth in the cone at 3 px and at 5 px exceeds A's by more than twice A's
     range over seeds, for every seed's `T` student (margin guessed, re-derived from Phase 0's spread before
     Phase 2);
   - (c) B' does not show the gain of (b);
   - (d) training and prediction take at most 1.3 times A's time (guessed).
   Otherwise B is rejected and the record goes to `rejected/`.
7. **Early stop.** Phase 2 runs B and B' on `T1000` and bgal first. If B's gain at 3 px is below A's range over
   seeds, the experiment stops there as a negative result.
8. **One training schedule.** All variants train with the base's schedule and shift augmentation (`shift_sigma`
   2 Å, as `validate_pose_student` sets it, `STUDENT_SHIFT_SIGMA`). A wider augmentation is a separate question, already measured (`learned_pose_state_assignment.md`).

Defaults taken without a question: the modulus, not its logarithm; the same block widths as the complex path; the
amplitude template initialised from the reprojections' moduli through the untrained path, as the complex templates
are.

## 4. Not in this plan

- The 2D class student: the amplitude path is tried on the pose student first; carrying it to the shared core's
  2D student is a follow-up if it is adopted.
- A learned centre estimator, a wider shift augmentation, a transformer encoder.
- The teacher, the probability tables and refine3D's use of the cones.
- Any file outside the table; the stream.

## 5. Change types and gates

| Phase | Kind | Gate |
| --- | --- | --- |
| 0 | measurement | the base's reads, three seeds, frozen |
| 1 | patch outside the tree | the variants' wiring: symmetries and thread identity in the HHCNN tester's checks, C's negative |
| 2 | experiment | early-stop rule (decision 7) |
| 3 | experiment | the verdict (decision 6) |
| 4 | result-changing for the pose student (if adopted) | fast gate, `unit_numerics` HHCNN and polar-student suites, `pose_student` workflow test, `validate_pose_student` on bgal, the decentring curve again |

## 6. Acceptance

- Spreads and floors come from Phase 0's three seeds per set (measured), recorded before Phase 2 runs. The 1.3 cost
  bound and the "twice the range" margin are guessed and say so.
- Reported, not gated: the maps of `validate_pose_student` on bgal for A and B (unseen particles' map from the
  student's poses and registered shifts), for Hans to inspect.

## 7. Common rules and stop conditions

- No floor or tolerance relaxed; a file outside the table gets a row before the edit; the thread policy (identity
  at one thread, a tester that forces a team for every changed parallel region).
- Stops: the 3D posterior weighting not landed; a teacher that cannot be rebuilt; C learning more than 0.05 within 6°
  on bgal (the path is not what it claims); a symmetry the amplitude path would have to break; a change needed in
  the teacher or the tables.
- Gates no ruling may waive: decision 6 (c), the base's off-state identity, thread identity of the training step.

## 8. Phases

- **Phase 0, baseline.** Taken from `student_inplane_posterior.md`'s Phase 0 when this plan's base is that plan's
  base; otherwise (a variant adopted there) A is the adopted network, retrained at three seeds on the same frozen
  teachers. Every read of decision 5, the decentring curve included; the spreads recorded.
  Exit: the reads table and the spreads frozen in Progress, the teachers' checksums.
- **Phase 1, the variants.** Patches B, B' and C on an export, with the HHCNN tester's symmetry and thread-identity
  checks extended to the amplitude path. C trained on bgal.
  Exit: the checks pass; C within 6° at most 0.05 on bgal.
- **Phase 2, first read.** B and B' on `T1000` and bgal, three seeds, all reads.
  Exit: decision 7 applied, the continue-or-stop recorded.
- **Phase 3, the matrix.** B and B' on every set; the decentring curve for `T`.
  Exit: decision 6 applied, the verdict recorded.
- **Phase 4, adoption or rejection.** Adopted: the patch into the tree, the model version raised, the tests, the
  algorithm page (`doc/algorithms/pose_student.md`). Rejected: the record to `rejected/` with its numbers.
  Exit: the gates of section 5, or the rejected note.

### File table

| File | What | Phases |
| --- | --- | --- |
| `src/main/learn/simple_hhcnn.f90`, `_tester.f90`, `_bench.f90` | the amplitude path, its checks, its benchmark | 1 (export), 4 |
| `src/main/learn/simple_polar_feats.f90`, `_tester.f90` | the amplitude channel of the features, if needed | 1 (export), 4 |
| `src/main/learn/simple_polar_student.f90`, `simple_pose_student.f90` | the model file, the version | 4 |
| `doc/algorithms/pose_student.md` | the network's description | 4 |
| this plan, its report, `doc/implementation_notes/README.md` | record | 0-4 |

## Progress

(one entry per phase, written by the run)
