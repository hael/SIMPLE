# The students learn the in-plane posterior: run report

Run `student_inplane`, 2026-11-15, on the Dell (Oracle Linux 8.10, Intel Core Ultra 9 285, 24 cores, GCC 15.2),
unattended. The plan is `planned/student_inplane_posterior.md`, and its Progress section holds the detail of every
phase. Paths below that start with `R/` are under `/data2/agent_runs/student_inplane`.

The run used three bases, and every number below says which:
- `86babf86b` with the launch's documentation: Phases 0 to 3 and R4's first audit (ruling R1);
- `4282d58a3`: the corrected R and the readings ruling R4 judges (ruling R5);
- `603c44490`: the class student's physical teacher P, the readings ruling R6 judges, and the adoption (rulings R6
  and R7).

## 1. What was learned

**Delivered: P adopted for the class student (ruling R7); the pose student unchanged; M and R rejected.**

SIMPLE's two polar students are small networks trained on the probability search's own posterior:
- the pose student proposes a projection direction with a cone and a confidence (3D);
- the class student proposes a class with a shortlist (2D).

Both use the Harmonic Hermitian convolutional network (HHCNN). The base teacher gives each direction or class the
loss at one in-plane rotation drawn from its posterior, tempered by `prob_beta` (0.3), and the readout takes the
best rotation. The run tried three variants:
- **M**: the in-plane marginal (the sum over every rotation) as the teacher, with a log-sum-exp readout;
- **R**: M plus a rotation loss against the teacher's in-plane posteriors at its 8 heaviest references;
- **P** (ruling R6, class student only): R's teacher and network with no `prob_beta`, so the class distribution
  and the in-plane posteriors are the physical posterior, the one a solve2D run ends at.

**The first verdict (Phase 3, decision 7, `86babf86b`): M and R rejected for both students.** Neither cleared the
bars of decision 7 (b), and each fell below the base's three-seed range on some reads.

**R4's audit found the teacher off the refinement's in-plane weighting, and that was fixed.**
- The learn programs run at `prob_beta` = 0.3, so the teacher weighed every in-plane row at 0.3 T.
- The in-plane posterior insertion weighs the same rows at the physical T, and in 3D adds the shift prior's term.
- On bgal the difference moved a reference's log partition sum by up to 1.94, about a factor of 7 in mass.
- The corrected teacher calls the insertion's own `inpl_quadrature`. On real rows it matches the insertion exactly,
  on every base. The learn program checks this on every run.

**R4 on `4282d58a3`: the corrected R is worse on three reads and better on seven** (eight paired seeds; R − A
and its standard error; worse means beyond twice the standard error in the worse direction):

| set | read | A | R | R − A (SE) |
| --- | --- | --- | --- | --- |
| bgal 2D | shortlist holds refine2D's drawn class (%) | 74.8 | 71.6 | −3.2 (1.1) worse |
| trpv1 | shortlist holds refine2D's drawn class (%) | 91.1 | 90.0 | −1.06 (0.36) worse |
| `T1000` | within 6° | 0.787 | 0.779 | −0.0082 (0.0037) worse |
| bgal 3D | in-plane error against the stored poses (°) | 0.897 | 0.846 | −0.051 (0.008) better |
| bgal 2D | held-out classes recovered / seeded | 0.974 / 0.974 | 0.985 / 0.988 | +0.011 / +0.014 better |
| bgal 2D | held-out in-plane error (°) | 0.901 | 0.841 | −0.060 (0.023) better |
| trpv1 | held-out classes recovered / seeded | 0.860 / 0.862 | 0.881 / 0.885 | +0.021 / +0.023 better |
| trpv1 | particles given a shortlist (%) | 69.4 | 71.9 | +2.5 (0.6) better |

The other 16 reads, among them the truth in the cone on `T1000` and bgal, are within noise. No implementation error
was found, and the run stopped.

**R6 narrowed the question to the class student.** Its teacher was made physical, and `validate_class_student`
judges against the physical posterior's best class: the class of largest in-plane marginal mass at T, replacing
one refine2D draw at `prob_beta` 100.

**R6 on `603c44490`: P against A, eight paired seeds, 24 threads:**

| set | read | A | P | P − A (SE) |
| --- | --- | --- | --- | --- |
| bgal | shortlist holds the reference's class (%) | 97.7 | 94.6 | −3.13 (0.45) worse |
| bgal | the same, perfectly centred (%) | 96.4 | 94.9 | −1.50 (0.58) worse |
| trpv1 | shortlist holds the reference's class (%) | 99.1 | 98.1 | −0.95 (0.11) worse |
| trpv1 | the same, perfectly centred (%) | 99.7 | 99.5 | −0.23 (0.07) worse |
| bgal / trpv1 | particles given a shortlist (%) | 70.6 / 70.0 | 80.7 / 75.6 | +10.1 / +5.6 better |
| bgal / trpv1 | shortlist length (median) | 5.25 / 7.13 | 4.38 / 6.50 | −0.88 / −0.63 better |
| bgal / trpv1 | held-out seeded | 0.974 / 0.862 | 0.984 / 0.906 | +0.010 / +0.044 better |
| trpv1 | held-out classes recovered | 0.860 | 0.900 | +0.040 (0.003) better |
| bgal | held-out in-plane error (°) | 0.901 | 0.842 | −0.058 (0.021) better |

The class FRC on both sets, bgal's classes recovered and trpv1's in-plane error are within noise.

**R7 accepted the four worse reads as the calibration's trade.**
- P's shortlists hold the reference's class at the 95% the threshold is calibrated for (94.6% on bgal), where A is
  over-cautious (97.7%).
- P shortlists 10.1 (bgal) and 5.6 (trpv1) more particles in 100: a correct shortlist for 7.3 and 4.8 more, at 2.8
  and 0.8 more wrong ones.
- On all unseen particles P ranks the reference's class as A does on bgal. On trpv1 it is 0.6 to 1.0 points lower
  from the third candidate on.

**What is now in the repository:**
- **The class student's teacher** is the physical posterior: the learn program sets `prob_beta` to 1, and each
  class's in-plane row is weighed by the 2D insertion's quadrature at the physical temperature.
- **Its network** reads out the log-sum-exp over the rotations and trains with the rotation loss. The class model
  file is at version 2, its readout recorded.
- **The learn program audits its teacher** on its first batch of real rows and stops on any difference.
- **`validate_class_student`** judges against the physical posterior's best class.
- **The pose student, its model files and the 3D teacher code are byte for byte the base's.**

## 2. What changed, phase by phase

| phase | diff | files changed in the repository |
| --- | --- | --- |
| 0, the baseline | `review/phase_0.diff` | the plan (Progress, base pin, decision 6's references, decision 7's bars) |
| 1, the teacher and the readout | `review/phase_1.diff` | the plan (decisions 2 to 4 made precise, three file-table rows, Progress) |
| 2, the first read | `review/phase_2.diff` | the plan (decision 6 revised before any variant read, Progress) |
| 3, the matrix | `review/phase_3.diff` | the plan (Progress) |
| 4, the record, R4 to R7, the adoption | `review/phase_4.diff` | 15 source files (below), `doc/algorithms/class_student.md`, the plan (decision 2's revision, R4 to R7, two file-table rows, status line, Progress), this report; the rejected record written and then deleted |

The phase diffs are under `/home/elmlundho/agent_runs/student_inplane/review/`. The driver re-based the repository
copy onto `4282d58a3` (R5) and then onto `603c44490` (R6), so the raw diff holds only the run's own work.

The adopted sources, all in the class student's path or the code it shares:
- `src/main/strategies/search/probabilistic/simple_eul_prob_marginal.f90` (new): the teacher's in-plane marginal,
  each entry drawn as the fill draws it and weighed as the insertion weighs it, and its audit.
- `simple_eul_prob_tab2D.f90`: the marginal kept beside the batch's rows, the row read at the marginal, the audit.
- `simple_eul_prob_tab_utils.f90`: `inpl_quadrature`'s optional posterior output, the row's marginal record in the
  rows file (version 2), `marginal_dist`.
- `simple_eul_prob_tab_tester.f90`: the marginal's tests.
- The network and the core:
  - `src/main/learn/simple_hhcnn.f90`: the log-sum-exp readout and the rotation loss, an option of the network;
  - `simple_polar_student.f90`: the readout the owner chooses, recorded only for the log-sum-exp; the rotation
    targets;
  - their testers and `simple_hhcnn_bench.f90`.
- `simple_class_student.f90`: the log-sum-exp readout, the rotation targets, the mapping check, model version 2.
- `src/main/strategies/parallelization/simple_strategy_class_student.f90`:
  - the physical teacher and its audits;
  - `run_class_posterior_mode`, the validation's reference, with its self-check under the near-tie rule.
- `simple_strategy_student_common.f90`: the audits' reports.
- `src/main/commanders/simple/simple_commanders_class_student.f90`: `prob_beta` 1 for the learn program; the
  validation's reference, its class averages and its report wording.
- `src/defs/simple_defs_hhcnn.f90`: `KMAX_INPL_ROW`, `CLASS_STUDENT_PROB_BETA`, the reference's file names;
  `CLASS_DENSE_PROB_BETA` removed.
- `src/main/commanders/test/simple_commanders_test_highlevel_class_student.f90`: the workflow test reads the
  reference's renamed timing line.

## 3. Decisions taken under delegation

Each is recorded in the decision's own text or in the Progress entry of its phase.

- **The teacher's marginal** (decision 2, ruling R2): made precise in Phase 1; revised in Phase 4 on R4's audit to
  the insertion's own quadrature, physical temperature and shift-prior term (no shift-prior term in 2D).
- **The rotation targets and loss** (decisions 3 and 4, R2, Phase 1):
  - the targets are the joint over reference and rotation;
  - the log-sum-exp temperature is the readout's own scale;
  - the loss weight is 1.
- **The in-plane reads taken against the stored poses** (decision 6, R2, Phase 2), fixed before any variant was read.
- **The audits made permanent:**
  - the learn program recomputes its first batch by the insertion's route, and checks the class masses at T and
    the targets' mapping on every row;
  - the reference checks every class it writes.
- **The core's readout chosen by its owner in code:** an optional last argument of `polar_student%new`, recorded as
  `readout=lse` only when used. The pose student's code and model files therefore stay the base's, and the weights
  file keeps its version 2 and format.
- **The reference is a mode of the learn routine** (`run_class_posterior_mode`), not a copy of its frame. It writes
  each class as refine2D writes a pose: the stored shift plus the refined increment, and e3 = 360 − rot of the
  in-plane mode.

## 4. Deviations from the plan

- **Files outside the original table**, each given its row before its first edit (rulings R3 and R6):
  - `simple_defs_hhcnn.f90`;
  - the new `simple_eul_prob_marginal.f90`;
  - `simple_commanders_class_student.f90`;
  - `simple_commanders_test_highlevel_class_student.f90` (R7);
  - `code_base_map.md`, whose row proved unnecessary: the file is git-ignored and the build regenerates it.
- **One student adopted, not both** (R6, R7). The pose student keeps A. The 3D teacher code written for R (the 3D
  table's marginal, the pose student's rotation targets and mirror mapping) is left in the run's patches,
  `R/scratch/keep/variants_n/`.
- **Three bases** (R5, R6).
- **16 threads for both arms on `4282d58a3`.** The maintainer's seed profile held about six cores beside the run.
- **The reference's check stopped one validation** (R6's A arm, trpv1, seed 6). It had compared the class taken from
  single-precision losses with an exact double-precision argmax.
  - It now applies the row's own rules: a near-tie within the losses' rounding, or the argmin for a particle without
    a temperature.
  - The re-run passed; its first cause is not pinned.
- **The reference runs in one process**, where the refine2D pass it replaced could run over `nparts`.
- **A planted readout key was not caught at first.** A test now pins that the maximum readout writes no key.
- **The `class_student` workflow test failed on its first Release run.** It reads the report's timing line by the old
  label, which R7's rewording changed. The test's label was updated (a file-table row added), and the gates were
  re-run on both trees from the fixed sources.
- **The phases' adoption list named `validate_pose_student`;** R7's gate list does not, the pose student being
  unchanged.
- **Thread identity:** training is reproducible at a fixed thread count but not between 1 and 24 threads, as on the
  base (R4 item 5).

## 5. Evidence against the baseline

- **The adoption's gates** (ruling R7), on the Release tree `repo/build` and the Debug tree `repo/build_debug`:

  | gate | Release | Debug |
  | --- | --- | --- |
  | fast gate | 15 of 15 | 15 of 15 |
  | `unit_numerics` | 1,900 assertions pass | 1,900 pass |
  | `unit_pftc_align2D3D` | 6,384 pass | 6,384 pass |
  | `pose_student` workflow test (CTest) | passed | passed |
  | `class_student` workflow test (CTest) | passed | passed |

- **`validate_class_student` on bgal at seed 1** (Release):
  - shortlists for 79.4% of the unseen particles;
  - they hold the physical posterior's best class for 95.9%, and 96.1% perfectly centred, both meeting the 95% goal.
- **Earlier builds:** the fast gate passed 15 of 15 on all three bases and every variant build.

- **Off-state identity** (`simulated_workflow` at `SIMPLE_SEED=7`, `test=identity`, no exclusions):

  | base | comparison | 1jxy files | 6vxx files |
  | --- | --- | --- | --- |
  | `86babf86b` | M and R against the base | 223 of 223 | 186 of 186 |
  | `4282d58a3` | corrected M and R against the base | 223 of 223 | 186 of 186 |
  | `603c44490` | R6's A and P arms against the base | 223 of 223 | 186 of 186 |
  | `603c44490` | the repository's Release tree against the base | 223 of 223 | 186 of 186 |

  Each base also matched itself.
- **The pose student is the base's.**
  - Its learn on bgal at seed 1 from the repository's tree gives the base's header and weights byte for byte (md5
    `8f0d176e…`), as on R6's arms.
  - Its sources, the pose strategy and the 3D table are unchanged in the repository.
- **The class student's teacher on real rows** (bgal 2D, seed 1):
  - the in-plane terms equal the insertion's exactly (45,950 entries);
  - the class distribution equals exp(lnz − T d_min) normalised, to 9.1e-6 (0.047 of the bound its single-precision
    losses allow);
  - all rows are read at the physical temperature.
- **Planted violations:**
  - Phase 1: 18, of which 17 were caught; the weak one was replaced and caught;
  - Phase 4: 9 for R4 and 7 for R6, all caught (`R/scratch/keep/plants/`).
- **The reads** are in `R/scratch/keep/`:
  - the tables `summary*.tsv`, `fixed_inpl.tsv`, `topk2d*.tsv`;
  - the paired tests `r5_paired.txt` (R4) and `r6_paired.txt` (R6);
  - the first verdict `verdict7.txt`.

## 6. Left open for the maintainer

- **A crash in base code** (`4282d58a3`). One `validate_class_student` died with SIGSEGV in `sp_project`'s segment
  writer inside refine2D. The same command then completed (`R/scratch/keep/crash_Rn_p5_bgal_s5_validate.log`).
- **Training between thread counts:** not bit-identical between 1 and 24 threads; pre-existing.
- **The pose student's in-plane posterior:** R's 3D code, corrected and audited, is in `R/scratch/keep/variants_n/`,
  should 3D's in-plane insertion come back on by default.
- **Calibration:** P's sharper posterior seeds more particles on unseen data centred by `center_ptcls` than on the
  teacher's centred rows. A calibration on centred-like held-out rows may be a follow-up.
- **Not adopted:** the solve2D-scale and stream uses of the class student were not re-measured with P.

## 7. Proposed skill lines (not applied)

- `simple-main-strategies-probabilistic`:
  - "A route that needs a table entry's in-plane row computes it once and draws with the same keyed variate; the
    entry is then the fill's, bit for bit."
  - "Anything that claims to weigh rows as the in-plane insertion does calls `inpl_quadrature` at the physical
    temperature, not at `prob_beta` × T."
- `simple-architecture` (testing): "A 24-thread job beside other work stalls at its barriers: time nothing beside
  other work. solve2D at more than one thread does not give the same teacher twice."

## 8. Changed and new files

`git status --porcelain` at the end of the run:

```text
 M doc/algorithms/class_student.md
 M doc/implementation_notes/planned/pose_student_amplitude_channel.md
 M doc/implementation_notes/planned/student_inplane_posterior.md
 M src/defs/simple_defs_hhcnn.f90
 M src/main/commanders/simple/simple_commanders_class_student.f90
 M src/main/commanders/test/simple_commanders_test_highlevel_class_student.f90
 M src/main/learn/simple_class_student.f90
 M src/main/learn/simple_hhcnn.f90
 M src/main/learn/simple_hhcnn_bench.f90
 M src/main/learn/simple_hhcnn_tester.f90
 M src/main/learn/simple_polar_student.f90
 M src/main/learn/simple_polar_student_tester.f90
 M src/main/strategies/parallelization/simple_strategy_class_student.f90
 M src/main/strategies/parallelization/simple_strategy_student_common.f90
 M src/main/strategies/search/probabilistic/simple_eul_prob_tab2D.f90
 M src/main/strategies/search/probabilistic/simple_eul_prob_tab_tester.f90
 M src/main/strategies/search/probabilistic/simple_eul_prob_tab_utils.f90
?? doc/implementation_notes/planned/student_inplane_posterior_report.md
?? src/main/strategies/search/probabilistic/simple_eul_prob_marginal.f90
```

## 9. Commit message

```text
class student: learn the physical posterior (in-plane marginal teacher,
log-sum-exp readout, rotation loss); validation against its best class

The class student's teacher becomes the physical posterior a solve2D run
ends at: prob_beta 1, and each class's in-plane row weighed by the 2D
insertion's own quadrature at the physical temperature, its in-plane
marginal the class mass (simple_eul_prob_marginal). The network reads out
the log-sum-exp over the rotations and trains with a rotation loss on the
8 heaviest classes' in-plane posteriors; the class model goes to version 2.
The learn program audits its teacher against the insertion's route on real
rows. validate_class_student now judges against the physical posterior's
best class instead of one refine2D draw. The pose student, its model files
and the 3D teacher are unchanged.

Over eight paired seeds on bgal and trpv1 the class student shortlists
10.1 and 5.6 more particles in 100 with shorter lists and seeds more
held-out rows; its shortlists hold the reference's class at the calibrated
95% (bgal 94.6% against 97.7%). Run student_inplane; plan
planned/student_inplane_posterior.md (M and R rejected).
```
