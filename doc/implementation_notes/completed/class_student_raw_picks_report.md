# Report: the class student on raw picks

Run `class_raw`, 2026-11-16, on the Dell, unattended. Plan: `class_student_raw_picks.md` (beside this report; its
Progress section holds the detail). Base: `1637eeb1c`. Outcome: **adopted** (ruling R4), uncommitted for review.

The class student no longer centres new particles. It predicts on the raw picks, then aligns each particle against
the teacher's class averages of its candidate list (every in-plane rotation, shifts over a 30 Å disk) and commits the
best alignment. Against the old route through `center_ptcls` it agrees more often with the full search's class and
gives sharper class averages on all three data sets, at a cost per particle 2.3 to 17 times lower. It costs less than
the full search on bgal and trpv1 and more on streptavidin. The plan's rule required "less" on every set, so the run
stopped once at that rule. The maintainer waived it for streptavidin, a 30-class set where the full search is
already cheap. The amplitude channel was tested and rejected, and the idea has been removed from the docs.

Terms. The class student is a small network trained on a solve2D run (the teacher). For each new particle it gives
a ranked list of likely classes: its candidate list, the classes holding 99% of its posterior, at most 16. A confident
particle's list is its shortlist. The reference throughout is the full search: each unseen particle's best class
under the physical posterior over all of the teacher's class averages, every rotation and the 30 Å disk. Particle
counts: streptavidin 11,402 (5,000 teacher particles, 30 classes), bgal 5,513 (2,500, 50 classes), trpv1 35,645
(2,500, 50 classes). Every comparison is eight seeds, paired by seed, with the paired standard error (SE) in
brackets. "Beyond 2 SE" is the plan's significance rule.

## 1. What was learned

**Dropping centring alone costs accuracy but keeps the candidate list** (Phase 0). C0 is the base with
`center_ptcls`; R0 is the same validation with the unseen particles at their picked positions and the old local
shift registration:

| R0 − C0 | streptavidin | bgal | trpv1 |
| --- | --- | --- | --- |
| candidate list holds the reference's class (points) | −1.9 (0.5) | −1.8 (0.7) | −0.7 (0.3) |
| most likely class is the reference's (points) | −5.6 (0.8) | −15.7 (0.9) | −19.3 (1.0) |
| route's cost (ms per particle) | 8.92 → 0.16 | 3.25 → 0.65 | 1.98 → 0.40 |

On raw picks the list still holds the reference's class for 89.5 to 96.4% of particles. The plan's early stop
(decision 8, below half) therefore did not apply, and aligning over the list was worth testing. The old local
registration does not centre raw picks: its shift misses the reference's by about the picking error itself (7 to 9 Å
median).

**Align against the candidate list (SK), not the most likely class alone (S1)** (Phase 2, decision 7 (a); one
validation measured both on the same prediction):

| SK − S1 | streptavidin | bgal | trpv1 |
| --- | --- | --- | --- |
| aligned class is the reference's (points) | +7.4 (0.4) | +19.3 (1.0) | +37.9 (0.5) |
| class averages, FRC 0.5 crossing (Å) | −0.43 (0.13) | −2.31 (0.74) | −4.26 (1.14) |
| alignment's cost (ms per particle) | +0.36 | +0.49 | +0.41 |

Where the aligned class is the reference's, the aligned pose equals the reference's to rounding. The alignment
reuses the same shift-map search from the same start.

**The amplitude channel** (Phase 3, decision 7 (b)). B adds an amplitude path to the network beside the complex
one; B' is the same with each particle's amplitudes taken from another particle (the capacity control); C is the
amplitude path alone (the wiring check). B was worse than the tree on reads of every set: streptavidin's aligned
class −4.7 (0.9) points, bgal's class averages 0.65 Å coarser at FRC 0.5, trpv1's share shortlisted −10.0 (3.5).
Its one large gain, bgal's aligned class at +8.0 (0.8) points, B' showed too, at +6.0 (0.7), so it comes from extra
capacity, not from the amplitudes. C recognised 2D classes (0.484 of bgal's held-out classes), as expected for an
amplitude path: shift-invariant and unable to resolve a half turn. The amplitude channel was tested on the class
student 2026-11-16 and rejected; not to be revisited.

**Dropping centring: better quality on every set, cost above the full search on streptavidin only** (Phase 4,
decision 7 (c)). The tree (SK, no amplitude path) against C0:

| SK on raw picks − C0 | streptavidin | bgal | trpv1 |
| --- | --- | --- | --- |
| committed class is the reference's (points) | +2.6 (0.8) | +4.4 (0.9) | +17.4 (0.9) |
| class averages, FRC 0.5 crossing (Å) | −0.75 (0.13) | −2.93 (0.73) | −4.80 (1.27) |
| class averages, mean FRC to 15 Å | +0.043 (0.012) | +0.052 (0.010) | +0.052 (0.021) |
| route's cost (ms per particle), C0's route → the tree's | 8.92 → 0.54 | 3.25 → 1.16 | 1.98 → 0.85 |
| route against the full search (ms per particle) | 0.537 vs 0.442 (above, 8 of 8 runs) | 1.162 vs 1.271 (below) | 0.853 vs 0.923 (below) |

Phase 2's runs give the same picture. Streptavidin fails on cost because the list cap is fixed: with 30 classes the
cap of 16 means SK aligns each particle against about half of all classes (15.98 on average). Ruling R4 waived the
cost condition for this set: the full search there costs 0.44 ms, and small data sets decide nothing.

**The final tree at seed 1** (the adoption's `validate_class_student` run, beside C0's seed 1). One seed only: the
24-thread solve2D teacher is not the same twice, so these numbers do not repeat Phase 3's seed 1.

| seed 1 | streptavidin: route | C0 | bgal: route | C0 | trpv1: route | C0 |
| --- | --- | --- | --- | --- | --- | --- |
| committed class is the reference's (%) | 28.8 | 28.8 | 53.2 | 52.4 | 72.9 | 52.1 |
| class averages, FRC 0.5 (Å); the reference's | 12.50; 13.05 | 13.16 | 13.65; 14.14 | 18.61 | 14.68; 21.13 | 24.23 |
| class averages, mean FRC to 15 Å | 0.897 | 0.884 | 0.824 | 0.778 | 0.788 | 0.696 |
| route's cost (ms per particle); the full search's | 0.484; 0.390 | 8.695 | 1.132; 1.235 | 3.372 | 0.792; 0.842 | 1.688 |

The student's class averages cross FRC 0.5 at a finer resolution than the reference's own on all three sets. The
shortlists hold the reference's class for 96.6% on bgal and 99.3% on trpv1. On streptavidin the student shortlists
9 of 6,402 particles, 8 of which hold the class (the report marks 88.9% as not met; the eight-seed mean is 94.0%).

## 2. What changed, phase by phase

The diffs are in `/home/elmlundho/agent_runs/class_raw/review/phase_N.diff`.

- **Phase 0, the baseline:** the plan only (the Progress entry and notes under decisions 4 and 8); 1 file.
- **Phase 1, the alignment after prediction:** 12 files (+967, −315).
  - A new class `class_student_align` (`src/main/strategies/parallelization/simple_class_student_align.f90`, 103
    lines) aligns one particle against a list of classes. It reuses the classic searches' polar shift-map object
    `shiftmap_srch` unchanged (`src/main/pftc` untouched).
  - `predict_classes` (`simple_strategy_class_student.f90`) aligns every particle at its picked position. Both S1
    and SK were timed separately; SK was committed.
  - `validate_class_student` (`simple_commanders_class_student.f90`) no longer calls `center_ptcls`, measures both
    alignments and reports decision 6's reads.
  - New tester `simple_class_student_align_tester.f90` (suite `class_student_alignment` in `unit_numerics`,
    registered in `simple_commanders_test_class.f90` and the test UI's help).
  - Also changed: the workflow test's reads, `simple_defs_hhcnn.f90`, the algorithm page, the test policy, and a note
    in `polar_student_future_directions.md`.
- **Phase 2, the alignment's read:** 9 files (+248, −224). S1 was removed from the prediction, its table columns,
  the validation, the defs, the workflow test, the tester's wording and the docs.
- **Phase 3, the amplitude channel:** 2 files (+112, −2): the plan, and a note in the future-directions document
  (removed again in Phase 4). The variants B, B' and C stayed outside the tree as patches in
  `/data2/agent_runs/class_raw/scratch/keep/variants/`.
- **Phase 4, the verdict, then the adoption under R4:**
  - the plan: status line, a file-table row, decisions 3 and 7 (c), section 4, Phase 3's entry condensed, and the
    Progress entry;
  - this report;
  - `doc/algorithms/class_student.md`: the raw-pick reads;
  - the README row;
  - the amplitude idea removed from 6 documents: `planned/pose_student_amplitude_channel.md` deleted (148 lines),
    `planned/polar_student_future_directions.md`, `completed/class_student_vs_reference_matching.md` (−6),
    `planned/selection_student.md` (+3, −3), `planned/student_inplane_posterior.md` (+14, −20) and its report
    (+2, −2).

The whole tree against the base:
- 14 modified files and 1 deleted (+485, −488). Of the modified files, 6 are Fortran (+390, −275):
  - `simple_commanders_class_student.f90`, +277 −207;
  - `simple_strategy_class_student.f90`, +89 −48;
  - `simple_commanders_test_highlevel_class_student.f90`, +17 −17;
  - `simple_defs_hhcnn.f90`, +4 −2;
  - `simple_commanders_test_class.f90`, +2;
  - `simple_test_ui_class.f90`, +1 −1.
- 4 new files:
  - the plan, 656 lines;
  - this report;
  - the alignment class, 103 lines;
  - its tester, 271 lines.

## 3. Decisions taken under delegation

Each is recorded in the Progress entry of its phase.

- **Phase 1, the search:** the alignment reuses `shiftmap_srch` and the commit rule of the classic 2D search. A class
  keeps its refined pose only when that pose scores above the class's score at the stored shift, and the class with
  the best score wins.
- **Phase 1, the centring state:** particles are marked not yet centred before the alignment. Afterwards their
  centring state records the move, as a 2D search or `center_ptcls` leaves it.
- **Phase 1, the settings:** when `predict_classes` is given class averages, it sets the prior width, the window
  and the physical temperature of the reference pass.
- **Phase 1, the times file:** the wall times of prediction and alignment go to a separate file
  (`class_student_times.txt`), because the prediction table must stay byte-identical between a run in parts and a
  shared-memory run.
- **Phase 1, the angle:** the recorded angle is wrapped into [0°, 360°), the student's own convention. The workflow
  test's record gate caught a value of 360.0.
- **Phase 2, the class averages' comparison:** only classes populated by both the student's and the reference's
  assignments are compared.
- **Phase 3, the variants' design and controls:**
  - B' takes its partner from a fixed keyed permutation in training and from the batch's next particle in
    prediction;
  - in the gradient check, an entry whose every finite-difference step switches a ReLU is reported, not judged,
    with a guessed bound of at most one in twenty;
  - the tree was re-run beside B and B' so all three arms share the same tools and machine conditions.
- **Phase 4, what R4's removal covers:**
  - also removed, as amplitude proposals: the "front end of Fourier moduli" among the shift-robust alternatives in
    `polar_student_future_directions.md` section 2, and its pitfall sentence;
  - kept, as a recorded command output: the `git status` listing in `student_inplane_posterior_report.md`, the
    one place that still names the deleted file;
  - in `student_inplane_posterior.md`, the launch's list of dirty files now reads "a planned note since deleted",
    so its count of five and its line counts stay true.
- **Phase 4, the gates' reuse:** the tree's source has not changed since Phase 2; the rebuild recompiled only the
  version string. The plants of Phase 1 were not re-planted: the alignment module is byte-identical to the copy
  they ran on, and the tester has only been renamed since.

## 4. Deviations from the plan

- **Decision 7 (c) was not met as written.** The run stopped and wrote a stop report. Ruling R4 waived the cost
  condition for streptavidin and adopted the route. The plan records this in decision 7 (c) and in Progress.
- **C's wiring check was revised under ruling R2 and confirmed by R4.** Decision 3 guessed a bound: at most 0.10 of
  the held-out classes recovered. C reached 0.484. The check became the path's own contract, with thresholds set
  after C's run:
  - shift-invariant reads, measured within 0.2 points;
  - the half-turn ambiguity of the angle, 56% at the opposite angle against B's 2.5%.
- **Files outside the original table**, each given a row before its edit:
  - the alignment class and its tester (R3);
  - the suite table and the test UI's help (R3);
  - the algorithm page and the test policy, which also gained Phase 2 (R3);
  - the five amplitude-removal documents, one of them deleted (R4).
- **An expected effect that did not occur:** the amplitude channel did not help.
- **A base fault seen once:** a segmentation fault in solve2D (`sp_project`'s `write_segment_inside`, called from the
  refine2D strategy) in 1 of about 150 validations (run C0, bgal, seed 1). Re-runs and a bounds-checked Debug run
  passed. It lies outside this run's files and was not pursued.
- **A tooling fix:** the interpolated FRC crossing, copied from the last run, returned the box size when the first
  shell was already below 0.5. It was fixed in Phase 0 and the affected rows were recomputed.

## 5. Evidence against the baseline

- **Off-state identity:** `simulated_workflow` 1jxy and 6vxx at one thread with `SIMPLE_SEED=7`, compared with
  `test=identity`.
  - Base against base: 223 and 186 files the same, so the exclusion lists are empty.
  - The tree against the base after Phase 1, after Phase 2 and on the final tree: each time 223 and 186 files the
    same, 0 differing (`/data2/agent_runs/class_raw/scratch/keep/identity/cmp_base1_tree*`).
- **The pose student:** byte-identical to the base on bgal at seed 1 and 24 threads (model header and weights)
  after Phases 1 and 2, and on variant B. CTest `pose_student` passes on the final tree.
- **Thread identity:**
  - the alignment tester runs a team of four against one thread and matches bit for bit;
  - the prediction with alignment, seeded, at 24 threads twice and at one thread: the same table and candidates;
  - unseeded runs differ in rare shifts, because the reused shift map plans its FFTs by measurement unless a seed is
    set.
- **Gates on the final source:**
  - fast gate 15 of 15 on Release and on Debug, including `unit_numerics` and `unit_pftc_align2D3D`;
  - CTest `class_student`: all ten gates pass, including the run in parts matching shared memory byte for byte;
  - the alignment tester: 172 checks pass, and its four planted violations were all caught.
- **The numbers:** the tables above. Per-run rows are in `/data2/agent_runs/class_raw/scratch/keep/`: `reads.tsv`
  (Phase 0), `reads2.tsv` (Phases 1 and 2) and `reads3.tsv` (Phases 3 and 4), each produced by its script in
  `tools/`.

## 6. Left open for the maintainer

- **A list cap relative to the class count** (R4: a possible follow-up, not done). With 30 classes the fixed cap of
  16 makes SK align each particle against half the classes, which puts the route above the full search on
  streptavidin.
- **To inspect:** the class averages of the final tree at seed 1. Each folder holds the teacher's, the student's and
  the reference's stacks:
  - final tree: `/data2/agent_runs/class_raw/scratch/keep/validate/P4/<set>_s1/keep_cavgs/`;
  - C0, the old centred route: `.../validate/C0/<set>_s1/keep_cavgs/`.
- **Data-dependent:** the costs come from 5,000 to 35,000 particles at 24 threads. As the maintainer ruled, small
  sets settle nothing; the costs of the route and of the full search at a million particles are unmeasured.
- **The base's intermittent solve2D crash** (section 4).
- **The pose student's post-prediction refinement** with its cone used as a neighbourhood, the way SK uses the
  candidate list: noted in `polar_student_future_directions.md` section 6.

**The maintainer's runs**, for each set:

```bash
SIMPLE_SEED=1 simple_exec prg=validate_class_student projfile=<copy of the import project> nsample=<n> ncls=<k> mskdiam=<m> nthr=24
```

The values for `<n> <k> <m>` and the import project for each set:

| set | `nsample ncls mskdiam` | import project |
| --- | --- | --- |
| streptavidin | `5000 30 85` | `/data/10335-streptav/tst/1_import_particles/tst.simple` |
| bgal | `2500 50 180` | `/data2/bgal_addon/1_import_particles/bgal_addon.simple` |
| trpv1 | `2500 50 175` | `/data2/trpv1_c2/1_import_particles/trpv1_c2.simple` |

What to read:
- `validate_class_student_report.txt`: the aligned class's agreement, the class averages' FRC beside the
  reference's, and the cost against the full search;
- `validate_class_student_details.txt`: the shifts, the picking error, and the costs of prediction and alignment;
- `cavgs/cavgs_student.mrc` against `dense/cavgs_reference.mrc` and the teacher's averages.

The workflow test: `simple_test_exec test=class_student`.

**Proposed skill lines** (not applied):
- `simple-main-strategies-parallelization`: "`predict_classes` given class averages aligns every particle at its
  picked position against its candidate list (`class_student_align` over `shiftmap_srch`, unchanged); without
  `SIMPLE_SEED` the shift map plans its transforms by measurement and a rerun may differ in the last digits of a
  shift."
- `simple-modern-fortran`, Testers: "A fixture's spectra assigned to a calculator by pointer
  (`assign_sigma2_noise`) must outlive the fixture routine: keep them at module scope in the tester."
- `simple-main-strategies-search`: "`shiftmap_srch` is the reusable disk search for one particle against a list of
  references (the class student's alignment uses it); its `prep` decides disk or window from `shsig`."

## 7. Files (`git status --porcelain`)

```
 M doc/algorithms/class_student.md
 M doc/implementation_notes/README.md
 M doc/implementation_notes/completed/class_student_vs_reference_matching.md
 M doc/implementation_notes/planned/polar_student_future_directions.md
 D doc/implementation_notes/planned/pose_student_amplitude_channel.md
 M doc/implementation_notes/planned/selection_student.md
 M doc/implementation_notes/planned/student_inplane_posterior.md
 M doc/implementation_notes/planned/student_inplane_posterior_report.md
 M doc/policies/test_environment_policy.md
 M src/defs/simple_defs_hhcnn.f90
 M src/main/commanders/simple/simple_commanders_class_student.f90
 M src/main/commanders/test/simple_commanders_test_class.f90
 M src/main/commanders/test/simple_commanders_test_highlevel_class_student.f90
 M src/main/strategies/parallelization/simple_strategy_class_student.f90
 M src/main/ui/simple_test/simple_test_ui_class.f90
?? doc/implementation_notes/planned/class_student_raw_picks.md
?? doc/implementation_notes/planned/class_student_raw_picks_report.md
?? src/main/strategies/parallelization/simple_class_student_align.f90
?? src/main/strategies/parallelization/simple_class_student_align_tester.f90
```

## 8. Commit message

```
class student on raw picks: align against the candidate list after prediction, no centring

predict_classes aligns each particle at its picked position against the
teacher's class averages of its candidate list: every rotation over the
30 A disk (class_student_align over shiftmap_srch). validate_class_student
drops center_ptcls. On streptavidin, bgal and trpv1 the class agrees with
the full search more often and the class averages are sharper than with
centring. The amplitude channel was tested and rejected, and its plan
removed.
```
