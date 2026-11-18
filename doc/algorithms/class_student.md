# The Class Student

## Problem

2D class assignment searches every particle against every class average at every in-plane rotation and over a
shift window. On data sets of millions of particles, and in the stream's pool, that search recurs. The class student
is the pose student's network (`pose_student.md`) trained on a 2D teacher instead of a 3D one. It learns, per data
set and against one frozen set of class averages, to propose a class, its in-plane rotation and a short list of
candidate classes for an unseen particle, at a fraction of a millisecond. It does not discover classes. Finding them
stays with solve2D on a subset, or with the pool's own iterations.

## Model

The model is the pose student's, through the shared core `polar_student`
(`src/main/learn/simple_polar_student.f90`):
- the features: the particle's polar Fourier half ring at the band, whitened by the frozen sigma2 shell weights;
- the network: the Hermitian HHCNN;
- the readout: a learned template per class, correlated over every rotation of the ring. A class's logit is the
  log-sum-exp of its scores over the rotations at the readout's learned scale, so the class posterior is the
  in-plane marginal of one softmax over class and rotation, and the read-out angle is that softmax's mode. The pose
  student keeps the core's other readout, the maximum over the rotations; the model file records the log-sum-exp
  (`readout=lse`) and the class model is at version 2.

Only the classes differ.

- **The classes are the teacher's class averages, all of them.** Rejected classes stay classes, so a prediction into
  one is a rejection signal; the model records each class's state.
- **No reflection.** The 2D search tests no mirrors, so the network scores the map alone (`hhcnn%new` with
  `reflect=.false.`; the pose student passes `.true.` for its mirror pairs).
- **The templates start from the class averages' polar transforms**, as the pose student's start from reprojections.
- **The candidates** are the classes holding 99% of the posterior (`CAND_MASS`), at most `KMAX_CLASS` (16), in
  descending mass. The read-out angle of the predicted class is e3 = 360 - the readout's rotation, the search's
  convention.

## Training

`learn_class_student` fills the dense 2D probability table (`eul_prob_tab2D`, `refine=prob`) of the selected
particles against the teacher's class averages, prepared as refine2D's matcher prepares them, one batch at a time.
`get_posterior_row` reads each row before the batch's assignment releases it.

**The teacher is the physical posterior**, the one a solve2D run ends at. The learn program sets `prob_beta` to 1
(`CLASS_STUDENT_PROB_BETA`), so the table neither draws nor reads at a tempered temperature. Each class's whole
in-plane row is weighed by the 2D in-plane insertion's own quadrature:
- `inpl_quadrature`, at the particle's physical temperature T;
- at the shift the insertion would take the class at: the refined entry's, else the seed's.

It gives ln Z = ln sum_r exp[-T (d_r - d_min)] and d_min, so the class's mass is the in-plane marginal
exp(ln Z - T d_min). The soft targets are those masses normalised over the classes
(`simple_eul_prob_marginal`, `eul_prob_tab2D%set_inpl_marginal`).

**The rotation targets.** For the 8 classes of largest mass (`KMAX_INPL_ROW`) the teacher keeps their in-plane
posteriors. The rotation loss is the cross-entropy of the network's joint softmax over class and rotation against
the joint targets: the class's mass times its in-plane posterior, binned onto the network's rotation grid. Its
weight is 1 per unit of the teacher's mass.

**The audit.** On its first batch the learn program recomputes every entry by the insertion's own route, and stops
(`<fname>_teacher_audit.txt`) unless:
- the shift, temperature, ln Z, d_min and in-plane posterior are exactly the insertion's;
- every row is read at the physical temperature, and the class distribution it trains on equals the normalised
  exp(ln Z - T d_min) within the rounding of its single-precision losses;
- on every row, the rotation targets read back the teacher's in-plane angle within half a rotation step.

Training is the core's, as for the pose student:
- every tenth micrograph is held out for early stopping and calibration;
- each training row is shifted by a keyed draw in a disk of three shift-prior widths;
- one design is trained, width 16 (`CLASS_STUDENT_NDESIGNS` = 1). Raised, the designs are tried cheapest first and
  a wider one is kept only when its resampled recovery gain is above zero (`choose_student_design`, shared with the
  pose student).

Two settings matter:
- **The shift disk is narrow** (`shift_sigma` 2 Å by default). The teacher's particles are centred. At refine2D's
  default of 10 Å the student learned almost nothing.
- **Width 16 alone** (2026-11-14, Hans). With the choice on, width 32 was kept in 2 of 6 trpv1 learns (gains of
  0.016 and 0.018 in held-out classes recovered; losses up to 0.048 in the others) and in none of the bgal ones,
  at 2 to 4 times width 16's training time. On trpv1 it was about three quarters of the learn program's time.

The confidence threshold is set on the held-out rows so that the teacher's class lies among the candidates of 95% of
the particles the student would seed.

## Prediction and the hand-off

`predict_classes` builds each selected particle's features under the model's band and weights and runs the
forward pass. From the posterior it takes the class of largest mass, its read-out angle and the candidates.

**The alignment after prediction.** Given the class averages, every particle is aligned to them at its picked
position (its stored shift; raw picks are not centred first): every in-plane rotation of the polar grid and shifts
over the disk of the validation's reference pass (three widths of the 10 Å prior, `CLASS_DENSE_SHIFT_SIGMA`, under
that prior), at the student's band and the physical temperature. The search is the classic searches' polar shift map
(`shiftmap_srch`): each class scored at the stored shift, a coarse map over the disk, the refinement from its best
cell; a class keeps the refined pose when it beats its stored-shift score, and the best class by score is kept, as
the classic 2D search commits (`simple_class_student_align.f90`). The particle is aligned against every class of its
candidate list (at most `KMAX_CLASS`), the best alignment kept: the network serves as a neighbourhood of classes,
not as the answer. Aligning against the most likely class alone was measured beside it and rejected: on raw picks
the list holds the reference's class for about nine particles in ten, the mode for one in three
(`doc/implementation_notes/completed/class_student_raw_picks.md`). The aligned class, angle and shift are committed;
the particle's centring state (`shsig`, `shincr`) records the move, as a 2D search leaves it.

It writes:
- each particle's class, e3, shift, confidence (`pconf`) and route (`route`: 1 seeded, 2 dense) to `ptcl2D`;
- a text table;
- the wall times of the prediction and of the alignment (`<fname>_times.txt`, a shared-memory run's);
- the candidate artifact `class_candidates`.

A particle is seeded when its candidates fit the cap and its confidence reaches the model's threshold.

**In parts.** With `nparts`, `predict_classes` is the master of that many parts. Given the class averages it
makes the even/odd halves and sigma2 first; the parts only check them. Each part predicts its range of rows and
writes their records, its table and a candidates part file. The master merges the records, joins the tables in
part order, publishes one artifact from the part files and writes the project once. The result is the
shared-memory run's, byte for byte (the `class_student` workflow test checks it at two parts). Training stays in
one process: its cost is set by the teacher's sample, not by the data set.

**The candidate artifact** (`src/main/project/simple_class_candidates.f90`) is a side file registered in the
project's `out` segment and tied to the particle rows by `ptcl_layout_digest`. It is refused when the rows, the
project or the class count differ, and the consumer then runs as without it. The consumer is the sparse 2D table
(`refine=prob_neigh`): a seeded particle scores its candidates alone and searches the shift window, not the disk.

**No centring.** Unseen particles enter at their picked positions, their shifts as imported; the alignment after
prediction centres them. `center_ptcls` (`center_ptcls.md`) stays a program of its own, no longer on this route
(`doc/implementation_notes/completed/class_student_raw_picks.md`).

## What the reads found

Against reference matching (one `refine2D` pass of the new particles against the teacher's averages) on bgal and
trpv1, see `doc/implementation_notes/completed/class_student_vs_reference_matching.md`: equal class averages, and
per-particle costs that small sets cannot settle for data sets of millions.

On the Dell, Release builds, three seeds unless stated (`doc/implementation_notes/completed/class_student_2D.md`):

- **Held-out rows, centred by the teacher:** the teacher's class among the candidates for 94 to 98% of particles
  on bgal (50 classes) and 82% on trpv1 (100 classes); the in-plane angle within about 1°.
- **Unseen particles** against one dense search of the teacher's averages, scored by its most likely class:
  - centred by `center_ptcls`, among a seeded particle's candidates for 0.90 (trpv1) and 0.73 (bgal);
  - at converged shifts, for 0.96 and 0.89.

  The centring's residual (3 to 4 Å on trpv1, 6 to 8 Å on bgal) against the 2 Å training disk is the gap.
- **Raw picks with the alignment after prediction** (2026-11-16, `completed/class_student_raw_picks.md`, eight
  paired seeds on streptavidin, 30 classes, and bgal and trpv1, 50 each, against the route centred by
  `center_ptcls`):
  - on raw picks the candidate list holds the reference's class for 91 to 97% of particles, the most likely class
    is it for 23 to 34%;
  - the aligned class agrees with the reference's 2.6, 4.4 and 17.4 points more often than the centred route's
    (31, 53 and 72%), and the class averages' FRC 0.5 crossing is 0.75, 2.93 and 4.80 Å finer;
  - the route (prediction and alignment) costs 0.54, 1.16 and 0.85 ms per particle against the full search's 0.44,
    1.27 and 0.92: above it on streptavidin, where a list of at most 16 classes covers half of the 30. Hans adopted
    the route (2026-11-16): small sets do not decide the cost.
- **The physical teacher** (2026-11-15, `completed/student_inplane_posterior.md`, eight paired seeds against the
  base's teacher, the reference being the physical posterior's best class):
  - it shortlists 10.1 (bgal) and 5.6 (trpv1) more particles in 100, with shorter lists;
  - its shortlists hold the reference's class for 94.6% on bgal, at the 95% the threshold is calibrated for,
    against the base's over-cautious 97.7%. On trpv1 they hold it for 98.1% against 99.1%;
  - it seeds more held-out rows (+0.010 on bgal, +0.044 on trpv1) and reads the angle better on bgal (0.84° against
    0.90°).
- **The exact-class read is too strict where the classes repeat views** (bgal's D2). Hans inspected the class
  averages of the student's assignments and judged them good enough (the plan's ruling 15).
- **Cost:**
  - training 205 to 240 s for 10,000 particles and 100 classes;
  - prediction 0.1 to 0.3 ms per particle, and 0.35 to 0.6 ms with reading and preparation;
  - a dense search 1.1 to 1.2 ms at the full box.
- **solve2D at scale** on 35,645 trpv1 particles, two runs per route:
  - route A, solve2D on all: 453 and 468 s;
  - route B, solve2D on 10,000, the student, `center_ptcls` and prediction on the rest, then five seeded
    `prob_neigh` iterations: 450 and 477 s.

  The class-quality model accepts about the same numbers (25 to 26 against 22 to 23 of 100). Seeding barely
  shortens an iteration (about 8 s of each is reading, preparing and summing particles). Route B's per-particle
  cost is about 3.4 ms against 13 ms, so a gain, if any, comes at larger scale; that is not measured.

## Rationale

- **Composition, not inheritance:** the pose and class students each hold a `polar_student`. A Fortran extension in
  another module cannot reach its parent's private state, and one copy of the training, the features, the model
  file and the calibration was the plan's contract.
- **Supervised, per data set:** the teacher is the probabilistic search's own posterior against frozen averages, as
  for the pose student. Unsupervised and generative directions are in
  `doc/implementation_notes/planned/polar_student_future_directions.md`.
- **No new command-line key:** the artifact's presence switches the consumer; without it every path is byte for
  byte the base's.

## Implementation

- `src/main/learn/simple_class_student.f90`: the class (`class_student`), its candidates and read-out, its
  model-file keys; tester `simple_class_student_tester.f90` (suite `class_student`).
- `src/main/learn/simple_polar_student.f90`: the shared core (suite `polar_student`).
- `src/main/strategies/parallelization/simple_strategy_class_student.f90`: `learn_class_student_exec` (and, in its
  reference mode, `run_class_posterior_mode`), `predict_classes_exec`; the shared helpers and the audit's reports in
  `simple_strategy_student_common.f90` (suite `student_strategy`).
- `src/main/strategies/parallelization/simple_class_student_align.f90`: the alignment after prediction
  (`class_student_align`, one per thread); tester `simple_class_student_align_tester.f90` (suite
  `class_student_alignment` in `unit_numerics`).
- `src/main/strategies/search/probabilistic/simple_eul_prob_marginal.f90`: the teacher's in-plane marginal, drawn
  and weighed as the fill and the insertion do (tested in `simple_eul_prob_tab_tester.f90`).
- `src/main/strategies/search/probabilistic/simple_eul_prob_tab2D.f90`: `keep_batch_rows` and `get_posterior_row`
  (the teacher's rows), and the candidates in the sparse fill.
- `src/main/project/simple_class_candidates.f90`: the artifact (suite `class_candidates` in `unit_project`).
- `src/defs/simple_defs_hhcnn.f90`: the constants the students' programs share or set (`KMAX_CLASS`, the training
  designs, the calibration, the validation's settings and file names).
- `src/main/commanders/simple/simple_commanders_class_student.f90`: the commanders of `learn_class_student` and
  `predict_classes` (private, `simple_private_exec`) and the controller `validate_class_student` (`simple_exec`,
  developer visibility).

## Validation

- `validate_class_student`:
  - a random sample of `nsample` particles and solve2D on it (the teacher);
  - the student trained on the result;
  - the unseen particles at their picked positions (shift zero, not yet centred), predicted and aligned,
    then assigned by the full search, the reference: the physical posterior's best class. That is each particle's
    class of largest in-plane marginal mass at the physical temperature, over all classes against the teacher's
    averages, every particle searching the shift disk at the full box from its picked position. It is the hard
    assignment a solve2D run ends at, and the posterior the student learns (`run_class_posterior_mode`, the learn
    program's table in its reference mode; each written class is checked against its row's marginals). The class
    averages of that assignment are the reference's;
  - the same particles predicted again at the reference's shifts (ideal centring, no alignment);
  - the class averages of the student's assignment (every unseen particle at its aligned class, angle and shift).

  Its report (`validate_class_student_report.txt`, also in the log) is written for someone who has not met the
  student. It says what the candidate list and the shortlist are (the candidates of a seeded particle, route 1),
  says what the alignment searches and calls the full search's class the reference. Each line carries a goal and whether
  it is met:
  - particles given a shortlist (as many as possible);
  - shortlists containing the reference's class (the target the model's threshold was calibrated for, 95%), and
    the same predicted at the reference's shifts;
  - the shortlist's length, against all classes;
  - the aligned class's agreement with the reference's;
  - the class averages' FRC 0.5 crossing against the teacher's (interpolated between shells), beside the
    reference's;
  - the cost per particle of the route, prediction and alignment, below the full search's.

  One or two sentences then say what failed and why (the picking error or the student; the cost), and a table
  gives the run time of each phase. The diagnostics go to `validate_class_student_details.txt`: the network's reads
  at the picked positions and at the reference's shifts (agreement on the most probable class, the lists); for
  the alignment the agreement, the shift's distance to the reference's (same class and all, median and 90th
  percentile) and the angle's; the picking error; the class averages' FRC 0.5 (whole shells and interpolated) and
  mean FRC to 15 Å on one common class set; the costs of the prediction and of the alignment.
  Agreement on the most probable class is not a verdict: near-duplicate classes count as misses there (ruling 15
  of the 2D plan).

- The workflow test `class_student` (`simple_test_exec test=class_student`) runs it in miniature on 2,000
  simulated beta-galactosidase particles with planted shifts and gates the programs' contract.
- The stream integration is planned in `doc/implementation_notes/planned/class_student_stream.md`.
