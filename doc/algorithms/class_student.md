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
- the readout: a learned template per class, correlated over every rotation of the ring.

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
`get_posterior_row` reads each row before the batch's assignment releases it, at the temperature the table draws
at. The soft targets are those rows' posteriors over the classes.

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
forward pass. From the posterior it takes the class of largest mass, its read-out angle and the candidates. Given
the class averages, it registers the shift against the predicted class at the predicted angle, with the pose
student's search object (coarse grid, then L-BFGS-B).

It writes:
- each particle's class, e3, shift, confidence (`pconf`) and route (`route`: 1 seeded, 2 dense) to `ptcl2D`;
- a text table;
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

**Unseen particles must be centred first.** Raw picks are centred by `center_ptcls` (`center_ptcls.md`) before
they are labelled. Particles that are already aligned, for example those the stream's sieve hands to the pool, need
no such step.

## What the reads found

On the Dell, Release builds, three seeds unless stated (`doc/implementation_notes/completed/class_student_2D.md`):

- **Held-out rows, centred by the teacher:** the teacher's class among the candidates for 94 to 98% of particles
  on bgal (50 classes) and 82% on trpv1 (100 classes); the in-plane angle within about 1°.
- **Unseen particles** against one dense search of the teacher's averages, scored by its most likely class:
  - centred by `center_ptcls`, among a seeded particle's candidates for 0.90 (trpv1) and 0.73 (bgal);
  - at converged shifts, for 0.96 and 0.89.

  The centring's residual (3 to 4 Å on trpv1, 6 to 8 Å on bgal) against the 2 Å training disk is the gap.
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
- `src/main/strategies/parallelization/simple_strategy_class_student.f90`: `learn_class_student_exec`,
  `predict_classes_exec`; the shared helpers in `simple_strategy_student_common.f90` (suite `student_strategy`).
- `src/main/strategies/search/probabilistic/simple_eul_prob_tab2D.f90`: `keep_batch_rows` and `get_posterior_row`
  (the teacher's rows), and the candidates in the sparse fill.
- `src/main/project/simple_class_candidates.f90`: the artifact (suite `class_candidates` in `unit_project`).
- `src/defs/simple_defs_hhcnn.f90`: the constants the students' programs share or set (`KMAX_CLASS`, the training
  designs, the calibration, the shift registration, the validation's settings and file names).
- `src/main/commanders/simple/simple_commanders_class_student.f90`: the commanders of `learn_class_student` and
  `predict_classes` (private, `simple_private_exec`) and the controller `validate_class_student` (`simple_exec`,
  developer visibility).

## Validation

- `validate_class_student`:
  - a random sample of `nsample` particles and solve2D on it (the teacher);
  - the student trained on the result;
  - `center_ptcls` on every particle;
  - the unseen particles predicted, then assigned by the full search, the reference: each particle's best-scoring
    class over all classes, in-plane rotations and shifts (`refine2D refine=prob` at `prob_beta` 100, one
    iteration against the teacher's averages). It is not `refine=greedy`: it scores as the table the student
    learned from;
  - the same particles predicted again at the full search's shifts (ideal centring);
  - the class averages of the student's assignments.

  Its report (`validate_class_student_report.txt`, also in the log) is written for someone who has not met the
  student. It calls the candidates of a seeded particle (route 1) its *shortlist*, and the class of the full
  search the reference. Each line carries a goal and whether it is met:
  - particles given a shortlist (as many as possible);
  - shortlists containing the full search's class (the target the model's threshold was calibrated for, 95%),
    and the same with the particles perfectly centred (at the full search's shifts);
  - the shortlist's length, against all classes;
  - the class averages' FRC 0.5 resolution against the teacher's, by the student and by the full search;
  - the cost per particle, centring and prediction, below the full search's.

  One or two sentences then say what failed and why (the centring or the student; the cost of centring), and a
  table gives the run time of each phase. The diagnostics (agreement on the most probable class, the angle and
  shift errors, the centring's frame offset, the per-step costs) go to `validate_class_student_details.txt`.
  Agreement on the most probable class is not a verdict: near-duplicate classes count as misses there (ruling 15
  of the 2D plan).

- The workflow test `class_student` (`simple_test_exec test=class_student`) runs it in miniature on 2,000
  simulated beta-galactosidase particles with planted shifts and gates the programs' contract.
- The stream integration is planned in `doc/implementation_notes/planned/class_student_stream.md`.
