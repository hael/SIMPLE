# Training speed of the pose student, and the review's corrections

Status: completed, 2026-11-10 (the run `student_speed`, phases 0 to 6, uncommitted for review; the
report is `pose_student_training_speed_report.md` beside this note). Planned 2026-11-09. Hans set
the priority ("everything we can to speed up training of the student; priority number one") and
asked for an autonomous run; this note is its plan. The same day a code review of the pose student
(`learned_pose_state_assignment_review.md`, at `519161ff3`) found seven open problems; Hans
asked for them to be addressed in the same run. The run is `student_speed` on the Dell (section 9).
Nothing in it is committed; the run hands over one working-tree diff.

The pose student (`doc/algorithms/pose_student.md`) is a small network trained per dataset on the
posterior rows of the probabilistic search. It proposes a direction, a cone, an in-plane angle and a
mirror for particles the search has not seen. Its training is now the slowest part of every use: on
bgal, 1,500 training particles take 180 to 540 s to train against 4 to 15 s to predict 2,756
particles, and on the review's macOS build the registered workflow test `pose_student` timed out at
2,400 s at eight threads while still training. The run first corrects what the review found, then
makes training as fast as it can be made on this machine without changing what the student learns
beyond the run-to-run spread of its own reads.

## 1. What is measured, and the yardsticks

**Speed** is the wall time of `learn_pose_student` on fixed inputs at 24 threads, split as its bench
record `POSE_STUDENT_LEARN` splits it (rows, then per design: batch assembly, gradients, reduction
and optimiser steps, held-out evaluation), and the kernel benchmark `test=pose_student_bench`
(per ring count and width: forward, readout, backward, a training step at one thread and at the
team, prediction throughput, the training step's identity across threads).

**Quality** is read on fixed teachers, so that a change in the student is not drowned by a change
in the teacher (solve3D is unseeded and takes a different stage path from run to run):

- **The truth sets.** The medium assembly of the Ribosembly benchmark, kept by the pose_student
  run under `/data2/agent_runs/pose_student/scratch/keep/medium`: the training half `T` (24,000
  particles) with its map from the truth poses, its subsets `T200`, `T1000`, `T5000`, and the
  centred test half `S` with the truth poses. The student trains on a subset at 15 Å and
  `shift_sigma=2`, predicts `S`, and `pred_stats.awk` reads the folded mode against the truth:
  within 6° and 12°, the truth in the cone, the cone sizes, the confidence. This is the
  generalisation read: `S` was never seen.
- **The bgal fixed teacher.** The training project and the teacher's map of validation run 5
  (`/data2/agent_runs/validate_ps/run5/1_validate_pose_student`): the held-out recovery from the
  model, and the training rows' agreement with the teacher from a prediction on a copy of the
  training project (its rows carry the teacher's poses); the seeded fraction of the unseen particles
  from a prediction on a copy of the run's all-particles project (`predict_poses` writes into the
  project it predicts, so the run's own copy no longer holds pure teacher poses).
- **The spread.** Training draws its batches and augmentations from keys of the run's seed
  (`SIMPLE_SEED`); three seeds per set measure the spread of every read. A change is judged against
  that spread, never against one pair of runs. The review's corrections change what is calibrated
  (the threshold over confidence ties, the held-out split), so the spread is measured twice: at the
  base (phase 0) and after the corrections (phase 1), and the speed phases are judged against the
  corrected one.
- **End to end.** `validate_pose_student` on bgal at the end (phase 6), with the timings table and
  the maps from the student's poses.

The training step must be bit-identical between one thread and any team, on every supported BLAS
configuration (the network's contract; the review's first finding shows it is not yet). Results may
change against the base where a change regroups floating-point sums; they are then judged on the
reads.

## 2. The baseline (2026-11-09, the Dell, 24 threads)

**`learn_pose_student`**, bgal, 1,500 training particles (150 held out), 500 classes, 21 rings,
220-angle half ring:

| run | widths trained | per width (s) | kept | training total (s) |
| --- | --- | --- | --- | --- |
| validation run 2 | 16, 32 | 61, 104 | 16 | 178 |
| validation run 5 | 16, 32, 64 | 48, 115, 351 | 32 | 537 |

Run 2's 16-channel design, split: gradients (the parallel region) 54.5 s of 60.6 (90%), reduction
and Adam steps 3.7 s, held-out evaluation 2.0 s, batch assembly 0.4 s; 391 steps over 23 epochs at
a batch of 84 (17 batches per epoch); 533 samples per second. The rows (the teacher's table fill)
take under 5 s in all.

**Early stopping** keeps the best epoch and stops 250 optimiser steps later (about 15 epochs at 17
batches per epoch). The best epochs of runs 2, 4 and 5 were 4 to 15 of 19 to 30: half or more of
every design's steps train past the model that is kept.

**The truth sets** (the pair student, before this run; wall times under a load of 10 to 20):
`T200` 1.5 min, `T1000` 5.4 min, `T5000` 5.2 min, `T` (24,000) 10 min.

**The workflow test** `test=pose_student` (2,000 simulated particles): passes on the Dell; on the
review's macOS build it passed at one thread in 1,449 s and timed out at eight threads at the
registered 2,400 s, training the 32-channel design at epoch 91.

**The kernels** (`test=pose_student_bench`, 160-angle half ring, 500 classes, 24 threads; two
runs agree within 20%):

| rings | width | forward at one thread (ms) | of it the readout (ms) | backward (ms) | step of 256 at the team (ms) | one thread / team |
| --- | --- | --- | --- | --- | --- | --- |
| 12 | 16 | 4.3 to 5.2 | 3.1 to 3.3 | 4.3 to 6.8 | 140 to 150 | 12 |
| 25 | 16 | 6.4 to 6.7 | 3.3 to 4.1 | 7.4 to 8.4 | 270 to 290 | 8 to 12 |
| 25 | 64 | 13.0 | 4.9 | 20.5 | 1,500 to 1,670 | 12 |
| 37 | 64 | 16.2 to 18.2 | 4.6 to 5.5 | 27.5 to 30.2 | 2,190 to 2,430 | 12 to 13 |

Phase 0 found that this benchmark's "one-thread" forward and backward passes ran with the team's
thread count set, so the OpenMP BLAS spread their matrix products over the team; the training
steps were measured correctly. The corrected single-thread figures and the arithmetic rates are in
Progress, phase 0.

## 3. Where the time goes

- **The schedule trains too much.** Three widths are trained to the end one after the other, and the
  widest costs as much as the other two together, or more; early stopping waits 250 steps (15
  epochs at 17 batches per epoch) after the best epoch.
- **Small batches starve the team.** A mini-batch is cut into chunks of 8 samples (`CNN_CHUNK`),
  computed in parallel and reduced in chunk order (the reduction that makes the step identical on
  any thread count). bgal's batch of 84 (a sixteenth of the training rows) gives 11 chunks: at most
  11 of 24 threads work. A batch of 256 gives 32 chunks, two rounds of which the second keeps 8
  threads busy. The kernel benchmark's one-thread to team ratio of 8 to 13 on 24 cores is the
  symptom. The chunk loop is already scheduled dynamically, so the eight performance cores take
  more chunks than the sixteen efficiency cores.
- **The reduction and the Adam step are serial.** Every chunk's gradient (the templates are most of
  the 160,000 to 920,000 weights) is summed in double on one thread, then one Adam update: 5% of a
  step now, a larger share once the gradients are faster.
- **The readout's rotation search.** Per sample, each of the 500 classes is scored against the
  template and against its reflection: the product of 17 angular harmonics over the readout's
  channels and ring bands, then the score at every rotation of the half ring's grid (`2 x pftsz`
  rotations, 440 at bgal's 220 angles). The rotation scores are about 30 million floating-point
  operations per sample, in hand-written loops at roughly 8 GFLOP/s (one core's single-precision
  peak is about 20 times that); they are 50 to 70% of the forward pass at 16 channels. Each class's
  score is a trigonometric polynomial of order 16 in the rotation angle, so it is smooth on the
  grid's scale.
- **The blocks** are im2col plus `sgemm` (OpenBLAS on the calling thread) where the build's BLAS
  stays on the calling thread inside OpenMP teams (`SIMPLE_BLAS_TEAM_SAFE`, set by
  `cmake/Dependencies.cmake`; on for the Dell's bundled OpenMP OpenBLAS), and plain triple loops
  inside a team otherwise (the review's macOS build with the system OpenBLAS). Their arithmetic
  rate is unmeasured. They dominate at 64 channels.

## 4. The changes, in the order of the phases

### 4.1 The review's findings (phase 1)

Each finding is corrected with the regression coverage the review asks for; its section in the
review states the evidence and the expected behaviour. The run records each finding's resolution
(the change, the test, the evidence) in a closing section of the review file.

- **PS01, the training arithmetic changes with the team** (P1). `blas_ok()` lets `sgemm` run
  outside an active team and selects the explicit loops inside one when the BLAS is not known to be
  team-safe; a one-thread region is not an active team, so one thread and a team of four take
  different arithmetic. The arithmetic path must depend on the build only: with a team-safe BLAS the
  products always go to `sgemm`, without one always to the same loops, whatever the team. The
  HHCNN tester's identity check keeps an actual team of four (it must not pass by running serially),
  and the run checks it on two builds: the Dell's default (team-safe) and a variant configured
  without the definition (section 9.1). A build without a team-safe BLAS then trains on the loops
  everywhere: their speed is part of phase 5.
- **PS02, a model can load another network** (P2). The weights file is resolved against the header's
  directory, never the working directory; the reader checks the binary's shape (and the seed and
  version the header records) against the header and stops on a mismatch. Tests: two model
  directories with the same default names, a foreign weights file in the working directory, a
  relocated model loading its own network.
- **PS03, prediction workers ignore their part** (P2). A worker of `predict_poses` derives its range
  from `part` and `nparts` through the machinery the distributed workers already use (how the qsys
  job descriptions hand `fromp` and `top` to a part), or the program's contract requires the range
  explicitly; a part whose range holds no selected particle writes the empty contribution the merge
  expects instead of stopping. Tests: disjoint coverage by the parts, uneven ranges, a part of only
  unselected rows.
- **PS04, prediction accepts other physical geometry** (P2). Prediction checks the model's
  recorded sampling (`smpd`), crop (`box_crop`) and mask (`mskdiam`) against the project's within
  stated serialisation tolerances, besides the half ring, the box and the CTF treatment it checks
  now; no resampling is implied. Tests: differing sampling, crop and mask with equal `box` and
  `pftsz`.
- **PS05, calibration and prediction select different particles** (P2). The threshold is taken only
  at the end of a group of equal confidences, so that `conf >= threshold` at prediction selects
  exactly the calibrated prefix; it is written with enough digits to round-trip (the `F8.5` form
  does not), and the accepted population is recomputed after reading the model back. Tests: tied
  confidences (the review's 18-and-2 example), values next to the serialisation's boundary, the
  round trip.
- **PS06, a training set can have no held-out cohort** (P2). The held-out split is a deterministic
  choice over the selected particles' groups (stacks, or particles of a single stack) that leaves
  both cohorts non-empty, not `stkind modulo 10`; a selection too small to split is refused before
  training with a clear message; the resampled design comparison handles zero groups. Tests: two to
  nine stacks, selections without a tenth stack, small single-stack selections.
- **PS07, wide-band models cannot be read back** (P2). The reader takes records of any length (or
  the shell weights move to a form without the limit, such as one value per line or the weights
  file); tests: exact round trips on both sides of 227 and 228 shells.

The model file may change for PS02, PS05 and PS07 (its version goes up); release 4 keeps no
compatibility with earlier models. After the corrections the base's three seeds are repeated on
every fixed set (the corrected spread of section 1).

### 4.2 Occupancy and the step (phase 2)

- The chunk size follows the batch alone, never the team (the step must not depend on the thread
  count): chunks of `max(1, ceiling(batch / CNN_NCHUNK))` samples with a named target of 48 to 96
  chunks, the target and the chunk buffers' memory measured (the buffers hold one gradient per
  chunk: at 920,000 weights and 64 chunks, 240 MB).
- The reduction in parallel over slices of the weights: each slice sums the chunks in chunk order,
  so the sums are bit-identical to the serial ones; the Adam update in the same pass.
- The chunk buffers zeroed in parallel, or not at all where the first sample writes them.

### 4.3 The readout (phase 3)

- **As matrix products.** The class products for all classes and both reflections as one complex
  product of the templates with the sample's harmonic spectrum, the rotation scores as one real
  matrix product with the rotation basis (`cos(m theta_r)`, `sin(m theta_r)`), both through
  `sgemm` (through the same build-only rule as PS01); the backward through the winning rotation
  likewise.
- **Coarse to fine.** The score's order-16 polynomial sampled on a coarse grid of a named size (32
  to 64 rotations), the best few coarse cells refined on the full grid around them. Measured: the
  speed, and how often the coarse search's maximum differs from the full grid's on real tensors
  (the bgal and truth-set features), per class and for the mode class. Kept if the reads stay
  within the spread; the full-grid maximum stays the reference in the tester.
- **The backward over negligible classes.** Most classes have a probability and a target near zero;
  their template gradients are near zero. Skipping classes with `|p - t|` below a named bound is an
  approximation: measured for speed and reads, kept only with a gain above 10% and the reads within
  the spread.

### 4.4 The schedule (phase 4)

- **Early stopping in epochs**, with a named patience of a few epochs, and the budget in epochs.
- **A learning-rate schedule** (a short warm-up and a cosine or step decay; named constants) so
  that training converges in fewer epochs; Adam's rate and decay as now unless a measurement says
  otherwise.
- **The batch size** larger where the training set allows it (the rate scaled with it), measured
  for epochs to the best held-out loss and for wall time.
- **The design choice.** Today every width trains to the end and the wider one is kept when its
  held-out recovery gain is resampled clearly above zero. Two candidates: successive halving (every
  width trained a short named budget, only the leader by held-out loss continued, the recovery test
  kept between the finalists), and a rule from the training-set size (the size curve: 64 channels
  never won below some thousands of rows). Each is measured on every set and seed against the full
  search: the width it picks, the reads, the time. The tree keeps one.

### 4.5 The blocks and the rest (phase 5)

What phase 0's arithmetic rates and phases 2 to 4's profiles leave as the largest share: the blocks'
`sgemm` sizes (several samples' positions in one product), im2col reuse between forward and
backward, the gating and pooling loops vectorised (`!$omp simd` where the compiler does not),
memory layout; the feature tensor and the batch assembly if they reach 5% of training. The explicit
loops of a build without a team-safe BLAS (PS01) are measured on the variant build; blocked and
vectorised if they leave training more than three times slower than the team-safe build. If
nothing remaining exceeds 10% of training, phase 5 records that and ends.

### 4.6 Not in this run

The teacher's rows (the table fill), the network's architecture, the prediction table's columns and
the route, the shift registration, warm starts from an earlier model (the streaming use; recorded
as outlook), a GPU, reduced precision, a macOS build (the run has none; the maintainer checks the
review's macOS configuration after the run).

## 5. Decisions already taken

- The review's findings are corrected first (phase 1), each with its regression coverage, before
  any speed change is measured.
- No new command-line key of any kind; every new number is a named constant with its reason in a
  comment.
- No second code path stays in the tree: alternatives are built and measured as variants in the
  run's scratch, and only the kept one is in the diff. The two arithmetic paths of PS01 are the
  build's BLAS configuration, not a choice at run time.
- The training step stays bit-identical between one thread and any team, on both BLAS
  configurations.
- The student may change against the base; it is judged on the reads of section 1 against the seed
  spread, never on one pair of runs.

## 6. Risks

- A correction that changes the reads (PS05's threshold, PS06's split): expected; the corrected
  spread of phase 1 is the reference for the speed phases.
- A faster schedule that stops too early: the reads at `T200` and `T1000` are the first to show it.
- The coarse rotation search that misses a maximum: measured directly, and visible in the in-plane
  read-out of the validation.
- Load from the maintainer's work on the Dell distorting timings: every timing records `uptime`'s
  load and is repeated under load.

## 7. How the run is judged

Every review finding closed with its test, the tests passing on both BLAS configurations. Per fixed
set and seed, the training wall time and the reads at the base, after the corrections and after
each speed phase; the end-to-end validation; the registered workflow test's wall time at its
registered thread count. Hans rules on whether the reads are equal; the run sets no floor. Its
working target is training at least five times faster on every fixed set at equal reads, the kernel
benchmark's training samples per second at least three times the base's.

## 8. Phases

Detailed in section 9.2; in short: 0 baselines and the harness, 1 the review's findings, 2 occupancy
and the step, 3 the readout, 4 the schedule, 5 the blocks and the rest, 6 the end-to-end validation
and the report.

## 9. The run on the Dell

The development runs as the autonomous run `student_speed` on the Dell, with the engine of
`~/agent_runs/ms_regression/engine` and the spec in `~/tasks/student_speed_spec/` (run.conf,
prompt.md, phases.md with the launch procedure, rulings.md), the run root
`/data2/agent_runs/student_speed` behind the link `~/agent_runs/student_speed`. It is launched by
hand on the Dell after Hans has committed this note and the review.

### 9.1 The machine and the data

**The Dell** (Oracle Linux 8.10, user elmlundho): Core Ultra 9 285, 24 cores (8 performance cores,
CPUs 0-7, and 16 efficiency cores; no hyperthreading), 62 GB; `/data2` about 900 GB free, `/home`
little: nothing of the run goes under `/home`. GCC 15.2 from gcc-toolset-15, CMake 3.26,
`python3.12` for CMake; `OPENBLAS_CORETYPE=Haswell` and `TMPDIR=$HOME/tmp` exported. The maintainer
uses the machine during the run: timings record the load and are repeated under load; his checkout
`~/src/SIMPLE` and its `build/` are read-only for the run.

**Builds:** `./compile_clean.sh` about 4 minutes with the machine free; increments by
`make -C build -j20`, `scripts/run_fast_gate.sh build`, `make -C build install`. The default build
links the bundled OpenMP OpenBLAS and defines `SIMPLE_BLAS_TEAM_SAFE`. The variant without it, for
PS01 and phase 5, is a second build tree of the same source configured so that
`cmake/Dependencies.cmake` does not define it (the system OpenBLAS without the OpenMP policy, as the
review's macOS build; the run finds the configure options and records them). The base HEAD is
exported with `git archive` into `scratch/keep/src_base` and built in `scratch/keep/build_base`; HEAD
is never moved.

**Data** (read-only; copied into the run's scratch where a program writes beside its input):

- The truth sets of section 1 under `/data2/agent_runs/pose_student/scratch/keep/medium`:
  `config.sh` names the training half's project and map (`TPROJ`, `TMAP`), the test half `SREF`
  (`S_project/M48noeo.simple`, which carries the truth shifts, so that prediction centres each
  particle: the "centred" test half of the study script; `S_project/1_zero_project_shifts/` holds
  a copy with the shifts removed, whose particles are off centre by up to 10 pixels and which
  phase 0 found unusable for the reads), the subsets under
  `subsets/T200`, `T1000`, `T5000` (one `*_select_particles/M48noeo.simple` each), and
  `pred_stats.awk`. The study script `/data2/agent_runs/mirror_tie/study_pairs.sh` shows the
  commands (`learn_pose_student ... lp=15 nspace=1000 shift_sigma=2 mskdiam=300 pgrp=c1`,
  `predict_poses` on the centred test half).
- The bgal fixed teacher in `/data2/agent_runs/validate_ps/run5/1_validate_pose_student`:
  `learn/training_set.simple` (the training rows at the teacher's poses), the teacher's map
  `1_solve3D/recvol_state01.mrc`, `predict/all_particles.simple` (for the unseen rows only);
  `mskdiam=180`, `pgrp=d2`. Every program runs on a copy.
- The bgal project for the end-to-end validation: `/data2/bgal_addon/3_selection/bgal_addon.simple`
  with `nsample=1500 mskdiam=180 pgrp=d2 pgrp_start=c1 nthr=24`.
- The review's reproductions are described in its sections; their macOS files are not on this
  machine. The run rebuilds each case as a test.

### 9.2 The phases of the run

**Phase 0: baselines and the harness.**

1. Build the tree and the base export; the variant build without `SIMPLE_BLAS_TEAM_SAFE`, and on it
   the HHCNN suite with an actual team of four (PS01 reproduced on Linux, or recorded as not
   reproduced with the configuration used).
2. The kernel benchmark's new columns: the arithmetic rate of the blocks, of the readout's class
   products and of its rotation scores (operation counts from the shapes), and per training step the
   chunks and the threads they can keep busy; the learning bench record's step occupancy (chunks per
   step against the team). Measurement only: the training's results do not change.
3. The fixed sets copied under `scratch/keep/sets/` and a harness script in `scratch/keep/tools/`
   that, for a build and a seed, trains and predicts each set (`T200`, `T1000`, `T5000`, `T`, bgal)
   with `/usr/bin/time -v`, keeps the bench records, and writes one line per set to
   `scratch/keep/summary.tsv` (build, phase, set, seed, wall time, the training split, the widths
   trained and kept, the epochs, the held-out recovery, the reads). The script drives SIMPLE's
   programs and reads their outputs; it replaces none of them.
4. The base: three seeds on every set, the kernel benchmark twice, the registered `pose_student`
   test's wall time at its registered thread count; the spread of every read.

Exit: the base's times, reads and spread in Progress; the kernel table with the arithmetic rates;
PS01 on the variant build recorded; the fast gate and `unit_numerics suite=hhcnn` pass on the
default build.

**Phase 1: the review's findings** (section 4.1). Exit: each of PS01 to PS07 corrected with its
tests, and its resolution recorded in the review file; the fast gate and `unit_numerics` pass on both
builds (the HHCNN identity with an actual team of four on each); `test=pose_student` passes; the
corrected spread (three seeds on every fixed set) in Progress beside the base's.

**Phase 2: occupancy and the step** (section 4.2). Exit: the step identical at one thread and at
24 (the kernel benchmark and the HHCNN tester at 16 threads, repeated) on both builds; the gain per
set and kernel; the reads within the corrected spread.

**Phase 3: the readout** (section 4.3). Exit: the HHCNN tester passes on both builds (the gradient
check, the Hermitian contract, the rotation and mirror tests, the thread identity); the coarse
search's agreement with the full grid measured; the gain; the reads within the corrected spread.

**Phase 4: the schedule** (section 4.4). Exit: every candidate's time and reads on every set and
seed in Progress; one schedule and one design choice kept; the pose student's tester passes.

**Phase 5: the blocks and the rest** (section 4.5). Exit: the profile after phases 2 to 4 on both
builds, each change's gain, the reads within the corrected spread; or the record that nothing
exceeds 10% of training.

**Phase 6: the end-to-end validation and the report.**

1. `validate_pose_student` on bgal twice with the run's build and twice with the base build: the
   learning phase's time and split, the reads, the maps from the student's poses (the teacher
   differs between runs; the comparison is the learning time and whether the reads sit inside the
   runs' spread).
2. The registered `pose_student` test through CTest at its registered thread count and deadline on
   both builds, with its wall time against the base's; the fast gate and `unit_numerics` on both
   builds.
3. The algorithm page (the training section: the schedule, the chunking, the readout's search; the
   model file's weights path, the geometry checks, the threshold over ties, the held-out split), the
   review file's resolution section complete, the report beside this note (what was corrected, what
   was gained where, the evidence, the measurements, a commit message), this note moved to
   `completed/`, its row in `doc/implementation_notes/README.md`.

Exit: all tests pass on both builds; the report and the moved note.

### 9.3 The rules of the run

- Nothing is committed; git is read-only; the base is a `git archive` export.
- No new command-line key, public or internal; no second code path kept (section 5).
- The teacher's rows, the architecture, the prediction table's columns and the route are unchanged
  (section 4.6); the model file changes only as PS02, PS05 and PS07 require, with its version raised.
- New or changed OpenMP code: a tester at 16 threads, repeated, compared bit for bit with one thread,
  on both BLAS configurations; no string temporaries inside parallel regions; work arrays sized by
  the team in use; the BLAS on the calling thread inside a region only where the build is team-safe.
- No `goto`; THROW macros never continued with `//&`; no `.and.` short circuit on allocatable
  arrays; `use ..., only:`; no large type as a plain component; a tester that needs private entities
  is a submodule of the module it tests.
- Timed work records the load; single-thread timings are pinned to the performance cores
  (`taskset -c 0-7`); a change's gain is the mean of three alternating runs, never one pair.
- The data are read-only; bulky outputs are deleted once their numbers are in Progress;
  `scratch/keep/` holds the base export and build, the variant build, the sets, the harness, the
  summary and the logs.
- Write prose for people: Progress entries and the report read without other documents, with
  SIMPLE terms explained; the review's finding IDs may be cited, no other item IDs in the report.

### File table

| file | why | phase |
| --- | --- | --- |
| `doc/implementation_notes/planned/pose_student_training_speed.md` | this plan and its progress (moved to `completed/` at the end of phase 6) | each |
| `doc/implementation_notes/completed/pose_student_training_speed.md` | this plan, completed | 6 |
| `doc/implementation_notes/completed/pose_student_training_speed_report.md` | the run's report | 6 |
| `doc/implementation_notes/completed/learned_pose_state_assignment_review.md` | the review: each finding's resolution | 1 6 |
| `doc/implementation_notes/README.md` | the note's row | 6 |
| `doc/algorithms/pose_student.md` | the model file and prediction contracts (phase 1); the training section | 1 4 6 |
| `doc/policies/test_environment_policy.md` | the kernel benchmark's description, if its columns change it | 0 6 |
| `src/main/learn/simple_hhcnn.f90` | the network: PS01's arithmetic path, the profile counters, the chunking, the reduction, the readout, the blocks | 0 1 2 3 5 |
| `src/main/learn/simple_hhcnn_tester.f90` | its tests | 1 2 3 5 |
| `src/main/learn/simple_hhcnn_bench.f90` | the kernel benchmark | 0 1 2 3 5 |
| `src/main/learn/simple_pose_student.f90` | the model file (PS02, PS05, PS07), the training loop (batch, early stopping, rate schedule), the profile | 0 1 4 5 |
| `src/main/learn/simple_pose_student_tester.f90` | its tests | 1 4 |
| `src/main/learn/simple_polar_feats.f90` | the feature tensor, if it reaches 5% of training | 5 |
| `src/main/strategies/parallelization/simple_strategy_pose_student.f90` | PS03 to PS06; the design choice; the learning bench record | 0 1 2 4 5 |
| `src/main/strategies/parallelization/simple_strategy_pose_student_tester.f90` | the tests of the strategy's private helpers (split, threshold, range, geometry), a submodule | 1 4 |
| `src/main/commanders/simple/simple_commanders_pose_student.f90` | the programs' wrappers (PS03's range) | 1 |
| `src/main/ui/simple_private_prgs.f90` | `predict_poses`' keys (PS03's contract) | 1 |
| `src/main/commanders/test/simple_commanders_test_class.f90` | the suite table, for a new tester | 1 |
| `src/main/ui/simple_test/simple_test_ui_class.f90` | the suites' help, for a new tester | 1 |
| `src/main/commanders/test/simple_commanders_test_highlevel_pose.f90` | the kernel benchmark program; the workflow test's program-level checks (PS02 to PS04) | 0 1 3 |
| `src/main/ui/simple_test/simple_test_ui_highlevel.f90` | the kernel benchmark's help | 0 3 |

## Progress

### Phase 0: baselines and the harness (2026-11-09 to 2026-11-10, done)

**What changed in the tree** (measurement only; the training's results are unchanged):
`src/main/learn/simple_hhcnn.f90` (the network's cost record counts each training step's chunks
and the chunk slots of the team that ran it; the readout's forward pass is split into the readout
map with its spectrum, `readout_spectrum`, and one class's spectral product, `class_product`,
called in the same order with the same arithmetic, so that the benchmark can time them alone; the
matrix-product kernels are bound to the type so that the benchmark submodule links them),
`src/main/learn/simple_hhcnn_bench.f90` (the single-thread measurements now run at one OpenMP
thread; the kernels timed alone; a second table of operation counts and arithmetic rates and the
step's occupancy), `src/main/strategies/parallelization/simple_strategy_pose_student.f90` (the
learning bench record's line "chunks per step / team / occupancy" per design),
`src/main/ui/simple_test/simple_test_ui_highlevel.f90` and
`doc/policies/test_environment_policy.md` (the benchmark's description), this note (section 9.1's
test half corrected, a note under section 2's kernel table, this entry).

**Builds.** GCC 15.2, CMake 3.26, Release with `-O3 -march=native`. The tree (`repo/build`) and
the base export of `6541633df` (`scratch/keep/build_base`) link the system OpenMP OpenBLAS 0.3.15
(`/usr/lib64/libopenblaso.so`; the BLAS policy prefers it to the bundled build) and print
`BLAS in OMP teams: ON`. The variant without `SIMPLE_BLAS_TEAM_SAFE`
(`scratch/keep/build_noteam`, the tree's source) is configured with
`-DSIMPLE_BLAS_OPENMP=OFF -DBLAS_LIBRARIES=/usr/lib64/libopenblasp.so
-DLAPACK_LIBRARIES=/usr/lib64/libopenblasp.so` (the system pthreads OpenBLAS; this machine's
`libopenblas.so` is the single-threaded build) and prints `BLAS in OMP teams: OFF`. Configure
lines and logs: `scratch/logs/cmake_*.log`, the journal.

**The review's first finding reproduced on Linux.** On the variant build,
`OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=4 OMP_WAIT_POLICY=PASSIVE OPENBLAS_NUM_THREADS=1
simple_test_exec test=unit_numerics suite=hhcnn` fails 1 of 31 checks: "the chunked reduction
gives the same weights at one thread and at a team of 4, bit for bit"
(`scratch/keep/logs/ps01_hhcnn_team4_noteam.log`); with `OMP_THREAD_LIMIT=1` it passes 31 of 31,
the serial diagnostic that proves nothing (`ps01_hhcnn_limit1_noteam.log`). The default build
passes 31 of 31 with the actual team of four (`ps01_hhcnn_team4_tree.log`).

**Identity of the change.** `learn_pose_student` on `T200` at seed 1 and 24 threads with the base
and with the tree: the weights files are identical byte for byte and the model headers differ only
in the wall time and the project path (`scratch/keep/logs/idcheck_p0.log`). The fast gate passes
on the tree (15 of 15, `scratch/keep/logs/fast_gate_tree_p0_final.log`), and
`unit_numerics suite=hhcnn` passes 31 of 31 at its default and with a team of four
(`hhcnn_tree_p0_final*.log`).

**A flaw of the base's benchmark.** Its forward and backward passes "at one thread" ran with the
team's thread count set, outside a parallel region, where the OpenMP BLAS spreads one product over
the team: the blocks' matrix products appeared at 330 to 740 GFLOP/s, several times one core's
peak (`scratch/keep/bench/kernels_tree_blasthreaded.txt`). The tree's benchmark sets one thread
around them. The training steps were always measured correctly (the one-thread step inside a
one-thread team, the team's step inside the team). Section 2's per-sample columns are therefore
too low at 32 and 64 channels; the corrected ones follow.

**The kernels** (`simple_test_exec test=pose_student_bench nthr=24`, 160-angle half ring, 500
classes, mean of two runs that agree within 5%; the single-thread figures on CPU 0, a performance
core, by `OMP_PLACES=cores OMP_PROC_BIND=close`, the main thread's affinity read as CPU 0 while it
ran; load 12 to 17, mostly the benchmark's own team; `scratch/keep/bench/kernels_tree_{1,2}.txt`):

| rings | width | forward (ms) | of it the readout | of it the blocks | backward (ms) | step of 256 at one thread (ms) | at 24 threads (ms) | ratio |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 12 | 16 | 3.3 | 2.5 | 0.8 | 1.7 | 1,326 | 142 | 9.3 |
| 12 | 32 | 4.7 | 2.6 | 2.2 | 4.1 | 2,336 | 258 | 9.1 |
| 12 | 64 | 9.7 | 2.7 | 7.0 | 14.2 | 6,117 | 643 | 9.5 |
| 18 | 16 | 4.1 | 2.9 | 1.2 | 2.6 | 1,727 | 188 | 9.2 |
| 18 | 32 | 6.0 | 3.0 | 3.1 | 6.4 | 3,172 | 362 | 8.8 |
| 18 | 64 | 12.9 | 3.2 | 9.7 | 20.9 | 8,684 | 910 | 9.5 |
| 25 | 16 | 4.7 | 3.2 | 1.5 | 3.8 | 2,190 | 247 | 8.9 |
| 25 | 32 | 7.6 | 3.3 | 4.2 | 9.7 | 4,387 | 510 | 8.6 |
| 25 | 64 | 17.9 | 3.7 | 14.2 | 30.0 | 12,265 | 1,308 | 9.4 |
| 37 | 16 | 5.7 | 3.5 | 2.2 | 5.7 | 2,911 | 350 | 8.3 |
| 37 | 32 | 10.0 | 3.7 | 6.3 | 14.9 | 6,338 | 779 | 8.1 |
| 37 | 64 | 26.2 | 4.2 | 22.0 | 44.0 | 17,993 | 1,968 | 9.1 |

The one-thread forward and backward per sample times 256 give the one-thread step (at 37 rings
and 64 channels, (26.2 + 44.0) ms x 256 = 18.0 s against 18.0 s measured): the per-sample and
the step figures agree.

**The arithmetic rates** (operations counted from the shapes: the blocks' products
`2 cout 9 cin 2 npos` per block, forward once and backward twice except the first block; the class
products 8 operations per complex multiply-accumulate over 17 orders, the bands, 8 channels, 500
classes and both reflections; the rotation scores 4 operations per order and rotation, 2 x 500
classes x 320 rotations x 16 orders = 20.5 million per sample; the gates, the pooling, im2col and
the softmax not counted). GFLOP/s is billions of floating-point operations per second; one
performance core's single-precision peak is about 170.

| rings | width | blocks' operations forward / backward (millions per sample) | readout (millions) | blocks' products alone at one thread | the same backward | the blocks with im2col, gate and pooling | class products | rotation scores | step at one thread | step at 24 threads |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 12 | 16 | 54 / 107 | 21.6 | 140 | 132 | 70 | 3.6 | 9.6 | 36 | 332 |
| 12 | 32 | 215 / 427 | 21.6 | 145 | 144 | 100 | 3.5 | 9.6 | 73 | 660 |
| 12 | 64 | 854 / 1,703 | 21.6 | 149 | 150 | 122 | 3.5 | 9.5 | 108 | 1,028 |
| 18 | 16 | 75 / 149 | 22.7 | 133 | 129 | 65 | 3.6 | 9.4 | 37 | 341 |
| 18 | 32 | 298 / 593 | 22.7 | 143 | 143 | 97 | 3.6 | 9.5 | 74 | 650 |
| 18 | 64 | 1,186 / 2,366 | 22.7 | 149 | 149 | 122 | 3.7 | 9.6 | 106 | 1,007 |
| 25 | 16 | 108 / 215 | 23.7 | 135 | 131 | 73 | 3.8 | 9.6 | 41 | 365 |
| 25 | 32 | 429 / 854 | 23.7 | 144 | 144 | 101 | 3.8 | 9.6 | 77 | 660 |
| 25 | 64 | 1,708 / 3,407 | 23.7 | 149 | 148 | 120 | 3.9 | 9.5 | 107 | 1,008 |
| 37 | 16 | 163 / 322 | 24.8 | 135 | 130 | 74 | 3.9 | 9.6 | 45 | 377 |
| 37 | 32 | 644 / 1,281 | 24.8 | 143 | 141 | 102 | 3.9 | 9.6 | 79 | 644 |
| 37 | 64 | 2,562 / 5,110 | 24.8 | 145 | 145 | 116 | 4.0 | 9.5 | 110 | 1,003 |

What the rates say: the blocks' products run at 75 to 90% of a core's peak, so the blocks gain
only from doing less (im2col, gate and pooling cost half of their time at 16 channels and a fifth
at 64). The readout's rotation scores run at about 6% of peak and the class products at about
2%; together they are most of the forward pass at 16 and 32 channels (2.5 to 4.2 ms of 3.3 to 10
ms). A step of 256 samples is 32 chunks of 8: on 24 threads two rounds, the second a third full
(occupancy 0.67), and the team runs the step 8 to 9.5 times faster than one thread.

**Occupancy of the real training steps.** The batch is a sixteenth of the training rows between
16 and 256: `T200` takes batches of 16 (2 chunks per step, 1.9 on average; occupancy 0.08 on 24
threads, from the tree's learning bench record of the identity run), `T1000` 56 (7 chunks, 0.29),
bgal 84 (11 chunks, 0.46), `T5000` and `T` 256 (32 chunks, 0.67).

**The harness** (`scratch/keep/tools/harness.sh <build> <bin> <phase> <seed> [sets]`): per set
it trains with `learn_pose_student` under `/usr/bin/time -v` at 24 threads and `SIMPLE_SEED`,
predicts with `predict_poses`, reads the prediction table (`pred_stats.awk`, and the seeded
fraction: the rows the route sends to the neighbourhood search), and appends one line to
`scratch/keep/summary.tsv` (wall time, memory, the learning bench record's split per design with
epochs, batch, steps and occupancy, the widths trained and kept, the kept model's epochs and best
epoch, its held-out recovery, threshold and seeded fraction, the reads, the load at start and
end). Model headers, bench records, logs, read summaries and the compressed prediction tables are
kept under `scratch/keep/runs/<build>_p<phase>/<set>_s<seed>/`; the projects, feature caches,
rows and weights are deleted. `scratch/keep/tools/spread.sh` prints the mean and range over
seeds. The sets are copied under `scratch/keep/sets/` (`MD5SUMS`).

**The test half.** The first pass predicted the copy that section 9.1 named as the centred test
half, `S_project/1_zero_project_shifts/M48noeo.simple`; that copy has the truth shifts removed,
its particles sit up to 10 pixels off centre, and the student found the truth within 6 degrees
for 3 to 5% of them (`scratch/keep/summary_wrongS.tsv`). The study script's "centred" half is
`S_project/M48noeo.simple`, which carries the truth shifts so that prediction centres the
particles; with it `T200` at seed 1 reads 0.690 within 6 degrees against the study's 0.696. The
truth sets were rerun on it (the models are identical to the first pass's: the same seed gives
the same held-out record), and section 9.1 now names it.

**The base, three seeds on every set** (24 threads; load 4 to 26 at the start of a run, most of
it the previous run's team; no other heavy process was seen). Reads on the truth sets are the
student's folded mode against the truth on the 24,000 particles of the test half; on bgal they are
the 1,500 training rows' mode against the teacher's direction, and the seeded fraction of the
1,256 particles the student never saw. "Seeded" is the share of particles the route sends to the
neighbourhood search (confidence at or above the threshold, cone within the cap and the bound);
"in cone" is the share whose truth (teacher) direction lies in the student's cone.

| set | learning wall (s) | kept widths | median error (deg) | within 6 deg | within 12 deg | in cone | seeded | seeded with the truth in the cone | held-out seeded |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `T200` | 89 [79, 100] | 16, 16, 16 | 3.8 [3.5, 4.1] | 0.743 [0.690, 0.792] | 0.839 [0.800, 0.870] | 0.968 [0.962, 0.973] | 0.204 [0.117, 0.305] | 0.982 [0.980, 0.983] | 0.55 [0.50, 0.65] |
| `T1000` | 192 [118, 333] | 16, 32, 16 | 3.6 [3.5, 3.8] | 0.797 [0.756, 0.823] | 0.895 [0.868, 0.913] | 0.982 [0.974, 0.987] | 0.241 [0.219, 0.260] | 0.987 [0.979, 0.993] | 0.50 [0.48, 0.54] |
| `T5000` | 468 [310, 784] | 16, 64, 16 | 3.3 [3.1, 3.5] | 0.881 [0.837, 0.914] | 0.957 [0.935, 0.974] | 0.986 [0.979, 0.992] | 0.443 [0.370, 0.509] | 0.986 [0.981, 0.990] | 0.61 [0.55, 0.66] |
| `T` | 606 [586, 630] | 16, 16, 16 | 3.1 [3.0, 3.2] | 0.930 [0.914, 0.946] | 0.993 [0.992, 0.995] | 0.995 [0.991, 0.997] | 0.532 [0.464, 0.586] | 0.992 [0.987, 0.995] | 0.54 [0.47, 0.59] |
| bgal | 351 [169, 457] | 16, 32, 32 | 2.5 [2.5, 2.5] | 0.874 [0.871, 0.877] | 0.929 [0.924, 0.931] | 0.960 [0.956, 0.966] | 0.424 [0.414, 0.437] | 0.996 [0.995, 0.997] | 0.36 [0.33, 0.39] |

Mean [minimum, maximum] over seeds 1, 2 and 3. bgal's unseen particles: seeded 0.058 [0.045,
0.071]. The held-out recovery is 1.000 on every set but `T5000` (0.997) and `T` (0.999); the
confidence threshold 0.990 to 0.991 everywhere.

Per seed (the widths trained and kept, the kept design's best epoch of the epochs run):

| set | seed | learning wall (s) | widths trained | kept | best / epochs | median error (deg) | within 6 deg | in cone | seeded (bgal: training / unseen) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `T200` | 1 | 79 | 16, 32 | 16 | 3 / 24 | 4.08 | 0.690 | 0.962 | 0.117 |
| `T200` | 2 | 100 | 16, 32 | 16 | 5 / 26 | 3.79 | 0.747 | 0.970 | 0.190 |
| `T200` | 3 | 88 | 16, 32 | 16 | 5 / 26 | 3.51 | 0.792 | 0.973 | 0.305 |
| `T1000` | 1 | 118 | 16, 32 | 16 | 3 / 18 | 3.78 | 0.756 | 0.974 | 0.219 |
| `T1000` | 2 | 333 | 16, 32, 64 | 32 | 4 / 19 | 3.51 | 0.813 | 0.984 | 0.245 |
| `T1000` | 3 | 126 | 16, 32 | 16 | 5 / 20 | 3.46 | 0.823 | 0.987 | 0.260 |
| `T5000` | 1 | 310 | 16, 32 | 16 | 6 / 20 | 3.47 | 0.837 | 0.979 | 0.370 |
| `T5000` | 2 | 784 | 16, 32, 64 | 64 | 4 / 18 | 3.13 | 0.913 | 0.992 | 0.509 |
| `T5000` | 3 | 311 | 16, 32 | 16 | 7 / 21 | 3.21 | 0.892 | 0.989 | 0.452 |
| `T` | 1 | 586 | 16, 32 | 16 | 4 / 7 | 3.09 | 0.930 | 0.996 | 0.464 |
| `T` | 2 | 630 | 16, 32 | 16 | 6 / 9 | 3.17 | 0.914 | 0.991 | 0.547 |
| `T` | 3 | 602 | 16, 32 | 16 | 5 / 8 | 3.01 | 0.946 | 0.997 | 0.586 |
| bgal | 1 | 169 | 16, 32 | 16 | 6 / 21 | 2.45 | 0.874 | 0.958 | 0.419 / 0.045 |
| bgal | 2 | 457 | 16, 32, 64 | 32 | 5 / 20 | 2.52 | 0.871 | 0.956 | 0.414 / 0.071 |
| bgal | 3 | 428 | 16, 32, 64 | 32 | 5 / 20 | 2.47 | 0.877 | 0.966 | 0.437 / 0.059 |

**Where the training time goes** (the learning bench records): per design, 16 channels take 24
to 29 s on `T200`, 35 to 40 s on `T1000`, 96 to 106 s on `T5000`, 170 to 216 s on `T` and 50 to
52 s on bgal; 32 channels about twice that, 64 channels about five times (205 s on `T1000`, 459 s
on `T5000`, 255 to 284 s on bgal). Every design trains 18 to 26 epochs on the subsets and bgal and
keeps an epoch between 3 and 7, so three quarters of its epochs train past the kept model; on `T`
(85 batches per epoch) the patience of 250 steps is three epochs and it trains 7 to 9 epochs,
keeping epoch 4 to 6. Of a design's time the parallel gradients are about 90%, the
reduction and Adam steps 4 to 8%, the held-out evaluation 3 to 5%, the batch assembly under 1%. The
teacher's rows take under 1 s on the subsets and bgal and 56 s on `T` (outside the scope). A third
width is trained when the second wins the resampled recovery test: one seed of three on `T1000`
and `T5000`, two of three on bgal, and it doubles or nearly triples the wall time; the seed spread of the wall
time is mostly that choice.

**The registered workflow test** `pose_student` on the base through CTest at its registered 8
threads and 2,400 s deadline: passed in 313.9 s, of it 296 s training
(`scratch/keep/logs/ctest_pose_student_base.log`, `pose_student_metrics_base.tsv`).

**Exit items.** The base's times, reads and spread: above. The kernel table with the arithmetic
rates: above. The review's first finding on the variant build: reproduced, above. The fast gate and
`unit_numerics suite=hhcnn` on the default build: pass. The run's chain logs:
`scratch/keep/logs/chain_p0.log`, `chain_p0b.log`.

### Phase 1: the review's findings (2026-11-10, done)

**What changed.** The seven findings of the code review are corrected, each with its tests; the review
file's closing section "Resolution in student_speed" records every change, test and log. In short:

- *The training arithmetic no longer depends on the thread team* (PS01): the network's matrix products
  take `sgemm` when the build's BLAS stays on the calling thread inside OpenMP teams and the explicit
  loops otherwise, whatever the team (`simple_hhcnn.f90`). The HHCNN tester's identity check now also
  asserts that its team of four formed. Reproduced before the change on the variant build (30 of 31),
  passing after it on both builds (32 of 32 with an actual team of four); a one-thread cap now fails the
  team check instead of passing silently.
- *A model loads its own network* (PS02): the header names the weights by file name, resolved in the
  header's directory, and the reader checks the network's shape, seed and weight count against the
  header (`simple_pose_student.f90`). The pose student's tester became a submodule of its module to
  reach these checks.
- *Prediction parts take an explicit range* (PS03): a part needs `nparts`, `fromp`, `top` and
  `outfile`, as the distributed masters hand them to every part; a part of only unselected rows writes
  its records unchanged; a part's tables carry its number (`simple_commanders_pose_student.f90`, the
  strategy, `simple_private_prgs.f90`: `fromp` and `top` added to the private program's optional keys,
  both existing keys).
- *Prediction checks the physical geometry* (PS04): sampling within 1e-4 Å, mask within 0.01 Å, crop
  exactly, besides the half ring, box and CTF.
- *Calibration and prediction select the same particles* (PS05): the threshold ends a group of equal
  confidences; the threshold and the cone bound are written to round-trip; learn_pose_student reads
  the model back and checks the seeded held-out rows.
- *Every training set has a held-out cohort or is refused* (PS06): every tenth of the selected
  particles' stacks when there are ten or more, else every tenth selected particle; fewer than ten
  particles are refused; the bootstrap of zero groups gives zero; training refuses a set without
  held-out rows.
- *Wide-band models read back* (PS07): the header reader takes records of any length.

The model file's version is 4. New tests: the strategy's helpers in
`simple_strategy_pose_student_tester.f90` (a submodule; suite `pose_student_strategy` of
`unit_numerics`, 41 checks; the helpers are public module procedures, as the test environment policy's
section 4.3 requires for a submodule's calls into a procedure module under GCC 15); the pose student's
suite grew from 38 to 64 checks; the workflow test `pose_student` gained a relocated-model check and
four checks of a distributed caller's parts, including `merge_algndocs` of their documents. A review
by a separate agent of the diff found that the read-back check used the in-memory cone bound while the
file kept three decimals, that parts sharing `fname` would overwrite each other's tables, and that an
empty part skipped the refusal of `vol1`; all three were corrected before the final runs.

**Files changed in this phase:** `src/main/learn/simple_hhcnn.f90`, `simple_hhcnn_tester.f90`,
`simple_pose_student.f90`, `simple_pose_student_tester.f90`,
`src/main/strategies/parallelization/simple_strategy_pose_student.f90`,
`simple_strategy_pose_student_tester.f90` (new), `src/main/commanders/simple/simple_commanders_pose_student.f90`,
`src/main/commanders/test/simple_commanders_test_class.f90`,
`src/main/commanders/test/simple_commanders_test_highlevel_pose.f90`, `src/main/ui/simple_private_prgs.f90`,
`src/main/ui/simple_test/simple_test_ui_class.f90`, `doc/algorithms/pose_student.md`, the review file,
this note.

**Exit items.**
- Each of PS01 to PS07 corrected with its tests, its resolution in the review file: done.
- The fast gate and `unit_numerics` pass on both builds: 15 of 15 on the default build
  (`scratch/keep/logs/p1_fast_gate_tree.log`) and on the variant (`p1_fast_gate_noteam.log`); the
  HHCNN identity with an actual team of four, 32 of 32 on each (`p1_ps01_hhcnn_team4_{tree,noteam}.log`).
- `test=pose_student` passes: through CTest at its registered eight threads in 317.7 s (the base 313.9
  s), 13 of 13 checks (`p1_ctest_pose_student_tree_final.log`, `p1_pose_student_metrics_tree.tsv`).
- The corrected spread, three seeds on every fixed set, beside the base's: below.

**The corrected spread** (`scratch/keep/tools/harness.sh tree ... 1 <seed>`, 24 threads, load 5 to 25
at the start of a run, mostly the previous run's team; `scratch/keep/summary.tsv`, build `tree`,
phase 1; logs under `scratch/keep/runs/tree_p1/`). The corrections change no trained network on these
sets: the held-out tables (each held-out particle's confidence, cone and recovery) are identical to the
base's for all fifteen runs, so the splits, the kept widths, the epochs and the truth-set reads are the
base's to every printed digit. What moves is the threshold, which now keeps its full precision: the
five-decimal form had rounded it both ways (0.9900014 down to 0.99000 on `T` seed 1, 0.9901377 up to
0.99014 on bgal seed 2), and the seeded fractions change by at most 0.0013.

| set | learning wall (s) | kept widths | within 6 deg | in cone | seeded, base | seeded, corrected | held-out seeded |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `T200` | 101 [99, 105] | 16, 16, 16 | 0.743 [0.690, 0.792] | 0.968 [0.962, 0.973] | 0.2039 [0.1169, 0.3050] | 0.2038 [0.1169, 0.3050] | 0.55 [0.50, 0.65] |
| `T1000` | 197 [123, 340] | 16, 32, 16 | 0.797 [0.756, 0.823] | 0.982 [0.974, 0.987] | 0.2414 [0.2188, 0.2601] | 0.2414 [0.2187, 0.2602] | 0.50 [0.48, 0.54] |
| `T5000` | 472 [312, 789] | 16, 64, 16 | 0.881 [0.837, 0.914] | 0.986 [0.979, 0.992] | 0.4434 [0.3700, 0.5087] | 0.4434 [0.3698, 0.5087] | 0.61 [0.55, 0.66] |
| `T` | 608 [588, 631] | 16, 16, 16 | 0.930 [0.914, 0.946] | 0.995 [0.991, 0.997] | 0.5323 [0.4640, 0.5857] | 0.5321 [0.4637, 0.5853] | 0.54 [0.47, 0.59] |
| bgal | 347 [160, 455] | 16, 32, 32 | 0.874 [0.871, 0.877] | 0.960 [0.956, 0.966] | 0.4235 [0.4140, 0.4373] | 0.4240 [0.4153, 0.4367] | 0.36 [0.33, 0.39] |

Mean [minimum, maximum] over seeds 1 to 3; bgal's reads are the training rows against the teacher.
bgal's unseen particles seeded: 0.0581 [0.0454, 0.0709] (base 0.0584 [0.0454, 0.0709]). The learning
wall times are within 1% of the base's on `T5000`, `T` and bgal, 2% longer on `T1000` and 13% on
`T200`, the shortest runs (seed by seed 79, 100, 88 s at the base against 99, 105, 99 s; a build of the
review's corrections ran once during seed 1's `T1000` and `T5000`).

Per seed (the corrected tree; widths trained and kept, the kept design's best epoch of the epochs run):

| set | seed | learning wall (s) | widths trained | kept | best / epochs | median error (deg) | within 6 deg | in cone | seeded (bgal: training / unseen) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `T200` | 1 | 99 | 16, 32 | 16 | 3 / 24 | 4.08 | 0.690 | 0.962 | 0.117 |
| `T200` | 2 | 105 | 16, 32 | 16 | 5 / 26 | 3.79 | 0.747 | 0.970 | 0.190 |
| `T200` | 3 | 99 | 16, 32 | 16 | 5 / 26 | 3.51 | 0.792 | 0.973 | 0.305 |
| `T1000` | 1 | 123 | 16, 32 | 16 | 3 / 18 | 3.78 | 0.756 | 0.974 | 0.219 |
| `T1000` | 2 | 340 | 16, 32, 64 | 32 | 4 / 19 | 3.51 | 0.813 | 0.984 | 0.245 |
| `T1000` | 3 | 128 | 16, 32 | 16 | 5 / 20 | 3.46 | 0.823 | 0.987 | 0.260 |
| `T5000` | 1 | 314 | 16, 32 | 16 | 6 / 20 | 3.47 | 0.837 | 0.979 | 0.370 |
| `T5000` | 2 | 789 | 16, 32, 64 | 64 | 4 / 18 | 3.13 | 0.913 | 0.992 | 0.509 |
| `T5000` | 3 | 312 | 16, 32 | 16 | 7 / 21 | 3.21 | 0.892 | 0.989 | 0.452 |
| `T` | 1 | 588 | 16, 32 | 16 | 4 / 7 | 3.09 | 0.930 | 0.996 | 0.464 |
| `T` | 2 | 631 | 16, 32 | 16 | 6 / 9 | 3.17 | 0.914 | 0.991 | 0.547 |
| `T` | 3 | 606 | 16, 32 | 16 | 5 / 8 | 3.01 | 0.946 | 0.997 | 0.585 |
| bgal | 1 | 160 | 16, 32 | 16 | 6 / 21 | 2.45 | 0.874 | 0.958 | 0.420 / 0.045 |
| bgal | 2 | 455 | 16, 32, 64 | 32 | 5 / 20 | 2.52 | 0.871 | 0.956 | 0.415 / 0.071 |
| bgal | 3 | 424 | 16, 32, 64 | 32 | 5 / 20 | 2.47 | 0.877 | 0.966 | 0.437 / 0.058 |

The corrected spread is the base's: the speed phases are judged against these numbers. Two things
for them: the held-out split of the fixed sets is unchanged (the truth sets lie on one stack, and bgal's
selected stacks are numbered consecutively, so the ranked rule picks the stacks the old one did); and a
build without a team-safe BLAS now trains on the explicit loops at every thread count (its HHCNN suite
takes 26 s against 7 s), which phase 5 measures.

### Phase 2: occupancy and the step (2026-11-10, done)

**What changed** (`src/main/learn/simple_hhcnn.f90`, its tester and the kernel benchmark):

- *Chunks follow the batch alone.* A mini-batch is cut into chunks of `max(1, ceiling(n / CNN_NCHUNK))`
  samples with `CNN_NCHUNK = 96`; the chunk count depends on the batch, never on the team, so a step is
  the same at any thread count. Batches of 16 (`T200`), 56 (`T1000`) and 84 (bgal) now give one sample
  per chunk (16, 56, 84 chunks, against 2, 7 and 11 before); a batch of 256 (`T5000`, `T`) gives 86
  chunks of three (32 of eight before). The learning bench records' occupancy (the share of the team's
  chunk slots the chunks fill) rose from 0.08 to 0.63 on `T200`, 0.29 to 0.77 on `T1000`, 0.43 to 0.87
  on bgal and 0.67 to 0.89-0.90 at 256.
- *The reduction in parallel.* The kernels' columns and the templates' slices of 4,096 weights
  (`CNN_SLICE`) are summed over the chunks in chunk order by one task each, which gives the serial sums
  bit for bit at any team; the templates' Adam step runs in the same pass, the kernels' after their
  symmetric offsets are tied. The Adam updates became elemental (the same formulas). The chunk buffers
  were already zeroed inside the parallel region, by the chunk that fills them.
- *The tests.* The HHCNN tester's identity check runs at a team of four (the fast gate's cap) and at
  sixteen where the thread limit allows it, each asserting that its team formed; the benchmark counts
  the chunks from the batch.

**The chunk target.** Measured on the kernel benchmark (the training step of 256 samples at 24 threads,
summed over its twelve shapes; `scratch/keep/bench/kernels_*`), each target built as a variant of the
tree (patches in `scratch/keep/variants/`): the base 7,670 and 7,659 ms; 48 chunks 8,014 ms; 64 chunks
7,375, 7,650, 7,641 ms; 96 chunks 7,246, 7,306, 7,260 ms (the 64 and 96 runs alternated). 96 is faster
in every run, by 3.8% on the means, and is kept. Its chunk buffers hold one gradient per chunk: at the
largest benchmark design (37 rings, 64 channels, 915,713 weights) 86 chunks are 315 MB, at bgal's
widest (21 rings, 64 channels) about 220 MB. At a batch of 256 the team's step is only 2 to 5% faster
than the base's: the base's second round of eight chunks already ran on the eight performance cores
(the threads are bound close, so the first threads sit on CPUs 0-7), and a step at 24 threads runs 9
to 10 times faster than at one, about what eight performance and sixteen efficiency cores give.

**Identity.** The training step is bit-identical between one thread and a team on both builds: the
HHCNN tester at a team of four (31 of 31) and, by hand with `OMP_THREAD_LIMIT=24`, at sixteen, three
times each (34 of 34; `scratch/keep/logs/p2_hhcnn_{tree,noteam}_{team4,t16_rep1-3}.log`); the kernel
benchmark's step identical at 24 threads for every shape on the default build (six runs,
`kernels_tree_p2_*`, `kernels_nchunk96_*`) and on the variant build for the three 12-ring shapes
(`kernels_noteam_p2_1_partial.txt`): its explicit loops made the remaining nine shapes a matter of
hours (the one-thread step at 12 rings and 64 channels took 217 s against 6.1 s with `sgemm`), so the run
was stopped there; the identity is a property of the code path, which the tester at four and sixteen
threads covers on that build. The fast gate passes on both builds (`p2_fast_gate_{tree,noteam}.log`).

**The gain per set** (the learning wall time, mean [minimum, maximum] over seeds 1 to 3, 24 threads,
load 2 to 25 at a run's start, mostly the previous run's team; per design the time of the 16-channel
design, which every run trains):

| set | base (s) | phase 2 (s) | gain | 16-channel design, base (s) | phase 2 (s) | its reduction and Adam steps, base / phase 2 (s) |
| --- | --- | --- | --- | --- | --- | --- |
| `T200` | 89 [79, 100] | 27 [26, 29] | 3.3 | 23.5-29.3 | 6.1-7.2 | 0.9-1.6 / 0.4-0.5 |
| `T1000` | 192 [118, 333] | 106 [69, 179] | 1.8 | 35.3-40.2 | 19.9-21.9 | 1.7-2.5 / 1.2-1.3 |
| bgal | 351 [169, 457] | 241 [118, 313] | 1.5 | 49.8-52.4 | 34.3-35.9 | 2.3-2.5 / 1.8-1.9 |
| `T5000` | 468 [310, 784] | 466 [304, 786] | 1.0 | 95.6-105.9 | 90.3-99.5 | 6.5-8.2 / 2.2 |
| `T` | 606 [586, 630] | 597 [574, 621] | 1.0 | 169.7-216.0 | 158.5-203.8 | 12.2-13.9 / 3.5-4.5 |

The parallel reduction makes the reduction and Adam steps three to four times faster at a batch of 256;
the gradients at 256 take as long as before. The rows (the teacher's table fill) are unchanged.

**The reads** (`scratch/keep/summary.tsv`, build `tree`, phase 2; mean [minimum, maximum]; the
corrected spread of phase 1 is the base's to the printed digits):

| set | within 6 deg | within 12 deg | in cone | seeded | held-out seeded | base: within 6 deg | base: in cone | base: seeded |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `T200` | 0.743 [0.690, 0.792] | 0.839 [0.800, 0.870] | 0.968 [0.962, 0.973] | 0.215 [0.117, 0.334] | 0.55 [0.50, 0.65] | 0.743 [0.690, 0.792] | 0.968 [0.962, 0.973] | 0.204 [0.117, 0.305] |
| `T1000` | 0.797 [0.756, 0.823] | 0.895 [0.869, 0.912] | 0.982 [0.974, 0.987] | 0.243 [0.223, 0.262] | 0.50 [0.48, 0.54] | 0.797 [0.756, 0.823] | 0.982 [0.974, 0.987] | 0.241 [0.219, 0.260] |
| `T5000` | 0.881 [0.836, 0.913] | 0.957 [0.934, 0.974] | 0.986 [0.978, 0.992] | 0.443 [0.370, 0.508] | 0.61 [0.55, 0.66] | 0.881 [0.837, 0.914] | 0.986 [0.979, 0.992] | 0.443 [0.370, 0.509] |
| `T` | 0.930 [0.914, 0.946] | 0.993 [0.992, 0.995] | 0.995 [0.991, 0.996] | 0.532 [0.464, 0.583] | 0.54 [0.47, 0.59] | 0.930 [0.914, 0.946] | 0.995 [0.991, 0.997] | 0.532 [0.464, 0.586] |
| bgal | 0.873 [0.871, 0.877] | 0.928 [0.921, 0.931] | 0.961 [0.957, 0.966] | 0.396 [0.377, 0.431] | 0.34 [0.31, 0.39] | 0.874 [0.871, 0.877] | 0.960 [0.956, 0.966] | 0.424 [0.414, 0.437] |

bgal's unseen particles seeded: 0.054 [0.041, 0.063] (base 0.058 [0.045, 0.071]); held-out recovery
1.000 on bgal and the subsets but `T5000` (0.997) and `T` (0.999), as at the base. The new summation
order trains slightly different networks (the sums of the chunks regroup); the agreement with the truth
and with the teacher, the truth in the cone, the held-out recovery and the unseen seeded fraction stay
within the base's spread by their means. One read moves on bgal: the share of the training rows the
route seeds, 0.396 against 0.424, two seeds of three below the base's lowest. It follows the threshold,
not the network's accuracy: the confidences of the training rows cluster just above 0.990 (median
0.9907), so a threshold that moves in its fourth decimal moves many rows. Seed 2 shows it: the same
held-out seeded share as the base (0.387) at a threshold of 0.99036 against 0.99014, and 0.377 of the
training rows seeded against 0.414. `T200`'s seed 3 moved the other way (threshold 0.99004 against
0.99066, 0.334 seeded against 0.305). The seeded rows' agreement is unchanged or better (truth in the
cone 0.997 against 0.996). For the maintainer's ruling: the threshold's sensitivity is the calibration's
(a quantile of 150 held-out rows), and every speed change that regroups sums will move it this much.

Per seed (phase 2):

| set | seed | learning wall (s) | widths trained | kept | best / epochs | median error (deg) | within 6 deg | in cone | seeded (bgal: training / unseen) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `T200` | 1 | 26 | 16, 32 | 16 | 3 / 24 | 4.08 | 0.690 | 0.962 | 0.117 |
| `T200` | 2 | 28 | 16, 32 | 16 | 5 / 26 | 3.79 | 0.747 | 0.970 | 0.193 |
| `T200` | 3 | 29 | 16, 32 | 16 | 5 / 26 | 3.51 | 0.792 | 0.973 | 0.334 |
| `T1000` | 1 | 69 | 16, 32 | 16 | 3 / 18 | 3.78 | 0.756 | 0.974 | 0.223 |
| `T1000` | 2 | 179 | 16, 32, 64 | 32 | 4 / 19 | 3.51 | 0.813 | 0.984 | 0.243 |
| `T1000` | 3 | 71 | 16, 32 | 16 | 5 / 20 | 3.46 | 0.823 | 0.987 | 0.262 |
| `T5000` | 1 | 309 | 16, 32 | 16 | 6 / 20 | 3.47 | 0.836 | 0.978 | 0.370 |
| `T5000` | 2 | 786 | 16, 32, 64 | 64 | 4 / 18 | 3.13 | 0.913 | 0.992 | 0.508 |
| `T5000` | 3 | 304 | 16, 32 | 16 | 7 / 21 | 3.21 | 0.892 | 0.989 | 0.451 |
| `T` | 1 | 574 | 16, 32 | 16 | 4 / 7 | 3.09 | 0.931 | 0.996 | 0.464 |
| `T` | 2 | 621 | 16, 32 | 16 | 6 / 9 | 3.18 | 0.914 | 0.991 | 0.548 |
| `T` | 3 | 598 | 16, 32 | 16 | 5 / 8 | 3.01 | 0.946 | 0.996 | 0.583 |
| bgal | 1 | 118 | 16, 32 | 16 | 5 / 20 | 2.51 | 0.871 | 0.960 | 0.381 / 0.041 |
| bgal | 2 | 313 | 16, 32, 64 | 32 | 5 / 20 | 2.52 | 0.871 | 0.957 | 0.377 / 0.063 |
| bgal | 3 | 293 | 16, 32, 64 | 32 | 5 / 20 | 2.47 | 0.877 | 0.966 | 0.431 / 0.057 |

**Files changed in this phase:** `src/main/learn/simple_hhcnn.f90`, `simple_hhcnn_tester.f90`,
`simple_hhcnn_bench.f90`, this note.

**Exit items.** The step identical at one thread and at 24 (the kernel benchmark) and at sixteen (the
HHCNN tester, repeated) on both builds: above. The gain per set and kernel: above (the subsets and bgal
gain 1.5 to 3.3 times; at a batch of 256 the step gains 2 to 5% and the sets about 1%). The reads within
the corrected spread: yes by their means, with bgal's training-row seeded fraction flagged above.
The variant build's benchmark also shows what phase 5 inherits: its explicit loops make the blocks
24 to 39 times slower than `sgemm` (12 rings: a one-thread step of 47 s against 2.3 s at 32 channels,
217 s against 6.1 s at 64; the blocks' forward pass per sample 271 ms against 7 ms at 64 channels).

### Phase 3: the readout (2026-11-10, done)

**What changed** (`src/main/learn/simple_hhcnn.f90`, its tester and the kernel benchmark). The readout
scores every direction class against the sample at every in-plane rotation of the polar grid (and
against the sample's reflection) and keeps the best; it was the larger part of the forward pass at 16
and 32 channels.

- *The class products as matrix products.* Per angular order the sample's spectrum enters as four
  columns (the spectrum, the spectrum turned by a quarter, and the two of its reflection, the conjugate
  spectrum), and one product with the templates of all classes gives every class's real and imaginary
  products for the map and its reflection (`readout_products`).
- *The rotation scores as one matrix product.* Every class's score at every rotation is the product of
  a fixed rotation basis (the grid's `cos(m theta)` and `sin(m theta)`, scaled, held by the network) with
  the classes' coefficients (`readout_rotations`); the maximum over the rotations and the reflection is
  then a scan, the map's columns first as the ties went before.
- *The backward through the winning rotation likewise.* Each class's gradient at its maximising rotation
  becomes four coefficients per order; the spectrum's gradient is one product of the templates with them
  per order. The templates' gradient leaves the per-chunk buffers: each sample keeps its spectrum columns
  and coefficients by its place in the batch, and after the batch one product per order and block of 64
  classes (`CNN_TCOLS`) sums them over the batch (`template_grads`), each block one call on fixed shapes,
  so the sums are the same at any team. The chunk buffers no longer hold a template gradient (1.09
  million weights at 500 classes, 4.4 MB per chunk, 86 chunks at a batch of 256).
- All products go through the build's rule of the first finding: `sgemm` with a team-safe BLAS, the
  explicit loops otherwise. The accumulations are now single precision (they were double in the loops);
  the scores and gradients agree to single-precision rounding.

**The kernels** (`test=pose_student_bench`, 24 threads; the one-thread figures on CPU 0; phase 2 against
phase 3, `scratch/keep/bench/kernels_tree_p2_3.txt`, `kernels_tree_p3b_{1,2,3}.txt`, which agree within 1%):

| rings | width | readout per sample (ms), phase 2 | phase 3 | backward (ms), phase 2 | phase 3 | step of 256 at one thread (ms), phase 2 | phase 3 | at 24 threads, phase 2 | phase 3 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 12 | 16 | 2.5 | 0.4 | 1.5 | 1.2 | 1,357 | 643 | 131 | 92 |
| 12 | 32 | 2.5 | 0.5 | 4.0 | 3.5 | 2,358 | 1,664 | 248 | 208 |
| 12 | 64 | 2.7 | 0.6 | 14.2 | 13.3 | 6,170 | 5,447 | 637 | 555 |
| 18 | 16 | 2.8 | 0.5 | 2.5 | 1.7 | 1,712 | 894 | 176 | 132 |
| 18 | 32 | 2.9 | 0.6 | 6.4 | 5.2 | 3,183 | 2,370 | 348 | 310 |
| 18 | 64 | 3.2 | 0.8 | 20.9 | 19.7 | 8,733 | 7,877 | 900 | 820 |
| 25 | 16 | 3.2 | 0.6 | 3.8 | 2.6 | 2,221 | 1,232 | 233 | 192 |
| 25 | 32 | 3.3 | 0.7 | 9.8 | 8.2 | 4,405 | 3,448 | 502 | 454 |
| 25 | 64 | 3.7 | 1.1 | 30.0 | 28.7 | 12,325 | 11,353 | 1,321 | 1,212 |
| 37 | 16 | 3.5 | 0.6 | 5.7 | 4.0 | 2,925 | 1,803 | 346 | 296 |
| 37 | 32 | 3.7 | 0.8 | 14.8 | 13.1 | 6,401 | 5,271 | 766 | 718 |
| 37 | 64 | 4.2 | 1.3 | 43.9 | 42.5 | 18,065 | 16,951 | 2,035 | 1,875 |

The class products now run at 17 to 35 GFLOP/s and the rotation scores at 70 to 78 (3.6 and 9.6
before). The step at 24 threads, summed over the twelve shapes, went from 7,641 ms to 7,121 ms with the
matrix products and to 6,869 ms (three runs, 6,863 to 6,874) once the templates' gradient left the
chunk buffers: 10% in all, 20% at 16 channels. At 64 channels the blocks' products, already at 75 to
90% of a core's peak, are nine tenths of the step (phase 5).

**The coarse-to-fine rotation search, measured and dropped.** A variant (`scratch/keep/variants/c2f_patch.py`,
`c2f_speed.patch`) scored every class on a coarse grid of 40 rotations by a matrix product and refined
the best two coarse cells on the full grid. *Agreement* (a probe variant, `c2f_probe.patch`, which kept
the full grid and ran the coarse search on its scores alongside; models trained by the tree on
`T5000` and bgal, seed 1, predicting real particles: the test half, 27,072 forward passes with the shift
check's, and bgal's training and all-particle projects, 4,572 and 5,828; `scratch/keep/logs/chain_p3a.log`):
with one refined cell the maximising rotation or reflection differs for 1.8 to 2.6% of the class
searches and the most probable class changes in 0.01 to 0.03% of the samples; with two cells for 0.14 to
0.22% of the class searches and never for the most probable class or its rotation. *Speed* (three
alternating benchmark runs each): the step at 24 threads 1.7% faster (5.6% at 16 channels), the
one-thread readout 8% slower (its scalar refinement costs more than the product it replaces); not worth
an approximation, so the full grid stays. *Skipping the backward over negligible classes* was not built:
with the readout's backward as matrix products it is about 0.1 ms of a 2.4 ms sample at 16 channels, so
even skipping every class would gain under 5%, below the plan's bar of 10%.

**Identity and tests.** The HHCNN tester passes on both builds: the gradient check (359 entries, none off
at 1e-2), the Hermitian contract, the rotation and mirror tests, and the thread identity at four and, by
hand with `OMP_THREAD_LIMIT=24`, at sixteen threads, three times each (34 of 34;
`scratch/keep/logs/p3_hhcnn_{tree,noteam}_t16_rep{1,2,3}.log`); the pose student's suite (64 of 64); the
kernel benchmark's step identical at 24 threads for every shape on the default build; the fast gate on
both builds (`p3_fast_gate_{tree,noteam}.log`).

**The gain per set and the reads** (mean [minimum, maximum] over seeds 1 to 3; `scratch/keep/summary.tsv`,
build `tree`, phase 3):

| set | base (s) | phase 2 (s) | phase 3 (s) | gain on the base | 16-channel design, phase 2 / phase 3 (s) | within 6 deg | in cone | seeded | base: within 6 deg / in cone / seeded |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `T200` | 89 [79, 100] | 27 | 20 [19, 21] | 4.4 | 6.1-7.2 / 4.4-5.0 | 0.743 [0.690, 0.793] | 0.968 [0.962, 0.973] | 0.206 [0.117, 0.307] | 0.743 / 0.968 / 0.204 |
| `T1000` | 192 [118, 333] | 106 | 90 [55, 158] | 2.1 | 19.9-21.9 / 14.2-15.8 | 0.797 [0.757, 0.823] | 0.981 [0.974, 0.986] | 0.240 [0.219, 0.258] | 0.797 / 0.982 / 0.241 |
| bgal | 351 [169, 457] | 241 | 217 [100, 285] | 1.6 | 34.3-35.9 / 26.8-26.9 | 0.873 [0.871, 0.877] | 0.961 [0.956, 0.965] | 0.407 [0.371, 0.431] | 0.874 / 0.960 / 0.424 |
| `T5000` | 468 [310, 784] | 466 | 428 [269, 739] | 1.1 | 90.3-99.5 / 73.9-81.5 | 0.880 [0.836, 0.913] | 0.986 [0.979, 0.992] | 0.441 [0.367, 0.506] | 0.881 / 0.986 / 0.443 |
| `T` | 606 [586, 630] | 597 | 538 [520, 559] | 1.1 | 158.5-203.8 / 127.5-166.5 | 0.930 [0.913, 0.946] | 0.995 [0.991, 0.996] | 0.531 [0.470, 0.585] | 0.930 / 0.995 / 0.532 |

bgal's unseen particles seeded 0.057 [0.037, 0.078] (base 0.058 [0.045, 0.071]); held-out recovery 1.000
on bgal, `T200` and `T1000`, 0.997 on `T5000`, 0.999 on `T`, as at the base. The reads stay within the
base's spread by their means; bgal's training-row seeded fraction, flagged in phase 2, is back near the
base's (0.407 against 0.424; its seeds move with the threshold, 0.371 to 0.431). The 16-channel design is
17 to 28% faster than in phase 2 (28% on the subsets, 17 to 19% at a batch of 256); the rows (the teacher's table fill, out of scope) are now
a larger share on `T` (56 s of 538).

Per seed (phase 3):

| set | seed | learning wall (s) | widths trained | kept | best / epochs | median error (deg) | within 6 deg | in cone | seeded (bgal: training / unseen) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `T200` | 1 | 19 | 16, 32 | 16 | 3 / 24 | 4.09 | 0.690 | 0.962 | 0.117 |
| `T200` | 2 | 21 | 16, 32 | 16 | 5 / 26 | 3.79 | 0.748 | 0.970 | 0.193 |
| `T200` | 3 | 20 | 16, 32 | 16 | 5 / 26 | 3.51 | 0.793 | 0.973 | 0.306 |
| `T1000` | 1 | 54 | 16, 32 | 16 | 3 / 18 | 3.78 | 0.756 | 0.974 | 0.219 |
| `T1000` | 2 | 158 | 16, 32, 64 | 32 | 4 / 19 | 3.52 | 0.812 | 0.984 | 0.243 |
| `T1000` | 3 | 57 | 16, 32 | 16 | 5 / 20 | 3.46 | 0.823 | 0.986 | 0.258 |
| `T5000` | 1 | 275 | 16, 32 | 16 | 6 / 20 | 3.47 | 0.836 | 0.979 | 0.367 |
| `T5000` | 2 | 738 | 16, 32, 64 | 64 | 4 / 18 | 3.14 | 0.913 | 0.992 | 0.506 |
| `T5000` | 3 | 269 | 16, 32 | 16 | 7 / 21 | 3.21 | 0.891 | 0.989 | 0.450 |
| `T` | 1 | 520 | 16, 32 | 16 | 4 / 7 | 3.09 | 0.931 | 0.996 | 0.469 |
| `T` | 2 | 558 | 16, 32 | 16 | 6 / 9 | 3.18 | 0.913 | 0.991 | 0.539 |
| `T` | 3 | 535 | 16, 32 | 16 | 5 / 8 | 3.01 | 0.946 | 0.996 | 0.585 |
| bgal | 1 | 100 | 16, 32 | 16 | 5 / 20 | 2.51 | 0.871 | 0.961 | 0.371 / 0.037 |
| bgal | 2 | 285 | 16, 32, 64 | 32 | 5 / 20 | 2.52 | 0.871 | 0.956 | 0.431 / 0.078 |
| bgal | 3 | 266 | 16, 32, 64 | 32 | 5 / 20 | 2.47 | 0.877 | 0.965 | 0.419 / 0.054 |

**Files changed in this phase:** `src/main/learn/simple_hhcnn.f90`, `simple_hhcnn_tester.f90`,
`simple_hhcnn_bench.f90`, this note.

**Exit items.** The HHCNN tester passes on both builds (gradient check, Hermitian contract, rotation and
mirror tests, thread identity): above. The coarse search's agreement with the full grid measured on
real tensors, per class and for the mode class: above (dropped for speed). The gain: above, on the
kernels and per set. The reads within the corrected spread: yes by their means.

### Phase 4: the schedule (2026-11-10, done)

**Kept.** *The training in epochs* (`simple_pose_student.f90`, `simple_strategy_pose_student.f90`): at most
30 epochs (`MAX_EPOCHS`, the base's 2,500 steps on the full training half), stopped after 3 epochs
without a held-out improvement (`PATIENCE_EPOCHS`). The base counted both in optimiser steps (2,500 and
250), which on the subsets and bgal, with 12 to 18 batches per epoch, meant 14 to 21 epochs of patience
after the best one. In the phase 3 runs (34 designs over the fixed sets and seeds) the best held-out
epoch came at 3 to 7 and its improvements were never more than 2 epochs apart, so a patience of 3 epochs
keeps the same best epoch everywhere: the fifteen runs of this schedule kept exactly the models of
phase 3 (their held-out tables identical, file for file), at a fraction of the time. *The widest design
only on a large set* (the design choice's size rule, `WIDEST_MIN_ROWS = 10,000` training rows): the
search tries 16, then 32 channels, and 64 only from 10,000 training rows; below, the 64-channel design
won one of the four searches that trained it on the fixed sets and costs two and a half times the
32-channel one. The rate, the decay and the batch stay as they were.

**The candidates** (each on every fixed set and seed, three seeds; mean training wall time and reads;
the schedule candidates with the full design search, the design candidates with the kept schedule):

| candidate | `T200` (s) | `T1000` (s) | bgal (s) | `T5000` (s) | `T` (s) | kept widths `T1000` / bgal / `T5000` | reads against the base |
| --- | --- | --- | --- | --- | --- | --- | --- |
| phase 3 (patience in steps) | 20 | 90 | 217 | 428 | 538 | 16,32,16 / 16,32,32 / 16,64,16 | the base's spread |
| S1: patience of 3 epochs, the full search | 8.0 | 36.5 | 87.8 | 184.1 | 538.5 | the same | identical models to phase 3 |
| S2: S1 with a one-epoch warm-up and the rate halved at every epoch without improvement | 8.5 | 40.3 | 122.8 | 275.4 | 927.9 | 16,32,16 / 32,64,32 / 32,64,16 | wider designs kept more often, no read gained; bgal's truth in cone 0.955, below the base's lowest |
| S3: S1 with twice the batch (an eighth of the rows within 32 and 512), the rate times sqrt(2) | 8.0 | 25.5 | 71.4 | 307.1 | 902.1 | 16,16,16 / 16,64,16 / 32,32,16 | worse: bgal within 6 deg 0.856 and unseen seeded 0.020 (base 0.871-0.877 and 0.058), held-out recovery 0.990 |
| D2: S1 with 64 channels only from 10,000 rows (**kept**) | 8.5 | 23.3 | 44.6 | 126.6 | 537.8 | 16,32,16 / 16,32,32 / 16,32,16 | the base's spread; the same models as S1 on 14 of 15 runs |
| D1: S1 with successive halving (every width 3 epochs, the leader by held-out loss completed, set against a completed 16-channel design by the recovery test) | 19.0 | 58.6 | 97.6 | 248.2 | 1,206.6 | 16,64,64 / 64,64,16 / 16,64,16 | the 64-channel design always leads after 3 epochs; bgal within 6 deg 0.865, below the base's lowest |

Per seed, the training wall time (s) / kept width / within 6 deg / truth in the cone (bgal: the training
rows against the teacher):

| set | seed | S1 | S2 | S3 | D2 (kept) | D1 |
| --- | --- | --- | --- | --- | --- | --- |
| `T200` | 1 | 7 / 16 / 0.690 / 0.962 | 7 / 16 / 0.667 / 0.952 | 7 / 16 / 0.648 / 0.939 | 8 / 16 / 0.690 / 0.962 | 19 / 64 / 0.759 / 0.953 |
| `T200` | 2 | 8 / 16 / 0.748 / 0.970 | 10 / 16 / 0.753 / 0.970 | 8 / 16 / 0.720 / 0.949 | 9 / 16 / 0.748 / 0.970 | 18 / 16 / 0.748 / 0.970 |
| `T200` | 3 | 8 / 16 / 0.793 / 0.973 | 8 / 16 / 0.787 / 0.976 | 9 / 16 / 0.779 / 0.977 | 8 / 16 / 0.793 / 0.973 | 20 / 16 / 0.793 / 0.973 |
| `T1000` | 1 | 22 / 16 / 0.756 / 0.974 | 24 / 16 / 0.740 / 0.971 | 28 / 16 / 0.730 / 0.952 | 22 / 16 / 0.756 / 0.974 | 56 / 16 / 0.756 / 0.974 |
| `T1000` | 2 | 64 / 32 / 0.812 / 0.984 | 71 / 32 / 0.796 / 0.976 | 26 / 16 / 0.763 / 0.970 | 23 / 32 / 0.812 / 0.984 | 57 / 64 / 0.825 / 0.984 |
| `T1000` | 3 | 24 / 16 / 0.823 / 0.986 | 26 / 16 / 0.818 / 0.986 | 23 / 16 / 0.810 / 0.978 | 24 / 16 / 0.823 / 0.986 | 63 / 64 / 0.804 / 0.968 |
| bgal | 1 | 43 / 16 / 0.871 / 0.961 | 130 / 32 / 0.893 / 0.949 | 40 / 16 / 0.852 / 0.953 | 45 / 16 / 0.871 / 0.961 | 104 / 64 / 0.879 / 0.960 |
| bgal | 2 | 119 / 32 / 0.871 / 0.956 | 119 / 64 / 0.879 / 0.959 | 126 / 64 / 0.877 / 0.959 | 44 / 32 / 0.871 / 0.956 | 104 / 64 / 0.875 / 0.960 |
| bgal | 3 | 102 / 32 / 0.877 / 0.965 | 119 / 32 / 0.865 / 0.955 | 48 / 16 / 0.840 / 0.947 | 45 / 32 / 0.877 / 0.965 | 85 / 16 / 0.841 / 0.947 |
| `T5000` | 1 | 118 / 16 / 0.836 / 0.979 | 350 / 32 / 0.867 / 0.985 | 371 / 32 / 0.852 / 0.977 | 126 / 16 / 0.836 / 0.979 | 247 / 16 / 0.836 / 0.979 |
| `T5000` | 2 | 313 / 64 / 0.913 / 0.992 | 341 / 64 / 0.903 / 0.987 | 399 / 32 / 0.883 / 0.984 | 132 / 32 / 0.883 / 0.985 | 245 / 64 / 0.913 / 0.992 |
| `T5000` | 3 | 121 / 16 / 0.891 / 0.989 | 135 / 16 / 0.881 / 0.984 | 152 / 16 / 0.871 / 0.977 | 121 / 16 / 0.891 / 0.989 | 252 / 16 / 0.891 / 0.989 |
| `T` | 1 | 518 / 16 / 0.931 / 0.996 | 583 / 16 / 0.921 / 0.996 | 1,547 / 32 / 0.911 / 0.988 | 518 / 16 / 0.931 / 0.996 | 1,190 / 16 / 0.931 / 0.996 |
| `T` | 2 | 557 / 16 / 0.913 / 0.991 | 601 / 16 / 0.913 / 0.993 | 580 / 16 / 0.924 / 0.996 | 558 / 16 / 0.913 / 0.991 | 1,227 / 16 / 0.913 / 0.991 |
| `T` | 3 | 540 / 16 / 0.946 / 0.996 | 1,600 / 32 / 0.936 / 0.994 | 579 / 16 / 0.950 / 0.996 | 537 / 16 / 0.946 / 0.996 | 1,203 / 16 / 0.946 / 0.996 |

**The kept schedule's reads** (D2; mean [minimum, maximum] over seeds; the base's in brackets after):
`T200` within 6 deg 0.743 [0.690, 0.793] (base 0.743), truth in cone 0.968 (0.968), seeded 0.206 (0.204);
`T1000` 0.797 [0.757, 0.823] (0.797), 0.981 (0.982), 0.240 (0.241); bgal 0.873 [0.871, 0.877] (0.874),
0.961 (0.960), training rows seeded 0.407 (0.424), unseen seeded 0.057 [0.037, 0.078] (0.058); `T5000`
0.870 [0.836, 0.891] (0.881 [0.837, 0.914]), 0.984 [0.979, 0.989] (0.986 [0.979, 0.992]), 0.434
(0.443); `T` 0.930 (0.930), 0.995 (0.995), 0.531 (0.532). Held-out recovery as at the base. The one
model that differs from the full search's is `T5000` seed 2, where the size rule keeps 32 channels in
place of the 64 the search kept (0.883 within 6 deg against 0.913), inside the base's spread. The
kept tree's runs are D2's: the variant was built from the same source with the same rule, and the tree
reproduces its held-out tables on `T200` and bgal at seed 1 (`scratch/keep/runs/tree_p4fin/`).

**Where the time now goes** (the kept schedule's learning bench records): on `T` the rows (the teacher's
table fill, out of scope) are 56 s of 538 and the two designs 7 to 9 epochs each; on bgal the 16- and
32-channel designs 8 epochs each, against 20 before.

**Tests.** The pose student's suite (64 of 64) and the strategy's (41 of 41) and the fast gate pass on
both builds (`scratch/keep/logs/p4_{pose_student,strategy,fast_gate}_{tree,noteam}.log`). The registered
workflow test `pose_student` passes through CTest at its eight threads in 38.7 s, its training 21.8 s
(the base 313.9 s and 296 s; `p4_ctest_pose_student.log`, `p4_pose_student_metrics.tsv`).

**Files changed in this phase:** `src/main/learn/simple_pose_student.f90` (the budget and patience in
epochs), `src/main/strategies/parallelization/simple_strategy_pose_student.f90` (`MAX_EPOCHS`,
`PATIENCE_EPOCHS`, `WIDEST_MIN_ROWS`), `doc/algorithms/pose_student.md` (the training section), this
note. The variants measured and dropped are kept as patches and scripts in `scratch/keep/variants/`
(`sched_patch.py`, `sched_s2.patch`, `sched_s3.patch`, `design_patch.py`, `design_d1.patch`,
`design_d2.patch`); the runs' logs in `scratch/keep/runs/{tree_p4s1,sched_s2_p4s2,sched_s3_p4s3,design_d2_p4d2,design_d1_p4d1}/`.

**Exit items.** Every candidate's time and reads on every set and seed: above. One schedule (S1) and one
design choice (D2) kept. The pose student's tester passes.

### Phase 5: the blocks and the rest (2026-11-10, done)

**The profile after phases 2 to 4** (the kernel bench, `scratch/keep/bench/kernels_tree_p5prof_1.txt`
and `kernels_tree_p3b_{1,2,3}.txt`, per sample at one thread; the learning bench records of phase 4's
runs). Within a design the gradients (the parallel region) are about 90% of the time, the held-out
evaluation 4 to 5%, the reduction and Adam about 1%, the batch assembly under 1% (`asm` 0.0 to 0.3 s
per design). Within a sample:

| width | the blocks' products (`sgemm`) | the blocks' other work (im2col, gate, pooling, their backward, col2im, the first block's normalisation) | the readout |
| --- | --- | --- | --- |
| 16 | 50-53% | 33-37% | 9-17% of the forward |
| 32 | 66-70% | 22-30% | 6-10% |
| 64 | 78-81% | 16-19% | 3-5% |

The blocks' other work held more than 10% of training at every width; its largest kernels at 25
rings and 16 channels were the gate's backward (469 us per sample), the gate's forward (253) and im2col
(239). The products run at 130 to 150 GFLOP/s at one thread, near the core's single-precision peak
with AVX2. The team-safe build's step scaled 7.0 times at 16 channels over 24 threads (649 ms at one
thread, 92 ms at the team). *The build without a team-safe BLAS* (`scratch/keep/build_noteam`, PS01's
explicit loops) trained 14 times slower: `T200` seed 1 in 121.0 s against 8.5 s (the 32-channel design
47.3 s against 2.6), bgal about 53 s per epoch of the 16-channel design against 2.5 (stopped at its
fifth epoch; `scratch/keep/runs/noteam_p5prof/`). That is above section 4.5's factor of three.

**Kept, in the order made** (`simple_hhcnn.f90`):

1. *The explicit loops of a build without a team-safe BLAS* run down columns: `gemm_nn` and
   `gemm_nt_acc` loop over the output's columns and the inner terms, with the column the inner vector
   loop. `gemm_tn` transposes `a` in blocks of `CNN_LOOP_ROWS = 64` rows of the product (64 kB on the
   stack) so that it loops the same way. The HHCNN suite on that build ran 26 s before and 7.8 s after.
2. *The gate without branches or divisions per element*: the forward is out of place (`pre` to
   `gated`, which removes a copy) with the factor `max(|z| + b, 0) / |z|`. The backward uses `1/|z|`
   once and a 0/1 mask for the rectifier. The first block's normalisation multiplies by a per-ring
   scale (`nscale`, the filter norm over the normaliser) in place of a division per element. At 25
   rings and 16 channels, at one thread: gate forward 253 to 58 us, gate backward 469 to 77 us, the
   forward's remainder 214 to 79 us. The one-thread step went 7 to 12% faster (sum over the bench's
   shapes 59,169 to 54,449 ms, `kernels_tree_p5_{1..4}.txt`), but the team step barely moved (6,872 to
   6,710 ms, -2%). At the team the step was limited by something other than arithmetic.
3. *Blocks 2 and 3 build im2col in tiles* of `CNN_TILE = 128` columns. The forward and the backward
   each build the tile, the backward again rather than storing the whole matrix. The weights' gradient
   and the input's gradient run per tile, with col2im adding onto the input's gradient. The block's
   im2col matrix (about 2.2 MB per sample at 32 channels) no longer passes through memory three times
   per sample: a tile is 0.6 MB at most and stays in a core's cache. Two observations showed the team
   step was memory-bound at small widths. Eight performance cores ran the 16-channel step at 96 ms
   against 92 ms at all 24 cores (`kernels_tree_p5_t8.txt`), while the 64-channel step gained 35% from
   the efficiency cores. And change 2 removed arithmetic without moving the team step. The tiled
   variant, built under scratch and measured in three alternating pairs against the tree: team step
   sum 5,070 / 5,079 / 5,098 ms against 6,723 / 6,717 / 6,710 (-24.4%). By width: 16 channels -40%,
   32 -33%, 64 -19%. At one thread -2%. Tiles of 64 and 256 columns ran within 2% of 128 (5,039 and
   4,949 ms, one run each, against 128's four runs of 5,063 to 5,098), so 128 is kept
   (`kernels_tile{64,128,256}_*.txt`). The first block keeps its whole matrix: it is reused by the
   normalisation and is only nine rows deep.

Not changed: the products' sizes (several samples in one product: the products already run near
peak at one thread, and samples are what the team divides), `!$omp simd` (not tried: the
branch-free forms gave the gains above), the feature tensor and the batch assembly (under 1%), the reduction and Adam
(about 1%). The readout stays as phase 3 left it, although it is the largest share left outside the
blocks: 14 to 18% of a one-thread sample at 16 channels, under 10% at 32 and above.

**The kept tree's kernels** (`kernels_tree_p5final_1.txt`, load 0.8 at the start; sums over the
bench's 12 shapes): team step 4,918 ms against 6,872 after phase 4 (-28%; by width: 16 channels 429
against 722, 32 1,077 against 1,692, 64 3,412 against 4,458). One-thread step 53,546 against 59,169
(-9.5%). The 16-channel step now scales 10.2 times over the team (583 / 57 ms at 12 rings), against
7.0. The training step is identical at one thread and at the team on all 12 shapes. (The bench's
im2col and col2im columns time one tile in blocks 2 and 3; the rebuilt tiles fall in the remainder
columns.)

**Training on the fixed sets** (the tree, three seeds, `scratch/keep/tools/harness.sh tree ... 5
<seed>`, 24 threads, load 15 to 24; mean training wall time [minimum, maximum] in seconds, against
phase 4's kept D2):

| set | phase 4 (D2) | phase 5 | gain | kept widths |
| --- | --- | --- | --- | --- |
| `T200` | 8.5 [8.0, 9.1] | 8.5 [8.2, 8.8] | 0% (fixed costs) | 16,16,16 (the same) |
| `T1000` | 23.3 [22.3, 24.4] | 18.3 [17.9, 18.8] | -21% | 16,32,16 (the same) |
| bgal | 44.6 [44.1, 44.9] | 34.7 [33.9, 35.8] | -22% | 16,32,32 (the same) |
| `T5000` | 126.6 [120.8, 132.4] | 86.1 [82.5, 90.3] | -32% | 16,32,16 (the same) |
| `T` | 537.8 [518.0, 558.3] | 369.3 [359.0, 379.8] | -31% | 16,16,16 (the same) |

On `T` the designs take 299 s of the 369 (468 of 538 before); the rest is mostly the teacher's rows
(56 s). The phase's gains come from the tiling at the team. The gate's gain shows only at one thread
and does not separate from the tiling in these runs.

**The reads** (`scratch/logs/p5_spread_tree5.txt`). Every set kept the same widths as phase 4 on every
seed. Every read differs from phase 4's in the third or fourth decimal, the size of the rounding of a
reordered sum: within 6 deg `T` 0.930 (0.930), `T5000` 0.870 (0.870), bgal 0.872 (0.873); truth in
the cone bgal 0.961 (0.961); held-out recovery 0.997 to 1.000. Against the corrected spread (phase 1's
three seeds) every read's mean lies inside the range except two of bgal's. Training rows seeded: 0.414
against the range 0.415 to 0.437 (phase 4's D2 0.407; the seed range 0.386 to 0.435 overlaps the
corrected one). Seeded rows in the cone: 0.998 against 0.995 to 0.997, above the range.

**The build without a team-safe BLAS after the phase** (`scratch/keep/runs/noteam_p5prof2/`, seed 1,
load 20): `T200` 13.0 s (121.0 before; the tree 8.2), bgal 60.6 s (the tree 35.8). Its 32-channel
design on bgal runs 4.3 s per epoch against the tree's 2.2. It now trains 1.6 to 1.75 times slower
than the team-safe build, within section 4.5's factor of three.

**Tests.** `test_tiled_blocks` is new in the HHCNN suite. It runs the network at 48 angles, where block
2 has 384 columns (three tiles, the real/imaginary boundary inside the second) and block 3 has 192 (a
tile and a part tile). The responses, the weights' gradient and the input's gradient must equal those
of the whole im2col matrix within 1e-5 of the largest entry. The suite's other shapes fit in one tile.
A mutation that shifts only the tiles after the first fails these four checks and no other
(checked by hand, then restored). The HHCNN suite passes 41 of 41 on both builds, by default and three
times under `OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=24`; the teams of 4 and 16 are bit-identical to one
thread each time (`scratch/keep/logs/p5_hhcnn_{tree,noteam}*.log`). Also passing on both builds: the
pose student's suite (64 of 64), the strategy's (41 of 41) and the fast gate
(`p5_{pose_student,pose_student_strategy}_{tree,noteam}.log`, `p5_gate_{tree,noteam}.log`). The
registered workflow test `pose_student` passes through CTest in 31.5 s (38.7 after phase 4;
`p5_ctest_pose_student.log`).

**Files changed in this phase:** `src/main/learn/simple_hhcnn.f90` (the loops, the gate, the
normalisation scale, the tiles), `src/main/learn/simple_hhcnn_tester.f90` (`test_tiled_blocks`),
`src/main/learn/simple_hhcnn_bench.f90` (the table of the blocks' other kernels, the tiled products),
this note. The tile variants are kept as `scratch/keep/variants/tile_patch.py` and `tile128.patch`,
and the bench records in `scratch/keep/bench/`.

**Exit items.** The profile after phases 2 to 4 on both builds: above. Each change's gain: the loops
(the build without a team-safe BLAS from 14 to under 1.75 times the team-safe build), the gate and
normalisation (-7 to -12% at one thread), the tiles (-24% of the team step, -21 to -32% of training
from `T1000` up). The reads lie inside the corrected spread, apart from bgal's seeded fraction, 0.0015
below its lowest (phase 4: 0.008 below), and its seeded-in-cone fraction, above its highest.

### Phase 6: the end-to-end validation and the report (2026-11-10, done)

**The end-to-end validation.** `validate_pose_student` ran on bgal four times, alternating the run's
build and the base build (`prg=validate_pose_student projfile=bgal_addon.simple nsample=1500
mskdiam=180 pgrp=d2 pgrp_start=c1 nthr=24`, each in its own directory with a copy of
`/data2/bgal_addon/3_selection/bgal_addon.simple`, load 7 to 15). The program runs solve3D on 1,500
class-balanced training particles as the teacher, trains the student on its poses, predicts all
2,756 selected particles and reconstructs maps from the student's poses. Its solve3D is unseeded, so
every run has its own teacher.

| run | solve3D (s) | training (s) | widths trained; kept | training rows seeded | within 6° of the teacher | unseen seeded | unseen map at the student's poses, FSC 0.143 against the teacher's (Å) | teacher's map between halves, 0.143 (Å) | total (s) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| base 1 | 321.2 | 474.8 | 16, 32, 64; 32 | 0.412 | 0.859 | 0.044 | 4.80 | 4.08 | 855 |
| run 1 | 655.3 | 50.8 | 16, 32; 32 | 0.545 | 0.843 | 0.143 | 6.40 | 4.47 | 761 |
| run 2 | 288.3 | 36.4 | 16, 32; 32 | 0.484 | 0.890 | 0.048 | 4.87 | 4.08 | 381 |
| base 2 | 320.6 | 185.9 | 16, 32; 16 | 0.625 | 0.915 | 0.043 | 5.63 | 4.08 | 564 |

The learning phase is 7.6 times faster on average (44 s against 330 s; the plan's baseline runs 2
and 5 took 178 and 537 s). Every held-out recovery is 1.000 at a threshold of 0.990. Run 2's reads
lie between the two base runs', apart from its unseen seeded fraction (0.048 against 0.043 and
0.044). Run 1 had the weakest teacher: its solve3D took twice as long, and its control map (the
training particles at the teacher's poses) reached 4.47 Å between halves against 4.08 in the other
runs. Against that teacher the student is below the base runs within 6° and in the unseen map,
above them in the unseen seeded fraction, and between them in the training rows seeded. With the teacher changing from run to run, two runs per build do
not separate the student from the teacher. The fixed-teacher reads of phases 0 to 5 are the
comparison. Evidence: `scratch/keep/validation_p6/val_{tree,base}_{1,2}/` (the logs,
`validate_time.log`, `1_validate_pose_student/validate_pose_student_report.txt` and the bench
records; the volumes and projects deleted), `scratch/keep/tools/val_stats.sh`,
`scratch/keep/logs/chain_p6v.log`.

**Tests on both builds** (`scratch/keep/logs/chain_p6t.log`):

- The workflow test `pose_student` through CTest, at its registered eight threads and 2,400 s
  deadline: 31.9 s on the default build and 51.5 s on the build without a team-safe BLAS, against
  307.2 s for the base build the same evening (313.9 s in phase 0;
  `p6_ctest_pose_student_{tree,noteam,base}.log`).
- The fast gate: 15 of 15 on both (`p6_gate_{tree,noteam}.log`).
- `unit_numerics`: 28 suites and 2,006 assertions on both, including the HHCNN suite (41), the pose
  student's (64) and the strategy's (41) (`p6_unit_numerics_{tree,noteam}.log`).

**The documents.**
- The algorithm page's training section: the schedule, the step (chunking, the parallel reduction,
  the tiled blocks), the readout's full-grid rotation search, the build's arithmetic rule; its
  benchmark paragraph names the new tables.
- The kernel benchmark's line in the test policy.
- The review file's resolution section, with a closing paragraph after the speed phases.
- The report beside this note, `pose_student_training_speed_report.md`.
- This note, moved to `completed/`, and its row in `doc/implementation_notes/README.md`, moved to the
  completed table.

**Files changed in this phase:** `doc/algorithms/pose_student.md`,
`doc/policies/test_environment_policy.md`,
`doc/implementation_notes/completed/learned_pose_state_assignment_review.md`,
`doc/implementation_notes/completed/pose_student_training_speed_report.md` (new),
`doc/implementation_notes/README.md`, this note (moved from `planned/` to `completed/`).

**Exit items.** All tests pass on both builds: the fast gate, `unit_numerics` and the workflow test
through CTest, above. The report is beside this note, and this note is moved.
