# Training speed of the pose student, and the review's corrections

Status: planned, 2026-11-09. Hans set the priority ("everything we can to speed up training of the
student; priority number one") and asked for an autonomous run; this note is its plan. The same day
a code review of the pose student (`../completed/learned_pose_state_assignment_review.md`, at
`519161ff3`) found seven open problems; Hans asked for them to be addressed in the same run. The
run is `student_speed` on the Dell (section 9). Nothing in it is committed; the run hands over one
working-tree diff.

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
  (its centred copy is `S_project/1_zero_project_shifts/M48noeo.simple`), the subsets under
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

(filled by the run, one entry per phase)
