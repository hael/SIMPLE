# Training speed of the pose student

Status: planned, 2026-11-09. Hans set the priority ("everything we can to speed up training of the
student; priority number one") and asked for an autonomous run; this note is its plan. The run is
`student_speed` on the Dell (section 9). Nothing in it is committed; the run hands over one
working-tree diff.

The pose student (`doc/algorithms/pose_student.md`) is a small network trained per dataset on the
posterior rows of the probabilistic search. It proposes a direction, a cone, an in-plane angle and a
mirror for particles the search has not seen. Its training is now the slowest part of every use: on
bgal, 1,500 training particles take 180 to 540 s to train against 4 to 15 s to predict 2,756
particles. This run makes training as fast as it can be made on this machine without changing what
the student learns beyond the run-to-run spread of its own reads.

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
  (`SIMPLE_SEED`); three seeds per set measure the spread of every read at the base. A change is
  judged against that spread, never against one pair of runs.
- **End to end.** `validate_pose_student` on bgal at the end (section 8, phase 5), with the timings
  table and the maps from the student's poses.

The training step must stay bit-identical between one thread and the full team (the network's
contract; the HHCNN tester and the kernel benchmark check it). Results may change against the base
where a change regroups floating-point sums; they are then judged on the reads.

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
- **The blocks** are already im2col plus `sgemm` (OpenBLAS on the calling thread); their arithmetic
  rate is unmeasured. They dominate at 64 channels.

## 4. The changes, in the order of the phases

### 4.1 Occupancy and the step (phase 1)

- The chunk size follows the batch alone, never the team (the step must not depend on the thread
  count): chunks of `max(1, ceiling(batch / CNN_NCHUNK))` samples with a named target of 48 to 96
  chunks, the target and the chunk buffers' memory measured (the buffers hold one gradient per
  chunk: at 920,000 weights and 64 chunks, 240 MB).
- The reduction in parallel over slices of the weights: each slice sums the chunks in chunk order,
  so the sums are bit-identical to the serial ones; the Adam update in the same pass.
- The chunk buffers zeroed in parallel, or not at all where the first sample writes them.

### 4.2 The readout (phase 2)

- **As matrix products.** The class products for all classes and both reflections as one complex
  product of the templates with the sample's harmonic spectrum, the rotation scores as one real
  matrix product with the rotation basis (`cos(m theta_r)`, `sin(m theta_r)`), both through
  `sgemm`; the backward through the winning rotation likewise.
- **Coarse to fine.** The score's order-16 polynomial sampled on a coarse grid of a named size (32
  to 64 rotations), the best few coarse cells refined on the full grid around them. Measured: the
  speed, and how often the coarse search's maximum differs from the full grid's on real tensors
  (the bgal and truth-set features), per class and for the mode class. Kept if the reads stay
  within the spread; the full-grid maximum stays the reference in the tester.
- **The backward over negligible classes.** Most classes have a probability and a target near zero;
  their template gradients are near zero. Skipping classes with `|p - t|` below a named bound is an
  approximation: measured for speed and reads, kept only with a gain above 10% and the reads within
  the spread.

### 4.3 The schedule (phase 3)

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

### 4.4 The blocks and the rest (phase 4)

What phase 0's arithmetic rates and phases 1 to 3's profiles leave as the largest share: the blocks'
`sgemm` sizes (several samples' positions in one product), im2col reuse between forward and
backward, the gating and pooling loops vectorised (`!$omp simd` where the compiler does not),
memory layout; the feature tensor and the batch assembly if they reach 5% of training. If nothing
remaining exceeds 10% of training, phase 4 records that and ends.

### 4.5 Not in this run

The teacher's rows (the table fill), the network's architecture and model file format, the
prediction contract (the route, the cone, the table's columns), the shift registration, warm starts
from an earlier model (the streaming use; recorded as outlook), a GPU, reduced precision.

## 5. Decisions already taken

- No new command-line key of any kind; every new number is a named constant with its reason in a
  comment.
- No second code path stays in the tree: alternatives are built and measured as variants in the
  run's scratch, and only the kept one is in the diff.
- The training step stays bit-identical between one thread and the full team.
- The student may change against the base; it is judged on the reads of section 1 against the
  base's seed spread, never on one pair of runs.

## 6. Risks

- A faster schedule that stops too early: the reads at `T200` and `T1000` are the first to show it.
- The coarse rotation search that misses a maximum: measured directly, and visible in the in-plane
  read-out of the validation.
- Load from the maintainer's work on the Dell distorting timings: every timing records `uptime`'s
  load and is repeated under load.

## 7. How the run is judged

The run reports, per fixed set and seed, the training wall time and the reads at the base and after
each phase, and the end-to-end validation. Hans rules on whether the reads are equal; the run sets
no floor. Its working target is training at least five times faster on every fixed set at equal
reads, the kernel benchmark's training samples per second at least three times the base's.

## 8. Phases

Detailed in section 9.2; in short: 0 baselines and the harness, 1 occupancy and the step, 2 the
readout, 3 the schedule, 4 the blocks and the rest, 5 the end-to-end validation and the report.

## 9. The run on the Dell

The development runs as the autonomous run `student_speed` on the Dell, with the engine of
`~/agent_runs/ms_regression/engine` and the spec in `~/tasks/student_speed_spec/` (run.conf,
prompt.md, phases.md with the launch procedure, rulings.md), the run root
`/data2/agent_runs/student_speed` behind the link `~/agent_runs/student_speed`. It is launched by
hand on the Dell after Hans has committed the tree (the validation program, the benchmarks, the
shift registration) and this note.

### 9.1 The machine and the data

**The Dell** (Oracle Linux 8.10, user elmlundho): Core Ultra 9 285, 24 cores (8 performance cores,
CPUs 0-7, and 16 efficiency cores; no hyperthreading), 62 GB; `/data2` about 900 GB free, `/home`
little: nothing of the run goes under `/home`. GCC 15.2 from gcc-toolset-15, CMake 3.26,
`python3.12` for CMake; `OPENBLAS_CORETYPE=Haswell` and `TMPDIR=$HOME/tmp` exported. The maintainer
uses the machine during the run: timings record the load and are repeated under load; his checkout
`~/src/SIMPLE` and its `build/` are read-only for the run.

**Builds:** `./compile_clean.sh` about 4 minutes with the machine free; increments by
`make -C build -j20`, `scripts/run_fast_gate.sh build`, `make -C build install`. The base HEAD is
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

### 9.2 The phases of the run

**Phase 0: baselines and the harness.**

1. Build the tree and the base export.
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
4. The base: three seeds on every set, the kernel benchmark twice; the spread of every read.

Exit: the base's times, reads and spread in Progress; the kernel table with the arithmetic rates;
the fast gate and `unit_numerics suite=hhcnn` pass.

**Phase 1: occupancy and the step** (section 4.1). Exit: the step identical at one thread and at
24 (the kernel benchmark and the HHCNN tester at 16 threads, repeated); the gain per set and kernel;
the reads within the spread.

**Phase 2: the readout** (section 4.2). Exit: the HHCNN tester passes (the gradient check, the
Hermitian contract, the rotation and mirror tests, the thread identity); the coarse search's
agreement with the full grid measured; the gain; the reads within the spread.

**Phase 3: the schedule** (section 4.3). Exit: every candidate's time and reads on every set and
seed in Progress; one schedule and one design choice kept; the pose student's tester passes.

**Phase 4: the blocks and the rest** (section 4.4). Exit: the profile after phases 1 to 3, each
change's gain, the reads within the spread; or the record that nothing exceeds 10% of training.

**Phase 5: the end-to-end validation and the report.**

1. `validate_pose_student` on bgal twice with the run's build and twice with the base build: the
   learning phase's time and split, the reads, the maps from the student's poses (the teacher
   differs between runs; the comparison is the learning time and whether the reads sit inside the
   runs' spread).
2. `test=pose_student` (the nightly route), the fast gate, `unit_numerics`.
3. The algorithm page's training section (the schedule, the chunking, the readout's search), the
   report beside this note (what was gained where, the evidence, the measurements, a commit
   message), this note moved to `completed/`, its row in `doc/implementation_notes/README.md`.

Exit: all tests pass; the report and the moved note.

### 9.3 The rules of the run

- Nothing is committed; git is read-only; the base is a `git archive` export.
- No new command-line key, public or internal; no second code path kept (section 5).
- The teacher's rows, the architecture, the model file format and the prediction contract are
  unchanged (section 4.5); a model trained by the base still predicts with the run's build.
- New or changed OpenMP code: a tester at 16 threads, repeated, compared bit for bit with one thread;
  no string temporaries inside parallel regions; work arrays sized by the team in use; the BLAS on
  the calling thread inside a region.
- No `goto`; THROW macros never continued with `//&`; no `.and.` short circuit on allocatable
  arrays; `use ..., only:`; no large type as a plain component.
- Timed work records the load; single-thread timings are pinned to the performance cores
  (`taskset -c 0-7`); a change's gain is the mean of three alternating runs, never one pair.
- The data are read-only; bulky outputs are deleted once their numbers are in Progress;
  `scratch/keep/` holds the base export and build, the sets, the harness, the summary and the logs.
- Write prose for people: Progress entries and the report read without other documents, with
  SIMPLE terms explained; no item IDs in the report.

### File table

| file | why | phase |
| --- | --- | --- |
| `doc/implementation_notes/planned/pose_student_training_speed.md` | this plan and its progress (moved to `completed/` at the end of phase 5) | each |
| `doc/implementation_notes/completed/pose_student_training_speed.md` | this plan, completed | 5 |
| `doc/implementation_notes/completed/pose_student_training_speed_report.md` | the review report | 5 |
| `doc/implementation_notes/README.md` | the note's row | 5 |
| `doc/algorithms/pose_student.md` | the training section: the schedule, the chunking, the readout's search | 3 5 |
| `doc/policies/test_environment_policy.md` | the kernel benchmark's description, if its columns change it | 0 5 |
| `src/main/learn/simple_hhcnn.f90` | the network: the profile counters, the chunking, the reduction, the readout, the blocks | 0 1 2 4 |
| `src/main/learn/simple_hhcnn_tester.f90` | its tests | 1 2 4 |
| `src/main/learn/simple_hhcnn_bench.f90` | the kernel benchmark | 0 1 2 4 |
| `src/main/learn/simple_pose_student.f90` | the training loop: batch, early stopping, rate schedule; the profile | 0 3 4 |
| `src/main/learn/simple_pose_student_tester.f90` | its tests | 3 |
| `src/main/learn/simple_polar_feats.f90` | the feature tensor, if it reaches 5% of training | 4 |
| `src/main/strategies/parallelization/simple_strategy_pose_student.f90` | the design choice, the learning bench record | 0 1 3 4 |
| `src/main/commanders/test/simple_commanders_test_highlevel_pose.f90` | the kernel benchmark program | 0 2 |
| `src/main/ui/simple_test/simple_test_ui_highlevel.f90` | the kernel benchmark's help | 0 2 |

## Progress

(filled by the run, one entry per phase)
