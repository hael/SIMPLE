# Pose student: the review's corrections and training speed, run report

The autonomous run `student_speed` (2026-11-09 to 2026-11-10, on the Dell: Oracle Linux 8.10, Core
Ultra 9 285, 24 cores, GCC 15.2) worked on the pose student. The pose student is a small network
trained per dataset on the posterior rows of the probabilistic 3D search. For particles the search
has not seen, it proposes a direction, a cone about it and a confidence
(`doc/algorithms/pose_student.md`). The network is HHCNN, the Harmonic Hermitian convolutional neural
network. The run first corrected the seven findings of the code review
(`learned_pose_state_assignment_review.md`, PS01 to PS07). It then made the training faster without
changing what the student learns beyond the run-to-run spread of its reads. The plan and its full
record are in `pose_student_training_speed.md` beside this report; its Progress section holds every
phase's detail. Nothing is committed. The change is one working-tree diff against `6541633df`, and
the driver wrote each phase's part to `/home/elmlundho/agent_runs/student_speed/review/phase_N.diff`
(N = 0 to 5, and 6 when this phase ends).

**Terms.**
- *The fixed teachers* are training inputs kept constant for the whole run, so that a change in the
  student is not hidden by a change in its teacher. *`T200`, `T1000`, `T5000` and `T`* are 200,
  1,000, 5,000 and 24,000 particles of a simulated ribosome assembly with known poses; the student
  predicts an unseen test half and is read against the truth. *bgal* is the 1,500 training particles
  of beta-galactosidase from an earlier validation run, read against that run's solve3D poses (the
  teacher).
- *The reads*: within 6° (the student's most probable direction within 6° of the truth or the
  teacher), in the cone (the truth or teacher inside the student's cone), and seeded (the fraction of
  particles whose confidence clears the calibrated threshold and is handed to the neighbourhood
  search). *Held-out recovery* is the share of held-out particles whose teacher direction falls in
  the cone.
- *The spread* is the range of a read over three training seeds. *The corrected spread* is that
  range after the review's corrections (phase 1), the reference for the speed phases.
- *Load* is the machine's one-minute load average during a timing; the maintainer used the machine
  during the run.
- *The blocks* are the network's three convolution blocks. Each is *im2col* (the 3 by 3
  neighbourhoods laid out as columns of a matrix) followed by a matrix product, a *gate* on the
  modulus, and pooling.
- *BLAS* is the matrix-product library (OpenBLAS here). A *team-safe* BLAS stays on the calling
  thread inside an OpenMP team (the build definition `SIMPLE_BLAS_TEAM_SAFE`). *The variant build*
  is a build of the same source without a team-safe BLAS, as on the review's macOS machine. *FSC* is
  the Fourier shell correlation, the resolution of a map.

## What was corrected

Each finding has a regression test. Unless stated otherwise the tests are in the pose student's
tester or the strategy's new tester (suite `pose_student_strategy` of `unit_numerics`). The review
file's closing section, "Resolution in student_speed", holds each change, its test and its evidence.

- **The training arithmetic changed with the thread team (PS01).** The network chose between the
  BLAS product and explicit loops at run time, by whether an OpenMP team was active. So on a build
  without a team-safe BLAS, one thread and a team took different arithmetic. The choice is now the
  build's alone. Test: the HHCNN tester's identity check asserts that its team actually formed (a
  serial run now fails instead of passing) and compares teams of 4 and 16 threads with one thread,
  bit for bit, on both builds.
- **A model could load another network (PS02).** The weights file is resolved against the model
  header's directory, never the working directory. The reader refuses weights whose shape, seed or
  weight count differ from the header's. Tests: two model directories with the same default names, a
  foreign weights file in the working directory, a relocated model. The workflow test predicts with a
  relocated model.
- **Prediction workers ignored their part (PS03).** A part of `predict_poses` now requires its range:
  `fromp` and `top`, with `part`, `nparts` and `outfile`. These are existing keys, handed to every
  part as the distributed masters already hand them. A range without a selected particle writes its
  records unchanged, and part tables carry the part's number. Tests: uneven ranges cover a selection
  exactly once; a part of only unselected rows; refused ranges. The workflow test runs four uneven
  parts and merges their documents.
- **Prediction accepted another physical geometry (PS04).** Prediction refuses a model whose sampling
  (within 1e-4 Å), mask (within 0.01 Å) or crop differ from the project's. Test: each is refused with
  the same box and half ring; differences inside the header's decimals pass.
- **Calibration and prediction could select different particles (PS05).** The confidence threshold
  ends only where the confidence drops, never inside a group of equal confidences. It is written so
  that it reads back bit for bit, and the trained model is read back and checked to seed exactly the
  calibrated held-out particles. Tests: the review's tie of eighteen recovered and two missed rows,
  thresholds next to the old five-decimal boundary, the round trip.
- **A training set could have no held-out cohort (PS06).** The held-out split takes every tenth of the
  selected particles' stacks, in order, when there are ten or more stacks; otherwise every tenth
  selected particle. Fewer than ten selected particles are refused before training. Tests: two to
  nine stacks, fifteen stacks none numbered a multiple of ten, small single-stack selections.
- **Wide-band models could not be read back (PS07).** The header reader takes records of any length.
  Test: models of 227 and 228 shells round-trip bit for bit.

The model file's version is now 4; release 4 keeps no compatibility with earlier models. The
corrections changed no trained network on the fixed teachers: the held-out tables were the base's on
all fifteen runs. The seeded fractions moved by at most 0.0013, because the threshold now keeps its
full precision.

## What was gained, and at what cost to the reads

**Training wall time** of `learn_pose_student` at 24 threads on the fixed teachers (mean of three
seeds in seconds; the base is the code before the run):

| set (training particles) | base | after the run | factor |
| --- | --- | --- | --- |
| `T200` (200) | 88.9 | 8.5 | 10.5 |
| `T1000` (1,000) | 192.5 | 18.3 | 10.5 |
| bgal (1,500) | 351.3 | 34.7 | 10.1 |
| `T5000` (5,000) | 468.4 | 86.1 | 5.4 |
| `T` (24,000) | 606.1 | 369.3 | 1.6 |

Where the gain came from:
- **The schedule** gave most of it on the small and middle sets. Training now counts in epochs, with
  a patience of three epochs; the base waited 250 optimiser steps after the best epoch, which is 14
  to 21 epochs on the small sets. The 64-channel design is tried only from 10,000 training rows.
- **The step itself** got faster:
  - the readout as matrix products;
  - the batch cut into chunks by its size alone (a target of 96; 15 to 86 chunks in practice), so
    that more threads work on small batches, with the reduction in parallel;
  - the blocks' im2col built in tiles that stay in a core's cache.
- **On `T`** the schedule changed little, because the base already stopped after 7 to 9 epochs. Of
  the 369 s, 56 s are the teacher's table fill, which the run kept out of scope.

**The kernel benchmark** (`test=pose_student_bench`, a training step of 256 samples at 24 threads)
reaches 1.4 to 2.4 times the base's samples per second: 2.0 to 2.4 times at 16 channels, 1.4 to 1.5
at 64.

**The registered workflow test** `pose_student` (2,000 simulated particles, eight threads, deadline
2,400 s) passes through CTest on every build. It takes 31.9 s on the default build and 51.5 s on the
variant build, against 307.2 s for the base build the same evening (313.9 s at the run's start). On
the review's macOS build the base had timed out at 2,400 s.

**End to end.** `validate_pose_student` ran on bgal twice per build: 1,500 training particles of
2,756. Its teacher is a fresh, unseeded solve3D in every run, so only the learning time compares
directly. The map column is the unseen particles reconstructed at the student's poses with the shifts
registered there, FSC 0.143 against the teacher's map, in Å.

| run | solve3D, the teacher (s) | training (s) | widths trained; kept | training rows seeded | within 6° of the teacher | unseen seeded | map (Å) | total (s) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| base 1 | 321.2 | 474.8 | 16, 32, 64; 32 | 0.412 | 0.859 | 0.044 | 4.80 | 855 |
| base 2 | 320.6 | 185.9 | 16, 32; 16 | 0.625 | 0.915 | 0.043 | 5.63 | 564 |
| run 1 | 655.3 | 50.8 | 16, 32; 32 | 0.545 | 0.843 | 0.143 | 6.40 | 761 |
| run 2 | 288.3 | 36.4 | 16, 32; 32 | 0.484 | 0.890 | 0.048 | 4.87 | 381 |

- **Learning time:** 44 s on average against 330 s, 7.6 times faster. The plan's baseline runs had
  taken 178 and 537 s. Every held-out recovery is 1.000.
- **Run 2** lies between the two base runs on the table's reads, apart from its unseen seeded
  fraction (0.048 against 0.043 and 0.044).
- **Run 1** had the weakest teacher of the four. Its solve3D took twice as long as the others', and
  its control map (the training particles at the teacher's poses) reached 4.47 Å between its halves,
  against 4.08 Å in the other three runs. Its student is below both base runs within 6° (0.843) and
  in the unseen map (6.40 Å), above them in unseen seeded (0.143), and between them in training rows
  seeded.

Two runs per build cannot separate the student from the teacher; the fixed-teacher reads are the
comparison.

**The cost to the reads**, on the fixed teachers:
- **The means.** Every read's three-seed mean lies inside the corrected spread, except two on bgal.
  Training rows seeded is 0.414, against the range 0.415 to 0.437 (0.396 after phase 2, 0.407 after
  phases 3 and 4). Seeded rows in the cone is 0.998, above the range 0.995 to 0.997.
- **Per seed**, two runs move by more than rounding:
  - bgal seed 1 seeds 0.386 of its training rows against 0.419 at the base. Those confidences cluster
    just above 0.99, and the threshold is a quantile of 150 held-out rows, so any change that regroups
    a floating-point sum moves it this much.
  - `T5000` seed 2 keeps 32 channels where the base kept 64, under the rule that tries 64 channels
    only from 10,000 rows: within 6° 0.883 against 0.913, seeded 0.481 against 0.509. This seed sets
    `T5000`'s mean apart from the base's in the second decimal (0.870 against 0.881 within 6°).
  - Every other read of every set and seed differs from the base's in the third or fourth decimal.

**The plan's working targets** were five times faster training on every fixed set and three times the
base's kernel samples per second. The first is met on four sets of five (`T`: 1.6 times); the second
is not (1.4 to 2.4 times).

## Decisions taken under delegation

Each is recorded with its numbers in the plan's Progress entry for the phase named.

- **The variant build** uses the system pthreads OpenBLAS (`-DSIMPLE_BLAS_OPENMP=OFF`, with
  `BLAS_LIBRARIES` and `LAPACK_LIBRARIES` set to `/usr/lib64/libopenblasp.so`), because this
  machine's plain `libopenblas.so` is single-threaded (phase 0).
- **The kernel benchmark's one-thread passes** had run with the team's thread count set, so the
  OpenMP BLAS spread them over the team. They now run at one thread (phase 0).
- **The truth sets' test half** is the copy that carries the truth shifts. The copy the plan named had
  the shifts removed; on it the student found the truth within 6° for only 3 to 5% of particles
  (phase 0; the plan's data section corrected).
- **The review's open choices:**
  - PS03: the part's range is required rather than derived.
  - PS04: the tolerances are the header's decimals.
  - PS06: every tenth stack in order of number (phase 1).
- **The chunk target is 96 chunks per batch**, against 64 by three alternating benchmark runs and
  against 48 by one. Reduction slices are 4,096 weights (phase 2).
- **The templates' gradient** moved out of the per-chunk buffers into one product over the batch per
  block of 64 classes, an extension the plan did not list (phase 3).
- **The coarse-to-fine rotation search** was built, measured and dropped. It was 1.7% faster at 24
  threads and 8% slower at one thread, and it changed the winning rotation in 0.14 to 0.22% of the
  class searches. The backward skip of negligible classes was not built, because it could gain under
  5% (phase 3).
- **The schedule** is at most 30 epochs with a patience of 3, and the 64-channel design only from
  10,000 training rows. A warm-up with rate halving, a doubled batch and successive halving were each
  measured on every set and seed, and dropped for their reads or their time (phase 4).
- **The blocks** (phase 5):
  - the explicit loops reordered to run down columns: the variant build trained 14 times slower than
    the default, now 1.6 to 1.75 times;
  - the gate without branches;
  - im2col in tiles of 128 columns.

## Deviations from the plan

- **The variant build's thread identity in the kernel benchmark** at 24 threads was checked on three
  of its twelve shapes; the loops would have taken hours for the rest. The HHCNN tester checks the
  identity on that build at 4 and 16 threads (phase 2).
- **The readout's sums are now single precision**, as matrix products; the explicit loops had summed
  in double (phase 3). Scores and gradients agree to single-precision rounding.
- **The tile sizes 64 and 256** were measured once each, not three times; both ran within 2% of 128
  (phase 5).
- **The default build's BLAS** is the system OpenMP OpenBLAS, not a bundled one as the plan assumed.
- **`T200` after the corrections:** its training took 13% longer than the base's on average, with all
  three seeds slower. The corrections did not change its networks, and the phase did not find the
  cause; the later phases made the question moot.

## Left open for the maintainer

- **Whether the reads are equal**, which the plan leaves to the maintainer. In particular: bgal's
  seeded fraction at the edge of the corrected spread, and the threshold's sensitivity behind it; and
  `T5000` seed 2's narrower model.
- **The review's macOS configuration** (the system OpenBLAS without the team-safe definition) was
  checked only through the Linux variant build. The workflow test's timeout at eight threads there is
  the case to rerun.
- **`T`**, the 24,000-particle set, trains 1.6 times faster, short of the five-times target. What
  remains there:
  - the designs, about 300 s; at 16 and 32 channels the blocks' matrix products, already near the
    core's peak, are half to two-thirds of a sample;
  - the teacher's rows, 56 s, out of scope.
- **The readout** is the largest share left outside the blocks: 14 to 18% of a sample at 16 channels.

## Evidence

All under `/data2/agent_runs/student_speed/scratch/keep/` (the run directory's scratch):

- **Tests at the end, both builds:**
  - the fast gate (`logs/p6_gate_{tree,noteam}.log`);
  - `unit_numerics`, 28 suites and 2,006 assertions (`logs/p6_unit_numerics_{tree,noteam}.log`);
  - the workflow test through CTest (`logs/p6_ctest_pose_student_{tree,noteam,base}.log`);
  - the HHCNN suite at teams of 4 and 16, three times (`logs/p5_hhcnn_*.log`).
- **The fixed teachers:** `summary.tsv` (one line per build, phase, set and seed), `runs/` (logs and
  bench records), `tools/harness.sh` and `tools/spread.sh`.
- **The kernels:** `bench/kernels_*.txt` (the base: `kernels_base_1.txt`; the final tree:
  `kernels_tree_p5final_1.txt`).
- **The validation:** `validation_p6/val_{tree,base}_{1,2}/` (logs, reports, bench records) and
  `tools/val_stats.sh`.
- **The variants measured and dropped:** `variants/`.

## Measurements

**Training wall time** (s) per set and seed (seeds 1 / 2 / 3), at 24 threads, under a load of 2 to
26. "Phase 4" is the schedule the run kept.

| set | base | after the corrections | phase 2 | phase 3 | phase 4 | phase 5 (final) |
| --- | --- | --- | --- | --- | --- | --- |
| `T200` | 79.2 / 100.0 / 87.5 | 98.7 / 105.1 / 99.4 | 25.5 / 27.7 / 28.6 | 18.9 / 20.8 / 20.3 | 8.0 / 9.1 / 8.5 | 8.2 / 8.8 / 8.5 |
| `T1000` | 118.3 / 333.1 / 126.0 | 122.9 / 339.6 / 127.5 | 69.1 / 178.8 / 70.9 | 54.5 / 158.2 / 57.0 | 22.3 / 23.2 / 24.4 | 17.9 / 18.3 / 18.8 |
| bgal | 168.9 / 456.7 / 428.3 | 160.4 / 455.1 / 424.3 | 117.6 / 312.8 / 292.7 | 100.0 / 284.8 / 265.6 | 44.7 / 44.1 / 44.9 | 35.8 / 33.9 / 34.4 |
| `T5000` | 310.3 / 784.0 / 310.9 | 314.1 / 788.6 / 312.0 | 309.0 / 785.6 / 303.9 | 275.2 / 738.5 / 268.9 | 126.5 / 132.4 / 120.8 | 85.4 / 90.3 / 82.5 |
| `T` | 586.0 / 630.3 / 601.9 | 588.0 / 630.9 / 605.5 | 573.5 / 620.7 / 598.2 | 520.1 / 558.5 / 534.9 | 518.0 / 558.3 / 537.1 | 359.0 / 379.8 / 369.1 |

**The reads per seed** (seeds 1 / 2 / 3; within 6°; in the cone; seeded) at the base, after the
corrections and at the end:

| set | base | after the corrections | final |
| --- | --- | --- | --- |
| `T200` | 0.690 / 0.747 / 0.792; 0.962 / 0.970 / 0.973; 0.117 / 0.190 / 0.305 | the same | 0.690 / 0.748 / 0.792; 0.962 / 0.970 / 0.973; 0.116 / 0.191 / 0.300 |
| `T1000` | 0.756 / 0.813 / 0.823; 0.974 / 0.984 / 0.987; 0.219 / 0.245 / 0.260 | the same | 0.756 / 0.812 / 0.823; 0.974 / 0.984 / 0.986; 0.227 / 0.248 / 0.261 |
| bgal | 0.874 / 0.871 / 0.877; 0.958 / 0.956 / 0.966; 0.419 / 0.414 / 0.437 | seeded 0.420 / 0.415 / 0.437, the rest the same | 0.869 / 0.871 / 0.877; 0.960 / 0.957 / 0.966; 0.386 / 0.421 / 0.435 |
| `T5000` | 0.837 / 0.913 / 0.892; 0.979 / 0.992 / 0.989; 0.370 / 0.509 / 0.452 | the same | 0.836 / 0.883 / 0.892; 0.978 / 0.985 / 0.989; 0.368 / 0.481 / 0.452 |
| `T` | 0.930 / 0.914 / 0.946; 0.996 / 0.991 / 0.997; 0.464 / 0.547 / 0.586 | seeded 0.464 / 0.547 / 0.585, the rest the same | 0.931 / 0.912 / 0.946; 0.996 / 0.991 / 0.996; 0.467 / 0.536 / 0.584 |

**The reads' means** by phase (within 6°; in the cone; seeded):

| set | base | phase 2 | phase 3 | phase 4 | final |
| --- | --- | --- | --- | --- | --- |
| `T200` | 0.743; 0.968; 0.204 | 0.743; 0.968; 0.215 | 0.743; 0.968; 0.206 | 0.743; 0.968; 0.206 | 0.743; 0.968; 0.202 |
| `T1000` | 0.797; 0.981; 0.241 | 0.797; 0.982; 0.242 | 0.797; 0.981; 0.240 | 0.797; 0.981; 0.240 | 0.797; 0.981; 0.245 |
| bgal | 0.874; 0.960; 0.424 | 0.873; 0.961; 0.396 | 0.873; 0.961; 0.407 | 0.873; 0.961; 0.407 | 0.872; 0.961; 0.414 |
| `T5000` | 0.881; 0.986; 0.443 | 0.881; 0.986; 0.443 | 0.880; 0.986; 0.441 | 0.870; 0.984; 0.434 | 0.870; 0.984; 0.434 |
| `T` | 0.930; 0.995; 0.532 | 0.930; 0.995; 0.532 | 0.930; 0.995; 0.531 | 0.930; 0.995; 0.531 | 0.930; 0.995; 0.529 |

The other phases' reads per seed are in `summary.tsv`.

**The kernel benchmark's training samples per second** at 24 threads, by rings / channels:

| rings / channels | base | final | rings / channels | base | final |
| --- | --- | --- | --- | --- | --- |
| 12 / 16 | 1,899 | 4,485 | 25 / 16 | 1,016 | 2,165 |
| 12 / 32 | 984 | 1,895 | 25 / 32 | 486 | 832 |
| 12 / 64 | 396 | 613 | 25 / 64 | 188 | 270 |
| 18 / 16 | 1,343 | 3,176 | 37 / 16 | 720 | 1,478 |
| 18 / 32 | 687 | 1,363 | 37 / 32 | 323 | 574 |
| 18 / 64 | 269 | 388 | 37 / 64 | 126 | 185 |

## Proposed skill text

The maintainer applies skill edits. For `simple-modern-fortran`, section OpenMP:

> - A tester's thread-identity check asserts that its team formed (`omp_get_num_threads()` inside a
>   region against the count asked for): under `OMP_THREAD_LIMIT=1` or a serial runtime the
>   comparison passes trivially. Check teams of 4 (the fast gate's cap) and, by hand, 16.
> - Whether the BLAS runs inside a team is the build's property (`SIMPLE_BLAS_TEAM_SAFE`), never a
>   run-time test of `omp_in_parallel()`: a one-thread region is not an active team, so a run-time
>   test gives one thread and a team different arithmetic.

For the guide of `src/main/learn` (or `simple-main-strategies-parallelization`), section Tests:

> - The pose student: `unit_numerics suite=hhcnn`, `suite=pose_student`, `suite=pose_student_strategy`;
>   the workflow `test=pose_student` (CTest, eight threads); the kernels `test=pose_student_bench`
>   (no CTest entry; run with `OMP_PLACES=cores OMP_PROC_BIND=close`).

## Commit message

```text
Pose student: the review's corrections and faster training

Correct the code review's seven findings, each with its test: the training
arithmetic is the build's alone (identical at any team), a model loads only
its own weights, prediction parts take their range, refuse another geometry
and select exactly the calibrated particles, every selection has a held-out
cohort, and wide-band models read back (model version 4). Then train faster:
epochs with a patience of three, the widest design only from 10,000 rows,
the batch in chunks by its size with a parallel reduction, the readout as
matrix products, the blocks' im2col in cache-sized tiles. Training on the
fixed teachers is 10 times faster below 2,000 rows, 5.4 times at 5,000 and
1.6 times at 24,000; the reads' means stay within the seed spread apart
from bgal's seeded fraction, at its edge.
```
