# Pose student code review

**Reviewed:** 2026-11-09, `master` at `519161ff395de08f72b4a97d09629c9a0b86264e`.
**Development:** `71c4a41aa`, compared with its parent `20d24017d` and current callers.
**Status:** completed review; seven open findings, no source fixes applied.

The development needs corrections before acceptance. The Release build succeeds,
but the new HHCNN thread-identity assertion fails the fast gate. Six further
findings concern model selection, partition ownership, physical geometry,
confidence calibration, holdout selection and model-file capacity.

The review covers the private CLI wiring, teacher-row export, polar features,
HHCNN forward and backward calculations, student training and model I/O,
per-particle cone handoff to `prob_neigh`, orientation-record storage and the
registered tests. The reviewed pose-student sources are unchanged between the
development commit and the reviewed HEAD.

Related records: [implementation and decisions](learned_pose_state_assignment.md),
[development run report](learned_pose_state_assignment_report.md) and
[algorithm description](../../algorithms/pose_student.md).

## Findings

P1 denotes an issue requiring prompt correction; P2 denotes a correctness or
robustness issue to address in the normal development cycle. All findings below
remain open at the reviewed revision.

| ID | Priority | Finding | Evidence |
| --- | --- | --- | --- |
| PS01 | P1 | Training arithmetic changes with the active thread team | Fast-gate and isolated-suite failures |
| PS02 | P2 | Model loading can silently select another network | CLI reproduction with valid foreign weights |
| PS03 | P2 | Prediction workers do not derive their partition range | CLI reproduction, including an empty selected range |
| PS04 | P2 | Prediction accepts incompatible physical geometry | CLI reproduction with an internally consistent model header |
| PS05 | P2 | Calibrated prefixes and deployed thresholds can select different rows | Source analysis and a counterexample emulation |
| PS06 | P2 | A valid input selection can have no holdout cohort | Source analysis of split, training and bootstrap paths |
| PS07 | P2 | Wide-band models exceed the model reader's line buffer | CLI reproduction at 228 shells |

### PS01 P1 Training arithmetic changes with the active thread team

**Location:** [`blas_ok` and the matrix products](../../../src/main/learn/simple_hhcnn.f90#L430),
lines 430–480.

When `SIMPLE_BLAS_TEAM_SAFE` is absent, `blas_ok()` permits SGEMM outside an active
OpenMP team and selects explicit loops inside one. A one-thread parallel region
therefore takes a different arithmetic path from a multi-thread region. Different
rounding occurs before the fixed-order chunk reduction, so that reduction cannot
provide the promised bit-for-bit weight identity.

On the reviewed macOS Release build, GNU Fortran 16 with system OpenBLAS and no
`SIMPLE_BLAS_TEAM_SAFE` definition, both the fast gate and the isolated HHCNN
suite fail this assertion:

```text
the chunked reduction gives the same weights at one thread and at a team of 4, bit for bit
```

The isolated suite passes 27/28 assertions. Capping the runtime with
`OMP_THREAD_LIMIT=1` makes it pass 28/28 because both requested team sizes execute
serially; that diagnostic does not validate an actual four-thread team.

**Correction and regression coverage:** choose a consistent arithmetic path across
team sizes on supported BLAS configurations. Retain the existing test with an
actual four-thread team and check configurations both with and without the
team-safe BLAS definition. The fast gate must pass before installation.

### PS02 P2 Model loading can silently select another network

**Location:** [`pose_student%read`](../../../src/main/learn/simple_pose_student.f90#L768),
lines 768–773.

The reader first interprets `weights_file` against the process working directory.
It falls back to the weights beside the requested header only when that first
path does not exist. Consequently, an absolute `infile` does not uniquely select
its associated network when the header stores a relative weights name.

An isolated probe used a valid model header and weights in one directory, then
placed a different valid network named `student_weights.bin` in the prediction
working directory. With the same requested header and particle, prediction
changed from `pmax=0.5`, `k99=2` to `pmax=1.0`, `k99=1`, and completed normally.
The header's calibration and provenance can therefore be paired with another
network without an error.

**Correction and regression coverage:** make the writer and reader agree on a
header-relative weights-path contract, and validate the binary shape against the
header. Test two model directories with the same default filenames and a foreign
weights file in the working directory. Relocating a complete model must preserve
which network is loaded.

### PS03 P2 Prediction workers do not derive their partition range

**Location:** [particle filtering](../../../src/main/strategies/parallelization/simple_strategy_pose_student.f90#L733),
lines 733–737, and [part-file writing](../../../src/main/strategies/parallelization/simple_strategy_pose_student.f90#L838),
lines 838–840. Parameter defaults are set in
[`simple_parameters_phases`](../../../src/main/params/simple_parameters_phases.f90#L659).

The program advertises `part` and `nparts`, but neither its wrapper nor its body
derives `fromp` and `top` from those keys. With their defaults, filtering and
part-file writing both cover the full project.

On an isolated 2,000-row project containing 1,000 selected particles,
`part=2 nparts=2000` predicted all 1,000 selected particles and wrote all 2,000
project records. Multiple workers using only the advertised partition keys thus
produce overlapping contributions. Supplying `fromp=2 top=2` explicitly for an
unselected row instead aborts at line 737 with `no selected particle`.

**Correction and regression coverage:** derive and validate the worker range
through the established partition machinery, or require explicit ranges in the
program contract. A valid range with no selected particles must still produce
the contribution required by the merger. Check disjoint coverage, uneven ranges
and parts containing only unselected rows.

### PS04 P2 Prediction accepts incompatible physical geometry

**Location:** [model checks and band assignment](../../../src/main/strategies/parallelization/simple_strategy_pose_student.f90#L721),
lines 721–729.

Prediction checks `pftsz`, native `box` and CTF treatment, but does not compare the
model's recorded `smpd`, `box_crop` or mask geometry. Equal array dimensions do
not establish equal physical Fourier sampling or equal image preparation.

An internally consistent synthetic model specified 4 Å/pixel, a 200 Å mask and
a 16 Å band. It was accepted for a project at 2 Å/pixel with a 100 Å mask and
matching pixel geometry. Prediction completed normally after adopting the
model's shell indices. At the project's sampling, those indices represent a
different physical band from the one recorded in the model.

**Correction and regression coverage:** validate the full preprocessing and
physical-geometry contract before applying frozen weights and templates, with
appropriate serialization tolerances. Any supported resampling must be explicit.
Cover differing sampling, crop and mask settings even when `box` and `pftsz`
remain equal.

### PS05 P2 Calibration and prediction can select different populations

**Location:** [threshold calibration](../../../src/main/strategies/parallelization/simple_strategy_pose_student.f90#L481),
lines 481–490; [prediction's threshold comparison](../../../src/main/strategies/parallelization/simple_strategy_pose_student.f90#L813);
and [threshold serialization](../../../src/main/learn/simple_pose_student.f90#L620).

Calibration can retain a qualifying prefix ending inside a confidence tie.
Prediction uses `conf >= conf_threshold`, which includes the entire tie rather
than the selected prefix.

A counterexample emulation uses 20 eligible rows with confidence 1, ordered as
18 recovered rows followed by two misses. The prefix of 18 qualifies, but the
stored threshold 1 routes all 20, whose recovery is 90%, below the declared 95%.
This establishes a possible failure of the selection rule; it is not a measured
frequency of ties in the workflow fixture. Independently, saving the threshold
with `F8.5` can broaden or narrow the selected population after model reload.

**Correction and regression coverage:** evaluate thresholds only after complete
confidence groups and preserve enough precision to retain the calibrated
membership. Recompute the accepted population after a model-file round trip.
Include tied confidences and values near a serialization boundary.

### PS06 P2 A selected training set can have no holdout cohort

**Location:** [holdout selection](../../../src/main/strategies/parallelization/simple_strategy_pose_student.f90#L132),
lines 132–144; [zero-group bootstrap division](../../../src/main/strategies/parallelization/simple_strategy_pose_student.f90#L439);
and [training without held-out rows](../../../src/main/learn/simple_pose_student.f90#L363).

For a project with multiple stacks, the split uses raw `stkind modulo 10`.
Projects whose stack IDs are 1–9 have no holdout rows. The same occurs when the
selected subset excludes every stack divisible by ten. The single-stack fallback
also has no holdout when fewer than ten particles are selected.

The source then permits training to substitute training loss for held-out loss.
The design comparison reaches a bootstrap mean with `ngroups=0`, dividing zero
by zero, while threshold calibration retains `1.00001` and seeds nobody. This
finding follows the source paths and split examples; a complete no-holdout
training run was not performed during the review.

**Correction and regression coverage:** form a deterministic split over eligible
groups that leaves both cohorts nonempty, or reject an insufficient cohort before
training. Explicitly handle zero bootstrap groups. Cover two through nine stacks,
selected subsets with no tenth stack and small single-stack selections.

### PS07 P2 Wide-band models exceed the model reader's line buffer

**Location:** [fixed-size input record](../../../src/main/learn/simple_pose_student.f90#L667),
[shell-weight reads](../../../src/main/learn/simple_pose_student.f90#L747) and
[shell-weight writes](../../../src/main/learn/simple_pose_student.f90#L628).

The writer emits each complete shell-weight vector on one line using `ES17.9`
values separated by commas. At 228 shells, those lines contain 4,111 and 4,110
characters. The reader's `character(len=4096)` buffer truncates them before the
internal list-directed read.

A valid synthetic 228-shell header in the writer's format aborted at line 748
with `Fortran runtime error: End of file`. Thus the model format can represent
bands that its reader cannot reload, including the reload of a candidate model
during training.

**Correction and regression coverage:** read complete variable-length records or
store the shell weights in a representation without this line-length limit.
Check exact model round trips on both sides of the 227/228-shell boundary.

## Observed validation

The following results belong to the reviewed macOS Release build, GNU Fortran 16
and system OpenBLAS. The review applied no source fixes. Installation was skipped
because the fast gate failed. No Linux/BOX validation was performed.

| Check | Observed result |
| --- | --- |
| `make -C build -j8` | Passed, including complete UI JSON validation; no compiler warnings in the build log. |
| `scripts/run_fast_gate.sh build` | 14/15 suites passed in 43.0 s. `unit_numerics` passed 1,925/1,926 assertions; its sole failed assertion was HHCNN thread weight identity. |
| New learn suites within the fast gate | Polar features 39/39; pose student 38/38; HHCNN 27/28. All 359 finite-difference gradient entries matched. |
| Other affected areas in the fast gate | Orientation, project and polar alignment suites passed. |
| Isolated `test=unit_numerics suite=hhcnn` | Reproduced PS01 with `OMP_THREAD_LIMIT=4`. The one-actual-thread diagnostic passed 28/28. |
| Direct `test=pose_student nthr=1` | Passed in 1,449 s: eight contract checks and all 33 metric rows. Prediction, search, fallback merge and both reconstructions completed. |
| Registered CTest `pose_student`, eight threads | Timed out at 2,400.03 s while training width 32, epoch 91; prediction, search and reconstruction were not reached. |

The serial workflow seeded 5/1,000 test particles, with truth inside 3/5 cones.
These quality values are ungated reports. Its successful contract checks do not
establish 95% truth coverage, useful production seeding rates or thread identity.
The eight-thread timeout is an observed runtime failure; this review did not
isolate its complete performance cause.

## Reproduction entry points and retained evidence

The build and gate entry points were:

```sh
make -C build -j8
scripts/run_fast_gate.sh build
env OMP_WAIT_POLICY=PASSIVE OPENBLAS_NUM_THREADS=1 \
    ctest --test-dir build -R '^pose_student$' --output-on-failure
```

The isolated HHCNN run used `OMP_NUM_THREADS=1`, `OMP_THREAD_LIMIT=4`,
`OMP_WAIT_POLICY=PASSIVE` and `OPENBLAS_NUM_THREADS=1`. Changing only the thread
limit to 1 produced the serial diagnostic. Direct test invocations used the
fresh executable in `build/production`, the build as `SIMPLE_PATH`, local queue
settings and separate temporary working directories. The one-thread workflow
used `test=pose_student nthr=1` with passive waiting and one OpenBLAS thread.

The CLI probes used isolated project copies and valid synthetic network files to
exercise file resolution, geometry and partition handling. Their observed
outcomes are included above so that the findings do not depend on retaining the
temporary fixtures. Local evidence from this review is available at:

- `/private/tmp/simple-pose-review-build.log`
- `/private/tmp/simple-pose-review-fast-gate.log`
- `/private/tmp/simple-pose-review-hhcnn.log`
- `/private/tmp/simple-pose-review-hhcnn-serial-diagnostic.log`
- `/private/tmp/simple-pose-review-highlevel.log`
- `/private/tmp/simple-pose-review-highlevel-serial.log`
- `/private/tmp/simple-pose-review-serial/metrics.tsv`
- `/private/tmp/simple-pose-review-model-path-valid.log`
- `/private/tmp/simple-pose-review-sampling.log`
- `/private/tmp/simple-pose-review-empty-part.log`
- `/private/tmp/simple-pose-review-empty-range.log`
- `/private/tmp/simple-pose-review-bigband.log`
- `/private/tmp/simple-pose-review-probes/`

These local files are transient and are not tracked repository artifacts.
After corrections, rerun the fast gate, the focused regression cases above and
the registered workflow. A serial pass cannot substitute for actual-team
identity or for the registered workflow's completion within its budget.

## Resolution in student_speed

The autonomous run `student_speed` corrected the seven findings in its phase 1 (2026-11-10, on the Dell:
Oracle Linux 8.10, GCC 15.2). Every finding has its regression test; the evidence below is under
`/data2/agent_runs/student_speed/scratch/keep/logs/`. Two builds were checked: the default, which links
the system OpenMP OpenBLAS and defines `SIMPLE_BLAS_TEAM_SAFE`, and a variant linking the system
pthreads OpenBLAS without it (`BLAS in OMP teams: OFF`), the configuration of this review's macOS build.
The model file's version is now 4.

**PS01, the training arithmetic changed with the team.** Reproduced on Linux before the correction: on
the variant build the isolated HHCNN suite with `OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=4
OMP_WAIT_POLICY=PASSIVE OPENBLAS_NUM_THREADS=1` failed the identity check, 30 of 31
(`ps01_hhcnn_team4_noteam.log`). Change (`simple_hhcnn.f90`): the matrix products take `sgemm` when the
build's BLAS is team-safe and the explicit loops otherwise, whatever the team; the run-time test of an
active team (`blas_ok`) is gone, so the arithmetic is a property of the build. Test
(`simple_hhcnn_tester.f90`): the identity check also asserts that its team of four formed, so a serial
run (`OMP_THREAD_LIMIT=1`) now fails instead of passing. Evidence: with the actual team of four, 32 of 32
on both builds (`p1_ps01_hhcnn_team4_tree.log`, `p1_ps01_hhcnn_team4_noteam.log`); the serial diagnostic
fails the team check (`p1_ps01_hhcnn_limit1_noteam.log`); the fast gate passes on both builds
(`p1_fast_gate_tree.log`, `p1_fast_gate_noteam.log`). The variant build's HHCNN suite takes 26 s against
7 s: its training runs on the loops at any thread count.

**PS02, a model could load another network.** Change (`simple_pose_student.f90`): the header names the
weights by their file name alone, resolved against the header's directory and never the working
directory; a name with a directory is refused; after reading, the network's half ring, rings, classes,
width, seed and weight count must equal the header's, else the reader stops. Tests (the pose student's
tester, now a submodule of the module): two model directories written with the same default names and a
foreign `student_weights.bin` in the working directory each load their own network bit for bit; the
relocated model loads its own; a weights name with a directory, a missing name, and networks of another
seed, width or weight count are refused. Workflow test: the model copied to another directory beside a
file of the weights' default name that is no network predicts the same table row for row.

**PS03, prediction workers ignored their part.** The distributed masters hand every part `fromp`, `top`,
`part` and `nparts` together (`qsys_ctrl%augment_partition_job_descr`); no master runs `predict_poses`
yet. Change: the program's contract requires the range. A part needs `nparts`, `fromp`, `top` and
`outfile` (the commander), `1 <= part <= nparts` and `1 <= fromp <= top <=` the particles (the strategy);
`fromp` and `top` are optional keys of the private program (existing keys). A part whose range holds no
selected particle writes the header-only table and its range's records unchanged instead of stopping.
A part's tables carry its number (`<fname>_part<n>_pred.txt`), so that parts given one `fname`, as a
master gives it, do not overwrite each other's; a part refuses `vol1` (the shift registration runs in
shared memory) before it looks at its range, empty or not.
Tests: the strategy's tester (`simple_strategy_pose_student_tester.f90`, suite `pose_student_strategy`
of `unit_numerics`) covers uneven ranges covering a selection once, a part of only unselected rows, and
refused ranges and parts; the workflow test runs `predict_poses` as four uneven parts with one shared
`fname` (1-700, 701-1998, one unselected row, one selected row), checks that each predicts only its range, that together they
cover every selected row once with the whole run's table rows, and that `merge_algndocs` of the four
documents (which requires every record covered) gives the whole run's records.

**PS04, prediction accepted other physical geometry.** Change: prediction compares the model's half ring,
box, crop (`box_crop`) and CTF treatment exactly and its sampling within 1e-4 Å and mask within 0.01 Å
(the header's five and three decimals) with the project's, and stops on a difference; nothing is
resampled. Test (strategy tester): another sampling, mask or crop with the same box and half ring is
refused, differences within the decimals pass.

**PS05, calibration and prediction could select different particles.** Change: the threshold is taken
only where the confidence drops (the end of a group of equal confidences), so that `conf >= threshold`
at prediction selects exactly the calibrated prefix; it is written as `ES17.9`, which reads back bit for
bit, and so is the route's cone bound it was calibrated with (formerly three decimals);
learn_pose_student reads the written model back and stops if its threshold and bound would seed other
held-out rows than the calibration. Tests: the review's eighteen recovered and two missed rows at confidence one give no
threshold; the same tie below a recovered row ends the prefix before the tie; a recovered tie is taken
whole; ineligible rows are passed over (strategy tester); thresholds next to the five-decimal boundary
read back bit for bit and select the same population (pose student tester).

**PS06, a training set could have no held-out cohort.** Change: the split ranks the selected particles'
stacks by number and holds out every tenth when there are ten or more; with fewer, every tenth selected
particle; fewer than ten selected particles are refused before the rows are made. The resampled design
comparison of zero groups gives zero gain, and `pose_student%train` refuses a set without held-out rows
(the substitution of the training loss is gone). Tests: two to nine stacks, fifteen stacks none of whose
numbers is a multiple of ten, nine particles on one or on nine stacks (refused), ten on one stack, and the
zero-group quantile (strategy tester). The bgal training set's split changes with this rule (its held-out
stacks are now every tenth selected stack, not every stack whose number is a multiple of ten).

**PS07, wide-band models could not be read back.** Change: the header reader reads records of any length
in pieces. Test (pose student tester): models of 227 and 228 shells round-trip with their shell weights
bit for bit and the same prediction.

**After the corrections.** The fast gate passes on both builds (15 of 15; `unit_numerics` 10.8 s on the
default build, 30.6 s on the variant, whose network runs on the loops), with the new suite
`pose_student_strategy` (41 checks) and the pose student's (64) and HHCNN's (32) suites. The registered
workflow test `pose_student` passes through CTest at its eight threads in 317.7 s against the base's
313.9 s, with its thirteen checks, four of them new: the relocated model, the parts' ranges, their
coverage and agreement with the whole run, and the merge of their documents
(`p1_ctest_pose_student_tree_final.log`, `p1_pose_student_metrics_tree.tsv`). On the fixed sets of the
run (three seeds on four truth-set subsets and on bgal), the trained networks and held-out splits are
the base's on every set and seed (the held-out tables are identical); the confidence thresholds now keep
their full precision, and the seeded fractions move by at most 0.0013 either way (the five-decimal form
had rounded the calibrated threshold both ways: 0.9900014 down to 0.99000, admitting a few more
particles, and 0.9901377 up to 0.99014, admitting fewer).

**After the run's speed phases** (2026-11-10). The identity check of PS01 now compares teams of 4 and 16
threads with one thread bit for bit, and a new check runs the blocks across several im2col tiles; the
HHCNN suite passes 41 of 41 on both builds, by default and three times under `OMP_THREAD_LIMIT=24`
(`p5_hhcnn_{tree,noteam}*.log`). The variant build's explicit loops were reordered: its training is now
1.6 to 1.75 times the team-safe build's, against 14 times before. At the end the fast gate and
`unit_numerics` (28 suites, 2,006 assertions) pass on both builds (`p6_gate_{tree,noteam}.log`,
`p6_unit_numerics_{tree,noteam}.log`), and the workflow test `pose_student` passes through CTest at its
eight threads in 31.9 s on the default build and 51.5 s on the variant, against 307.2 s for the base
build (`p6_ctest_pose_student_{tree,noteam,base}.log`). The macOS configuration itself was not
available to the run.
