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
