# The 2D student: amortized class assignment for the stream and for solve2D at scale

Status, 2026-11-12: approved by Hans ("yes to all"), carried out interactively phase by phase, a report and his
go-ahead after each. Base: `f2e907c24`. Supersedes: none.
Depends on: the 3D pose student with the angle cap and the shared-memory shift registration
(`doc/algorithms/pose_student.md`; validated: yes, five fixed sets, three seeds, `validate_pose_student`
on bgal); the 2D probability table `eul_prob_tab2D` (in production in solve2D's later stages and
`refine2D refine=prob|prob_neigh`; validated: yes). Ideas outside this plan's supervised scope are in
`planned/polar_student_future_directions.md`.

## 1. Why

Clustering particles in 2D is the stream's main recurring cost, and the cost of solve2D on the data
sets of millions of particles the maintainer has to process. Three places pay for it:
- the sieve (p05) runs a whole `solve2D` per coarse and per fine chunk;
- the pool (p06) runs a distributed `refine2D` iteration over a sample of up to 500,000 particles
  (`STREAM_NPTCLS_MAX`) against 150 classes, and gives each new particle a random populated class
  before its first, greedy alignment (`draw_new_classes`);
- solve2D at scale searches every particle against every class at every rotation in its dense stages.

The 3D student shows what an amortized search costs instead. On bgal (1,500 training particles,
24 threads) it trains in 22 s, predicts at 670 particles per second (1.5 ms each), and its poses alone,
with registered shifts, reconstruct a map within the seed spread of the base's (`validate_pose_student`,
2026-11-11). Its readout already scores every class at every in-plane rotation: a 2D alignment against
fixed references is the same computation with class averages for templates. Phase 0 measures where the
2D time goes; no saving is claimed before Phase 4.

## 2. Contract

- **Supervised and per data set.** A 2D student learns to assign particles to one frozen set of class
  averages, the teacher's, from the 2D probability table's rows against those averages, as the 3D student
  learns from the 3D table's rows against a frozen map. It discovers no class. Finding the classes stays
  with solve2D on a subset, or with the pool's own iterations.
- **Off state.** Without a 2D model, and without its candidate artifact, the sieve, the pool, `refine2D`
  and `solve2D` behave byte for byte as on the base at a fixed seed.
- **The 3D student is unchanged by the refactor.** Phase 1 moves code. The 3D programs' models, weights,
  rows, prediction tables and project records stay byte-identical at a fixed seed, at one thread and at
  24.
- **One copy of everything shared.** The network, the features, the training loop, the held-out split,
  the confidence calibration, the design choice, the model file's common part, the provenance, the part
  ranges and the benchmark records exist once. The 2D and 3D students differ only in what their classes
  are, how a candidate set is formed, and how a pose is read out.
- **The library's object rules hold.** Every new stateful module is an encapsulated class: private
  components, `new` and `kill`, behaviour as type-bound procedures, a type named after its module. Large
  components are allocatable. Commanders are thin `execute(cline)` wrappers of strategy routines
  `run_<x>(cline)` and `<x>_exec(params, build, ...)` (`doc/policies/layering_and_project_ownership_policy.md`).

## 3. Decisions and rulings

1. (2026-11-12, Hans) Supervised only. Embeddings, self-supervision, shift-invariant front ends and
   generative models are recorded in `planned/polar_student_future_directions.md` and are not part of
   this plan.
2. (2026-11-12, Hans) The stream is the target application; solve2D at scale on data sets of millions
   of particles is the second.
3. (2026-11-12, Hans) Mirror the 3D path as closely as possible: a learn program, a predict program, a
   validation program, a candidate set handed to a sparse search, fixed-teacher yardsticks over seeds.
4. **The shared core is composed, not inherited.** A new class `polar_student`
   (`src/main/learn/simple_polar_student.f90`) owns the features (`polar_feats`), the network (`hhcnn`),
   the frozen shell weights, the geometry (box, crop, sampling, mask, band, CTF treatment), the training
   on soft class targets with shift and rotation augmentation and early stopping, prediction, the
   readout's template initialisation, the confidence threshold and its record, provenance, the training
   record, and the model file's common header and weights. `pose_student` (3D) and the new
   `class_student` (2D) each hold one `polar_student` as an allocatable component and delegate to it.
   Reason: a Fortran extension in another module cannot reach its parent's private components. A shared
   parent would have to make its state public or grow accessors for every field, and the record shows
   parallel routes beside an owner being rewritten wholesale within weeks (`how2write_agent_plans.md`,
   section 2).
5. **Candidates are class lists in the core.** The core's calibration works on a candidate list of
   classes per particle and an eligibility mask. In 3D the list is the folded cone's directions, and
   eligibility adds the cone-radius bound (ruling R1 of the pose student). In 2D the list is the classes
   holding 99% of the posterior, at most `KMAX_CLASS` of them, with no geometric bound. The held-out
   "recovered" test, the threshold and the design choice are then one routine for both.
6. **The network's reflection branch is a construction property.** `hhcnn%new` gains a logical argument
   for the reflected scores: 3D passes true (mirror pairs), 2D false, because the 2D table searches no
   mirrors. The arithmetic stays a property of the build and the shape, never of the team.
7. **Classes and templates.** A 2D student's classes are the teacher's class averages, all of them,
   rejected classes included. The readout's templates start from the class averages' polar transforms, as
   3D starts them from reprojections (ruling R3 of the pose student). The model records each class's
   state, so a prediction into a rejected class is a rejection signal.
8. **Teacher rows.** `eul_prob_tab2D` gains the posterior-row read the 3D table has
   (`get_posterior_row`): the class posterior marginalised over in-plane rotation and shift at the
   particle's physical temperature, written with the existing row file (`prob_row_file`). The 2D learn
   program fills the dense table against the frozen class averages in one pass, as `learn_pose_student`
   fills the 3D table against `vol1`.
9. **Shifts.** Prediction reads the class and the in-plane angle, then registers the shift against the
   predicted class average at the predicted rotation, with the 3D route's search object (`pftc_shsrch_grad`
   `new_fixed`, coarse grid and L-BFGS-B, shared memory). Training augments shifts over the picking-error
   width Phase 0 measures on real sets, never over the polar shift map's 3-sigma disk (the first 3D
   validation run learned nothing with it).
10. **The candidate artifact.** `predict_classes` writes, per particle, the top class, its rotation, its
    registered shift, the confidence, the route (seeded or dense) and the candidate list. The particle
    record carries the first five in its existing slots. The candidate lists go into a side file that is
    registered in the project's `out` segment with the particle-layout digest (`ptcl_layout_digest`), as
    other per-particle side files are. The sparse 2D table reads it when the artifact is present and its
    digest matches; otherwise it runs as on the base. No command-line key reaches the consumer.
11. **Names.** The class is `class_student`. The programs are `learn_class_student` and `predict_classes`
    (private, in `simple_private_exec`, like the 3D pair) and `validate_class_student` (public, in
    `simple_exec`, like `validate_pose_student`).
12. **No new parameter** beyond what decision 10's artifact and the programs' required keys need. Every new
    number is a named constant with its reason. No second code path stays in the tree; alternatives are
    measured as variants outside it.

13. (2026-11-12, Hans) The development uses the data on the Dell: bgal, trpv1 and the simulated medium
    assembly. trpv1 (35,645 particles) stands in for solve2D at scale in Phase 4. The stream phases are
    measured on the registered stream tests, since neither set has micrographs to replay. The real scale
    (the maintainer's data sets of millions of particles and a real stream) is a set of runs for the
    maintainer, reported, not gated (section 6).
14. (2026-11-12, Hans) Decisions 7 and 10 as written; the pool learns its student at iteration
    `FIRST_EXPORT_ITER` (10); Phases 0 to 4 first, Phases 5 to 7 after their merge.

Defaults taken without a question:
- the angle cap `NET_ANGLE_OVERSAMPLING = 2` applies unchanged;
- `KMAX_CLASS = 16`;
- the 3D student's `CONF_FRACTION` (0.95), held-out rule (every tenth micrograph) and design widths
  (16, 32, 64) apply unchanged;
- the feature cache's format is shared.

## 4. Not in this plan

- Unsupervised clustering, embeddings, self-supervision, shift-invariant or generative front ends
  (decision 1; the future-directions note).
- The 3D student's behaviour, the refine3D and solve3D searches, and the continuous machinery.
- The sieve's chunk lifecycle beyond Phase 6's conditional use, the stream's restart and IPC contracts,
  and the 3D ingestion policy.
- Multi-state or heterogeneity use of class assignments (`refine3D_states_cavgs`).
- The Mac. Every timing is on the Dell; the build without a team-safe BLAS runs the unit suites only.

## 5. Change types and gates

| Phase | Kind | Gate |
| --- | --- | --- |
| 0 | measurement | builds, fast gate, the 3D base profile, the 2D teachers, the picking-error distribution, the time breakdowns |
| 1 | identical (the refactor) | 3D byte identity (models, weights, rows, tables, records) at one and 24 threads; all student suites; CTest `pose_student` |
| 2 | identical for every existing path; tests-only for the new reads | 2D table byte identity on the 2D workflow tests; the new row read's unit test; HHCNN suites on both BLAS builds |
| 3 | new programs | unit suites, a registered workflow test, the reads on the fixed 2D teachers, the feasibility verdict (section 6) |
| 4 | result-changing where the artifact exists; identical otherwise | off-state identity; solve2D at scale, route against route |
| 5 | result-changing in the stream once a model exists; identical before | off-state identity; stream workflow tests; a real-set stream comparison |
| 6 | conditional on Phase 0's breakdown | as Phase 5 |
| 7 | end to end and close | docs, report |

## 6. Acceptance

**Fixed 2D teachers** (frozen in Phase 0, never moved). For each set, a solve2D on a subset of
`nsample` particles gives the frozen class averages, FRCs and assignments. Held-out micrographs judge the
student against that teacher:
- the medium-assembly simulated set of the pose student (truth poses, so truth view clusters) at three
  subset sizes;
- trpv1 (`/data2/trpv1_c2`, its `2_solve2D` as one teacher);
- bgal (`/data2/bgal_addon`);
- trpv1 in full (35,645 particles, 256 px at 1.2156 Å) for the route comparison at scale (decision 13).

**Reads**, per set and seed (seeds 1, 2, 3):
- top-1 agreement with the teacher's class;
- the teacher's class in the candidate list, overall and among seeded particles;
- the median candidate count, the seeded fraction, the in-plane angle error against the teacher's;
- the registered shift's error against the teacher's;
- class averages rebuilt from the student's assignments alone (no search), compared with the teacher's
  by FRC and by `cavg_quality`'s acceptance.

**The feasibility verdict** (guessed, frozen here before Phase 3; it decides only whether Phase 4
starts). On every fixed set, the three-seed mean of "teacher's class among the candidates of a seeded
held-out particle" is at least 0.95 (the 3D calibration fraction), the mean seeded fraction is at least
0.30, and the median candidate count is at most `KMAX_CLASS`. If any set fails, stop with the table.

**Comparisons across seeds.** Each set has three seeds per arm. A read is outside when the arm's
three-seed mean differs from the reference arm's three-seed mean by more than the larger of the two arms'
seed ranges and a floor of 0.01 for fractions (1° for angles). Phase 0 applies the rule to the base
against itself by splitting six base seeds into two groups of three, both ways. The rule is satisfiable
only if no read is outside. This replaces the leave-one-out self-check that stopped `student_crop`.

**Speed** is reported, not gated:
- the learn and predict wall times and their phase splits;
- solve2D at scale, route against route (Phase 4);
- the stream's time to its first and later publications and the pool's time per iteration (Phase 5).

**The maintainer's runs at real scale** (decision 13; reported, not gated): Phase 4's two routes on a data
set of millions of particles, and Phase 5's stream on a real stream, with the commands and scripts this plan
leaves in Progress.

**Reported, not gated (the maintainer's eye):**
- the class-average montages of each route;
- the stream's published class averages, with their paths in Progress.

## 7. Common rules and stop conditions

- Variants are built and measured outside the tree, with their patches kept; only the kept change enters
  the diff. A file outside the table gets its table row before the edit.
- OpenMP: a new or changed parallel region gets a tester that forces a team and compares bit for bit with
  one thread, repeated at 16 threads by hand; no string temporaries in a region; work arrays sized by the
  team in use.
- Timings on an idle machine with the load recorded; a timing taken under load is repeated. The
  maintainer runs other work on this machine.
- Thread identity of training (the HHCNN suites' team checks) on both BLAS builds.

Stop when:
- Phase 0's data preconditions fail (a set's checksums differ, a teacher run fails);
- the base fails Phase 0's self-check;
- the refactor cannot keep the 3D outputs byte-identical;
- decision 4's composition cannot express a needed piece without making the core's state public;
- the feasibility verdict fails;
- the student cannot learn under the measured shift augmentation (record the evidence; the options are
  a centring step, a narrower band, or the future note's shift-robust front end, and the choice is
  Hans's);
- the consumer needs a key beyond decision 10;
- a unit suite of the base fails on the tree and the cause is outside the table.

Gates no ruling may waive:
- the fast gate;
- the student suites (`polar_features`, `hhcnn`, `pose_student`, `pose_student_strategy`, the new 2D
  suites) and the 3D probability table suite;
- CTest `pose_student` and the new 2D workflow test;
- the 3D byte identity of Phase 1;
- the off-state identities of Phases 4 to 6.

## 8. Phases

- **Phase 0, baseline and measurements.**
  - Record the base hash, the machine and the load. Build the base export and the tree (Release and the
    variant without `SIMPLE_BLAS_TEAM_SAFE`) and run the fast gate.
  - The 3D base profile for Phase 1's identity: models, weights, rows and prediction tables of `T200`
    and bgal at seeds 1 to 3, CTest `pose_student`, one `validate_pose_student`.
  - The 2D teachers of section 6: one solve2D per set and subset, frozen under `scratch/keep/teachers/`
    with checksums. Six seeds of the reads for the self-check need a student, so the self-check runs in
    Phase 3 on its first student; Phase 0 records the teachers' own seed spread instead (two solve2D runs
    of one subset).
  - **The picking-error distribution:** the shifts that converged solve2D runs end with, per set, as a
    radial histogram in Å. Phase 3's augmentation width is its 90th percentile.
  - **The 3D student's shift sensitivity at those magnitudes:** `predict_poses`' shift check
    (`SHIFT_CHECK_MAGS`) rerun at the measured percentiles.
  - **Time breakdowns:**
    - a stream run of the registered stream tests (`mini_stream_*`, `stream_chains`): the sieve's
      chunk solve2D time against the pool's iteration time;
    - solve2D on trpv1 in full: stage times and the dense stages' share.
  - Exit: the teachers frozen with checksums; the picking-error table; the time breakdowns; all in
    Progress. If the breakdown shows the sieve, not the pool, dominates the stream, Phase 6 runs before
    Phase 5 (the reorder clause).

- **Phase 1, the shared core (identical).**
  - Create `polar_student` (decision 4) and move into it what `pose_student` holds that is not 3D
    geometry.
  - Move the shared strategy helpers into one stateless module,
    `src/main/strategies/parallelization/simple_strategy_student_common.f90`:
    - the held-out split, the bootstrap quantile, the confidence threshold and the part ranges;
    - the geometry check, the design-choice loop and the learn and predict benchmark records;
    - the feature cache's contract.
  - `learn_pose_student_exec` and `predict_poses_exec` call them.
  - The candidate lists of decision 5 replace the 3D-only calibration inputs.
  - Tests move with their code; the core gets its own tester submodule.
  - Exit:
    - byte identity of the 3D outputs against Phase 0's profile, at one and 24 threads, timing fields
      excluded;
    - every student suite passes;
    - CTest `pose_student` passes;
    - no routine exists in two places (`rg` for the moved names);
    - `scripts/check_test_registry.py` passes.

- **Phase 2, the 2D table's rows and the network without reflection.**
  - `eul_prob_tab2D%get_posterior_row` (decision 8), with a unit test on a synthetic batch: marginals sum
    to one, the temperature matches the 3D read's convention, and ties are drawn as the table draws them.
  - The reflection argument of `hhcnn%new` (decision 6), with HHCNN tests at false: no reflected scores,
    the gradients by finite differences, thread identity at teams of 4 and 16.
  - Exit:
    - the 2D workflow tests (`simulated_workflow suite=1jxy`'s solve2D, `mini_stream`) are byte-identical
      to the base;
    - the new unit tests pass;
    - the HHCNN suites pass on both BLAS builds.

- **Phase 3, the 2D student and its programs.**
  - `class_student` (`src/main/learn/simple_class_student.f90`): one `polar_student`, the classes, their
    states, the candidate lists, the pose read-out (class and in-plane angle).
  - Strategy routines `run_learn_class_student` / `learn_class_student_exec` and `run_predict_classes` /
    `predict_classes_exec` in `simple_strategy_class_student.f90`; commanders as thin wrappers; the UI,
    exec and private-program chains.
  - `validate_class_student` mirrors `validate_pose_student`: a class-balanced subset, a solve2D on it,
    learn, predict on the rest, class averages from the student's assignments, timed phases and a report.
  - The augmentation width from Phase 0.
  - A registered workflow test in miniature (a blob or simulated 2D fixture), with mechanics gates and the
    reads reported.
  - Exit:
    - the suites and the workflow test pass;
    - the reads of section 6 on every fixed set, three seeds (six on one set for the self-check);
    - the feasibility verdict recorded. On failure, STOP.

- **Phase 4, the 2D consumer and solve2D at scale.**
  - The candidate artifact (decision 10): written by `predict_classes`, registered with its layout digest,
    read by the sparse 2D table (`prob_neigh`, the `snhc` sparse fill) for seeded particles. Seeded
    particles score only their candidates; the others run as on the base.
  - solve2D at scale on trpv1 in full (decision 13):
    - route A: the full solve2D;
    - route B: solve2D on `nsample`, learn, predict the rest, then the final stage's sparse refine2D
      seeded by the artifact.
  - Compare the time, the class averages' FRC and `cavg_quality`, and the montages (reported).
  - Exit:
    - the off-state identity (no artifact: byte-identical);
    - route B's numbers beside route A's, two runs of each;
    - an artifact with a stale digest is ignored, and the test for that is mutation-tested.

- **Phase 5, the stream pool.**
  - The pool learns a `class_student` as an asynchronous job, like its iterations, once iteration
    `FIRST_EXPORT_ITER` is reached. It learns on a sample of the pool's particles against its current
    class averages.
  - New particles (`append_sets`) take the student's class, rotation and shift instead of
    `draw_new_classes`' random class, and seeded particles search only their candidates in the pool's
    iterations.
  - **Keeping the student current as the classes change:** re-template (new class averages into the
    readout, a short fine-tune of the existing network) or retrain. Both are measured as variants, and
    the kept one is chosen on the reads and the time.
  - Exit:
    - the off-state identity (before the first model, and in a stream without one, byte-identical);
    - the stream workflow tests pass;
    - on the registered stream tests, the time to each publication, the pool's time per iteration and the
      published class averages' `cavg_quality`, beside the base's, two runs each;
    - the commands for the maintainer's real-stream run in Progress (decision 13).

- **Phase 6, the sieve (conditional).** Runs only if Phase 0's breakdown shows that the sieve's chunk
  solve2D holds at least a third of the stream's 2D time.
  - Once the pool has a student, fine chunks are assigned by it and judged by `cavg_quality`, in place of
    a chunk solve2D.
  - Exit as Phase 5.

- **Phase 7, close.**
  - `doc/algorithms/class_student.md` (the algorithm page, with the reads and the measurements).
  - Updates to `doc/policies/stream/pool2D_policy.md`, `doc/policies/2D/solve2D_policy.md` and, if Phase 6
    ran, `ptcl_sieve_policy.md`.
  - The report `completed/class_student_2D_report.md`, and this plan moved to `completed/`.
  - The fast gate, `unit_numerics`, the stream and 2D workflow tests, and the CTest entries on both
    builds.

Planning note: Phases 0 to 4 and Phases 5 to 7 suit two runs. The stream phases start from the merge of the
first run, so that master's 2D and stream changes are not outrun (`how2write_agent_plans.md`, section 7).

### File table

| File | What | Phases |
| --- | --- | --- |
| `src/main/learn/simple_polar_student.f90` (new) | the shared core class | 1, 3 |
| `src/main/learn/simple_polar_student_tester.f90` (new) | its tester submodule | 1, 3 |
| `src/main/learn/simple_pose_student.f90` | 3D geometry around the core | 1 |
| `src/main/learn/simple_pose_student_tester.f90` | tests that stay 3D | 1 |
| `src/main/learn/simple_class_student.f90` (new) | the 2D class | 3, 5 |
| `src/main/learn/simple_class_student_tester.f90` (new) | its tests | 3, 5 |
| `src/main/learn/simple_hhcnn.f90`, `simple_hhcnn_tester.f90` | the reflection argument | 2 |
| `src/main/learn/simple_polar_feats.f90`, `simple_polar_feats_tester.f90` | the cache contract, if it moves | 1 (if needed) |
| `src/main/strategies/parallelization/simple_strategy_student_common.f90` (new) | the shared stateless helpers | 1, 3 |
| `src/main/strategies/parallelization/simple_strategy_student_common_tester.f90` (new) | their tests | 1, 3 |
| `src/main/strategies/parallelization/simple_strategy_pose_student.f90`, `_tester.f90` | 3D routines on the shared helpers | 1 |
| `src/main/strategies/parallelization/simple_strategy_class_student.f90` (new), `_tester.f90` (new) | the 2D routines | 3, 4 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab2D.f90`, `simple_eul_prob_tab2D_tester.f90` | the row read; the candidate-restricted sparse fill | 2, 4 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_utils.f90` | row types, if shared fields are needed | 2 (if needed) |
| `src/main/commanders/simple/simple_commanders_class_student.f90` (new) | the thin commanders, validation's controller | 3 |
| `src/main/commanders/simple/simple_commanders_pose_student.f90` | shared validation steps, if any move | 1 (if needed) |
| `src/main/ui/simple_private_prgs.f90`, `src/main/ui/simple/simple_ui_*.f90` (the 2D program's group) | program registration | 3 |
| `src/main/exec/simple_exec_*.f90`, `src/main/apis/simple_private_exec_api.f90`, `simple_exec_api.f90` | routing | 3 |
| `src/main/commanders/test/simple_commanders_test_highlevel_class_student.f90` (new), `simple_test_registry.f90`, `simple_commanders_test_class.f90`, `src/main/ui/simple_test/*.f90`, `production/CMakeLists.txt` | the workflow test and suites | 1, 3 |
| `src/main/project/` (the `out` segment's artifact registration only) | the candidate artifact | 4 (if needed) |
| `src/main/stream/pool2D/simple_stream_pool2D.f90`, `simple_stream_refine2D_utils.f90`, `simple_stream_pool2D_tester.f90` | the pool's learn job and seeding | 5 |
| `src/main/stream/stages/simple_stream_stage_pool2D.f90` | the stage's hand-off of the model | 5 (if needed) |
| `src/main/sieve/` | fine-chunk assignment | 6 (conditional) |
| `doc/algorithms/class_student.md` (new), `doc/algorithms/pose_student.md` | algorithm pages | 1, 7 |
| `doc/policies/stream/pool2D_policy.md`, `doc/policies/2D/solve2D_policy.md`, `doc/policies/sieving_and_rejection/ptcl_sieve_policy.md` | the living rules | 4-7 |
| `doc/policies/test_environment_policy.md` | the new tests' rows | 1, 3 |
| `doc/implementation_notes/planned/class_student_2D.md` | this plan, Progress | 0-7 |
| `doc/implementation_notes/completed/class_student_2D.md`, `..._report.md`, `doc/implementation_notes/README.md` | the close | 7 |

## Progress

### Phase 0, baseline and measurements (2026-11-12, the Dell, base `f2e907c24`)

**Builds and reference outputs.** `/data2/agent_runs/class_student/build_base`: a Release build of the base
(GCC 15.2, `-O3 -march=native`, BLAS in OpenMP teams on). The 3D reference outputs for Phase 1's byte
identity are in `keep/p3d/`:
- `learn_pose_student` and `predict_poses` on `T200` and bgal at seeds 1 to 3, every output kept (model
  header, weights, rows, feature cache, prediction table); training took 7.6 to 8.1 s on `T200` and 20.9 to
  22.1 s on bgal;
- the `pose_student` workflow test (36 s, metrics in `keep/p3d/workflow/`);
- `validate_pose_student` on bgal (404 s).

The script is `tools/chain_p0.sh`, its log `logs/chain_p0.log`.

**The 2D teachers** (`teachers/`, each a copy of the import project, a random subset by
`selection nran=... prune=no`, whose other particles stay as state 0 for testing, then `solve2D` on the
subset at 24 threads):

| teacher | particles (subset of) | classes | mask | solve2D | class FRC 0.5 (50% of particles) |
| --- | --- | --- | --- | --- | --- |
| `bgal_s1` (seed 1) | 2,756 (5,513) | 50 | 180 Å | 35 s | 14.8 Å |
| `bgal_s2` (seed 2, another subset) | 2,756 (5,513) | 50 | 180 Å | 35 s | 15.5 Å |
| `trpv1_10k` | 10,000 (35,645) | 100 | 175 Å | 139 s | 13.5 Å |
| `trpv1_full` | 35,645 | 100 | 175 Å | 468 s | 8.2 Å |
| `sim_T5000` | 5,000 | 50 | 300 Å | 81 s | 24.0 Å |

The simulated medium assembly is a poor 2D teacher. Its class averages reach only 24 Å, and its shifts drift
from the truth's (median 24 Å, at most 42) to a median of 41 Å and a maximum of 926 Å, in a 384 Å box. It
stays a 3D yardstick and is reported, not judged, in 2D.

**The picking error**, measured as the radius of the shift each converged solve2D ends with (imported
particles start at zero):

| set | median | 75% | 90% | 95% | 99% |
| --- | --- | --- | --- | --- | --- |
| trpv1, all 35,645 | 8.0 Å | 11.6 Å | 15.9 Å | 20.5 Å | 40.0 Å |
| trpv1, 10,000 | 8.1 Å | 11.8 Å | 17.0 Å | 24.3 Å | 45.6 Å |
| bgal, 2,756 | 10.1 Å | 15.7 Å | 22.1 Å | 27.5 Å | 41.2 Å |

**The 3D student's shift sensitivity at that error**, from `predict_poses`' shift check on bgal (1.275 Å per
pixel, trained with `shift_sigma` 2 Å): the folded mode is kept for 73% of particles at 2 px (2.6 Å), 21% at
5 px (6.4 Å) and 0.4% at 10 px (12.8 Å, near the median picking error). Its cone still holds the mode for
96% at 5 px but only 28% at 10 px.

**Augmentation does not close the gap** (`runs/shiftprobe_s*`, bgal, seed 1, the student trained with wider
shift disks and read on the centred particles of the full set):

| training `shift_sigma` (disk 3 sigma) | within 6° | truth in cone | held-out seeded | cone holds the mode at 10 px |
| --- | --- | --- | --- | --- |
| 2 Å (6 Å, as now) | 0.501 | 0.557 | 0.313 | 0.28 |
| 4 Å (12 Å) | 0.387 | 0.542 | 0 | 0.63 |
| 7 Å (21 Å, about the 90th percentile) | 0.176 | 0.485 | 0 | 0.94 |

Robustness to shifts costs accuracy, and at the picking error's scale the student no longer seeds anyone. On
these data and this network, a student is only useful on particles centred to within about 3 Å. The stream's
pool meets that: its particles arrive with the sieve's chunk 2D shifts (`append_sets` keeps them). solve2D at
scale and the sieve do not.

**solve2D's time** on trpv1 in full (35,645 particles, 100 classes, 24 threads, 468 s):

| stage | mode | time |
| --- | --- | --- |
| 1 | `snhc_smpl` | 33 s |
| 2 | dense `prob` | 98 s |
| 3, 4, 5 | sparse `prob_neigh` | 61, 60, 63 s |
| 6 | dense `prob` | 139 s |

- The two dense stages hold 52% of the time.
- Per iteration, the dense table fill takes 8 to 15 s of 16 to 23 s, and the sparse fill 3 to 4 s of 11 to
  13 s.
- The rest, about 8 s per iteration, is reading 9.3 GB of particles, preparing them and assembling the class
  sums. That part does not shrink with a student; at millions of particles it is the floor.

**The stream's breakdown could not be measured here.** `mini_stream_1jxy` simulates 6 micrographs and
20 particles (48 s), with no sieve-to-pool hand-off of any size. Phase 6's condition (the sieve's share)
therefore waits for the maintainer's real stream; Phase 5 (the pool) stays first.

**What it means, in short.**
- The teachers are fixed and fast to remake.
- The real sets' picking error is 8 to 10 Å at the median and 16 to 22 Å at the 90th percentile.
- The student needs particles centred to about 3 Å, which the pool supplies and raw picks do not. Decision 9
  (train over the picking-error width) does not survive this measurement, and its replacement is for Hans to
  choose.
- In solve2D at scale the student can replace at most the dense table fills, about half the time, because
  reading and preparing particles stays.

