# The 2D student: amortized class assignment for the stream and for solve2D at scale

Status, 2026-11-12: approved by Hans ("yes to all"), carried out interactively phase by phase, a report and his
go-ahead after each. Paused 2026-11-13 after Phase 0 for `completed/center_ptcls.md`, which centres the particles the
student needs, and resumed the same day with `center_ptcls` in the working tree: unseen particles are centred by
it before they are labelled (decision 9). Stopped 2026-11-13 at Phase 3's feasibility verdict, and
resumed 2026-11-14 by ruling 15 (the verdict judged on the class averages by Hans's inspection). **Closed
2026-11-14 after Phase 4** (decision 17): Phases 5 to 7 are replaced by `planned/class_student_stream.md`, owned
by Joe and Cyril; the algorithm page is `doc/algorithms/class_student.md`. Base: `f2e907c24`. Supersedes: none.
Depends on: the 3D pose student with the angle cap and the shared-memory shift registration
(`doc/algorithms/pose_student.md`; validated: yes, five fixed sets, three seeds, `validate_pose_student`
on bgal); the 2D probability table `eul_prob_tab2D` (in production in solve2D's later stages and
`refine2D refine=prob|prob_neigh`; validated: yes); `center_ptcls` (`completed/center_ptcls.md`; validated: yes on
trpv1 and bgal and by its workflow test, uncommitted). Ideas outside this plan's supervised scope are in
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
9. **Shifts and centring** (2026-11-13, Hans: mirror the pose student, and centre the unseen data with
   `center_ptcls`). The student works on centred particles, as the pose student always did.
   - **Training particles** are the teacher's, centred by the teacher's solve2D. The student trains on them
     with the 3D student's narrow augmentation (`shift_sigma` 2 Å).
   - **Unseen particles** are centred by `center_ptcls` before they are labelled: in `validate_class_student`,
     in Phase 4's route at scale, and wherever raw picks meet the student. The stream's pool needs no such
     step, because its particles arrive centred by the sieve's chunk 2D.
   - **After prediction**, the shift is registered against the predicted class average at the predicted
     rotation, with the 3D route's search object (`pftc_shsrch_grad` `new_fixed`, coarse grid and L-BFGS-B,
     shared memory).
   - **The two frames differ.** `center_ptcls`'s centres and the teacher's differ by a common translation
     and, per view, by each set of class averages' own centring. The registration absorbs both. Phase 3
     measures them (section 6).
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

15. (2026-11-14, Hans) **The feasibility verdict is judged on the class averages.** The frozen verdict (section 6)
    failed on the exact-class read (Progress, Phase 3). Hans inspected the class averages of the student's
    assignments on bgal (seeds 1 to 3) and judged them good enough: bgal's D2 symmetry gives the teacher
    near-duplicate views, between which the exact-class read counts a miss. The ruling replaces the gate by his
    inspection; the exact-class reads and the class averages' FRC are reported, not gated. A numeric FRC rule (the
    student's median within one shell of the dense pass's) was proposed first and fails on the runs (1.5 shells
    behind on bgal, 1.8 on trpv1, three-seed means): it is not the rule. trpv1's averages are to be inspected
    too. Phase 3 passes; Phase 4 starts.
16. (2026-11-14, Hans) The CTest budget of 46 (the `class_student` entry) is confirmed.
17. (2026-11-14, Hans) The plan closes after Phase 4. The stream phases (5 to 7) go to their own plan,
    `planned/class_student_stream.md`, for Joe and Cyril, written from what this run measured.
18. (2026-11-14, Hans, after the close) The class student trains width 16 alone: `CLASS_STUDENT_NDESIGNS` = 1 caps
    the designs tried; raising it restores the choice. Replaces the default below for the class student's widths.

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

**Fixed 2D teachers** (frozen in Phase 0, never moved; `/data2/agent_runs/class_student/teachers/`). Each is a
solve2D on a random subset, whose other particles are the unseen ones:
- trpv1: 10,000 of 35,645 particles, 100 classes (`trpv1_10k`);
- bgal: 2,756 of 5,513 particles, 50 classes, two subsets (`bgal_s1`, `bgal_s2`);
- trpv1 in full (`trpv1_full`, 100 classes), for the route comparison at scale (decision 13) and as the
  converged shifts of every trpv1 particle;
- the medium-assembly simulated set: reported, not judged (Phase 0 found it a poor 2D teacher).

**Unseen test particles in two states**, so that the centring's cost is visible:
- raw picks centred by `center_ptcls` (the route's case);
- the same particles with the full run's converged shifts (the centring's ceiling; trpv1 only).

**The frame offset:** `center_ptcls`'s shifts against the teacher's on the training particles, as their common
translation and the residual after it.

**Reads**, per set and seed (seeds 1, 2, 3):
- top-1 agreement with the teacher's class;
- the teacher's class in the candidate list, overall and among seeded particles;
- the median candidate count, the seeded fraction, the in-plane angle error against the teacher's;
- the registered shift's error against the teacher's;
- class averages rebuilt from the student's assignments alone (no search), compared with the teacher's
  by FRC and by `cavg_quality`'s acceptance.

**The feasibility verdict** (guessed, frozen here before Phase 3; it decides only whether Phase 4
starts), judged on the unseen particles centred by `center_ptcls`. On every fixed set, the three-seed mean of
"teacher's class among the candidates of a seeded held-out particle" is at least 0.95 (the 3D calibration fraction), the mean seeded fraction is at least
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
- the student passes the feasibility verdict on the unseen particles with converged shifts but fails it on
  the same particles centred by `center_ptcls` (record both tables; the options are a better centring or the
  future note's shift-robust front end, and the choice is Hans's);
- the consumer needs a key beyond decision 10;
- a unit suite of the base fails on the tree and the cause is outside the table.

Gates no ruling may waive:
- the fast gate;
- the student suites (`polar_features`, `hhcnn`, `polar_student`, `pose_student`, `student_strategy`, the new 2D
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
    radial histogram in Å. (Planned as Phase 3's augmentation width; decision 9 as revised replaced that use.)
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

- **Phase 3, the 2D student and its programs.** (Unseen particles centred by `center_ptcls`, decision 9.)
  - `class_student` (`src/main/learn/simple_class_student.f90`): one `polar_student`, the classes, their
    states, the candidate lists, the pose read-out (class and in-plane angle).
  - Strategy routines `run_learn_class_student` / `learn_class_student_exec` and `run_predict_classes` /
    `predict_classes_exec` in `simple_strategy_class_student.f90`; commanders as thin wrappers; the UI,
    exec and private-program chains.
  - `validate_class_student` mirrors `validate_pose_student`:
    - a class-balanced subset and a solve2D on it (the teacher, which centres the training particles);
    - learn;
    - `center_ptcls` on the unseen particles, called as a child command line;
    - predict them;
    - class averages from the student's assignments, timed phases and a report.
  - The augmentation: the 3D student's (decision 9).
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
    - route B: solve2D on `nsample`, learn, `center_ptcls` on the rest, predict them, then the final
      stage's sparse refine2D seeded by the artifact.
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
| `src/main/project/simple_class_candidates.f90` (new), its tester (new) | the candidate artifact's class: the side file, its layout digest, publish and open (row added 2026-11-14, before the edit) | 4 |
| `src/main/stream/pool2D/simple_stream_pool2D.f90`, `simple_stream_refine2D_utils.f90`, `simple_stream_pool2D_tester.f90` | the pool's learn job and seeding | 5 |
| `src/main/stream/stages/simple_stream_stage_pool2D.f90` | the stage's hand-off of the model | 5 (if needed) |
| `src/main/sieve/` | fine-chunk assignment | 6 (conditional) |
| `doc/algorithms/class_student.md` (new), `doc/algorithms/pose_student.md` | algorithm pages | 1, 7 |
| `doc/policies/stream/pool2D_policy.md`, `doc/policies/2D/solve2D_policy.md`, `doc/policies/sieving_and_rejection/ptcl_sieve_policy.md` | the living rules | 4-7 |
| `doc/policies/test_environment_policy.md` | the new tests' rows | 1, 3 |
| `doc/implementation_notes/completed/class_student_2D.md` | this plan, Progress | 0-7 |
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
- The student needs centred particles, as the pose student always had; the pool supplies them and raw
  picks do not. Decision 9
  (train over the picking-error width) does not survive this measurement (replaced 2026-11-13: decision 9), and its replacement is for Hans to
  choose.
- In solve2D at scale the student can replace at most the dense table fills, about half the time, because
  reading and preparing particles stays.

### Probes after Phase 0 (2026-11-13, the Dell): the shift question, closed by ruling

Three probes asked whether the 2D student could work on raw picks (`runs/p1_trpv1`, `runs/p1_floor`,
`runs/p2`, `runs/p3_search`).

- **Classical centring** (a shift search against the rotationally averaged mean of the class averages, all
  35,645 trpv1 picks): within 3 Å of the full run's converged shift for 19% of particles (raw picks 9%;
  median residual 5.7 Å).
- **Power-spectrum features** (|c|^2 in place of (Re, Im), shift-invariant by construction): the student
  learned almost nothing (within 6° 0.006 on T1000 against 0.80 with phases; 0.05 on bgal against 0.50).
  The phases carry the information.
- **The full search from raw picks** (refine2D against the teacher's 100 class averages): 30%, 32% and 36%
  within 3 Å after one, two and three passes. Two converged solve2D runs agree for 69% of the same
  particles (1.9 Å median), and bgal's two runs for 76% (1.5 Å).

So a single particle at the student's band does not fix its centre to 3 Å in a few passes, whoever estimates
it. Ruling (Hans, 2026-11-13): mirror the pose student, which always worked on centred particles; decision
9 as revised. Raw picks are then centred by `center_ptcls` (`completed/center_ptcls.md`, built the same day) before
the student labels them.

### Phase 1, the shared core (2026-11-13, the Dell, HEAD `38e033cbf` with `center_ptcls` in the tree)

**What moved.**
- `src/main/learn/simple_polar_student.f90` (new) holds the core class `polar_student`. It owns the features,
  the network, the band and the shell weights; the training on soft class targets; prediction; the readout's
  template initialisation; the threshold and training records; provenance; and the model file's common
  header and weights, with the header reader `polar_student_header`.
- `simple_pose_student.f90` is now the 3D wrapper around one core: the grid, its mirror pairs, the folded cone,
  the read-out pose, the cone bound and the `vol1` hash. Its public interface is unchanged.
- Its training computes each row's targets folded onto the mirror pairs, as before but once per row before
  training instead of once per batch, then calls the core's training.
- `src/main/strategies/parallelization/simple_strategy_student_common.f90` (new) holds the stateless strategy
  helpers, moved from the pose student's strategy:
  - the held-out cohort, the bootstrap quantile, the confidence threshold, the geometry check, a part's range,
    `file_hash`;
  - their constants (`CONF_FRACTION` and `BOOT_QUANTILE` are public).
- **Tests moved with the code.**
  - New suite `polar_student`, 23 checks: the targets, the weights-path and weights-mismatch checks moved
    from the pose student's tester, and a new round trip of the core's header with its owner's keys.
  - The pose student's suite keeps 50 checks.
  - The helpers' tester moved to `simple_strategy_student_common_tester.f90`; its suite is now
    `student_strategy` (41 checks).

**Deviations.**
- **Header line order.** The model header's lines are now the core's block, the shell weights, then the
  owner's keys. Readers are key-based, so the format and the version (4) are unchanged. The identity check
  therefore compares the header as a set of lines.
- **What waits for Phase 3.** The design-choice loop and the benchmark records stay inside
  `learn_pose_student_exec`, because they are internal procedures on its local state. They move to the common
  module in Phase 3, when the class student's learn routine is their second user.

**Identity** (`/data2/agent_runs/class_student/keep/p3d_ref` against `keep/p3d_tree`). Both are Release builds,
24 threads: the snapshot of the tree before the refactor, and the refactored tree.
- On `T200` and bgal at seeds 1 to 3, these are byte-identical: the weights, the rows (`.bin` and `.txt`) and
  the prediction tables.
- The model headers are identical as sets of lines, excluding `time_train_s` and `teacher_project`, the run
  folder's path.
- The feature caches are identical beyond their first 32 bytes, the contract with the project's path and
  modification time.
- The `pose_student` workflow test's metrics are identical, timings excluded.
- At one thread, the weights and rows of `T200` seed 1 are byte-identical between the two builds.

**A pre-existing finding, not from this phase.**
- In both builds, `T200` seed 1 trains to different weights at one thread than at 24 (held-out loss 1.61094
  against 1.60973). The rows and features are identical.
- So the program-level training is not thread-identical, although the HHCNN suite's team checks pass. The
  likely cause is a BLAS product called outside a parallel region, which then runs threaded. Left for Hans.

**Gates.**
- the fast gate, 15 of 15;
- the suites: `polar_student` 23, `pose_student` 50, `student_strategy` 41, `polar_features` 48, `hhcnn` 41;
- `check_test_registry.py` consistent.

**A process note.** The tester's move was made with `git mv`, which staged the rename in the index. It was
reported to Hans and left as it is.

### Phase 2, the 2D table's rows and the network without reflection (2026-11-13, the Dell)

**The teacher's rows.** `eul_prob_tab2D` gains two procedures:
- `keep_batch_rows` collects the dense batch's rows from the fill's thread buffers before the assignment. When
  rows are kept, `assign_batch` assigns from them (checking that they are the batch's); otherwise it collects
  them as before.
- `get_posterior_row`, the 3D table's read for the 2D table: the class posterior at beta times the particle's
  temperature over the classes at their best in-plane rotations, the kept classes in descending mass, the
  tail, the row statistics, the mode's rotation and refined shift, and optionally the whole row's losses.

Dense table only; the sparse table refuses.

New test `test_dense_posterior_rows`:
- the masses, their ratio exp(-T dL), kmax and the tail;
- the mode's class, loss, rotation and refined shift;
- the argmin read when a particle has no temperature;
- the batch assigned from the kept rows as from the buffers.

The 2D table suite passes 112 of 112.

**The reflection switch.**
- `hhcnn%new` takes a required `reflect`. Without it, the score product covers only the map's columns and the
  rotation search skips the reflected ones, so no class is scored by the reflected map and its backward pass
  never reaches those columns.
- `read_weights` keeps the setting the network was built with.
- `polar_student%new` passes it through: the pose student `.true.`, the class student `.false.`.
- **Tests:** the gradient check and the thread identity now run at both settings. A new
  `test_no_reflection` checks two networks with the same weights: without the reflection, no class is scored
  by the reflected map; each class scores at most what it scores with it, and the same where the reflection
  did not win the class. The comparison is at single-precision rounding, because the score product's width
  changes OpenBLAS's blocking of its last column (differences near 1e-7 relative were seen).
- The HHCNN suite passes 52 of 52 on the default build and on the build without a team-safe BLAS, and 46 of 46
  with a team limit of 4 (its 16-thread checks skip there).

**Existing paths unchanged.**
- `center_ptcls` on bgal (the dense 2D table, references and pass) gives byte-identical particle records and
  all 24 class-average files byte-identical, between the snapshot before Phase 1 and the tree.
- `mini_stream_1jxy` has identical metrics (timings excluded) and identical 2D progress lines.
- The 3D student on `T200` and bgal, seed 1, has byte-identical weights, rows and predictions against Phase 1's
  reference.
- The fast gate passes 15 of 15. The suites: `polar_student` 23, `pose_student` 50, `student_strategy` 41,
  `polar_features` 48.

### Phase 3, the 2D student and its programs (2026-11-13, the Dell): stopped at the verdict, passed by ruling 15

**What was built** (uncommitted):
- `simple_class_student.f90` and its tester (suite `class_student`, 12 checks): the classes, their states, the
  candidates (the classes holding 99% of the posterior, at most `KMAX_CLASS`), the read-out angle, the model
  file's class keys.
- `simple_strategy_class_student.f90`: `learn_class_student_exec` (the dense 2D table's rows against frozen class
  averages, prepared as refine2D's matcher prepares them, then the designs) and `predict_classes_exec` (the class,
  angle, candidates and confidence; given the class averages, the shift registered against the predicted class
  at the predicted angle with the 3D route's search object; the records class, e3, shift, `pconf`, `route`).
- `simple_commanders_class_student.f90`: the two thin commanders (private programs) and the controller
  `validate_class_student` (public, developer visibility, in the refine2D UI group).
- The design-choice loop and the training constants (`STUDENT_WIDTHS`, `WIDEST_MIN_ROWS`, `MAX_EPOCHS`,
  `PATIENCE_EPOCHS`, `SHIFT_AUG_NSIG`, `KMAX_ROW`) and `fname_base` now live once, in
  `simple_strategy_student_common.f90`. The 3D student stays byte-identical (`T200` and bgal, seed 1: weights,
  rows and prediction tables).
- The workflow test `class_student` (`simple_commanders_test_highlevel_class_student.f90`): `validate_class_student`
  on the `center_ptcls` entry's fixture, its contract gated, the reads reported. 36 s at 8 threads (Release), 1:49
  on the Debug build, both passing. The CTest budget is now 46 (confirmed, decision 16).
- Not built: a tester for the class strategy (its helpers are the common module's, tested there); the benchmark
  records stay per program; `predict_classes` has no distributed parts.

**Two corrections on the way.**
- **The training disk.** At refine2D's default `shift_sigma` (10 Å, a 23.5 px augmentation on bgal) the student
  learned almost nothing (bgal_s1: 6% of held-out rows recovered). At the plan's 2 Å (4.7 px) it recovered 96%.
  `learn_class_student` now defaults to 2 Å (decision 9's width).
- **The yardstick.** The dense pass draws each class from the tempered posterior (beta 0.3): 11% of its classes were
  not the row's mode. The validation now runs it at `prob_beta` 100 (`DENSE_PROB_BETA`), so its class is the row's
  mode, as the held-out rows are judged. The dense pass also runs at the full box (the teacher's averages are
  full-box, and the class averager does not crop).

**The held-out reads** (the frozen teachers, seed 1; the times of these runs are the Debug build's and are not
reported):

| teacher | kept width | held-out recovered | top-1 | in-plane error (median, agreeing) |
| --- | --- | --- | --- | --- |
| `bgal_s1` | 16 | 0.98 | 0.66 | 0.9° |
| `bgal_s2` | 16 | 0.94 | 0.64 | 0.8° |
| `trpv1_10k` | 16 | 0.82 (seeded 0.825, 0.995 of them recovered) | 0.44 | 1.1° |

Width 32 never beat 16 by the resampled rule.

**The feasibility verdict** (`validate_class_student`, Release copy of the tree, 24 threads, seeds 1 to 3; bgal:
2,756 training and 2,757 unseen particles, 50 classes; trpv1: 10,000 and 25,645, 100 classes). Each cell is the
three-seed mean of the unseen particles against the dense pass's class (each row's mode); "converged" is the same
student on the same particles at the dense pass's shifts:

| set, shifts | top-1 | seeded | dense class among a seeded particle's candidates | median candidates |
| --- | --- | --- | --- | --- |
| trpv1, `center_ptcls` | 0.370 | 0.704 | **0.904** (0.917, 0.894, 0.900) | 9 to 10 |
| trpv1, converged | 0.540 | 0.771 | **0.964** (0.969, 0.953, 0.969) | 8 |
| bgal, `center_ptcls` | 0.311 | 0.760 | **0.726** (0.714, 0.716, 0.750) | 6 to 8 |
| bgal, converged | 0.551 | 0.833 | **0.886** (0.901, 0.891, 0.868) | 5 |

- The verdict (0.95, 0.30, at most `KMAX_CLASS`) fails on both sets with `center_ptcls`. **trpv1 meets the stop
  condition of section 7 exactly**: it passes at converged shifts and fails on the same particles centred by
  `center_ptcls`, at every seed.
- bgal fails at converged shifts too against the dense pass's mode. Against the teacher's own rows at those shifts
  (the learn program's fill with the 2 Å prior, as the held-out rows are made; seed 1 only, `runs/p3_validate_bgal_s1`
  before the rerun) the same student reached 0.99 converged and 0.83 centred. The dense pass's 30 Å disk lets a
  poorly centred particle trade class against shift, so its mode is a harder yardstick than the teacher's rows.
- The gap is the centring: `center_ptcls` leaves the unseen particles 2.8 to 3.6 Å (trpv1) and 5.8 to 7.7 Å (bgal)
  from the dense pass's shifts (median), against the student's 2 Å training disk. Against the teacher's shifts on the
  training particles its common translation is below 0.4 Å and the residual after it 3.0 to 3.7 Å (trpv1) and 5.2 to
  7.5 Å (bgal).
- Where the student's class agrees with the dense pass, its read-out angle is within 1.05 to 1.19° and its
  registered shift within 0.22 to 0.32 Å (medians).
- Class averages from the student's assignments alone against the teacher's (median FRC 0.5 over the classes with
  at least 10 particles on both routes): trpv1 22.2 Å against the dense pass's 18.3 to 20.8 Å; bgal 14.8 to 20.4 Å
  against 15.5 to 19.2 Å.

**Times** (Release, 24 threads; the machine carried other work, load 15 to 20, so these are indicative):

| phase | bgal | trpv1 |
| --- | --- | --- |
| solve2D on the training set | 32 to 35 s | 138 to 139 s |
| `learn_class_student` (rows and two designs) | 43 to 49 s | 185 to 264 s |
| `center_ptcls` on every particle | 11 s | 34 to 45 s |
| `predict_classes` on the unseen (with registration) | 1.6 s (0.56 ms each) | 10.6 to 15.4 s (0.41 to 0.60 ms each) |
| the dense pass on the unseen (full box) | 3.3 s (1.18 ms each) | 29 s (1.14 ms each) |

At the full box the dense pass against 50 to 100 classes costs only about twice the student's prediction.

**Not run, because the verdict stopped the phase:** the six-seed self-check of the comparison rule, the
`trpv1_full` converged-shift arm (the validation's own converged arm stands in for it), the reads on `bgal_s2`.

**The choice for Hans** (section 7): a better centring (for example a second `center_ptcls` pass, or one against
the teacher's own class averages before prediction), or the future note's shift-robust front end. Run directories:
`/data2/agent_runs/class_student/runs/p3_learn/`, `runs/p3_validate_{bgal,trpv1}_s{1,2,3}/`, the workflow test's
`runs/p3_wftest_{release,debug}/`.

**Ruling 15 (2026-11-14).** Hans inspected the bgal class averages of the student's assignments and judged them good
enough; Phase 3 passes on his inspection (decision 15). The proposed numeric rule failed: the student's median
FRC 0.5 against the teacher's lies 1.00, -0.99 and 4.49 Fourier shells behind the dense pass's on bgal (mean 1.50) and
1.48, 3.00 and 1.00 on trpv1 (mean 1.83). With about 55 unseen particles per class on bgal, one seed moves that read
by 1.5 shells.

### Phase 4, the 2D consumer and solve2D at scale (2026-11-14, the Dell)

**What was built** (uncommitted):
- `src/main/project/simple_class_candidates.f90` (new; its file-table row added before the edit) and its tester
  (suite `class_candidates` in `unit_project`, 10 checks). The class `class_candidates` writes the candidates of the
  predicted rows (route, count, classes in descending mass) to a stream file under the particles'
  `ptcl_layout_digest`, registers it in the `out` segment (imgkind `class_candidates`, through `add_entry2os_out`, as
  the flex gateway does; `sp_project` is unchanged), and opens it only when the digest, the row count and the class
  count match the consumer's (statuses none, unreadable, stale).
- `predict_classes` publishes it beside the project (`<fname>_cands.bin`).
- `eul_prob_tab2D`: the sparse table (`refine=prob_neigh`) opens the artifact at construction. A seeded particle
  scores its candidates alone, in descending mass, and searches the window instead of the disk, because its shift
  was registered against its predicted class. A refused artifact is reported and the table runs as without it.

**Gates.**
- The stale digest: the test of another particle layout is mutation-tested (with the digest comparison removed, 2 of
  the suite's 10 checks fail; restored, 10 pass).
- **Off-state identity** (no artifact): solve2D on a 1,000-particle bgal selection, 20 classes, at one thread, the
  snapshot before Phase 1 (`build_ref`) against the tree. All 82 class-average files are byte-identical, and all 27
  cluster tables identical apart from the run directory's name, through the `prob_neigh` stages 3 to 5. At 24
  threads the two builds differ from iteration 1, as solve2D's first stage differs between runs of one build (its
  draws are per thread); the test policy judges identity at one thread.

**solve2D at scale on trpv1** (35,645 particles, 100 classes, Release copy of the tree, 24 threads, seeds 1 and 2;
the machine carried other work at a load of 20 to 25 throughout, so the times are indicative; `tools/routes2d.sh`).
Route A is solve2D on every particle. Route B is solve2D on a random 10,000 (the teacher), `learn_class_student`,
`center_ptcls` and `predict_classes` on the other 25,645, then the last sparse stage's settings (`prob_neigh`, lp
8 Å, box 120, `trs` 4.6, `ml_reg`, `extr_iter` past the limit) for 5 iterations on every particle, seeded by the
artifact, from the teacher's last averages at the crop (its iteration 29; its final restoration is full-size).

| | route A, seed 1 | route A, seed 2 | route B, seed 1 | route B, seed 2 |
| --- | --- | --- | --- | --- |
| solve2D (A: all; B: the 10,000) | 453 s | 468 s | 133 s | 140 s |
| `learn_class_student` | | | 240 s | 205 s |
| `center_ptcls` (25,645) | | | 37 s | 37 s |
| `predict_classes` (25,645) | | | 9 s | 9 s |
| the seeded sparse stage (5 iterations, 35,645) | | | 57 s | 58 s |
| **total** | **453 s** | **468 s** | **477 s** | **450 s** |
| seeded by the artifact | | | 17,327 (49%) | 17,979 (50%) |
| classes accepted by `model_cavgs_rejection` (defaults) | 22 | 23 | 26 | 25 |
| mean class FSC 0.143 of the last iteration | 9.7 Å | 9.6 Å | 8.7 Å | 8.3 Å |

- **At 35,645 particles the routes cost the same.** Training is route B's largest step (205 to 240 s); it does not
  grow with the data set. Corrected 2026-11-14: the teacher's row fill takes 3.0 s (0.30 ms per particle) and the
  learn program about 7 s before training (`runs/p4_routeB_s1/learn.log`); the designs take the rest, width 32 about
  three quarters. The 45 s written here first has no source in the logs.
- **Seeding barely shortens an iteration of the sparse stage:** 11.4 s, against 11 to 13 s for solve2D's own sparse
  stages (Phase 0). Phase 0 found about 8 s of each iteration in reading, preparing and summing particles, which the
  student does not touch.
- **What grows with the particles** in route B: centring 1.4 ms, prediction 0.35 ms and the sparse stage 1.6 ms (5
  iterations) per particle, about 3.4 ms, against route A's 13 ms per particle. A straight extrapolation (not measured,
  and blind to solve2D's own sampling at scale) puts a million particles at about 1 h for route B against 3.6 h for
  route A. The maintainer's runs at real scale (decision 13) are what can show it.
- **Quality:** route B's class averages pass the quality model in about the same numbers (25 to 26 of 100, route A 22
  to 23). Their mean FSC is better, but the iterations differ (route A ends with the dense full-size stage), so that
  read is not like for like. Montages for inspection: route A `runs/p4_routeA_s{1,2}/cavgs_iter030.mrc`, route B
  `runs/p4_routeB_s{1,2}/sparse/cavgs_iter006.mrc`; the accepted classes `runs/p4_quality_{A,B}{1,2}/
  quality_selected_cavgs.mrc`.
- Prediction seeded 68 and 70% of the unseen particles (mean 10.7 candidates). The other half of all particles is
  not seeded: the teacher's 10,000 carry no artifact entry, and unseen particles outside the cap or under the
  threshold search as before.

**Faults on the way (the run script's, not the programs').** The first route B runs centred in the run directory,
where `center_ptcls` read solve2D's full-box `frcs.bin` and stopped, and gave the sparse stage the full-size
averages, which the class averager cannot crop. Both steps now run in their own directories.

**Exit:** the off-state identity holds; route B's numbers stand beside route A's, two runs each; the stale-digest test
is mutation-tested. Not done: route B with a dense last stage (as route A ends), and the runs at real scale.

### Close (2026-11-14, decision 17)

Phases 0 to 4 are delivered, uncommitted for Hans's review. Phases 5 to 7 are not started; they are replaced by
`planned/class_student_stream.md`, whose section 3 records two things the code check found for them: the pool runs
`refine=snhc_smpl`, which the candidate artifact does not reach (only `prob_neigh` reads it), and the artifact's
layout digest changes whenever the pool appends particles. The algorithm page `doc/algorithms/class_student.md` is
written (Phase 7's first item); the policy updates of Phase 7 belong to the stream plan.

### After the close (2026-11-14)

- **Ruling 18 (Hans): the class student trains width 16 alone.** `CLASS_STUDENT_NDESIGNS` = 1
  (`src/defs/simple_defs_hhcnn.f90`) caps the designs `learn_class_student` tries; raising it restores the resampled
  choice unchanged. The record behind it, every real-teacher learn of this run: width 32 kept in 2 of 6 trpv1 learns
  (`p3_validate_trpv1_s1`, `_s3`: +0.018 and +0.016 held-out classes recovered), never on bgal (6 learns), at 2 to 4
  times width 16's training time. The pose student keeps its choice of three designs.
- **The validation report answers one question** (Hans): can the student replace the full search? Five numbers
  (the seeded share; the full-search class among a seeded particle's candidates, against the calibration target;
  the candidates per seeded particle; the class averages against the teacher's, by both routes; the cost per
  particle, by both routes), the same safety at ideal shifts, how to read them, and the run time. Agreement on
  the single most probable class, the angle and shift errors and the frame offset moved to
  `validate_class_student_details.txt`. In what users read, "dense pass" became "full search" and a seeded
  particle's candidates its "shortlist", each line with its goal and a sentence on what failed and why.
- **Training stays in one process** (Hans): the teacher's rows are 3% of the 2D learn and 15% of the 3D one, and the
  designs' step is too short to exchange gradients between nodes through files. Prediction in parts is delivered in
  `completed/student_predict_parts.md`.
