# The class student in the stream: the pool's learn job, seeding and candidates

Status, 2026-11-14: planned. The owners are Joe and Cyril; Hans approves the plan before any code. It continues
`completed/class_student_2D.md` (Phases 0 to 4 delivered there; its Phases 5 to 7 are replaced by this plan).
Base: the commit that merges the 2D student (HEAD `38e033cbf` plus that plan's working-tree diff); pin its hash
here when it lands. Depends on:
- the class student and its programs: `learn_class_student`, `predict_classes`, `validate_class_student`, the
  candidate artifact `class_candidates` (`doc/algorithms/class_student.md`; validated: yes on bgal and trpv1,
  three seeds, Hans's ruling 15 of the 2D plan);
- the stream pool (p06) and its policy `doc/policies/stream/pool2D_policy.md` (validated: yes, in production).

## 1. Why

The pool (p06) runs a distributed refine2D iteration over a sample of up to 500,000 particles
(`STREAM_NPTCLS_MAX`) against its classes, and gives each new particle a random populated class before its first
search (`draw_new_classes`). A class student trained on the pool's own class averages could give new particles a
predicted class, rotation and shift instead, and restrict their search to a few candidate classes. The 2D plan
built the student and measured it (section 2). Whether it pays in the pool is not measured: the pool's time per
iteration has not been broken down, and the stream test is too small to show it (2D plan, Phase 0).

## 2. What the 2D run established (the facts this plan builds on)

All on the Dell, Release builds, 24 threads; the full tables are in the 2D plan's Progress.

1. **The student needs centred particles.** Trained over a 2 Å shift disk (`learn_class_student`'s default), it
   recovered the teacher's class among its candidates for 94 to 98% of held-out bgal particles and 82% on trpv1,
   with the in-plane angle within about 1°. At refine2D's default 10 Å disk it learned almost nothing (6%).
2. **Centring decides the result.** On unseen particles centred by `center_ptcls` (residual 3 to 4 Å on trpv1, 6 to
   8 Å on bgal) the dense search's class was among a seeded particle's candidates for 0.90 (trpv1) and 0.73 (bgal);
   at converged shifts for 0.96 and 0.89 (three seeds each). The pool's particles arrive with the sieve chunks' 2D
   shifts, which should be the good case; it is not measured (Phase 0).
3. **The exact-class read is too strict on redundant classes.** bgal's D2 views give near-duplicate classes, between
   which the read counts a miss. Hans judged the student's class averages good enough by eye (ruling 15). A numeric
   rule on the class averages' FRC failed and is noisy: about 55 particles per class move it by 1.5 shells between
   seeds.
4. **Training is the fixed cost.** 205 to 240 s for 10,000 trpv1 particles and 100 classes with two designs of the
   network; the teacher's row fill is 3 s of it, width 32 about three quarters. Width 32 was kept in 2 of 6 trpv1
   learns (gains under 0.02) and never on bgal, so the class student now trains width 16 alone
   (`CLASS_STUDENT_NDESIGNS`, ruling 18 of the 2D plan).
5. **Seeding barely shortens an iteration.** In refine2D `prob_neigh` on trpv1 an iteration took 11.4 s with half
   the particles seeded, against 11 to 13 s without; about 8 s of it is reading, preparing and summing particles.
   The gain at scale is the per-particle cost (about 3.4 ms against 13 ms in solve2D), not the iteration.
6. **Prediction is cheap:** the features, the forward pass and the shift registration take 0.1 to 0.3 ms per
   particle; the whole program, with reading and preparing the particles, 0.35 to 0.6 ms.
7. **solve2D at scale, route against route** (35,645 trpv1 particles): the same time (450 to 477 s against 453 to
   468 s), class averages accepted by the quality model in about the same numbers (25 to 26 against 22 to 23 of 100).

## 3. What the code check changed (at `38e033cbf` with the 2D diff)

- **The pool does not run the search the artifact feeds.** `stream_pool2D%init_state` accepts only `refine=snhc` or
  `snhc_smpl`, and p06 sets `snhc_smpl`. The candidate artifact is read by `eul_prob_tab2D`'s sparse fill
  (`refine=prob_neigh`) alone. Seeding the pool's search therefore needs either the `snhc_smpl` strategy to read the
  candidates, or the pool to run `prob_neigh` (decision 4).
- **The artifact goes stale as the pool grows.** `class_candidates` is tied to the project's rows by
  `ptcl_layout_digest`; every `append_sets` adds rows, so an artifact written before an append is refused (by
  design) after it. It has to be republished for the iteration's project, or its tie changed (decision 5).
- **The pool runs distributed** (`refine2D_distr` through `simple_private_exec`, workers `refine2D`), with
  `ml_reg=no`. The workers read the project file, so the artifact must be registered in the iteration's project
  before the job starts.
- `draw_new_classes(spproj, nptcls, iter, ncls)` in `simple_stream_refine2D_utils.f90` sets the class of every
  selected row whose `updatecnt` is 0 to a populated class drawn with a seed of the iteration. `FIRST_EXPORT_ITER`
  is 10 (`MSKDIAM_SWITCH_ITER`).
- The class averager pads smaller averages up but does not crop: references must be at the consumer's box.
  `learn_class_student` and `predict_classes` work at the full box with full-size averages; the pool works cropped.

## 4. Contract

- **Off state.** Without a model, and without an artifact, the stream behaves byte for byte as on the base at a
  fixed seed (one thread for the identity runs; the test policy's section 1.5).
- **The pool owns its project.** The learn and predict jobs read the pool's files and write their own; the pool
  writes the project (`doc/policies/layering_and_project_ownership_policy.md`, section 5; the stream's IPC and
  restart policies).
- **No new command-line key**, as in the 2D plan (its decision 12): a model or an artifact switches the route.
- **The 2D student's programs are the building blocks.** Extend `learn_class_student_exec` and
  `predict_classes_exec`; do not copy their frames.

## 5. Decisions for the owners to close before a run

Each needs an answer (or a measured Phase 0 read) before the phase that uses it. Defaults are proposals.

1. **When the pool learns.** Proposal: at iteration `FIRST_EXPORT_ITER` (10), as the 2D plan's ruling 14 set, on a
   random sample of the pool's selected particles (proposal 10,000, the trpv1 teacher's size), against the pool's
   current class averages at the full box.
2. **Width.** Proposal: train width 16 only (fact 4), halving the cost; keep the design loop for a measured case.
3. **How a new particle is seeded.** Proposal: when a model exists, `append_sets` predicts the new rows (class,
   rotation, registered shift) instead of `draw_new_classes`; rows the student does not seed keep the random draw.
4. **Candidates in the pool's search.** (a) the `snhc_smpl` strategy scores a seeded particle's candidates first
   and stops there, or (b) the pool moves to `prob_neigh`. (b) is a change to the pool's search on its own and
   needs its own measurement; (a) keeps the pool's search. Proposal: (a), or no candidates at all if Phase 0
   shows the search is a small share of the iteration (fact 5).
5. **The artifact's tie to the rows.** Republish per iteration from the stored predictions (cheap), or tie it to
   the stack index and `indstk` of each row instead of the whole layout. Proposal: republish.
6. **Keeping the student current** as the classes change: re-template (the new averages into the readout, a short
   fine-tune) or retrain. Both measured as variants; the cheaper one at equal reads is kept.
7. **The acceptance rule**, frozen in Phase 0 from measured spread (fact 3 shows the exact-class read is the wrong
   one): proposal, the published class averages' `cavg_quality` acceptance and the time to each publication, beside
   the base's, over at least three stream runs each.
8. **The sieve (the 2D plan's Phase 6).** Only if Phase 0 shows the sieve's chunk solve2D holds at least a third of
   the stream's 2D time.

## 6. Not in this plan

- The student itself (its network, features, training, model file): changes go through `doc/algorithms/
  class_student.md` and the 2D student's suites.
- solve2D at scale (the 2D plan's route B is measured and stays a script, `tools/routes2d.sh` in the run folder).
- The stream's restart and IPC contracts beyond what the learn job needs; the 3D ingestion policy.
- Unsupervised or generative 2D (`planned/polar_student_future_directions.md`).

## 7. Phases

| Phase | Kind | Gate |
| --- | --- | --- |
| 0 | measurement | the stream's 2D time breakdown on a real stream; the pool particles' shift residual; base stream tests twice |
| 1 | identical until a model exists | the pool's learn job; off-state identity; the model's held-out reads on pool data |
| 2 | result-changing once a model exists | seeding by prediction; the acceptance rule of decision 7 |
| 3 | result-changing once a model exists | candidates in the pool's search (decision 4); as Phase 2 |
| 4 | result-changing | keeping the student current (decision 6); as Phase 2 |
| 5 | conditional | the sieve (decision 8) |
| 6 | close | docs, policies, report |

- **Phase 0, the base.**
  - A real stream (Hans's data, the 2D plan's decision 13), and the registered stream tests (`mini_stream_*`,
    `stream_chains`) twice on the base.
  - The breakdown of a pool iteration: read and preparation, the search, the class sums, the job overhead; and the
    sieve's chunk solve2D time against the pool's.
  - The pool particles' shifts at entry against their converged shifts after N pool iterations (the residual of
    fact 2, measured where it applies).
  - Offline: `validate_class_student`'s learn and predict steps on a snapshot of a pool's project (its class
    averages as the teacher), for the held-out reads on pool data.
  - Exit: the breakdown, the residual and the offline reads in Progress; decisions 4, 7 and 8 closed with them.
- **Phase 1, the learn job.** The pool starts `learn_class_student` as an asynchronous job at the iteration of
  decision 1, on its files, writing the model in the pool's folder; the stage records the model in its state and
  restart files. Exit: the off-state identity; the stream tests pass; the model's held-out reads beside Phase 0's.
- **Phase 2, seeding.** Decision 3. Exit: the acceptance rule against the base; the off-state identity.
- **Phase 3, candidates.** Decisions 4 and 5. Exit: as Phase 2, and the stale-artifact case tested.
- **Phase 4, keeping current.** Decision 6. Exit: as Phase 2, both variants measured.
- **Phase 5, the sieve.** Only under decision 8.
- **Phase 6, close.** `doc/algorithms/class_student.md` and `pool2D_policy.md` updated, the report, this plan to
  `completed/`, the gates on both BLAS builds.

## 8. Stop conditions and gates

Stop when:
- Phase 0 shows the pool particles' residual well beyond the student's 2 Å disk (fact 2) and no centring step fits
  the stream's budget: the student would not seed reliably;
- the search is a small share of a pool iteration and seeding by prediction (Phase 2) does not change the published
  averages: there is nothing left to gain;
- a phase needs a command-line key, or a write of the pool's project outside the pool;
- the off-state identity fails.

Gates no ruling may waive: the fast gate; the student suites (`polar_features`, `hhcnn`, `polar_student`,
`class_student`, `student_strategy`, `class_candidates`); the stream workflow tests; CTest `class_student`; the
off-state identity.

## 9. Practical notes from the 2D run

- `build/` on the Dell was configured Debug for a while; take every time on a Release build (a separate build tree
  is safest).
- `center_ptcls` and refine2D runs write `cavgs_iter*.mrc` and `frcs.bin` into their working directory: give each
  step its own directory when several share a project, or one reads another's files.
- The dense search draws classes from the posterior; for a yardstick take each row's mode (`prob_beta` large, as
  `validate_class_student` does), or the read counts the draws' noise.
- The student's training is not thread-identical at the program level (a pre-existing finding of the 2D plan's
  Phase 1, likely a BLAS call outside a parallel region); identity runs are at one thread.
- The runs, scripts and logs of the 2D plan are in `/data2/agent_runs/class_student/` on the Dell (`tools/`,
  `runs/`).

## 10. File table

| File | What | Phases |
| --- | --- | --- |
| `src/main/stream/pool2D/simple_stream_pool2D.f90`, `simple_stream_refine2D_utils.f90`, `simple_stream_pool2D_tester.f90` | the learn job, seeding, the artifact per iteration | 1-4 |
| `src/main/stream/stages/simple_stream_stage_pool2D.f90` | the model in the stage's state and restart | 1-4 |
| `src/main/strategies/parallelization/simple_strategy_class_student.f90` | learn and predict on the pool's files, if needed | 1-3 |
| `src/main/strategies/search/simple_strategy2D_snhc_smpl.f90` (or the `prob_neigh` route) | candidates in the pool's search (decision 4) | 3 |
| `src/main/project/simple_class_candidates.f90` | the artifact's tie (decision 5), if needed | 3 |
| `src/main/sieve/` | fine-chunk assignment | 5 (conditional) |
| `doc/algorithms/class_student.md`, `doc/policies/stream/pool2D_policy.md`, `doc/policies/sieving_and_rejection/ptcl_sieve_policy.md` | the living rules | 6 |
| this plan, `doc/implementation_notes/README.md` | Progress, the close | 0-6 |

## Progress

(none yet)
