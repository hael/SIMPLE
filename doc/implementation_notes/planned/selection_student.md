# The selection student: a particle-level challenger to the class-average rejection

Status, 2026-11-15: planned, awaiting Hans's approval. Base: `a4881e6c4` (the working tree also holds uncommitted
refine3D/solve3D parameter edits that this plan does not touch). Supersedes: none. Depends on: the polar students'
core, `simple_polar_student` and `simple_hhcnn` (validated: yes, per data set; `doc/algorithms/class_student.md`,
`doc/algorithms/pose_student.md`); the class-average rejection `model_cavgs_rejection` with its `pool` preset
(validated: yes as a class-level model on its training tables; as a reproducible per-particle label, no: Phase 0
measures it); solve2D (validated: yes).

## In one paragraph

The class student already sorts particles into class averages, and it keeps rejected classes as classes. The
selection student is the same network asked a simpler question: will this particle end up in a selected class? Its
classes are the class averages of several training data sets, each marked selected or rejected. Its target is the
mark of the class each training particle ended up in. Its answer, `psel`, is the posterior mass on the selected
classes. The existing machinery lacks two things: a loss on the mass of the target's group instead of on a
distribution over classes (one option in the network), and a frame that is the same on every data set (the same
band and half ring, each data set whitened by its own spectrum). The challenge is pre-registered. On a data set it
has never seen, does `psel` predict the class-average machinery's verdict nearly as well as the machinery reproduces
its own verdict from another seed?

## 1. Why

Rejection in SIMPLE is decided on class averages. The sieve's and the pool's `model_cavgs_rejection` score each
class, and every particle inherits its class's state (`map_cavgs_selection`). This has two costs.

- **A class-level rule is capped by the purity of the classes.** The nanoparticle plan defines that cap as what a
  perfect class gate can do at a given class population (`planned/trajectory_cavgs_rejection.md`, section 3).
  Junk inside accepted classes, and good particles inside rejected ones, are out of its reach.
- **The verdict needs a 2D run first.** The sieve runs a whole solve2D per chunk before the pool sees a particle.

The 2D student showed that the HHCNN learns a class assignment from a teacher in 205 to 240 s on 10,000 particles
and predicts at 0.1 to 0.3 ms per particle (`doc/algorithms/class_student.md`). The per-data-set students do not
transfer by construction, because their classes are one data set's averages. Whether one network can learn the
retention criterion across data sets is open. This plan measures it against the class machinery's own
reproducibility before any integration is considered.

## 2. Contract

- **A challenger only.** No program, preset or default of the class-average rejection changes. The following are
  out of bounds: `model_cavgs_rejection`, `src/main/cavg_quality/`, `src/main/quality/`, the sieve, the stream,
  solve2D and refine2D. `predict_selection` reads a project and writes a table. It never changes a state, a record
  or the project.
- **Off state.** The network's new loss option is off by default. With it off, `hhcnn`, `polar_student`,
  `class_student` and `pose_student` give byte-identical models, weights, tables and records at one thread and a
  fixed seed (Phase 1's gate).
- **Selection confidence is its own number.** `psel` never enters `pconf`, `route` or the class student's
  candidates.
- **One copy of everything shared.** The network, the training loop, the held-out split, the shift augmentation,
  the model file's common part, the part ranges and the distributed master stay where they are: in the core and in
  the students' common strategy module. The selection student composes a `polar_student`, as the class and pose
  students do (`completed/class_student_2D.md`, decision 4).
- **The library's object rules and the layering policy hold.** The student is an encapsulated class with `new`
  and `kill`. Its commanders are thin `execute(cline)` wrappers of strategy routines, `run_<x>(cline)` and
  `<x>_exec(params, build, ...)`. No new command-line parameter is added.

## 3. What the code check found (at `a4881e6c4`)

- **The core trains on soft class targets only.** `hhcnn%backward_sample` sets `dlogits = prob - t`, and
  `sample_loss` is the cross-entropy against `t`. A target that names a group G of classes needs the loss
  −log Σ_{c∈G} p_c. Its logit gradient is p_c − p_c·[c∈G] / Σ_{c∈G} p_c: the same backward pass, with `t` replaced
  by the posterior renormalised on the group. This is the one change to the network (decision 3).
- **`polar_feats` whitens a tensor by shell weights frozen in the model** (wk = k/σ²(k), `feats%new(..., wk)`).
  The class and pose students take them from the teacher's calculator (`build%pftc%get_sigma2_weights`). A
  cross-dataset model cannot hold one data set's weights. The selection student's tensors therefore carry their own
  data set's whitening, and its core holds flat weights (decision 6).
- **Pinning `pftsz` and `kfromto` pins the network's shape.**
  - The calculator's half ring comes from the mask (`magic_pftsz` in `simple_parameters_phases`) unless `pftsz` is
    given on the command line. `box_crop` is a command-line key too.
  - A Fourier index k means the resolution box·smpd/k at any crop, so a fixed `kfromto` is the same fraction of the
    box on every data set.
  - The network's half ring is min(pftsz, 2⌈π k_max⌉) (`polar_feats_net_pftsz`).

  See decision 5.
- **Particles are prepared at their stored 2D shifts** (`prepimg4align`, `shvec`). A project copy with the shifts
  zeroed gives the particles at their picks, unless the stacks were re-extracted at the shifts. Phase 0 records
  which.
- **The prediction can reuse `predict_classes`' path without references.** That path builds the calculator with one
  reference slot (`build%pftc%new(params, 1, ...)`) and needs no sigma2.
  - Reused unchanged from `simple_strategy_student_common`: `student_holdout_split`, `student_part_pinds`,
    `student_predict_parts`, `concat_part_tables`, `fname_base` and `file_hash`.
  - Not reused: `student_geometry_mismatch`. The geometry varies by design, so the frame is checked instead.
- **Particle records are a static enumeration** (`simple_defs_ori`: slots 1 to 59 named, 60 to 64 spare). A `psel`
  slot would change the record layout. In this plan the prediction is a table (decision 10).
- **The HHCNN, polar-student and class-student suites run in `unit_numerics`** (`suites_numerics` in
  `simple_commanders_test_class`). The highlevel gates are routed in `simple_test_exec_highlevel` and declared in
  `simple_commanders_test_highlevel`.
- **Two planned notes also edit `hhcnn` and `polar_student`:** `planned/student_inplane_posterior.md` and
  `planned/pose_student_amplitude_channel.md`. Both keep their variants outside the tree until a verdict. Whichever
  of them or this plan lands its `hhcnn` edit second rebases on the first. The marginal loss acts on the logits,
  whatever readout computes them. Phase 0 depends on neither.

## 4. Decisions and rulings

1. (2026-11-15, Hans) **A particle-level challenger to the class-average rejection:** a network trained to
   recognise the particle images that end up in selected classes. Reuse as much of the students' code as possible
   and keep the architecture clean. It is simpler than the same day's proposals of a solve2D posterior readout and
   label tiers, which this decision replaces.
2. **The label is the state of the particle's class.** A training particle's label is 1 when the `cls2D` state of
   its `ptcl2D` class is positive, and 0 when it is zero. Particles without a class are not used. The programs read
   whatever set the states. The experiment sets them by one pipeline (decision 9), which gives four properties: the
   label is a defined function of the data, it is available at scale, it is free of selection by intent, and its
   reproducibility can be measured.
3. **The loss is the mass on the label's group.**
   - `hhcnn%new` and `polar_student%new` gain an optional logical `marginal`, default false.
   - With it on, the targets are a 0/1 support over the classes. A sample's loss is −log Σ_c t_c p_c (in training
     and in `eval_batch`), and its logit gradient is p − p∘t / Σ(p∘t).
   - The mode count counts a sample whose largest class lies in its support.
   - Prediction is unchanged. The weights file does not record the option, since it does not change the network.
4. **The classes are a bank of the training sets' classes.**
   - Membership: every class of every training set with at least `SELECTION_BANK_MINPOP` = 10 training particles
     (guessed; the value of `CLASS_MINPOP_FRC`) is a class of the student, in the group of its state.
   - Starting templates: each member's tensor is turned into its class's frame by its stored in-plane angle, at the
     nearest grid rotation, with `polar_feats%rotate` (the inverse of the class student's read-out
     e3 = 360 − ang). The turned tensors are averaged, and `init_templates` takes the averages, as the class
     student's takes its class averages.
   - `psel` is the posterior mass on the selected group. There is no reflection, since this is 2D.
5. **The frame** (all three values guessed). Every data set is prepared with `box_crop` = `SELECTION_BOX_CROP` = 64
   (or the box when smaller), `pftsz` = `SELECTION_PFTSZ` = 96 and `kfromto` = `SELECTION_KFROMTO` = (2, 15).
   - The band runs from half the box to a fifteenth of it.
   - The half ring is the network's cap at k = 15 (2⌈15π⌉ = 96, an FFT-friendly size), so nothing is resampled.
   - The commanders put `box_crop` and `pftsz` on the command line before `params%new`. The strategy sets
     `kfromto`, as `predict_classes` sets the model's.
   - Each data set's own mask diameter and CTF are used. `with_ctf` must agree between the model and the data.
6. **Each data set's whitening is its own.** wk(k) = k / p(k), where p(k) is the mean per-sample power at shell k of
   flat-weighted, per-particle-normalised tensors.
   - The sample: `SELECTION_WK_NSAMPLE` = 2,000 particles (guessed), at a fixed stride through the data set's
     selection.
   - wk is normalised to mean 1 over the band and is the same for both halves.
   - One routine serves learning and prediction. Every part of a distributed prediction computes wk from the same
     particles.
   - The model's core holds flat weights.
7. **Training.**
   - At most `SELECTION_NROWS_DATASET` = 10,000 labelled particles per training set (guessed; the 2D student's
     teacher size), drawn at a fixed stride.
   - Every tenth stack of each set is held out (`student_holdout_split`), for early stopping and the thresholds.
   - Width 16 alone, as for the class student (`completed/class_student_2D.md`, ruling 18). The core's budget
     applies (`STUDENT_MAX_EPOCHS`, `STUDENT_PATIENCE_EPOCHS`).
   - The shift augmentation's disk has a radius of `SELECTION_SHIFT_AUG_FRAC` × `SELECTION_BOX_CROP` frame pixels.
     The core's features take the frame box, so a shift is in frame pixels. The fraction is measured and frozen in
     Phase 0 (decision 9).
   - Sets and labels are not re-weighted beyond the cap.
8. **One probability, two thresholds, three outcomes.** Both thresholds are set on the held-out rows.
   - t_rej is the largest threshold at which at least `SELECTION_RECALL` = 0.95 of the selected particles of every
     training set have psel ≥ t_rej (the recall-first rule of the class-average learners).
   - t_acc is the smallest threshold at which at least `SELECTION_PRECISION` = 0.95 (guessed) of the held-out
     particles at or above it are selected, on every training set. When no threshold qualifies, t_acc is set above
     1, so that no particle is accepted outright.
   - A particle is accepted at psel ≥ t_acc and rejected below t_rej. Otherwise it is sent to collective analysis
     (a 2D of the uncertain particles; a later plan). When t_acc ≤ t_rej, the middle band is empty.
9. **The labels and the control come from the class machinery itself.**
   - Each data set (at most 100,000 particles, with a keyed random subset above that; guessed, for cost) goes through
     solve2D at its defaults and then `model_cavgs_rejection quality_model=pool prune=no`, twice, at `SIMPLE_SEED` 1
     and 2.
   - solve2D runs with `ncls` = 100, or 50 below 10,000 particles (guessed).
   - The preset is `pool`, the pre-3D cleanup of pooled class averages, because that is the closest context to a
     whole data set's 2D.
   - Seed 1's states are the labels. Seed 2's states, compared with seed 1's on the particles both runs labelled,
     measure the pipeline's own reproducibility: the control of the verdict (section 7).
   - Phase 0 also sets `SELECTION_SHIFT_AUG_FRAC`: the largest, over the data sets extracted at their picks, of the
     90th percentile of |stored shift| / box over the seed-1 particles.
10. **Programs.**
    - `learn_selection_student` (private): a `filetab` of labelled projects in; the model (`fname`) and its learn
      report out.
    - `predict_selection` (private): a project and the model (`infile`) in; `<fname>_pred.txt` out, one row per
      particle with iptcl, eo, psel, the outcome, and the label (the state of the particle's class, −1 without one).
      The project is never written. With `nparts`, it runs through the 2D student's distributed master
      (`student_predict_parts`); the parts write their records unchanged, because the master's merge expects them.
    - `validate_selection_student` (`simple_exec`, developer visibility): the leave-one-data-set-out read of
      section 7, and its report.
    - The training set is assembled in memory, one project at a time, each project with its own parameters and
      builder (as `validate_class_student` chains its programs). There is no feature cache.
11. **Data.**
    - At least `SELECTION_MIN_DATASETS` = 6 distinct specimens (guessed), each with particle stacks, CTF parameters
      and at least 5,000 particles, copied to the Dell.
    - At least four of them must have at least 5% of their particles rejected by seed 1's pipeline; otherwise the
      student sees few negatives.
    - The list is Hans's (open item A). bgal and trpv1 on the Dell may be two of them. The simulated sets are not
      eligible.

**Open items for Hans before launch.** The proposals above are in force until he rules.

- A: the data sets (decision 11).
- B: the preset that labels: `pool` or `chunk100mics` (decision 9).
- C: the frame (decision 5).
- D: whether the workflow test also becomes a CTest entry (budget 46 → 47). The default is no.

**Defaults taken without a question:**

- the names: `selection_student`, `psel`, the three programs, and the `SELECTION_` constants in
  `simple_defs_hhcnn`;
- the table's outcome codes: 1 accept, 0 reject, 2 collective analysis;
- the validation's data-set directories, each holding `labels_seed1.simple` and `labels_seed2.simple`;
- one learn seed per fold (`SIMPLE_SEED` 1) for the verdict.

## 5. Not in this plan

- Any change to the class-average rejection (`model_cavgs_rejection`, `cavg_quality`, `quality`), the sieve, the
  stream, solve2D, refine2D, or the behaviour of the class or pose student.
- Integration as a pre-filter (its natural owner is the particle sieve, `ptcl_sieve_policy.md`), the
  collective-analysis lane, and a `psel` slot in the particle record. A positive verdict opens that plan; this plan
  ends at the read.
- Manual selections as labels. The programs read the states as they are, so an experiment with manual labels needs
  no code; it is a later run. Also excluded: a class posterior read out from solve2D, label tiers, and the
  posterior-weighted baselines.
- Re-boxing data sets to their particle size. The ratio of box to particle is a variable the network sees; Phase 0
  records it per data set.
- Runs on the Mac. Every run is on the Dell; the Mac runs the unit suites only.

## 6. Change types and gates

| Phase | Kind | Gate |
| --- | --- | --- |
| 0 | measurement | builds, fast gate, the identity references, the labelled data sets and the control |
| 1 | identical for the existing students | the class and pose students' outputs byte-identical at one thread; the HHCNN suite with the new gradient, toy and thread checks |
| 2 | new code; tests-only for existing paths | builds, fast gate, `unit_numerics`, `unit_ui`, the registry and description checks, the workflow test twice, parts byte-identical to shared memory |
| 3 | scientific read | the pre-registered verdict (section 7) |
| 4 | report | the disagreement averages collected; the report |

## 7. Acceptance

The read is leave-one-data-set-out. For each data set H in turn, the student learns on the others and predicts
every labelled particle of H, not just the training cap.

**Definitions.** Per H, write v for the seed-1 label, u for the seed-2 label, s for psel and t for the fold's t_rej.

- The student's recall R is the share of particles with v = 1 that have s ≥ t. Its removal Q is the share with
  v = 0 that have s < t. Its balanced accuracy is BA = (R + Q)/2.
- The control comes from the pipeline: R_pipe is the share of v = 1 with u = 1, Q_pipe the share of v = 0 with
  u = 0, and BA_pipe = (R_pipe + Q_pipe)/2 (measured in Phase 0).
- Reported beside: the AUC of s against v; the shares of the three outcomes and the labels within each; and BA_in,
  the fold's BA on its own held-out stacks. The cost of transfer is BA_in − BA.

A data set enters the verdict only if seed 1 labelled at least 200 particles of each label (guessed). The other
data sets are used for training only.

**Two reads, each with its own verdict:**

- **V1, at the seed-1 shifts:** the particles are centred and their classes unknown to the student. Can the network
  carry the class-level judgement down to particles across data sets?
- **V2, at zeroed shifts:** the particles at their picks, on the data sets whose stacks were extracted at the picks.
  Is the student usable before any 2D?

**Verdict, written now.** The margins are guessed. The control is measured in Phase 0 and frozen before Phase 3.

- Positive: the median over H of BA − BA_pipe is at least −0.05, its minimum is at least −0.15, and R ≥ 0.90 on
  every H.
- Negative: the median over H of BA − BA_pipe is below −0.15.
- Otherwise partial: a lead, not a finding.

**Seed spread.** Fold 1 is also learnt at `SIMPLE_SEED` 2 and 3. If fold 1's BA ranges over more than 0.05 across
the three seeds, every verdict is reported as a lead, not a finding, and Hans rules.

**Reported, not gated (Hans's eye): the disagreement averages of Phase 4.** They answer what the verdict cannot:
where the student and the classes disagree, is the student right?

- They are made per H, at V1, for two sets of particles: (a) the particles of selected classes that the student
  rejects, and (b) the particles of rejected classes that it accepts.
- Each set is averaged by its seed-1 class, for classes with at least 10 such particles: `make_cavgs` on a project
  copy that holds only those particles.
- Each average is scored by the class-average model: `model_cavgs_rejection quality_mode=apply quality_model=pool
  prune=no` on the copy.
- The stacks and tables go under `scratch/keep/disagreement/<H>/`.

## 8. Common rules and stop conditions

- **The repository's rules** (`CLAUDE.md`, as in the 2D plan):
  - no new parameter;
  - every number is a named constant in `simple_defs_hhcnn` with its reason;
  - `use ..., only:` imports; large types allocatable;
  - `THROW_HARD` never continued with `//&`;
  - new OpenMP code tested at 1, 11 and 24 threads;
  - the arithmetic of training and prediction is a property of the build and the shape, never of the team;
  - `check_test_registry.py` and `check_descr.py` after any change to the registry, a router, the UI or a header;
  - `git diff --summary` checked for mode changes;
  - nothing committed.
- **Frozen numbers.**
  - No margin or floor of section 7 moves after Phase 0.
  - No constant of section 4 moves to improve a held-out read. The constants are set by their rules on training
    rows only, and the held-out data sets never choose anything.
- **Fix as you go,** with a file-table row added before the edit.
- **Stop conditions:**
  - Phase 0: fewer than `SELECTION_MIN_DATASETS` data sets meet decision 11. Stop and say so.
  - Phase 0: BA_pipe is below 0.70 on more than a third of the data sets (guessed). The machinery's verdict is then
    not stable enough to be a target. Stop and report.
  - Phase 0: no data set was extracted at its picks, so `SELECTION_SHIFT_AUG_FRAC` has no measurement. Stop for a
    ruling.
  - Phase 1: an output of the class or pose student differs at one thread and cannot be made identical. Restore it
    and stop.
  - Phase 2: a data set on which the frame cannot be built leaves the list, recorded. Stop if fewer than
    `SELECTION_MIN_DATASETS` remain.
  - Phase 3: the learn on all data sets reaches a BA_in below 0.70 (guessed). The student then does not learn the
    criterion even within its training sets. Stop before the held-out reads and report.
  - Phase 3: a fold takes more than 2 hours. Cap every fold at 5,000 rows per data set (declared now), record it and
    continue.
  - Any phase: a file outside the table, a decision the code contradicts, or a floor that would have to move.
- **Gates no ruling may waive:** Phase 1's byte identity; the fast gate; `unit_ui`; the registry check; section 7's
  verdict as written.

## 9. Phases

**Phase 0, baseline, labels and control** (the Dell, the base, no source edit).
- Release and Debug builds and the fast gate; the machine recorded; a base build kept.
- Identity references for Phase 1, at one thread and `SIMPLE_SEED` 1, under `scratch/keep/identity_base/`:
  - `learn_class_student` and `predict_classes` on the 2D plan's fixed bgal teacher
    (`/data2/agent_runs/class_student/teachers/`);
  - `learn_pose_student` and `predict_poses` on `T200`.

  Each runs twice, and the outputs and checksums are kept. The files that differ between the two runs are the
  declared exclusions.
- The data sets of decision 11 inventoried: specimen, particles, box, smpd, mskdiam, box/mskdiam, CTF, and whether
  the stacks were extracted at the picks.
- The pipeline of decision 9 at seeds 1 and 2 on each data set. Each data set gets a directory
  `scratch/keep/datasets/<id>/` with `labels_seed1.simple` and `labels_seed2.simple`.
- Per data set: the labelled particles, the share rejected, R_pipe, Q_pipe and BA_pipe. From the data sets extracted
  at their picks: `SELECTION_SHIFT_AUG_FRAC`.

Exit: every data set has its two labelled projects; the control table is complete; `SELECTION_SHIFT_AUG_FRAC` and
the control are frozen in Progress; the identity references are kept. Or a stop of section 8.

**Phase 1, the marginal loss** (identical for the existing students).
- Decision 3 in `simple_hhcnn`: the component, `backward_sample`, `sample_loss`, the mode count and `eval_batch`.
- Decision 3 in `simple_polar_student`: `new`'s optional argument, passed to the network.
- Tests in the HHCNN suite:
  - the gradient check with a support target, over every layer and the readout, in the manner of
    `test_gradient_check`;
  - a one-class support gives the cross-entropy against that class and its gradient, within one unit in the last
    place;
  - a toy problem with two groups of classes, on which the marginal training reaches a group accuracy above a floor
    set from the finished code's profile;
  - thread identity with `marginal` on, at 1 and 11 threads.
- The core's paragraph in `doc/algorithms/pose_student.md`.

Exit: `unit_numerics` passes (HHCNN, polar student, class student, pose student); Phase 0's identity cases reproduce
byte for byte, apart from the declared exclusions; the fast gate passes.

**Phase 2, the student and its programs.**
- `simple_selection_student` (domain):
  - the class, composing a `polar_student` with `marginal` on: training on labels, `psel`, the two thresholds, the
    outcome;
  - the model file: the core's header, plus the bank's groups, the thresholds, the frame, and the training sets'
    names and project hashes;
  - pure rules: the whitening of decision 6, the thresholds of decision 8, section 7's numbers (R, Q, BA, AUC) and
    the bank's starting templates (decision 4).

  Its tester:
  - each rule checked in closed form on small inputs;
  - the model's round trip;
  - the bank's rotation convention: a tensor turned by a known angle, read out as the class student reads e3, comes
    back to its template.
- `simple_strategy_selection_student`:
  - the frame (decision 5) and the whitening pass;
  - a project's labelled rows (`run_selection_rows`, `selection_rows_exec`);
  - `learn_selection_student`: the rows of every project in the `filetab`, the bank, the training, the thresholds,
    the model and its learn report (per training set: rows, held-out rows, label shares, BA_in, AUC; the bank per
    group; epochs and time);
  - `predict_selection`: shared memory, its parts and its master.
- `simple_commanders_selection_student`: the commanders of learn and predict, with their UI entries and routers.
- The workflow test `selection_student`, with its registry entry, test UI entry, router and test policy entry.
  - Fixture: three simulated sets, each with 30% of its particles replaced by noise of the same power, labelled by
    construction (class 1 selected, class 2 rejected).
  - Run: learn on two sets, predict on the third, in shared memory and in two parts.
  - Gates:
    - the model reads back with its bank and thresholds;
    - every labelled particle of the third set is predicted;
    - the two-part table equals the shared-memory one byte for byte;
    - the AUC of psel against the planted labels is above a floor set from the finished code's profile (test policy,
      section 6.2);
    - the AUC gate fails when the labels are permuted (the plant is removed afterwards).
- `doc/algorithms/selection_student.md`.

Exit:
- Release and Debug builds pass, and so does the fast gate.
- `unit_numerics` (including its new sub-suite `selection student`) and `unit_ui` pass.
- `check_test_registry.py` and `check_descr.py` pass.
- The workflow test passes twice.
- On every Phase 0 data set, the frame is built and the whitening and the rows are written.
- One learn on all data sets runs, and its time is recorded: the estimate of a fold's cost, to be tested against
  section 8's two-hour rule.

**Phase 3, the leave-one-out read.**
- `validate_selection_student`: its commander, its UI entry in `simple_ui_cavgproc` beside `model_cavgs_rejection`,
  and its router case in `simple_exec_cavgproc`. It does the following:
  - takes a `filetab` of data-set directories and reads each set's rows once;
  - for each H: learns the fold on the other sets, then runs `predict_selection` on H at the seed-1 shifts and on a
    copy with the shifts zeroed;
  - computes section 7's numbers and both verdicts;
  - writes `validate_selection_student_report.txt`, in plain words for a reader who has not met the student, with a
    details file beside it.
- The run, in order: the learn on all data sets (BA_in, for the stop); then the folds; then fold 1 at seeds 2 and 3.

Exit: the report is complete for every data set; V1 and V2 are written by section 7's rule and checked by hand from
the table in Progress; every model, table and log is under `scratch/keep/loo/`.

**Phase 4, the disagreement averages and the report.**
- Section 7's reported-not-gated averages, for every H that enters the verdict.
- `selection_student_report.md`, in the order of `doc/how2s/how2write_agent_plans.md`, section 8.

Exit: the stacks and tables exist for every H with a verdict; the report is written; the algorithm page states the
outcome.

### File table

Testers of a listed Fortran file are covered by its row.

| File | What changes | Phases |
| --- | --- | --- |
| `doc/implementation_notes/planned/selection_student.md` | this note: Progress, file-table rows | each |
| `doc/implementation_notes/planned/selection_student_report.md` | the report | 4 |
| `doc/implementation_notes/README.md` | this plan's row | each |
| `doc/algorithms/selection_student.md` | the algorithm page (new) | 2 3 4 |
| `doc/algorithms/pose_student.md` | the core's marginal loss, one paragraph | 1 |
| `doc/policies/test_environment_policy.md` | the workflow test and the sub-suite | 2 |
| `src/defs/simple_defs_hhcnn.f90` | the `SELECTION_` constants and file names | 2 3 |
| `src/main/learn/simple_hhcnn.f90` | decision 3 | 1 |
| `src/main/learn/simple_polar_student.f90` | decision 3's optional argument | 1 |
| `src/main/learn/simple_selection_student.f90` | the domain (new) | 2 |
| `src/main/strategies/parallelization/simple_strategy_selection_student.f90` | the strategy routines (new) | 2 3 |
| `src/main/strategies/parallelization/simple_strategy_student_common.f90` | if needed: a shared helper generalised, never a behaviour change | 2 |
| `src/main/commanders/simple/simple_commanders_selection_student.f90` | the commanders (new) | 2 3 |
| `src/main/ui/simple_private_prgs.f90` | `learn_selection_student`, `predict_selection` | 2 |
| `src/main/apis/simple_private_exec_api.f90` | their commanders' imports | 2 |
| `production/simple_private_exec_driver.f90` | their router cases | 2 |
| `src/main/ui/simple/simple_ui_cavgproc.f90` | `validate_selection_student` | 3 |
| `src/main/exec/simple_exec_cavgproc.f90` | its router case | 3 |
| `src/main/commanders/test/simple_commanders_test_class.f90` | the `selection student` sub-suite of `unit_numerics` | 2 |
| `src/main/commanders/test/simple_commanders_test_highlevel.f90` | the workflow test's commander type | 2 |
| `src/main/commanders/test/simple_commanders_test_highlevel_selection_student.f90` | the workflow test (new submodule) | 2 |
| `src/main/exec/simple_test_exec_highlevel.f90` | its router case | 2 |
| `src/main/commanders/test/simple_test_registry.f90` | its entry | 2 |
| `src/main/ui/simple_test/simple_test_ui_highlevel.f90` | its description | 2 |
| `production/CMakeLists.txt` | only if Hans rules open item D | 2 |

## Progress

(one entry per phase, written by the run)
