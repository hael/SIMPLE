# particle_sieving_2: a seed solve2D, the class-average CNN's labels and the class student for the rest

Status, 2026-11-24: delivered, uncommitted, for review (decisions 25 to 29 added 2026-11-24; the key `retrain` of
decision 27 approved by Joe the same day); report
[particle_sieving_2_report.md](particle_sieving_2_report.md). Base: `638d99ae3` (the working tree also holds
untracked scripts under `scripts/`, which this plan does not touch). Supersedes: none. Owner: Joe. Depends on:
- the class student (hhcnn) and its programs `learn_class_student` and `predict_classes`
  (`doc/algorithms/class_student.md`; validated: yes, per data set, on bgal, trpv1 and streptavidin);
- the class-average CNN (sscnn) in Fortran, `simple_cavg_cnn` and `simple_cavg_cnn_predict`
  (`completed/cavg_cnn_fortran_inference_report.md`; validated: yes as a copy of the Python network, to 2e-7;
  **as a rejection criterion: no**, the real-data read of its predictions is still open, report section 6);
- `solve2D` and `refine2D` with `refine=greedy` and `update_frac` (validated: yes, in production).

## 1. Why

`particle_sieving` runs a `solve2D` job on every chunk of the data set, in two tiers. The class student's run
measured `solve2D` at about 13 ms per particle on trpv1, and the student's prediction with its alignment at 0.54
to 1.16 ms per particle (`doc/algorithms/class_student.md`, "What the reads found"). `particle_sieving_2` pays for
the search once, on the particles of at most 100 micrographs or 10,000 particles (decision 25), and sorts the rest
of the data set into those classes with the student, with a greedy refine2D iteration every 50,000 particles to keep
the classes current. Classes are
rejected by the class-average CNN instead of the `sieve` preset of `model_cavgs_rejection`. Whether this selects
as well as the existing sieve is not judged in this plan (decision 14); the plan delivers the program and shows
that it runs and is reproducible.

## 2. Contract

- **`particle_sieving` does not change.** `ptcl_sieve`, its policy, the stream's p03 and p05 and their defaults
  are out of bounds. No line of the new module runs unless the program `particle_sieving_2` is called.
- **One project owner.** The `ptcl_sieve2` object owns the work project. The child programs (`solve2D`,
  `learn_class_student`, `predict_classes`, `refine2D`) run on copies in their own folders; the owner reads their
  records back and writes the user's project once, at the end, as `particle_sieving` does (a temporary file,
  renamed).
- **The children are whole programs run as stages**, through their commanders, as `validate_class_student` runs
  them. The module calls no part of a stage and copies no frame of one.
- **Nothing is a default.** The program is registered at developer visibility until a real-data comparison with
  `particle_sieving` has been read by the maintainer (a later plan).
- **Reproducible:** at one thread and a fixed seed two runs give the same accepted set, class labels and class
  averages (section 6 says where this is gated).

## 3. Decisions and rulings

All dated 2026-11-23. 1 to 14 are Joe's; the request gave 1 to 6, the questions put to him gave 7 to 14.

1. **The seed.** `solve2D` on the selected particles of the first `nmics` (default 200) selected micrographs of
   the project, in project order, into `ncls` (default 100) classes. (Its size revised by decision 25.)
2. **The mask** is the box's default everywhere: `mskdiam = (box - COSMSKHALFWIDTH) * smpd` (the stream's
   `mskdiam_box`), passed to every child. (Replaced by decisions 23 and 28: `mskdiam` given, and required.)
3. **Bad classes.** After the seed, `cavg_cnn` scores every class average; a class is bad when `p_junk > 0.9`
   (strictly). A class with no particles is neither scored nor bad.
4. **The rest of the data set comes in batches of 10,000 particles**, cut at stack boundaries in project order:
   stacks are added until the batch holds at least 10,000 selected particles; the tail is the last batch. Each
   batch is classified by `predict_classes` (the student's prediction and its alignment against the candidate
   classes).
5. **Every 5 batches a greedy round:** one `refine2D` iteration with `refine=greedy` on the particles added by
   those batches and the same number of earlier particles (decision 10).
6. **A new commander beside the existing one, the logic in one new module** (decision 11).
7. **No search after a batch.** The "local refinement" of a batch is `predict_classes`' own alignment after
   prediction. No `refine2D` iteration follows a batch.
8. **Bad classes stay classes and are relabelled after every greedy round** (revised the same day, Joe: rejected
   particles are not left out of refinement). Every ingested particle stays active in the work project and takes
   part in the greedy rounds whatever its class's label, so a bad class keeps its particles and its average, and
   a particle can move between bad and good classes. After each round `cavg_cnn` scores the averages again; a
   class may turn bad or good. Rejection is a label during the run: only the final project sets a particle's
   state, 0 when its class's last label is bad.
9. **The class averages are rebuilt only by the greedy rounds.** Between rounds the batches are predicted against
   the averages and the model of the last round.
10. **The earlier particles of a greedy round are those of lowest `updatecnt`**, by the existing
    `sample4update_cnt` through `update_frac = min(1, 2 * n_new / n_active)`, where `n_new` is the number of
    particles added since the last round and `n_active` all ingested particles. The new particles have
    `updatecnt` 0 and fill the first tier; the lowest remaining tiers give the other half. No new bookkeeping
    field. With no search between rounds, `updatecnt` counts greedy rounds only.
11. **Code home.** All workflow logic in `src/main/sieve/simple_ptcl_sieve2.f90` (type `ptcl_sieve2`, `new`,
    `run`, `kill`) with its tester. A thin `commander_sieve_ptcls2` in `simple_commanders_sieve.f90`. Elsewhere
    only registration (exec router, UI), the constants, and one extraction in `simple_cavg_cnn_predict.f90`
    (Phase 1).
12. **The student is trained after the seed and again after every greedy round**, on a random sample of 10,000
    ingested particles against the current class averages. (Replaced by decision 26: once, after the seed; again
    after the rounds `retrain` lists, decision 27.)
13. **Command-line keys:** `nmics` (200; 100 since decision 25) and `ncls` (100) from the existing keys; the class-average CNN's weights
    file through the existing `infile`, required. The batch size (10,000), the round interval (5), the junk
    threshold (0.9) and the training sample (10,000) are named constants in `src/defs/simple_defs_sieve2.f90`. No
    new key (until decision 27's `retrain`).
14. **Acceptance is mechanical:** a miniature workflow test, a run to completion on a real set, and identical
    results twice at one thread. No comparison with `particle_sieving` and no judgement of the selection in this
    plan.
15. (2026-11-23, Joe, on stop condition 1) **The greedy rounds run at the full box** (`autoscale=no`). `solve2D`
    leaves full-box class averages and FRCs (its final `make_cavgs` overwrites the last stage's cropped ones) beside
    carried class sums at its crop (`cavg_state.bin`; MotAB: 320 px averages, 160 px sums), and `refine2D` neither
    crops references nor blends sums of another box. Round 1 therefore runs as a full update over every ingested
    particle and publishes full-box carried sums; decision 10 applies from round 2. Replaces default f's "startit
    above 1" with decision 18.
16. (2026-11-23, Joe, on stop condition 3) **The workflow test's junk threshold is frozen after the Phase 0 dry
    run**, between the miniature's `p_junk` values, so that both labels occur. The fixture network has random
    weights (`p_junk` 0.48 to 0.52 on MotAB's 100 averages); the CNN is not judged, the mechanics are. If a later
    change moves the split, the test reports which label is missing.
17. (2026-11-23, Joe) **Phase 4 runs on thunder** with the trained `cavg_cnn_v3.bin`
    (`~/Data/ss_sim/fortran_ref/` here). The project is to be found there; Phases 0 to 3 and 5 do not wait for it.
18. (2026-11-23, under delegation) Round r runs `refine2D` iteration `SIEVE2_ROUND_ITER0 + r - 1`, with
    `SIEVE2_ROUND_ITER0 = LPLIM3ITERBOUND + 1` (8): past `refine2D`'s onsets of the in-plane shift search (iteration
    3) and of the band from the class FRCs (`LPLIM1ITERBOUND`, 5), as `solve2D`'s terminal pass is.
19. (2026-11-23, under delegation; confirmed by Joe 2026-11-24, decision 29) **The workflow test is a CTest entry of its own**, which
    raises `SIMPLE_CTEST_BUDGET` from 45 to 46 (an owner decision by the budget's rule): Phase 3 registers it in the
    high-level table, `test=identity` finds only registered entries, and no existing entry runs this program. The
    alternative is a manual test program outside CTest, compared by hand.
20. (2026-11-23, under delegation) **The result carries no sigma2 store**, as `particle_sieving`'s combined result
    has none (`merge_chunk_projfiles_without_sigma2`). A round rewrites a particle's sigma2 record over its band
    only; the batches' particles were inactive when the seed's `calc_pspec` filled every shell, so their records
    stay invalid and every round's commit skips them (the run's spectra are the seed particles'). Handed on, the
    store would stop a later commit on a subset without seed particles (an empty half, seen in the mutation test).
    Seeding the new rows from the global curves (`sigma2_state_extend`, the add-on's route) is left to Joe.
21. (2026-11-23, Joe) **The batch particles' sigma2 records are seeded with `sigma2_state_extend`.** Replaces the
    reason of decision 20 (the result still carries no store, as `particle_sieving`'s has none). After each batch's
    prediction, in the batch's folder, the batch's particles join a pending set; when it holds both halves, their
    records start from the store's global curve of their half, every other record is kept, and the candidate is
    published over the folder's copy of the store, which the next steps take. The extension takes each half's mean
    over the seeded rows, so a pending set of one half (halves follow micrographs: a one-micrograph tail) waits for
    the next batch, and one left at the end is never seeded (logged). The rounds then measure every ingested
    particle's noise over their band.
22. (2026-11-23, under delegation) **The program makes the halves of every selected particle before the seed**
    (`ensure_field_halves` on the work project, the owner's one maker). Found by decision 21's first real-data run:
    the seed's `solve2D` labelled the halves with only the seed's particles active, and the transaction gives every
    unmarked stack a half by micrograph group, each to the half with fewer active particles; the stacks without
    active particles add none, so all of them went to the lighter half and were marked. Every batch particle was
    even in all three earlier thunder runs (so every round summed them into the even class sums) and no batch could
    be seeded. The miniature hid it: the fixture is labelled by index and has no micrograph rows. Its fixture now has
    a micrograph per stack and no halves, and two gates (every batch holds both halves; the halves balanced within
    the partition's 5 % warning) fail without the fix.
23. (2026-11-23, Joe) **`mskdiam` on the command line**, an existing key: given, every step takes it; without it,
    the box's default (decision 2; required since decision 28). On the betagal stream set, `solve2D` collapses many classes into rotationally
    symmetric ring averages at stage 2 (its first centring stage); under the box's mask (283 A, or the automatic
    mode, whose stage 2 runs under the box's mask) they persist to the end, about two thirds of 100 classes, while
    under 174.72 A they recover (three builds, one run each). Run 5's seed had the ring classes.
24. (2026-11-23, Joe) **`mskdiam=0` is automatic:** the seed's `solve2D` runs in its automatic mode (it estimates
    the mask from its class averages) and every later step takes the mask `solve2D` registers with its final class
    averages (`get_mskdiam('cavg')`). The commander drops the zero from the command line before `params%new`, which
    would replace it with a default; without `mskdiam` the box's default stays (decision 2; since decision 28 the
    program stops without it).
25. (2026-11-24, Joe) **The seed ends at `nmics` micrographs (default now 100) or 10,000 selected particles,
    whichever comes first** (`SIEVE2_SEED_NPTCLS`, a constant, decision 13; the stack that reaches 10,000 is the
    seed's last). Replaces decision 1's 200 micrographs. Run 9's seed (20 micrographs, 4,177 particles, 42 per class)
    was thin: its first round relabelled 29 classes and moved 24 % of the seed's particles from rejected to
    accepted. On the betagal stream set the rule takes 45 micrographs and 10,083 particles (runs 10 and 11).
26. (2026-11-24, Joe) **The student is trained once, after the seed, on the seed's particles.** Replaces decision 12's
    retraining after every greedy round. Measured on betagal at `mskdiam=180`, batches of 1,000 (a variant build)
    and a round every 5 batches: training took 3,424 of run 9's 3,918 s (nine trainings, 98 to 603 s each). The share
    of new particles that a round leaves in the student's class was 45 to 65 % trained once (run 10) and 45 to 63 %
    retrained (run 11), on the same batches, 0.3 points apart on average over rounds 2 to 8; replaying run 9's
    batches of rounds 2, 5 and 9 with its first student lost 1.7, 2.6 and 5.5 points against the retrained ones.
    The next round searches every new particle against all classes, so the student only places a particle until
    then. Run 10 took 1,109 s, run 11 4,008 s. Their selections differ (42,918 and 38,266 accepted; the accepted
    sets' Jaccard 0.73, against 0.71 and 0.75 with run 9) as independent runs do: without `SIMPLE_SEED` the seed's
    `solve2D` draws from `/dev/urandom`, so the two seeds differed (75 of 10,083 particles in the same class).
    Joe's first proposal, retraining when the share of particles changing class falls below 5, 2 and 1 %, was not
    taken: over the particles a round searches, that share levels off at 8 to 11 % (run 9), and over all
    particles it falls with the share a round searches (56 % to 10 %), not with the classes' convergence.
27. (2026-11-24, Joe, who approved the new key the same day) **`retrain` retrains the
    student after chosen greedy rounds:** a comma-separated list of rounds (`retrain=2,5`), none by default (decision
    26's one training). After each listed round, once its labels are written and while batches remain, the student
    is trained on a random sample of at most 10,000 ingested particles against the round's class averages
    (`learn_<rr>/`); a listed round no batch follows is ignored with a warning, and an entry below 1 (or not an
    integer, which reads as 0) stops the commander before `params%new`. A new key against decision 13: a
    `character(len=STDLEN)` field like `clustinds`, parsed with `list_of_ints2arr`.
28. (2026-11-24, Joe) **`mskdiam` is required.** Replaces decision 2's box default on the command line: under the
    box's mask (283 A on betagal) the seed's `solve2D` left about two thirds of its classes as rings (decision 23).
    `mskdiam=0` stays the automatic mode (decision 24). The module has no box default any more (a path only the
    miniature reached): its settings take a diameter or 0, and the miniature passes the box's mask itself.
29. (2026-11-24, Joe) **Kept as they are:** batches of 10,000 and a round every 5 batches (on betagal, 5 batches and
    one round, after the last; with the defaults `retrain` serves larger sets); `p_junk` > 0.9 (decision 3; the
    CNN's most probable group was cleaner on four runs, 62.0 to 68.8 % of the accepted particles in protein-view
    classes against 47.6 to 53.5 %); classes the first round empties stay empty (17 to 23 of 100 on betagal; the
    policy says so); every step keeps its copies of the project and the sigma2 store. **Decision 19 confirmed:** the
    miniature stays a CTest entry of its own (budget 46).

Defaults taken without a question (veto cheaply):

- a. Names: program `particle_sieving_2`, type `ptcl_sieve2`, commander `commander_sieve_ptcls2`, suite
  `particle_sieve2`, workflow test `particle_sieving_2_mini` (not the program's name: in `simple_test_exec`,
  `params%new` would find the test's descriptor and not copy the project).
- b. The training sample of decision 12 is drawn from all ingested particles, those of bad classes included, so
  that the student keeps predicting into bad classes (since decision 26: the seed's particles).
- c. After the seed every seed particle's `updatecnt` is set to 1, so that decision 10 counts greedy rounds and
  not the seed's iterations.
- d. A last greedy round runs at the end of intake when batches have arrived since the last one.
- e. The work project is a copy of the user's project with the same rows. Particles not yet ingested are state 0
  there and ingested ones keep their state; the object keeps in memory which rows are ingested and each class's
  label, and the final write sets the states from them (decision 8). A batch's `predict_classes` runs on a copy with only the batch's rows selected.
- f. The greedy round runs with `center=no`, `startit` above 1 and `maxits=1`, as the stream's pool runs its
  fractional iterations; its other settings are those of the seed's last `solve2D` stage (read in Phase 0).
- g. `nparts`, when given, is passed to `solve2D`, `predict_classes` and `refine2D`; training is one process.
- h. A project whose selected particles all fit in the seed runs the seed and the labelling only.
- i. The settings are a record with the constants as defaults, so that the tester and the workflow test can run
  small batches. The command line cannot set them.
- j. No restart. A rerun needs a fresh folder (`mkdir=yes`).

## 4. Not in this plan

- `particle_sieving`, `ptcl_sieve`, the stream (p03, p05, p06) and NICE.
- The student (`src/main/learn/`, `simple_strategy_class_student.f90`), `solve2D`, `refine2D`, the class
  averager, `model_cavgs_rejection`.
- Training or shipping a class-average CNN model: the user exports their own weights file
  (`doc/how2s/how2train_cavg_cnn.md`).
- Whether the CNN's `p_junk` is a good rejection criterion, and whether the selection matches the existing
  sieve's: a later plan, with the maintainer's data.
- Restart, a GUI panel, a stream stage.
- `doc/implementation_notes/planned/class_student_stream.md` plans the student in the pool; its files do not
  overlap with this plan's except the sieve policy, which it edits only in its conditional Phase 5.

## 5. Change types and gates

| Phase | Kind | Gate |
| --- | --- | --- |
| 0 | measurement | builds, fast gate, the children twice at one thread, the data precondition |
| 1 | identical | `cavg_cnn_predict`'s table byte-identical on the fixture and one real stack; suites `cavg_cnn`, `cavg_cnn_input` |
| 2 | new code, off unless called | suite `particle_sieve2`; `unit_ui`; the fast gate |
| 3 | new code | the workflow test, twice at one thread, identical |
| 4 | measurement | the real set to completion, twice |
| 5 | close | docs, report |

## 6. Acceptance

- **Identity (gated), on the miniature.** The workflow test at `nthr=1` and a fixed seed, twice: the final
  project's `ptcl2D` states and classes, the `cls2D` states and the class-average stack are byte-identical.
  Exclusions are declared in Phase 0 from what the children's two base runs show, and none is added later. The
  student's training is known not to be identical across thread counts at the program level
  (`class_student_stream.md`, section 9), hence one thread. Identity at one thread on a real set is not asked:
  the run would take days.
- **Counts that must hold (gated), miniature and real set:** every selected particle of the input ends accepted or
  rejected, none pending; the accepted particles are exactly those whose class is not bad at the end; every greedy
  round's work project holds every ingested particle active, those of bad classes included; the number
  of greedy rounds is `ceiling(nbatches / 5)`; every batch but the tail holds at least the batch size.
- **Negative controls (gated), in the unit suite:** a threshold above 1 rejects nothing; a threshold below 0
  rejects every populated class and the run ends with the "no class left" stop; a batch-status array with one
  pending row fails the final check.
- **The real set (reported, not gated).** Two runs at the machine's full thread count: wall time per phase,
  particles accepted, the number of bad classes after each round, the overlap of the two runs' accepted sets, the
  final class-average montage and its `cavg_cnn` table. The montage and the selection are the maintainer's to
  judge. Artifacts under `<run folder>/report/`.
- No numeric floor is set on the selection. The only numbers in the gates are the request's (200, 100, 0.9,
  10,000, 5; the seed's limits since decision 25) and the miniature's sizes (guessed, frozen in Phase 0 after the
  dry run; the seed's particle cap of 450 added with decision 25, reached by the same fifth stack).

## 7. Common rules and stop conditions

- No new command-line key (decision 13; decision 27 adds `retrain`). No edit outside the file table without a row added first. No
  per-file compile options; large types (`parameters`, `cmdline`, `sp_project` copies) are allocatable components
  released in `kill` (`doc/policies/compile_time_policy.md`).
- Constants in `simple_defs_sieve2.f90`, not in the commander or the module.
- The commander normalizes the command line before `params%new`; the module receives `params` and its settings.
- Never commit or push. Long runs in the background; the real set on thunder if the development machine's disk
  cannot hold it.

Stop and report when:
- `refine2D` cannot run a fractional greedy iteration from the seed's `solve2D` output without an edit outside
  the file table (the carried class sums, `refine2D_blends_carryover`);
- `predict_classes` or `learn_class_student` cannot be run on a copy with a subset of rows selected, as
  `validate_class_student` does, without an edit to the student;
- the miniature cannot produce both a bad and a good class with the fixture's weights file by any threshold-free
  means (the mechanics would be untested);
- the real set or a weights file is not in hand at Phase 4 (Phases 0 to 3 still run).

Gates no ruling may waive: the fast gate; `unit_project suite=particle_sieve` (the existing sieve, unchanged);
`particle_sieve2`; `cavg_cnn`, `cavg_cnn_input`; `class_student`, `student_strategy`; `unit_ui`;
`scripts/check_cli_chain.py`; the identity of section 6.

## 8. Phases

- **Phase 0, the base.**
  - Build Release (and Debug where the disk allows; the CNN run used a second machine), run the fast gate, record
    machine, compiler and load.
  - On the workflow test's fixture (the 2,000 simulated beta-galactosidase particles of `test=class_student`), at
    `nthr=1` and a fixed seed, twice each: `solve2D`; `learn_class_student`; `predict_classes` on a copy with a
    subset of rows selected; one `refine2D` iteration with `refine=greedy`, `update_frac=0.5`, `startit` above 1
    from the `solve2D` output. Record which outputs are byte-identical; these fix the exclusions.
  - Read the seed's last-stage settings a greedy round must repeat (objective, low-pass, box and crop, `trs`).
  - Check the required keys of every child command line against the UI.
  - Dry-run the miniature by hand: seed size, batch size and interval that give at least 3 batches and 2 greedy
    rounds on 2,000 particles; `cavg_cnn` with `testdata/cavg_cnn/weights.raw` on the seed's averages.
  - Data precondition: ask Joe for the real project and the weights file; record their paths.
  - Exit: the table of identical outputs, the settings, the miniature's sizes and the two paths in
    Progress; no stop condition met.
- **Phase 1, the CNN's probabilities as a routine.** In `simple_cavg_cnn_predict.f90`, move the body that reads
  the averages and returns each class's three probabilities and its empty flag into a public routine;
  `cavg_cnn_predict_exec` calls it and writes its table as before.
  Exit: `cavg_cnn_predict`'s table byte-identical to the base's on the fixture and on one real stack; both suites
  pass.
- **Phase 2, the module and the program.** `ptcl_sieve2` with `new(params, settings)`, `run`, `kill` and private
  steps: `run_seed`, `label_classes`, `apply_labels`, `next_batch`, `predict_batch`, `greedy_round`,
  `train_student`, `write_result`. The thin commander, the exec route, the UI descriptor (developer visibility,
  required `projfile` and `infile`), the constants. A unit suite for the steps that need no child: batch cutting
  at stack boundaries with a tail, the label and state rule, the `update_frac` of decision 10, the final check,
  the negative controls of section 6.
  Exit: `test=unit_project suite=particle_sieve2` passes; `unit_ui` and `check_cli_chain.py` pass; the fast gate
  passes; `rg` finds no use of `ptcl_sieve2` outside its module, tester and commander.
- **Phase 3, the workflow test.** `test=particle_sieving_2` on the fixture with the sizes frozen in Phase 0,
  registered in the high-level test table.
  Exit: it passes; run twice at `nthr=1`, the identity of section 6 holds; the counts hold.
- **Phase 4, the real set.** Two runs to completion at full threads.
  Exit: both end with a normal stop and the counts hold; the report folder holds the artifacts of section 6.
- **Phase 5, close.** A section on `particle_sieving_2` in `ptcl_sieve_policy.md` (what it is, that it shares
  nothing with `ptcl_sieve`, its folder layout); the notes index; `<plan>_report.md` in the order of the how-to's
  section 8, with a proposed addition to the `simple-main-sieve` guide for Joe's approval (not applied).
  Exit: `scripts/check_guide_imports.py` and the fast gate pass; the report exists.

### File table

| File | What | Phases |
| --- | --- | --- |
| `src/main/sieve/simple_ptcl_sieve2.f90` (new) | the workflow | 2-4 |
| `src/main/sieve/simple_ptcl_sieve2_tester.f90` (new) | suite `particle_sieve2` | 2-3 |
| `src/defs/simple_defs_sieve2.f90` (new) | the constants and file names | 2 |
| `src/main/commanders/simple/simple_commanders_sieve.f90` | `commander_sieve_ptcls2` | 2 |
| `src/main/exec/simple_exec_sieve.f90` | the route | 2 |
| `src/main/ui/simple/simple_ui_preproc.f90`, `src/main/ui/simple_ui.f90` | the descriptor | 2 |
| `src/main/cavg_quality/simple_cavg_cnn_predict.f90` | the probabilities as a routine | 1 |
| the unit-suite table and the high-level test registry (`src/main/commanders/test/`), `src/main/ui/simple_test/` | the suite and the workflow test, a new submodule for the test | 2-3 |
| `src/main/sieve/sieve.inf` or the folder's source list | the new files, if needed | 2 |
| `src/main/exec/simple_test_exec_highlevel.f90`, `production/CMakeLists.txt` (entry and `SIMPLE_CTEST_BUDGET`), `doc/policies/test_environment_policy.md` (the budget's number) | the workflow test's route and CTest entry (decision 19); added 2026-11-23 | 3 |
| `doc/policies/sieving_and_rejection/ptcl_sieve_policy.md` | the new program's section | 5 |
| `src/main/params/simple_parameters.f90`, `src/main/params/simple_parameters_parse.f90` | the key `retrain` (decision 27); added 2026-11-24 | - |
| this plan, its report, `doc/implementation_notes/README.md` | Progress, the close | 0-5 |

## Progress

**Phase 0, 2026-11-23**, the development machine (Ubuntu, gfortran 15.2, 16 cores, load 0.1), base `638d99ae3`.
- Release build of the base and its fast gate: 15 of 15 (19.7 s). Debug: see Phase 3.
- The code check changed three premises (stop conditions met, ruled the same day): `solve2D` leaves full-box class
  averages and FRCs beside cropped carried sums (MotAB's finished run: `cavgs_iter030` 320 px, `cavgs_iter029` and
  `cavg_state.bin` 160 px), so a greedy round cannot start fractionally from it (ruling 15); the CNN fixture's
  random weights give `p_junk` 0.478 to 0.523 on MotAB's 100 averages, none above 0.9 (ruling 16); the trained
  `cavg_cnn_v3.bin` puts 86 of those 100 classes above 0.9 (reported, not judged).
- Child command lines against the code: standalone `refine2D` defaults to `autoscale=yes` and `mkdir=yes` (both set
  explicitly), `objfun=euclid`, `ml_reg=yes`, `restore_cavgs=yes`, measured sigma2; `refine=greedy` takes
  `trs=MINSHIFT`; the band comes from the class FRCs from iteration `LPLIM1ITERBOUND` (decision 18).
- The children twice at one thread: measured inside the program's own one-thread runs (Phase 3), not on their own;
  no child changed in this run but `cavg_cnn_predict` (Phase 1, byte-identical). No exclusion was needed.
- Dry run of the miniature (eight threads, no junk threshold): the class_student test's 2 000 particles cut into 20
  stacks of 100; a seed of 5 stacks (500 particles), 10 classes, batches of 300 (5), a round every 2 batches (3
  rounds), a training sample of 400. Frozen: these sizes and the threshold 0.503 (measured: the last labels'
  `p_junk` 0.4907 to 0.5164 over 9 populated classes; 4 junk, 5 not).
- Data: thunder, `/mnt/data/test_22Jun26/.workspace_1/404_simple_stream/reference_based_picking/
  reference_based_picking.simple` (233 micrographs and stacks, 53 050 particles, box 224), the source of the
  `mic_series_rbp` study; weights `cavg_cnn_v3.bin` (md5 `1f375006...`, the same file on both machines).
- What it means, in short: the program can be built from existing pieces, with the greedy rounds at the full box.

**Phase 1, 2026-11-23**, development machine. `cavg_cnn_class_probs` (the averages and halves read, prepared and
scored) moved out of `cavg_cnn_predict_exec`, which keeps the table and the log. Exit: the table byte-identical to
the base's on MotAB's 100 averages with the fixture weights and with `cavg_cnn_v3.bin`, the log's lines too; suites
`cavg_cnn` 41 of 41, `cavg_cnn_input` 12 of 12.

**Phase 2, 2026-11-23**, development machine. Files: `simple_ptcl_sieve2.f90` and its tester (new),
`simple_defs_sieve2.f90` (new), the commander, the route, the UI descriptor, the suite row and `suite=` help.
- Exit: `unit_project suite=particle_sieve2` 22 of 22 (both stops in forked children, each with its passing
  control); `unit_ui` 1 042 of 1 042; `check_cli_chain.py`; the fast gate 15 of 15. `ptcl_sieve2` is used by its
  module, tester and commander, and by the suite table and the Phase 3 test.
- The suite's first run caught a defect: the threshold was single precision, so 0.9 compared as 0.89999998 and a
  `p_junk` of exactly 0.9 was junk (decision 3 says strictly above). The threshold is double precision now.
- Found while building the module: a step's copy must keep the work project's name, since the sigma2 store's
  layout digest includes `projname` (a copy named after its step would have its store refused and rebuilt).

**Phase 3, 2026-11-23**, development machine. The test `particle_sieving_2_mini` (submodule
`simple_commanders_test_highlevel_sieve2.f90`), registry row and CTest entry (decision 19).
- Two deviations from the first draft: the test is not named after the program (in `simple_test_exec`,
  `params%new` found the test's descriptor, which needs no project, and skipped the project copy), and in
  `simple_test_exec` the program runs in the test's folder (no execution directory), so the outputs are declared
  relative to it.
- Exit: passes at eight threads (29 to 130 s, 0.25 GB, 136 MB) and twice at one thread (92 to 128 s); 1 167 of 2 000
  accepted at one thread (5, 4, 4 junk classes after the rounds), 1 315 and 1 337 in two eight-thread runs (the
  children are not thread-identical). `test=identity`: 14 of 14 declared outputs byte-identical across the two
  one-thread runs (the result project, five labels tables, three rounds' averages and halves), twice (before and
  after decision 20).
- Mutation test: two plants at once (round copies without batch 1's particles; every particle accepted) failed the
  round, active-particle and accepted-particle gates and the unit suite; the round gate then got its own flag.
  A first plant (rounds without the seed's particles) stopped `refine2D`'s sigma2 commit on an empty half, which
  showed that every round skips the batch particles' sigma2 records (decision 20).

**Phase 4, 2026-11-23**, thunder (gcc-toolset-15 Release build of this working tree as of decision 20, 72 cores,
load about 70 with other work), the set and weights of Phase 0, `nmics=200 ncls=100 nthr=28`, two runs side by
side in `~/ps2_phase4/run1` and `run2`.
- Exit: both ended with a normal stop (5 355 s and 5 800 s; 3.4 GB peak each; 1.1 GB of run directory each). The
  counts hold in both: 53 050 particles decided, the accepted exactly those whose class the last labels keep (0
  mismatches), the round's copy holding all 53 050 ingested particles active, one batch (the tail, 7 089 particles)
  and one round (`ceiling(1 / 5)`). The artifacts are in `~/ps2_phase4/run<n>/report/` (the result's and the
  copies' particle fields, both labels tables, the run report, the round's class averages, halves and montage
  `cavgs_iter008.jpg`).
- Reported, not gated. Wall time per step (run 1 / run 2): seed `solve2D` 4 595 / 4 775 s, labels 2 s, student
  359 / 674 s, the batch 17 / 14 s (prediction 0.67 / 0.64 ms and alignment 1.23 / 0.94 ms per particle; 22 / 19 %
  of particles seeded, 31 % kept the mode), round 1 (a full update of 53 050 particles at 224 px) 381 / 334 s.
  Junk classes (`cavg_cnn_v3`, `p_junk` > 0.9) after the seed 8 / 14, after the round 17 / 27 of 100. Accepted
  45 580 / 45 411 of 53 050; the accepted sets share 40 088 particles (Jaccard 0.79; 5 492 and 5 323 in one run
  only). The round's commit skipped the 7 089 batch particles' sigma2 records, as decision 20 expects.
- The seed ran at about 2.5 min per iteration against about 33 s for the `mic_series_rbp` study's `solve2D` of the
  same 200 micrographs (alone, 16 threads, mask from its class averages): two runs reading the same stacks from
  `/mnt/data` beside other work, and the box's mask (283 A against the study's estimated 145 A).
- Not exercised on real data: a fractional round (the set holds one batch after a 200-micrograph seed) and `nparts`.
- What it means, in short: the program runs to completion on a real set and keeps its contract; the selection
  differs between runs by about a fifth of the accepted particles, which follows the multi-threaded seed.

**Phase 5, 2026-11-23**, development machine. Section 15 of `ptcl_sieve_policy.md`; the report
`particle_sieving_2_report.md` beside this plan; the notes index; a proposed addition to the `simple-main-sieve` guide
in the report (not applied). Exit: `check_guide_imports.py` and the fast gate (15 of 15) pass; the report exists.

**Decisions 21 and 22, 2026-11-23** (Joe asked for the seeding; its first real run found the halves' defect).
Files: `simple_ptcl_sieve2.f90` (`seed_sigma2`, `both_halves`, the halves before the seed), its tester, the
workflow test (fixture with micrographs and no halves; three gates), section 15 of the policy, this plan and the
report.
- Development machine: `particle_sieve2` 26 of 26; `particle_sieving_2_mini` 12 of 12 gates at eight and at one
  thread (942 of 2 000 accepted at one thread; 5, 5 and 5 junk classes after the rounds at eight, 7, 5, 5 at one);
  `test=identity` 14 of 14 byte-identical; the fast gate 15 of 15. Mutation: the seeding removed fails the sigma2
  gate (the rounds skip 600 / 1200 / 1500 records again); the halves left to the seed fail the halves and sigma2
  gates (the seed's partition, active 200 / 300, sent the 15 other stacks to the even half).
- thunder, `run3_eo_one_half` (`nmics=10 ncls=10`, seeding without decision 22, stopped after round 1): the batch
  particles were all even (10 136 even, 0 odd in batch 1) and every batch waited on one half; round 1 skipped all
  50 601 batch records. The two Phase 4 runs above had the same halves.
- thunder, `run4_mics10` (`nmics=10 ncls=10 nthr=32`, both decisions): normal stop in 18:45 (3.4 GB). The
  program's partition: 233 micrograph groups, 26 569 / 26 481. Six batches (10 004 to 10 251 particles, the tail
  656), each with both halves (4 311 / 5 825 to 398 / 258) and seeded; no commit skipped a record. Round 1 a full
  update of 52 394 particles (92 s), round 2 fractional at 0.025 (656 new and 656 of the lowest update count,
  blended with the carried sums, mass ratio 1.009; 5 s). 32 412 of 53 050 accepted, 0 mismatches with the last
  labels, 5 of 10 classes junk. Step times: seed 223 s, a batch about 19 s, the student 113 s on 1 793 particles
  and 595 s on 10 000.
- thunder, `run5_defaults` (`nthr=24`, both decisions, beside `run4_mics10`), replacing the Phase 4 runs' numbers:
  normal stop in 4 509 s (2.9 GB, 1.3 GB on disk); the same partition; the batch (7 089) with both halves (2 388 /
  4 701) and seeded; no record skipped. 45 359 of 53 050 accepted, 0 mismatches; junk classes 18 after the seed,
  25 after the round. Step times: seed 3 733 s, student 401 s, batch 15 s, round 1 357 s. Artifacts in
  `~/ps2_phase4/run5_defaults/report/` (and `run4_mics10/report/`).

**Decisions 23 to 26, 2026-11-23 and 24** (the mask, then the seed's size and the student's training), thunder
(`nthr=64` unless said; the share of accepted particles in protein-view classes judged by eye from the final
montages, `p_junk` > 0.9).
- Masks, 200-micrograph seed: `run6_mask174` (174.72 A) 36 720 accepted, 50.0 %; `run7_mask174` (174 A) 36 614,
  48.8 % (the CNN's most probable group: 28 826, 62.0 %); `run8_maskauto` (`mskdiam=0`, 140 A) 50 044, 37.9 %;
  run 5 (the box's 283 A) 45 359, 41.5 %. No run rejected a protein-view particle.
- `run9_mask180_mics20_b1000` (180 A, 20 micrographs, a variant build with batches of 1 000, the student retrained
  after every round): 3 918 s, training 3 424 s (nine, 98 to 603 s); 45 062 accepted, 50.0 % (most probable group
  35 551, 63.4 %). Earlier particles changing class per round, of those the round searched: 18.4 % in round 2, 8 to
  11 % from round 4; new particles left in the student's class 42 to 62 %. Replaying the batches of rounds 2, 5 and
  9 with the first student (the run's sigma2 store copied in; replays deterministic) kept 41.6, 47.4 and 41.2 %
  against 43.3, 50.0 and 46.7 % with the retrained ones.
- `run10_seed10k_train_once` and `run11_seed10k_retrain` (180 A, the seed rule of decision 25: 45 micrographs,
  10 083 particles; 39 batches of about 1 000, 8 rounds; variant builds; run sequentially): 1 109 s (training
  504 s) and 4 008 s (training 3 398 s, eight). New particles left in the student's class 45 to 65 % and 45 to 63 %;
  accepted 42 918 (47.6 %; most probable group 32 225, 63.3 %) and 38 266 (53.5 %; 29 745, 68.8 %); protein-view
  particles 20 409 and 20 478, 17 823 shared. Classes emptied by round 1: 17 and 23 (17 in run 9).
- The module (decisions 25 and 26, this tree), development machine: `particle_sieve2` 36 of 36 (`test_plan_seed_cap`
  new); `particle_sieving_2_mini` 12 of 12 gates at eight threads (21 s) and twice at one (81 s; 933 of 2 000
  accepted; 7, 5, 5 junk classes after the rounds), its seed of five stacks ended by the cap of 450 before `nmics`
  10, one training; `test=identity` 14 of 14 byte-identical; `unit_ui` 1 042 of 1 042; `check_cli_chain.py`,
  `check_test_registry.py`, `check_submodule_linkage.py`, `check_guide_imports.py`, `check_descr.py`; the fast gate
  15 of 15. Not rerun on thunder: run 10 took this path through its variant.
- What it means, in short: one training on a 10,000-particle seed gives the same placement at a quarter of the
  time; which rule and mask select best is still open (report, section 6).

**Decisions 27 to 29, 2026-11-24**, development machine. Files: the key `retrain` (`simple_parameters.f90`,
`_parse.f90`), the commander (`retrain` parsed and refused before `params%new`, `mskdiam` required), the module
(`retrain_rounds`, `retrains_after`, the folders `learn_<nn>` again), its tester, the workflow test, the UI
descriptor, section 15 of the policy, this plan and the report.
- Exit: `particle_sieve2` 39 of 39 (`test_retrains_after` new); `particle_sieving_2_mini` 13 of 13 gates (the
  training folders' gate new) at eight threads (27 s) and twice at one (82 s; 935 of 2 000 accepted): the student
  trained after the seed (400 of 500) and after round 1 (400 of 1 100), round 3 of the list ignored with its
  warning; `test=identity` 14 of 14 byte-identical; `unit_ui` 1 042 of 1 042; `unit_core suite=parameters` 66 of
  66; `check_cli_chain.py` (595 accepted keys, 556 bindings, 1 823 UI inputs); the fast gate 15 of 15.
- By hand: without `mskdiam` the parse stops on the missing key; `retrain=a`, `retrain=0,2` and `retrain=,` stop in
  the commander before `params%new`, and no execution directory is made.

**Review against the repository's rules, 2026-11-24**, development machine, the tree now at `a4788cb44` (eight
commits pulled after the base; none touches this plan's files).
- Fixed: the box default of the mask left the module (`mask_diameter` and its test deleted: only the miniature
  reached it, which now passes the box's mask itself); three procedure headers cut to three lines; `int2str`
  imported at module level; the report's path built with `filepath`; run measurements removed from a test comment.
- Exit: `particle_sieve2` 36 of 36; `particle_sieving_2_mini` 13 of 13 at eight threads and twice at one,
  `test=identity` 14 of 14 between them; the fast gate 15 of 15. Against the one-thread run before the pull only
  `shsig` differs, in its last digits (1 039 particles; classes, states and averages identical): commit
  `0a8b0903e` changed the polar shift maps.
- Then, on Joe's word: `simple_ptcl_sieve2` named among the controllers that import commander modules (layering
  policy, sections 3 and 8); a row for `particle_sieving_2_mini` in `test_environment_policy.md`; the paragraph
  on `particle_sieving_2` added to the `simple-main-sieve` guide.

