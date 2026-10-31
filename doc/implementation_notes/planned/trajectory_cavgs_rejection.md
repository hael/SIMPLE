# `trajectory_cavgs_rejection`: the nanoparticle class-average rejection as its own program, scored jointly and trained on the maintainer's selections

An unattended run on the Dell (`trajectory_cavgs_rejection`), 2026-10-31. It takes the rejection of nanoparticle
class averages out of `model_cavgs_rejection`, gives it the alignment of each class average to a simulated
particle as a second measurement beside the lattice features, scores the two together, and sets that score on
the maintainer's manual selections. It follows the `trajectory_solve3D` run of 2026-10-29/30, whose two gates
(the lattice gate at an absolute grade and a model gate cut relative to the survivors' median) kept, by the
maintainer's eye on 71-17/traj01 and traj02, far too many poor class averages: in situ they removed 6 to 25
percent of the classes made mostly of frames he had rejected.

## 1. What the maintainer fixed (2026-10-31)

- The nanoparticle rejection is separated completely from the `model_cavgs_rejection` workflow: its own
  program, its own module, its own policy note and tests. `model_cavgs_rejection` goes back to the cryo-EM
  contexts (`chunk`, `pool`, `sieve`) and nothing of the nano context stays in it.
- The alignment of the class averages to the simulated particle (the start model `simulate_nanoparticle` builds
  for the estimated diameter; `refine3D_nano` on the `cls3D` field, one iteration, a global search) is a
  measurement of quality in the new program, beside the lattice features.
- The decision is one joint score of the lattice features and the model correlation, trained on the
  maintainer's manual selections, with its threshold set by the recall-first rule as before; it replaces the
  two coupled gates and their constants (`LAT_PEAK_MIN`, `LAT_PEAK_MIN_RAW`, `K_MODEL`, `model_gate_states`).
- As before: the evidence is independent of the data being judged (the lattice is intrinsic to the class
  average, the model is simulated; no map made from the frames enters the score); the constants are set once
  on labelled data and never moved to let a trajectory through; bugs are fixed in the phase that finds them.
- Added 2026-10-31, 01:40: the training data are new. The maintainer will select by eye on the class averages
  of the workflow's own 2D (the raw frames, centred, `refine2D_nano`), not on his earlier classes, so this run
  produces those class averages for every trajectory, and at several class populations: the raw averages of
  about twenty frames are likely too noisy to judge or to measure, and populations at least three times higher
  are to be tried. The experiments come first and are done by the morning of 2026-11-01; the training on his
  selections is the next run. The redesign (the program of section 2, with the model correlation in the feature
  bank) follows the experiments in this run.

What the tables of the previous run say, and this plan builds on (its `result/keep/` on the maintainer's Mac,
read in the night of 2026-10-30; the detail is in section 3): on the three labelled trajectories rebuilt on raw
frames, the correlation to the simulated particle separates his rejected classes (median 0.15 to 0.22 on every
set, 90th percentile at most 0.30) from his selected ones (median 0.47 to 0.61), and the level transfers between
start models. The relative cut failed because, where he rejected most of a trajectory, the survivors' median is
the junk level. At the 0.98-recall rule nothing removes much on 71-6/traj01 (a few of his selected classes sit at
the junk level); at 0.95 a joint logistic score on peak and correlation removes 96, 57 and 88 percent of his
rejected classes while keeping 97, 98 and 94 percent of his selected frames, and it is robust where an absolute
correlation floor is not (71-17/traj06, the largest particle, has a lower correlation level throughout). By eye,
the lattice peak alone and the correlation alone are each fooled on some classes (a granular blob with a strong
Fourier peak; a fringed class the sphere does not match), which is why the score is joint.

### Open decisions (the maintainer rules before launch; the plan's proposals are in force until he does)

1. The program's name: `trajectory_cavgs_rejection` (a `single_exec` program of the `trajectory_*` category, after
   the program it replaces).
2. The lattice bank (`simple_cavg_quality_lattice.f90`, with its tester) moves to the SINGLE domain as
   `src/main/nano/single_cavg_lattice.f90`, so that `src/main/cavg_quality/` is cryo-EM only. (Alternative: it
   stays where it is as shared code; then only its nano-context plumbing leaves `cavg_quality`.)
3. The recall target of the threshold rule: 0.95 (the plan's proposal, from the table above); the learn report
   prints 0.98, 0.95 and 0.90 beside each other so that the choice is visible.
4. The features of the joint score: the six lattice features and the model correlation. The class population
   is a validity gate only (empty classes), not a feature: his labelled classes and the workflow's time chunks
   have different population semantics, and a learned population weight would not transfer. The learn report
   gives the single-feature baselines (lattice peak alone, correlation alone, the two together) beside the joint
   score, by the same rule.
5. Training data (settled by the ruling of 01:40): the maintainer's selections on the class averages of the
   workflow's own 2D, made after phase 0 at the population he chooses; the earlier labelled sets may join them in
   the next run if their populations are comparable. The class populations tried: 20, 40, 60 and 100 frames per
   class.
6. The alignment of class averages to a model exists three times (`autorefine3D_nano`'s final pass,
   `validate_cavgs_vs_model` and `cavgsproc_nano` in `single_commanders_experimental.f90`, the model gate of
   `trajectory_solve3D`). This run makes it one routine and uses it in the new program, in `trajectory_solve3D`
   and in `autorefine3D_nano`; the two experimental programs keep their own copies unless the maintainer rules
   that they take the routine too.

Terms. A trajectory is one tracked particle: a stack of graphene-subtracted frames (`subtr_frames.mrc`), 0.358 Å.
`refine2D_nano` groups the frames into class averages that are time chunks of about twenty frames. The lattice
features are Fourier measurements of one class average (its strongest Bragg peak over its ring mean, the spot
count, arc, diameter, signal-to-noise, and the band presence). The start model is the face-centred cubic
platinum particle `simulate_nanoparticle` builds for the estimated diameter. Aligning a class average to the
model means finding the model projection it matches best and its correlation to it in the band 3 to 1 Å.
Recall is the fraction of the maintainer's selected classes a rule keeps; removal is the fraction of his
rejected classes it removes. The labelled trajectories are 71-17/traj04, 71-6/traj01 and 71-6/traj04, whose
class averages he sorted by eye in October 2026 (`cls2D_field.txt`, `ptcl2D_field.txt` beside the frames).

## 2. The program

`single_exec prg=trajectory_cavgs_rejection projfile=<p> vol1=<start model> mskdiam=<A> nthr=24`, on a project
whose class averages are registered (`refine2D_nano` or `make_cavgs`), in its own execution directory. Its modes
are the `quality_mode` key the cryo-EM program uses, with the same meanings:

- `apply` (default): measures every class average, aligns the valid ones to the model, scores them, sets the
  class states (`cls2D`, `cls3D`) and the particle states through `map_cavgs_selection`, writes the project;
  `prune=yes` drops the rejected frames. `thres=<probability>` moves the decision threshold for one run, as the
  lattice gate's `thres` did. The model is the built-in preset unless `infile=<model file>` is given.
- `analyze`: the same measurements and score, with the project's `cls2D` state as the manual reference; writes
  the training table and the products below, leaves the project unchanged.
- `learn`: fits the joint score from the training tables of `filetab`, writes the model file (`fname`) and the
  learn report (section 3). `evaluate` scores a model on training tables; `promote` emits the Fortran of a
  built-in preset from a model file.

The measurements, per class average, in `apply` and `analyze`:

1. **Validity.** A class is rejected outright when it is empty, holds invalid pixels, or has no measurable
   spectrum (the nano validity gates as they are now, flatness on the unfiltered average).
2. **The lattice bank**, the six features of `single_cavg_lattice` (the measurement as it is now; `mskdiam` sets
   the soft mask of the spectrum).
3. **The model alignment.** The valid classes are aligned to `vol1` (`refine3D_nano oritype=cls3D maxits=1 lp=1.0
   trs=10 mskdiam=<mskdiam> volrec=no`, the pass `autorefine3D_nano` runs against its final map); the correlation
   of each to its best reprojection is `corr` of `cls3D`, with the orientation and shift. One routine,
   `align_cavgs_to_model` in `single_commanders_nano3D.f90`, serves the new program, `trajectory_solve3D` (its
   commander is `single_commanders_trajectory_solve3D.f90`, a module of its own since 2026-10-31: a workflow
   commander of that size gets its own module) and `autorefine3D_nano`'s final pass.
4. **The score.** A logistic score over the seven features: the lattice peak enters raw (an absolute log10
   ratio), the correlation raw (its level transfers, section 3), the other five by median and MAD over the valid
   classes, clipped, as now. The preset carries the weights, the intercept, the normalization choice and the
   threshold; a class is accepted at probability at least the threshold.

Products in the execution directory: `cavgs_rejection_classes.tsv` (class, population, the seven raw features,
score, state, manual state when there is one, best orientation and shift); `cavgs_ranked.mrc` with
`cavgs_ranked_reprojs.mrc`, every valid class average and its model reprojection sorted by score, lowest first,
and `cavgs_ranked.txt` naming the class of each slot (the maintainer judges from these; the twenty-pair sample
of the previous run is gone); `cavgs_selected.mrc` and `cavgs_rejected.mrc`; in `analyze`, the training table
`cavgs_rejection_training.txt` (the same `# key=value` header and one row per class as the cryo-EM table, so
that the same reader serves). `apply` writes `diam` (the spot diameter) and `lattice_peak` on every scored
`cls2D` row and the diameter's median over the accepted classes on the class-average entry of `os_out`, as now.

Ownership: the program is `ui -> exec -> commander -> domain`: the UI entry in `single_ui_trajectory.f90` with
every key named here, the router case in `single_exec_trajectory.f90`, the commander
`commander_trajectory_cavgs_rejection` in `single_commanders_trajectory.f90`, the domain in
`src/main/nano/single_cavg_rejection.f90` (feature assembly and normalization, the score and the decision, the
model file and the built-in preset, the training table, the learner on `simple_quality_logistic`, the ranked
products) with `single_cavg_rejection_tester.f90`, and the lattice bank in `single_cavg_lattice.f90` (decision 2).
The program's checks (`mskdiam` required by `apply` and `analyze`, `vol1` by `apply` and `analyze`, `filetab` by
`learn` and `evaluate`, `thres` a probability in (0,1)) are the commander's, on the command line before
`params%new`, beside its defaults: the parameter class holds no program-specific checks (the maintainer's ruling
of 2026-10-31; the block `trajectory_solve3D` had put there was moved into its commander the same night).
Follow `doc/policies/layering_and_project_ownership_policy.md`: the commander reads the project once, the domain
never reads a project, the commander writes the changed segments once.

What leaves `model_cavgs_rejection` and `src/main/cavg_quality/`: the `nano` context (`quality_context=nano`, the
nano validity gates and local-variance measurement, the nano bank size and names, the nano relational band),
the `nano` preset (`quality_model=nano`, `nano_model_spec`, `set_nano_threshold`, `thres` on that program), the
`lattice` feature policy of `learn`, the nano annotations of `apply` (`diam`, `lattice_peak`), the lattice
comments of the training table, and the `cavg quality nano` sub-suite of `unit_numerics`. The cryo-EM program's
`analyze`, `learn` and `evaluate` outputs are byte-identical before and after (section 5).

`trajectory_solve3D` changes with it: its stages 6 and 7 (the lattice gate and the model gate) become one
stage, `6_trajectory_cavgs_rejection`, run with `prune=yes`; the stage count becomes 10; its contained model-gate
code, its `model_gate_nsd` key and its sampled pairs go; `thres` passes through as the score threshold. Its
records and resume are unchanged otherwise.

## 3. The class averages for the maintainer's selections, and the class-population experiment (phase 0)

For every trajectory the run makes the class averages he will select from: the workflow's own 2D on the raw
frames, exactly as `trajectory_solve3D` runs it up to its 2D stage (`trajectory_import` of `subtr_frames.mrc`;
`center2D_nano ncls=10`; `estimate_diam` on the centred chunk averages, `mskdiam = 1.5 x max_diam`, `moldiam =
min_diam`; `simulate_nanoparticle` for the start model; `refine2D_nano mskdiam=<mskdiam> center=no`), with one
difference: the class count. `refine2D_nano` at its default gives time chunks of about twenty frames (649 classes
for the 71-17 trajectories, 684 for 71-6: the ceiling of frames over twenty). The experiment runs it at four
populations, P = 20, 40, 60 and 100 frames per class (`ncls = ceiling(nframes / P)`; the centering, the diameter and the start model are made once per
trajectory and shared), so that the maintainer can judge at which population the raw averages show their lattice
and choose the one to select on, and so that the numbers below say what each population buys.

Per trajectory and population, in `scratch/keep/popstudy/<dataset>/<trajectory>/P<frames>/`: the project after
the 2D and its class averages (`cavgs_iter*.mrc` of the last iteration; kept, these are what he selects on);
the features of every non-empty class: the six lattice features (`model_cavgs_rejection quality_model=nano
mskdiam=<mskdiam> prune=no thres=0.01`, which writes `cavgs_quality_features.txt` and rejects only the invalid
classes) and the correlation, orientation and shift from the alignment to the start model (`refine3D_nano
oritype=cls3D maxits=1 lp=1.0 trs=10 volrec=no`, as the previous run's model gate did; its scripts under
`gate_raw/` and the model-gate code in `single_commanders_nano3D.f90` are the templates), merged into
`classes.tsv` (class, population, the six features, correlation, e1 e2 e3 x y); the class average / reprojection
pairs of every non-empty class as one stack, `pairs_by_corr.mrc`, sorted by correlation lowest first, and the
class averages alone sorted by correlation and by lattice peak (`cavgs_by_corr.mrc`, `cavgs_by_peak.mrc`), each
with a `.txt` naming the class of every slot (pure-Python MRC writing as the previous run's `mrcpairs.py`; the
contact sheets are made on the maintainer's Mac after the fetch). The maintainer's selection on one of these
stacks becomes a `cls2D` state of that project (`selection` with the kept classes), the training set of the next
run.

The numbers, in `scratch/keep/popstudy/summary.tsv`, per trajectory and population: classes with frames, the
population's median and range, the wall time of the 2D, the fraction of consecutive frames sharing a class, the
distribution of the lattice peak (median, 10th and 90th percentiles, the fraction at or above 0.52 and 0.65)
and of the correlation (the same percentiles; the gap between its lower and upper modes, if the histogram shows
two). On the three labelled trajectories, where every frame carries the maintainer's earlier state
(`ptcl2D_field.txt`, through the frame's stack index), each class gets the majority label of its frames and a
purity, and the table adds, per population: the fraction of frames in classes of purity at least 0.8; what a
perfect class gate at that granularity could do (reject every majority-rejected class: the recall of his selected
frames and the removal of his rejected frames, the ceiling any class-level rule has at that population); the AUC
of the lattice peak and of the correlation against the majority label; and, by the recall-first rule at 0.95 on
the three sets together, the removal of the lattice peak alone, the correlation alone and both. The Progress
entry says, in plain words, at which population the features separate his labels best and what it costs in
granularity.

### The training of the score, in the next run

When the maintainer's selections exist (the stacks above, one population per trajectory, `cls2D` states on the
projects), `trajectory_cavgs_rejection quality_mode=analyze` on each gives the training tables, and the three
labelled trajectories' earlier classes, rebuilt on the raw frames as the previous run's raw study did (its
scripts under `~/agent_runs/trajectory_solve3D/scratch/keep/gate_raw/`, read only), may join them as further sets
if their populations are comparable.

`learn` fits the joint score on the tables of `filetab`: the logistic problem of `simple_quality_logistic`
(L-BFGS-B, the dataset roles and class weights the cryo-EM learner uses, its lambda grid), features as in
section 2. The threshold is the recall-first rule of the previous runs: the strictest threshold at which the
recall of the maintainer's selected classes is at least the target on every training set (decision 3), the
frame-weighted recall reported beside. The learn report (`cavgs_rejection_learn_report.txt`) holds, for the
target and for 0.98, 0.95 and 0.90: per set and pooled, recall, removal, frame-weighted recall, the threshold;
the same for the baselines of decision 4 by the same rule; the leave-one-trajectory-out check (the score
refitted without each set, its recall and removal on that set at the threshold the others give); the feature
weights and the AUC of every feature. The score is promoted as the built-in preset (`quality_mode=promote`, the
snippet pasted into `single_cavg_rejection.f90`), and `K_MODEL`, `LAT_PEAK_MIN`, `LAT_PEAK_MIN_RAW` and
`model_gate_states` are deleted with the old gates.

A stop condition for that run: if at the target recall the joint score removes no more of the rejected classes
than the lattice peak alone on every set (within 0.02), the alignment carries no evidence beyond the lattice;
stop and report, the maintainer decides. The tables of the previous run say it will not come to that.

## 4. What is kept, and the comparison

Phase 0 keeps everything of section 3 under `scratch/keep/popstudy/` (the projects and class-average stacks of
every trajectory and population are the deliverable; the 2D's earlier iterations and any copied frame stacks are
dropped), `summary.tsv`, the scripts, and the logs. Phase 1 keeps its test logs, the byte gate's before/after
files and the `fcc` run's text products under `scratch/keep/`.

For the next run, when the ten trajectories go through `trajectory_solve3D` with the trained score: per
trajectory, under `scratch/keep/pipeline/<dataset>/<trajectory>/`, as the previous run kept it, plus the new
program's ranked stacks and class table; the comparison table beside the previous run's row
(`~/agent_runs/trajectory_solve3D/scratch/keep/pipeline/table.tsv`, read only) and the hand chain's: frames in
and kept, classes rejected (validity, the score), the score's threshold, the diameter and mask, the 3D's FSC at
0.143 at the program's mask, the strong atoms, the largest shifts, the wall times; and the in-situ check on the
labelled trajectories (`insitu.py`, read only) repeated for the new score.

## 5. Rules of the run

- The repository's rules as in the previous run: release 4 keeps no backwards compatibility; `THROW_HARD` and
  `THROW_WARN` are macros, never continue their argument with `//&`, build the message in a string variable;
  Fortran does not short-circuit, never test a flag and an element of an allocatable array in one condition; no
  goto; no reliance on integer overflow; dead code met during the work is deleted; `use ..., only:` imports;
  large types never as plain components; new OpenMP code tested at 11 and 24 threads; `git diff --summary`
  checked for mode changes before a phase ends; no new CTest entries, `SIMPLE_CTEST_BUDGET` unchanged;
  `python3 scripts/check_test_registry.py . --verbose` and `python3 scripts/check_descr.py .` whenever a suite
  table, a tester registration, a router, a test-UI help or a file header changes; a new program needs its UI
  entry, router case and commander together, and its checks are the commander's, on the command line before
  `params%new` — nothing program-specific goes into the parameter class.
- Nothing is rejected on evidence the denoiser has touched, and nothing in the score comes from a map made of
  the frames being judged.
- The constants come from the labelled data by the rule of section 3, set once in the next run, and are not
  revisited on the downstream numbers; no threshold or default moves to let a trajectory through. The
  provisional preset of phase 1 is a placeholder that lets the program and its test run; it is not a result and
  the policy note says so.
- Tests: `single_cavg_rejection_tester` (the score in closed form on a small table; the threshold rule on a toy
  training set with known recall and removal; the model file written and read back; the ranked products of a
  synthetic stack) and the lattice tester, moved, as sub-suites `cavg lattice` and `cavg rejection` of
  `unit_single`; `unit_numerics` loses `cavg quality nano`; `unit_ui`, `unit_project`; the `single_workflow`
  `fcc` suite runs `trajectory_solve3D` with the new stage and keeps its assertion that the gate loses at most
  5 percent of the simulated trajectory's frames beside the map's floors and the pose checks; `wurtzite` as it
  is. The cryo-EM byte gate: before any source change, phase 1 runs `model_cavgs_rejection quality_mode=analyze` on a
  cryo-EM project with class averages (a `refine2D` of the bgal set under `/data`, its `cls2D` states as they
  are, `quality_context=chunk` and `pool`) and `learn` on the resulting table, and keeps the outputs under
  `scratch/keep/byte_gate/`; after the separation the same commands give byte-identical files (`cmp`), and
  the fast gate, `unit_numerics` and `unit_single` pass.
- The data under `/data2/NanoX_trajs` and `/data/NanoX/zenodo` are read only; the previous run's directory
  `~/agent_runs/trajectory_solve3D` is read only; python without numpy (`python3.12`).

## 6. Phases

**Phase 0, the class averages and the population experiment** (section 3), on the build as it is, before any
source change. The ten trajectories in the previous run's inventory order, each at the four populations; the
features, the ranked stacks, `summary.tsv`; the Progress entry with the plain learning. Exit: every trajectory
has its four stacks and tables, the summary is complete, and the labelled three carry the purity, ceiling and
AUC numbers. This phase is done first so that the maintainer has the stacks in the morning whatever the next
phase does.

**Phase 1, the separation and the program.** The modules, the program (`apply`, `analyze`, `learn`, `evaluate`,
`promote`) with a provisional preset (the lattice peak and the correlation at the weights the table of section 1
gave, threshold at the 0.95 rule), the shared alignment routine, `trajectory_solve3D` on the new stage with a
key for the class population of its gating 2D (`ncls`, passed through to `refine2D_nano`; the default as now
until the next run sets it from phase 0's result), the nano context removed from `model_cavgs_rejection` and
`cavg_quality`, the tests moved and written, the policy note split (the nano sections of
`model_cavgs_rejection.md` become `trajectory_cavgs_rejection.md`), the skills' two pointers updated. Exit: the
build, the fast gate, `unit_single`, `unit_numerics`, `unit_ui`, `unit_project`, the registry and description
checks pass; the byte gate is identical; `single_workflow` `fcc` and `wurtzite` pass; `trajectory_cavgs_rejection
quality_mode=analyze` on one phase-0 project (its states as they are) writes a training table the learner reads.

**Phase 2, the report.** `trajectory_cavgs_rejection_report.md` beside this note: the population experiment's
numbers and the plain learning, where the stacks for the maintainer's selections are, the program as built and
what left the cryo-EM program, every fix, the decisions taken under delegation, what is open (the training, the
ten trajectories: the next run).

**The next run**, after the maintainer's selections: the training of the score (section 3), its promotion, the
`fcc` test with it, the ten trajectories through `trajectory_solve3D` at the chosen population, the comparison
table and the in-situ check (section 4), and its report.

## 7. File table

A phase edits only its files; a file the work turns out to need is added here, with its phase and a one-line
reason, before it is edited. Testers (`*_tester.f90`) of a listed Fortran file are covered by its row.

| File | What changes | Phases |
| --- | --- | --- |
| `doc/implementation_notes/planned/trajectory_cavgs_rejection.md` | this note: Progress, file-table rows | each |
| `doc/implementation_notes/planned/trajectory_cavgs_rejection_report.md` | the report | 2 |
| `doc/policies/sieving_and_rejection/trajectory_cavgs_rejection.md` | the new policy note (the nano sections of `model_cavgs_rejection.md`, rewritten for the program) | 1 2 |
| `doc/policies/sieving_and_rejection/model_cavgs_rejection.md` | the nano context and the two gates leave it | 1 |
| `doc/policies/test_environment_policy.md` | the `single_workflow` claim and the moved sub-suites | 1 |
| `.github/skills/simple-cluster-cavgs-quality/SKILL.md` | the nano pointers leave; one line names the new program | 1 |
| `.github/skills/simple-main-nano/SKILL.md` | the pointer to the new program and modules | 1 |
| `src/main/nano/single_cavg_lattice.f90` | the lattice bank, moved from `cavg_quality` (decision 2), with its tester | 1 |
| `src/main/nano/nano.inf` | the folder description names the class-average modules | 1 |
| `src/main/nano/single_cavg_rejection.f90` | the domain: features, score, model file, preset, training table, learner, ranked products; tester | 1 |
| `src/main/commanders/single/single_commanders_trajectory.f90` | the commander `trajectory_cavgs_rejection` | 1 |
| `src/main/commanders/single/single_commanders_nano3D.f90` | `align_cavgs_to_model`; `autorefine3D_nano`'s final pass on the routine | 1 |
| `src/main/commanders/single/single_commanders_trajectory_solve3D.f90` | `trajectory_solve3D` (its own module since 2026-10-31) on the new stage, with the `ncls` key | 1 |
| `src/main/ui/single/single_ui_trajectory.f90` | the UI entry; `trajectory_solve3D`'s keys | 1 |
| `src/main/exec/single_exec_trajectory.f90` | the router case | 1 |
| `src/main/params/simple_parameters.f90` | `model_gate_nsd` goes | 1 |
| `src/main/params/simple_parameters_parse.f90` | the same | 1 |
| `src/main/params/simple_parameters_phases.f90` | only removals: the nano checks of `model_cavgs_rejection` go; no program-specific check is added | 1 |
| `src/main/ui/simple/simple_ui_cavgproc.f90` | `model_cavgs_rejection` loses `thres` and the nano choices | 1 |
| `src/main/ui/simple_ui_params_common.f90` | `quality_context` and `quality_model` lose `nano` | 1 |
| `src/main/cavg_quality/` | the nano context, preset, policy and annotations leave every module; the lattice module and tester move out | 1 |
| `src/main/commanders/simple/simple_commanders_cavgs.f90` | `exec_model_cavgs_rejection` loses the nano branches | 1 |
| `src/main/commanders/single/single_commanders_experimental.f90` | the two experimental programs on the shared routine, if decision 6 says so | 1 |
| `src/main/commanders/test/simple_commanders_test_class.f90` | the sub-suites moved and added | 1 |
| `src/main/ui/simple_test/simple_test_ui_class.f90` | the `unit_single` and `unit_numerics` suite help | 1 |
| `src/main/commanders/test/simple_commanders_test_single.f90` | `single_workflow`: the stage directories of the new layout | 1 |
| `src/main/ui/simple_test/simple_test_ui_highlevel.f90` | the `single_workflow` description | 1 |
| `production/CMakeLists.txt` | the `single_workflow_fcc` comment | 1 |
| `src/main/quality/simple_quality_logistic.f90` | the core the learner calls: fixes only (its `learn` output stays byte-identical) | each |
| `src/main/commanders/single/single_commanders_nano2D.f90` | center2D_nano, refine2D_nano: fixes | each |
| `src/main/commanders/simple/simple_commanders_sim.f90` | simulate_nanoparticle: fixes | each |
| `src/main/commanders/simple/simple_commanders_imgproc.f90` | estimate_diam: fixes | each |
| `src/main/strategies/parallelization/simple_refine3D_strategy.f90` | refine3D on `cls3D`: fixes | each |
| `src/main/strategies/parallelization/simple_make_cavgs_strategy.f90` | the rebuilt classes: fixes | each |
| `src/main/project/simple_sp_project_ptcl.f90` | map_cavgs_selection, prune: fixes | each |


## 8. Progress

(one entry per phase, written by the run)
