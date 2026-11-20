# Multi-state workflows under scoring-noise annealing: validation on the SNR 0.05 benchmark

Status, 2026-11-20: planned, for the autonomous run `ms_anneal` on dell2. Base: `f9fd7ca21`
(the commit that adds this plan sits on top; the launch copies the checkout at that commit and
Phase 0 records the base it got). Supersedes: none. Depends on: `5709265cb` (`sigma2_anneal=yes`
the default of solve3D, refine3D_auto, refine3D_local and refine3D_states; validated: by the
`multistate_workflows` entry once, on 3,000 particles at SNR 1; never on the benchmark the
multi-state work is judged on).

## 1. Why

Every multi-state control so far scored the particles against measured noise. The de novo runs
that set the standing success criterion (`multistate_ab_initio_recovery.md`, ruling 12) ran on
`1c12bb4eb` and separated the three states of the SNR 0.05 benchmark F05 in 4 of 5 runs (final
accuracy 0.996, 0.863 partial, 0.969, 0.969, 0.969); the validation of `4282d58a3`
(2026-11-15) separated its one finished run (0.927 at iteration 82, recall 0.89 / 0.97 / 0.79). In
all of them solve3D ran without annealing: its default was `no` until `5709265cb` (2026-11-20), and
the log of the 2026-11-15 run has no `SIGMA2 ANNEAL` line. Since then the multi-state five-stage
schedule scores stage 1 (shc, 30 A) from the white ceiling toward the smoothed bootstrap
(weight -1 to -0.5), stage 2 (shc, 24 A) on to the bootstrap (-0.5 to 0), stage 3 (the first
`prob` stage, where the states are gathered: in run 1 of 2026-11-11 the largest assembly went from
0.69 to 0.95 during it) from the bootstrap to measured noise (0 to 1), and stages 4-5 at measured
noise (`set_refine3D_sigma2_schedule`, `simple_solve3D_controller.f90`, with
`NSTAGES_INDEPENDENT = 5`, `PROB_REFINE_STAGE = 3`, `STICKY_LAST_STAGE = 4`). refine3D_states,
the route Hans and Cyril use on real data (`consensus`: single-state solve3D, then
refine3D_states), now smooths the spectra it starts from, scores its whole frequency march at
that anchor and ramps to measured noise over the march's last three planned iterations
(`run_frequency_march`, `FREQUENCY_BLOCK_NITS = 3`). Whether the states still separate under the
new scoring is unmeasured. This run measures it, on the same commit, with and without the
switch, and only if annealing loses states does it change code, confined to `nstates > 1`.

The 23 commits between `4282d58a3` and the base also touch solve3D's path (the sticky subset of
`8363af240`, the equal class quotas of `11f832a8e`, the removal of the 3D in-plane posterior
insertion in `fd9af0add`, `nsample_ramp` given up in `6668fb842`, the PorF low-pass fix
`a4881e6c4`). `8363af240` says multi-state is unchanged by the sticky subset. The arm without
annealing is therefore also the regression check of the base for multi-state since
`4282d58a3`.

## 2. Contract

- A validation on a frozen export of the base. Phases 0-3 change no production code. Phase 1
  changes test code only: the `multistate_workflows` routes that run solve3D or refine3D_states
  pass `sigma2_anneal=` on when the test is given it, and `test=profile` forwards the key as it
  forwards `suite=`. With the key absent the routes run exactly as before (the entry's command
  lines are unchanged: the off-state invariant of the test change).
- Phase 4 runs only if the verdict of decision 5 or 7 names a loss. A candidate there changes
  the schedule or the default for `nstates > 1` only; `nstates = 1` stays byte-identical at one
  thread (`simulated_workflow_6vxx` identity, base export against the tree); no new parameter;
  no floor, tolerance or timeout of an existing test relaxed; no edit of the anchor smoothing,
  the blend, or the 2D schedule.
- The two arms differ in one key on the command line: arm Y runs the default
  (`sigma2_anneal=yes`), arm N adds `sigma2_anneal=no`. With `no`, solve3D's stages run refine3D
  without the key (refine3D's own default is `no`), and refine3D_states prepares no anchor and
  sets no weights: arm N is the scoring every earlier control had.

## 3. Decisions and rulings

1. (2026-11-20, Hans) Pick the multi-state controls on dell2 up again now that `sigma2_anneal=yes`
   is the default of solve3D and the refine3D controllers, as an autonomous testing run.
2. (2026-11-20, the plan) Design: a same-commit A/B. Arm Y, the default; arm N, `sigma2_anneal=no`
   on the command line; nothing else differs. `SIMPLE_SEED` stays unset for every benchmark run
   (each has its own random start), so the arms are two profiles, not pairs.
3. (standing: `multistate_ab_initio_recovery.md` ruling 12 and section 4) The judge is the SNR 0.05
   benchmark F05 on dell2: `simple_exec prg=solve3D pgrp=c1 mskdiam=300 nstates=3 nparts=2 nthr=8
   print_states=yes projfile=P2_05.simple`, from a fresh copy of the 2D solution `P2_05` per run.
   A run's result is its last `states_iter*.txt` (the final assignment at `prob_beta` 1), scored by
   `acc3.py --truth labels_f05.txt`: *separated* when accuracy is at least 0.90 and every true
   state's recall at least 0.75; *partial* when accuracy is at least 0.75 or exactly two states are
   recovered; *failed* otherwise (origin: the recovery plan, section 4; measured there: the
   separated runs sit at 0.969-0.996, the partial one at 0.863). The stopping rule stays: a run
   whose stage-4 accuracy is below 0.60 is stopped and counts as *failed* (measured: separated runs
   are at 0.965 or more after stage 4, the partial one at 0.765).
4. (2026-11-20, the plan) Three F05 de novo runs per arm, one at a time, in the order Y N Y N Y N,
   so that whatever drifts on the machine touches both arms alike. About 6-7.5 h each (measured,
   `nparts=2 nthr=8`, load about 16).
5. (2026-11-20, the plan; frozen before Phase 2) The reading of Phase 2, written before any run:
   - Y separates at least 2 of 3: the standing criterion holds under the default. No code change.
     The difference between the arms is reported as a lead (three runs per arm cannot measure a
     modest shift: at a per-run separation probability of 0.8, 2 of 3 or better comes up 9 times
     in 10).
   - Y separates at most 1 of 3 and N at least 2 of 3: annealing loses states in multi-state de
     novo. Phase 4 runs candidate A1, then A2, each judged by the same rule.
   - Both arms separate at most 1 of 3: the base itself no longer meets the criterion (at 0.8 per
     run this happens by chance about once in a hundred). STOP with the evidence: locating which of
     the 23 commits since `4282d58a3` did it is not in this plan.
   - Y separates 3 of 3 and N at most 1 of 3: annealing helps. No code change; reported.
6. (2026-11-20, the plan; frozen before Phase 1) The cheap fixtures, as profiles of 20 seeds per arm
   (`test=profile target=multistate_workflows suite=<route> nseeds=20 ncunits=16`, two executions at
   a time at eight threads each): the `shc` route (de novo solve3D on 3,000 particles; last profile
   15 of 20 resolved every state, 2026-11-05, dell2; floor 6 of 20), the `states` route
   (refine3D_states from the true poses and the consensus map; 20 of 20 in its profile) and the
   `consensus` route (single-state solve3D then refine3D_states; 40 of 40 in October, borderline on
   2026-11-04: 1 of 3 eight-thread runs lost a state). Per route the count of runs that resolve
   every true state and the medians of accuracy, adjusted Rand index and the smallest recall are
   recorded for both arms. A difference of 5 or more runs of 20 between the arms is a finding; less
   is a lead (guessed from the binomial spread: at 15 of 20 the standard deviation is 1.9 runs, so 5
   is about 2.6 of them). The `shc` floor of 6 of 20 is a gate for arm Y (it is the production
   default); below it, Phase 4 candidate A2 is tried whatever Phase 2 says. The `labels` route runs
   plain refine3D, whose default is `no`: it is not profiled.
7. (2026-11-20, the plan; frozen before Phase 3) refine3D_states on F05, the real-data route on the
   benchmark, has never been measured. Phase 3 runs it from the true poses and the consensus map
   (as the `states` route does: `reconstruct3D` of the true-pose project `PT05` with `nstates=1`,
   then `simple_exec prg=refine3D_states nstates=3 pgrp=c1 mskdiam=300 nthr=16 [sigma2_anneal=no]`),
   two runs per arm in the order Y N Y N, judged by decision 3's verdict on the project's `ptcl3D`
   state column (`print_project_field`), with a wall-clock budget of 10 h per run (guessed; a run
   over it is stopped and recorded as *unfinished* with its last labels). Two per arm is a lead,
   not a finding: a loss here (Y 0 of 2 resolved, N 2 of 2) names candidate B1 for Phase 4; any
   other outcome is reported.
8. (2026-11-20, the plan) Candidates of Phase 4, in this order, each confined to `nstates > 1`,
   each judged on F05 by decision 3's verdict over 3 runs against arm N's count, with the gates of
   section 5, and tried only when a decision above names it:
   - A1: the five-stage multi-state schedule reaches measured noise at the first `prob` stage: the
     white-to-bootstrap ramp spans the shc stages 1-2, stage 3 onward scores at weight 1 (in
     `set_refine3D_sigma2_schedule`: for `nstates > 1`, `last_ramp = PROB_REFINE_STAGE - 1`).
   - A2: solve3D's commander sets `sigma2_anneal=no` for `nstates > 1` unless the user gave the key
     (`exec_solve3D`, where it sets `yes` today).
   - B1: refine3D_states scores its march at measured noise (its commander's default `no`), the
     user's explicit `yes` honoured.
   At most two candidates run (about 21 h each). A kept candidate needs: separation at least 2 of 3
   and not below arm N's count; `nstates = 1` identity; the fast gate; the entry; the `shc` profile
   at or above 6 of 20 if its code path changed.
9. (2026-11-20, the plan) Sync rule: the run works on the export it started from and never merges
   master. The report names the commit it judged; Hans applies a kept candidate to the current
   tree.
10. (2026-11-20, the plan) Threads and jobs: one F05 de novo at a time (`nparts=2 nthr=8`); one
    refine3D_states at a time (`nthr=16`); profiles at `ncunits=16`; the nightly at `ncunits=16`;
    builds beside a job only when nothing is being timed for its wall time (nothing in this plan
    is judged by time, but every run's wall time and load are recorded for the next plan).

Defaults taken without a question: run name `ms_anneal`, root `dell2:~/agent_runs/ms_anneal`;
F05 inputs read-only from `~/agent_runs/ms_abinitio/scratch/keep/F05` on dell2 (`P2_05/`, `PT05/`,
`labels_f05.txt`, `truth_f05.txt`, `truthcheck/` maps); the tools and scripts of
`~/agent_runs/ms_validate_4282d58/keep/{tools,scripts}` copied and repointed (acc3.py,
run_metrics3.py, largest_traj.py, dock_eval.sh, dn_run_f05.sh, dn_f05_one_p3.sh, dn_analyze.sh,
stats_watch.sh, stage_proj_watch.sh, stop_rule_watch.sh); run names `dn_F05_Y_{1,2,3}`,
`dn_F05_N_{1,2,3}`, `st_F05_Y_{1,2}`, `st_F05_N_{1,2}`.

## 4. Not in this plan

- Locating a base regression among the 23 commits since `4282d58a3` (a STOP, decision 5).
- The 2D schedule (`solve2D`, default `yes` since 2026-11-14, judged in
  `sigma2_anchor_anneal_2D.md`), refine3D_auto and refine3D_local (single-state), the anchor
  smoothing, the blend, the stage tables of `nstates = 1`.
- The SNR 0.01 set (`/data2/multistate_bench`, the 20k subset S and the balanced subset B): a
  stress case of the earlier plan, not run here.
- The stream, the GUI, the 2D tables, the skills (skill text is proposed in the report).
- Any new parameter, any floor or timeout change, any edit of the multistate test's gates.

## 5. Change types and gates

| Phase | Kind | Gate |
| --- | --- | --- |
| 0 | measurement, no edit | Release build of the export, fast gate 15/15, `unit_pftc_align2D3D suite=sigma2_anchor`; the entry and the nightly recorded (known failures classified, not judged) |
| 1 | tests-only | fast gate (`check_test_registry.py` included); `multistate_workflows suite=states` once without the key and once with `sigma2_anneal=no`, both passing; the six profiles complete |
| 2 | measurement | six F05 verdicts or stops by the rule; decision 5 read |
| 3 | measurement | four refine3D_states verdicts or budget stops; decision 7 read |
| 4 | result-changing, `nstates > 1` only | `simulated_workflow_6vxx` identity at one thread (base export against the tree, the test environment's seed); fast gate; the touched unit suites (`unit_reconstruction`, `unit_pftc_align2D3D suite=sigma2_anchor`); the entry; three F05 runs by decision 3 |
| 5 | docs-only | the report exists; the plan is in `completed/`; the policies say what the code does |

## 6. Acceptance

Numbers and their origins (all frozen at the date of this plan):

- Separated: accuracy at least 0.90 and every recall at least 0.75; partial: accuracy at least 0.75
  or exactly two states recovered (measured: the recovery plan's runs; `acc3.py` defaults).
- Standing criterion: at least 2 of 3 F05 de novo runs separated (ruling 12 of the recovery plan).
- Stopping rule: stage-4 accuracy below 0.60 (measured: 0.765-0.995 after stage 4 in the five
  runs of 2026-11-11/13).
- Profile margin: 5 of 20 runs (guessed, binomial spread, decision 6). `shc` floor: 6 of 20
  (standing, `test_environment_policy.md`).
- refine3D_states on F05: 10 h wall per run (guessed budget, decision 7).
- Reported, not gated: the final maps docked on the truth maps in both hands (`dock_eval.sh`,
  `keep/runs/<run>/dock_eval.txt`), the accuracy at each stage end and each iteration
  (`metrics.txt`), the `SIGMA2 ANNEAL` readout of every iteration of arm Y (weight, beta, band,
  the mean log ratio of measured to anchor noise) tabulated per stage, the wall time and load per
  run. Hans reads the maps.
- No identity case in Phases 0-3 (no production code changes). Phase 4's identity case:
  `simulated_workflow_6vxx` at one thread with the test environment's seed, base export against
  the tree, compared with `test=identity target=simulated_workflow_6vxx dir1= dir2=` (about 13 min
  per build on dell2); no exclusions.

## 7. Common rules and stop conditions

- No new parameter. Experimental variants, if any, are patches on `git archive` exports under
  `scratch/keep/variants/`, never keys.
- No floor, tolerance, timeout or assertion of an existing test is relaxed. A floor never comes
  from the first run of the code being accepted.
- Fix as you go: a defect met on the way is fixed in the run only if it sits in a file of the
  table; otherwise it is recorded with its evidence and the run goes on. A file table row is
  added before any edit outside the table.
- Every benchmark run keeps its labels, log, statistics snapshots, metrics, docking and final maps
  under `scratch/keep/runs/<run>/` with one line in `scratch/keep/runs/summary.tsv` (run, arm,
  build, wall, accuracy, recall per state, ARI, populations, per-stage accuracy, verdict, notes);
  the rest of the run directory is deleted once that line is written. `SIMPLE_SEED` is unset for
  every benchmark run.
- Thread and disk policy: decision 10; the engine's disk floor (`DISK_MIN_GB`).

Stops (written up as a question for Hans with the evidence, never worked around):

- both arms fail the standing criterion (decision 5, third case);
- the F05 inputs or the truth are missing or do not verify (Phase 0);
- the export does not build or the fast gate fails on the base;
- a candidate cannot keep `nstates = 1` identical, or would need a new key;
- the machine cannot run what a phase needs (disk, a job the budget cannot hold).

Gates no ruling may waive: `nstates = 1` identity for a kept change; the `shc` floor of 6 of 20
for a kept change; the fast gate; no relaxed floor.

## 8. Phases

- **Phase 0, baseline** (about 4 h). Export the base (`git archive HEAD | tar -x -C
  scratch/keep/src_base`), build Release with tests into `scratch/keep/build_base`, install there,
  fast gate; record compiler, cores, load, build time. `unit_pftc_align2D3D suite=sigma2_anchor`.
  The entry `multistate_workflows nthr=8` once (default scoring), its `metrics.tsv` kept. The
  nightly (`test=nightly ncunits=16 dir=scratch/nightly`), its report kept, every failure classified
  as pre-existing (named in `nightly_report.tsv` and the memory of 2026-11-15: `pose_student`) or
  new, recorded, not judged. Verify the F05 inputs: `P2_05/P2_05.simple` with 80,000 `ptcl3D` rows,
  `labels_f05.txt` with 24,000 / 48,000 / 8,000 of states 1 / 2 / 3, `PT05/`, the truth maps under
  `truthcheck/`; copy and repoint the scripts and tools; run `acc3.py` on the kept labels of
  `dn_F05_428_1` (`~/agent_runs/ms_validate_4282d58/keep/runs/`) to show the tools reproduce
  0.927. Write down the annealing schedule the five-stage run will get (section 1) to check
  against the first Y run's `SIGMA2 ANNEAL` lines. Check `simple_exec prg=solve3D` and
  `prg=refine3D_states` accept `sigma2_anneal` (the UI's inputs). Write the time plan.
  Exit: build and gate recorded; the entry's `metrics.tsv` and the nightly report under
  `scratch/keep/phase0/`; the inputs verified by counts; the tools reproduce 0.927; the schedule
  written; the time plan written.
- **Phase 1, the cheap profiles** (about 14 h). Tests-only: the `multistate_workflows` routes that
  run `solve3D` (`run_solve3D`) and `refine3D_states` (`run_refine3D_states`) set `sigma2_anneal`
  on their command lines when the test's command line has it (as `prob_athres` is passed on to
  the `cavgs` route); `thread_args` of the profile runner forwards `sigma2_anneal=` as it does
  `suite=`; the policy row of `multistate_workflows` in `test_environment_policy.md` names the
  key. Rebuild, fast gate, then `suite=states` once without the key and once with `no`. Then the
  profiles of decision 6: `shc`, `states`, `consensus`, arm Y then arm N each, 20 seeds, under
  `scratch/profiles/<route>_<arm>/`, their `profile` tables kept. Table in Progress: per route and
  arm, the count resolving every state, the medians, the shuffled-label control count, and the
  reading by decision 6.
  Exit: the test change builds and passes the gate; six profile directories with their tables;
  the Progress table and reading.
- **Phase 2, F05 de novo** (about 42 h). Six runs by decision 4 with the watchers (stats, stage
  projects, the stopping rule), `dn_analyze.sh` after each, `dock_eval.sh` on the final maps, the
  `SIGMA2 ANNEAL` lines of the Y runs tabulated per stage (`grep` of `solve3D.log`). Check the first
  Y run's weights against Phase 0's schedule and record a mismatch as a finding (the run goes on).
  Exit: six lines in `summary.tsv`, each with its verdict; the reading by decision 5 in Progress,
  naming Phase 4's candidate or none.
- **Phase 3, refine3D_states on F05** (up to 40 h, likely much less). The consensus map once from
  `PT05` (`reconstruct3D nstates=1`; the `states` route's recipe in
  `simple_commanders_test_highlevel_msworkflows.f90`, `run_states_route`), then four runs by
  decision 7, each from a fresh copy, timed, the state column printed after each
  (`print_project_field`), scored with `acc3.py`, the maps docked.
  Exit: four lines in `summary.tsv` (verdict or *unfinished*), the reading by decision 7.
- **Phase 4, candidates** (only if Phase 1, 2 or 3 names one; about 21 h per candidate). Implement
  the candidate for `nstates > 1`, in the file of its row; the identity case, the fast gate, the
  unit suites, the entry; then three F05 runs (de novo for A1/A2, refine3D_states for B1), judged
  by decision 8. A candidate that fails is restored (`git show <base>:<path> > <path>`) and
  recorded with its mechanism before the next. A kept candidate updates its policy text in the
  same phase (`solve3D_policy.md` "Scoring-noise annealing" for A1/A2, `refine3D_states_policy.md`
  section 4 for B1, `refine3D_policy.md` 5.1 for either).
  Exit: a kept candidate with its evidence and policy text, or every named candidate judged and
  restored, in Progress.
- **Phase 5, report** (about 1 h). `multistate_sigma2_anneal_validation_report.md` beside the
  completed plan: what was learned first (does annealing keep, lose or gain states, per route, with
  the counts and the spreads), what changed, decisions under delegation, deviations, evidence
  against the baseline (the 1c12bb4eb and 4282d58a3 numbers), what is left for Hans (the maps;
  a regression bisect if decision 5's third case; the refine3D_states lead), skill text proposed,
  a commit message. Move the plan to `completed/` with `mv`.
  Exit: the report exists; the plan is in `completed/`; `git status` lists only files of the table.

### File table

| File | What | Phases |
| --- | --- | --- |
| `doc/implementation_notes/planned/multistate_sigma2_anneal_validation.md` | this plan and its Progress | each |
| `doc/implementation_notes/completed/multistate_sigma2_anneal_validation.md` | the plan, moved at the end | 5 |
| `doc/implementation_notes/completed/multistate_sigma2_anneal_validation_report.md` | the report | 5 |
| `src/main/commanders/test/simple_commanders_test_highlevel_msworkflows.f90` | pass `sigma2_anneal=` on to the solve3D and refine3D_states routes | 1 |
| `src/main/commanders/test/simple_commanders_test_highlevel_suite.f90` | `thread_args` forwards `sigma2_anneal=` to a profile's executions | 1 |
| `src/main/ui/simple_test/simple_test_ui_highlevel.f90` | the key as an input of the test program, if the parser needs it | 1 if needed |
| `doc/policies/test_environment_policy.md` | the `multistate_workflows` row: the key passed on | 1 |
| `src/main/solve/simple_solve3D_controller.f90` | candidate A1: the multi-state schedule reaches measured noise at the first prob stage | 4 |
| `src/main/solve/simple_solve3D_reseed_tester.f90` | the stage-schedule test of A1 (or the tester that owns `set_refine3D_sigma2_schedule`'s checks) | 4 |
| `src/main/commanders/simple/simple_commanders_solve3D.f90` | candidate A2: the default `no` for `nstates > 1` | 4 |
| `src/main/commanders/simple/simple_commanders_refine3D.f90` | candidate B1: refine3D_states' default | 4 |
| `src/main/simple_refine3D_states_controller.f90` | candidate B1, if the march needs it | 4 if needed |
| `src/main/ui/simple/simple_ui_solve3D.f90`, `src/main/ui/simple/simple_ui_heterogeneity.f90` | the UI default text of a kept A2 or B1 | 4 if needed |
| `doc/policies/3D/solve3D_policy.md` | a kept A1/A2: "Scoring-noise annealing" | 4, 5 |
| `doc/policies/3D/refine3D_policy.md` | a kept candidate: section 5.1 | 4, 5 |
| `doc/policies/heterogeneity/refine3D_states_policy.md` | a kept B1: section 4 | 4, 5 |

## Progress

(one entry per phase, written by the run: date, machine, base, files changed, each exit item and
how it was met, log paths under the run root, decisions taken under delegation, what it means in
short)
