# Why the anchor's early stages end with better maps: a multi-state diagnosis

Planned 2026-10-21 for a short Dell run (arun run `ms_stage3diag`), base f8195504b. It follows the
regression work in `multistate_solve3D_regression.md`, whose runs on the Dell (`ms_regression`,
stopped, and `ms_fixcheck`) hold the data, the tools and the evidence this note starts from.

## What is known

Cyril's simulated set: 30 000 particles of two RNA polymerase conformations (17 000 of EMD-51275,
13 000 of EMD-51273), box 180 at 2.055 Å, SNR 0.1, CTF, random poses. The two true maps agree
closely, FSC about 0.8 out to 10 Å after docking, so the information that tells the states apart
lies at high resolution. Particles classified against the two true maps keep 99.5 % of their
states, so the set is separable in principle.

Runs of `solve3D nstates=2` from the same 2D project, distributed (`nparts=3 nthr=8`):

- **The anchor build** (1bbbdcc24, 2026-10-15) separated the states (accuracy 0.971). At the end
  of stage 3 one state's half-maps agreed at FSC=0.5 out to the crop's edge (5.6 Å) while the
  states were still unseparated (accuracy 0.51). Because of that, the stage-boundary FSC=0.5 rule
  ran stage 4 at 6.0 Å instead of its planned 7.0 Å, and the states split during stage 4
  (iterations 54-61), finishing in stage 5.
- **The fixed build** (stage skipping removed for multi-state, f8195504b): every run so far
  ended stage 3 with maps whose FSC=0.5 stayed coarser than stage 4's planned 7.0 Å (11.2 Å in
  the first run, with FSC=0.143 at the 7.4 Å matching limit), so stage 4 stayed at 7.0 Å. The states began to split only in stage 5, the first stage at 6 Å, and
  either ran out of iterations while splitting (accuracy 0.674) or one state drained (0.572).
  Times: about 87 minutes per run against 4.3 hours for the anchor.

The FSC code is the same in both builds (checked in `ms_fixcheck`). While the states are
unseparated each state's map is a mixture of two nearly identical structures, which both
half-sets reconstruct alike, so a state's half-map FSC then says how consistent the poses are, not
whether the states are separated. The anchor's stage-3 maps therefore point to better poses after
stage 3, which this run measures directly against the true poses.

The maintainer will not change resolution policy on the evidence of simulated data. This run
diagnoses what in the code makes stages 1-3 end with worse poses; it changes no ladder,
promotion, exit or iteration-budget rule.

## Terms

- **Stage**: one of the five refine3D executions of a multi-state `solve3D` run: stages 1-2 `shc`,
  stages 3-4 `prob`, stage 5 `prob_neigh` in state mode, each with a planned low-pass limit.
- **Accuracy**: the fraction of particles whose state equals the true state, under the better of
  the two matchings of found to true states (0.567 if all particles are in one state).
- **Pairwise pose error**: for pairs of particles in the same true state and the same found state,
  the difference between the angle of their relative rotation in the found poses and in the true
  poses; free of the frame and hand each found map chose (`pair_pose_error` in
  `src/main/commanders/test/simple_test_truth_metrics.f90`). Reported as the median in degrees and
  the fraction of pairs within 5°.
- **Stage snapshot**: the project file copied the moment `solve3D` writes `states_stage<N>.txt`
  after stage N, which holds the poses and states of that moment.

## Rulings that hold

From 2026-10-21: the speed changes stay (a diagnostic variant may undo one, a fix may not); NU
filtering at stage 5 stays; resolution-limit rules stay as they are (this run varies none);
distributed and shared memory must behave the same; stages 1-2 keep `shc` with sticky sampling.
No commits; bugs found are fixed in the run when they are code defects, as the phases say.

From 2026-10-21 ~16:20 (the maintainer's ruling during Phase 0): the fixed build had failed in all
four full runs of `ms_fixcheck` (0.572, 0.674, 0.515, 0.589), and its fourth run ended stage 3 the
way the anchor's first run did (one state's half-maps to 5.6 Å, stages 3 and 4 promoted, stage 4
at 6 Å) and still failed; only one anchor run had ever finished. So Phase 0 first asks whether the
anchor separates reliably. The three anchor runs of Phase 0 go to the end (snapshots after stages 2
and 3 still taken; final labels, final maps and per-stage reprojection images kept); the three
fixed-build runs stop after stage 3 as planned, and the fixed build's separation count comes from
the five full runs of `ms_fixcheck` (four on P2, one on P2H). Order FS, AS, FS, AS, FS, AS. The
anchor counts as *clearly better* when at least three of its four full runs (its first run and the
three new ones) separated (accuracy 0.90 or more). Phase 1 then depends on both answers:
anchor clearly better and the pose rule holds, Phase 1 as written; anchor clearly better but the
pose rule does not hold, probes and variants run to the end and are judged by separation (at most
four probes, two runs each, each probe's time recorded); anchor not clearly better, no probes: the
report says there is no clear regression on this set and describes, from the per-stage labels,
poses and images of all runs, where and how runs of both builds fail (including the fixed build's
third run, which never found the structure).

## How runs are made and judged

The command is the one of the earlier runs:
`solve3D pgrp=c1 mskdiam=270 nparts=3 nthr=8 nstates=2 print_states=yes` on a copy of the 2D
project made by the anchor (P2), one run at a time on the whole machine, without `SIMPLE_SEED`.
Runs here stop after stage 3: the run's process group is ended once `states_stage3.txt` exists.
Stage snapshots are taken after stages 2 and 3; each snapshot's states must equal the
`states_stage<N>.txt` it goes with.

Measured after stages 2 and 3 of every run: the pairwise pose error per true state; accuracy and
populations; each state's half-map FSC=0.5 and FSC=0.143 (from the log); the stage's snapshot maps
(`recvol_state0<s>_stage0<N>_lp.mrc`) docked onto the true maps in both hands, correlation; the
iterations and wall time of each stage.

A build *ends stage 3 better* when its pairwise pose error after stage 3 is lower in every one of
its runs than in every run of the other build (three against three); anything less is *not
reproducible* and is reported as such.

## Phases

### Phase 0: is the anchor's stage 3 reproducible?

1. Build f8195504b (the run's tree). The anchor's build from `ms_regression` may be used in place
   (read only).
2. Three anchor runs and three fixed-build runs through stage 3, alternating, with snapshots.
3. The same measurements for the stage-2 and stage-3 snapshots, plus, for reference, what the logs
   of the full runs already made say about stages 2 and 3 (the anchor's first run, the fixed
   runs).

Exit: the table of both builds after stages 2 and 3, and the decision: *the anchor ends stage 3
better* (Phase 1 localises it) or *not reproducible* (Phase 1 is skipped: write the report).
Under the ruling of 2026-10-21 (see "Rulings that hold") the anchor runs go to the end and the
exit also answers whether the anchor separates reliably; what Phase 1 does follows both answers.

### Phase 1: which change makes the difference

1. Diagnostic variants first, on the fixed build, two runs each through stage 3, only for the
   suspects the Phase 0 table points at:
   - **stage 2 without the scan exit** (a scratch patch restoring stage 2's earlier convergence
     rule): the anchor's stage 2 ran 20 iterations, the fixed build's 13-15. Diagnostic only; the
     scan exit is a speed change and stays.
   - **the prob table before 2026-10-15**: the commit pair e2673a5d7 against its parent (stage 3
     is the first prob stage).
2. If neither explains the difference: at most five bisection probes, two runs each through stage
   3, over the first-parent commits between 1bbbdcc24 and f8195504b that touch the path, in
   order: b57538b7b, e2673a5d7, 65bee3073, 57a447151, 13795c1d9, 04e520e78, 189e0a355, fb4e653ac,
   cd13708eb, 8c8d8b427 (these five built with 43de9b340's fix of the table-read race applied),
   43de9b340, 9531dee62, f4d6accff, aea5b1846, af13bc828, 0c9542e50, f078cf683, dafd4d046,
   24d50aabf, 497cc974b, 68278fa9a, 72930a3a4, f8195504b. Builds from cd13708eb up to 72930a3a4
   can omit stage 3 or end it early (stage omission, the mid-stage resolution exit); they get
   f8195504b's stage-skipping fix applied, so that every probe runs stage 3 in full and only the
   other changes differ.
3. For the change found: what it does to stages 1-3 (iterations, scan, pose error per iteration
   where the logs allow).
4. If it is a code defect (the code does not do what its commit intended), fix it at its cause in
   the run's tree, add the files to the file table, and verify: three runs through stage 3 (pose
   error as good as the anchor's) and two full runs (accuracy). If it is a speed change or a
   resolution rule, describe it with the evidence and stop there: the maintainer decides.

Exit: the change named and its effect shown (or "no single change" with the evidence), any defect
fixed and verified, the report written.

## Common rules

- Nothing is committed; git is read-only; other commits are built from `git archive` exports.
- One SIMPLE job at a time on the whole machine; record the load when timing.
- Never rerun a run to change a verdict; the decision rule above is fixed before the runs.
- New or changed OpenMP code: a tester at 24 threads, repeated, compared bit for bit with one
  thread. THROW macros never continued with `//&`; no `.and.` short-circuit on allocatable arrays;
  no goto.
- Bulky outputs are deleted once their numbers are in the journal; snapshots, logs and the stage
  maps used go to `scratch/keep/`.

### File table

| file | why | phase |
| --- | --- | --- |
| `doc/implementation_notes/planned/multistate_stage3_diagnosis.md` | this plan and its progress | each |
| `doc/implementation_notes/completed/multistate_stage3_diagnosis_report.md` | the report | 1 |
| `doc/implementation_notes/completed/multistate_stage3_diagnosis.md` | this plan, moved to the completed directory at the end of the run | 1 |

A fix in Phase 1 adds its files here, with the reason.

## Progress

### Phase 0 (done 2026-10-22)

**Answer.** Neither build is reliably better on this set. The anchor separated the two states in
two of its four full runs (accuracy 0.971 and 0.931; then 0.768, partly separated, and 0.512,
failed); the fixed build separated in one of five (0.912 on the anchor's 2D project rebuilt by the
current 2D code, called P2H; its four runs on P2 failed at 0.572, 0.674, 0.515 and 0.589). Under
the maintainer's ruling two of four is *not clearly better*. The pose question is *not
reproducible*: two anchor runs ended stage 3 with far better poses than any fixed-build run, but
the third ended it with worse poses than all three. Under the ruling Phase 1 therefore makes no
probes and no variants: it writes the report, which describes where and how the runs of both builds
fail.

**What decides a run, as far as these runs show.** In every run, both builds, stages 1 and 2 end
with poor poses (pairwise pose error after stage 2 18-32°, stage-2 maps docking onto the true maps
at correlation 0.1-0.6) and no state separation (accuracy 0.50-0.52). Stage 3 is where runs part ways. In the two
anchor runs that went on to separate, one found state's map became right during stage 3 (docked
correlation 0.93-0.97 onto the truth, its particles' pose error 2-6°, its half-maps agreeing to the
edge of the 132-pixel crop, 5.6 Å), which promoted stage 4 to 6 Å; the states then began to split
late in stage 4 and finished splitting in stage 5 (AS 2 ran out of iterations before it was done). The other found state still had near-random poses
(pose error 18-34°). In the three fixed-build runs and the third anchor run no map became right in
stage 3 (correlation 0.75 at best, pose error 13-29°), stage 4 stayed at 7 Å, and the third anchor
run failed with one map matching neither true structure. So a run's fate is largely settled by
whether one state's map locks onto the structure during stage 3, and both builds sometimes do and
sometimes don't. The fixed build's fourth run in `ms_fixcheck` (no snapshots, so no pose error) did
show the anchor-like stage-3 maps and the promotion and still failed, so locking on in stage 3 is
not sufficient either.

The table below lists each run after stage 2 and after stage 3. Pose error is given as the
median for true state 1 / true state 2, with the fraction of pairs within 5° in brackets.
Docked correlation is the better hand's correlation of each found map (m) onto each true map (t).
FS is the fixed build, stopped after stage 3; AS is the anchor, run to the end.

| run | stage | pose error | accuracy | populations | FSC=0.5 / 0.143 per state (Å) | docked correlation | iterations | wall (s) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| FS 1 | 2 | 28.8° (0.14) / 25.8° (0.17) | 0.503 | 17161 / 12839 | 10.9/9.7, 18.5/16.1 | m1t1 0.31 m1t2 0.36 m2t1 0.19 m2t2 0.31 | 19 | 822 |
| FS 1 | 3 | 26.5° (0.17) / 22.1° (0.20) | 0.501 | 17417 / 12583 | 9.0/7.7, 19.5/13.2 | m1t1 0.48 m1t2 0.50 m2t1 0.16 m2t2 0.15 | 17 | 1167 |
| FS 2 | 2 | 28.3° (0.15) / 23.6° (0.18) | 0.502 | 14287 / 15713 | 14.8/13.7, 12.8/11.9 | m1t1 0.32 m1t2 0.39 m2t1 0.26 m2t2 0.44 | 20 | 865 |
| FS 2 | 3 | 22.9° (0.21) / 12.8° (0.30) | 0.515 | 14684 / 15316 | 9.7/7.7, 20.6/8.8 | m1t1 0.73 m1t2 0.69 m2t1 0.62 m2t2 0.60 | 17 | 1295 |
| FS 3 | 2 | 31.6° (0.12) / 29.8° (0.13) | 0.510 | 15071 / 14929 | 12.3/10.6, 11.9/9.3 | m1t1 0.12 m1t2 0.12 m2t1 0.35 m2t2 0.16 | 12 | 490 |
| FS 3 | 3 | 24.0° (0.17) / 18.3° (0.21) | 0.514 | 16548 / 13452 | 8.8/7.4, 16.1/10.6 | m1t1 0.20 m1t2 0.63 m2t1 0.33 m2t2 0.14 | 17 | 1347 |
| AS 1 | 2 | 23.8° (0.18) / 18.1° (0.22) | 0.505 | 16544 / 13456 | 8.6/7.7, 9.0/7.4 | m1t1 0.46 m1t2 0.61 m2t1 0.10 m2t2 0.37 | 20 | 933 |
| AS 1 | 3 | 9.0° (0.45) / 2.7° (0.63) | 0.512 | 23546 / 6454 | 5.6/5.6, 18.5/12.8 | m1t1 0.95 m1t2 0.97 m2t1 0.16 m2t2 0.16 | 17 | 2389 |
| AS 2 | 2 | 25.7° (0.16) / 20.3° (0.19) | 0.510 | 13278 / 16722 | 13.7/13.2, 10.6/7.7 | m1t1 0.36 m1t2 0.32 m2t1 0.24 m2t2 0.40 | 20 | 930 |
| AS 2 | 3 | 13.9° (0.36) / 4.7° (0.51) | 0.510 | 19997 / 10003 | 5.6/5.6, 16.8/9.0 | m1t1 0.93 m1t2 0.97 m2t1 0.50 m2t2 0.62 | 17 | 2447 |
| AS 3 | 2 | 29.8° (0.14) / 26.6° (0.16) | 0.519 | 13255 / 16745 | 17.6/11.6, 9.3/7.6 | m1t1 0.14 m1t2 0.09 m2t1 0.07 m2t2 0.13 | 20 | 937 |
| AS 3 | 3 | 28.7° (0.15) / 24.0° (0.18) | 0.525 | 13128 / 16872 | 9.3/8.8, 7.9/7.6 | m1t1 0.30 m1t2 0.33 m2t1 0.08 m2t2 0.20 | 15 | 1743 |

Stage 4's low-pass limit: promoted to 6.0 Å in AS 1 and AS 2; the planned 7.0 Å in AS 3 and all
three FS runs. Every snapshot's states equal the `states_stage<N>.txt` it goes with.

The anchor runs to the end:

| run | final accuracy | verdict | populations | stage 4 / stage 5 (iterations, wall s, accuracy at the end) | final maps docked onto the truth | wall |
| --- | --- | --- | --- | --- | --- | --- |
| A_1 (`ms_regression`) | 0.971 | separated | 13414 / 16586 | 15, 4915, 0.86 / 12, 6699, 0.97 | 0.994, 0.991 | 4.27 h |
| AS 1 | 0.931 | separated | 13511 / 16489 | 15, 4621, 0.62 / 12, 8250, 0.93 | map 1 onto truth 2 0.997, map 2 onto truth 1 0.989 | 4.53 h |
| AS 2 | 0.768 | partial (still rising when stage 5 ended) | 16550 / 13450 | 15, 4133, 0.58 / 12, 6002, 0.77 | 0.992, 0.993 | 3.79 h |
| AS 3 | 0.512 | failed | 17713 / 12287 | 15, 3228, 0.50 / 12, 6353, 0.51 | map 1 onto both truths 0.79/0.82 (a mixture), map 2 onto neither (0.51/0.27) | 3.44 h |

For reference, the logs of the full runs made before this phase say the same about stages 2 and 3.
The anchor's first run (A_1) ended stage 3 with one state's half-maps at 5.6 Å and was promoted at
stage 4, like AS 1 and AS 2. Of the fixed build's runs in `ms_fixcheck`, HF_1 and HF_2 ended stage 3
coarser (FSC=0.5 at 10-11 Å for the better-resolved state) and were not promoted, like FS 1 to FS 3.
HF_3's better state reached 7.0 Å, but that run never found the structure, and it was not promoted.
HF_4 ended stage 3 at 5.6 Å, was promoted to 6.0 Å, and still failed (0.589). EF_1 separated
(0.912), splitting in stage 4. Accuracy after stages 2 and 3 was 0.50-0.52 in all of them.

**Exit items.**
- The table of both builds after stages 2 and 3: above. The pose error was measured for every
  stage-2 and stage-3 snapshot (six runs, twelve snapshots, all consistent with their states files).
  Accuracy, populations, FSC, docked correlation, iterations and wall time are given with it, plus
  the logs of the earlier full runs.
- The decision on the pose rule: *not reproducible*. For true state 1 the anchor's stage-3 medians
  are 9.0, 13.9 and 28.7° against the fixed build's 26.5, 22.9 and 24.0°; for true state 2 they are
  2.7, 4.7 and 24.0° against 22.1, 12.8 and 18.3°. The third anchor run lies above the fixed build's
  range for both states. At stage 2 the ranges also overlap (state 1: 23.8-29.8° against 28.3-31.6°;
  state 2: 18.1-26.6° against 23.6-29.8°).
- The separation question added by the ruling: the anchor separated in 2 of 4 runs, the fixed build
  in 1 of 5. The anchor is not clearly better, so Phase 1 makes no probes and writes the report.

**Evidence** (under `/home/elmlundho/agent_runs/ms_stage3diag/scratch/`):
- Per run: `keep/runs/<run>/`, holding `solve3D.log`, `run.log`, `stage_log.txt`, and for each of
  `snap2/` and `snap3/` the project snapshot, states, `pose_err.txt` (per true state and per
  found-state cell), `acc.tsv`, `dock.txt` and the gzipped particle dump. The anchor runs also hold
  `metrics.txt` (per-stage and per-iteration accuracy), `final_states.txt`, the final maps,
  `dock_summary.txt`, `pose_err_final.txt`, `stage_log_full.txt` and the per-stage images in `jpg/`.
- The table: `keep/phase0_table.md`.
- Logs: `logs/run_<run>.log`, `logs/collect_<run>.log`, `logs/queue_phase0.log`, and the build log
  `logs/build_head_full.log`.
- Tools: `keep/tools/` (`pose_err.py` with its self-test, `run_stage3.sh`, `collect_stage3.sh`,
  `collect_final.sh`, `stage_log.py`, `table.py`, `acc.py`, `run_metrics.py`, `queue_phase0.sh`).
- The runs started at 16:19 on 2026-10-21 and ended at 05:56 on 2026-10-22, one at a time with the
  machine otherwise idle. Two starts of the first anchor run do not count (described under
  "Order" below): `logs/run_AS_1_aborted_*`.

Files changed in this phase: this plan only (the Phase 0 Progress entry and the ruling under
"Rulings that hold" and in Phase 0's exit). No code changed.

Decisions taken before the first run, so that no number of this phase can move them:

- **How the decision rule reads with two true states.** Pairwise pose error is measured per true
  state, so each run gives two medians after stage 3. A build *ends stage 3 better* only when, for
  each of the two true states, the median of every one of its three runs is lower than the median
  of every run of the other build. If the two true states disagree, or any run of one build falls
  inside the range of the other, the result is *not reproducible*. The fraction of pairs within 5°
  and the other measurements are reported beside it but do not enter the verdict.
- **Pairs.** For a true state, a pair is two different particles of that true state in the same
  found state (found state 0 excluded): the first particle is drawn uniformly over the true state,
  the second uniformly over the particles of the same true state in the first one's found state;
  200 000 draws with a fixed seed, so the same snapshot always gives the same number. The rotation
  matrices are SIMPLE's (`euler2m`). The tool's self-test gives zero error for the true poses, for
  the true poses in a rotated frame and in the other hand, about 1.8° for 2° of noise on each Euler
  angle and about 36° for random poses.
- **Where a run ends.** After `states_stage3.txt` appears the run is let go on until its log shows
  the start of stage 4 (seconds later), so that the low-pass limit chosen for stage 4 (planned, or
  promoted by the FSC=0.5 rule) is in the log; then the run's processes are ended.
- **Stage maps.** The snapshot maps are at the stage's crop box (112 px at stage 2, 132 px at
  stage 3). Before docking they are padded in Fourier space to the truth's 180 px box (`scale
  newbox=180`), then docked onto each true map in both hands (`dock_volpair hp=40 lp=10`, the band
  the earlier runs used); the better hand's correlation is reported.
- **Builds.** The fixed build is this run's tree (f8195504b plus this plan), built Release with
  `compile_clean.sh` like the builds of the earlier runs (fast gate 15 of 15 passed); the anchor is
  the earlier run's Release build of 1bbbdcc24, used in place.
- **Order.** First planned as AS 1, FS 1, AS 2, FS 2, AS 3, FS 3; under the maintainer's ruling of
  2026-10-21 it is FS 1, AS 1, FS 2, AS 2, FS 3, AS 3, with the anchor runs going to the end. One
  run at a time on the whole machine, starting only once the `ms_fixcheck` run that was using the
  machine at the start of this phase had ended (16:12). Two starts of the first anchor run do not
  count and gave no numbers: the first ran the wrong binaries for a minute (the anchor's
  environment file was missing; the run script now refuses to start unless the build's own
  `simple_exec` is the one found), the second was ended in stage 2 together with the session that
  had started it.

### Phase 1 (done 2026-10-22)

**What this phase did.** Phase 0 found the anchor not clearly better (two of four runs separated)
and the pose rule not reproducible. That is the third case of the maintainer's ruling, so this
phase made no probes and no diagnostic variants, and ran no SIMPLE job and no build. It describes,
from the kept per-stage labels, poses, maps, logs and images of all twelve runs, where and how runs
of both builds fail. The twelve runs are the anchor's four (A_1 from `ms_regression`, AS 1 to 3),
the fixed build's five full runs from `ms_fixcheck` (HF_1 to HF_4 on P2, EF_1 on P2H) and its three
runs through stage 3 (FS 1 to 3). Then it wrote the report and moved this plan to the completed
directory. No code defect was found, so nothing was fixed.

**Where runs part ways: whether one state locks on during stage 3.** In every run stages 1 and 2
end alike: no state separation (accuracy 0.50-0.52), poor poses (pairwise pose error 18-32° in the
six runs with snapshots), and stage-2 maps that dock onto the true maps at correlation 0.6 at best.
The images of the stage-1 maps of a run that separated (AS 1) and one that failed (AS 3) look
alike. During stage 3 a run either *locks on* or it does not. Locking on means that one found
state's map becomes right: in the runs with snapshots it docks onto the truth at 0.93-0.97, its
particles' pose error falls to 2-6°, and its population grows (to 18 000-25 000 of the 30 000
particles) while the other state keeps near-random poses. In the log it shows as that state's
half-maps agreeing out to the edge of the 132-pixel crop (5.61 Å at FSC=0.5), which in turn makes
the stage-boundary rule promote stage 4 to 6 Å. A run that does not lock on ends stage 3 with
neither map right (correlation 0.75 at best, pose error 13-29°), populations between 12 500
and 17 500, FSC=0.5 of the better state at 7-11 Å, and stage 4 at its planned 7 Å. The reprojections show the
difference too: after stage 3, AS 1's locked state shows the molecule's outline and domains, while
AS 3's maps are streaky blobs, and in AS 3 they stay that way to the end. The twelve runs:

| run | build, input | stage-2 iterations | locked on in stage 3 | final accuracy |
| --- | --- | --- | --- | --- |
| A_1 | anchor, P2 | 20 (budget) | yes | 0.971 separated |
| AS 1 | anchor, P2 | 20 (budget) | yes | 0.931 separated |
| AS 2 | anchor, P2 | 20 (budget) | yes | 0.768 partial |
| AS 3 | anchor, P2 | 20 (budget) | no | 0.512 failed |
| HF_1 | fixed, P2 | 13 (scan exit) | no | 0.572 failed |
| HF_2 | fixed, P2 | 15 (scan exit) | no | 0.674 failed |
| HF_3 | fixed, P2 | 11 (scan exit) | no | 0.515 failed |
| HF_4 | fixed, P2 | 17 (scan exit) | yes (stage 3 itself promoted to 7.0 Å) | 0.589 failed |
| EF_1 | fixed, P2H | 20 (budget) | yes | 0.912 separated |
| FS 1 | fixed, P2 | 19 (scan exit) | no | stopped after stage 3 |
| FS 2 | fixed, P2 | 20 (budget) | no | stopped after stage 3 |
| FS 3 | fixed, P2 | 12 (scan exit) | no | stopped after stage 3 |

Three of the five runs that locked on separated, one partly. None of the four full runs that did
not lock on separated. The anchor locked on in three of its four runs, the fixed build in one of
its seven runs on P2 (and in its one run on P2H). With these counts that difference could still be
chance: if lock-on were equally likely in both builds, a split at least this uneven would turn up
about one time in eleven. It is a lead, not a finding.

The six snapshotted runs give one more hint. Ranked by pose error after stage 2, the order is the
same for both true states: AS 1, AS 2, FS 2, FS 1, AS 3, FS 3. The two runs that locked on are the
two whose stage-2 poses were best (24/18° and 26/20° against 28-32° and 24-30° for the others).
But the gaps are small, and FS 2 ran its full 20 iterations of stage 2 without locking on.

**How runs fail.**
- *The structure is never found well* (HF_3 and AS 3, one in each build). No map is right after
  stage 3 (AS 3: docked correlation 0.33 at best, pose error 29/24°). At the end one map is a
  mixture of the two states and the other matches neither true map (AS 3: correlation 0.51 and
  0.27; HF_3: at most 0.39). This is a de novo failure, not a state failure, and it happens in the
  anchor as well.
- *No lock-on, then a state drains or splits too late* (HF_1, HF_2). Stage 4 at 7 Å leaves the
  states unseparated (accuracy 0.51-0.52). In stage 5, at 6 Å, one state drains in HF_1 (to 2 934
  particles, 0.572). In HF_2 the states begin to split at the fifth iteration of stage 5 and are
  still splitting when its twelve iterations end (0.674).
- *Lock-on, but the second state does not rebuild in stage 4* (HF_4; partly AS 2). In the runs
  that separated, the minority state grew back from about 4 000-5 000 to 12 000-13 000 particles in the
  first iterations of the promoted stage 4, and the split followed (AS 1 and EF_1). In HF_4 the
  minority state stayed near 7 500 throughout stage 4 (accuracy 0.52), and stage 5 got only to
  0.589. AS 2's minority state reached about 10 700 and its split began late in stage 4; it was
  still rising (0.768) when stage 5's iteration budget ran out.

**What this suggests, for the maintainer** (no change made: these are speed and schedule rules
the ruling keeps):
- The one systematic difference visible before stage 3 is stage 2's length. The anchor ran stage 2
  to its 20-iteration budget in all four runs. In the fixed build the scan exit ended it after 11-19
  iterations in six of its eight runs. Lock-on is a yes/no marker visible about 35 minutes into a
  fixed-build run, so the plan's diagnostic variant (stage 2 without the scan exit) could be judged
  cheaply on lock-on. Telling a three-in-four rate from a one-in-seven rate reliably would take
  about eight to ten runs per arm, roughly five to six hours for the variant.
- Because the anchor itself fails one run in two, the regression question needs more than four
  runs per build to settle, whatever variant is tried.
- The second failure point, the minority state not rebuilding after lock-on, is not tied to either
  build in these runs.

**Exit items.**
- *The change named, or "no single change", with the evidence:* there is no clear regression to
  localise. The anchor separates in two of four runs and the fixed build in one of five, so
  bisection was not run (per the ruling). The evidence is the Phase 0 tables and the lock-on table
  above.
- *Any defect fixed and verified:* no code defect was found; nothing was fixed.
- *The report written:* `doc/implementation_notes/completed/multistate_stage3_diagnosis_report.md`.

**Evidence:** Phase 0's run directories `scratch/keep/runs/{AS,FS}_{1,2,3}/` (snapshots, pose
errors, per-iteration accuracy in the anchor runs' `metrics.txt`, the images in `jpg/`). The full
runs of the earlier runs are in `~/agent_runs/ms_fixcheck/scratch/keep/runs/` (HF_1 to HF_4, EF_1)
and `~/agent_runs/ms_regression/scratch/keep/runs/A_1`. The analysis commands and outputs are in
the Phase 1 journal (`state/journal_phase_1.md` of the run directory).

Files changed in this phase: this plan (the Phase 1 entry, a file-table row, moved to
`doc/implementation_notes/completed/`) and the new report.
