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

A fix in Phase 1 adds its files here, with the reason.

## Progress

Nothing run yet.
