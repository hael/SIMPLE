# Review report, run `ms_fixcheck` (2026-10-21): the stage-skipping fix on the two-state set

Plan: `multistate_solve3D_regression.md` beside this report. This run made one check, Phase 6 of
that plan's progress; the section "Phase 6, 2026-10-21: the stage-skipping fix on the two-state set"
holds every number. The journal is `/home/elmlundho/agent_runs/ms_fixcheck/state/journal_phase_6.md`.

## The answer

**The fix works as specified, but by itself it does not restore the anchor's behaviour. The fixed
build separated the two states in 1 of 5 runs, at about 1.5 hours per run.**

The fix stops a multi-state `solve3D` from omitting a stage or ending one early on resolution
grounds. None of the five runs omitted a stage or took the mid-stage exit; stages 3 and 4 ran until
their own convergence rule or iteration limit. On Cyril's simulated set:

- **On P2** (the 2D project made by the anchor build, the input of the earlier pair): 0 of 4
  separated. Final accuracy 0.572, 0.674, 0.515 and 0.589; all *failed* against the threshold fixed
  in Phase 0 (separated at 0.90 or more, partial from 0.75 to 0.90).
- **On P2H** (the 2D project made by the current build, i.e. Cyril's whole pipeline): the one run
  *separated*, accuracy 0.912. Each of its two maps docks onto its own true map at correlation 0.99.

"Accuracy" is the fraction of particles whose final state is their true state. Putting every
particle in one state scores 0.567; a random split scores about 0.51.

For comparison: the current build before the fix (H_1) failed at 0.573 in 48 minutes, and the anchor
build of 2026-10-15 (A_1) separated at 0.971 in 4.3 hours. A fixed run takes 83 to 100 minutes. That
is about 40 minutes more than H_1, because stages 3 and 4 now run in full. It is about a third of the
anchor's time, because the speed changes made after the anchor remain.

**Where the states split.** When a run separates, it does so the way the anchor did. In A_1 and
EF_1 the labels stay at chance through stages 1 to 3, and in stage 4, at 6.0 Å, the smaller state
first shrinks and then the states split from about the fifth iteration. Accuracy then gains 0.025 to
0.04 per iteration, and stage 5 completes the separation. The failed fixed runs show nothing in
stage 4:

- HF_2 starts the same kind of split only in the fifth iteration of stage 5, and runs out of
  iterations at 0.674 while still gaining.
- HF_4 gains a little in stage 5.
- HF_1 loses its second state in stage 5, as H_1 did.
- HF_3 never finds the structure (see "Unexpected" below).

## What changed

One phase, diff in `/home/elmlundho/agent_runs/ms_fixcheck/review/phase_6.diff`.

- **The fix itself**, in `simple_solve3D_controller.f90`, `simple_solve3D_utils.f90`,
  `simple_commanders_solve3D.f90`, `simple_solve3D_reseed_tester.f90` and
  `doc/policies/3D/solve3D_policy.md`, is the maintainer's, received uncommitted. It built as
  received and its unit test passed, so the run did not change these files.
- **The plan:**
  - the Progress section of Phase 0 of the stopped run `ms_regression`, carried over unchanged with
    its correction under "The simulated set";
  - the new Phase 6 section;
  - two rows of the file table: the fix's files and this report, both phase 6.
- **This report.**

## Decisions taken under delegation

All are recorded in the journal.

1. **Phase 0's Progress carried into this plan.** This copy of the plan still said "Nothing run
   yet". The Phase 0 text existed only in the stopped run's working tree. It was applied unchanged,
   with a line saying where it comes from, so that Phase 6 has its baseline in the same note.
2. **Inputs read in place.**
   - P2, P2H and the particle stack stay in the stopped run's directory. They are read only and
     copied per run, as the tools do.
   - The two truth maps and the true labels were copied into this run's `scratch/keep`.
   - The tools were copied and their paths adapted. To check them, H_1's metrics were recomputed
     from its kept labels, and came out identical.
3. **Run names.** HF_1 to HF_4 are the fixed build on P2; EF_1 is the fixed build on P2H.
4. **Per-stage images kept from HF_4 on.** HF_3's maps turned out to match nothing. From then on the
   collector also kept each run's small per-stage reprojection images. HF_1 to HF_3 have only their
   final maps.
5. **Two read-only diagnoses that change no verdict.**
   - A reading of the code: why the anchor ran stage 4 at 6.0 Å and HF_1 to HF_3 at 7.0 Å.
   - A docking and projection check of HF_3's maps. The docking ran at two threads while HF_4 was
     running.

No run was repeated, and no threshold moved.

## Deviations from the plan

- The plan's file table puts the report under `completed/`. The run's instructions put it here,
  beside the plan, which stays under `planned/` because its phases 1-5 have not been carried out.
- Otherwise none. The five runs used arm H's command exactly, one at a time, alone on the machine.

## Evidence against the baseline

Each cell gives accuracy after the stage, iterations and wall time.

| stage | H_1 (before the fix) | A_1 (anchor) | HF_1 | HF_2 | HF_3 | HF_4 | EF_1 (P2H) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 0.52, 12 it, 1 min | 0.51, 12 it, 2 min | 0.51, 12 it, 1 min | 0.52, 12 it, 1 min | 0.53, 12 it, 1 min | 0.52, 12 it, 1 min | 0.53, 12 it, 1 min |
| 2 | 0.51, 20 it, 14 min | 0.51, 20 it, 16 min | 0.51, 13 it, 9 min | 0.51, 15 it, 10 min | 0.52, 11 it, 7 min | 0.50, 17 it, 12 min | 0.50, 20 it, 14 min |
| 3 | omitted | 0.51, 17 it, 45 min | 0.51, 17 it, 21 min | 0.51, 17 it, 21 min | 0.52, 15 it, 18 min | 0.51, 17 it, 27 min | 0.50, 17 it, 21 min |
| 4 | 0.53, 3 it, 6 min (exit) | 0.86, 15 it, 82 min | 0.52, 15 it, 29 min | 0.51, 15 it, 29 min | 0.51, 15 it, 30 min | 0.52, 14 it, 34 min | 0.78, 15 it, 37 min |
| 5 | 0.57, 12 it, 25 min | 0.97, 12 it, 112 min | 0.57, 12 it, 26 min | 0.67, 12 it, 26 min | 0.51, 12 it, 26 min | 0.59, 12 it, 26 min | 0.91, 12 it, 26 min |
| stage 4 limit | 7.0 Å | 6.0 Å | 7.0 Å | 7.0 Å | 7.0 Å | 6.0 Å | 6.0 Å |
| run | **failed** 0.573, 48 min | **separated** 0.971, 256 min | **failed** 0.572, 86 min | **failed** 0.674, 88 min | **failed** 0.515, 83 min | **failed** 0.589, 100 min | **separated** 0.912, 100 min |

Where the evidence is:

- Labels after every stage and iteration, logs, final maps, confusion matrices and `summary.tsv`:
  `/home/elmlundho/agent_runs/ms_fixcheck/scratch/keep/runs/`.
- Build: `scratch/logs/build_head_full.log`. Release build, fast gate 15 of 15 passed.
- Unit test: `scratch/unit/unit_reconstruction.log`. 356 of 356 assertions passed, including the
  stage-rules suite with `test_stage_skipping`.

## Unexpected

- **HF_3 never found the structure.**
  - Its two final maps match neither truth map, nor each other. Docked correlation is at most 0.39,
    also at 20 Å and in both hands; the two maps correlate at 0.12.
  - Their projections are streaky blobs.
  - Its half maps still agreed to 5.4 Å in stage 4: the model is consistent but wrong.

  This is a failure of the poses, not only of the states. It happened in one of the six runs of the
  current code on this set. Where it began is unknown, because its intermediate maps were not kept.
  Stages 1 and 2, where a de novo model forms, are not touched by the fix.
- **The smaller state's map can be a striped blob after stage 3.** In HF_4 and in EF_1 it held 5 000
  to 6 000 particles whose labels barely moved: in HF_4's stage 3, 18 particles changed state over
  the last 14 iterations. In both runs the map recovered into a proper structure in stage 4. HF_4
  failed and EF_1 separated, so this alone does not decide the outcome.
- **Stage 4 ran at 6.0 Å in the anchor run, HF_4 and EF_1, and at 7.0 Å in HF_1 to HF_3.** The limit
  depends on whether the FSC=0.5 promotion fires. FSC=0.5 is the resolution at which the two
  half-maps' correlation falls to 0.5.
  - HF_1 to HF_3 reached only about 11 Å after stage 3, so no version of the rule promotes.
  - A reading of the code found no change since the anchor in how the half-map FSC is measured.
  - The rule itself changed on 10-18, from the measured resolution to one rung.
- No relabel or empty-state message, no crash and no NaN. Every run printed the two known
  "NSTATES is not given" warnings of the power-spectrum step at stage 1, which the stopped run
  traced as harmless.

## Left open for the maintainer

- **Keep the fix on its own merits.** It does what the ruling asks and costs about 40 minutes a run
  on this set.
- **It is not the whole story.** On the anchor's 2D project the fixed build failed 4 of 4, so the
  stage omission and the early exit were not the only difference from the anchor. The new paragraph
  in `solve3D_policy.md` section 6 presents the two rules as the cause of the 2026-10-20 failure. That
  is still a plausible part of it, but these runs do not show they were sufficient.
- **What separates a run that splits in stage 4 from one that does not is open.** Five runs do not
  settle it. The candidates the evidence leaves are:
  - **Stage 4 at 6.0 Å.** Both separated runs had it, but so did HF_4, which failed.
  - **The stage-2 scan exit of 10-16.** The two separated runs, and H_1, ran stage 2 for its full 20
    iterations. The four fixed P2 runs stopped after 11 to 17.
  - **The greedy state assignment of stages 3-4** (the plan's cause 4). The stages run their full
    length and gain nothing.
  - **The 2D project.** The one separated fixed run came from P2H, but one run cannot tell this from
    chance.

  The plan's Phase 2 diagnostic variants would test these directly: the schedule before 10-16, which
  restores a 20-iteration stage 2, and the prob table before 10-15. So would Phase 4, a temperature on
  the stage 3-4 state assignment (ruling 6). A few more EF runs would say whether P2H matters.
- **The anchor's own rate rests on one run.** The stopped run planned four anchor runs; only A_1 was
  made. Two or three more, at 4.3 hours each, would say how reliably the anchor separates this set.
  Without them, the anchor's 1 of 1 cannot be compared firmly with the fixed build's 1 of 5.
- **HF_3's pose failure** deserves a look if it recurs. A de novo run that converges to a consistent
  wrong model is a different failure from a state that drains.
