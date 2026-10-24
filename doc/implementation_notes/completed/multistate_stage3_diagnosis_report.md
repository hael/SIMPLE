# Multi-state stage-3 diagnosis: report (arun run `ms_stage3diag`, 2026-10-21 to 2026-10-22)

## The answer

There is no clear regression to find on Cyril's simulated two-state set. The anchor build
(1bbbdcc24, 2026-10-15) separated the two states in two of its four full `solve3D` runs (0.971,
0.931; then 0.768, partly separated, and 0.512, failed). The current build (f8195504b, stage
skipping removed for multi-state runs) separated in one of five (0.912 on P2H, its four runs on P2
failed). The anchor's first run was the lucky half of a coin that comes up heads about half the
time.

The anchor's early stages do not reliably end with better poses either. Two anchor runs ended
stage 3 with far better poses than any run of the current build (median pairwise pose error 9.0°
and 13.9° for one true state, 2.7° and 4.7° for the other, against 13-27°). The third ended stage 3
worse than all three runs of the current build (28.7° and 24.0°). By the decision rule fixed
before the runs, that is *not reproducible*.

What does decide a run, in both builds, is whether one state *locks on* during stage 3:

- **What locking on looks like.** One found state's map becomes right (docked onto the true map at
  correlation 0.93-0.97, its particles' pose error 2-6°). Its half-maps then agree out to the edge
  of the crop (5.6 Å), which promotes stage 4 to 6 Å.
- **What follows.** Of the five runs that locked on, three separated and one partly. Of the four
  full runs that did not, none separated.
- **How often each build locks on.** The anchor locked on in three of its four runs. The current
  build locked on in one of its seven runs on P2. That is a lead, not a finding: if both builds were
  equally likely to lock on, a split this uneven would still turn up about one time in eleven.

The anchor ran stage 2 to its 20-iteration budget every time, while the current build's scan exit
ended stage 2 after 11-19 iterations in six of its eight runs. That is the one systematic
difference visible before stage 3, and it is untested here.

No code defect was found and no code was changed.

## What changed, phase by phase

- **Phase 0** (`/home/elmlundho/agent_runs/ms_stage3diag/review/phase_0.diff`): the plan only.
  - The maintainer's ruling of 2026-10-21, under "Rulings that hold" and in Phase 0's exit.
  - Phase 0's Progress entry: the decisions fixed before the runs, the table of six runs after
    stages 2 and 3, the anchor runs to the end, the reference from the earlier full runs, and the
    exit items.
- **Phase 1** (`/home/elmlundho/agent_runs/ms_stage3diag/review/phase_1.diff`): the plan's Phase 1
  Progress entry (how runs of both builds fail, the lock-on table), a file-table row, this report,
  and the plan moved from `planned/` to `completed/`.

## Evidence against the baseline

- **Runs.** Six new runs on the same 2D project (P2) with the same command as the earlier runs,
  one at a time on an otherwise idle machine:
  - three of the current build, stopped after stage 3;
  - three of the anchor, run to the end.
- **What was measured.** The project was snapshotted after stages 2 and 3, and each snapshot's
  states were checked against the run's own states file. For each snapshot: the pairwise pose
  error against the true poses, accuracy, populations, half-map FSC, and each stage map docked onto
  both true maps.
- **The baseline.** The full runs already on this machine: the anchor's first run (`ms_regression`)
  and the current build's five (`ms_fixcheck`).
- **The tables** are in the plan's Progress: Phase 0 (both builds after stages 2 and 3, the anchor
  runs to the end) and Phase 1 (all twelve runs, lock-on and outcome).
- **Raw outputs** are under `/home/elmlundho/agent_runs/ms_stage3diag/scratch/keep/` (`runs/`,
  `phase0_table.md`, `tools/`), and the logs under `scratch/logs/`.
- **Times.** An anchor run took 3.4-4.5 hours to the end, and a current-build run 32-37 minutes
  through stage 3.

## Decisions taken under delegation

All are recorded in the plan's Progress, under Phase 0 unless noted.

- **Verdict with two true states.** The pose rule had to hold for each of the two true states for
  a build to count as better. It was fixed before the first run.
- **How pose error is measured.** Pairs are two particles of the same true state in the same
  found state, 200 000 draws with a fixed seed, using SIMPLE's rotation convention. The tool's
  self-test covers the true poses, a rotated frame, the other hand, small noise and random poses.
- **Where a stopped run ends.** It was let go on until the log showed the start of stage 4, so
  that stage 4's low-pass decision is recorded.
- **Stage maps.** They were padded to the true maps' 180-pixel box before docking, using the
  docking band of the earlier runs.
- **Builds.** The current build was built Release with `compile_clean.sh` (the fast gate passed
  15 of 15), to match the earlier runs. The anchor build was used in place.
- **Waiting for the machine.** No run started until the `ms_fixcheck` jobs that were using the
  machine had ended.
- **Phase 1 itself.** It ran no job, because the ruling's third case allows no probes or variants.

## Deviations from the plan

- **The maintainer's ruling** (2026-10-21, applied in Phase 0):
  - anchor runs went to the end;
  - the order was changed to current, anchor, alternating;
  - Phase 0 also answered whether the anchor separates reliably;
  - Phase 1's content followed from that answer, so no diagnostic variants, no bisection and no
    fix.
- **Two discarded starts.** The first anchor run was started twice without counting.
  - The first start ran the wrong binaries for a minute: the anchor's environment file was missing.
    The run script now checks which `simple_exec` it uses.
  - The second start ended with the session that had started it.
  - Neither produced numbers.
- **The plan's "commits"** are read as phases: nothing was committed.

## Left open for the maintainer

- **Whether to test the stage-2 scan exit.** The plan's diagnostic variant (stage 2 without the
  scan exit) could be judged cheaply by the lock-on marker after stage 3. Telling a three-in-four
  lock-on rate from one-in-seven would take about eight to ten variant runs, roughly five to six
  hours. The scan exit is a speed change the ruling keeps, so this would only be diagnostic.
- **How many runs a regression question needs.** Even the anchor fails one run in two on this set,
  so any future regression question here needs more than four runs per build.
- **The second failure point.** After locking on, the minority state sometimes does not rebuild in
  stage 4: the current build's fourth run, and partly the anchor's second. It is not tied to either
  build in these runs.
- **Runs that never find the structure.** One in each build (HF_3 and AS 3): no map becomes right
  even by the end. That is a de novo failure on this low-SNR set, not a state-separation failure.
