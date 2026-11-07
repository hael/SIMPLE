# Multi-state posterior state draw: report of the run `ms_posterior2`

Report, 2026-11-06. The run carried out the plan
`doc/implementation_notes/completed/multistate_posterior_state_draw.md` on
dell2: the home build machine, with 16 cores of three kinds, Ubuntu 26.04 and
gfortran 15.2. It ran unattended, as one uncommitted diff against HEAD
ea4c1c20a, the tree the batch alignment pipeline refactoring left. Each phase's
diff is `/home/hael/agent_runs/ms_posterior2/review/phase_N.diff` (written by
the run's driver at the end of each phase); `R` below is
`/home/hael/agent_runs/ms_posterior2`. The plan's Progress section holds the
detail and the evidence paths; this report is the summary.

The change concerns multi-state `solve3D` (de novo 3D reconstruction with
`nstates > 1`) and the probabilistic refinement it runs:
- A particle's state is now one draw from its tempered posterior over state and
  direction together, in place of the closest state followed by a draw within
  it.
- The labels stay hard, and the weighted insertion's extra poses come from the
  drawn state.
- The table fill splits its shift-refinement budget over the states.
- The neighbourhood table draws over every state's candidates.
- The annealing tail waits until no state gains resolution.
- A final assignment at the physical posterior labels every particle against
  the final maps.

No parameter was added.

## What we learned

- **The change works as designed, and single-state runs are untouched.**
  Single-state results are bit for bit identical to the base tree: at one and
  eight threads, at two parts, and on a full single-state workflow. The unit
  suites, the fill harness at 16 threads, the multi-state test routes and the
  nightly suite (19 of 19 entries) pass.
- **On small simulated data it behaves like the rule it replaced.** On true
  poses with random labels (the new `labels` test route, four minutes) the
  joint draw recovers the three beta-galactosidase states along the same path
  as the closest-state rule. The state posterior there is almost as sharp as an
  argmin from the second iteration. Over 20 random starts it resolves every
  state 13 times. In the other 7 a state drains as the maps of a random start
  fall behind. At the three draining seeds checked the old rule drains too,
  and so do hotter draws (`prob_beta` 0.1 at all three, 0.03 at two). The de
  novo `shc` route resolves every state in 15 of 20 seeds,
  against the floor of 6 (13 of 20 in October, on older code).
- **On the three-state benchmark subset it does not separate the states.**
  The subset is 20,000 simulated ribosome-assembly particles: 30 % the smallest
  assembly, 60 % the medium, 10 % the largest.
  - With the true maps and poses every particle goes to its true state (the
    state signal is there at 10 Å).
  - From one shared random labelling on true poses, two runs end at accuracy
    0.69 and 0.65 (accuracy: the share of particles in the found state matched
    to their true state).
    The smallest assembly separates; the medium is split over two states; the
    largest never gets a state and merges into one of them.
  - Three de novo runs end at accuracy 0.59, 0.42 and 0.58, all failed by the
    thresholds fixed beforehand. Two give a clean map of the medium assembly
    (correlation 0.97 to the truth) holding half to 60 % of its particles. The
    largest assembly is never isolated.
- **The readout places the failure in the low-resolution probabilistic
  stages.** In stage 3 (22.6 Å) the state posterior at the search's factor 0.3
  spreads over 2.2-2.7 of the three states (1.6-2.3 in stage 4, at 20.2 Å).
  Accuracy stays near chance while about half the labels change each
  iteration, and a quarter to over half of each state's posterior mass rests
  on entries the fill did not shift-refine. The labels start to organise at the
  end of stage 4 (it adds 0.09-0.10 of accuracy in two runs) and through the
  last stage at 6 Å, which adds 0.03-0.10.
  - Stage 3 meets the plan's test for "too hot".
  - But even at the physical factor the state posterior is broad (1.3-1.7
    states). That fits maps that cannot yet tell the states apart at that
    resolution as well as it fits a draw that is too hot.
  - The existing annealing tail at the end of the last stage moved accuracy by
    less than 0.01.
  - So the plan's fallback, an annealing tail at the end of every
    probabilistic stage, is the named next step, but these runs do not show it
    would help.
- **The final assignment works and is cheap.** It labels all 20,000 particles
  at the physical posterior against the final maps in 8-12.5 minutes, changes
  1-3 % of the labels, and the coverage check after it finds nothing missing.

## What changed, phase by phase

- **Phase 0** (2026-11-04; `review/phase_0.diff`), the builds, the test route
  and the tools:
  - the new `labels` route of `multistate_workflows`, in the default route
    list, with its variants `labels_identical` and `hard_labels`;
  - the registry comment and the test policy;
  - the subset's 2D solution, P2S;
  - the tools, copied and adapted.

  The phase ran before the batch alignment pipeline refactoring, which the
  maintainer ruled should go first (ruling 12 of the plan's section 10); its
  diff was re-applied on the tree that refactoring left (ruling 13).
- **Phase 1** (2026-11-05; `review/phase_1.diff`), the joint draw and its
  readout:
  - the dense table's joint state-and-direction draw in the per-row selection
    `assign_rows`;
  - the state readout: `PROB_TAB_STATE_DRAW`, `PROB_TAB_STATE_MASS` and four
    statistics-file keys, carried in the parts' readout records and stored on
    the convergence object in both execution modes;
  - the unit tests, the single-state identity, and the `labels` reports.
- **Phase 2** (2026-11-05; `review/phase_2.diff`), the rest of the design:
  - the per-state split of the fill's refined-reference budget on both shift
    paths;
  - the neighbourhood table's pooled draw, with the extras gathered from the
    committed state;
  - the per-state stall test of the annealing tail;
  - the final assignment, replacing the greedy missing-update pass;
  - the fill harness's per-state polish check;
  - the identity again.

  An independent review found two defects, both fixed. First, the distributed
  master re-seeded and checked the multi-state labels at every `refine3D`
  start; a state drained in one stage would then stop the next stage, and the
  final assignment. It now does that only at a fresh start, as the
  shared-memory strategy does. Second, a neighbourhood particle whose
  candidates all carry a non-finite distance got no reference. That case
  crashed the single-state `simulated_workflow_6vxx` on the binary before the
  fix; such a particle is kept again, as before.
- **Phase 3** (2026-11-05; `review/phase_3.diff`), robustness on small data:
  - the `labels` profile and its floors, written in the suite;
  - the `shc` profile;
  - the single routes, the full entry and the nightly suite;
  - the entry's new time in the registry comment and the test policy.
- **Phase 4** (2026-11-05/06; `review/phase_4.diff`, the plan only),
  validation on the subset:
  - the truth-frame check and the known-reference control;
  - two label-recovery runs and three de novo runs, judged against the truth.
- **Phase 5** (2026-11-06; `review/phase_5.diff`), documents and report:
  - `solve3D_policy.md`, `refine3D_policy.md` and `refine3D_states_policy.md`
    updated;
  - two path references follow the plan's move;
  - this report, and the plan moved to `completed/`.

## Decisions taken under delegation

Each is recorded in the plan's Progress section of its phase and in that
phase's journal (`R/state/journal_phase_N.md`); the journal holds the reasons
in more detail.

1. *How `labels` starts* (Phase 0). Shared-memory `refine3D` stops without
   starting maps, so the route draws the uniform labelling itself, reconstructs
   the maps, and names them as the distributed master's fresh start does. Its
   settings mirror `solve3D`'s multi-state probabilistic stages.
2. *The `consensus` route's failure* in Phase 0's default entry (a state merged
   and emptied) was recorded, not fixed. It is independent of Phase 0's diff,
   which changed only test code: at one thread the metrics are byte-identical
   on base and tree, and other eight-thread runs pass with a partial merge. On
   the state draw's code the route passed in Phases 2 and 3.
3. *No per-stage tail in Phases 1-2.* On `labels` the readout read "too cold"
   (the state posterior near-sharp from the second iteration), not "too hot",
   so the fallback was not indicated.
4. *The single-state identity on `simulated_workflow_6vxx`* instead of
   `_1jxy`: at this HEAD `_1jxy` stops before any probabilistic stage
   (ruling 13).
5. *The polish rule kept* (Phase 2). The polar shift map records a polish only
   where it beats the stored shift's distance. Changing that would change
   single-state results, so the harness checks each state's share and the
   budget instead (a deviation, below).
6. *The file table followed the work*:
   - the 3D matcher added under ruling 13 (the readout's route to the caller
     now that the pre-step is gone);
   - the refine3D strategy's row extended to Phase 2 for the distributed label
     start;
   - the fill harness and its driver added (Phase 2);
   - the speed-week closeout added for its link (Phase 5).

   The rows of the retired pre-step files were removed.
7. *The `labels` drains are the plan's risk 4, not a stop* (Phase 3). Risk 4
   is a state that drains because its map is garbage at any temperature. At the
   three draining seeds checked the closest-state rule drains too, and so do
   hotter draws (`prob_beta` 0.1 at all three, 0.03 at two). No allowed setting
   cures it, and reseeding is outside the plan. The route's floors were
   therefore set from its profile by the floor rule. For a measured value the
   floor sits no higher than the worst run less one robust spread, and no
   higher than the median less three. For a success count it is the count less
   three binomial standard deviations:
   - every state resolved in 6 of 20 seeds;
   - adjusted Rand index 0.70, accuracy 0.75, population error 0.30;
   - rotation error 2.0°, shift error 0.10 px.
8. *`hard_shc`'s collapse checked against the base build* (Phase 3): the base
   collapses the same way at the same seed.
9. *The entry's time updated* (Phase 3), and the build reinstalled before the
   nightly ran, so the nightly ran on the final tree.
10. *The subset runs' settings* (Phase 4). The label-recovery start was drawn
    once (seed 20261106) and applied with `selection`. The `refine3D` runs added
    `iterstats=yes` and `mkdir=no`, which change no result. The verdict
    thresholds kept their defaults: the known-reference ceiling was 1.0000, and
    they were fixed before the first de novo run.
11. *The per-stage tail not tried in Phase 4*, although stage 3 meets the "too
    hot" test:
    - ruling 10: the subset validates and does not drive development;
    - the fast tests do not call for it;
    - a controller change would reopen the `shc` profile and the nightly;
    - the evidence does not show it would help.

## Deviations from the plan

- **The design landed in the batch pipeline.** Sections 2-4 of the plan
  describe the pre-step `ref_assign` / `prob_align_exec`; under ruling 13 the
  draw went into the per-row selections `assign_rows` and `assign_rows_neigh`,
  every draw keyed, and the readout into the parts' records. Ruling 13 also
  widened the identity to eight threads and two parts.
- **The fill's polish coverage.** The plan asks that "every active block
  carries a polished entry per particle". The code guarantees each block its
  share of screen and polish attempts; 1 of 64 harness particles ended without
  a recorded polish in one state.
- **Phase 1's fixture known-reference run** was not the known-label control it
  was meant to be: the particle import resets states, so it started from random
  labels. Corrected in Phase 2's Progress.
- **The `shc` route was run in Phase 2** beyond the plan's list, because it is
  the one route that reaches the final assignment.

## Evidence against the baseline

Ruling 10 forbids a comparison arm of the current code on the subset, so the
comparisons are on small data and are named as such.
- **Single-state identity**, base build (HEAD exported with `git archive`)
  against the tree, same seed and layout: every map, FSC, sigma2 store,
  statistics file and particle record is identical. That holds for a small
  `refine3D` at one thread, at eight, and at two parts, and for
  `simulated_workflow_6vxx` at one thread (181 of 181 declared outputs), in
  Phases 1 and 2.
- **`labels` on the unchanged pipeline tree against the joint draw** (Phase 1):
  final accuracy 0.999 against 0.999, along the same per-iteration path. The
  conformational pair (`hard_labels`) stays unseparated (accuracy 0.67 before,
  0.69 after, below the 0.80 of one state for all), and identical states still
  drain a state.
- **The budget split on `labels`** (Phase 2): the unrefined share of each
  state block's mass, 0.64 / 0.40 / 0.25 at the last iteration of Phase 1's
  annealing tail, now stays at or below 0.061, with the labels unchanged.
- **The fill's cost** (Phase 2): +1.7 % table-fill time against Phase 1's
  binary, within the tenth the plan allows.
- **`shc` profile**: 15 of 20 against 13 of 20 on 2026-10-19 (older code, not
  paired). The new `labels` profile: 13 of 20.
- **Drains and collapses** found in Phase 3 occur on the base build too: the
  closest-state rule drains at the three of the seven draining `labels` seeds
  checked, and `hard_shc` collapses at its seed.
- **The benchmark**: the first run's de novo run of the code before this work
  on the full 80,000 particles (the plan's section 9) failed in the shape of
  two of Phase 4's runs: one medium-locked state and mixtures, though without
  the smallest-rich state those two had. It is not paired (other data size,
  older code).

## Open for the maintainer

1. **The per-stage annealing tail** (the plan's section 6 fallback; the
   controller would set the existing key `prob_beta_end=1` on stages 3-5). Stage
   3 on the benchmark meets the test that calls for it; the evidence does not
   show it would help. Trying it means a controller change, then the `labels`
   and `shc` routes, the `shc` profile and the nightly, then new de novo runs.
   If it does not suffice, a separate state factor would be a new parameter,
   which ruling 4 leaves to you.
2. **The polish budget per state.** In stages 3-4 a quarter to over half of
   each state block's posterior mass rests on stored-shift losses. One polish
   count for every block (the plan's section 4.2) would cost up to 2.5-2.8
   times the table time.
3. **The rare state on true poses.** From random labels, the partition that
   forms at the flat factor splits the dominant class. Whether the row's
   in-plane bias plays a part is not settled by these runs. That is the plan's
   risk 3: each row entry holds the loss at one drawn in-plane angle, not the
   loss marginalised over the angle. Its remedy, the in-plane marginal, changes
   what the row holds, which is your decision.
4. **The `labels` drains** (7 of 20 random starts lose a state). Reseeding an
   emptied state is outside the plan, and nothing map-driven is allowed
   (ruling 7).
5. **`consensus` is borderline at eight threads.** It failed once in Phase 0
   (a merged and emptied state) and passed in Phases 2 and 3; a profile of the
   current route is pending (test policy).
6. **The stream and the GUI** need no edit. The stream's multi-state `solve3D`
   stage runs `solve3D` and inherits the change, and nothing there reads the
   deleted pass or the new statistics keys. The GUI could show the state
   readout; nothing was done.
7. **The skill text** below.

## Proposed skill text (not edited)

- `.github/skills/simple-refine3d/references/bayesian-model.md`, the
  `ref_assign` bullet: "`assign_rows` (the per-row selection of the batch
  pipeline) draws each particle's state and projection together: one draw from
  the tempered posterior over its whole row, every active state's block,
  `p_r ∝ exp(-prob_beta T_i (L_r - L_min))` (`draw_posterior_row`); the state is
  the drawn reference's. With one state the row is the block. The
  neighbourhood table draws the same way over its pooled candidates. The
  extras of the weighted insertion come from the drawn state alone. With
  several states the refined-reference budget is split over the states'
  blocks."
- `.github/skills/simple-solve3d-importance-sampling/references/coupling-map.md`:
  replace "the state label is the deterministic argmin over states" with "the
  state is drawn with the projection from the whole row"; and add "a
  multi-state `solve3D` ends with one assignment of every active particle at
  `prob_beta=1` against the last stage's maps (`run_multistate_final_assignment`;
  solve3D's greedy missing-update pass is gone, while `refine3D_states` and
  `classify3D_refs` keep their own)".
- `.github/skills/simple-main-root/references/solve_controllers.md`, item 9:
  replace "(a multi-state run first assigns the particles no stage updated)"
  with "(a multi-state run first assigns every active particle at
  `prob_beta=1` against the last stage's maps)".
- `.github/skills/simple-main-strategies-probabilistic/SKILL.md`. It still
  describes the two-phase pre-step (table workers, `prob_align*_exec`, part
  files, the master's `ref_assign`), which the batch pipeline retired. It needs
  a rewrite to the pipeline (`fill_tab_range`, `assign_batch`, `assign_rows`,
  `assign_rows_neigh`, the readout records), with one line on the state rule:
  "state and direction are one keyed draw over the whole row; the state readout
  (`prob_state_readout`) travels in the readout record."
- `.github/skills/simple-refine3d/SKILL.md` and
  `.github/skills/simple-solve3d-importance-sampling/SKILL.md`: drop
  `prob_align_exec`, `prob_align_neigh_exec`, `simple_prob_align_strategy.f90`
  and `simple_commanders_prob.f90`, which the pipeline retired. Name the
  matcher's batch pipeline in their place.

## Commit message

```
multistate solve3D: draw the state from the posterior, final assignment at beta 1

A particle's state is now one draw from its tempered posterior over state and
direction together (dense and neighbourhood tables); labels stay hard and the
weighted insertion's extras come from the drawn state. The fill splits its
refined-reference budget over the states, the annealing tail waits until no
state gains, and a final assignment at prob_beta=1 replaces the greedy
missing-update pass. A state readout (PROB_TAB_STATE_DRAW/_MASS, STATS_FILE
keys) shows the draw's sharpness. New multistate_workflows route `labels`
(true poses, random labels), gated by its 20-seed profile. nstates=1 is bit for
bit unchanged. On the three-state benchmark subset the states are still not
separated de novo; see the report.
```
