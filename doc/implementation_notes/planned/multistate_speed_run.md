# Multi-state speed: the run `ms_speed`

Implementation note, 2026-10-15. **Plan; nothing implemented.** The
execution plan of the unattended run `ms_speed` on the Oracle Linux test
machine. The design is in `single_state_speed_week.md` (beside this note):
its section "Multi-state" (the rulings and items 1 to 3) and items 3 and 4 of
its section "Day 3". This note holds the run's work, cases, rules, phases and
progress.

## Goal

Make multi-state runs faster for the release; multi-state matters even more
than single state (maintainer, 2026-10-15). The baseline is the three-state
msp1 run (141,838 particles, analysis in
`~/Processing/msp1/bench/bench_summary/analysis.md`): 10,415 s before the run
`prob_speed`. That run is estimated to bring it to about 7,000 s, and the
work of this run to about 5,000 s; these are estimates from msp1's reports,
not measurements. This run measures the gains on beta-galactosidase; the
maintainer measures msp1 afterwards.

## The work

1. **The multi-state schedule** (rulings of 2026-10-15; results change):
   - A multi-state `solve3D` switches from `prob` to `prob_neigh` with
     `prob_neigh_mode=state` at the first stage whose directions go from
     1,000 to 2,500. In the default schedule that is stage 5; the switch is
     derived from the stage table of directions (`NSPACE`), not from a second
     stage constant. `NSTAGES_INDEPENDENT` keeps the default stop after
     stage 5 with a definition of its own.
   - Directions are capped at 2,500 in multi-state `solve3D` runs (stages 7
     and 8, 5,000 today, when a user runs them) and in `refine3D_states`,
     whose `prob_neigh` stage ran at 5,000 (`STAGE2_NSPACE`). Its 500
     subspace directions (`STAGE2_NSPACE_SUB`) stay.
   - Stage 2 of a multi-state `solve3D` ends once `shc` scans more than 95%
     of the search space on average, whatever the state overlap. The present
     rule (state overlap and scan together) still ends it earlier if it is
     met. Stage 1 and every other user of the multi-state rule stay as they
     are.
     - The measure is the mean over the particles of the percentage of
       references scanned, which the search stores per particle as `frac`
       and the log prints as "% SEARCH SPACE SCANNED".
     - The controller tells the rule through an internal key on stage 2's
       `refine3D` command line, set together with `fracsrch=95`. The key is
       declared like `cohort_sampling`: parameters and parser, marked
       internal, no interface entry.
   - The strategy's benchmark report (`REFINE3D_STRATEGY`) records the mean
     percentage scanned, and the summary script shows it.
2. **`ref_assign` without rescans** (identical results), in the dense table's
   assignment:
   - Single state: count the available particles instead of testing all of
     them (`any(ptcl_avail)`) for every assignment. Advance only the
     references whose head was the assigned particle, as the neighbourhood
     assignment already does.
   - Several states: the state labels come from a global greedy frontier, and
     a particle leaves it at its smallest distance over all references. So
     its label is the state of its best reference. Compute that per particle,
     in parallel.
   - The two agree except when two references of different states tie
     exactly for a particle's minimum, for example a flat particle whose
     distances are all equal. The tie rule is the lowest reference index,
     written in the code and tested.
   - Within a state: loop over that state's references only, count the
     available particles, and advance only the heads that pointed at the
     assigned particle. The random draws and their order are unchanged.
3. **Table reads** (identical results):
   - The `prob` driver reads the part files in parallel; each part fills
     only its own particles' columns of the master table. The master table
     is not initialised in full, since every entry is overwritten and the
     read already stops on an incomplete part.
   - The neighbourhood table (`read_tabs_to_glob` of the `prob_neigh`
     table) reads each part file once, counting and filling in one pass,
     instead of twice.

The neighbourhood assignment (`ref_assign_neigh`), which stage 5 of a
multi-state run uses after item 1, already works on a sparse graph with
counters and is not changed.

## Common rules

- **Data.** Beta-galactosidase: `/data2/bgal_addon/control/bgal_addon.simple`
  (5,513 particles, box 256 at 1.275 Å per pixel, the 2D classification,
  `mskdiam=180 pgrp=d2`). It is copied into a fresh directory under the run's
  scratch directory for each run; the originals are never written.
- **Builds.** Timed runs use an out-of-tree Release build with tests built.
  The gates use the Debug build of `./compile_debug.sh`, which runs the fast
  gate.
- **Execution settings.** `nparts=4 nthr=6` for distributed runs and
  `nthr=24` for shared-memory runs, identical in every phase. One case at a
  time on the otherwise idle machine.
- **Timing.** Gains are reported as numbers from this machine.
  - Every timed comparison is the same case before and after, at least twice
    each way, interleaved, read with `scripts/plot_refine3d_bench.py`: wall
    time, the critical path per iteration, the table parts' thread-seconds,
    and the pre-step driver's table reads and `ref_assign`.
  - The long case MS-shc runs once each way.
  - A change is kept unless it is slower than the build without it beyond
    the run-to-run spread.

**Cases.** Exact command lines are recorded in Phase 0 and repeated
unchanged.

| case | starts from | command line | what it measures |
| --- | --- | --- | --- |
| S1 | the control project | `simple_exec prg=solve3D_cavgs projfile=bgal.simple mskdiam=180 pgrp=d2 nparts=4 nthr=6` | the common single-state start, kept |
| S3 | the control project | S1's line with `nstates=3` | the common three-state start, kept |
| A | S1's final project | `simple_exec prg=solve3D projfile=bgal.simple mskdiam=180 pgrp=d2 nsample=2500 cavg_ini_ext=yes nparts=4 nthr=6` | single state: the `prob` and `prob_neigh` tables, reads and assignment |
| MS | S3's final project | A's line with `nstates=3` | three states from the class averages: the particle stages that follow (stages 4 and 5 by default) |
| MS-shc | the control project | `simple_exec prg=solve3D projfile=bgal.simple mskdiam=180 pgrp=d2 nstates=3 nsample=2500 nparts=4 nthr=6` | three states from the particles, as msp1 ran: the `shc` stages 1 and 2, then `prob` |
| M | the final project of MS's first Phase 0 run | `simple_exec prg=refine3D projfile=bgal.simple nstates=3 refine=prob objfun=euclid nspace=1000 maxits=2` with `vol1` to `vol3` its final state maps, `box_crop`, `lp` and `trs` of MS's first particle stage, `mskdiam=180 pgrp=d2 nparts=4 nthr=6` | the dense three-state table at 3,000 references: its fill, the reads and `ref_assign` |

- **S3.** If S3 leaves a state empty, MS cannot start, since
  `cavg_ini_ext=yes` needs every state populated. Rerun it once with
  `reseed_states=yes` and record that; if a state is still empty, stop.
- **MS-shc.** With `nsample=2500` and three states the sample exceeds 90% of
  the particles, so every particle is searched in every iteration. msp1
  searched a fraction; the stage pattern is the same.
- **M.** It starts from the same kept project in every phase, so its runs
  stay comparable after the schedule changes MS.

**Acceptance.**

- **Single-state paths** in every phase:
  - The one-thread `simulated_workflow_1jxy` (both steps) and
    `simulated_workflow_6vxx` stay byte-identical to Phase 0, under the
    comparison rules of the run `prob_speed`. The run identifier in
    `solve3D_manifest.txt` and the layout digest in the header of
    `sigma2_state.bin` hash absolute paths, so they are ignored.
  - A's final FSC 0.143 lies within the larger of its Phase 0 spread and 3%
    of its Phase 0 value.
  - `solve3D_addon` and `cont_refine3D_1jxy` pass.
- **Identical changes** (phases 2 and 3):
  - The 3D table tester compares the new assignment with the present one,
    kept in the tester as the reference, on synthetic dense tables of one
    and of three states with a fixed seed. The assignments are identical,
    and the tie case follows the written rule.
  - M's state populations and last-iteration mean score stay within the
    spread of the before-runs; multi-threaded runs are not bit-identical.
- **The schedule change** (Phase 1):
  - MS and MS-shc run to the end, and every state keeps particles and a map.
  - Every state's final FSC 0.143 is no worse than 1.05 times the worst
    state's value over the Phase 0 runs of the same case.
  - MS-shc's stage 2 ends by the 95% rule: record the iteration and the
    percentage, against Phase 0's stage-2 length.
  - `flex_pca_blobs` and `solve3D_pyramid_blobs` pass their checks (or do
    no worse than on the base, if the base already fails one).
  - If a multi-state check fails, stop and describe it: the maintainer
    decides with msp1.
- **Gates.** The fast gate (`./compile_debug.sh`), the CTest unit suites of
  the touched areas, and `python3 scripts/check_test_registry.py .` whenever
  tests or their registration change.

**Not in this run.**

- The compact table record, the states' volume assembly in parallel and the
  worker fold (items 4 to 6 of the week plan's section "Multi-state").
- Day 3 items 1, 2, 5 and 6, and Day 4.
- The cap on refined references across states (open).
- Skill files: proposed text goes into the report.
- The generated indexes under `doc/code_overview/`.

The plans stay in `planned/`.

## Phases

### Phase 0: baseline

1. Build the tree as it is: the Debug build with the fast gate, and the
   Release build for timed runs, kept for the before-runs of Phase 1.
2. Run the cases on the Release build:
   - S1 and S3, kept;
   - A twice, MS twice, M twice, MS-shc once.
   - Record for every case:
     - the exact command line and the wall time;
     - the final FSC 0.5 and 0.143 per state, and the state populations;
     - the summary script's report; keep its benchmark reports.
   - Record MS's list of stages (refine mode, directions, iterations), and
     MS-shc's stage-2 length and the percentage scanned at its last
     iteration.
3. Run the tests on the Release build, keeping the outputs of the two
   `simulated_workflow` suites for the byte comparisons:
   - `simulated_workflow_1jxy` (both steps) and `simulated_workflow_6vxx`;
   - `solve3D_addon` and `cont_refine3D_1jxy`;
   - `flex_pca_blobs` and `solve3D_pyramid_blobs`.
   Record their results.

Exit criteria: the build and the fast gate pass; every case and test is
recorded with its spread where it ran twice; S1, S3 and M's start are kept.

### Phase 1: the multi-state schedule

Item 1 of "The work".

- Tests:
  - the convergence tester (`simple_convergence_tester.f90`) gets cases for
    the stage-2 exit: it ends at a scan above 95% with a low state overlap,
    and stage 1 does not;
  - the single-state acceptance;
  - `flex_pca_blobs` and `solve3D_pyramid_blobs`.
- Timed against Phase 0:
  - MS twice and MS-shc once;
  - S3 once, recording which of its stages changed (the class-average route
    shares the mode policy);
  - A once, as a check that single state did not move.
- Documents: `solve3D_policy.md`, sections 3 and 6 (stage 5 in `prob_neigh`,
  the cap, the stage-2 exit), the `refine3D_states` policy (the cap), and
  the controllers' comments.

Exit criteria: the gates pass; the acceptance of the schedule change holds;
the timings and MS-shc's stage-2 length are recorded against Phase 0.

### Phase 2: `ref_assign` without rescans

Item 2 of "The work".

- Tests:
  - the tester's reference comparisons;
  - the byte comparisons and the rest of the single-state acceptance.
- Timed against Phase 1: A and M, twice each way. Record the `ref_assign`
  seconds of the pre-step reports.
- If an assignment cannot be made identical apart from the written tie rule,
  stop and describe it.

Exit criteria: the gates pass; the assignments are identical; the timings are
recorded against Phase 1.

### Phase 3: table reads

Item 3 of "The work".

- Tests: the byte comparisons and the rest of the single-state acceptance.
- Timed against Phase 2: A and M, twice each way. Record the table-read
  seconds of the `PROB_ALIGN` and `PROB_ALIGN_NEIGH` reports.

Exit criteria: the gates pass; the outputs are identical; the timings are
recorded against Phase 2.

### Phase 4: documents and close

- `refine3D_policy.md`: the pre-step's reads and assignment.
- The week plan's section "Multi-state": mark items 1 to 3 done, with the
  measured numbers.
- Final runs:
  - A, MS and M once each, against Phase 0;
  - every simulated test, `flex_pca_blobs` and `solve3D_pyramid_blobs`.
- The progress section, then the report beside this note, with a table of
  every case: wall time, per-state resolution and populations in Phase 0 and
  now.

Exit criteria: the gates pass; the report is written.

### File table

| file | why | phases |
| --- | --- | --- |
| `doc/implementation_notes/planned/single_state_speed_week.md` | the status of the multi-state items | 4 |
| `scripts/plot_refine3d_bench.py` | the summary script: the scan percentage, fixes | each |
| `src/main/solve/simple_solve3D_controller.f90` | the multi-state schedule: the switch to the neighbourhood table, the cap, the stage-2 key | 1 |
| `src/main/simple_refine3D_states_controller.f90` | the direction cap of `refine3D_states` | 1 |
| `src/main/simple_convergence.f90` | the stage-2 exit | 1 |
| `src/main/simple_convergence_tester.f90` | tests of the stage-2 exit | 1 |
| `src/main/params/simple_parameters.f90`, `src/main/params/simple_parameters_parse.f90`, `src/main/params/simple_parameters_phases.f90` | the internal stage-2 key | 1 |
| `src/main/strategies/parallelization/simple_refine3D_strategy.f90` | the scan percentage in the strategy report | 1 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90` | `ref_assign`; the read and the master table's initialisation | 2 3 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_tester.f90` | the reference assignment and the comparisons | 2 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_neigh.f90` | the neighbourhood table read once | 3 |
| `src/main/commanders/simple/simple_commanders_prob.f90` | the parallel read of the parts | 3 |
| `doc/policies/3D/solve3D_policy.md` | the multi-state schedule | 1 4 |
| `doc/policies/heterogeneity/refine3D_states_policy.md` | the direction cap | 1 4 |
| `doc/policies/3D/refine3D_policy.md` | the assignment and the reads | 2 3 4 |

## Progress

Nothing yet.
