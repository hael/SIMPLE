# Multi-state speed: the run `ms_speed`

Implementation note, 2026-10-15. The execution plan of the unattended run
`ms_speed` on the Oracle Linux test machine. The design is in
`single_state_speed_week_completed.md` (beside this note): its section
"Multi-state" (the rulings and items 1 to 3) and items 3 and 4 of its section
"Day 3". This note holds the run's work, cases, rules, phases and progress.

**Status, 2026-10-20: split.** Done: Phase 0 (baseline), Phase 1 (the
multi-state schedule, accepted by the maintainer's ruling of 2026-10-16),
Phase 2 (`ref_assign` without rescans: implemented and identical on the tester
and on one-thread runs; it stopped at two noisy multi-threaded checks, and the
code is applied in the tree), and Phase 3 (table reads, carried out by the run
`speed_mac`, see `identical_speed_run.md`). Open: the ruling on the Phase 2
stop, and Phase 4 (documents, final runs, the report). They are in
[`../planned/multistate_speed_run_unfinished.md`](../planned/multistate_speed_run_unfinished.md).
Where the text below says "plan" or "nothing implemented", read it with this
status.

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
  - Each state's final FSC 0.143 and population are recorded against
    Phase 0. The earlier rule (no worse than 1.05 times the worst state of
    Phase 0) was dropped by the maintainer's ruling of 2026-10-16, after
    Phase 1: on beta-galactosidase the states that missed it are noise
    classes of a dataset without real states, and the schedule is judged on
    msp1. A multi-state case is accepted when it runs to the end and every
    state keeps particles and a map.
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

The plans stayed in `planned/` while work remained; the open part is now the `_unfinished` note.

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

Item 3 of "The work". Carried out instead by the run `speed_mac` on
2026-10-17, with identical results (`identical_speed_run.md`, Phase 1).

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
| `doc/implementation_notes/completed/single_state_speed_week_completed.md` | the status of the multi-state items | 4 |
| `scripts/plot_refine3d_bench.py` | the summary script: the scan percentage, fixes | each |
| `src/main/solve/simple_solve3D_controller.f90` | the multi-state schedule: the switch to the neighbourhood table, the cap, the stage-2 key | 1 |
| `src/main/simple_refine3D_states_controller.f90` | the direction cap of `refine3D_states` | 1 |
| `src/main/commanders/simple/simple_commanders_solve3D.f90` | the comment on the default multi-state stop, which said it comes before `prob_neigh` | 1 |
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

### Phase 0, 2026-10-15/16: baseline

No source file changed; this section is the only edit. The base is HEAD
75eab7f90 with a clean tree. Paths below are relative to the run
directory `/home/elmlundho/agent_runs/ms_speed` unless absolute.

**Machine and builds.**

- The machine is dell: Oracle Linux 8.10, 24 cores, no GPU. It uses GCC
  15.2 from gcc-toolset-15, with no environment modules. The login shell
  points `SIMPLE_PATH` at another checkout, so every run sets
  `SIMPLE_PATH`, `PATH` and `SIMPLE_QSYS=local` itself.
  `/usr/bin/python3` is 3.6, so a `python3` link to 3.12 goes first on
  `PATH` for builds.
- The Release build is `scratch/keep/build_phase0`
  (`scratch/build_release.sh`):
  `cmake -S repo -B <dir> -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS=ON -DPython3_EXECUTABLE=/usr/bin/python3.12`,
  then `make -j24` and `make install` (the install prefix is the build
  directory). It is kept for Phase 1's before-runs.
- The Debug build is `./compile_debug.sh`, run by
  `scratch/run_debug_build.sh`. The fast gate passed: 15 of 15 unit
  suites, return code 0 (`scratch/logs/debug_build_phase0.log`).

**How a case runs.**

- `scratch/run_cases.sh <case>...` runs the named cases one after another
  on the Release build. For each case it:
  - makes a fresh `scratch/cases/<case>/` holding a copy of the input
    project;
  - runs the command line below with `projfile=bgal_addon.simple mkdir=yes`
    (the plan writes `bgal.simple`; the copy keeps its file name), under
    `/usr/bin/time -v`;
  - keeps the case's `BENCH/`, the summary of
    `scripts/plot_refine3d_bench.py`, the final project, the
    `solve3D_manifest.txt` (the stage record that `solve3D` writes) and
    the log, under `scratch/keep/bench/<phase>/<case>/`.
  - `scratch/run_case.sh` does the single run. Phase 0 ran the same
    script under the name `run_phase0_cases.sh`; the copy adds only the
    S3 restore for M.
- A kept project names its maps by absolute paths inside the directory
  that wrote them. So before A, MS and M, the script recreates
  `scratch/cases/S1`, `scratch/cases/S3` or `scratch/cases/MS_r1` from
  `scratch/keep/`.
- The kept starts (`keep/S1`, `keep/S3`, `keep/M_start`) were trimmed
  after use. Removed: per-part and density volumes, reprojection models,
  `_unfil` and post-processed copies, and per-iteration plots. Every path
  the kept projects reference still exists. The logs show that A and MS
  read only the project's maps (their noise spectra are seeded afresh),
  so the trimmed starts are equivalent inputs.
- Metrics are read with three scripts in `scratch/keep/`:
  - `ms_metrics.sh` gives each state's last post-processing FSC 0.5 and
    FSC 0.143 from the log, and the state populations of the final
    project;
  - `prestep_sums.py` sums the pre-step reports;
  - `plot_refine3d_bench.py` (the repository's summary script) gives the
    rest.

**Command lines** (exact; later phases repeat them unchanged, each with
`projfile=bgal_addon.simple mkdir=yes` appended):

| case | input project | command line |
| --- | --- | --- |
| S1 | copy of `/data2/bgal_addon/control/bgal_addon.simple` | `simple_exec prg=solve3D_cavgs mskdiam=180 pgrp=d2 nparts=4 nthr=6` |
| S3 | the same | `simple_exec prg=solve3D_cavgs mskdiam=180 pgrp=d2 nparts=4 nthr=6 nstates=3` |
| A | `scratch/cases/S1/1_solve3D_cavgs/bgal_addon.simple` | `simple_exec prg=solve3D mskdiam=180 pgrp=d2 nsample=2500 cavg_ini_ext=yes nparts=4 nthr=6` |
| MS | `scratch/cases/S3/1_solve3D_cavgs/bgal_addon.simple` | `simple_exec prg=solve3D mskdiam=180 pgrp=d2 nsample=2500 cavg_ini_ext=yes nparts=4 nthr=6 nstates=3` |
| M | `scratch/cases/MS_r1/1_solve3D/bgal_addon.simple` | `simple_exec prg=refine3D nstates=3 refine=prob objfun=euclid nspace=1000 maxits=2 vol1=/home/elmlundho/agent_runs/ms_speed/scratch/cases/MS_r1/1_solve3D/recvol_state01.mrc vol2=…/recvol_state02.mrc vol3=…/recvol_state03.mrc box_crop=144 lp=7.96 trs=5.29 mskdiam=180 pgrp=d2 nparts=4 nthr=6` (full line in `scratch/M_args.txt`) |
| MS-shc | copy of the control project | `simple_exec prg=solve3D mskdiam=180 pgrp=d2 nstates=3 nsample=2500 nparts=4 nthr=6` |

M's values come from the first particle stage (stage 4) of the kept MS
run, as its `solve3D_manifest.txt` records it: `box_crop` 144, low-pass
limit 7.960975647 Å (written as `lp=7.96`) and shift limit 5.294117451
pixels (written as `trs=5.29`). Its maps are the run's final state maps
(`recvol_state0N.mrc`, the maps the project records). S3 left no state
empty, so it did not need a rerun with `reseed_states=yes`.

**Results** (Release build; one case at a time on an otherwise idle
machine, in the order S1, S3, A, A, MS, MS, M, M, MS-shc). FSC is the
Fourier shell correlation of the half-maps: the 0.5 and 0.143 values are
the resolutions in Å where it crosses those thresholds. Populations are
the particles in each state of the final project.

| case | wall time (s) | peak memory of the master | FSC 0.5 / 0.143 per state (Å) | populations |
| --- | --- | --- | --- | --- |
| S1 | 182.2 | 3.07 GB | 7.77 / 6.80 | 5,513 |
| S3 | 422.7 | 5.26 GB | 8.82 / 7.96, 27.20 / 23.31, 6.94 / 5.44 | 2,699, 917, 1,897 |
| A, run 1 | 920.0 | 3.24 GB | 4.35 / 3.84 | 5,513 |
| A, run 2 | 926.1 | 3.22 GB | 4.35 / 3.84 | 5,513 |
| MS, run 1 | 1,504.0 | 5.35 GB | 5.93 / 4.08, 17.18 / 9.89, 9.07 / 5.63 | 3,017, 783, 1,713 |
| MS, run 2 | 1,501.3 | 5.36 GB | 6.04 / 4.08, 16.32 / 10.88, 9.07 / 6.04 | 2,994, 794, 1,725 |
| M, run 1 | 63.6 | 0.94 GB | 6.40 / 4.53, 17.18 / 10.20, 9.33 / 6.80 | 3,394, 816, 1,303 |
| M, run 2 | 61.0 | 0.94 GB | 6.40 / 4.53, 17.18 / 10.20, 9.60 / 7.42 | 3,378, 816, 1,319 |
| MS-shc | 2,045.2 | 5.41 GB | 4.47 / 3.84, 12.55 / 8.59, 20.40 / 13.60 | 4,566, 672, 275 |

- **Spreads.**
  - A's two runs differ by 6.1 s (0.7%) and not at all in resolution.
  - MS: by 2.7 s; in resolution by up to 0.99 Å (state 2); in the
    populations by up to 23 particles.
  - M: by 2.6 s; by 16 particles; and in state 3's FSC 0.143 by 0.62 Å.
    The last iteration's mean score, on the [0,1] scale the log prints,
    is 0.410 in both runs.
- **Thresholds the later phases use.**
  - A's FSC 0.143 must stay within 3% of 3.84 Å (3% is larger than the
    zero spread), so at most 3.96 Å.
  - The worst state's FSC 0.143 over the Phase 0 runs is 10.88 Å for MS
    and 13.60 Å for MS-shc. Every state must therefore stay at or better
    than 11.42 Å in MS and 14.28 Å in MS-shc.
- **Where the time goes.** These are critical-path seconds from the
  summary script, summed over the run (run 1 of each case).
  - In MS, of 1,450 s of refinement, the probability tables take 1,158 s
    (80%), the search 111 s and volume assembly 42 s. The parts' table
    thread-seconds are 22,568 (22,513 in run 2).
  - In M, of 62.8 s, the tables take 44.5 s; table thread-seconds are 773
    (805 in run 2).
  - In MS-shc, of 1,987 s, the tables take 1,425 s; table thread-seconds
    are 27,189.
  - In A, of 870 s, the tables take 310 s, the search 196 s and volume
    assembly 177 s. Table thread-seconds are 2,461 for `prob` (the dense
    table) and 2,171 for `prob_neigh` (the neighbourhood table); run 2 has
    2,502 and 2,178.
- **The pre-step driver's table reads and `ref_assign`** (the assignment
  of particles to references), in seconds per iteration. They come from
  the `PROB_ALIGN` and `PROB_ALIGN_NEIGH` reports, summed by
  `prestep_sums.py`, because the summary script does not report them yet.
  - Three states × 1,000 directions: reads 0.61, `ref_assign` 0.35
    (MS, both runs, and MS-shc); M 0.69 / 0.40 and 0.66 / 0.37.
  - Three states × 2,500 directions: reads 1.50, `ref_assign` 0.85 to
    0.86 (MS); MS-shc 1.49 / 0.85.
  - Single state, A: `prob` at 1,000 directions 0.105 / 0.050; at 2,500
    directions 0.27 / 0.12. `prob_neigh` at 2,500 directions 0.021 /
    0.053; at 5,000 directions 0.037 / 0.092.
- **MS's stages** (refine mode, directions, iterations): stage 4 `prob`
  at 1,000 directions, iterations 1-15; stage 5 `prob` at 2,500
  directions, iterations 16-27, about 82 s per iteration, of which 69 s
  is the tables. Then comes the final noise-spectrum iteration (29) and
  the final reconstruction.
- **MS-shc's stages.** Stage 1 is `shc` (stochastic hill climbing, which
  stops at the first better reference) at 500 directions, iterations
  1-12. Stage 2 is `shc` at 1,000 directions, iterations 13-32. Stages 3
  and 4 are `prob` at 1,000 directions, iterations 33-49 and 50-64.
  Stage 5 is `prob` at 2,500 directions, iterations 65-76.
- **MS-shc's stage 2** ran its full 20 iterations (its iteration limit)
  without converging.
  - Its last iteration had a state overlap of 0.985 and scanned 95.79% of
    the search space on average.
  - The scan by iteration: 28.7, 30.0, 35.9, 47.1, 56.7, 67.0, 75.6,
    82.2, 87.0, 90.2, 92.6, 94.0, 94.8, 94.6, 94.7, 95.0, 95.2, 95.4,
    95.5, 95.8%.
  - It first exceeded 95% at iteration 29 (95.23%), the 17th iteration of
    the stage.
- **S1 and S3's stages** (class averages, used by Phase 1 to see which
  stages change): `prob_neigh` at 500 directions (iterations 1-30),
  `prob` at 500, 1,000 and 2,500 directions, then `prob_neigh` at 2,500
  and at 5,000 directions (iterations 85-96).

**High-level tests** (CTest entries on the Release build, after the
timed cases; logs `scratch/logs/ctest_<test>_build_phase0*.log`):

- `simulated_workflow_1jxy` passed in 1,008 s and
  `simulated_workflow_6vxx` in 733 s.
  - Kept for the byte comparisons under
    `scratch/keep/sim/phase0/{1jxy/7_solve3D,1jxy/8_solve3D_distr,6vxx/7_solve3D}`
    (`scratch/keep/keep_sim.sh`): the final maps, half-maps, FSC,
    noise-spectrum state, a text dump of the `ptcl3D` orientations and
    the manifest.
  - `scratch/keep/compare_sim.sh` compares two such directories. It
    ignores the manifest's run identifier, layout digest, checksum and
    artifact lines (they hash absolute paths) and the layout digest in
    the header of `sigma2_state.bin`.
  - Both suites ran a second time on the same build (849 s and 613 s;
    `scratch/run_ctests.sh` now runs `ctest -V`, so their checks are in
    the logs). The second run is byte-identical to the first under these
    rules, for 1jxy's two steps and for 6vxx, so the byte comparison is
    usable as the acceptance.
  - The checks of the second run were:
    - 1jxy shared-memory: docking correlation 0.9653 (mirrored hand),
      whole-volume correlation 0.8863, FSC 0.143 at 4.62 Å;
    - 1jxy `nparts=2`: 0.9652, 0.8837, 4.71 Å;
    - 6vxx: 0.9797 (direct hand), 0.8763, 4.62 Å.
- `solve3D_addon` passed in 454 s and `cont_refine3D_1jxy` in 447 s.
- **`solve3D_pyramid_blobs` fails on the base**, the same way in two runs
  (50 s).
  - Its step `solve3D_cavgs nstates=3` on 64-pixel particles builds a
    resolution ladder with a cropped box of 88 pixels. After "generating
    diverse labelling", `reconstruct3D` stops: "box_crop cannot exceed
    the native box" (`simple_parameters_phases.f90`, line 825).
  - This happens before any `refine3D` stage. Phase 1 has to do no worse
    than this.
- **`flex_pca_blobs` fails on the base**, the same way in two runs
  (180-196 s).
  - In `refine3D_states_cavgs`, after its `prob_neigh` refinements, flex
    PCA (the principal-component analysis of conformational variability
    that labels the states) puts all 360 particles in one state.
  - The distributed PCG reconstruction (preconditioned conjugate
    gradients) then stops: "state 2 has no particles and no previous maps
    in this directory to carry forward" (`simple_rec3D_pcg_strategy.f90`,
    line 611).
  - Phase 1 has to do no worse than this.

**Findings for later phases.**

- `scripts/plot_refine3d_bench.py` reports the pre-step's table reads and
  `ref_assign` only as means over all iterations, mixing the numbers of
  states and directions. The timing rule needs them per number of states
  and directions. (Corrected in Phase 1, which also adds that breakdown to
  the script; this entry first said the script did not report them at all.)
- The script also warns that one family/iteration/partition key occurs
  twice, and that the final noise-spectrum iteration has no strategy
  record.
- M's log warns twice that "your input doc has multiple states but
  NSTATES is not given" (`simple_builder.f90`, line 258). The warning
  comes from a program that the three-state `refine3D` starts.

**Exit criteria.**

- The build and the fast gate pass: Release and Debug both built, and
  the fast gate passed 15 of 15.
- Every case and test is recorded, with its spread where it ran twice:
  see the tables above.
- S1, S3 and M's start are kept as `scratch/keep/S1`, `scratch/keep/S3`
  and `scratch/keep/M_start`. Benchmark reports and summaries are under
  `scratch/keep/bench/phase0/<case>/`.

### Phase 1, 2026-10-16: the multi-state schedule

**Status: done, accepted by the maintainer's ruling of 2026-10-16.** The
phase first stopped at its acceptance, as the plan required: the multi-state
resolution check failed on beta-galactosidase (below). The ruling accepts
Phase 1 as implemented. The states that failed are noise classes on a
dataset without real states, the dominant states held or improved, and the
maintainer judges the schedule on msp1. The ruling also drops that
resolution rule for the rest of the run (see "Acceptance"). The account
below is the one written at the stop. The gates pass, single state is unchanged (byte-identical
simulated outputs), and the schedule is markedly faster. What fails is the
acceptance rule on the final FSC 0.143 (the Fourier shell correlation
resolution) of the weakest state, details below. Paths are relative to
`/home/elmlundho/agent_runs/ms_speed`.

**Files changed.**

- `src/main/solve/simple_solve3D_controller.f90`:
  - The new constant `NSPACE_MULTI_MAX = 2500` is the multi-state
    direction cap.
  - The function `multistate_prob_neigh_stage` returns the first stage of
    the `NSPACE` table that reaches the cap (stage 5).
  - A multi-state run uses `prob` from stage 3 up to that stage, then
    `prob_neigh` with `prob_neigh_mode=state` (in
    `set_refine3D_mode_policy`). The particle route and the class-average
    route both switch there.
  - `set_refine3D_stage_controls` caps `nspace` at 2,500 when
    `nstates > 1`. It also sets `fracsrch=95` and the internal key
    `scan_exit=yes` in stage 2 when the stage runs `shc` with several
    states. `emit_refine3D_stage_cfg` writes the key, or deletes it in
    every other stage.
  - The `solve3D_addon` early-stopping branch and its registration pass
    clear the key.
  - `NSTAGES_INDEPENDENT` keeps the value 5 with a definition of its own.
    `PROB_NEIGH_REFINE_STAGE` now concerns single state only.
- `src/main/params/simple_parameters.f90`, `simple_parameters_parse.f90`,
  `simple_parameters_phases.f90`: the internal key `scan_exit`
  (`yes|no`, default `no`) and its flag `l_scan_exit`. It is declared like
  `cohort_sampling`: the field, the parser registration, the derived flag
  and the yes/no validation, with no interface entry. The argument list of
  the generated sources picks it up, so distributed workers parse it.
- `src/main/simple_convergence.f90`: in the multi-state branch of the 3D
  convergence check, `scan_exit` adds "or the mean percentage scanned
  exceeds `fracsrch`" to the rule of every state's orientation overlap
  above 0.96 together with the scan.
  - With `fracsrch=95` the scan alone therefore decides stage 2. The old
    rule can only end the stage at the same iteration, never earlier.
  - Without the key (stage 1, every later stage, every other user) the
    rule is unchanged.
- `src/main/simple_convergence_tester.f90`: the case
  `test_multistate_scan_exit` (it runs in the unit suite
  `unit_cart_align3D`). It builds three states of two particles, all
  updated:
  - a 96% scan with a joint overlap of 0.5 ends stage 2 (key, 95);
  - the same field does not end stage 1 (no key, 99), nor the rule without
    the key at 95;
  - a 94% scan does not end stage 2;
  - an overlap of 0.98 with a 96% scan converges with or without the key.
- `src/main/strategies/parallelization/simple_refine3D_strategy.f90`: the
  strategy report (`REFINE3D_STRATEGY`) records "refine3D % search space
  scanned" after each convergence check, on both the shared-memory and the
  distributed path.
- `scripts/plot_refine3d_bench.py`:
  - a `scan%` column in the table by stage (the mean percentage scanned
    at the stage's last iteration);
  - `search_space_scanned_pct` in `bench_iterations.csv`;
  - a new table of the pre-step driver per table, number of states and
    directions: the parts' scheduling, the master's table reads, its
    `ref_assign` and the total, in seconds per iteration.
- `src/main/simple_refine3D_states_controller.f90`: `STAGE2_NSPACE`
  5,000 → 2,500. `STAGE2_NSPACE_SUB` stays 500.
- `src/main/commanders/simple/simple_commanders_solve3D.f90` (row added
  to the file table): only the comment on the default multi-state stop,
  which said it came before `prob_neigh`.
- `doc/policies/3D/solve3D_policy.md` sections 3 and 6, and
  `doc/policies/heterogeneity/refine3D_states_policy.md` section 3: the
  switch, the cap and the stage-2 exit.

A review of the diff by a separate reading pass found no defect. The key is
set or deleted in every stage, and nothing else consumes the multi-state
mode of stage 5.

**Gates.**

- Debug build with the fast gate: 15 of 15 suites passed
  (`scratch/logs/debug_build_phase1.log`), with the same seven linker notes
  as the base.
- The convergence tester's new case passes.
- The Release build is kept as `scratch/keep/build_phase1`.

**Single-state acceptance (holds).**

- The one-thread `simulated_workflow_1jxy` (both steps) and
  `simulated_workflow_6vxx` outputs are byte-identical to Phase 0
  (`scratch/logs/compare_phase1_*.log`).
  - The comparison now also ignores the manifest's `stack_digest`. That
    line hashes the stacks' absolute paths, which differ between the two
    build trees' test directories, like the run identifier and the layout
    digest already ignored.
  - Their checks equal Phase 0's (1jxy 0.9653 / 0.8863 / 4.62 Å, and with
    two partitions 0.9652 / 0.8837 / 4.71 Å; 6vxx 0.9797 / 0.8763 /
    4.62 Å).
  - Outputs are in `scratch/keep/sim/phase1/`.
- A: 928.4 s (Phase 0: 920.0 and 926.1 s). FSC 0.143 is 3.80 Å, within
  the 3.96 Å limit, and the stage pattern is unchanged.
- `solve3D_addon` passed (470 s) and `cont_refine3D_1jxy` passed (448 s).
- `flex_pca_blobs` and `solve3D_pyramid_blobs` fail exactly where they
  fail on the base, so they do no worse.
  - flex PCA leaves state 2 empty and the distributed PCG reconstruction
    stops.
  - `reconstruct3D` refuses a cropped box larger than the native one.
  - Neither reaches the changed paths: `refine3D_states_cavgs` stops
    before its `prob_neigh` march. So no test runs the `refine3D_states`
    cap.

**Timings on this machine** (Release build, one case at a time; reports
and summaries in `scratch/keep/bench/phase1/<case>/`):

| case | Phase 0 wall time (s) | Phase 1 wall time (s) | change |
| --- | --- | --- | --- |
| MS | 1,504.0, 1,501.3 | 884.9, 892.4 | −41% |
| MS-shc | 2,045.2 | 1,415.1 | −31% |
| S3 | 422.7 | 389.9 | −8% |
| A | 920.0, 926.1 | 928.4 | within the spread |

- **MS, stage 5:** `prob_neigh` over 2,500 directions takes 31.7 s per
  iteration against 82.8 s for `prob`. The tables' critical path falls
  from 1,158 s to 557 s, and their thread-seconds from 22,568 to 10,477
  (5,824 dense + 4,653 neighbourhood). The pre-step's table reads at that
  stage drop from 1.50 to 0.14 s per iteration, and `ref_assign` from 0.85
  to 0.33 s.
- **MS-shc:**
  - Stage 2 ended by the 95% rule at iteration 28 with 95.38% scanned and
    a state overlap of 0.996. That is 16 iterations against Phase 0's 20
    (its iteration limit; there it first exceeded 95% at its 17th
    iteration, and the runs diverge from stage 1).
  - Stage 5 runs `prob_neigh` as in MS.
  - The tables' thread-seconds fall from 27,189 to 16,237.
- **S3** (the class-average route shares the mode policy):
  - its stage 5 (iterations 61-72) changed from `prob` to `prob_neigh` at
    2,500 directions;
  - its stages 7 and 8 (iterations 85-96) from 5,000 to 2,500 directions;
  - its earlier stages are unchanged.

**The multi-state acceptance (fails).** Each row gives the final FSC 0.143
per state in Å, then the particles per state.

| run | state 1 | state 2 | state 3 | populations |
| --- | --- | --- | --- | --- |
| MS Phase 0, runs 1 / 2 | 4.08 / 4.08 | 9.89 / 10.88 | 5.63 / 6.04 | 3,017 / 783 / 1,713; 2,994 / 794 / 1,725 |
| MS Phase 1, run 1 | 4.18 | **14.84** | 5.26 | 2,754 / 687 / 2,072 |
| MS Phase 1, run 2 | 4.08 | 11.26 | 5.02 | 2,750 / 688 / 2,075 |
| MS-shc Phase 0 | 3.84 | 8.59 | 13.60 | 4,566 / 672 / 275 |
| MS-shc Phase 1 | 3.89 | 9.07 | **14.84** | 4,428 / 523 / 562 |
| S3 Phase 1 (no acceptance; for the record) | 23.31 | 6.66 | 15.54 | 398 / 4,120 / 995 |

- Every run ran to the end, and every state kept particles and a map.
- **The rule fails twice:** every state's FSC 0.143 must be no worse than
  1.05 times the worst state's value over the Phase 0 runs of the case.
  - MS run 1's state 2 is at 14.84 Å against a limit of 11.42 Å; run 2's
    state 2 (11.26 Å) passes.
  - MS-shc's state 3 is at 14.84 Å against a limit of 14.28 Å.
- **What the logs show** (iteration by iteration, MS):
  - The weak state 2 sits at 15-16 Å through stage 4 and most of stage 5
    in the Phase 0 runs too.
  - It reaches 9.9-10.9 Å only in the last five or six iterations of the
    `prob` stage, as its population grows from about 570 to 790.
  - With `prob_neigh` in stage 5 it grows less (from about 570 to 690) and
    stays at 14.8-15.5 Å in run 1, reaching 11.3 Å in run 2.
  - The strong states are unchanged or better: state 1 4.08-4.18 Å, state 3
    5.02-5.26 Å against 5.63-6.04 Å.
  - In MS-shc the two weak states exchange particles: state 2 went from 672
    to 523, state 3 from 275 to 562. The weakest map is 14.84 Å against
    13.60 Å.
- Beta-galactosidase has no real states, so the weak states are mostly
  noise classes, and their resolution moves by about 1 Å between identical
  Phase 0 runs.

**For the maintainer to decide (with msp1).** Whether the schedule, which
cuts MS by 41% and MS-shc by 31% on this machine, is accepted although the
weakest state of a dataset without real states ends coarser:

- in one of two MS runs, 14.84 Å against a limit of 11.42 Å;
- in MS-shc, 14.84 Å against 14.28 Å.

Nothing was tuned to meet the rule.

**Exit criteria.**

- The gates pass: yes.
- The acceptance of the schedule change holds:
  - under the rule as written, no (as above);
  - under the maintainer's ruling of 2026-10-16, yes: MS and MS-shc ran
    to the end, and every state kept particles and a map;
  - the per-state FSC 0.143 and populations are recorded against Phase 0;
  - MS-shc's stage 2 ended by the 95% rule;
  - `flex_pca_blobs` and `solve3D_pyramid_blobs` do no worse than on the
    base.
- The timings and MS-shc's stage-2 length are recorded against Phase 0:
  yes.
- Phase 1 is closed under that ruling. The Release build is kept as
  `scratch/keep/build_phase1` for the before-runs of Phase 2.

### Phase 2, 2026-10-16: `ref_assign` without rescans (stopped at two multi-threaded checks)

**Status: implemented. The assignment is identical, shown by the tester
and by one-thread runs on real data. Stopped because two acceptance checks
on multi-threaded runs fail by small margins, for the maintainer to
decide.** Paths are relative to `/home/elmlundho/agent_runs/ms_speed`.
"Before" is the Phase 1 build (`scratch/keep/build_phase1`). "After" is
this phase's Release build, kept as `scratch/keep/build_phase2`.

**Files changed.**

- `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90`,
  `ref_assign` (the assignment of particles to references from the dense
  table of the `prob` refinement). The distance preparation and each
  reference's sort of the particles are unchanged. A reference's head is
  its closest particle not yet assigned.
  - The routine keeps a count of the unassigned particles instead of
    testing them all (`any(ptcl_avail)`) for every assignment.
  - For every particle it keeps the list of references whose head it is
    (`head_first`, `ref_next`). After an assignment it advances only that
    particle's list. The earlier code ran the head update on every
    reference, which changes only those heads.
  - **State labels with several states.** A particle's state is now the
    state of its closest reference over all states, computed per particle
    in parallel; an exact tie goes to the lowest reference index (written
    in the code).
    - The earlier global greedy frontier removed the closest available
      pair of reference and head particle first, and among equal distances
      the lowest reference index (`minloc`).
    - When a particle leaves that frontier, the pair's distance is the
      particle's minimum over all references. The chosen reference is the
      lowest index whose head is at that distance. The lowest-index
      reference at the particle's minimum has the particle or an equal
      distance as its head, so it is that reference.
    - The labels are therefore identical, ties included. The plan expected
      a difference in the tie case; there is none.
  - **Within a state.** A state's references are a contiguous block (the
    state's rank times the number of directions). Each of them keeps a head
    while a particle of the state is left, so the earlier gathering of the
    state's active references always returned the whole block, in order.
    The loop now runs over the block only, with the count and the lists.
  - The random draws (the calls of the likelihood sampler, their arguments
    and their order) are unchanged.
- `src/main/strategies/search/probabilistic/simple_eul_prob_tab_tester.f90`:
  the test `test_ref_assign_against_reference` in the unit suite
  `unit_pftc_align2D3D`.
  - It keeps the earlier `ref_assign` as `reference_assign`, a copy of its
    algorithm working on the same distance array.
  - The tables have 150 particles: one state of 200 directions, three
    states of 200, and two states of one direction (each state's draw has
    a single reference).
  - The particles' states take turns, their angular distances are uniform
    over 0 to 20 degrees, and the distances are uniform over 0 to 5. The
    table's random values come from seed 20261016, and both assignments'
    draws from seed 20261017.
  - Every table holds the tie cases:
    - particle 1 is flat (every distance equal);
    - particle 2 ties its minimum between a reference of state 2 and one
      of state 3;
    - particle 3 between state 1 and state 3;
    - particle 4 has only invalid distances;
    - particles 5 and 6 are identical, which ties them in every
      reference's order.
  - The test asserts that:
    - every particle gets the same reference from both assignments;
    - the tie cases get the state of the lowest-index reference;
    - the reference assignment changes under another draw seed, so the
      agreement is not that of a draw-free assignment.
  - 24 checks pass.
- `doc/policies/3D/refine3D_policy.md`, section 8: a subsection on the
  dense-table assignment.

**Identity on real data.**

- **M on one thread.** The three-state `refine3D` case M ran with
  `nparts=4 nthr=1` and `SIMPLE_SEED=20261016` on both builds (290 and
  288 s). The orientation dump, the three state maps and their FSC files
  are byte-identical.
- **`solve3D_addon` on one thread.** The test ran on both builds through
  `simple_test_exec test=solve3D_addon nthr=1` with one OpenMP thread
  (`simple_test_exec` fixes the seed). The logs are identical apart from
  timings and run identifiers, and so are all 45 check values.
  `cohort_pose_error_excess_deg` is 4.409 in both.

**Gates and single-state acceptance.**

- Debug build with the fast gate: 15 of 15 passed, both after the first
  build and after the tester's final edit
  (`scratch/logs/debug_build_phase2.log`,
  `scratch/logs/fast_gate_phase2_final.log`).
- The one-thread `simulated_workflow_1jxy` (both steps) and
  `simulated_workflow_6vxx` are byte-identical to Phase 0
  (`scratch/logs/compare_phase2_*.log`; outputs in
  `scratch/keep/sim/phase2/`).
- A's final FSC 0.143 is 3.80 and 3.84 Å after the change, as before.
- `cont_refine3D_1jxy` passes.
- `flex_pca_blobs` and `solve3D_pyramid_blobs` fail where they fail on
  the base.

**Timings on this machine** (Release build, interleaved before/after;
reports in `scratch/keep/bench/phase2_before/<case>/` and
`scratch/keep/bench/phase2/<case>/`):

| case | before: wall time (s) | after: wall time (s) | `ref_assign`, s per iteration, before → after |
| --- | --- | --- | --- |
| A | 923.7, 932.6 | 921.8, 928.4 | one state × 1,000 directions 0.050, 0.053 → 0.045, 0.045; × 2,500: 0.121, 0.120 → 0.101, 0.099 |
| M | 63.6, 65.6, 58.8, 59.6, 60.3 | 62.5, 62.6, 62.6, 62.3, 63.6 | three states × 1,000 directions: 0.375-0.385 → 0.225-0.240 |

- On beta-galactosidase `ref_assign` was already small: about 2,500
  sampled particles in A, 5,513 in M.
- The change cuts it by 15-20% in single state and by 37% with three
  states. Neither case's wall time moves beyond its spread, which the
  table fill dominates.
- The `prob_neigh` assignment is not changed (A's `prob_neigh` rows are
  equal before and after).

**What fails.**

1. **M's state populations.** The plan's rule: they stay within the spread
   of the before-runs. With two runs each way that spread was very narrow
   (state 1: 3,388 and 3,391), and both after-runs fell outside it (3,341
   and 3,437).
   - Before running more, I decided to run three more of each way,
     interleaved, to judge on all five against all five, and to add no
     further runs.
   - Before, state 1/2/3: 3,391/816/1,306, 3,388/821/1,304, 3,409/812/1,292,
     3,456/819/1,238, 3,438/815/1,260.
   - After: 3,341/820/1,352, 3,437/814/1,262, 3,377/824/1,312,
     3,389/819/1,305, 3,395/815/1,303.
   - Two after-runs fall outside the before range (state 1 3,388-3,456,
     state 2 812-821, state 3 1,238-1,306): run 1 in states 1 and 3, run 3
     in all three states.
   - The last iteration's mean score is 0.410 in every run (0.411 in one
     before-run).
   - The one-thread M above shows the assignment is identical. The
     multi-threaded table fill is not reproducible, and the one-thread
     run's own state-1 population is 3,425.
2. **`solve3D_addon` in its CTest registration (8 threads).**
   `cohort_pose_error_excess_deg` must stay at or below 5.0 degrees.
   - After: 5.343 and 5.195 (fail).
   - Before: 2.614 (Phase 1) and 2.957 (pass).
   - The one-thread runs give 4.409 on both builds, with identical logs.
     The test's floor is therefore sensitive to the multi-threaded
     nondeterminism of the run.
   - A spread of three more eight-thread runs per build is recorded below.

Neither check can be passed without a tolerance that neither the plan nor
a ruling gives, or without repeating runs until they pass. So the phase
stops. **The maintainer decides:**

- whether the identity evidence (tester and one-thread runs) settles the
  identical change, and how M's population check and `solve3D_addon`'s
  5-degree floor should treat multi-threaded runs;
- whether Phase 2 is accepted.

**`solve3D_addon` spread** (eight threads, `scratch/run_addon_spread.sh`
and the CTest runs; logs `scratch/logs/addon8t_*.log` and
`scratch/logs/ctest_solve3D_addon_build_phase{1,2}.log`). The check is
`cohort_pose_error_excess_deg`, with a 5.0-degree limit.

- Before: 2.614, 2.957, 2.716, 4.281, 2.857 (all pass).
- After: 5.343, 5.195, 5.970, 1.507, 2.177 (three of five fail).

The two sets overlap, and five runs each way cannot separate them: a Fisher
exact test on three failures of five against none of five gives p ≈ 0.17.

The two builds differ only in `ref_assign` and the tester. `ref_assign`
returns the same assignment for the same table: the tester shows it, and so
do the one-thread runs of M and of `solve3D_addon`. The run-to-run variation
of the eight-thread runs comes from the table fill, whose dynamic
scheduling decides which thread's random draws a particle meets, and that
fill is the same code in both builds. The test's 5-degree floor lies inside
the spread of its eight-thread runs.

**Exit criteria.**

- The gates pass: yes.
- The assignments are identical:
  - yes in the tester (every particle, ties included) and on one thread
    with real data (M and `solve3D_addon`);
  - the multi-threaded acceptance checks (M's populations, the eight-thread
    `solve3D_addon`) fail as described.
- The timings are recorded against Phase 1: yes.
