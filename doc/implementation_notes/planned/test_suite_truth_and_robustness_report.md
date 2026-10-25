# A test suite that knows the truth: report of the runs `test_truth` and `test_truth_sync`

Report, 2026-10-19, on the plan `test_suite_truth_and_robustness.md` (beside this file), whose
"Progress" section holds the detail of every phase. The run `test_truth` carried out Phases 0 to 6
on 25d256886 (16-18 October 2026) and was stopped in Phase 6 so that its work could be merged onto
master 43de9b340; `test_truth_sync` carried out Phases 7 to 9 on the commit that holds the merge
(a67d487d9), unattended on the 24-core Oracle Linux test machine. Nothing is committed: the whole
change is one working-tree diff. The diffs of this run's phases are in
`~/agent_runs/test_truth_sync/review/phase_7.diff`, `phase_8.diff` and `phase_9.diff`; the earlier
run's work arrived as part of the base commit.

## What the suite verifies now

The suite keeps four kinds of test apart (test policy, section 2.3):

- **identity**: a refactoring changed nothing, judged by one-thread runs compared byte for byte
  (`simple_test_exec test=identity`, over the outputs each registry row declares);
- **correctness**: one run recovers a known truth, held to floors that sit far from both the
  observed scatter and a negative control (a shuffled label, the other molecule, copies of one map);
- **robustness**: how reliably, measured by profiles over 20 seed offsets and thread settings
  (`test=profile`), never a single-run gate;
- **performance**: timings, never in a pass/fail test.

Every high-level test now writes its numbers to a table (`metrics.tsv`) and has a card in the test
policy (section 6.5): its claim, the path it asserts ran, its truth metrics, each floor with the
profile it came from, its negative control, its determinism and its cost. Every floor comes from a
profile of the merged code by the floor rule (at most the worst run minus one robust spread, at most
the median minus three, at least the control plus three), or the metric is reported. The suite runs
as one command, `simple_test_exec test=nightly`, which packs the entries on the machine's cores,
memory and disk; `ctest -L highlevel -j 24` packs the same entries by their CPU reservations.

What the suite could not say before and says now:

- **Multi-state results are judged against known state labels**, on simulated data only.
  `simulate_particles` makes mixtures of several volumes with given populations; a fixture builds
  compositional beta-galactosidase states (the tetramer, the tetramer less one chain, a dimer),
  blob phantoms and their truth; shared metrics map estimated to true states once (labels, maps and
  poses under one mapping) and score labels by the adjusted Rand index and an accuracy judged
  against its own null. `multistate_assignment` tests `refine3D nstates=3`; `multistate_workflows`
  tests the class-average route, `solve3D` from particles, `refine3D_states` and the route used on
  real data (single-state `solve3D` then `refine3D_states`).
- **The de novo multi-state routes are a robustness measurement** (ruling of 2026-10-18): they
  resolve every true state in 13 (particles) and 9 (class averages) of 20 seeds, a limitation of the
  outright state assignment from random labels. The profile holds those counts to binomial floors
  (6 and 2 of 20); the nightly keeps their path, checks and pose ceilings, which held in every run.
- **Fragile floors are gone.** The add-on's pose gate is relative to its own base; floors set from
  one run were replaced (the mini-stream picking floor, the simulated workflow's correlation); the
  SINGLE workflows' pose metric is symmetry-aware.

## What changed, phase by phase

- **Phases 0 to 6 (`test_truth`, merged into the base):** the tooling (the registry, the test
  dispatcher with its failure and deadline contract, `test=nightly`, `test=profile`,
  `test=identity`, the seed offset), the first profile of every entry, `simulate_particles` with
  several states, the fixture and the metrics, the multi-state tests, the repairs (seven production
  defects fixed: among them `flex_pca` parking a real state in state 0, shared-memory `refine3D`
  not seeding its states, an empty or one-sided state stopping the distributed PCG, `autoscale`
  planning a crop above the box), the review list of Phases 0 to 3, packing by cores, memory and
  disk. Its Progress entries list each with its evidence.
- **Phase 7 (`phase_7.diff`): the merged code, and what still holds.** The gates pass on the merged
  code. The speed week (04e520e78) is byte-identical to 25d256886 at one thread on x86, on the
  simulated workflow, the three-state, correlation and PCG cases; its part files, whose layout it
  changed on purpose, compared through lossless decoders; its one intended result change (keeping
  back-filled projection labels) is reached by no run (details in `identical_speed_run_report.md`).
  Entries without `solve3D` or `refine3D` (and `single_workflow_fcc`) are byte-identical between the
  old run's final build and the merged build, so their profiles stand; `test=identity` gained the
  declared outputs to show it. Nine entries went on the list to profile again.
- **Phase 8 (`phase_8.diff`): profiles and floors on the merged code.** The nine entries profiled
  (20 seeds; one thread as well for the two the speed week checks at one thread, for
  `multistate_assignment`, and for the new route), every floor rechecked and its origin written. One
  floor moved by the rule (the `prob` adjusted Rand index of `multistate_assignment`, 0.90 to 0.89,
  from the one-thread profile); two gates became reports for want of an admissible bound (the
  `simulated_workflow_1jxy` FSC ceiling, which failed 2 of 20, and the class-average state route's
  map-match check in `flex_pca_blobs`, which failed 1 of 20); the two rulings that had reached no
  phase were applied. Tools: `test=profile keep_outputs=no` (outputs deleted as executions end, so
  long profiles fit on the disk inside the profile commander) and `success` rows with their profile
  verdict.
- **Phase 9 (`phase_9.diff`): timing, documents and close.** The wall time of the suite and the
  speed week's timing on x86 (below); the test cards in the policy; the stale text listed in
  `testing_and_benchmarking_plan.md` section 3 (unit-suite tables, the seed-offset wording, the
  four-thread fast tier that nothing runs, `test=nightly` as the tier of tests too expensive for the
  build, the CMake comments); the real-data layer's design in the policy (section 6.6); this report.

## Evidence against the baseline

- **Every floor of the suite comes from a 20-seed profile of the merged code** (Phase 8, kept under
  `~/agent_runs/test_truth_sync/scratch/keep/profiles/phase8/`), or from Phases 1 and 5 for the
  entries whose outputs Phase 7 showed byte-identical on the merged code. The numbers are in each
  test's source beside its floor and in the cards.
- **The whole suite on the final code** (the otherwise idle 24-core machine, Phase 9):

  | run | Phase 5 (25d256886 and that phase's changes) | Phase 9 (final code) |
  | --- | --- | --- |
  | `test=nightly ncunits=24` | 3 724 s, 2 entries failed | 2 633 s, 1 failed |
  | `ctest -L highlevel -j 24` | 2 721 s, 2 failed | 2 775 s, 1 failed |

  The one failure is `solve3D_pyramid_blobs`, known since Phase 1 (below, under "Open for the maintainer"). The nightly is
  29 % faster despite the fourth multi-state route: master's smaller simulated workflows no longer
  set its wall time.
- **The speed week on x86.** Byte-identical to its base at one thread (Phase 7), and faster where it
  aimed (Phase 9; twice each way, interleaved): the three-state `solve3D` (MS) −11.7 %, the dense
  three-state `refine3D` (M) −15.2 %, distributed single-state `solve3D` (A) −3.7 %, `solve2D` (E)
  −4.0 %, their shared-memory forms −0.5 % and −0.8 %. The table fills and reads are faster; the
  nonuniform filter is slower on this machine (about 20 s per run of A), see
  `identical_speed_run_report.md`, "Timing on x86".

## Decisions taken under delegation

Each is recorded where it applies; the Progress entries of the plan give the reason.

- **Phase 7:** 13795c1d9 (the axis-search move, contained in 04e520e78) was not built separately for
  the comparison, since no case searches an axis; three exclusions were declared on x86 beyond the
  Mac's (PDF plots' dates and ids, the part files through decoders, the PCG diagnostics' timing
  tables); entries working in a directory named by process id are compared through that directory.
- **Phase 8:** the profiles ran inside `test=profile` with outputs deleted as they end (a tool
  change proposed in the Phase 5 list for the maintainer) rather than in shell loops; a success
  counts only in an execution that passed; the consensus route resolves every state in all 40
  runs, so its resolution is checked per run and its labels gated by the floor rule (as the ruling
  asks for what holds run to run), while its per-state recall is reported; the class-average route is the one the ruling of 2026-10-17
  defined (imported images), so its success count measures `solve3D_cavgs` alone; de novo routes
  were profiled at eight threads only (the ruling names one thread for the new route alone);
  `multistate_workflows`' deadline rose to 7 200 s for its four routes.
- **Phase 9:** the timing inputs are copies of `prob_speed`'s and `ms_speed`'s kept projects with
  their paths remapped by `update_project` (the originals point into those runs' deleted case
  directories); E-shm was timed beside A, A-shm, E, MS and M because the testing plan names it.

## Deviations from the plan

- The speed week's identity on x86 covers the cases the continuation's plan names (the 1jxy
  simulated workflow, the three-state, correlation and PCG cases); the 6vxx workflow, the filter's
  standalone programs, the add-on and the continuous-refinement cases, which the Mac also compared,
  were not compared at one thread on x86.
- The class-average multi-state route does not run `solve3D cavg_ini_ext=yes` after
  `solve3D_cavgs`, as the ruling of 2026-10-18 describes the route: its images are imported class
  averages without particles (ruling of 2026-10-17); the plan's Phase 8 entry proposes, for the
  maintainer, a second class-average route on particles.
- The hard (conformational) profiles were not repeated on the merged code; they inform no floor.

## Open for the maintainer

The plan's Phase 5 and Phase 8 entries list fourteen items, each with a diagnosis and a proposal.
The ones that decide the suite's next step:

1. **De novo multi-state assignment** loses or merges states in about half the runs (particles 13,
   class averages 9 of 20 resolved); the planned annealed state draw is the fix, and the per-seed
   metrics of both routes under `~/agent_runs/test_truth_sync/scratch/keep/profiles/phase8/`
   (`msw_t8_shc`, `msw_t8_cavgs`, `msw_t8_consensus`, `msw_t1_consensus`, built from this run's
   `build_p8`) are its baseline, run for run on the same seed offsets.
2. **`solve3D_pyramid_blobs` fails 20 of 20**, mostly on production's own stops after its state step merges
   blob states that differ in 6 % of the mass: move it to the compositional states.
3. **`flex_pca_blobs` part 3** (`refine3D_states_cavgs` on pose-less sub-averages) stays at chance:
   carry the members' poses into the sub-averages.
4. **The simulated workflows have no map gate** (2 of 20 and 1-2 of 20 missed maps): hold "the map
   is found" as a success count, or raise the signal.
5. **`simple_cmdline` exits 0 on missing keys**, so a CTest entry missing a key passes.
6. **The nonuniform filter is slower on x86 since the speed week** (A 52.7 against 73.4 s of filter
   per run): time `nu_filt3D` on this machine with the candidates filtered again and with the scratch
   files before that choice stands.
7. **The Mac has not profiled the changed cases.** Commands, on a Release build with tests, from its
   build directory with `SIMPLE_PATH` set to it:

   ```text
   simple_test_exec test=profile target=solve3D_addon,cont_refine3D_1jxy nthr=1,8 nseeds=20 keep_outputs=no dir=mac_profile_a
   simple_test_exec test=profile target=simulated_workflow_6vxx,simulated_workflow_1jxy,single_workflow_wurtzite nseeds=20 keep_outputs=no dir=mac_profile_b
   simple_test_exec test=profile target=multistate_assignment nthr=1,8 nseeds=20 keep_outputs=no dir=mac_profile_c
   simple_test_exec test=profile target=multistate_workflows suite=consensus nthr=1,8 nseeds=20 keep_outputs=no dir=mac_profile_d
   simple_test_exec test=profile target=multistate_workflows suite=shc nseeds=20 keep_outputs=no dir=mac_profile_e
   ```

   Each writes `profile_<entry>_t<n>.tsv` with each metric's median, robust spread and worst run;
   the floor rule against the floors in the sources tells whether they hold on the Mac too. The Mac's values of 2026-10-17 (the add-on's absolute
   excess 6.00° against 5°, the 6vxx correlation 0.780 against 0.80) were measured on cases that
   have changed since.
8. **The real-data layer** (ruling of 2026-10-18, late; policy section 6.6): public datasets with
   previously solved maps, as Zenodo subsets fetched by name with checksums, each judged against
   its earlier map by FSC and correlation with floors from a profile and another dataset's map as
   the control, never an input. Open: the datasets and their licences, the design of each case
   (starting point, downsampling, number of images, budget), how the fetch is done, and whether
   `SIMPLE_data_testing`'s datasets are the same data. Nothing was chosen, downloaded or uploaded.
