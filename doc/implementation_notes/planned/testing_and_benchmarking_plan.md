# Testing and benchmarking: what is left, and in what order

Joint plan, 2026-10-18. **Partly done.** Of the order in section 8, the `speed2d_mac` run (step 1) and `test_truth_sync` Phases 7 to 9 (step 2, commit `34f314fe5`) have landed; `scripts/plot_refine3d_bench.py` exists; the matching-band test is dropped (ruling 1). The rest is open. Rulings of 2026-10-18 in section 2. It gathers the
open testing and benchmarking items of the planned notes into one order of
work. The notes it draws on are listed at the end; each item names where its
detail lives.

Terms used throughout:

- An **identity test** runs a case at one thread before and after a change and
  compares the outputs byte for byte. It is the only test of a change that
  claims to leave results unchanged.
- A **correctness test** recovers a known truth from simulated data in one run,
  against a floor.
- A **robustness profile** runs a case over 20 seed offsets (and thread
  settings) and reports the spread. It is never a single-run gate.
- The **floor rule** sets a floor from a profile, far from both the observed
  scatter and a negative control (a run that must fail, such as shuffled
  labels or the wrong molecule).
- `BENCH/` is the directory of appended timing reports that `refine3D`,
  `solve3D` and their parts write. `scripts/plot_refine3d_bench.py` summarises
  one run's reports.

## 1. Where things stand

**The test suite.** The test-suite overhaul (the run `test_truth`, Phases 0 to
5 on the Dell) built:

- a test registry;
- `simple_test_exec test=nightly`, `test=profile` and `test=identity`;
- seed offsets;
- a `metrics.tsv` file of gates and reports per high-level test;
- the floor rule;
- multi-state tests with known labels.

It was merged onto master 43de9b340 on 2026-10-18 and is not yet committed. The
continuation run `test_truth_sync` (the test-suite plan's Phases 7 to 9) is
specified and waits for that commit. Its work:

- the speed week's identity check on x86;
- profiles and floors again on the merged code;
- de novo multi-state measured as a success rate;
- the route the maintainer uses on real data;
- profiles at one thread;
- the suite's timing and the report.

**Benchmarks.**

- The speed week is on master. Its runs were `prob_speed` and `ms_speed` on the
  Dell, and `speed_mac` on the Mac.
- Timing reports are appended to `BENCH/` and summarised one run at a time.
- The run `speed2d_mac` is planned for the Mac. It brings the 2D timing records
  and the 2D speed items, and a new testing rule: every changed parallel region
  is stress-tested at the machine's full thread count, 20 repetitions, because
  identity at one thread cannot show thread safety. A race in commit 43de9b340
  passed identity.

**Real data.** Real data is handled by three separate mechanisms today:

- **`SIMPLE_data_testing`:** GitHub Actions on datasets trimmed to 50
  micrographs, with an FSC pass/fail and a web page (test policy section 6.4).
- **The hardware resource planner:** Python under `scripts/`, with its own
  manifest, fixtures and a NumPy re-implementation of the simulated workflows.
  It runs on the Xeon desktop and the CSB cluster
  (`hardware_resource_planning_completed.md`).
- **The native dataset layer:** planned for the Dell (test policy section 6.6),
  not built.

The reference benchmarks are still run by hand:

- PfCRT and msp1, with reports on the Mac;
- the beta-galactosidase cases S, A, A-shm, E, E-shm, MS and M, with inputs on
  the Dell.

## 2. Rulings of 2026-10-18

1. **No matching-band test.** The test proposed in
   `matching_band_contract_and_test.md` section 4 is not pursued. What the
   latest `solve3D` changes leave testable is the stage rules themselves
   (section 3, "Stage-shift rules").
2. **The stage shift is judged on real data.** The maintainer tests the
   FSC=0.5 stage shift (cd13708eb) and the nonuniform filtering from stage 5
   (8c8d8b427) on real datasets; no paired truth profile.
3. **Case A** is `prob_speed`'s: `solve3D` from S in shared memory (660 s,
   then 479 s, 4.03 Å). A, A-shm, E and E-shm come from `prob_speed`, MS and M
   from `ms_speed`, each judged against its own run's numbers. `ms_speed`'s
   variant from S1 (`nparts=4 nthr=6`, 920 s, 3.84 Å) is dropped.
4. **One real-data mechanism per purpose** (refined by ruling 9):
   - the native dataset layer on the Dell for the full-size benchmarks and
     their profiles, with public cases from Zenodo in `test=nightly`;
   - `SIMPLE_data_testing` kept for the trimmed public datasets on GitHub;
   - the hardware planner kept for resource prediction only. Its workflow
     re-implementation and its testing role are retired in favour of
     `simple_test_exec`.
5. **`test=nightly` is the tier of tests too expensive for the build.** Today
   that is the simulated tests. It is not a nightly of record with archives
   and notifications, and it does not take in the library, platform or fast
   tiers. Real datasets may join it, run from within SIMPLE like any other
   entry, once it is settled how experimental data reaches the machine
   (section 5).
6. **The seed key is approved.** `rnd_seed_parameter.md` goes ahead. Its file
   table gains the seed-offset code, the dispatcher's child seeding and the
   hardware planner's use of `SIMPLE_SEED`.
7. **Order of the method projects:** the annealed state draw first, then the
   placement of pose-less sub-averages in `refine3D_states_cavgs`.
8. **The stage after a class-average start is never skipped, but may end
   early.** In `cavg_ini` and `cavg_ini_ext` runs the particles' entry stage is
   excluded from the FSC=0.5 stage shift (as the code does today) and keeps the
   mid-stage exit (as the code does today). The eligibility goes into one
   shared function with unit tests (section 3).
9. **Real test data on Zenodo, judged against previously solved maps.** The
   real-data cases of `test=nightly` use public datasets whose maps have been
   solved before. Subsets go up on Zenodo, and the earlier map serves as the
   reference ("ground truth" for the map). Each case needs careful design:
   starting point, downsampling and number of images, so that it tests what
   it is meant to test in the time it is given (section 5).

## 3. Now: the suite on the merged code (`test_truth_sync`, Dell)

Specified in `test_suite_truth_and_robustness.md`, Phases 7 to 9. Phase 9
times A, A-shm, E and E-shm as `prob_speed` defined them and MS and M as
`ms_speed` did. It also closes these documentation gaps:

- The test cards go into the test policy (section 6.5 is still
  `CARDS_PLACEHOLDER`; drafts are on the Dell in
  `~/agent_runs/test_truth/scratch/p6_cards.md`).
- Stale text in the test policy:
  - the unit-suite tables of sections 1.1 and 1.2 miss about ten sub-suites;
  - section 1.4 says a seed offset "adds n" (it multiplies by a large odd
    stride) and that seeds come through `parameters` (they do not yet);
  - section 4.5 says the nightly repeats the fast tier at four threads
    (nothing does, and by the ruling it will not);
  - section 6 presents `test=nightly` as the nightly suite; by the ruling it is
    the tier of simulated tests too expensive for the build.
- `production/CMakeLists.txt` comments: "30 s", "RUN_SERIAL", the heartbeat's
  seven children, and `multistate_workflows` running one route.
- The test plan's header.

**Stage-shift rules.** The latest `solve3D` changes are tested on real data
(ruling 2). What is cheap to test in the fast tier is the rules:

- **Eligibility, written twice.** The stage shift and the mid-stage exit are
  meant to share their eligibility. The shift
  (`simple_commanders_solve3D.f90`) excludes:
  - stages 1 and 2 (inside `solve3D_stage_behind_fsc05`);
  - the entry stage of the loop;
  - the last stage;
  - the axis-search stage when an axis is searched.

  The exit (`lpexit`, emitted in `simple_solve3D_controller.f90`) excludes the
  same except the entry stage. The two agree for a run from stage 1. They
  differ for runs that enter later (`cavg_ini`, `cavg_ini_ext`, input volumes):
  there the exit may end the entry stage early, while the shift never skips it.
  Ruling 8 keeps that asymmetry for the class-average starts: the entry stage
  is never skipped, and may end early. One pure function used by both, with
  the two answers (may skip, may exit early), is unit-tested over the entry
  routes: from stage 1, `cavg_ini`, `cavg_ini_ext`, input volumes, the add-on,
  state continuation, with and without an axis search. Where comments and the
  policy say the entry stage of a `cavg_ini` run "runs its full budget", they
  gain "unless the mid-stage exit ends it".
- **The thresholds.**
  - The exit compares the FSC=0.5 crossing with the next rung, never with the
    band in use. The first version did the latter and was reverted after PfCRT
    converged in 2 of 6 restarts.
  - The advance moves one rung, is bounded by the cap and only sharpens. There
    is no advance from the last rung.

  These are a few assertions once the rung arithmetic takes the crossing as an
  argument instead of reading the project.
- **The ladder is strictly decreasing.** The advance relies on it, and so does
  the reasoning that a skipped stage is always followed by a finer one.
  `simple_lpstages_tester` has the cases (from an FRC, a flat FRC with its
  linear fallback, class averages against particles) but does not assert
  monotonicity.

## 4. Next: runner and benchmark tooling (one run, after `test_truth_sync`)

The runner:

- **Critical path.** Order the `test=nightly` queue by expected core-seconds
  (measured wall time times reservation), not by reservation alone. Today the
  longest one-thread workflow starts last and alone sets the wall time: 3,724
  s, against a target of under 45 minutes (2,721 s through CTest).
- **Disk budget.** Make the dispatcher's disk budget count only running
  executions when their outputs are deleted, so a profile no longer runs one
  execution at a time.
- **Missing keys.** `simple_cmdline` exits with status 0 when a required key is
  missing, so a CTest entry that lacks a key passes. It must exit non-zero.
  This is a production fix.
- **Stage-shift rules** (ruling 8, section 3): the shared eligibility
  function and its unit tests, the rung arithmetic taking the crossing as an
  argument, and the monotonic-ladder assertion in `simple_lpstages_tester`.
- **Seeds** (ruling 6):
  - the `rnd_seed` key;
  - no fixed seeds in production code: `flex_pca`'s class expansion seeds
    itself, so part 3 of `flex_pca_blobs` ignores the seed offset;
  - the five testers that call `seed_rnd_fixed` directly.
- **Mac limits** stated, not fixed: no CPU or memory samples, since the
  dispatcher uses `/proc`.

Benchmark reporting:

- **Run identity.** Each `BENCH/` record names its run, so a harness no longer
  has to clear the directory between runs.
- **Before/after summaries.** `plot_refine3d_bench.py` compares two runs (today
  it summarises one) and shows the 2D critical path once `speed2d_mac` adds the
  2D records. Fix its two known warnings: a duplicate `VOLASSEMBLE` record at
  iteration 1, and a missing strategy record at the closing sigma iteration.
- **The matching band per stage.** Log `res05`, the planned rung and the band
  used per stage, so the band's trajectory can be compared run to run
  (`matching_band_contract_and_test.md` section 3, "not folded in").
- **Memory.** Measure the master's peak memory in a large distributed
  `refine3D`: about 1 GB per handed-on project copy at 750,000 particles,
  unmeasured. Also the nonuniform filter's setup memory in volume assembly under
  ML regularisation (+217 MB at box 300, measured only standalone).
- **The prob_tab objective-evaluation line** in the `refine3D` reports changed
  meaning in the speed week. Mark it, so old and new reports are not compared.

## 5. Then: real data in the suite

Design in test policy section 6.6, changed by rulings 5 and 9: real data
reaches the suite as cases of `test=nightly`, on public datasets kept on
Zenodo, each judged against a previously solved map.

**The data.**

- **Public datasets with a solved map.** EMPIAR datasets whose maps are
  deposited in EMDB, for example beta-galactosidase, apoferritin, TRPV1 or the
  Plasmodium 80S ribosome. Each must have a licence that allows
  redistribution of a subset; check every one.
- **On Zenodo,** as versioned records with a DOI. A record holds:
  - the subset;
  - the reference map;
  - the case's command lines;
  - a checksum per file.

  Each record must stay within Zenodo's per-record size limit (50 GB at the
  time of writing).
- **Fetched by dataset name** into a local cache, checked against the
  checksums. A missing or corrupt dataset skips its case with a visible status,
  never a pass. How the fetch is done is part of the design: a SIMPLE program,
  or a one-time setup step per machine.
- **The full-size benchmarks stay where they are.** PfCRT, msp1 and the
  beta-galactosidase speed cases, which are not public, stay as local copies
  behind a per-machine manifest, for profiles and timings.
- **`SIMPLE_data_testing`.** Agree with its owner whether its trimmed datasets
  and the Zenodo records are the same data, so the two do not drift apart.

**What a reference map can and cannot judge.**

- **The map.** After docking in both hands, compare the run's map with the
  reference: the FSC between the two and their correlation. Use a band that
  stops short of the reference's own resolution, with floors by the floor rule
  from a 20-seed profile. The negative control is another dataset's reference.
- **Never an input.** The reference is never the start or a reference of the
  run, or the comparison measures reference bias.
- **Not independent truth.** The reference came from the same images, often
  with more of them and other software. It judges the map, not the poses or
  the state labels, which stay unknown. A multi-state dataset is compared state
  by state only where per-state maps were deposited, and its labels are never
  judged.
- **Robustness and performance, as before:**
  - does the workflow complete;
  - does it hold its resolution limit;
  - how long does it take;
  - how much do the results scatter across seeds.

  cFAR (the conical FSC area ratio, a measure of anisotropy) is a stable
  metric beside the FSC=0.143 resolution. FSC=0.5 is never a pass/fail gate
  on a dataset whose FSC runs along 0.5. On PfCRT a 0.03 change moves the
  crossing by 1.5 Å (`matching_band_contract_and_test.md`).

**Designing a case.** This is the work. Each case fixes, and writes in its
test card:

- **the claim:** which workflow or stage it tests, and what a regression would
  look like;
- **the starting point:**
  - movies, micrographs, extracted particles with CTF, or a 2D solution;
  - expensive steps that another case covers are skipped;
  - a start such as `cavg_ini_ext` follows the `solve3D` test protocol
    (`nsample` 2500 and the common `solve3D_cavgs` start);
- **the downsampling:** pixel size and box chosen so the target resolution
  sits well inside Nyquist and the run fits its budget;
- **the number of images:** enough to reach the target resolution with margin,
  as a fixed subset (named micrographs or particle indices) in the record;
- **the budget:** wall time and CPU reservation within the nightly's target. A
  case that needs longer belongs to the profiles or the timings, not the
  nightly;
- **the reference comparison:** the reference's resolution, symmetry and hand,
  and the comparison band;
- **its profile:** 20 seed offsets before its floors are set; one thread where
  a refactoring run checks it at one thread.

**Components.**

- the fetch and cache, with checksums;
- a SIMPLE test program that runs each case through the same dispatcher;
- `metrics.tsv` rows:
  - map-to-reference FSC and correlation as floors;
  - FSC=0.143 resolution and cFAR as ceilings and floors;
  - wall time and per-stage `BENCH` times as reports;
- seed profiles through `test=profile`.

**The speed week's real-data checks** stay with the local full-size data:

- PfCRT on the merged code with the split timers. The week plan wanted one
  PfCRT run per block; one profile replaces them.
- The msp1 timing (10,415 s at baseline, estimated 4,000 to 4,700 s,
  unmeasured).
- The multi-state schedule judged on msp1 by populations and per-state maps
  across seeds.

**Later.**

- The 66 real platinum maps under `/data/NanoX/zenodo`: plain detection must
  reproduce the stored coordinates.
- FePt on the Mac.
- The `_unfil` half-map pairs of the nonuniform evidence-mask diagnostic.
- The "dataset with one stack per micrograph" that the week plan asks for is
  defined nowhere: define it or drop it.

## 6. Test gaps to close (small items, a day run or by hand)

Simulated workflows and streams:

- **`simulated_workflow`:**
  - `simulate_movie` saves particle positions, so picking gets precision and
    recall;
  - per-movie defocus and drift are compared with the truth files the simulator
    already writes;
  - the map correlation is gated again at master's sizes if the profile allows
    (in `test_truth_sync`).
- **`mini_stream`:** exercise motion correction, and use its eight CTest
  threads or lower its reservation.

Negative controls:

- `pcg_recon` against a wrong phantom;
- `cont_refine3D_1jxy` against shuffled poses.

SINGLE:

- **`single_workflow`:**
  - an absolute pose ceiling from the existing profile (medians 4.9° fcc, 1.5°
    wurtzite; only the +0.5° per stage relative ceiling is gated);
  - per-stage counts;
  - Ruben to confirm `smpd=0.358`;
  - confirm `single_workflow_fcc` passes on the merged code, now that the clamp
    is reverted.
- **`single_atoms_stats`:**
  - a noisy variant (today its spread is zero, so it is a regression check, not
    a calibrated gate);
  - one-to-one atom matching;
  - the nearest-neighbour peak for platinum.
- **The C-alpha finder:**
  - one-to-one recovery;
  - N- and C-terminal exclusion;
  - candidate count and score;
  - a negative control;
  - rotation.

  Today it has two assertions (`single_area_tests_handover.md`).
- **The SINGLE atomic-model comparison** end to end.

Pose and picking paths without a test:

- **`solve3D` with external volumes** copies 2D shifts onto random orientations,
  the pattern that cost the class-average route most of its poses. Test it, then
  fix it with the pose-less rule. This is a production defect candidate.
- **The second picker (`picker=new`):** register it or defer it (a ruling on the
  process budget).

Test helpers of the next-release refactor:

- The staged refactor's test helpers (`staged_domain_driven_design_refactor.md`)
  are mostly met by `test=identity` and `test=profile`. Its blessed baselines
  conflict with the policy's "no blessed baselines"; settle that at its M0
  ruling.

**Stream (Joe and Ruben):**

- the sieve's full rejection reason and one class per hard gate;
- optics: a one-cluster control, the threshold, restart;
- picking references;
- pick and extract (the cap equals the asserted count; positions and a
  negative control are missing);
- the heartbeat's failure paths;
- a recorded passing run of the two chained library tests (`lib_stream`, which
  the nightly does not run).

## 7. Method work the suite now makes possible

Each is a separate project that changes results. Each is judged by a paired
20-seed truth profile against the baselines `test_truth_sync` keeps.

- **The annealed state draw** for de novo multi-state (ruled 2026-10-18, after
  the suite). Today every true state is resolved in 11 of 20 runs from particles
  and 9 of 20 from class averages.
- **`refine3D_states_cavgs` placing pose-less sub-averages.** Carry the members'
  poses into the class expansion, copy them to both halves, and refine them
  locally. Sub-average state accuracy is at chance today (0.42 against a
  shuffled 0.37). (`solve3D_pyramid_blobs`, which waited on it, was removed on
  2026-10-20.)
- **The objective in multi-state class-average stages.** `solve3D_cavgs` forces
  `euclid`. `cc` resolved 2 of 3 runs against 1 of 3, with poses four times
  worse.
- **`refine3D` started from given maps under `euclid`** merges states in its
  first iteration (no noise model yet). The test runs `cc` meanwhile.
- **The single-state class-average start at `prob_athres` 10** (ruled): map
  correlation median 0.91 against 0.98 at 90, while FSC and rotation error are
  slightly better. Check the map's scale and mask before the default settles.
- **`flex_pca`'s population-floor path** gives unsupported particles random
  labels instead of state 0.
- **The FSC=0.5 promotion** takes the best state's crossing (its optimistic
  tail); see `matching_band_contract_and_test.md` section 3.
- **The FSC=0.5 stage shift and mid-stage exit** use the best state's crossing
  (`minval` over states), so in a multi-state run the best state can end a stage
  for all of them.
- **The speed items that change results** (`speed_architecture_proposals.md`):
  - the master draws the stage 1-2 sample;
  - 2D partitions balanced by active particles;
  - the `maxits`/`minits`, `fillin` and `prob_snhc` renames, which also need
    the stream tests;
  - single-precision weights;
  - a cap on refined references across states.

Validation tests added when their features land. These are not now, but are
budgeted here so they are not forgotten:

- the `refine_tomo` phantom tilt series with known angle errors;
- gate 2 of the PCG atomic prior (a synthetic tilt series, true and damaged
  models);
- the B-factor `single_workflow` case, and identity with `fit_bfactors=no`;
- a two-species `autorefine3D_nano` run;
- the C-alpha evidence-mask fixture.

One-off validation runs stay with their notes: species phase 9 on the 66 maps,
FePt with two species, the nonuniform evidence-mask diagnostic.

## 8. Order

1. **Tonight, Mac:** `speed2d_mac` as planned (2D timing records, identical 2D
   items, thread stress tests).
2. **Dell, as soon as the merged code is committed:** `test_truth_sync`, Phases
   7 to 9.
3. **Dell, one run:** runner and benchmark tooling (section 4), with the seed
   key and the stage-shift rules.
4. **The real-data cases (section 5).**
   - First, a design note: datasets, starting points, downsampling, numbers of
     images, budgets, references.
   - Then the Zenodo records.
   - Then one or two Dell runs: the fetch, the test program, the profiles.
   - The local full-size benchmarks (PfCRT, msp1, the beta-galactosidase speed
     cases) go alongside.
5. **Alongside, by hand or in day runs:** the test gaps of section 6.
6. **Then the method projects of section 7,** each its own run (ruling 7): the
   annealed state draw first, then the placement of pose-less sub-averages.

## Sources

Under `doc/implementation_notes/planned/`:

- `test_suite_truth_and_robustness.md`
- `single_state_speed_week_unfinished.md`
- `multistate_speed_run_unfinished.md`
- `speed_architecture_proposals.md` (the maintainer's tree holds the newer copy)
- `hardware_resource_planning_unfinished.md`
- `species_discovery_status.md`
- `atom_bfactor_model.md`
- `single_pcg_integration_unfinished.md`
- `nu_evidence_mask_cross_validated_threshold.md`
- `next_release/calpha_nu_evidence_mask.md`

Under `doc/implementation_notes/completed/`:

- `matching_band_contract_and_test.md`
- `single_state_speed_week_completed.md` and its report
- `multistate_speed_run_completed.md`
- `identical_speed_run.md`, its report and review
- `hardware_resource_planning_completed.md`
- `single_pcg_integration_completed.md`

Under `doc/refactoring_notes/completed/`:

- `solve3D_pyramid_plan_2026-10-15.md`

Under `doc/refactoring_notes/planned/`:

- `phase5_workflow_gates_and_nightly_runner_handover.md`
- `single_area_tests_handover.md`
- `stream_area_tests_handover.md`
- `rnd_seed_parameter.md`
- `refine_tomo.md`
- `next_release/staged_domain_driven_design_refactor.md`

Elsewhere:

- `doc/policies/test_environment_policy.md`
- `doc/policies/3D/refine3D_policy.md` (section 11.1)
- the run kit `~/agent_runs/runs/speed2d_mac/plan.md`
