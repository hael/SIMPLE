# A test suite that knows the truth: correctness, robustness and multi-state

Implementation note, 2026-10-16. **Phases 0 to 2 done (see Progress); the
later phases are planned.** Written
after the run `ms_speed`, which could not tell whether its multi-state
schedule change was correct: no test knows the state labels, and the checks it
used on beta-galactosidase measured run-to-run scatter, not correctness. The
maintainer's points: multi-state results must be validated against known
state labels, so on simulated data only; every high-level test must be shown to
test what it claims; and robustness testing must be separated from correctness
testing. A weekend on the Oracle Linux test machine is available, and a suite
that is both robust and fast is the goal.

## Rulings of 2026-10-16

- **Heterogeneity models as proposed.** Compositional beta-galactosidase
  states (chains removed) for the correctness tier, a conformational state (a
  rigid chain rotation) for the hard profile, blob phantoms for fast logic
  tests.
- **`simulate_particles` learns several volumes and populations.** It no
  longer forces state 1; the mixture is made by the production program, not
  by a test-side workaround (section "The simulator").
- **Repair what can be repaired**: the two failing multi-state tests and the
  weaknesses the audit finds, each with the profile that justifies the
  change.
- Threads in correctness gates: the CTest thread count, with floors
  calibrated across seeds and thread settings (proposed, not objected to).
- **No outside driver.**
  - Every high-level test is a SIMPLE commander, so the runner is one too:
    `simple_test_exec test=nightly` runs the whole suite.
  - A robustness profile is a commander that executes a test commander
    under command lines it generates itself (seed offsets, thread counts).
  - Nothing that runs, profiles or compares tests is a Python script. The
    existing registry lint (`scripts/check_test_registry.py`) stays a lint.

- **A review of this plan (2026-10-16)** found eight design faults, all
  corrected below:
  - one state mapping for labels, maps and poses;
  - a discrimination metric with a duplicated-map control;
  - declared point groups and one frame and hand per state;
  - explicit fixture modes;
  - null distributions for unequal populations;
  - CPU reservations instead of thread counts;
  - a dispatcher with a failure and deadline contract;
  - paired production variants for the schedule comparison.

## The plain learning

The policy (`doc/policies/test_environment_policy.md`, sections 2 and 6)
already asks for the right things: truth instead of completion, floors with a
stated origin, a negative control where cheap, reproducible draws. The suite
does not yet live up to it in four ways:

1. **Multi-state has no working truth test.**
   - Two tests involve several states, `flex_pca_blobs` and
     `solve3D_pyramid_blobs`, and both fail on the current code before the
     state-separating stages.
   - No test runs multi-state `solve3D`, multi-state `refine3D` or
     `refine3D_states` against known labels.
   - The state metrics that exist use many-to-one majority recoding, which
     also passes identical maps, rather than a one-to-one matching.
2. **Correctness and robustness are confounded.**
   - Multi-threaded runs are not reproducible, so a single run's metric is a
     draw from a distribution.
   - Floors were set from one or two runs, and several sit inside that
     distribution:
     - `solve3D_addon`: the 5° pose floor against values of 1.5-6°;
     - `simulated_workflow_6vxx`: about 80 particles, 6 of 8 seeds pass;
     - `simulated_workflow_1jxy`: a correlation floor of 0.80 against
       recorded values of 0.80-0.90.
   - So a failure does not say whether the code is wrong or the draw was
     unlucky. That is the stop pattern of the last two runs.
3. **Claims are not checked against what runs.**
   - Each claim was found unexercised or unchecked:
     - `flex_pca_blobs` claims the class-average state refinement, but its
       `refine3D_states_cavgs` stops before the neighbourhood stage, so the
       2,500-direction cap of `refine3D_states` was never run by any test;
     - `simulated_workflow` claims orientations are compared, and only the
       map is;
     - `mini_stream` claims picking and 2D, and checks a particle count and
       internal consistency.
   - No test asserts that the path it claims actually ran.
4. **Some gates cannot fail, or cannot fail for the right reason.**
   - Map floors of 40 Å against 4-6 Å results.
   - Pose checks only relative to another run, with no absolute floor.
   - A pyramid purity floor of 0.45, barely above what random labels reach.
   - A docking correlation that turns a non-finite value into 1.

## Four kinds of test, kept apart

| kind | question | how it is run | how it passes |
| --- | --- | --- | --- |
| **identity** | did a refactoring change any result? | one thread, fixed seed, the same case on two builds | declared outputs byte-identical (run identifiers and path hashes excluded by a declared list) |
| **correctness** | does the code recover the truth it was given? | a single run on simulated data with a large, known signal | truth metrics against floors that sit far from both the observed distribution and the negative control |
| **robustness** | how reliably, and how well, across draws? | the same case over many seeds and thread settings | a distribution (median, spread, worst, failure rate against the correctness floor) recorded and compared over time; never a single-run gate |
| **performance** | how fast? | benchmark reports, as today | timings, never in a pass/fail test |

Consequences:

- **A correctness gate must be insensitive to the draw.** Its case is built so
  that the right algorithm passes on every seed and thread count, by a wide
  margin. Its floor is set from a robustness profile of at least 20 runs:
  below the worst run by a stated margin, and well above what a broken
  algorithm or the negative control gives. If no floor satisfies both, the
  case is not a correctness gate; it becomes a robustness profile, or its
  signal is raised until it can be one.
- **Identity is the only test of a refactoring that claims identical
  results.** It never uses multi-threaded runs. The run `ms_speed` settled
  its assignment change this way (one-thread runs byte-identical) after its
  multi-threaded checks had stopped it.
- **Robustness profiles inform decisions; they do not stop runs.** A
  results-changing change, like the multi-state schedule, is judged by
  comparing robustness profiles before and after, on cases with known
  truth.

## How every high-level test is validated

Each high-level test gets a test card, kept in the policy, and the card is
checked by running the test, not by reading it.

1. **Claim.** One sentence: the guarantee it checks.
2. **Path assertions.** The test asserts, from logs or outputs, that the
   claimed path ran. Examples:
   - the refine mode and number of states of each stage;
   - that a distributed step used several worker processes;
   - that the stage the claim concerns was reached.
3. **Truth metrics.** Every gate compares with the generating truth (poses,
   shifts, defocus, positions, maps, state labels). Completion and
   consistency checks stay, but as checks, not as the gate.
4. **Calibrated floors.** A robustness profile of the test (20 seeds at its
   CTest thread count; one thread as well where the test is deterministic)
   gives each metric's distribution. The card records the median, the worst
   run, the floor and the margin, in units of the observed spread.
5. **Sensitivity.** A negative control shows the gate can fail for the right
   reason: run the test once with the truth deliberately broken (shuffled
   state labels, randomised poses, the wrong truth map) and confirm that the
   gate fails. A gate that passes a broken input is not a gate.
6. **Determinism class and cost.** Reproducible at one thread or not; wall
   time and threads.

What the inventory of 2026-10-16 already shows, test by test, is in the
appendix. Every one of these becomes a card. The two broken multi-state tests
are either repaired or replaced by the new multi-state tests below; that is a
decision for the maintainer.

## Multi-state tests with known labels

### The simulator

`simulate_particles` (`simple_commanders_sim.f90`, `exec_simulate_particles`;
interface `simple_ui_sim.f90`) today reads one volume and sets every
particle's state to 1. It learns:

- **Several volumes.** `nstates=K` with `vol1` to `volK`. Every particle is
  projected from the volume of its state. The particle loop keeps its index
  order, so the draws keep their order; the K padded volumes are held
  together, a cost stated in the documentation.
- **States from the orientation table, or drawn.** A `state` field in
  `oritab` is respected; values outside 1 to K are refused.
  - Without one, the states are assigned from given population fractions.
  - The counts are exact: the fraction times the number of particles,
    rounded, with the remainder to the largest fractions. The assignment is
    shuffled with the seeded generator.
  - The fractions default to equal. Their key follows an existing pattern
    for a list of reals in the parameters, or a small text file if there is
    none; the run records which.
- **The truth written out.** The output orientation document carries every
  particle's state with its pose, shift and CTF, as it does today for one
  state.
- **One state is unchanged.** With `nstates` absent or 1 the program's
  outputs are byte-identical to today's at one thread.
- Its documentation (`doc/user_guide/simulate_particles.md` and the
  interface text) describes the mixture.

### The fixture

A test-side fixture builds a heterogeneous dataset whose labels, poses,
shifts and CTF are all known, by calling `simulate_particles` once with all
states. It has three explicit modes, because existing tests rely on
contracts a pure search fixture would break:

| mode | poses | in the project | for |
| --- | --- | --- | --- |
| search | random, with shifts | no poses, no states (the truth beside it) | the new multi-state tests: the programs must find poses and states |
| known poses | random or spiral | the true poses, no states | tests of state separation given the poses, as `flex_pca_blobs` parts 1 and 2 run today |
| exact views | a number of views, each repeated exactly | the true poses, and the 2D class equal to the view | class-average tests, as `flex_pca_blobs` part 3 and `solve3D_pyramid_blobs` run today |

The planted-state fixture of `flex_pca_blobs`
(`create_flex_pca_phantom_fixture`) simulates one stack per state and
interleaves them, with equal populations. The flex and pyramid tests move
onto the new fixture only with the mode whose contract they need. The move
keeps their contracts, not their bytes, and the moved tests are profiled
again.

- **States in one frame.** Every state's coordinates are edited from one
  model and then translated by one shared vector, so the truth maps share a
  frame.
  - The present `pdb2mrc` centres each model on its own centre, and it
    forces centring whenever a coordinate is negative, as 1JYX's are.
  - So the shared translation (the full model's centre to the box centre)
    is applied to every state before `pdb2mrc` runs without centring.
- **Two families.**
  - **Blob phantom** (box 64, seconds per simulation), from the existing
    `write_phantom_truths`: a common body and a state-specific lobe. For
    fast logic tests.
  - **Beta-galactosidase states** (1JYX, embedded; waters and other het
    atoms dropped; box 96 at about 2.2 Å per pixel):
    - compositional states: the tetramer, the tetramer without one chain,
      and a dimer;
    - a conformational state: one chain rotated rigidly about the axis
      through its centre of mass parallel to the tetramer's long axis, by a
      stated angle.
- **Every truth declares its point group.** The tetramer is D2, the
  tetramer without one chain is C1, and the dimer has the group of its two
  chains. The multi-state programs run in C1, so a symmetric state's
  particles have several equivalent true poses. A fixture check confirms
  each declared group: the truth map correlates with itself under every
  operator of its group (above 0.99), and not under the operators of the
  next larger group in the same family.
- **Mixture.** Made by the extended `simulate_particles` in one call:
  - given populations, deliberately unequal (for example 50, 30 and 20%);
  - in search mode, random poses (not a spiral) and shifts within ±2-3
    pixels;
  - a defocus range and a stated SNR.
- **Truth kept apart.** The true state, pose and shift of every particle go
  into a truth document beside the project. In search mode the project the
  program sees carries no poses and no states. Starting volumes, where a
  test needs them, are the truths low-passed to a stated resolution.
- **Cached.** A fixture is generated once per build from its parameters and
  mode, and reused by every test and every robustness run that asks for the
  same one.

### The metrics

Shared, in `simple_test_truth_metrics.f90`, so every multi-state test uses the
same definitions. Every metric has a unit test on synthetic inputs with a
known answer.

**One state mapping for everything.**

- The estimated states are matched to the true states once: the one-to-one
  assignment that maximises label agreement over all permutations (up to
  six states).
- That single mapping is used for:
  - the populations;
  - the maps;
  - the poses and shifts;
  - the discrimination.
- The best map matching is also derived on its own, from the docked
  correlations, and the test asserts that the two mappings agree.
  Otherwise correct labels with swapped maps could score perfectly.

**Labels.**

- **Primary: the adjusted Rand index.** It is permutation-invariant, 0 in
  expectation for independent partitions whatever the populations, and 1
  for a perfect partition.
- **Secondary: the accuracy under the mapping**, judged against its own
  null distribution, never against 1/K:
  - **The null.** The true labels are shuffled across particles with the
    population counts kept. The matched accuracy is recomputed for each of
    1,000 shuffles, giving the null's median and 99.9th percentile. With
    populations of 50, 30 and 20% the null sits near 0.38.
  - **The collapsed predictor.** Every particle in one state scores the
    largest population (0.50 there), and an adjusted Rand index of 0.
  - **The floor.** A floor on accuracy must clear both the null's 99.9th
    percentile and the collapsed predictor by the margin of the floor rule.
- **Unit tests.**
  - A global renaming of the estimated state identifiers leaves the
    accuracy and the adjusted Rand index unchanged.
  - Known partitions give their published adjusted Rand values.

**Populations.** The error per state, under the mapping.

**Maps, and whether they discriminate the states.**

- **Docking.** Each estimated map is docked (both hands) onto its mapped
  truth. That gives one rigid transformation per estimated state: rotation,
  shift and hand.
- **Fit.** Per state: the band correlation and the masked FSC 0.143 with its
  mapped truth.
- **Discrimination**, the guard against duplicated maps:
  - The difference mask is where the truths differ (the absolute difference
    of each truth from the others above a stated fraction of its maximum,
    dilated).
  - In the docked frame and inside that mask, each state's margin is its
    correlation with its own truth minus its highest correlation with
    another truth.
  - For every pair of states, the difference of the two docked estimated
    maps is also correlated with the difference of their truths.
  - Both must clear floors. Copies of one map have margins near 0 and a
    difference near 0 (or not finite, which fails).
- **The duplicated-maps negative control.** Every estimated map is replaced
  by one of them (and once by their mean), with the truths kept distinct.
  The fit metrics may still pass; the discrimination metrics must fail.

**Poses and shifts.**

- **Frame.** Under the state's docking transformation, each particle's pose
  and shift are carried into the truth frame by one composition rule, with
  the hand's mirror applied to poses and shifts alike.
  - A unit test pins the convention: a truth rotated, shifted and mirrored
    by a known transformation, with poses moved accordingly, must give zero
    error.
- **Error.** The rotation error is the smallest geodesic distance over the
  operators of that truth's declared point group, and only those: no group
  is learnt from the data.
  - So `sym_pose_error`, which today learns rotation clusters from the
    errors, is replaced, not moved: a wrong frame could become a learnt
    "symmetry" cluster.
  - The error is computed over the particles whose state is correct under
    the mapping.
- **Shifts.** The error after the same transformation.

A non-finite value fails its metric; it never maps to a pass.

### The tests

| test | what it runs | signal | kind |
| --- | --- | --- | --- |
| multi-state assignment | `refine3D nstates=3 refine=prob`, and once `prob_neigh`, from the blob truths low-passed to 20 Å, one and several threads | large lobes, SNR 1-2 | correctness: adjusted Rand index near 1 and accuracy far above its null |
| multi-state `solve3D` | `solve3D_cavgs nstates=3` then `solve3D nstates=3 cavg_ini_ext=yes`, and `solve3D nstates=3` from the particles (the `shc` route), on the compositional states; by the ruling of 2026-10-17 the class-average step is tested on simulated images of the states imported with `import_cavgs` as its class averages (known state and view per image, unequal populations), scored per imported image against its known state (the single state mapping, the adjusted Rand index, per-state recall) and per-state maps, in place of `make_cavgs` on the views fixture; by the rulings of 2026-10-18 the class averages are about 600, each with even and odd halves that are independent noise realisations of the same state, view and shift (no CTF, shifts within about a pixel), multi-state `solve3D_cavgs` runs direct `shc` in stages 1 and 2 as particle `solve3D` does, its `prob_athres` default is 10 as in `solve3D` (the taper deleted), and the collapse is localised by controlled experiments (particle `solve3D` on the same images, one imposed setting at a time) before it may be called a design question | missing chains | correctness with calibrated floors; robustness profile |
| state refinement from a consensus | `refine3D_states` on a single-state solution of the compositional states, through its neighbourhood stage; by the ruling of 2026-10-17 state 0 means exclusion, so the gate requires every true state to be resolved (no true state left in state 0) | missing chains | correctness, with a path assertion that the neighbourhood stage ran |
| hard heterogeneity | the same three programs on the conformational states, small rotation, low SNR, one small state | a few degrees of domain motion | robustness profile only, the basis for decisions such as the stage-5 schedule |

Negative controls, each run once at calibration:

- **Shuffled labels:** the truth labels are shuffled across particles with
  the counts kept. The accuracy must fall inside its null distribution and
  the adjusted Rand index to about 0.
- **Duplicated maps:** the estimated maps are replaced by copies of one of
  them, and once by their mean. The discrimination metrics must fail.
- **Identical states:** the fixture is built with all states equal. The
  adjusted Rand index must sit near 0, and the discrimination metrics must
  not report distinct states.

## Fast

- **Small where signal allows.** Correctness cases are sized to the
  smallest box and particle count at which their signal still separates
  them from the negative control by a wide margin. The blob tier runs in
  minutes, the beta-galactosidase tier in tens of minutes.
- **Fixtures cached** per build instead of simulated in every test.
- **Packed by CPU reservation, not by thread count.**
  - A test's thread count understates what it uses. `simulated_workflow_1jxy`
    declares one thread and launches two worker processes in its
    distributed step; `flex_pca_blobs` takes four worker threads and uses
    eight.
  - So every registry entry carries its peak CPU reservation: the cores in
    use at the busiest moment, master and workers together. That number is
    first derived from the test's settings and then confirmed by the
    profile, which samples the CPU use of each execution's process group.
  - `test=nightly` and `test=profile` pack executions so that the sum of
    their reservations stays within the core budget.
  - The CTest entries, kept for continuous integration, declare the
    reservation as their processor count instead of `RUN_SERIAL`, so
    `ctest -L highlevel -j 24` packs them by the same numbers.
- **A correctness tier and a profile tier.**
  - The correctness gates are what `test=nightly` runs, aimed at under 45
    minutes on 24 cores.
  - The robustness profiles run on the weekend or on demand through
    `test=profile`.

## Tooling

All of it in SIMPLE, as test programs of `simple_test_exec`. Each is a
commander, and every test execution it starts is a command-line object it
builds, run in a child process.

- **The registry.** One table in SIMPLE lists every high-level entry with:
  - its sub-suites;
  - its master thread count and its worker threads and worker processes;
  - its peak CPU reservation (section "Fast");
  - its deadline (today's CTest timeout);
  - its kind (correctness gate or robustness profile).

  `test=nightly`, `test=profile` and the CTest registration all read the
  same entries. The registry check (`scripts/check_test_registry.py`, an
  existing lint) keeps `production/CMakeLists.txt` in line with it.
- **The seed offset.**
  - `simple_test_exec ... seed_offset=<n>` adds n to every fixed seed of the
    test, with 0 giving today's draws.
  - The offset enters through `fix_test_seed` and `set_fixed_seed`. The two
    high-level commanders that call `seed_rnd_fixed` directly (the
    simulated workflow and the SINGLE workflow) are moved onto
    `set_fixed_seed`.
- **The dispatcher**, shared by `test=nightly` and `test=profile`.
  - **Not `simple_qsys_ctrl`.** Its `schedule_jobs` cancels the remaining
    jobs and stops the master when one part fails, which is right for a
    distributed workflow and wrong for a test suite. Its local queue also
    sets no deadline.
  - **Built on `simple_forked_process`.** That module forks a child that
    executes a commander under a given command line, polls its status
    without blocking, and can send it SIGTERM and SIGKILL. The stream
    pipeline already uses it.
  - Each child:
    - becomes the leader of its own process session (`c_setsid`, already
      bound in `unix_unistd.F90`), so the workers it launches belong to
      it;
    - changes into its own directory;
    - sets its seed offset and threads;
    - executes the test commander.
  - **Contract.** A failed child is recorded with its exit status and the
    others go on. A child past its deadline gets SIGTERM, then after a grace
    period SIGKILL, sent to its whole process group, so no worker outlives
    it; it is recorded as timed out. When the dispatcher is itself
    interrupted, it stops every child group before exiting. Every execution
    ends as passed, failed or timed out, with its wall time, its peak CPU
    use and its `metrics.tsv` if one was written.
  - **Its own tests.** Test programs that fail at once, hang, or start a
    worker process that hangs are dispatched together with a passing one.
    The dispatcher must report each correctly, leave no process behind
    (checked by process group), and finish within the deadlines plus the
    grace period.
- **`simple_test_exec test=nightly`.**
  - It runs every correctness entry of the registry once through the
    dispatcher, each in its own directory, packed by reservation.
  - It collects every entry's exit status, wall time and `metrics.tsv`.
  - It writes a nightly report: one row per entry and metric, with value,
    floor and pass, plus the wall times.
  - It exits non-zero if any gate fails.
- **`simple_test_exec test=profile target=<entry|all> [suite=] nseeds=<n> [nthr=<list>]`.**
  - It generates one command line per seed offset and thread setting, runs
    them through the dispatcher as `test=nightly` does, and collects their
    `metrics.tsv`. Failed and timed-out executions count in the failure
    rate.
  - It writes each metric's distribution: runs, median, robust spread,
    worst, and the failure rate against the floor. The statistics use
    `simple_stat`.
  - With two profiles given, it writes their comparison.
- **`simple_test_exec test=identity target=<entry> dir1=<run> dir2=<run>`.**
  - It byte-compares two run directories over the entry's declared outputs
    (in its registry row; hence `target=`, added in Phase 0).
  - It ignores the declared run-specific fields: the run identifier in
    `solve3D_manifest.txt`, the layout digest of `sigma2_state.bin` and the
    manifest's `stack_digest`; Phase 0 found that two directories also differ
    in the manifest's own `layout_digest`, its `checksum` and the paths of its
    `artifact` lines, and that project files hold absolute paths, so projects
    are compared by the records of their particle and class segments.
  - It replaces the ad hoc scripts of the last two runs.
- **Test cards** in the policy, one per high-level test, with the
  calibrated numbers.

## Programs and runs that need inputs

**In the repository**, every high-level CTest entry is hermetic: it makes its
own data from embedded models or phantoms. Five programs reachable through
`simple_test_exec` need inputs a user supplies:

| program | what it is | proposal |
| --- | --- | --- |
| `pcg_frac_update` | real assertions (additivity of the reconstruction sums, the update weighting, the continuation replay) on a user's `reconstruct3D` project | run it on a project made by the fixture, so it becomes a hermetic test with known inputs |
| `rec3D_backends` | gridding against PCG on a user's project; measurement only, no thresholds | a developer tool (policy, section 3.4), listed as such, not a test |
| `solve3D_addon_snapshots`, `write_state_weights_labels`, `write_state_weights_mixed` | generators of inputs for manual runs | tools, listed as such |

**Real-data runs** (beta-galactosidase, PfCRT, msp1, the ribosome, the
benchmark datasets of `SIMPLE_data_testing`) need datasets and have no truth.

- They are robustness and performance evidence, never correctness:
  - does the workflow still complete;
  - does the resolution hold its dataset limit;
  - how long does it take;
  - how much do the results scatter across seeds.
- Multi-state real data can show stability (populations and per-state
  resolutions across seeds), never correctness: the labels are unknown.
- **Automated** as a dataset layer on the test machine:
  - a manifest per machine (the dataset's path, the exact command line, the
    resolution limit, the reference timing), kept outside the repository
    because the paths are the machine's;
  - run by the same SIMPLE runner as the simulated tests, as a test program
    that reads the manifest and executes each dataset's command line.
- The beta-galactosidase and PfCRT cases of the speed runs are its first
  entries. The layer is designed in this run's report and built after it.

## Execution: the run `test_truth`

The weekend on the Oracle Linux test machine, as an unattended run in seven
phases. Nothing is committed; the maintainer reviews the complete diff. Every
phase records itself in the progress section.

### Common rules

- **Base.** The HEAD the copy starts from. The run's own tests are the
  subject, so a test that fails on the base is a finding, not a stop.
- **Builds.**
  - An out-of-tree Release build with tests (`BUILD_TESTS=ON`) for the
    profiles and the high-level runs.
  - The Debug build of `./compile_debug.sh`, which runs the fast gate, for
    the gates.
  - Kept per phase as `scratch/keep/build_phaseN`.
- **Packing.**
  - Profile runs are packed on the machine's 24 cores: the sum of their
    peak CPU reservations at most 24, with free memory watched.
  - Nothing in this run is a timing comparison, so the machine need not be
    idle. Wall times are recorded as the tests' cost.
- **Seeds.** A profile uses seed offsets 0 to 19 (0 is the CTest draw), at the
  entry's CTest thread count. One-thread entries are also run twice at
  offset 0 and compared byte for byte.
- **Floors from profiles.** Metrics where higher is better are written here;
  the opposite case is the mirror image. Let `s` be the robust spread
  (1.4826 times the median absolute deviation, or the observed range over 4
  if that is zero). A floor is admissible when it lies:
  - at most the worst of the 20 runs minus `s`;
  - at most the median minus `3s`;
  - at least the negative control's value plus `3s`.

  If no value satisfies all three, the case is not a correctness gate. It is
  recorded as a robustness profile, or its signal is raised until it can be
  a gate. Every floor is written next to its check with its origin (the
  profile, its date, its median and spread).
- **Identity.** A change the plan calls identical keeps the one-thread
  outputs of the affected tests byte-identical. `test=identity` checks this
  (until Phase 0 builds it, a byte comparison by hand), with its declared
  list of run-specific fields:
  - the run identifier in `solve3D_manifest.txt`;
  - the layout digest of `sigma2_state.bin`;
  - the manifest's `stack_digest`;
  - since Phase 0, also the manifest's `layout_digest`, `checksum` and
    artifact paths, with projects compared by their particle and class
    records (section "Tooling").

  Multi-threaded runs never decide identity.
- **Fix as you go** (ruling of 2026-10-17, replacing "findings, not
  stops"). Every bug and issue found is fixed in the phase that finds it,
  not deferred and not only recorded:
  - a defect in this run's own code, or in production code that a test
    exposes, is fixed at once, its file added to the file table for the
    current phase with a one-line reason; the fix is verified by the
    identity check where the change is meant to be identical, otherwise by
    the test that found it and a profile before and after;
  - where the cause is a question of algorithm or design (a workflow that
    merges or loses states), it is diagnosed to its root cause first: a
    defect is fixed; a design decision is written, with the diagnosis and a
    proposed fix, into the progress section under "For the maintainer", and
    the work goes on (a stop only if the open question blocks the rest of
    the phase);
  - a floor never moves except by the floor rule.

  The run stops only when:
  - a change the plan calls identical cannot be made identical;
  - a repair needs a design the plan does not give;
  - a production fix would change the results of workflows beyond the
    defect.
- **Gates.**
  - The fast gate.
  - The CTest unit suites of the touched areas.
  - `python3 scripts/check_test_registry.py .` whenever tests or their
    registration change; a new CTest entry raises `SIMPLE_CTEST_BUDGET` with
    its reason.
- **Release 4 rules.**
  - No backwards compatibility.
  - `THROW_HARD` and `THROW_WARN` arguments are never continued with `//&`.
  - No element test of an allocatable on the line that tests its allocation.
  - No `goto`.
  - Dead code is deleted.
  - Executable bits are checked with `git diff --summary`.
  - The Windows build is not a concern.

### Phase 0: tooling, and the existing tests report their numbers

1. **The seed offset.**
   - `simple_test_exec ... seed_offset=<n>` adds n to every fixed seed of a
     test, through `fix_test_seed` and `set_fixed_seed` in
     `simple_test_utils.f90`.
   - The commanders that call `seed_rnd_fixed` directly move onto
     `set_fixed_seed`: `simple_commanders_test_highlevel_workflow.f90` and
     `simple_commanders_test_single.f90` on 2026-10-16. A search confirms no
     other test draws from a fixed seed outside the two routines.
   - It is a key of the test front door, not of `parameters`. The default 0
     gives today's draws, byte for byte (identity).
2. **Every high-level test reports its gate values** to `metrics.tsv`,
   through `test_gate`'s `metric` and `report`.
   - No floor, check or outcome changes.
   - The one-thread suites stay byte-identical (identity).
3. **The registry** of high-level entries (section "Tooling"), and the
   CTest registration and the registry check made to agree with it.
4. **The dispatcher, `test=nightly`, `test=profile` and `test=identity`**
   (section "Tooling"), as commanders of `simple_test_exec`.
   - The dispatcher's own tests: failing, hanging and worker-hanging
     children beside a passing one.
   - Unit tests for the pure parts: the statistics of a profile, the packing
     of executions into the core budget by reservation, and the identity
     comparison with its ignore rules.

Exit criteria: the gates pass; the identity of items 1 and 2 holds;
`test=profile` and `test=identity` work on one entry, and `test=nightly` on
the whole suite; their use is described in the policy.

### Phase 1: the profile of the existing suite

- `simple_test_exec test=profile target=all nseeds=20`: every high-level
  entry, 20 seed offsets at its CTest thread count. The one-thread entries
  are also run twice at offset 0 and compared with `test=identity`.
- Each entry's peak CPU reservation confirmed from the sampled CPU use of
  its executions, and corrected in the registry where the derived value
  was wrong.
- A test card per entry, from the profile and the logs:
  - its claim;
  - which path assertions the logs support;
  - its truth metrics;
  - each floor against its distribution, with its margin in units of `s`;
  - whether a negative control exists or could;
  - determinism and cost.
- Nothing is changed in this phase; the cards record findings.

Exit criteria: every entry profiled, or recorded as failing deterministically;
every card drafted in the progress section.

### Phase 2: `simulate_particles` with several states, the fixture, the metrics

- **The simulator** (section "The simulator"):
  - identity for one state;
  - a test that the mixture is right: the state counts are exact, every
    state's volume is used, and every image correlates better with the
    projection of its own state's volume at its true pose than with the
    other states'.
- **The fixture**, both families, in its three modes, cached per build, with
  the check of each truth's declared point group.
- **The shared metrics** in `simple_test_truth_metrics.f90` (section "The
  metrics"), with unit tests:
  - the matching on known confusion matrices, and its invariance under a
    global renaming of the estimated states;
  - the adjusted Rand index on known partitions;
  - the null distribution of the accuracy for unequal populations, and the
    collapsed predictor;
  - the assertion that the label and map mappings agree, failing on correct
    labels with swapped maps;
  - the discrimination metrics failing on duplicated maps and on their
    mean;
  - the frame and hand convention for poses and shifts on a known
    transformation;
  - the pose error under a declared group, against a wrong frame that a
    learnt cluster would have accepted;
  - a non-finite value failing.
- The pose error under declared point groups replaces `sym_pose_error` in
  `cont_refine3D_1jxy`, with that test's D2 truth declared.

Exit criteria: the gates pass; identity holds for one state; the new tests
pass.

### Phase 3: the multi-state assignment test

- `refine3D nstates=3` with `refine=prob`, and once with `refine=prob_neigh`,
  on the blob family from the truths low-passed to 20 Å.
- Profile over 20 seeds at one thread and at the entry's thread count.
- The negative controls: shuffled labels, and identical states.
- Floors from the profile by the rule above; registered as a high-level
  entry.

Exit criteria: the gates pass; the floors are admissible, or the case is
recorded as a profile with the reason.

### Phase 4: the multi-state workflow tests

- **Tests on the compositional beta-galactosidase states**, with path
  assertions on each stage's refine mode, directions and number of states:
  - `solve3D_cavgs nstates=3` then `solve3D nstates=3 cavg_ini_ext=yes`;
  - `solve3D nstates=3` from the particles (the `shc` route);
  - `refine3D_states` from a single-state solution, through its
    neighbourhood stage.
- Their profiles, negative controls and floors, as in Phase 3.
- **The hard profile** on the conformational states (robustness only).
- **The schedule comparison**, with paired production variants:
  - `75eab7f90` differs from the present code in 88 files, reconstruction
    among them, so it cannot isolate the schedule.
  - Instead, a scratch copy of the present tree with only the schedule's
    production changes reversed:
    - the switch to `prob_neigh` at 2,500 directions, the direction cap and
      the stage-2 exit in `simple_solve3D_controller.f90`;
    - the `scan_exit` rule in `simple_convergence.f90` and its key in the
      parameters;
    - `STAGE2_NSPACE` in `simple_refine3D_states_controller.f90`.
  - Its exact diff against the present tree is recorded. The variant is not
    part of the delivered change.
  - Both builds run the same test code, the same cached fixtures and the
    same seed offsets.
  - Reported, not gated: the adjusted Rand index and accuracy, the per-state
    maps and their discrimination, and poses, before and after.

Exit criteria: the gates pass; every test's floors admissible or the case
recorded as a profile; the comparison reported.

### Phase 5: repairs and the audit's actions

- **The flex and pyramid tests onto the new fixture**, each in the mode
  whose contract it needs: known poses for `flex_pca_blobs` parts 1 and 2,
  exact views for its part 3 and for `solve3D_pyramid_blobs`. Profiled
  again after the move.
- **The two failing multi-state tests:**
  - `solve3D_pyramid_blobs`: its `solve3D_cavgs` step plans a cropped box
    larger than the native one;
  - `flex_pca_blobs` part 3: an empty state stops the distributed PCG;
  - and the pyramid test's metrics: the single state mapping, the
    discrimination metrics, a floor from a profile, a non-finite
    correlation failing.
  - Production defects behind them are fixed where the fix is local.
  - If not, the finding is recorded and the run stops for a ruling.
- **Every other entry**, by its card:
  - path assertions added;
  - floors moved to admissible values with their origin written;
  - fragile cases redesigned and profiled before and after: the pose floor
    of `solve3D_addon`, and the 80 picked particles of
    `simulated_workflow_6vxx`;
  - claims corrected;
  - `pcg_recon` relabelled library;
  - `simulate_particles` checking images against projections.
- **The programs that need inputs**, as in the section above.
- **Packing.** `test=nightly` packs its executions within the core budget.
  The CTest entries replace `RUN_SERIAL` by their processor count, with a
  resource lock only where entries share a resource. Measure the wall time
  of `test=nightly` and of `ctest -L highlevel -j 24`.

Exit criteria: the gates pass; every change carries its profile; the packed
suite's wall time recorded.

### Phase 6: documents and close

- `test_environment_policy.md`:
  - the four kinds of test;
  - the floor rule;
  - the cards;
  - the programs that need inputs and the dataset layer;
  - the registry, `test=nightly`, `test=profile` and `test=identity`.
- `simulate_particles.md`.
- A last `simple_test_exec test=nightly` on the final build.
- The report beside this note: what every test verifies now, what it costs,
  the profiles, and the findings.

Exit criteria: the gates pass; the report is written.

### Phases 7 to 9: the continuation run `test_truth_sync`

Phases 0 to 6 ran as `test_truth` on 25d256886. Its work was merged onto
master 43de9b340 on 2026-10-18 (Progress, Phase 6), and the continuation
finishes the plan on the commit that holds the merge. The old run's builds,
profiles and fixture caches under `~/agent_runs/test_truth/scratch/keep/` are
read, never written. The rulings of `test_truth` stay in force; two of them
reached no phase of it (floors at one thread and on the maintainer's Mac;
how the multi-state workflows are tested) and are carried out here.

### Phase 7: the merged code, and what still holds

- **Gates** on the merged code: the Release build with tests, the Debug fast
  gate, the registry lint, `check_descr.py`.
- **The speed week on x86.** The speed week (04e520e78) was judged
  byte-identical at one thread on the maintainer's Mac only. Compare a build
  of 25d256886 (the old run's `build_base`) with a build of 04e520e78 (from
  an export of that commit in scratch) at one thread, byte for byte, on the
  cases of `identical_speed_run.md`, section "The identity test":
  `simulated_workflow_1jxy` with both steps, the three-state, correlation
  and PCG cases (MSI, CCI, PCI there), with the base build's
  `simulated_workflow_1jxy` run directory as their fixture.
  - The speed week keeps back-filled projection labels in volume assembly,
    a fix that changes results where labels were back-filled. Find which
    callers and cases reach it; on those, every difference must be the kept
    labels or follow from them.
  - Any other difference is a defect of the speed week on x86: find the item
    that causes it (`identical_speed_run_report.md` lists them) and fix it in
    the merged code.
- **What still holds.** Entries that run `solve3D`, `refine3D` or their
  class-average forms are profiled again in Phase 8 without a check: master
  changed their results. Each other entry runs once at seed offset 0 at one
  thread on the old run's final build and on the merged build, compared with
  `test=identity`. An entry whose outputs are byte-identical keeps its
  profile; the others are profiled again.

Exit criteria: the gates pass; the speed week's verdict on x86 with its
evidence; the list of entries to profile.

### Phase 8: profiles and floors on the merged code

- Every entry on Phase 7's list: 20 seed offsets at its CTest thread count,
  and at one thread where a refactoring run checks it at one thread. Floors
  by the floor rule from both profiles, with their origin written. The
  profiles go under `scratch/keep/profiles/phase8/`.
- `simulated_workflow` at master's sizes: its correlations gated again if
  the floor rule allows.
- The two rulings not yet applied, as they say.
- Defects found are fixed as you go.

Exit criteria: every floor comes from profiles of the merged code; the gates
pass.

### Phase 9: timing, documents and close

- On the otherwise idle machine:
  - the wall time of `test=nightly` and of `ctest -L highlevel -j 24`,
    beside Phase 5's numbers, saying which of master's changes lie between
    them;
  - the speed week's timing. A, A-shm, MS, M and E as the runs `prob_speed`
    and `ms_speed` last ran them, on builds of 25d256886 and 04e520e78,
    twice each way, interleaved. The results go into
    `identical_speed_run_report.md`.
- Phase 6's documents, finished with the new floors and cards.
- The report beside this note.

Exit criteria: the gates pass; the report is written.

### File table

| file | why | phases |
| --- | --- | --- |
| `production/simple_test_exec.f90` | the seed offset | 0 |
| `src/main/exec/simple_test_exec_highlevel.f90`, `src/main/apis/simple_test_exec_api.f90` | dispatch of the new test programs | 0 3 4 5 |
| `scripts/check_test_registry.py` | the registry check reads the SIMPLE registry | 0 5 |
| `src/utils/simple_test_utils.f90` | the seed offset | 0 |
| `src/utils/simple_forked_process.f90` | the dispatcher's children: own session, group signals, deadlines | 0 |
| `src/utils/simple_forked_process_tester.f90` | unit tests of the own session and the group signals | 0 |
| `scripts/simple_args_generator.pl` | the command line accepts the test front door's keys (`seed_offset`, `target`, `nseeds`, `dir1`, `dir2`) without a field in `parameters` | 0 |
| `src/main/ui/simple_test/` | the test front door's key; new and corrected test descriptions | 0 2 3 4 5 |
| `src/main/commanders/test/` | the tests: reporting, fixture, metrics, new tests, repairs | 0 2 3 4 5 |
| `src/main/commanders/test/simple_test_registry.f90` | Phase 1 corrects the peak CPU reservations from the profile | 1 |
| `src/main/commanders/test/simple_test_profile_stats.f90` | the floor of a name written both as a report and against a floor is the floor's, not the median with the report's placeholder | 1 |
| `src/main/commanders/test/simple_test_tools_tester.f90` | its unit test | 1 |
| `src/main/commanders/test/simple_commanders_test_highlevel_suite.f90` | `test=profile` takes a list of entries, so that Phase 1 profiles every entry except the deterministic failures in one packed run | 1 |
| `src/main/ui/simple_test/simple_test_ui_highlevel.f90` | the help of `target=` names the list | 1 |
| `doc/policies/test_environment_policy.md` | the list of entries of `test=profile` | 1 |
| `doc/policies/test_environment_policy.md` | the new truth-metric sub-suites and what the `simulate_particles` entry now runs | 2 |
| `doc/policies/test_environment_policy.md` | the new high-level entry `multistate_assignment` in the list of entries | 3 |
| `doc/policies/test_environment_policy.md` | the new high-level entry `multistate_workflows` in the list of entries | 4 |
| `src/main/flex/run/simple_flex_pca_application.f90` | the population-floor path of flex_pca prunes states below the smaller of its two floors, so refine3D_states keeps a 20 % state of 3 000 particles (ruling of 2026-10-17: state 0 means exclusion) | 5 |
| `src/main/flex/run/simple_flex_pca_state_service.f90` | the prune fallback no longer keeps a state without mass (an empty state stopped the distributed PCG of flex_pca_blobs part 3) | 5 |
| `src/main/strategies/search/simple_matcher_smpl_and_lplims.f90` | the FSC=0.143 clamp of the matching low-pass leaves the band non-empty (single_workflow_fcc stopped on an empty Fourier range); dropped in the merge of 2026-10-18, since master reverted the clamp | 5 |
| `src/main/strategies/parallelization/simple_refine3D_strategy.f90` | shared-memory refine3D seeds the states of a multi-state start as the distributed path does | 5 |
| `src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90` | the distributed PCG carries a state's previous maps forward when its members are all in one halfset, as it does for an empty state, instead of stopping (refine3D_states_cavgs of flex_pca_blobs part 3 stopped) | 5 |
| `src/main/solve/simple_solve3D_controller.f90` | multi-state `solve3D_cavgs` runs direct `shc` in stages 1 and 2, as particle `solve3D` (ruling of 2026-10-18); the angular-window taper deleted | 5 |
| `src/main/commanders/simple/simple_commanders_solve3D.f90` | `solve3D_cavgs`' default `prob_athres` 10, as `solve3D` (ruling of 2026-10-18) | 5 |
| `src/main/params/simple_parameters.f90`, `src/main/params/simple_parameters_parse.f90`, `src/main/ui/simple/simple_ui_solve3D.f90`, `src/main/solve/simple_solve3D_utils.f90`, `src/main/solve/simple_solve3D_reseed_tester.f90`, `src/main/solve/simple_state_assignment.f90`, `src/main/solve/simple_state_assignment_tester.f90` | the key `prob_athres_taper`, its interface entry, its function, its test and its caller deleted (ruling of 2026-10-18) | 5 |
| `doc/policies/3D/solve3D_cavgs_policy.md`, `doc/policies/3D/solve3D_pyramid_policy.md`, `doc/policies/stream/stream_3D_ingestion_policy.md` | the class-average route's stages and angular window after the ruling of 2026-10-18 | 5 |
| `src/main/simple_external_reference_pose_initialization.f90`, `src/main/commanders/simple/simple_commanders_refine3D.f90` | the pose-initialization pass searches no shift for pose-less rows, and `refine3D_states_cavgs` says its rows are pose-less: the shift fitted first against an arbitrary projection biased every orientation (`flex_pca_blobs` part 3: random poses of the sub-averages) | 5 |
| `src/utils/simple_magic_boxes.f90` | autoscale never plans a crop at or above the input box (solve3D_pyramid_blobs stopped on a crop of 88 for a box of 64) | 5 |
| `src/main/commanders/simple/simple_commanders_sim.f90` | simulate_particles refuses state_pops when the orientation table carries states, an empty state, and the ctf modes it cannot apply | 5 |
| `src/utils/simple_test_utils.f90`, `src/utils/math/simple_rnd.f90` | seed offsets times a large odd stride in 64-bit integers, so offsets of one test cannot land on another test's base | 5 |
| `src/utils/simple_forked_process.f90`, `src/utils/simple_forked_process_tester.f90` | the child's signal handling is async-signal-safe, SIGHUP is reset in the child, and an exit can be detected without reaping | 5 |
| `src/utils/text/simple_string.f90`, `src/utils/text/simple_string_tester.f90` | one overflow-free FNV-1a hash, which the fixture's cache key uses, and its published test vectors | 5 |
| `production/simple_test_exec.f90` | the seed offset survives a queue resubmission | 5 |
| `src/main/ui/simple/simple_ui_sim.f90` | the interface text of the refused simulate_particles inputs | 5 |
| `doc/user_guide/simulate_particles.md` | the refused inputs (state_pops beside oritab states, an empty state, ctf modes other than yes and no) | 5 |
| `src/main/commanders/simple/simple_commanders_sim.f90` | `simulate_particles` with several states | 2 |
| `src/main/ui/simple/simple_ui_sim.f90` | its interface | 2 |
| `src/main/params/simple_parameters.f90`, `src/main/params/simple_parameters_parse.f90`, `src/main/params/simple_parameters_phases.f90` | the population key, if one is needed | 2 |
| `doc/user_guide/simulate_particles.md` | the mixture | 2 6 |
| `src/main/flex/run/simple_flex_pca_application_tester.f90` | the flex tests moved onto the new fixture; in Phase 2 it exports `write_phantom_truths`, the blob truths the fixture reuses; in Phase 5 that writes the A-B difference only when there is a B (a one-state fixture crashed it) | 2 5 |
| `production/CMakeLists.txt` | registrations, labels, processor counts, the budget | 0 3 4 5 |
| `doc/policies/test_environment_policy.md` | the tools, the cards, the rules | 0 5 6 |
| `doc/implementation_notes/planned/test_suite_truth_and_robustness_report.md` | the review report of the run, beside this note | 6 |
| `src/main/commanders/test/simple_commanders_test_highlevel_msworkflows.f90` | its header comment still described the class-average route before the rulings of 2026-10-18 (300 images); comment only | 6 |
| `src/main/commanders/test/`, `src/main/ui/simple_test/`, `src/main/exec/simple_test_exec_highlevel.f90`, `production/CMakeLists.txt`, `scripts/check_test_registry.py` | the tests, their registration and their floors on the merged code | 7 8 9 |
| `doc/policies/test_environment_policy.md`, `doc/user_guide/simulate_particles.md` | the cards, floors and tools after the merge | 8 9 |
| `doc/implementation_notes/planned/identical_speed_run_report.md` | the speed week's identity on x86 and its timing | 7 9 |
| `doc/implementation_notes/planned/test_suite_truth_and_robustness_report.md` | the review report of the run, beside this note | 9 |

## Appendix: the high-level suite on 2026-10-16

From reading the code of every entry. The profile of the weekend replaces the
"recorded" values with measured distributions.

- **`mini_stream_6vxx`, `mini_stream_1jxy`.**
  - Defocus is checked against the truth (≤ 0.1 µm).
  - Picking is checked by count only (picked/48 in 0.40-1.50; first run
    0.458).
  - The 2D result is checked only for consistency.
  - The input is the simulator's optimal average, so motion correction is
    not exercised.
- **`simulated_workflow_6vxx`, `simulated_workflow_1jxy`.**
  - One thread, reproducible.
  - Only the final map is gated: correlation ≥ 0.80 against recorded
    0.80-0.90; FSC 0.143 ≤ 40 Å against 4-6 Å.
  - CTF, motion, picking and poses are not checked, though the policy says
    orientations are.
  - 6vxx rests on about 80 picked particles and passed 6 of 8 seeds.
- **`single_workflow_fcc`, `single_workflow_wurtzite`.**
  - The map is gated well (band correlation ≥ 0.90, FSC ≤ 5 Å).
  - Poses are only compared with another run (+0.5°), with no absolute
    floor.
  - Input from a noise-free nanoparticle.
- **`pcg_recon`.** An in-process operator test of the reconstructor with
  tight numerical checks. Sound, but it is a library-tier test registered as
  high level.
- **`simulate_particles`.** Metadata and bounds only. Images are never
  compared with projections, and the applied SNR and CTF are not verified.
- **`solve3D_addon`.**
  - Eight threads.
  - Ground-truth poses: median ≤ 15°, excess ≤ 5° against observed values of
    1.5-6°, which is flaky.
  - The union map is gated (≥ 0.90).
  - Its description says one add-on; the code runs two.
- **`cont_refine3D_1jxy`.** Strong: ground-truth poses and shifts with tight
  floors, a distributed step compared with the shared-memory one, and
  symmetry-aware errors.
- **`single_atoms_stats` and the CdSe variants.** Atom recall, precision and
  RMS against the model, but on noise-free maps.
- **`species_discovery`.** Strong: every label and recall against the
  generating model, and a fixture sanity check.
- **`stream_preproc`.** Strong: defocus, drift and the integrated micrograph
  against the truth, with a wrong-movie control.
- **`flex_pca_blobs`.**
  - Planted states; parts 1 and 2 run on true poses with the true mean as
    reference, so they test state separation, not search.
  - Label accuracy floor 0.60 against 0.998 observed.
  - Map preference is relative only.
  - Part 3 fails on the current code (an empty state stops the distributed
    PCG), so `refine3D_states_cavgs` never reaches its neighbourhood stage.
- **`solve3D_pyramid_blobs`.**
  - Fails on the current code (a cropped box of 88 on a 64-pixel phantom).
  - Floors provisional; purity 0.45 against chance 0.33.
  - References are matched to their closest truth with no distinctness
    requirement, and a non-finite correlation counts as 1.
- **Across the suite.**
  - No test simulates non-zero shifts and checks their recovery.
  - No test simulates unequal state populations.
  - Particle positions from `simulate_movie` are not saved, so picking
    cannot be checked against them.

## Progress

### Phase 0, 2026-10-16: tooling, and the existing tests report their numbers

Done on the Oracle Linux test machine (24 cores) in one session. Base: HEAD
`25d256886`, unchanged at the start. Builds: gfortran 15.2.1, CMake 3.26.5,
Release with tests (`BUILD_TESTS=ON`) out of tree, and the Debug build of
`./compile_debug.sh`, which ran the fast gate.

**What was built.**

1. **The seed offset.** `simple_test_exec ... seed_offset=<n>` adds n to every
   fixed seed of a test: `fix_test_seed` (the process seed, and so the
   `SIMPLE_SEED` that commanders and workers inherit) and `set_fixed_seed`
   (every tester and high-level test). The front door reads the key and
   deletes it from the command line, so no commander sees it. The command-line
   check only accepts keys that are fields of `parameters`, so the key is added
   to the short list of keys `scripts/simple_args_generator.pl` accepts without
   a field (with `target`, `nseeds`, `dir1`, `dir2` of the new programs); the
   file was added to the file table for that reason.
   - The simulated workflow and the SINGLE workflow called the library's
     `seed_rnd_fixed` directly; they now call `set_fixed_seed` and print the
     seed they use.
   - The search for other fixed seeds found no other high-level test that
     seeds outside the two routines. Outside them: the production code of the
     class expansion of `flex_pca` (`simple_cls_expansion_strategy`,
     `simple_flex_cls_kernel`) seeds itself with fixed numbers, so those draws
     in `flex_pca_blobs` part 3 do not move with the offset (a finding for
     Phase 1); and five unit or library testers (`simple_flex_pca_tester`,
     `simple_flex_pca_pcg_tester`, `simple_flex_cls_kernel_tester`,
     `simple_eul_prob_tab_tester`, `simple_rnd_tester`) call `seed_rnd_fixed`
     directly, outside the high-level suite.
2. **Every high-level test reports its gate values** to `metrics.tsv` in the
   directory it was started in, through `test_gate`. No floor, check or
   outcome changed: each test still decides as before, and every value it
   decides on is written as a row first, so a failing run reports it too.
   - Four tests wrote their table inside a fixture directory they delete when
     they pass (`solve3D_addon`, `cont_refine3D_1jxy`, `flex_pca_blobs`,
     `solve3D_pyramid_blobs`); their table now goes to the run directory.
   - Seven tests had no table and now have one: `mini_stream` (defocus error
     per micrograph, pick ratio against its floor and its ceiling),
     `simulated_workflow` (extracted particles; for each `solve3D` step the
     map's box and sampling, its correlation and FSC 0.143 against the truth,
     the two docking correlations), `single_workflow` (the same for the map,
     the pair pose errors and their limits), `simulate_particles` (each check
     group, the e2 range, the largest defocus, astigmatism and shift
     deviations), `single_atoms_stats` (recall, precision, RMS position error,
     density correlation, the lattice constant and interior coordination of the
     binary crystals), `stream_preproc` (defocus, drift and image correlation
     per movie, and the margin over the wrong movie), `pcg_recon` (each stage's
     tolerance value).
   - The table has a fifth column, the sense of the floor: `min`, `max` (the
     floor is a ceiling), `check` or `report`; values are written in exponent
     form so that tolerances of 1e-6 keep their digits.
     `validate_reconstructed_volume` writes its decisions to the gate it is
     given.
3. **The registry** (`simple_test_registry.f90`): one row per high-level CTest
   entry with its command line, CTest threads, worker processes and threads,
   peak CPU reservation (derived from the settings: the master's threads or
   its workers', whichever is larger, since a master waits for its workers),
   deadline (the CTest `TIMEOUT`), kind (all 18 are correctness gates until
   Phase 1 profiles them), and the outputs an identity check compares (the four
   one-thread entries). `scripts/check_test_registry.py` now also holds the
   highlevel registrations of `production/CMakeLists.txt` to the registry
   (names, command lines, timeouts, threads), and reports a registration it
   cannot read. The CTest registrations already agreed; nothing changed in
   them but a comment.
4. **The dispatcher, `test=nightly`, `test=profile`, `test=identity`.** Their
   use is in the policy, section 1.5. Decisions taken on the way:
   - Each execution is a forked child that runs the test program in the
     process, through the same routers as the front door. The children are
     started in queue order while the reservations of those running fit the
     core budget, the key `ncunits=` (by default the machine's processors).
   - **The local queue starts its workers with `setsid nohup`**, in sessions
     of their own, so killing a child's process group or session would leave
     its workers running. Each child therefore exports a marker,
     `SIMPLE_TEST_EXECUTION=<dispatcher pid>.<execution>`, which every worker
     inherits; the processes of an execution are its session and every
     process carrying its marker (read from `/proc/<pid>/environ`). Deadlines
     (SIGTERM, then SIGKILL after 60 s), the end of an execution and an
     interrupt (SIGINT, SIGTERM, SIGHUP) signal all of them, and the CPU and
     memory samples (every second, from `/proc`) sum over them.
   - The children get the build's environment: `SIMPLE_PATH` is the build the
     running `simple_test_exec` belongs to (from `/proc/self/exe`, or the given
     `SIMPLE_PATH` when it runs from the build tree, as under CTest), with its
     `bin/` and `scripts/` first on `PATH`, `SIMPLE_QSYS=local`.
   - `test=identity` takes `target=<entry>` beside `dir1=` and `dir2=`: the
     declared outputs belong to the entry. Its ignore list had to be completed
     against real runs: besides the run identifier and `stack_digest` of
     `solve3D_manifest.txt` and the layout digest of `sigma2_state.bin`, the
     manifest's `layout_digest` (also a digest of the stack paths), its
     `checksum` (of the whole manifest) and the path at the end of each
     `artifact` line differ between two directories. Project files hold
     absolute paths, so they are compared by the orientation records of their
     particle and class segments. The CTF text files (absolute paths) and the
     `BENCH` timing files are not declared outputs.
   - `test=profile` with `dir1=` and `dir2=` compares two profiles. A thread
     count other than the entry's own scales the `nthr=` the entry passes, its
     reservation and its deadline.
   - Tests: the sub-suite `test tools` of `unit_parallel` (statistics of a
     profile, packing, the identity rules, the registry; 102 checks), the
     sub-suite `test dispatcher` of the platform entry `forked_process`
     (children that pass, fail with exit code 3, hang, ignore SIGTERM and hang,
     wait on a worker that hangs in its own session, leave a worker behind in
     their group or in a session of their own; and an interrupted dispatcher;
     30 checks), and a test of the new session and group signals of
     `simple_forked_process`.

**Evidence** (all under `/home/elmlundho/agent_runs/test_truth/scratch`; the bulky outputs were deleted once their
numbers were recorded, the logs, reports and gate tables kept).

- Builds: the base (HEAD) in `keep/build_base`, Phase 0 in
  `keep/build_phase0` (both Release with tests, about 85 s each); the Debug build of
  `./compile_debug.sh` with the fast gate, 15 of 15 suites passed in 9.7 s
  (`logs/compile_debug_full.log`). The registry check, `check_descr.py`, the
  platform entry `forked_process` (forked process, stream heartbeat, test
  dispatcher) and `unit_parallel` in reverse order pass. `git diff --summary`
  shows no mode change.
- **Identity of items 1 and 2.** One-thread runs of the base build and of the
  Phase 0 build at seed offset 0, compared with `test=identity` over the
  declared outputs: identical, both for runs through the front door and for
  the runs of `test=nightly`.
  - `simulated_workflow_1jxy`, both `solve3D` steps (the shared-memory and the
    two-worker one): 205 files the same, none different;
  - `simulated_workflow_6vxx`: 158 files;
  - `single_workflow_wurtzite`: 133 files;
  - `single_workflow_fcc`: 71 files, up to the point where both builds fail
    the same way (finding below).
  - Reports `logs/identity_<entry>_base_vs_{final,nightly}.txt`. A control:
    `simulated_workflow_6vxx` at `seed_offset=1` passes (seed 20261004, 61
    extracted particles instead of 51, correlation 0.817, FSC 5.7 Å) and
    differs from the base in 154 of its 159 declared files
    (`logs/identity_simulated_workflow_6vxx_base_vs_offset1.txt`), so the
    offset reaches the workflow and the check sees a different draw.
- **`test=nightly` on the whole suite**: 18 entries in 970 s of wall time at a
  core budget of 19 (the four identity runs used the rest of the machine);
  report `nightly_p0/nightly_report.tsv`, log `logs/nightly_p0.log`. 14
  entries passed. Four failed, each the same way on the base build (findings
  below): `mini_stream_6vxx`, `single_workflow_fcc`, `flex_pca_blobs`,
  `solve3D_pyramid_blobs`. Wall times (s) and peak CPU use (cores) against
  the derived reservation, for Phase 1 to confirm: `simulated_workflow_1jxy`
  970, 2.3 of 2; `simulated_workflow_6vxx` 539, 1.2 of 1;
  `single_workflow_wurtzite` 557, 1.2 of 1; `cont_refine3D_1jxy` 461, 8.8 of
  8; `solve3D_addon` 417, 8.1 of 8; `flex_pca_blobs` 226 (failed), 8.4 of 8;
  `stream_preproc` 81, 1.0 of 1; `mini_stream` 52-62, 1.7-1.9 of 1;
  `species_discovery` 42, 7.7 of 8; `single_atoms_stats` and its CdSe
  variants 12-21, 2.2-2.6 of 8; `simulate_particles` 7, 7.2 of 8;
  `pcg_recon` 1 s (too short for a sample).
- **`test=profile` on one entry**: `simulate_particles` with 3 seed offsets
  (`trial/prof_final/`) and with 2 seed offsets at 2 threads (peak CPU 2.0 of
  2 reserved; the offsets change the draws, the largest defocus deviation
  0.1987 against 0.1994 µm); `pcg_recon` at 1 and 4 threads; and
  `cont_refine3D_1jxy` with 1 seed (`envcheck/prof/`), started from a shell
  whose `SIMPLE_PATH` and `PATH` pointed at another checkout: every worker of
  the test ran `keep/build_phase0/bin/simple_private_exec` with that build's
  `SIMPLE_PATH` and carried its marker (`envcheck/workers.txt`). The
  comparison of two profiles (`dir1=`, `dir2=`) was run on the `pcg_recon`
  tables (a profile against itself, every shift 0).
- **`test=identity`**: as above, on the four one-thread entries.
- A review of the diff by a separate reader found six defects in the first
  version of the dispatcher and the profile, all fixed before the final runs:
  the thread axis of a profile did not reach entries that pass `nthr=`; a
  SIGKILL could reach a reused pid after the child was reaped; zombies counted
  as live processes; SIGHUP and the death of the dispatcher were not handled;
  a timed-out child lost its buffered output; the registry check skipped a
  registration it could not read.

**Exit criteria.**

- The gates pass: the fast gate (Debug, 15 suites), the unit suites of the
  touched areas (`unit_parallel` with the new `test tools`, the platform entry
  `forked_process` with the new `test dispatcher`), and the registry check.
  No CTest entry was added, so `SIMPLE_CTEST_BUDGET` is unchanged (39).
- The identity of items 1 and 2 holds (above).
- `test=profile` and `test=identity` work on one entry, and `test=nightly` on
  the whole suite (above).
- Their use is described in the policy, section 1.5, with the seed offset in
  section 1.4 and the table format in sections 5.2 and 6.2.

**Findings** (recorded, not acted on in this phase).

- `single_workflow_fcc` fails on the base build, deterministically at one
  thread, after 60-80 s: in the second iteration of `autorefine3D_nano` the
  matching low-pass is clamped from 1.5 Å to the FSC 0.143 resolution of
  4.77 Å, the resulting Fourier index range is empty (from 18 to 12), and the
  polar Fourier setup stops with "invalid kfromto" (`simple_polarft_core`).
  Its card in Phase 1 starts from this.
- `mini_stream_6vxx` fails on the base build: 19 of the 48 placed particles
  are picked, a ratio of 0.396 against the floor 0.40; the floor's comment
  records a first run of 22 of 48 (0.458). The same 19 by the front door and
  in the nightly, at the CTest thread count.
- `flex_pca_blobs` fails as the plan says (state 2 of the class-average part
  has no particles and stops the distributed PCG), and so does
  `solve3D_pyramid_blobs` (its reconstruction asks for a cropped box larger
  than the native one).
- The derived reservations are wrong for several entries, as the samples
  above show (`single_atoms_stats` uses about 2.5 cores of 8 reserved,
  `mini_stream` about 1.8 of 1); Phase 1 corrects them from its profile.
- The class expansion of `flex_pca` seeds itself with fixed numbers in
  production code, so `flex_pca_blobs` part 3 does not change its draws with
  the seed offset; five unit and library testers seed outside
  `set_fixed_seed` (item 1).
- `species_discovery` writes several values twice under one name (reported,
  then held to a floor); the profile takes such a name as one metric with the
  floor's sense.

Files changed: `production/simple_test_exec.f90`,
`production/CMakeLists.txt` (a comment), `scripts/check_test_registry.py`,
`scripts/simple_args_generator.pl`, `src/utils/simple_test_utils.f90`,
`src/utils/simple_forked_process.f90`,
`src/utils/simple_forked_process_tester.f90`,
`src/main/exec/simple_test_exec_highlevel.f90`,
`src/main/ui/simple_test/simple_test_ui_class.f90`,
`src/main/ui/simple_test/simple_test_ui_highlevel.f90`, in
`src/main/commanders/test/`: `simple_commanders_test_class.f90`,
`simple_commanders_test_highlevel.f90` and its submodules `_addon`, `_cont`,
`_flex`, `_pcg`, `_pyramid`, `_stream`, `_workflow`,
`simple_commanders_test_single.f90`, `simple_commanders_test_stream.f90`,
`simple_test_gate.f90`, `simple_test_truth_metrics.f90`; new there:
`simple_commanders_test_highlevel_suite.f90`, `simple_test_registry.f90`,
`simple_test_dispatch.f90`, `simple_test_profile_stats.f90`,
`simple_test_identity.f90`, `simple_test_tools_tester.f90`,
`simple_test_dispatch_tester.f90`; `doc/policies/test_environment_policy.md`
and this note.

### Phase 1, 2026-10-16/17: the profile of the existing suite

**What ran.** Every high-level entry was profiled with the Phase 0 tools on the Release build with
tests of this phase (`keep/build_phase1`), at its CTest thread count, at seed offsets 0 to 19
(offset 0 draws what CTest draws), packed on the 24 cores by CPU reservation. All paths below are
under `/home/elmlundho/agent_runs/test_truth/scratch`. Each run's gate table (`metrics.tsv`), full
log, last 40 log lines and exit status (`exit_status.txt`, with wall time, peak CPU and peak
memory) are kept in `keep/profiles/phase1/<entry>/<entry>_t<threads>_s<offset>/`, with the
entry's profile table `profile_<entry>_t<threads>.tsv` beside them. The bulky outputs were deleted
as each run ended.

- **Decisions.**
  - `test=profile` learnt a comma-separated list of entries (file-table rows added), so that one
    packed run could profile every entry except the three that fail deterministically on the base,
    which the run's instructions exclude from the profile; those were run twice instead.
  - The profile was split into groups by disk use. Started as one run, it began all 20
    `simulated_workflow_6vxx` runs at once (one core each, about 3 GB each), and the free disk fell
    to 40 GB. I stopped it with SIGTERM: it stopped every child and its workers and left no
    process, as its contract says (`logs/profile_main_interrupted.log`). The workflows then ran
    at a budget of 8 cores (at most 8 at once) beside the other entries at 14, later 16 cores. That
    second group was stopped and restarted once for packing (only one 8-core entry fitted in 14
    cores). For a while the groups reserved 26 cores together, a slight oversubscription that
    lowers the sampled CPU use a little.
  - A mistake while sorting the results overwrote the three 20-run tables of the one-thread
    workflows with one-run tables. They were rebuilt from the kept gate tables by
    `scratch/regen_profile.sh` (the same statistics, in awk) and equal the numbers printed before
    the mistake (`TABLE_NOTE.txt` beside them).
  - The profile tool took the floor of a name written both as a report and against a floor
    (`species_discovery` does this) as the median with the report's placeholder. Fixed in
    `simple_test_profile_stats.f90` with a unit check; the species table was rebuilt.
- **Identity at offset 0.** Two runs of each one-thread entry at seed offset 0 are byte-identical
  over its declared outputs (`test=identity`): `simulated_workflow_1jxy` 205 files,
  `simulated_workflow_6vxx` 158, `single_workflow_wurtzite` 133, `single_workflow_fcc` 71 (up to
  its failure). Reports `keep/profiles/phase1/identity_<entry>.txt`.
- **Peak CPU reservations, confirmed and corrected.** The registry now holds the largest CPU use
  sampled in the profile, less 0.1 for sampling jitter, rounded up (the earlier value in
  brackets): `mini_stream_6vxx` 2 (1), `mini_stream_1jxy` 3 (1), `simulated_workflow_6vxx` 2 (1),
  `simulated_workflow_1jxy` 3 (2), `single_workflow_fcc` and `_wurtzite` 2 (1), `pcg_recon` 8,
  `simulate_particles` 8, `solve3D_addon` 11 (8), `cont_refine3D_1jxy` 10 (8),
  `single_atoms_stats` 3 (8), `_cdsew` 4 (8), `_cdsez` and `_cdser` 3 (8), `species_discovery` 8,
  `stream_preproc` 2 (1), `flex_pca_blobs` 9 (8), `solve3D_pyramid_blobs` 9 (8). The one-thread
  workflows peak above one core because the jobs they submit with one part briefly run beside the
  master; the 8-thread add-on and refinement peak above 8 for the same reason. Peak memory: at
  most 8.8 GiB (`cont_refine3D_1jxy`). The CTest processor counts follow in Phase 5.
- **Outcomes over the profile** (runs passed of runs made): `mini_stream_6vxx` 3/20,
  `mini_stream_1jxy` 18/20, `simulated_workflow_6vxx` 16/20, `simulated_workflow_1jxy` 17/20,
  `single_workflow_wurtzite` 20/20, `single_workflow_fcc` 0/2, `pcg_recon` 20/20,
  `simulate_particles` 20/20, `solve3D_addon` 10/20, `cont_refine3D_1jxy` 20/20,
  `single_atoms_stats` and its three variants 80/80, `species_discovery` 11/20, `stream_preproc`
  17/20, `flex_pca_blobs` 0/2, `solve3D_pyramid_blobs` 0/2. Ten of the eighteen entries fail on
  at least one draw at their CTest settings, three of them on every draw: today's nightly result
  is partly a matter of luck.

**The test cards.** Each card states the claim, what the test actually checks, which path the logs
show and what the test asserts, the profile of its gate values against their floors (median, the
robust spread s = 1.4826 times the median absolute deviation, the worst run, margins in units of
s), whether the floor is admissible by the plan's floor rule (at most the worst run minus s and at
most the median minus 3 s; the third condition, the negative control, cannot be judged where no
control exists), the negative control, determinism and cost. Wall times are from a loaded machine.

#### Card: `mini_stream_6vxx` and `mini_stream_1jxy`

- **Claim.** The mini stream (micrograph import, CTF estimation, picking, particle extraction, 2D
  classification by `solve2D` and shape ranking of the class averages) works on six synthetic
  micrographs of eight particles each (48 placed) of the 6VXX spike or the 1JYX
  beta-galactosidase.
- **What it checks.** Only the defocus is compared with the truth (each micrograph within
  0.10 µm). Picking is checked by count only (picked/48 between 0.40 and 1.50), because the
  simulator does not save the particle positions; extraction, classes and shape ranking only for
  consistency (a run with no ranked class passes). The micrographs are the simulator's
  motion-free average, so motion correction is not exercised.
- **Path assertions.** None in the test. The log shows the path (`CTF_ESTIMATE: 1 PART(S)
  COMPLETED`, `EXTRACT: 1 PART(S) COMPLETED`, the `solve2D` refinement modes,
  `SIMPLE_SHAPE_RANK_CAVGS NORMAL STOP`).
- **Profile** (20 seed offsets at 8 threads, `keep/profiles/phase1/rest/`):
  - 6vxx fails 17 of 20 runs, all on the pick ratio: median 0.354, robust spread s 0.062, worst
    0.25 against the floor 0.40. The floor lies above the median: the gate is mostly failing, and
    its value was set from one run (22/48). One run also failed a defocus row (0.118 µm on
    micrograph 4 against 0.10, where the median of that row is 0.002 µm: a rare outlier).
  - 1jxy fails 2 of 20, on the pick ratio: median 0.458, s 0.031, worst 0.354; the floor sits 1.9
    s below the median, inside the distribution.
  - Defocus rows: median errors 0.001-0.003 µm against 0.10, margins of 20 to 140 s below the
    ceiling, apart from the 6vxx outlier above.
  - The ceiling of the pick ratio (1.50) is never approached (worst 0.65).
- **Floor rule.** The pick-ratio floor is not admissible for either suite (it does not lie a
  spread below the worst run nor three spreads below the median). The defocus ceiling is, except
  for its rare outliers.
- **Negative control.** None. Cheap: the defocus of each micrograph against the truth of its
  neighbour (the nominal defocus steps are 0.25 µm, above the limit), which must fail every row.
  Picking cannot be controlled until the simulator saves the positions: a picker that returns 19
  random boxes passes today.
- **Determinism and cost.** The test hands one thread and one part to every commander, so its
  CTest thread count of 8 is unused; not reproducible across runs with threads (the 2D search).
  Wall time 83-90 s (median) under a loaded machine; peak CPU 1.3-2.4 cores against a reservation
  of 1.

#### Card: `simulated_workflow_6vxx` and `simulated_workflow_1jxy`

- **Claim.** Simulated movies of the 6VXX spike (point group C3) or 1JYX beta-galactosidase (D2)
  go through import, motion correction, CTF estimation, picking, extraction, `solve2D` and
  `solve3D`, and the final map matches the map made from the atomic model; 1jxy also runs `solve3D`
  distributed over two worker processes from the same project.
- **What it checks.** Only the final map(s): docked in both hands onto the truth, whole-volume
  correlation at least 0.80 and FSC 0.143 at most 40 Å. Defocus and frame drift are not compared
  with the truth, although the simulator writes both; particle positions are not saved, so picking
  is not checked; orientations are not compared, although the policy says they are.
- **Path assertions.** None. The log shows every step (`MOTION_CORRECT: 1 PART(S) COMPLETED`, the
  last `solve3D` stage `STAGE 8 WITH LP 6.0 A`, and for 1jxy `REFINE3D: 2 PART(S) COMPLETED`), but
  the test does not assert that the distributed step used two workers.
- **Profile** (20 seed offsets at one thread, `keep/profiles/phase1/workflows/`):
  - 6vxx fails 4 of 20 runs, 1jxy 3 of 20 (in each of them both `solve3D` steps fail, though
    with different values). The map correlation is bimodal: the runs that find the
    structure give 0.84-0.90 (median 0.871 for 6vxx, 0.887 for 1jxy, robust spread s 0.02-0.04),
    the others 0.54-0.80. The floor 0.80 lies 1.9 s (6vxx) and 4.7 s (1jxy) below the median, and
    the failed runs are far below it: these are runs in which `solve3D` did not find the map, not
    values scattered around the floor.
  - The FSC ceiling of 40 Å is far from what a right map gives (4.2-6.6 Å in the runs whose
    correlation passes) and lets wrong maps through: maps that fail the correlation reach FSC
    values of 6.4, 20.8 and 35.7 Å and pass this row; it failed only at 41.6 and 124.8 Å. The
    correlation alone decides.
  - Extracted particles vary from 45 to 84 with the seed (of 100 placed): the workflow rests on
    few particles, as the plan says.
- **Floor rule.** The correlation floor is not admissible: the worst runs are below it. This is a
  case where the right algorithm does not succeed on every draw (a robustness question for
  `solve3D` on so few particles), not one where a different floor would help; the plan's Phase 5
  redesigns it.
- **Identity.** Both are byte-reproducible at one thread: two runs at seed offset 0 give the same
  declared outputs (1jxy: 205 files, both `solve3D` steps; 6vxx: see below).
- **Negative control.** None. Cheap: the final map against the other suite's truth, or against a
  phase-randomised truth beyond 20 Å, which shows whether a wrong but globular map clears 0.80;
  the defocus and drift against the truth files already written.
- **Determinism and cost.** One thread; 6vxx 821 s (median, under load), 1jxy 1 272 s. Peak CPU
  1.3-1.4 cores (6vxx: its jobs with one part briefly run beside the master) and 2.4-2.6 (1jxy,
  two workers), against reservations of 1 and 2.

#### Card: `single_workflow_wurtzite` and `single_workflow_fcc`

- **Claim.** A simulated nanoparticle (face-centred cubic platinum, or wurtzite cadmium selenide)
  goes through reprojection along a trajectory, noise, denoising, 2D analysis and `autorefine3D_nano`,
  and the final map matches the generating model; continuous pose refinement (`refine=cont`) and
  the pose polish do not worsen the poses against the trajectory truth.
- **What it checks.** The map: band correlation to 5 Å at least 0.90, FSC 0.143 at most 5 Å. The
  poses: only against another run (each pose error at most the previous one plus 0.5°).
- **Profile of wurtzite** (20 seeds, one thread): passes 20 of 20. Map correlation median 0.980,
  s 0.0015, worst 0.974 (50 s above the floor); FSC median 0.95 Å, worst 3.0 Å against 5 Å.
  **Finding:** every pair pose error, in every run and every stage, is 42.7-45.6° (median 44.3°),
  the level of random poses, while the map is excellent. The pair error measures relative rotation
  angles and is not symmetry-aware, so either the poses are right only up to the lattice's
  symmetry, or the particle order of the project does not match the truth table. Either way the
  two pose gates compare two chance-level numbers and pass whatever the poses are.
- **fcc** (run twice at offset 0): fails both times, byte-identically up to the failure: in the
  second iteration of `autorefine3D_nano` the matching low-pass is clamped from 1.5 Å to the FSC
  0.143 resolution (4.77 Å), the Fourier index range becomes empty (18 to 12) and the polar setup
  stops ("invalid kfromto"). Not profiled further.
- **Identity.** wurtzite: two offset-0 runs identical over 133 files; fcc: 71 files up to the
  failure.
- **Negative control.** None. Cheap: shuffled truth poses, which would leave the pose errors near
  45° and show that the pose gates cannot fail; the fcc truth against the wurtzite map.
- **Determinism and cost.** One thread; wurtzite 655 s (median under load), fcc 86-91 s to its
  failure; peak CPU 1.3 against 1.

#### Card: `pcg_recon`

- **Claim.** The preconditioned conjugate-gradient (PCG) reconstruction operator and solver are
  right: adjoint identities, symmetry and positive definiteness of the normal operator, recovery of
  a known phantom, agreement of the kernel and matrix-free operators, streaming, symmetry
  replication, cropping, the maximum-likelihood prior, observation parity and support semantics,
  in fourteen fail-fast stages (its description says nine).
- **What it checks.** Every stage against an analytic identity or a second code path, plus the
  phantom recovery (an inverse crime: the same operator simulates and solves) and deapodization
  against envelope-free data. It is an in-process numerical test of a library operator, registered
  as high level.
- **Path assertions.** In-process; `every_stage_passed` asserts that all stages ran.
- **Profile** (20 seeds): passes 20 of 20. The tolerances sit far from the values: adjoint errors
  1e-9 against 1e-5, kernel interior error 0.030 against 0.05 (26 s), the kernel's change under a
  changed CTF 0.44 against a floor of 1e-4. Most rows are identical in every run (s = 0): the seed
  changes only the random probes.
- **Floor rule.** Admissible wherever the value varies; the constant rows are exact properties
  checked against round-off tolerances.
- **Negative control.** Partly: the kernel must change when the CTF changes, and deapodization on
  must beat off. Missing: the recovered map against a wrong phantom.
- **Determinism and cost.** About 2 s; up to 8 cores.

#### Card: `simulate_particles`

- **Claim.** `reproject` and `simulate_particles` on the embedded 6VXX map produce the stacks
  and orientation tables asked for, with the requested CTF, defocus spread, astigmatism and shifts.
- **What it checks.** Metadata only: counts, box, sampling, finite non-constant images,
  orientation spread, and that the drawn defocus, astigmatism and shifts lie within the requested
  ranges. The images are never compared with the projections at their recorded poses, and the
  applied signal-to-noise ratio and CTF are not verified.
- **Path assertions.** Outputs and counts.
- **Profile** (20 seeds): passes 20 of 20. The bound rows are the largest of 200 uniform draws
  against the bound of the draw, so they sit just inside it by construction (defocus 0.1992 µm
  against 0.2001, shift 1.995 against 2.0001 pixels); the floor rule does not apply to them: they
  are range checks.
- **Negative control.** None. Cheap and needed: each image against the reprojection at its
  recorded pose and at a shuffled pose (Phase 5 adds the image check).
- **Determinism and cost.** 12 s; 6.7-7.6 cores of 8.

#### Card: `solve3D_addon`

- **Claim.** `solve3D` on a seeded selection of 2 000 simulated particles of a 6VXX map whose
  symmetry is broken by an added blob (the frozen base solution), then `solve3D_addon` on 2 500 and
  on 3 000 particles (two add-ons; the description says one), against the simulation truth and the
  base run: the new particles (the cohort) get poses as good as the frozen ones, and the union map
  is at least as good as the base map.
- **What it checks.** Provenance and integrity in depth (manifests, the carried sums, the frozen
  rows restored, sigma2 isolation; these are path assertions too: the first add-on reads the frozen
  particles once, the second reads none). Against the truth: the cohort's pair pose error at most
  15° and at most 5° above the frozen particles'; the union map's correlation with the truth at
  least 0.90; relative to the base: correlation loss at most 0.03, FSC at most 1 Å worse. The
  particles carry no shifts.
- **Profile** (20 seeds at 8 threads, `keep/profiles/phase1/solve3D_addon/`): **10 of 20 runs
  fail**, 9 on the cohort pose error (15.6-18.4° against 15) and 1 on the excess (5.02 against
  5). The cohort error has median 14.9°, robust spread s 2.0°: the floor sits at the median. The
  cause is the base: the frozen particles' own pair pose error ranges 4.6-16.6° (median 11.3°,
  only 38-51 % of pairs within 5°), although the base map correlates with the truth at 0.93-0.95
  in every run, and the cohort follows the frozen particles (excess median 3.2°, s 1.1°). The map
  rows pass with margins of 5 to 10 s: union correlation median 0.941 (worst 0.924, floor 0.90),
  FSC median 4.8 Å.
- **Floor rule.** The two pose floors are not admissible (the absolute one sits at the median, the
  excess one at the worst run); the map floors are. The absolute pose floor fails for a reason the
  add-on does not cause (a poorly posed base), so it does not test the add-on; Phase 5 decides
  between a floor relative to the base alone and a better-posed base.
- **Negative control.** None. Cheap: shuffled truth poses (cohort error near 40°), and the truth
  without its blob for the map rows.
- **Determinism and cost.** 8 threads, not reproducible run to run. 549 s (median, loaded
  machine); peak CPU 9.0 cores, up to 10.6 (its commanders' threads overlap with the reading
  threads), against the old reservation of 8.

#### Card: `cont_refine3D_1jxy`

- **Claim.** Continuous Cartesian pose refinement of 1 000 simulated 1JYX particles from poses
  perturbed by 15° and 2 pixels recovers the true poses: `refine3D refine=cont` under two
  objective functions and distributed over two workers, `refine3D_auto pose_cont=only`, the
  polish, and the continuation of a polar `refine3D_auto`.
- **What it checks.** Absolute, truth-based, symmetry-aware pose floors (median rotation error at
  most 1.5°, shift at most 0.2 pixels, at least 95 % of particles within the 8 Å basin), the
  distributed result against the shared-memory one, and route assertions (no polar model, no
  discrete pass assigned a projection, the polish refines exactly its sample). Most stages refine
  against the truth map itself, so they test local pose optimisation, not reference building; only
  `refine3D_auto` builds its own reference. The particles carry no shifts; the 2 pixels are put
  into the starting poses. That the distributed stage used two workers is not asserted.
- **Profile** (20 seeds): passes 20 of 20, every floor admissible by a wide margin: rotation
  errors median 0.94° (euclid) and 0.23° (cc), worst 1.01° against 1.5° (14.5 s); shifts 0.06-0.07
  pixels against 0.2; basin 0.984 worst against 0.95 (6.6 s); `refine3D_auto` 1.1° median, worst
  2.0° against 7.5°; the polish and continuation rows with margins of 1.9 to 300 s. The
  distributed run equals the shared-memory one exactly in every run (difference 0).
- **Floor rule.** Admissible throughout; the closest is the polish row (worst 1.9 s inside its
  relative floor).
- **Negative control.** Partly: the starting errors (15°, 2 pixels) are reported and must be
  beaten. Cheap: the floors against shuffled truth poses, and stage a with a wrong reference.
- **Determinism and cost.** 8 threads; 504 s; peak CPU 8.9-9.2 cores against 8; peak memory
  8.8 GiB, the largest of the suite.

#### Card: `single_atoms_stats`, `_cdsew`, `_cdsez`, `_cdser`

- **Claim.** Atom detection and statistics recover the atoms of a simulated nanoparticle
  (platinum; cadmium selenide in the wurtzite, zincblende and rocksalt lattices), and for the
  binary crystals the interior coordination and the fitted lattice constant.
- **What it checks.** Recall, precision and RMS position error against the model, the correlation
  of the detected-atom density with the simulated map, and for the binary crystals the lattice
  constant within 2 % and no interior atom with the wrong coordination; output files and their row
  counts.
- **Profile** (20 seeds): passes 80 of 80, and every value is the same in every run (s = 0): the
  simulated map is noise-free and does not depend on the seed, so the seed offset changes nothing.
  Recall 1.0 (rocksalt 0.970), precision 1.0, RMS error 0.013-0.015 Å against a ceiling of 0.75 Å,
  density correlation 0.98-0.9997 against 0.80, lattice error 0.4-0.5 % against 2 %.
- **Floor rule.** With s = 0 the rule cannot place a floor; the floors are far from the values and
  detect only gross breakage. The test is a deterministic regression check, not a correctness gate
  in the plan's sense. Matching is nearest-neighbour within 1 Å, not one-to-one.
- **Negative control.** None. Cheap: detected coordinates shifted by half a bond length (recall
  near 0); noise at a measured level (as `species_discovery` simulates) would give the floors a
  distribution.
- **Determinism and cost.** 9-19 s; peak CPU 2.2-3.1 cores against the old reservation of 8.

#### Card: `species_discovery`

- **Claim.** `detect_atoms` with a species list and with species discovery, on platinum,
  platinum-nickel (alloy and core-shell) and platinum-aluminium particles rendered with element
  kernels, surface B-factors and measured-spectrum noise, labels every atom with its element
  against the generating model, and leaves the existing products unchanged.
- **What it checks.** Recall and false atoms per run, the recall of the second element (at least
  0.90), half-map label agreement, class order, every label right, radial shells, product identity;
  negative controls built in (a pure particle must give one class; plain platinum detection must
  miss the aluminium, which proves the light case tests recovery).
- **Profile** (20 seeds): **9 of 20 runs fail**, all on one check: in the run on the pure platinum
  particle with the list platinum-nickel, the platinum class has fewer than two radial shells, so
  the radial-shell analysis writes its failure. The check depends on the draw. All recall rows
  pass: platinum recall at least 0.995 (12 s above 0.98); nickel recall median 0.985, s 0.022,
  worst 0.940 against 0.90 (1.8 s); aluminium recall median 0.985, worst 0.908 against 0.90 (0.4 s,
  inside the distribution).
- **Floor rule.** The second-element recall floors are not admissible for aluminium (worst less
  than one spread above the floor) and only just for nickel; the others are.
- **Determinism and cost.** 8 threads; 34 s; peak CPU 7.3-8.0 cores.

#### Card: `stream_preproc`

- **Claim.** Five simulated movies through the stream's preprocessing stage (motion correction
  and CTF estimation in its worker jobs) recover the simulated defocus, frame drift and image.
- **What it checks.** Against the truth: defocus within 0.10 µm, frame drift within 0.5 pixels,
  the integrated micrograph's correlation with the optimal average at least 0.50 and at least
  0.05 above its correlation with another movie's average (a negative control). Also the outputs
  and the number of micrographs.
- **Path assertions.** Outputs only; that the work ran as a worker job is not asserted (the log
  shows `qsys_name local`, `# MOVIES PROCESSED & IMPORTED : 5`).
- **Profile** (20 seeds, `keep/profiles/phase1/rest/`): 3 of 20 runs fail, all on a defocus row:
  0.164, 0.106 and 0.127 µm against 0.10, where the medians of the rows are 0.007-0.016 µm and s
  0.005-0.006 µm. The CTF fit sometimes lands on a wrong defocus: a heavy tail, not a broad
  distribution. Drift: medians 0.06-0.11 pixels, worst 0.36 against 0.5 (2.2 to 8.5 s below the
  ceiling). Image correlation: medians 0.85-0.92, worst 0.82, against 0.50; the margin over the
  wrong movie 0.84-0.91 against 0.05. (The failure rates of the later movies include the runs
  that stopped at an earlier movie.)
- **Floor rule.** The image floors are admissible by a wide margin. The drift ceiling is
  admissible. The defocus ceiling is not, because of the outliers: 15 % of runs fail.
- **Negative control.** Yes for the image (the wrong movie: 0.001-0.010). Missing and cheap for
  defocus (the neighbouring movie's truth) and drift (a zero-drift estimate, which would show how
  far the floor is from no correction at all).
- **Determinism and cost.** Seeded by the front door's process seed only; one one-thread worker
  job at a time. Wall time 92 s; peak CPU 1.2-1.3 against 1.

#### Card: `flex_pca_blobs` (failing; run twice)

- **Claim.** `flex_pca` separates the two states of a planted phantom in shared memory and with two
  local workers, in agreement with the truth and with each other; the per-state refinement of its
  states works; and the class-average state refinement (`cls_expansion`, then
  `refine3D_states_cavgs`) recovers the three states of a second phantom of 3 600 particles in 60
  exact views.
- **What it checks.** Parts 1 and 2 run on the true poses with the true mean as reference, so they
  test state separation, not search; the label accuracy floor is 0.60 against observed 0.998 and a
  chance level of 0.5. Map preferences are relative. Part 3: committed sub-averages, sub-average
  state accuracy at least 0.60, the median pair pose error of the sub-averages at most 15°.
- **Runs** (seed offsets 0 and 1, `keep/profiles/phase1/side_fail/`):
  - Offset 0 fails in part 3 before its gate: FLEX inside `refine3D_states_cavgs` publishes an
    empty state 2 and the distributed PCG stops ("state 2 has no particles and no previous maps in
    this directory to carry forward").
  - Offset 1 runs part 3 to its end, through the neighbourhood stages at 2 500 directions, and
    fails its gate: the sub-average state accuracy is 0.367 and the particle state accuracy 0.362,
    the chance level of three states, while the sub-averages are pure (0.95) and well posed
    (4.3°). So the failure depends on the draw, and when part 3 runs it does not separate the
    states.
  - Parts 1 and 2 pass in both.
- **Negative control.** None. Cheap: shuffled labels into the accuracy (expect about 0.5 in parts
  1 and 2, which the 0.60 floor fails by little), and a one-state phantom.
- **Determinism and cost.** Not reproducible (threads); `cls_expansion` seeds itself with fixed
  numbers in production code, so that part does not move with the seed offset. 470-625 s under
  load; 8.4 cores at peak against 8.

#### Card: `solve3D_pyramid_blobs` (failing; run twice)

- **Claim.** `solve3D_pyramid` with three states, end to end on the three-state phantom of part 3
  above, assigns the planted states (purity) and builds references and refined maps that dock onto
  the truths.
- **What it checks.** Purity at least 0.45 (chance 0.33); each reference docked onto its closest
  truth and each refined map onto the majority truth of its particles, at least 0.70, with no
  requirement that the truths be distinct; a non-finite correlation is not guarded. Floors
  provisional.
- **Runs** (offsets 0 and 1): both stop before any metric, in the pyramid's `solve3D_cavgs` step,
  which plans a cropped box of 88 pixels on the 64-pixel phantom ("reconstruct3D box_crop cannot
  exceed the native box").
- **Negative control.** None. Cheap: random state labels for the purity; a one-to-one matching of
  references to truths, so that three copies of one truth fail.
- **Cost.** 65 s to the failure.

**Findings across the suite** (for Phase 5; no test was changed in this phase).

- Floors inside or above their distribution: the pick ratio of `mini_stream` (6vxx fails 85 % of
  runs, 1jxy 10 %), the map correlation of `simulated_workflow` (bimodal: the runs that do not find
  the map), the pose floors of `solve3D_addon` (45 % and 5 % of runs fail, because the base run's
  poses vary), the defocus ceiling of `stream_preproc` (15 %, rare wrong CTF fits), the aluminium
  recall of `species_discovery`, and a draw-dependent radial-shell check in it (45 %).
- Gates that cannot fail for the right reason: the pose gates of `single_workflow` (both sides at
  the chance level of about 44°), the FSC ceiling of 40 Å in `simulated_workflow` (wrong maps of
  35.7 Å pass it), the metadata bounds of `simulate_particles`, and the noise-free floors of
  `single_atoms_stats` (identical in every run).
- Multi-state: `flex_pca_blobs` part 3 fails with the draw in one of two ways: an empty state stops
  the distributed PCG, or `refine3D_states_cavgs` runs through its 2 500-direction neighbourhood
  stages and reaches only the chance level of state accuracy (0.37 for three states).
  `solve3D_pyramid_blobs` stops on a cropped box larger than the native one at every draw.
- `single_workflow_fcc` stops on an empty Fourier range in `autorefine3D_nano` at every run.
- No test asserts that its distributed step used several workers, although every log shows it.

**Exit criteria.**

- Every entry profiled (15 entries, 20 seeds each, 300 runs) or recorded as failing
  deterministically (`single_workflow_fcc`, `flex_pca_blobs`, `solve3D_pyramid_blobs`, run twice
  each; `flex_pca_blobs` turned out to fail in two ways, as its card records).
- Every card drafted above.
- Gates after this phase's changes: the fast gate passes (Debug, 15 suites, 7.6 s,
  `logs/compile_debug_full_p1.log`), the platform entry `forked_process` passes, the registry
  check and `check_descr.py` pass, and `test tools` passes with its new check (105 checks).

Files changed in this phase: `src/main/commanders/test/simple_commanders_test_highlevel_suite.f90`
(entry lists), `src/main/commanders/test/simple_test_registry.f90` (reservations),
`src/main/commanders/test/simple_test_profile_stats.f90` (the floor of reported and gated names),
`src/main/commanders/test/simple_test_tools_tester.f90` (its check),
`src/main/ui/simple_test/simple_test_ui_highlevel.f90` (the help of `target=`),
`doc/policies/test_environment_policy.md` (the entry list), and this note.

### Phase 2, 2026-10-17: `simulate_particles` with several states, the fixture, the metrics

All paths below are under `/home/elmlundho/agent_runs/test_truth/scratch`. Build: Release with
tests in `keep/build_phase2`; the Debug build of `./compile_debug.sh` for the fast gate.

**The simulator.** `simulate_particles` takes `nstates=K` with `vol1` to `volK` and projects every
particle from the volume of its state; the K padded volumes are held together. A missing `volK`
or a volume of another box is refused.

- The states come from the state field of `oritab` when it has one (in a text orientation
  document, a `state=` key on a line; a project file always has states), and a state outside 1 to
  K is refused. Otherwise they are drawn with the population fractions of the new key
  `state_pops` (equal by default) in exact counts: each fraction times the number of particles,
  rounded, the remainder given to (or taken from) the states in order of decreasing fraction, the
  labels shuffled with the seeded generator after the poses and CTF are drawn.
- **Decision on the key.** `parameters` has no list of reals; its lists are comma-separated
  character keys that the consumer parses (`states`, `clustinds`), so `state_pops` follows that
  pattern (a field in `simple_parameters.f90`, its binding in `simple_parameters_parse.f90`), not
  a text file.
- With `nstates` absent or 1, every particle has state 1 as before. **Identity holds:** the
  `simulate_particles` entry at one thread on the Phase 1 build and on this build gives
  byte-identical `simulated_particles.mrc`, `simulated_oris.txt`, `reprojs.mrcs` and
  `reproject_oris.txt` (`logs/p2_simpart_p1.log`, `logs/p2_simpart_p2.log`, compared with `cmp`
  after the last change of the phase).
- Interface text (`simple_ui_sim.f90`) and `doc/user_guide/simulate_particles.md` describe the
  mixture, the rule for the counts and the memory cost.

**The fixture** (`simple_test_state_fixture.f90`). One `simulate_particles` call with all states
makes each dataset; the truth (state, pose, shift, CTF of every particle) is the simulator's
orientation document `truth.txt`; the project of the mode is written into the caller's working
directory.

- Families: blob phantoms (box 64 at 2 Å, the common body and per-state lobe of the existing
  `write_phantom_truths`, now exported by the FLEX tester, file-table row added; also an
  `identical` kind for the negative control) and beta-galactosidase 1JYX without its hetero atoms
  (box 96 at 2.2 Å): compositional states (the tetramer, the tetramer without chain D, the A-B
  dimer) and conformational states (the tetramer, and the tetramer with chain D rotated rigidly by
  a stated angle about the axis through its centre parallel to the tetramer's long axis).
- **One frame.** The D2 axes of 1JYX are found from its chains: a 2-fold axis through the centre
  carries the midpoint of a point and its image, so it runs along the sum of two related chain
  centres (C-alpha atoms of the residues all four chains hold). The coordinates are rotated so
  that the axes are SIMPLE's `d2` frame (z relating chains A and B, x relating A and C), translated
  by one shared vector that puts the tetramer's density centre on the box centre, and every state
  is cut from these coordinates and rendered by `pdb2mrc` without centring. The declared groups:
  tetramer `d2`, tetramer less a chain `c1`, dimer `c2` (its 2-fold is z), rotated chain `c1`.
- **Point-group check:** every truth correlates with itself under every operator of its group
  above 0.99 and not under the extra operator of the next larger group of its family (c1 to c2,
  c2 to c4, d2 to d4, about z). Values: tetramer 0.995 under D2, 0.49 under the 90° turn; less a
  chain 0.64 under C2; dimer 0.996 under C2, 0.37 under the 90° turn; rotated-chain state 0.956
  under C2 (a 10° rotation); blob truths 0.28-0.29 under C2.
- **Modes:** search (random poses, shifts within ±2.5 pixels; the project carries no poses and
  every particle in state 1), known (the same dataset; the project carries the true poses and
  shifts), views (a number of spiral directions each repeated, no shift; the project carries the
  true poses and the 2D class equal to the view).
- **The cache key rule.** A dataset is cached in
  `<SIMPLE_PATH>/test_runs/fixture_cache/<family>_<hash>`, where the hash is the 64-bit FNV-1a
  hash (16 hexadecimal digits) of a key text holding a fixture version number, the family and
  kind, the number of states, particles and (views mode) views, SNR, defocus, defocus spread,
  astigmatism, shift range (0 for views), the rotation angle, the seed plus the seed offset, the box
  and the population fractions. The search and known modes share a dataset. The directory holds
  the key as `key.txt` and a `COMPLETE` marker; it is built in a temporary directory renamed into
  place, so concurrent runs cannot read half a dataset. The fixture version is raised when the
  fixture or the simulator changes its output. After a build or a cache hit the process seed is
  reset to the front door's, so a test draws the same either way.
- Low-passed truths for starting volumes: `write_lowpassed`.

**The metrics** (`simple_test_truth_metrics.f90`), each with a unit test (`truth metrics` in
`unit_numerics`, 40 checks, 0.07 s Release, 0.6 s Debug; `truth map docking` in
`lib_reconstruction`, 0.4 s):

- the matching of estimated onto true states by exhaustive search (up to eight states), invariant
  under a global renaming; the adjusted Rand index against values from an independent pair count
  (0.2424, -0.5, 0.3911, 1 for a renamed identical partition, 0 for a collapsed one); the accuracy
  under the mapping; its null distribution from 1 000 shuffles of the true labels with the counts
  kept (for populations of 50, 30 and 20 % the median 0.383 and 99.9th percentile 0.426 of a
  Python emulation are reproduced) and the collapsed predictor (0.50, adjusted Rand index 0);
  population errors;
- the map mapping from the docked correlations, which must agree with the label mapping and does
  not when correct labels come with swapped maps (a non-finite correlation counts as no
  agreement);
- discrimination inside the mask where the truths differ (by more than a fraction of their
  maximum, dilated): each map's margin over the other truths, and the correlation of map
  differences with truth differences; the truths themselves pass (margins about 0.18 on the
  blobs), copies of one map and of their mean fail (margins not all positive, differences without
  a finite correlation);
- **the frame and hand rule**, found by projection and pinned by the test: a docked map is
  `rotvol` of the estimate (or of its x-mirror) by a rotation A and shift s; a particle's pose r
  and shift sh carry into the truth's frame as r·Aᵀ and sh − (r·s)₁,₂, after, for the mirrored
  hand, r → Mz·r·Mx and, for an even box, the one-voxel move of `mirror('x')`. The test checks
  that each image of the estimate equals the image of the truth at the carried pose (correlation
  above 0.995) and that a wrong composition does not; the library test does the same through a
  real `dock_both_hands`, which now also returns its rotation, shift and hand;
- the pose error under a declared point group only: min over its operators g of the angle between
  the estimated pose and r_true·g (the composition under which a symmetric map projects the same,
  pinned by projection); poses equal up to D2 score zero, and a global 30° frame error, which a
  cluster learnt from the data would accept, scores above 25°;
- a non-finite value fails its gate row.

**`cont_refine3D_1jxy`** now measures its symmetric pose errors under the declared D2 of 1JYX (the
identity and the 180° turns about the axes `betagal_frame` finds, in the truth map's frame, which
the refinements share); `sym_pose_error`, which learnt rotation clusters from the errors, is
deleted. Run on this build (`logs/p2_cont.log`): passes, seed 15.0°, polar 1.39°, continued 0.20°,
polish 0.41° against 0.84° without, `refine3D_auto` with polish 0.18°. A check of the truth map's
own D2 invariance was tried and dropped: that map keeps the waters and is centred on its geometric
centre, not its symmetry centre (0.85; `logs/p2_cont_first.log`).

**The new tests of the mixture**, in the `simulate_particles` entry: two mixtures of three blob
states (populations 50/30/20 % of 90 particles: counts 45, 27, 18; 0.333/0.333/0.334 of 100: counts
33, 33, 34, the remainder rule), every image closest to the projection of its own state's volume
at its true pose (no CTF, no shift, SNR 2); and the fixture checks above (blob modes, each built
twice with the second from the cache, the projects' contents, the betagal point groups). The entry
passes at one thread in 37 s, building its datasets (`logs/p2_simpart_p2.log`), and at its CTest
setting of 8 threads in 10 s, reading them from the cache (`logs/p2_simpart_new.log`).

**Review.** A separate reader reviewed the phase's diff; fixed before the final runs: one state
again sets state 1 before any `oritab` is consulted; the fixture's project moved out of the shared
cache; the process seed reset after a build; `match_maps` initialised and robust to non-finite
correlations; missing `volK` and unequal boxes refused; guards in `betagal_frame` and
`discrimination`; heavy test locals made allocatable; dead code removed.

**Exit criteria.**

- The gates pass: the fast gate (Debug, 15 suites, 9.7 s, `logs/compile_debug_full_p2.log`), the
  touched suites (`unit_numerics` also in reverse order, `lib_reconstruction` truth map docking,
  the `simulate_particles` and `cont_refine3D_1jxy` entries), the registry check and
  `check_descr.py`.
- Identity holds for one state (above).
- The new tests pass (above).

Files changed in this phase: `src/main/commanders/simple/simple_commanders_sim.f90`,
`src/main/params/simple_parameters.f90`, `src/main/params/simple_parameters_parse.f90`,
`src/main/ui/simple/simple_ui_sim.f90`, `doc/user_guide/simulate_particles.md`,
`src/main/flex/run/simple_flex_pca_application_tester.f90` (exports `write_phantom_truths`), in
`src/main/commanders/test/`: `simple_test_truth_metrics.f90`, new
`simple_test_truth_metrics_tester.f90` and `simple_test_state_fixture.f90`,
`simple_commanders_test_highlevel_stream.f90`, `simple_commanders_test_highlevel_cont.f90`,
`simple_commanders_test_class.f90`; `src/main/ui/simple_test/simple_test_ui_class.f90`,
`src/main/ui/simple_test/simple_test_ui_highlevel.f90`, `doc/policies/test_environment_policy.md`
and this note.

### Phase 3, 2026-10-17: the multi-state assignment test

All paths below are under `/home/elmlundho/agent_runs/test_truth/scratch`. Build: Release with
tests in `keep/build_phase3`; the Debug build of `./compile_debug.sh` for the fast gate.

**The new entry `multistate_assignment`** (`simple_commanders_test_highlevel_multistate.f90`,
registered in `production/CMakeLists.txt` with 8 threads and a 3 600 s timeout, its registry row
reserving 11 cores; the process budget rises from 39 to 40 for it, as this phase of the plan
decides). It runs `refine3D nstates=3` twice, with `refine=prob` and with `refine=prob_neigh`
(`nspace=2500`, `nspace_sub=500`), on the blob family of the state fixture in its search mode:
1 500 particles, populations 50, 30 and 20 %, SNR 2 with CTF, random poses and shifts of up to
2.5 pixels; starting volumes the three truths low-passed to 20 Å; `lp=8`, `trs=4`, 6 iterations,
`pgrp=c1`, in shared memory. `suite=prob`, `prob_neigh` or `identical` runs one case;
`suite=identical` is the identical-states control, run by hand.

- **The initial states are the test's.** A finding: `refine3D` in shared memory (no `nparts`)
  with `nstates` above 1, on a project whose particles are all in state 1, never seeds the states.
  Only the distributed path calls `seed_multistate_startup_labels`, and the probabilistic table
  searches only states with more than 5 particles, so the run silently refines one state. The test
  therefore draws uniform random initial states with the production `gen_labelling(..., 'uniform')`
  the distributed path uses, before `refine3D` runs.
- **The frame.** The starting volumes fix frame and hand, so maps, poses and shifts are compared
  in the truths' frame without docking (the docking transformation is the identity).
- **One mapping** of estimated onto true states (`match_states`) serves the populations, the maps,
  the poses and the discrimination. The map mapping on its own uses the correlations inside the
  difference mask of the truths (`discrimination` now also returns them): whole-volume
  correlations of the blob states, which differ only in a small lobe, barely discriminate.
- **Rows:** the adjusted Rand index; the accuracy, its null (median and 99.9th percentile of
  1 000 shuffles) and the collapsed predictor; the largest population error; the shuffled-label
  control (its adjusted Rand index; a check that its accuracy lies within the null and that it
  fails the label floors); every state has a map; the map mapping's agreement with the label
  mapping; per state the correlation and FSC 0.143 against its truth; the discrimination (smallest
  margin and difference correlation); the duplicated-maps control (copies of one map, and of their
  mean, must not count as distinct); the median rotation and shift errors of the correctly
  assigned particles. The log prints the confusion matrix and the masked map correlations.

**Signal and the merging of states.** The first runs (SNR 1) and four diagnostics (one iteration;
SNR 2; 15 iterations; SNR 3 with 15 iterations; `logs/p3_x_*.log`) show the same pattern: the
estimated state started from truth 2 shrinks to 30-185 particles and the 30 % state is absorbed
into the state of truth 3 (confusion such as 740/7/4, 9/40/7, 1/403/289), whatever the signal or
the number of iterations. The adjusted Rand index stays near 0.7, not near 1. A finding for
Phase 5: `refine3D` with three blob states, from the true references and random labels, merges
two of them in most runs. SNR 2 (the plan's range) was kept because it gives positive
discrimination margins where the states do separate.

#### Card: `multistate_assignment`

- **Claim.** `refine3D` with three states (`refine=prob`, `refine=prob_neigh`), from the true
  references low-passed to 20 Å and random labels, assigns the particles of a blob mixture far
  above chance, and finds the poses of the correctly assigned particles.
- **Path assertions.** Both refinements run with three states and leave a map for every state;
  the log carries the refine mode and the populations per iteration.
- **Profile** (`keep/profiles/phase3/multistate_assignment/`, 20 seed offsets at 1 and at 8
  threads, 40 executions, 2 089 s packed; tables `profile_multistate_assignment_t1.tsv` and
  `_t8.tsv`). Medians at 1 / 8 threads, robust spread s, worst run:
  - `prob`: adjusted Rand index 0.658 / 0.653 (s 0.039 / 0.058, worst 0.16 / 0.21), accuracy
    0.741 / 0.740 (worst 0.574 / 0.595), largest population error 0.172 / 0.172 (worst 0.30);
  - `prob_neigh`: adjusted Rand index 0.663 / 0.637 (s 0.034 / 0.047, worst 0.497 / 0.480),
    accuracy 0.760 / 0.745 (s 0.034 / 0.029, worst 0.701 / 0.689), largest population error
    0.167 / 0.165 (worst 0.235 / 0.268); runs where the states separate reach 0.90-0.93;
  - the null's median 0.40 and 99.9th percentile 0.43 (at most 0.48); the collapsed predictor 0.50;
  - median rotation error 4.3-4.7° (worst 5.56°), median shift error 0.22-0.27 px (worst 0.34);
  - state 1's map correlates 0.975 with its truth (FSC 0.143 at 8-8.5 Å); the other two vary with
    the merging (0.53-0.97);
  - the map mapping disagrees with the label mapping in 17 (`prob`) and 16 (`prob_neigh`) of 40
    runs, every one a run with two true states merged; every run with an adjusted Rand index of
    0.72 or more agrees.
- **Floors** (written in the source with their origin): `prob_neigh` adjusted Rand index at
  least 0.30 (admissible 0.146-0.434), accuracy at least 0.63 and above the null's 99.9th
  percentile (admissible 0.602-0.658), largest population error at most 0.30 (admissible
  0.295-0.386); both modes median rotation error at most 6.0° (admissible from 5.93°) and median
  shift error at most 0.40 px (admissible from 0.380). **Recorded as a profile, with the reason:**
  the labels and population error of `prob` (its worst runs lose two states; no floor clears the
  controls by 3 s), the map-mapping agreement (35-50 % of runs fail it, by the merging), the
  discrimination (spread 0.2-0.65, worst runs below the duplicated-maps control) and the map fits.
  These are reported; a non-finite reported value still fails its row.
- **Negative controls.** Shuffled labels, every run: adjusted Rand index -0.004 to 0.006,
  accuracy inside the null in all 40 runs, and below both label floors. Duplicated maps, every
  run: copies of one map (smallest margin about -0.46, difference correlation not finite) and of
  their mean (about -0.1 to -0.17) never count as distinct. Identical states (`suite=identical`,
  run at calibration): every particle ends in one state (adjusted Rand index -0.002 to 0, accuracy
  0.50 = the collapsed predictor), below the floor; the identical truths leave no difference mask,
  so the discrimination has nothing to judge there.
- **Determinism and cost.** Reproducible at one thread: two runs at offset 0 give identical
  `metrics.tsv`, equal to the profile's (`logs/p3_det_a.log`, `_b.log`, 1 142-1 169 s each).
  At 8 threads it is not (the values move in the third or fourth digit). 91-106 s at 8 threads
  alone, peak 10.15 cores, 0.92 GiB; CTest 96 s.

**Review.** A separate reader reviewed the phase's diff; applied: an unknown `suite` refused
before any directory is made, the working directory's changes checked and `CWD_GLOB` kept, the
mean control divided by the number of states, the floor comment quoting the binding profile of
each floor, the identical-states control's limit stated, `nthr` optional so that `suite=list`
runs alone.

**Exit criteria.**

- The gates pass: the fast gate (Debug, 15 suites, 9.7 s, `logs/compile_debug_full_p3.log`), the
  `truth metrics` sub-suite of `unit_numerics` (42 checks, with the new masked-correlation output
  of `discrimination`), the new entry through CTest (96 s, `logs/p3_ctest_multistate.log`) and by
  hand in each sub-suite, the registry check and `check_descr.py`.
- Floors: admissible for the labels and populations of `prob_neigh` and for the poses and shifts
  of both modes; the labels of `prob`, the map-mapping agreement, the discrimination and the map
  fits are recorded as a profile, with the reason (the merging of two states) in the card and the
  source.

Files changed in this phase: in `src/main/commanders/test/`: new
`simple_commanders_test_highlevel_multistate.f90`, `simple_commanders_test_highlevel.f90` (the
commander type), `simple_test_registry.f90` (the entry), `simple_test_truth_metrics.f90` and
`simple_test_truth_metrics_tester.f90` (the masked correlations of `discrimination`);
`src/main/exec/simple_test_exec_highlevel.f90`, `src/main/ui/simple_test/simple_test_ui_highlevel.f90`,
`production/CMakeLists.txt` (the entry, the budget), `doc/policies/test_environment_policy.md`
(the entry, the counts; file-table row added) and this note.

### Phase 4, 2026-10-17/18: the multi-state workflow tests

All paths below are under `/home/elmlundho/agent_runs/test_truth/scratch`. Builds: Release with
tests in `keep/build_phase4` (the one the profiles ran), `keep/build_phase4b` (the final code, for
the confirmation runs, with the same fixture cache), the paired variant built from a copy of the
tree, and the Debug build of `./compile_debug.sh` for the fast gate. After the runs the profiled
build, the variant's tree and build were deleted and the cache kept only the datasets of seed
offset 0 (in `keep/build_phase4b`); the variant's patch and diff stay in `keep/variant_schedule/`.

**The new entry `multistate_workflows`** (`simple_commanders_test_highlevel_msworkflows.f90`,
registered with 8 threads and a 3 600 s timeout, reserving 11 cores; the process budget rises from
40 to 41 for it, as this phase of the plan decides). Its dataset is the state fixture's
compositional beta-galactosidase states: 3 000 particles of the tetramer, the tetramer less chain D
and the A-B dimer (populations 50, 30 and 20 %), SNR 1, with the contrast transfer function. It runs
three workflows, each in its own sub-suite:

- `cavgs`: `solve3D_cavgs nstates=3` on the class averages of an oracle 2D classification, then
  `solve3D nstates=3 cavg_ini_ext=yes` from the labels and poses it hands to the particles. The
  fixture's views mode learnt a `state_classes` option: 50 views of 60 particles, each 2D class one
  view and one true state (150 classes), so the classes are pure in state and the route is judged on
  its 3D state assignment alone. The views mode now also scores its 2D solution (a particle with no
  angles and no score counts as never searched, and `solve3D` refused `cavg_ini_ext` on such a
  project, "no 2D solution").
- `shc`: `solve3D nstates=3` from the particles (random poses, no states), whose stages 1-2 run
  `shc`.
- `states`: `refine3D_states nstates=3` from a single-state solution. Decision: the solution is the
  true poses with the consensus map reconstructed from every particle in one state
  (`reconstruct3D`, registered as the project's state 1 map), so the test isolates the state
  refinement from the errors and the cost of a single-state `solve3D`.
- **Metrics** are those of the multi-state assignment test, now with every estimated map docked in
  both hands onto every truth: the docked correlations give the map mapping, the docking onto its
  mapped truth the frame in which the fit, the discrimination, the poses (carried by the docking's
  rotation, shift and hand, up to each truth's declared point group) and the shifts are judged. The
  fraction of particles a workflow leaves inactive (state 0) is reported; they count as wrong in the
  accuracy, while the adjusted Rand index sees them as one more cluster. A truth without a map gets
  the lowest discrimination values (margin -2, difference correlation -1), so an ungated row cannot
  fail on a missing map.
- **Path assertions**, from the record the refine3D strategy writes for every iteration
  (`BENCH/REFINE3D_STRATEGY_BENCH.txt`: refine mode, directions, states): the blocks of
  iterations must contain the stage ladder in order, every iteration must run the requested states,
  and none may exceed the multi-state cap of 2 500 directions. `shc`: `shc` over 500 and 1 000
  directions, `prob` over 1 000, `prob_neigh` over 2 500. The class-average stage: `prob_neigh`
  over 500, `prob` over 500 and 1 000, `prob_neigh` over 2 500; its `solve3D`: `prob` over 1 000,
  `prob_neigh` over 2 500. `refine3D_states`: the neighbourhood stage, `prob_neigh` over 2 500,
  with the states `flex_pca` delivered (at least two).
- The helpers `finite_min`, `offdiag` and `maps_distinct` moved from the multi-state assignment
  test into `simple_test_truth_metrics.f90`, which both tests use (unit-tested in `truth metrics`,
  46 checks).
- Controls: shuffled labels and duplicated maps in every run; `suite=<workflow>_identical` runs a
  workflow on identical states, `suite=hard_<workflow>` on the conformational states.

**Findings** (for Phase 5; the tests record them, nothing was tuned around them):

- **The class-average route never separates the states.** `solve3D_cavgs` on 150 state-pure class
  averages ends in every run with the particles' labels at chance (adjusted Rand index median 0,
  at most 0.07; accuracy 0.50, the collapsed predictor) and at least one state empty, so
  `solve3D cavg_ini_ext=yes`, which refuses an empty state, never runs (40 of 40 runs, 1 and 8
  threads; the same with the earlier schedule).
- **`refine3D_states` drops the 20 % dimer state.** Its `flex_pca` initialization delivers two
  states and leaves exactly the 600 dimer particles inactive in every run (40 of 40), so the
  accuracy is 0.80 and the dimer has no map. The two states it keeps are right (labels exact, maps
  0.993 and 0.977, rotation error 1.44°).
- **`solve3D` from the particles merges two states in a third of the runs**: in 11 (8 threads) and
  12 (1 thread) of 20 runs it separates all three (adjusted Rand index 0.92-1.0), in the others it
  merges two (0.46-0.75) or empties one.

#### Card: `multistate_workflows`

- **Claim.** `solve3D nstates=3` from the particles finds three compositional states of
  beta-galactosidase, their maps and the particles' poses, far above chance; the class-average
  route and `refine3D_states` are profiled against the same truth.
- **Path assertions.** The stage ladder of each workflow (above), from the refine3D strategy's
  record; they held in every run of every profile on this build.
- **Profile of the gated route** (`keep/profiles/phase4/current_t8_shc` and `current_t1_shc`, 20
  seed offsets each, at 8 and 1 threads). Medians at 8 / 1 threads, robust spread s, worst run:
  adjusted Rand index 0.970 / 0.992 (s 0.042 / 0.012, worst 0.467 / 0.462), accuracy 0.987 / 0.998
  (worst 0.551 / 0.546), largest population error 0.011 / 0.001 (worst 0.243 / 0.249), rotation
  error 1.81° / 1.79° (worst 2.05° / 2.13°), shift error 0.070 / 0.077 pixel (worst 0.098 / 0.104),
  map correlations 0.95-0.98 with FSC 0.143 at 4.5-4.7 Å, discrimination margin 0.43 (worst -0.79,
  the merged runs), map mapping agreeing with the label mapping in 17 of 20 runs.
- **Floors** (written in the source with their origin): adjusted Rand index at least 0.30
  (admissible 0.13-0.425 at 8 threads, 0.04-0.450 at 1), largest population error at most 0.30
  (admissible 0.259-0.452 and 0.251-0.494), median rotation error at most 2.5° (admissible from
  2.12° and 2.20°), median shift error at most 0.15 pixel (admissible from 0.107 and 0.115). The
  label floor clears the controls but not the merging: a run that merges two states still passes it.
  **Recorded as a profile, with the reason:** the accuracy (worst 0.551 against the collapsed
  predictor's 0.50 plus 3 s, 0.556), the discrimination and the map fits (bimodal with the merging),
  and the whole of the `cavgs` and `states` routes (their failures are deterministic: no floor can
  clear the controls), which the entry therefore does not run; they run by suite.
- **Profiles of the other routes** (`current_t{1,8}_cavgs`, `current_t{1,8}_states`): as in the
  findings. Their runs at 8 threads take 408 s (`cavgs`) and 154 s (`states`) on average.
- **Hard profile** (`current_t8_hard_*`: the tetramer and the tetramer with chain D rotated by 10°,
  populations 80 and 20 %, SNR 0.5): `solve3D` stays at chance in every run (adjusted Rand index at
  most 0.006) while still finding the poses (1.74°); the class-average stage collapses in every
  run; `refine3D_states` leaves 10-26 % inactive and reaches an adjusted Rand index of 0.42 at the
  median (from -0.007 to 1.0), accuracy 0.79 against the collapsed predictor's 0.80.
- **Negative controls.** Shuffled labels in every run: adjusted Rand index -0.003 to 0.005, accuracy
  inside the null in every run, below the label floor. Duplicated maps in every run: never distinct.
  Identical states, once per workflow (`logs/p4_control_*.log`): adjusted Rand index -0.0004
  (`shc`), 0.003 (`cavgs`, whose stage again empties a state) and 0.003 (`states`, one state and
  14 % inactive), all below the floor; on the final build `shc_identical` passes its checks
  (`logs/p4_final_shc_identical.log`).
- **Determinism and cost.** Reproducible at one thread: two runs of `suite=states` at offset 0 give
  identical `metrics.tsv`, equal to the profile's (`logs/p4_det_states_a.log`, `_b.log`). The entry
  (the `shc` route) takes 1 650 s at 8 threads with two others beside it, 1 443 s through CTest on
  the final build (`logs/p4_ctest_multistate_workflows.log`), peak 10.8 cores, 2.6 GiB; at one
  thread 1.9-4.4 h.

**The schedule comparison.** The paired variant is a copy of the tree with only the schedule's
production consumers reversed (`keep/variant_schedule/`: `schedule_consumers_revert.patch`, its
result against the tree `variant_vs_tree.diff`, three files): in `simple_solve3D_controller.f90`
the switch to `prob_neigh` at the first stage of 2 500 directions, the 2 500-direction cap and the
stage-2 exit on the scan; in `simple_convergence.f90` the `scan_exit` rule; in
`simple_refine3D_states_controller.f90` `STAGE2_NSPACE` back to 5 000. Decision: the inert
`scan_exit` key stays declared in the parameters, so both builds compile byte-identical test code
(the convergence tester's case still sets it). Both builds ran the same test code, the same cached
datasets and the same seed offsets, at 8 threads (`keep/profiles/phase4/variant_t8_*`). Reported,
not gated (variant against the present schedule):

- `solve3D` from the particles: all three states separated (adjusted Rand index above 0.9) in 16 of
  20 runs against 11 of 20; by seed, 6 runs separate only with the earlier schedule and 1 only with
  the present one, 3 merge with both. Medians: adjusted Rand index 0.999 against 0.970, accuracy
  0.9997 against 0.987, map correlations 0.985 / 0.968 / 0.972 against 0.982 / 0.954 / 0.961,
  discrimination margin 0.446 against 0.432, rotation error 1.61° against 1.81°, shift error 0.073
  against 0.070 pixel. Its runs take 2 480 s against 1 658 s on average (the variant runs `prob` in
  stage 5 and has no stage-2 exit).
- The class-average route: the same collapse with both (every run leaves a state empty); the
  variant's runs take 665 s against 408 s.
- `refine3D_states`: the same labels (the dimer dropped by `flex_pca` with both), the same maps
  (0.993 / 0.977), rotation error 1.15° with 5 000 directions against 1.44° with 2 500; 167 s
  against 154 s.
- So on these known states the schedule trades about a third of the time of `solve3D` from the
  particles for fewer runs that separate all three states (16 to 11 of 20, a paired sign test of
  6 against 1, one-sided p = 0.06) and slightly coarser poses; it changes neither the
  class-average route's collapse nor `refine3D_states`' labels.
- **Ruling of 2026-10-17: the multi-state schedule stays.** The maintainer accepted the comparison
  as it is ("0.970 is good enough"): the present schedule is kept, neither reverted nor split.

**Review.** A separate reader reviewed the test code; applied: maps required only for populated
states (an emptied state keeps its number and loses its map) with the number of empty states
reported, the volume name from `refine3D_state_vol_fname`, comments on why the consensus map is
registered by the test (`reconstruct3D` with `mkdir=no` registers nothing) and on the class-average
maps the next step does not read, a released string, a redundant absolute path removed.

**Exit criteria.**

- The gates pass: the fast gate (Debug, 15 suites, 9.1 s, `logs/compile_debug_full_p4.log`),
  `unit_numerics` (1 079 checks, `runs/p4_unit_numerics_final/run.log`), the new entry through CTest
  on the final build (1 443 s), the registry check and `check_descr.py`; no mode changes.
- Every test's floors admissible or the case recorded as a profile: the `shc` route's labels,
  populations, poses and shifts have admissible floors; its accuracy, discrimination and fits, and
  the `cavgs`, `states` and hard cases are recorded as profiles, with the reasons above.
- The comparison reported (above).

Files changed in this phase: in `src/main/commanders/test/`: new
`simple_commanders_test_highlevel_msworkflows.f90`, `simple_commanders_test_highlevel.f90` (the
commander type), `simple_commanders_test_highlevel_multistate.f90` (the helpers moved out),
`simple_test_state_fixture.f90` (state-pure classes, the scored 2D solution),
`simple_test_truth_metrics.f90` and `simple_test_truth_metrics_tester.f90` (the helpers),
`simple_test_registry.f90` (the entry); `src/main/exec/simple_test_exec_highlevel.f90`,
`src/main/ui/simple_test/simple_test_ui_highlevel.f90`, `production/CMakeLists.txt` (the entry, the
budget), `doc/policies/test_environment_policy.md` (the entry, the counts; file-table row added)
and this note.

### Phase 5, 2026-10-18: repairs and the audit's actions

All paths below are under `/home/elmlundho/agent_runs/test_truth/scratch` unless they start with
`src/`, `doc/`, `production/` or `scripts/`. Release builds with tests: `keep/build_phase5b` (the
first after-profiles), `keep/build_phase5e` (after the class-average rulings: most after-profiles)
and `keep/build_phase5f` (the final code: the flex profile, the last verifications and the
after-timing); the intermediate builds were deleted for disk, their evidence is in the logs. The
experiments ran in `keep/exp_tree`, a copy of the tree with switches read from the environment, and
its build `keep/build_exp`; nothing of it is in the repository. The journal of the phase is
`state/journal_phase_5.md` beside this run's scratch directory. The maintainer stopped the phase once
to deliver rulings (2026-10-18); the profiles then running were restarted as the ruling asked.

**Rulings applied in this phase** (`spec/rulings.md`; recorded where they belong: the common rules,
the test table, the Phase 4 entry and here):

- 2026-10-17, fix as you go: defects are fixed at once and listed in the file table; a design
  question gets a diagnosis and a proposed fix under "For the maintainer" below; floors move only by
  the floor rule.
- 2026-10-17, `refine3D_states`: state 0 means exclusion; the dimer at state 0 was a defect, and the
  gate requires every true state resolved.
- 2026-10-17, `solve3D_cavgs`: tested on simulated images imported with `import_cavgs`, scored per
  image (single mapping, adjusted Rand index, per-state recall) and per state map.
- 2026-10-17, the multi-state schedule stays (Phase 4 entry).
- 2026-10-17, the review list of Phases 0 to 3 and the production defects (both handled below).
- 2026-10-18, the class-average collapse: find whether the test or the code is wrong before calling
  it a design question (inputs checked, then controlled experiments); `flex_pca_blobs` part 3 on
  states `flex_pca` can resolve, its floors from a profile by the floor rule.
- 2026-10-18 (later): multi-state `solve3D_cavgs` runs direct `shc` in stages 1 and 2, as particle
  `solve3D`; `prob_athres` defaults to 10 in every `solve3D_cavgs` run and the taper is deleted; the
  test's class averages are about 600 even/odd pairs of independent noise, without CTF, with small
  shifts; the profiles the rulings touch are redone, finished profiles of untouched entries kept.

**Production defects fixed** (each with its file-table row):

- `single_workflow_fcc` stopped on an empty Fourier range in `autorefine3D_nano`: the FSC=0.143
  clamp of the matching low-pass (commit 65bee3073) never compared with the high-pass limit, so a
  3 Å high-pass and a 12 Å clamp left no band. The clamp now applies only when the band stays
  valid (`simple_matcher_smpl_and_lplims.f90`); the entry passes (583 s, all 20 runs of its profile).
- Shared-memory `refine3D` did not draw random states for a multi-state fresh start, as the
  distributed path does: it now seeds them (`simple_refine3D_strategy.f90`, fresh start with one
  populated state only).
- `flex_pca` in `refine3D_states` set the dimer (20 %) to state 0: its population-floor path then
  pruned every state below `min_neff`=2000, which no state of 3 000 particles reaches. It prunes
  below the smaller of the two floors now (`simple_flex_pca_application.f90`). The prune fallback
  also kept a state without mass, which stopped the distributed PCG (`simple_flex_pca_state_service.f90`).
- The distributed PCG stopped on a state whose members were all in one halfset (class-average
  halves are assigned independently); it carries the state's previous maps forward, as for an empty
  state (`simple_rec3D_pcg_strategy.f90`).
- `autoscale` planned a crop above the input box (88 for 64) and `solve3D_pyramid_blobs` stopped on
  it; the crop is clamped to the box (`simple_magic_boxes.f90`).
- The pose-initialization pass of `refine3D_states_cavgs` searched shifts first on rows without a
  pose (the shift fitted against an arbitrary projection then biased every orientation: on
  noise-free reprojections of the exact reference 45 % of the poses within 5° against 97 % without
  the shift search). Pose-less rows now skip the shift search in that pass
  (`simple_external_reference_pose_initialization.f90`, `simple_commanders_refine3D.f90`); the
  sub-averages' pose error in `flex_pca_blobs` part 3 fell from 24.5° to 4.7°.
- `simulate_particles` ignored `state_pops` beside an orientation table holding states, accepted an
  empty state and accepted CTF modes it cannot apply; it refuses all three
  (`simple_commanders_sim.f90`, its interface and the user guide).
- By the rulings of 2026-10-18: multi-state `solve3D_cavgs` runs direct `shc` in stages 1 and 2,
  `prob_athres` defaults to 10 in `solve3D_cavgs` and the taper (`prob_athres_taper`) is deleted
  with its key, interface entry, function, test, caller and documentation.
- Identity: the one-thread `single_workflow_wurtzite` gives the same 140 production files (maps,
  projects over every segment, text) on the Phase 4 build and the Phase 5 build (`logs/
  p5_identity_report_wz.txt`); only the test's new pose-docking maps differ. The fixes above leave
  a run they do not concern byte-identical.

**The multi-state routes after the fixes** (20 seed offsets each at 8 threads):

- `multistate_assignment`: the merging of two blob states (Phase 3: 16 to 17 of 40 runs) came from
  the test, which ran `refine3D` under `euclid` from starting maps that carry no noise model;
  production's state assignment scores by correlation. Experiments (offsets 0-2): adjusted Rand
  index under `euclid` 0.58-0.90 with merges, under `cc` 0.88-0.97 with every state found;
  normalised starts helped only partly. The test runs `objfun=cc`; its profile resolves every state
  in 20 of 20 runs (smallest recall 0.953 for `prob`, 0.893 for `prob_neigh`), and its labels are
  gated again (ruling 5): adjusted Rand index 0.90 / 0.80, accuracy 0.95 / 0.92, smallest recall
  0.93 / 0.85, population error 0.03 / 0.04, no estimated state empty, map mapping agreeing, a
  discrimination margin of 0.30 (`keep/profiles/phase5/after`).
- `multistate_workflows`, `states` route: every state resolved in 20 of 20 runs (adjusted Rand index
  at least 0.9989, nothing in state 0); labels, populations, discrimination and poses gated by the
  floor rule (`keep/profiles/phase5/ms_states`).
- `shc` route (particle `solve3D`, its schedule kept by ruling): unchanged code; the Phase 4 profile
  stands (a state lost in 9 of 20 runs at 8 threads, 7 of 20 at one); its checks that every state is
  resolved stay gates and fail in those runs, its poses meet their ceilings.
- Class-average route on the ruled inputs (600 images, independent halves at SNR 2.5, no CTF,
  shifts within a pixel): every state resolved in 9 of 20 runs (Phase 4: 0 of 20), adjusted Rand
  index median 0.76 (0.40-1.00); poses gated by the floor rule (4.5°, 0.40 px)
  (`keep/profiles/phase5/ms_cavgs`). The localisation the ruling asked for (offsets 0-2, `runs/
  p5_loc`): particle `solve3D` on the same 1 200 half-images resolved all states in 0 of 3 runs
  (the tetramer merged with the tetramer less a chain, as in the `shc` route), the class-average
  route in 1 of 3, with `objfun=cc` in 2 of 3 (poses 7.5-9.6° against 1.6-3.1°), with
  `noise_norm=yes` in 1 of 3. The remaining failures are therefore not in the class-average path but
  in the multi-state assignment both share; the forced `euclid` trades labels for poses (below).
- The single-state class-average start (`suite=cavgs_single`, 20 offsets each,
  `keep/profiles/phase5/cavgs_single`): with `prob_athres` 90 (before) rotation error median 2.88°
  (2.66-3.12), map correlation median 0.979 (0.976-0.981), FSC 5.56 Å; with 10 (after) rotation
  error median 2.60° (1.23-3.35), map correlation median 0.910 (0.875-0.987), FSC 5.28 Å, shifts
  up to 0.34 px in one run. Recorded as the ruling asks; the lower map correlation beside a better
  FSC is listed for the maintainer.
- `flex_pca_blobs` part 3 on the compositional beta-galactosidase states (ruling; blob states with
  a large distinguishing blob were tried first and give pose-less sub-averages no orientation at the
  12 Å fit band, 43° pair error; an 8 Å band gave noisier poses): with the seeding fix the poses
  hold (median pair error 11.7°), but the states stay at the chance level in its profile (sub-average
  accuracy median 0.42 against a shuffled-label null of 0.37, 99.9th percentile 0.43; three of 20
  runs deliver fewer than two states), so no floor is admissible: the accuracy, the pose error and
  the state count are reported, the route's contract stays gated. Parts 1 and 2 pass 20 of 20
  (label accuracy median 0.989) (`keep/profiles/phase5/flex_f`).
- `solve3D_pyramid_blobs` (on the blob states differing in a lobe of 6 % of the mass): 20 of 20
  runs fail, 16 on production's own stops (the reference check, or a state with fewer than 20
  sub-class averages after the state step collapsed), 4 on the assignment accuracy
  (`keep/profiles/phase5/after_e`). Listed for the maintainer.

**The review list of Phases 0 to 3** (ruling 5), all applied: the label rows of
`multistate_assignment` (above); the profile's failure rate combines each execution's outcome with
each metric (an execution that failed after passing rows fails), with a tester case; a stale
`metrics.tsv` is removed when an execution starts; an inherited SIG_IGN is kept; the child's SIGTERM
handler `_exit`s and the gate flushes each row; the child's process group is killed before it is
reaped (`waitid` with `WNOWAIT`); the marker search reads `/proc/<pid>/environ` to the end; the
`/proc` guarantees are stated Linux-only and their probes skipped elsewhere; an interrupted execution
has its own outcome, outside the failure rate (seen in this phase: the stopped profile recorded its
running execution as interrupted); packing starts the largest reservation first within cores,
memory and disk, with memory and disk in the registry and scaled reservations from master and worker
threads; `test=identity` compares every project segment and the text files without the run-directory
prefix and fails on a pattern matching nothing (which found a stale `single_workflow` pattern, now
corrected); run-directory tables of up to 1 024 characters; `nseeds` defaults to 20, `seed_offset`
survives `script_exec`, `suite=` needs one entry; seed offsets times the odd stride 2654435761 in
64-bit integers; one overflow-free FNV-1a hash with published test vectors; the fixture key holds the
build's commit; the point-group check is relative (the declared group against the next larger); no
negative coordinate before `pdb2mrc`; the fixture checks are their own entry `state_fixture`, nested
under the record counts; checks fail on a non-finite value; report-only gates became metrics
(`simulated_workflow` extracted particles and classes, `mini_stream` classes, `pcg_recon`'s stage
decisions); the stray shell text left Phase 0 item 4.

**Every other entry, by its card** (before: the Phase 1 profile `keep/profiles/phase1/<entry>`;
after: `keep/profiles/phase5/after_e`, `after_singlewf`, `after`; 20 seed offsets each):

- `mini_stream_6vxx`, `_1jxy`: the 6vxx picking (19 of 48) had not regressed: the 0.40 floor came
  from one unseeded run. Floors by the floor rule from the after-profile: pick ratio 0.05 (6vxx:
  median 0.33, worst 0.19) and 0.25 (1jxy: median 0.49, worst 0.40); defocus ceiling 0.13 µm (worst
  0.125); classes a metric; path assertions on the `solve2D` and `refine2D` records; consistency
  gates as rows; the negative control is the fit against the farthest other micrograph's truth (the
  neighbour's, 0.25 µm away nominally, came within the ceiling in some runs; verified 85 of 85).
- `simulated_workflow_6vxx`, `_1jxy`: 200 particles in 20 movies; a wrong-molecule control; the
  distributed step asserts its second part's record. `solve3D` still misses the map in 7 of 20
  (6vxx) and 1 of 20 (1jxy) runs, at the control's correlation (median 0.53), so no correlation
  floor is admissible: the correlations and the 6vxx FSC are reported, the 1jxy FSC ceiling (40 Å)
  stays. Listed for the maintainer.
- `single_workflow_fcc`, `_wurtzite`: the pose metric was not symmetry-aware (the pair error read
  correct poses of these near-symmetric clusters as random, 42.7-45.6° in every run). It is now the
  median rotation error in the truth's frame after one docking of the `autorefine3D_nano` map (to
  2 Å, where the lattice fringes fix the frame), up to the lattice's group (O, D6) and the in-plane
  rotation by 180° (both lattices are nearly centrosymmetric): fcc median 4.9°, wurtzite 1.5°. The
  relative ceilings (a stage loses at most 0.5°) are admissible (worst excess +0.19°).
- `pcg_recon`: relabelled library, its stage decisions as rows; deterministic (two runs, identical
  values), so one run is its profile.
- `simulate_particles`: images checked against the reprojections at their recorded poses (50
  images, every one closer to its own pose than to another's); the interval rows of the CTF and
  shift draws are checks of the simulator's contract (the largest of many uniform draws sits just
  inside its bound, so no floor is admissible), their maxima reported.
- `solve3D_addon`: the pose gate relative to the base (the cohort's excess over the frozen
  particles) with a shuffled-cohort control; the after-profile is wider than Phase 1's (median
  3.38°, s 2.14, worst 5.89), so the ceiling moves from 6.5 to 10.0° by the floor rule (the control
  sits at 23.9°, s 3.2).
- `cont_refine3D_1jxy`: the distributed stage asserts two parts; only check rows were added, its
  Phase 1 profile stands.
- `stream_preproc`: the defocus ceiling 0.20 µm (Phase 1 profile) holds (worst 0.134); worker jobs
  asserted (`JOB_FINISHED_*`, every recorded exit status 0); the farthest-truth control.
- `species_discovery`: the radial-shell check skipped where the two-class fit is not admissible;
  aluminium recall 0.88 (Phase 1); the second element's recall floor moves from 0.90 to 0.89 by the
  floor rule (Ni: median 0.977, worst 0.908).
- `single_atoms_stats` and its three CdSe variants: a shifted-atoms control; floors admissible.
- `flex_pca_blobs`: parts 1-2 on the fixture's known-pose mode, part 3 above; the distributed run
  asserts its worker jobs. `solve3D_pyramid_blobs` on the views mode with one-to-one mapping,
  matched and docked maps, discrimination and a shuffled control.
- `state_fixture` (new): passes 20 of 20 (11 s from the cache).

**Programs that need inputs.** `pcg_frac_update` builds its own project from the fixture (600 blob
particles at their true poses) when run without `projfile=` and is the library entry of that name;
it had printed its usage and exited 0 (its interface still required `pgrp` and `mskdiam`), and a
one-state fixture crashed the truth writer: both fixed, it passes in 3 s. `rec3D_backends` and the
input generators stay tools.

**Packing.** `test=nightly` packs by cores, memory and disk; every high-level CTest entry carries
`PROCESSORS` equal to its reservation instead of `RUN_SERIAL`, which the registry lint checks. Wall
time on an otherwise idle machine, before (Phase 4 code, `RUN_SERIAL`; `logs/p5_timing_before.log`)
and after (final code; `logs/p5_timing_after.log`):

| run | before: Phase 4 code, `RUN_SERIAL` | after: final code, `PROCESSORS` |
|---|---|---|
| `test=nightly ncunits=24` | 1 927 s (3 entries stopped early) | 3 724 s (2 failed) |
| `ctest -L highlevel -j 24` | 4 122 s (4 failed) | 2 721 s (2 failed) |

CTest now packs and is a third faster. The nightly is slower for two reasons: entries that stopped
early before now run through (`single_workflow_fcc`, `flex_pca_blobs`, the pyramid's later steps)
and the workflows simulate 200 particles; and largest-first packing, as ruled, started the longest
entry, `simulated_workflow_1jxy` (3 cores, 2 239 s), only after about 1 450 s, behind the 9- to
11-core entries, so it alone set the wall time (item 12 below).

A profile with outputs deleted as executions end cannot use the dispatcher's disk budget, which
counts every started execution's disk (outputs stay in a nightly); the after-profiles therefore ran
partly in bounded bash loops (`p5_profile_rest.sh`, `p5_cavgs_single_par.sh`).

**For the maintainer** (diagnosis and a proposed fix each):

1. Multi-state assignment from random labels loses or merges states: particle `solve3D` (the
   `shc` route, 9 of 20 runs; the same 1 200 class-average halves imported as particles, 3 of 3)
   and class-average `solve3D_cavgs` (11 of 20) alike; the pyramid's state step collapses the same
   way. Proposal: reseed a state that falls below a population floor during stages 1-2 (not only an
   empty one, and in particle `solve3D` too), or a population-balanced or annealed state assignment
   in the prob stages.
2. `solve3D_cavgs` forces `objfun=euclid`: under `cc` the class-average route resolved more runs
   (2 of 3 against 1 of 3) but its poses were four times worse; production's state assignment
   (`simple_state_assignment`) already uses `cc`. Proposal: `cc` in the multi-state state stages, or
   a `cc` stage before `euclid`.
3. `refine3D` under `euclid` from given starting maps without a noise model merges states in its
   first iterations (the multistate test now runs `cc`). Proposal: a `cc` first iteration, or a
   sigma2 estimate, when vol1..volN are given.
4. `refine3D_states_cavgs` aligns pose-less sub-averages by a global c1 search against the
   consensus; on beta-galactosidase (a consensus nearly two-fold symmetric) and at the sub-averages'
   12 Å band the states stay at chance. Proposal: carry the members' poses into `cls_expansion`'s
   class records, copy them to both halves in `cavgs_work_project`, refine locally
   (`pose_policy=local`), and cap the route's `lpstop` at the sub-averages' band.
5. The single-state class-average start at `prob_athres` 10 gives maps that correlate less with the
   truth (median 0.910 against 0.979) at a slightly better FSC; worth a look at the map's scale or
   mask before the default settles.
6. `simulated_workflow`: `solve3D` on 200 picked particles still misses the map in 7 of 20 runs
   (6vxx). Proposal: more particles (or movies) for 6vxx, or a cleaner initial model step; then
   restore the correlation floor from a profile.
7. `solve3D_pyramid_blobs` runs on the blob states that differ in a lobe of 6 % of the mass, the case
   the ruling called unfair for `flex_pca_blobs`; it fails 20 of 20. Proposal: the compositional
   beta-galactosidase states, after item 4.
8. `flex_pca`'s population-floor path gives particles it does not support random labels, not state
   0 (`simple_flex_pca_state_service.f90`, around line 633). Proposal: state 0, or the nearest state.
9. The external-reference pose pass of `solve3D` copies 2D shifts onto random orientations and then
   searches shifts first, as `refine3D_states_cavgs` did; no test covers it. Proposal: the same
   pose-less rule there.
10. `simple_cmdline` stops with exit status 0 when required keys are missing (it prints the usage),
    so a CTest entry missing a key passes; `pcg_frac_update` did. Proposal: a non-zero status for
    `simple_test_exec` (or all executables) on missing keys.
11. The dispatcher's disk budget cannot know about outputs deleted during a profile. Proposal: an
    option that counts only running executions' disk when outputs are removed as they end.
12. Largest-first packing (as ruled) starts the long one-thread workflows last, which doubles the
    nightly's wall time. Proposal: order by expected core-seconds (the registry's measured wall time
    times the reservation), largest first, which keeps large reservations from starving and starts
    the longest entries early.

**Exit criteria.**

- The gates pass: the Debug fast gate (15 suites, `logs/compile_debug_full_p5_final.log`),
  `forked_process` (70 checks), the library entries `pcg_recon` and `pcg_frac_update` through CTest,
  the registry lint, `check_descr.py`; no mode changes. Of the high-level entries, those failing
  in the final nightly fail on the behaviours listed for the maintainer (items 1 and 7), as the
  rulings ask (a lost state is a defect to keep visible): `multistate_workflows` and `solve3D_pyramid_blobs`, in the nightly and in CTest alike (`timing/after/`).
- Every change carries its profile: before from Phases 1, 3 and 4, after from this phase (above,
  20 seed offsets each), the identity check where the code was meant to be unchanged, and one run
  each for the deterministic library entries.
- The packed suite's wall time is recorded (above).

Files changed in this phase: production code `src/main/flex/run/simple_flex_pca_application.f90`,
`src/main/flex/run/simple_flex_pca_state_service.f90`,
`src/main/strategies/search/simple_matcher_smpl_and_lplims.f90`,
`src/main/strategies/parallelization/simple_refine3D_strategy.f90`,
`src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90`, `src/utils/simple_magic_boxes.f90`,
`src/main/simple_external_reference_pose_initialization.f90`,
`src/main/commanders/simple/simple_commanders_refine3D.f90`,
`src/main/solve/simple_solve3D_controller.f90`, `src/main/commanders/simple/simple_commanders_solve3D.f90`,
`src/main/params/simple_parameters.f90`, `src/main/params/simple_parameters_parse.f90`,
`src/main/ui/simple/simple_ui_solve3D.f90`, `src/main/solve/simple_solve3D_utils.f90`,
`src/main/solve/simple_solve3D_reseed_tester.f90`, `src/main/solve/simple_state_assignment.f90`,
`src/main/solve/simple_state_assignment_tester.f90`, `src/main/commanders/simple/simple_commanders_sim.f90`,
`src/main/ui/simple/simple_ui_sim.f90`; the runner `src/utils/simple_forked_process.f90` and its
tester, `src/utils/simple_test_utils.f90`, `src/utils/math/simple_rnd.f90`,
`src/utils/text/simple_string.f90` and its tester, `production/simple_test_exec.f90`,
`production/CMakeLists.txt`, `scripts/check_test_registry.py`; the tests in
`src/main/commanders/test/` (the dispatcher, suite runner, gate, identity, profile statistics,
registry and their testers; the state fixture; the truth metrics and their tester; the new entry
`simple_commanders_test_highlevel_fixture.f90` and the commander type in
`simple_commanders_test_highlevel.f90`; the multistate, msworkflows, flex, pyramid, stream,
test_stream, workflow, addon, single, cont, frac, pcg and backends test modules),
`src/main/flex/run/simple_flex_pca_application_tester.f90`,
`src/main/exec/simple_test_exec_highlevel.f90`, `src/main/ui/simple_test/simple_test_ui_highlevel.f90`;
documents `doc/policies/test_environment_policy.md`, `doc/policies/3D/solve3D_cavgs_policy.md`,
`doc/policies/3D/solve3D_pyramid_policy.md`, `doc/policies/stream/stream_3D_ingestion_policy.md`,
`doc/user_guide/simulate_particles.md` (file-table rows added for each) and this note.

### Phase 6, 2026-10-18: interrupted for the synchronisation with master

Phase 6 started at 19:34 and the maintainer stopped it at about 19:43. Master had moved nine commits
past this run's base (25d256886 to 43de9b340), several of them in files this run changed, so the
run's work was merged onto master before going further. Before the stop the phase had edited the test
policy (four kinds of test, the floor rule, packing, `test=nightly`, the cards, the dataset layer)
and the `simulate_particles` user guide; those edits are in the merged tree. The final build, the
last nightly and the report were not finished.

**The merge** (by hand, 2026-10-18): the run's complete working-tree diff, applied three-way onto
43de9b340. Eleven files had changed on both sides; the decisions:

- **The FSC=0.143 matching clamp.** Master reverted the clamp (fb4e653ac: it clamped the `solve3D`
  stage ladder), so this phase's repair of it goes with it. The empty Fourier band that stopped
  `single_workflow_fcc` came from the clamp itself.
- **The `prob_athres` taper.** Master had moved the symmetry axis search to stage 4 (13795c1d9) and
  adapted the taper to it. The taper is deleted as ruled; master's stage constants
  (`SYMSRCH_STAGE` 4, `NOSYMSRCH_STAGE` 3) and its FSC=0.5 stage rule (cd13708eb) stay. The taper's
  unit test gives way to a test of `solve3D_symsrch_stage` (master's two assertions), and
  `solve3D_policy.md` no longer mentions the taper.
- **`simulated_workflow`.** Master's sizes (61695ce7a: 6vxx 160 particles from 100 reprojections in 16
  movies, seed 20260923; 1jxy 100 particles in 10 movies) replace this run's 200 in 20. This run's
  gate rows, wrong-truth control and distributed-parts check stay. Its correlations stay reported,
  not gated, as the profile of the 200-particle version required; master's sizes are not yet
  profiled.

**What no longer holds.** The profiles and floors of Phases 1 to 5 were measured on 25d256886 with
this run's changes. Master changed results in `solve3D` and `refine3D` since: the clamp revert, the
FSC=0.5 stage omission, the nonuniform filter from stage 5 with its mid-stage exit (8c8d8b427), the axis
search at stage 4, and the speed week (04e520e78), which is meant to be identical apart from keeping
back-filled projection labels in volume assembly. So every entry that runs either needs its profile
again. An entry whose outputs stay byte-identical at one thread between the old and the merged
code keeps its profile.

**The continuation** is a new run on the merged code: the x86 identity check of the speed week
(25d256886 against 04e520e78), the identity check that sorts the entries, the profiles and floors
again where needed, the rulings not yet delivered to this run, the suite's timing, the documents
and the close.
