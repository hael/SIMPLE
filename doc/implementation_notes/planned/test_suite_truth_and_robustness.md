# A test suite that knows the truth: correctness, robustness and multi-state

Implementation note, 2026-10-16. **Proposal; nothing implemented.** Written
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
| multi-state `solve3D` | `solve3D_cavgs nstates=3` then `solve3D nstates=3 cavg_ini_ext=yes`, and `solve3D nstates=3` from the particles (the `shc` route), on the compositional states | missing chains | correctness with calibrated floors; robustness profile |
| state refinement from a consensus | `refine3D_states` on a single-state solution of the compositional states, through its neighbourhood stage | missing chains | correctness, with a path assertion that the neighbourhood stage ran |
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
- **`simple_test_exec test=identity dir1=<run> dir2=<run>`.**
  - It byte-compares two run directories over the entry's declared outputs.
  - It ignores the declared run-specific fields: the run identifier in
    `solve3D_manifest.txt`, the layout digest of `sigma2_state.bin` and the
    manifest's `stack_digest`.
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
  - the manifest's `stack_digest`.

  Multi-threaded runs never decide identity.
- **Findings, not stops.** The following are findings, recorded with their
  evidence, after which the run goes on:
  - a test that fails;
  - a distribution that overlaps its floor;
  - a multi-state workflow that cannot separate its states;
  - a production defect found by a new test.

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

### File table

| file | why | phases |
| --- | --- | --- |
| `production/simple_test_exec.f90` | the seed offset | 0 |
| `src/main/exec/simple_test_exec_highlevel.f90`, `src/main/apis/simple_test_exec_api.f90` | dispatch of the new test programs | 0 3 4 5 |
| `scripts/check_test_registry.py` | the registry check reads the SIMPLE registry | 0 5 |
| `src/utils/simple_test_utils.f90` | the seed offset | 0 |
| `src/utils/simple_forked_process.f90` | the dispatcher's children: own session, group signals, deadlines | 0 |
| `src/main/ui/simple_test/` | the test front door's key; new and corrected test descriptions | 0 2 3 4 5 |
| `src/main/commanders/test/` | the tests: reporting, fixture, metrics, new tests, repairs | 0 2 3 4 5 |
| `src/main/commanders/simple/simple_commanders_sim.f90` | `simulate_particles` with several states | 2 |
| `src/main/ui/simple/simple_ui_sim.f90` | its interface | 2 |
| `src/main/params/simple_parameters.f90`, `src/main/params/simple_parameters_parse.f90`, `src/main/params/simple_parameters_phases.f90` | the population key, if one is needed | 2 |
| `doc/user_guide/simulate_particles.md` | the mixture | 2 6 |
| `src/main/flex/run/simple_flex_pca_application_tester.f90` | the flex tests moved onto the new fixture | 5 |
| `production/CMakeLists.txt` | registrations, labels, processor counts, the budget | 0 3 4 5 |
| `doc/policies/test_environment_policy.md` | the tools, the cards, the rules | 0 5 6 |

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

Nothing yet.

