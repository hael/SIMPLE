# Test Environment Policy

## Purpose and scope

This policy says what a test is in SIMPLE, where each test goes, how to write
one, and what runs when. It applies to every Fortran test in the repository
and to everyone who adds, changes or removes one.

The refactoring itself is documented in
`doc/refactoring_notes/completed/uniform_test_environment_refactoring.md` (the plan and
its batch records, section 9.7). Every test that existed before it is listed,
with its review verdict, in `doc/refactoring_notes/completed/test_review_record.md`; the
full inventory, `doc/code_overview/test_inventory.md`, is generated from
the sources and that record on every build and is not committed. This
policy is the short version for day-to-day work; when the two disagree, fix
the one that is wrong.

## 1. The environment

SIMPLE has one test executable, `simple_test_exec`, built by every
`compile_*.sh` unless it is given `--exclude-tests` (CMake option
`BUILD_TESTS`, default ON). Every test is
a procedure inside the library, not a program of its own, and every test is
reached through `simple_test_exec test=<name>`. `simple_test_exec test=list`
lists what exists.

Tests are grouped into suites, and a suite is what CTest registers: one CTest
entry per suite, never one per test. The entries fall into two tiers and four
labels.

| label | entries | when | what |
|---|---|---|---|
| `fast` | 15 area suites `unit_<area>` | every `compile_*.sh` build (unless `--exclude-tests`), before installation | unit tests of the library: hermetic, in-process, one thread, seconds |
| `library` | 5 library suites `lib_<area>` and the PCG entries `pcg_recon` and `pcg_frac_update` | nightly | longer numerical tests on generated data: realistic sizes, minutes |
| `highlevel` | 19 high-level gates | explicit CTest command only | long simulated pipelines and commander integrations, including independent molecular and nanoparticle suites |
| `platform` | `forked_process` (with the test dispatcher's own tests), plus `coarrays`, `openmp_offload` when CMake finds the capability | by hand, nightly where the machine has the capability, and `coarrays` during `compile_coarrays.sh` | tests that need child processes, a launcher or a device; `coarrays` is the capability-gated two-image synchronization suite |

The fast tier is the build-time gate. Library and supported platform tests may
run overnight; high-level tests run only when explicitly selected with CTest.

**The fast gate is part of the build.** Every `compile_*.sh` build
runs `scripts/run_fast_gate.sh` between `make` and `make install`. It first
runs `scripts/check_test_registry.py`, `scripts/check_descr.py` and
`scripts/check_flex_dag.py` (section 4.5), then
`ctest -L fast` with half the cores rounded up, then
`scripts/ctest_budget.py`, which
fails the build when an entry fails and only warns when the gate takes more
than 60 s of real time (`scripts/run_fast_gate.sh`, `scripts/ctest_budget.py`).
A budget warning allows installation to continue. A failed gate installs nothing. The per-entry timings are kept in
`build/test_runs/ctest_fast.log.timing.txt`. Tests are on by default in every compile script;
`--exclude-tests` (`BUILD_TESTS=OFF`) builds the library and executables only
and skips the gate, for when only the executables are needed.

Each fast suite has a 600 s CTest deadline to allow slower or shared hosts
to finish. This is separate from the warning-only 60 s total budget. CTest's
explicit per-test `TIMEOUT` takes precedence over `ctest --timeout`, so the
registration and gate script must agree. A suite that reaches its deadline
still fails the gate.

Fast suites use `OMP_WAIT_POLICY=PASSIVE` so idle OpenMP workers sleep instead
of consuming a shared host's CPU quota. Thread-consistency checks can request
larger teams than the default `OMP_NUM_THREADS=1`; active waiting in those
teams can starve the other suites on CPU-limited login nodes.

The fast-suite environment also sets `OMP_THREAD_LIMIT=4` and
`OPENBLAS_NUM_THREADS=1`. The OpenMP limit bounds actual
teams even when a tester calls `omp_set_num_threads(omp_get_num_procs())`;
the BLAS limit prevents independent numerical-library worker pools from
multiplying across concurrent suites. Thread-consistency checks still run
with multiple workers. Full-machine thread stress can be run directly through
`simple_test_exec` without these CTest limits.

`compile_coarrays.sh` runs the same fast gate as the other build scripts and
then runs the capability-gated two-image synchronization suite before
installation. A failure in either gate prevents installation. The coarray
suite is registered only when `USE_COARRAYS=ON` and an OpenCoarrays launcher is
found. It can also be rerun explicitly after a coarray build with
`cd build && ctest -R '^coarrays$' --no-tests=error --output-on-failure`.

**The process budget.** The number of CTest entries is fixed in
`SIMPLE_CTEST_BUDGET` (`production/CMakeLists.txt`, currently 42: 15 fast,
7 library, 19 highlevel and 1 platform) and configuration fails when it does not
match. A CTest entry is an isolation unit, not a place for one more check:
checks are added inside existing suites. A new entry needs a stated reason and
the owner's agreement, and is recorded in the plan.

**Every entry** runs in its own directory `build/test_runs/<entry>`, with its
own `TIMEOUT` and an explicit `OMP_NUM_THREADS` (1 for the fast tier), and
`simple_test_exec` fixes the seed of its process (below). Each suite writes a report,
`SIMPLE_TEST_<suite>_<date>/simple_test_<suite>_report.txt`, with every check,
every failure by name and the time of each sub-suite.

### 1.1 The fast area suites

Each area suite runs its sub-suites in one process. One sub-suite is one
tester module (section 4.1).

| suite | sub-suites |
|---|---|
| `unit_core` | ANSI formatting, string, syslib, fileio, stack I/O, class sample I/O, character hash, hash, value-reference hash, linked list, record list, command line |
| `unit_ori` | orientation, orientation collection, symmetry, Euler shift |
| `unit_image` | image, mrc2jpeg, mrc validate, image header, Fourier iterator, B-spline smoother, masks, nano mask, volume shape, binary image, segmentation, trailing-reconstruction blend, CTF, image serialisation |
| `unit_numerics` | online variance, random draws, affinity propagation, hierarchical clustering, k-means, k-center, Gaussian mixture, extreme deconvolution, statistics, linear algebra, Kaiser-Bessel kernel, search/sort/locate, decay schedules, PCA, cavg quality relations, cavg quality selection, diffusion-map graphs, optimisers, low-pass stages, shift search, truth metrics (state matching, adjusted Rand index, the null accuracy, map matching, discrimination, the frame, hand and point-group conventions of poses) |
| `unit_project` | STAR file, STAR project, STAR stream export, binoris, project records, project merge, class compatibility, particle sieve, micrograph import, micrograph selection, micrograph subset, particle subset, optics groups, optics maps, solve3D manifest, project superset, solve3D addon report, solve3D chunks plan, state weight set, cavgs work project, eo partition |
| `unit_motion` | motion gain, motion model, particle extractor |
| `unit_ui` | UI JSON, GUI metadata, GUI assembler, UI hash, UI visibility |
| `unit_ipc` | IPC TCP socket, HTTP POST, persistent worker message, persistent worker server, stream pipe (localhost only) |
| `unit_reconstruction` | rec3D backend, observation noise, class-average carry-over, trailing chain identity, frozen accumulator, weighted insertion, volume pair metrics, nonuniform filtering, solve3D_cavgs stage rules, state assignment rules, atomic prior |
| `unit_pftc_align2D3D` | polar correlation, continuous in-plane, refine3D in-plane state, 2D probability table I/O, 3D probability table, sigma2 state, sigma2 bootstrap, cavg registration, cavgs state match |
| `unit_cart_align3D` | Cartesian Fourier, cart calculator, pose optimizer, pose strategy, pose statistics |
| `unit_heterogeneity` | flex PCA, flex PCG operator, flex cls_expansion, flex cls_expansion kernel |
| `unit_parallel` | qsys control, qsys environment, test tools (the profile statistics, packing, identity rules and the registry of the suite runner) |
| `unit_single` | atoms, cif2mrc, C-alpha finder, species |
| `unit_stream` | preprocessing, optics assignment, initial analysis, reference picking, particle sieving, pool 2D, pool 2D object, solve 3D (the steps of the stream stages, without a queue or waits), stream master, job sets, meta plots, stream watcher |

`simple_test_exec test=units` runs all fifteen in one process. It is a
convenience and deliberately not a CTest entry.

### 1.2 Long-running CTest entries

| entry | label | what it runs |
|---|---|---|
| `lib_reconstruction` | library | PCG half-set: independent half-set PCG solves against gridding; PCG sparse sampling: an atomic phantom from a tilt series, a preferred orientation and a dense control, PCG at two iterations against gridding; addon report docking: the solve3D_addon report docks a map rotated by 90 degrees and recovers the rotation; truth map docking: a map docked in both hands onto a truth made from its mirror, and every particle's pose and shift carried into the truth's frame; fractional reconstruction |
| `lib_cart_align3D` | library | pose 1JYX recovery: 1 000 simulated 1JYX particles perturbed by 15 degrees and 2 pixels, refined by the Cartesian 3D strategy (`strategy3D_cont`) under `objfun=cc` and `euclid` at 8 and 4 A, gated on the ground-truth orientations |
| `lib_heterogeneity` | library | flex PCA on a two-state phantom, flex PCA deconvolution of 20 000 particles, the PCG operator at box 64, the PCG solve sweep |
| `lib_single` | library | pdb2mrc coverage of the built-in molecular models |
| `lib_stream` | library | optics assignment, picking references, pick and extract, solve 3D reference docking, and the stream chains sieve to 3D and movies to 3D (they submit jobs to the local queue through the installed executables and take tens of minutes) |
| `pcg_recon` | library | gated stages of the PCG reconstruction operator on generated data, each stage's decision a row of its table; in process, so a library entry (it was a high-level one until 2026-10-18) |
| `pcg_frac_update` | library | the PCG fractional-update gates (additivity of the reconstruction sums, the continuation weighting, full-mass preservation, replay) on a project of 600 blob particles at their true poses that the multi-state fixture builds |
| `mini_stream_6vxx`, `mini_stream_1jxy` | highlevel | independent embedded-model mini-stream validations |
| `simulated_workflow_6vxx`, `simulated_workflow_1jxy` | highlevel | simulated particles in movies (6vxx: 160 particles in 16 movies from 100 reprojections; 1jxy: 100 particles in 10 movies) through import, motion correction, CTF, picking, extraction, `solve2D`, `solve3D`; 1jxy then runs `solve3D` a second time from the same `solve2D` project with `nparts=2` (two local worker processes, the distributed strategy), in a directory of its own and checked against the same truth, so that the master's job scheduling and assembly between parts are covered (owner decision 2026-10-15); the test asserts that the distributed run wrote its second part's record. A control scores the result against the other molecule's truth. The map's correlation with its truth and its FSC are reported as a robustness profile: in the profile of 2026-10-19 `solve3D` missed the map, at the correlation of the other molecule, in 2 of 20 runs for 6vxx and in 1 of 20 for 1jxy in each of its two steps, so no floor is admissible. Both suites reseed the whole workflow with the root seed 20260923, printed in the log |
| `single_workflow_fcc`, `single_workflow_wurtzite` | highlevel | independent SINGLE pipelines on simulated FCC Pt and wurtzite CdSe nanoparticles; the poses are judged by their median rotation error against the truth after one docking of the `autorefine3D_nano` map, up to the rotations of the lattice's group (O for FCC, D6 for wurtzite), since these near-symmetric particles are posed only up to them |
| `simulate_particles` | highlevel | `reproject` and `simulate_particles` on the embedded 6VXX volume; `simulate_particles` with three blob volumes and given populations (exact state counts, every image closest to the projection of its own state's volume at its true pose); 50 images without CTF, each correlating better with the reprojection at its own recorded pose than with one at another image's pose |
| `state_fixture` | highlevel | the multi-state test fixture (`simple_test_state_fixture`): blob phantoms in its search, known-pose and exact-view modes (built twice, the second from the cache) and the compositional and conformational beta-galactosidase states, each truth holding its declared point group and no larger one (its correlation over the group at least 0.99 and at least 0.01 above its correlation over the next larger group) |
| `solve3D_addon` | highlevel | `solve3D` on a seeded selection of a first set of simulated particles of a symmetry-broken 6VXX map (a 2000-row frozen project), then `solve3D_addon` on a 3000-row current project that appends a second set, checked against the simulation truth (poses, map) and the base run; its own entry by owner decision (2026-09-26) |
| `cont_refine3D_1jxy` | highlevel | continuous Cartesian pose refinement on 1 000 simulated 1JYX particles from perturbed poses: `refine3D refine=cont` (euclid, cc, distributed), `refine3D_auto pose_cont=only`, the polish `pose_cont=yes` in `refine3D` and `refine3D_auto`, and the continuation of a polar `refine3D_auto`, gated on the ground-truth orientations |
| `multistate_assignment` | highlevel | `refine3D nstates=3` with `refine=prob` and with `refine=prob_neigh` on 1 500 simulated particles of the three blob phantoms of the state fixture (populations 50, 30 and 20 %, SNR 2), from the truths low-passed to 20 A and random initial states, scored by correlation (`objfun=cc`, as production's state assignment), gated on the known state of every particle (adjusted Rand index, accuracy, each state's recall, the populations, no estimated state empty: a merge of two states fails), the maps (their mapping agreeing with the label mapping, their fit, their discrimination inside the difference mask) and the poses of the correctly assigned particles; the shuffled-label and duplicated-maps controls in both refinements; `suite=identical` runs the identical-states control by hand (labels only: identical truths leave nothing to discriminate) |
| `multistate_workflows` | highlevel | the multi-state workflows on 3 000 simulated particles of three beta-galactosidase states (the tetramer, the tetramer less chain D, the A-B dimer; populations 50, 30 and 20 %, SNR 1), in four routes: `cavgs`, `solve3D_cavgs` on 600 simulated images of the three states (class-average stand-ins: even and odd halves that are independent noise realisations of the same state, view and shift, each at SNR 2.5, no CTF, shifts within a pixel) imported with `import_cavgs` (owner ruling 2026-10-17); `shc`, `solve3D` from the particles; `states`, `refine3D_states` from the true poses and the consensus map; `consensus`, single-state `solve3D` in C1 followed by `refine3D_states`, the route used on real data (owner ruling 2026-10-18). Each is judged against the simulation truth: the path each refinement took (refine mode, directions and states per iteration, from the refine3D strategy's record), the state labels, the maps docked in both hands, their discrimination and the poses, with the shuffled-label and duplicated-maps controls. The de novo routes `shc` and `cavgs` resolve every true state in only part of the runs (13 and 9 of 20 in the profile of 2026-10-19), a limitation of the method: whether they do is a `success` row, held by `test=profile` to 6 and 2 of 20 seeds, and their labels are reported; their path, checks and pose ceilings are gates. The `states` and `consensus` routes resolve every state in every run: that is checked, and their labels, populations, discrimination and poses are held to floors. The entry runs all four (1 913 s at eight threads alone); `suite=<route>` runs one, `suite=cavgs_single` (the single-state class-average start, `prob_athres=` passed on), `suite=cavgs_as_particles` (the same simulated halves through particle `solve3D`), `suite=<route>_identical` (identical states) and `suite=hard_<route>` (a chain rotated by 10 degrees, SNR 0.5, populations 80 and 20 %) are run by hand |
| `single_atoms_stats` | highlevel | simulated Pt nanoparticle atom detection and statistics |
| `single_atoms_stats_cdsew`, `single_atoms_stats_cdsez`, `single_atoms_stats_cdser` | highlevel | the same test on simulated wurtzite, zincblende and rocksalt CdSe nanoparticles (`element=CdSeW`, `CdSeZ`, `CdSeR`), whose lattice lookups use the compound selector; also asserts the coordination number of the interior atoms (4, 4, 6) and the fitted lattice constant within 2% of the table |
| `species_discovery` | highlevel | `detect_atoms` with a species list, with `discover_species=yes` and without an element, on simulated Pt, Pt/Ni and Pt/Al particles rendered with the five-Gaussian element kernels and surface B factors (with and without half maps), judged against the generating models through `test_gate` (`metrics.tsv`); the present products must not change |
| `stream_preproc` | highlevel | five simulated movies through the stream's preprocessing stage and its worker jobs |
| `flex_pca_blobs` | highlevel | `flex_pca` on a two-state blob dataset of the multi-state fixture (2 000 particles at their true poses, SNR 1) in shared memory and with two local workers, checked against the truth and against each other; `refine3D_auto state=X` on the truth-matched states, with the frozen weights and with the hard labels; `cls_expansion` and `refine3D_states_cavgs` on the fixture's three compositional beta-galactosidase states (3 600 particles in 60 exact views; ruling of 2026-10-18: not blob states that differ in a lobe of 6 % of the mass): the route's contract gated (projects, hand-back, maps), its state accuracy against the planted states (one-to-one mapping, with the shuffled-label null), the number of planted states its maps match and the poses of the sub-averages reported as a robustness profile, since its profile admits no floor (accuracy at the chance level) |

### 1.3 Programs that are not registered

A few programs are reachable through `simple_test_exec` but are not CTest
entries: `rec3D_backends` (gridding against PCG on a user's project; a
measurement without thresholds, a developer tool in the sense of section 3.4),
`cont_refine3D_bench` (the `refine=cont` benchmark on simulated 1JYX particles:
timings, work counts and pose errors, no thresholds; owner decision
2026-10-19, `doc/implementation_notes/completed/cont_lm_acceleration.md`),
`polar_matching_bench` (the layouts of the dense probability table's euclid
scoring and production's cc scoring on synthetic polar data at production
sizes: timings and identity checks, including a dump of the production kernels'
distance vectors for comparing builds, no thresholds;
`doc/implementation_notes/planned/polar_matching_layout_bench.md`),
and the input generators `solve3D_addon_snapshots`, `write_state_weights_labels`
and `write_state_weights_mixed` (tools for manual runs). `pcg_frac_update`,
which needed a user's project, builds its own from the multi-state test fixture
when run without `projfile=` and is the library entry of that name. No new
program joins this list without a reason: a diagnostic that runs on a user's
data is a developer program (section 3.4), not a test.

For any suite-based entry, print its accepted sub-suite identifiers without
running tests using `simple_test_exec test=<entry> suite=list`; for example,
`simple_test_exec test=unit_image suite=list`.

### 1.4 Running tests

```text
./compile_debug.sh                                      # build, run the fast gate, install
simple_test_exec test=list                              # every test program
simple_test_exec test=unit_image                        # one area suite
simple_test_exec test=unit_image suite=masks            # one sub-suite of it
SIMPLE_UNIT_ORDER=reverse simple_test_exec test=units   # the sub-suites in reverse order
cd build && ctest -L fast --output-on-failure           # the gate by hand
cd build && ctest -L library                            # a nightly tier by hand
cd build && ctest -L highlevel --output-on-failure      # long tests, explicit only
./run_tests.sh [fast|library|highlevel|platform|all|list|<entry>]   # the same from an existing build
simple_test_exec test=simulated_workflow suite=6vxx seed_offset=3   # the same test on other draws
simple_test_exec test=nightly                           # the whole high-level suite (section 1.5)
```

`./run_tests.sh` reruns tests from an existing build without configuring or
compiling; it runs the fast gate when given no mode, and takes
`--build-dir <dir>`. `make check` in the build directory runs the same fast gate
as the compile scripts. The day-to-day loop in an existing `build/` is
`make -C build -j<n>`, `scripts/run_fast_gate.sh build`, then
`make -C build -j<n> install`; the compile scripts run that sequence from a
fresh directory. `simple_test_exec` parses its command line against the
generated argument vocabulary, so a key that only the front door reads
(`seed_offset`, `target`, `nseeds`, `dir1`, `dir2`, `keep_outputs`) is also
listed in `scripts/simple_args_generator.pl`; a new such key not listed there
is rejected by the parser.

The `suite=` selector is the sub-suite name in lower case, with blanks and
hyphens written as underscores and commas and slashes dropped
(`stack I/O` is `stack_io`, `search, sort, locate` is `search_sort_locate`).
The selectors of each suite are listed in the description of its `suite=`
input in the test UI (`src/main/ui/simple_test/simple_test_ui_class.f90`). A
run by hand draws what a CTest run draws: `simple_test_exec` sets
`SIMPLE_SEED=20260923` for its own process first (`fix_test_seed`), whatever
the environment holds. The seed is not a user setting yet: production programs
read `SIMPLE_SEED` from the environment, and a seed key of `parameters` is
planned (`doc/refactoring_notes/planned/rnd_seed_parameter.md`).

`seed_offset=<n>` moves every fixed seed of the test by n times a large odd
stride (2654435761, in 64-bit integers), so that the offsets of one test cannot
land on another test's base. It moves the process seed of `fix_test_seed` and
every `set_fixed_seed` of a tester or high-level test (and so the `SIMPLE_SEED`
its commanders and workers inherit). 0, the
default, draws what CTest draws. It is a key of the front door only: it is read
there and removed from the command line, so no commander or `parameters` sees
it. A test must draw from fixed seeds only through these two routines, or its
draws do not move with the offset. Production code that seeds itself with a
fixed number (the class-expansion strategy and kernel of `flex_pca`) does not
move with it.

`SIMPLE_UNIT_ORDER=reverse` is a diagnostic for the fused, in-process suites:
it reverses the sub-suite table so a result that differs from the normal run
exposes state leaked by one sub-suite into another. It does not reverse CTest
entries and is not part of the normal build gate; use it after changing shared
state, lifecycle or test-runner code.

### 1.5 The suite runner: the registry, `test=nightly`, `test=profile`, `test=identity`

The high-level entries are listed once, in the test registry
(`src/main/commanders/test/simple_test_registry.f90`): per entry its command
line, its `OMP_NUM_THREADS`, the worker processes and threads it starts, its
peak CPU reservation (the cores in use at its busiest moment, master and
workers together), its deadline (the CTest `TIMEOUT`), its kind (a correctness
gate or a robustness profile only), its peak resident memory and the disk its
run directory grows to (each from a measurement, with a margin), whether its
master's threads follow a profile's thread count, and the outputs an identity
check compares: declared for the one-thread workflows and for every entry that
runs no `solve3D` or `refine3D`. `scripts/check_test_registry.py`
holds the `highlevel` registrations of `production/CMakeLists.txt` to the
registry: the same names, command lines, timeouts and thread counts, and a
CTest `PROCESSORS` equal to the CPU reservation, so that `ctest -L highlevel -j
24` packs the entries as `test=nightly` does. A new high-level entry gets its
registry row and its CTest registration together.

Three test programs run and compare the suite; none of them is a CTest entry.
They run every test execution as a child process of the test dispatcher
(`simple_test_dispatch`), which:

- forks the child (`simple_forked_process`), which leads its own session,
  enters its own run directory, writes its output to `run.log` there, sets its
  `OMP_NUM_THREADS` and seed offset and runs the test program in the process;
- gives every child the environment CTest gives an entry: `SIMPLE_PATH` the
  build the running `simple_test_exec` belongs to, its `bin/` and `scripts/`
  first on `PATH`, `SIMPLE_QSYS=local`, so that the child and the workers it
  starts run the same build;
- starts the largest reservation queued (the first in queue order among
  equals) when its cores, memory and disk fit in what is free: the cores of
  the budget (`ncunits=`, by default the processors of the machine) less those
  of the running executions, 90 % of the memory available when the run starts
  less the running executions' memory, the disk free on the run's file system
  when it starts less the disk of every execution started (outputs stay).
  Smaller executions wait for every larger one to start, so they cannot starve
  it; one that does not fit even alone runs alone;
- removes a `metrics.tsv` left in a reused run directory before the child
  starts, so a table is always the execution's own;
- records a failed child with its exit status and goes on; sends a child past
  its deadline SIGTERM and, after a grace period (60 s), SIGKILL, and records
  it as timed out; kills every process left of a child when it ends, the
  child's process group before the child is reaped (so that its pid and group
  cannot belong to another process yet). The processes of a child are its
  session and every process that inherited its marker `SIMPLE_TEST_EXECUTION`
  from the environment, since the local queue starts workers in sessions of
  their own. A child ends at once on SIGTERM (`_exit`, async-signal-safe); its
  gate rows and log lines are flushed as they are written, so its table holds
  every row written before the signal;
- samples the CPU use and resident memory of each child's processes every
  second and reports the peak;
- when interrupted (SIGINT, SIGTERM, SIGHUP), stops every running child,
  records it as interrupted (an outcome without a verdict, outside the failure
  rates), starts no other and reports the rest as not run. A disposition the
  dispatcher inherited as ignored (`nohup`, a background job) stays ignored. A
  child dies with a dispatcher killed outright (SIGKILL); the workers of such
  a child may then outlive it.
- **Linux only**: the markers of processes outside a session, the CPU and
  memory samples and the count of processes left read `/proc`, and the exit
  of a child is detected without reaping it by `waitid` with `WNOWAIT`.
  Elsewhere the dispatcher says so when it starts, signals the session alone
  and kills after reaping; the dispatcher's own tests skip the probes that
  need `/proc`.
- **Known limits.** Largest-first packing starts the long one-thread
  workflows (`simulated_workflow_*`, reservations of 2 and 3) after every
  larger reservation, so on 24 cores they set the nightly's wall time (about
  an hour on 2026-10-18 against 45 minutes of CTest packing); ordering by
  expected core-seconds would start them early. The disk budget counts every
  started execution's disk, since outputs stay; a profile that deletes
  the outputs of passed executions as they end therefore runs one execution at
  a time once the budget is spent.

```text
simple_test_exec test=nightly [ncunits=24] [dir=<run>]
simple_test_exec test=profile target=<entry[,entry...]|all> nseeds=20 [nthr=1,8] [suite=<sub-suite>] [ncunits=24] [dir=<profile>]
simple_test_exec test=profile dir1=<profile> dir2=<profile>
simple_test_exec test=identity target=<entry> dir1=<run directory> dir2=<run directory>
```

- **`test=nightly`** runs every correctness entry of the registry once, each
  in `<run>/<entry>`, and writes `nightly_report.tsv`: one row per entry and
  gate metric with the outcome (passed, failed, timed out), the exit status,
  the wall time, the peak CPU use beside the reservation, the peak memory, and
  the metric's value, floor, verdict and sense. It fails when an entry or a
  gate failed. It runs at the entries' CTest thread counts and seed offset 0,
  so an entry run by it draws what CTest draws.
- **`test=profile`** runs one entry, a comma-separated list of them or all at seed offsets 0 to nseeds-1 and
  each thread count of `nthr=` (by default the entry's own), each execution in
  `<profile>/<entry>_t<threads>_s<offset>`, and writes `executions.tsv` (as the
  nightly report) and per entry and thread count `profile_<entry>_t<n>.tsv`:
  for every gate metric the runs that wrote it, the median, the robust spread
  `s` (1.4826 times the median absolute deviation, or the range over 4 when
  that is zero), the range, the worst run in the sense of its floor, and the
  failure rate: an execution fails a metric when one of its rows of it fails,
  when it wrote none, or when the execution itself failed or timed out (even
  after writing passing rows); interrupted and unstarted executions are left
  out. `nseeds=` is 20 when not given. `suite=` names one entry's sub-suite, so
  it needs a single `target=`. A thread count other than the entry's own scales
  the `nthr=` the entry passes to its commanders and, below the entry's own,
  its deadline (inversely); the reservation follows the master's share (the
  reservation less its workers' threads) and each worker's threads in
  proportion, while an entry whose threads are its own (the one-thread
  workflows, the mini stream, the stream preprocessing) keeps its reservation;
  memory and disk do not scale. Seed offsets move every base seed of the test
  by offset times a large odd stride (2654435761, in 64-bit integers), so the
  offsets of one test cannot land on another test's base. A profile informs;
  it does not fail on a metric, with one exception: a success metric (sense
  `success`, below) counts the executions that succeeded (the column
  `successes`), and the profile fails when that count is below its floor times
  the executions. `keep_outputs=no` deletes the subdirectories of the run
  directory of each passed execution (its outputs) as soon as it ends, keeping
  `run.log` and `metrics.tsv`, and returns its disk to the budget, so that a
  long profile counts only the running executions' disk; the outputs of a
  failed or timed out execution stay. With `dir1=` and
  `dir2=` it compares two profiles instead: for each table both hold,
  `comparison_profile_<entry>_t<n>.tsv` gives both medians, the shift of the
  median in units of the larger spread (positive is better) and both failure
  rates.
- **`test=identity`** compares two run directories of one entry byte for byte
  over the outputs its registry row declares, and writes `identity_report.txt`
  (SAME, DIFF or MISSING per file). It ignores only the declared run-specific
  fields, which hold the run's paths, run identifier or time:
  - in `solve3D_manifest.txt`, the lines `run_id`, `layout_digest` and
    `stack_digest` (digests of the absolute stack paths), `checksum` (of the
    whole manifest), and the path that ends each `artifact` line;
  - in `sigma2_state.bin`, the layout digest (header word 14, bytes 121-128);
  - project files (`*.simple`) are compared record by record over every
    segment, and text files (`*.txt`) line by line, each with its run
    directory's own path removed, since both hold absolute paths.
  A declared pattern that selects no file in either directory fails. The
  `BENCH` timing files are not declared outputs. Identity is decided only by
  one-thread runs. An entry that works in a directory named by its process id
  (`test_<name>_<pid>`: the mini stream, `simulate_particles`, the atom
  statistics, the stream preprocessing) declares its outputs relative to that
  directory, and `dir1=` and `dir2=` name it; `state_fixture` and
  `species_discovery` keep no files, so their gate table (`metrics.tsv`) is
  what they declare.

## 2. What a test is

A test states something the code guarantees and fails when that guarantee is
broken. If it cannot fail, it is not a test, whatever it is called and
wherever it lives.

### 2.1 A test

- **Makes assertions whose failure reaches the exit status.** Checks go
  through `simple_test_utils` (section 5.1), one assertion per guarantee, with
  a message that states the guarantee ("an integer shift is the circular shift
  out(x) = in(x + s)", not "shift test 3").
- **Derives the expected value independently of the code under test**: a
  closed form, a brute-force computation in the test, an exact array
  operation, or an emulation of the algorithm (for example in numpy) whose
  numbers are written into the test with a comment naming how they were
  obtained. Running the routine once and pasting its output into the test
  pins the bug along with the behaviour.
- **Names its tolerances** and derives them from the arithmetic: single
  precision, an interpolation order, a statistical standard error. A
  tolerance that was widened until the test passed is a finding, not a
  tolerance.
- **Pins what the routine guarantees**, not what it happened to produce: not
  the outcome of a random start, an unconverged iteration or a timing.
- **Includes a negative control where one is cheap**: an unrelated image does
  not correlate, a constant sample has no foreground, an atom's density does
  not correlate with a window beside it.
- **Pins conventions.** Sign, handedness, origin and index conventions are
  where SIMPLE's code breaks silently. When a test finds out which convention
  a routine follows, it asserts that convention.
- **Is reproducible**: every draw it makes comes from a fixed seed.

The review of September 2026 wrote tests to this standard for most of the
library, and they found and fixed more than twenty-five production defects,
none of which the programs they replaced could have caught. Examples of the
standard are the graphene-mask test (the excluded shells found by a
brute-force search in the test), the `rotate_ref_8` test (the explicit
rotation of `calc_frc` must peak where the FFT path of `gen_objfun_vals`
does, for probes that exercise every branch) and the `atom_validate` test (an
atom must correlate with its own simulated density; a numpy emulation put the
expected value at 0.99 and the defect at 0.4).

### 2.2 Not a test

- **A program that only runs.** It prints numbers, writes files or reports
  "completed". Of the 130 test identities that existed before the
  refactoring, 63 could not fail.
- **A wrapper of commanders.** A routine that calls `new_project`,
  `import_movies` and `preprocess` in sequence is a high-level workflow. It
  belongs with the production programs, with developer visibility, or in a
  local worktree (section 3.4). It becomes a workflow gate (section 6) only when it runs on
  simulated data and compares its result with the model the data was
  simulated from, against declared floors.
- **A benchmark or a timing loop.** Timings are not assertions; a benchmark
  is a developer program or stays in a local worktree.
- **A conversion or utility invocation that only produces output**: a program. The underlying
  conversion can be a unit test when it generates its own input and quantitatively verifies the output.
- **An experiment on downloaded or user data**: a local worktree, or a
  developer program if the team needs it again.
- **A demonstration** (sampling pictures, gnuplot windows, printed tables to
  eyeball): delete it or keep it in a local worktree.

### 2.3 Four kinds of test, kept apart

A run of SIMPLE with more than one thread is not reproducible, and a workflow
on simulated data is a draw: its metric scatters from seed to seed. A test
that mixes "is it right" with "how reliably" cannot say which of the two
failed. So there are four kinds of test, and each answers one question:

| kind | question | how it is run | how it passes |
|---|---|---|---|
| **identity** | did a refactoring change any result? | one thread, fixed seed, the same case on two builds (`test=identity`) | the declared outputs byte-identical, run-specific fields excluded by a declared list |
| **correctness** | does the code recover the truth it was given? | one run on simulated data with a large, known signal (CTest, `test=nightly`) | truth metrics against floors that sit far from both the observed distribution and the negative control |
| **robustness** | how reliably, and how well, across draws? | the same case over many seeds and thread settings (`test=profile`) | a distribution (median, spread, worst, failure rate) recorded and compared over time; never a single-run gate |
| **performance** | how fast? | benchmark reports | timings, never in a pass/fail test |

Consequences:

- **A correctness gate is insensitive to the draw.** Its case is built so that
  the right algorithm passes on every seed and thread count by a wide margin,
  and its floors come from a robustness profile by the floor rule (section
  6.2). If no floor is admissible, the case is not a correctness gate: its
  metric is reported as a robustness profile, or its signal is raised until it
  can be a gate.
- **Identity is the only test of a refactoring that claims identical
  results.** It never uses multi-threaded runs.
- **Robustness profiles inform decisions; they do not stop work.** A change
  meant to change results is judged by comparing profiles before and after,
  on cases with known truth.
- **Multi-state correctness is judged only against known labels on simulated
  data.** Real data without labels shows stability, never correctness
  (section 6.6).

## 3. Where a test goes

### 3.1 The decision

Answer these in order.

1. **Does it assert something (section 2.1)?** If not, it does not go in the
   test environment at all; see section 3.4.
2. **Which area does it test?** A test goes to the area of the code it
   tests, not the area of its author or of the workflow that uses the code.
   The areas are named after the machinery and kept short: `core`, `ori`,
   `image`, `numerics`, `project`, `motion`, `ui`, `ipc`, `reconstruction`,
   `pftc_align2D3D` (everything on the polar Fourier transform, 2D and 3D),
   `cart_align3D` (Cartesian continuous registration), `heterogeneity`,
   `parallel`, `single`, and `stream` for the stream stages.
3. **Which tier?** The first that fits:
   - **fast** (`unit_<area>`) when it meets all of the admission rules of
     section 3.2. This is the default: most tests of the library belong here.
   - **library** (`lib_<area>`) when it is hermetic and deterministic and runs
     in-process, but needs realistic sizes and minutes.
   - **highlevel** when it runs a pipeline that starts `simple_exec` or
     distributed workers, on simulated data, and is gated against the
     simulation truth (section 6).
   - **platform** when it needs child processes, a launcher (`cafrun`,
     `mpirun`) or a device.

### 3.2 Admission to the fast gate

All of them:

1. It can fail (section 2.1).
2. It is hermetic: no network, no download, no user-supplied file, no
   dependence on core count, load, wall clock or another test's files.
   Fixtures are generated from a fixed seed, or committed and small.
3. It runs in-process on one OpenMP thread, returns on success, restores the
   working directory, and leaves no global state behind (the random generator
   is reseeded by the runner; module memos such as the mask coordinates must
   be released: `unmemoize_mask_coords`). A sub-suite whose subject is a
   threaded path opens its own small team with `num_threads` (three is the
   convention); the entry's `OMP_NUM_THREADS` stays 1.
4. It is cheap. The whole gate must stay under 60 s of real time; keep every
   sub-suite well under a second, with the smallest fixture that exercises
   the path (a 64-pixel box tests an FFT as well as a 256-pixel one).
5. It uses no sleeps, polls or real timers.

### 3.3 Is it part of the build?

The fast tier is, and nothing else is. A test that does not meet section 3.2
must not be forced into the gate by shrinking it until it no longer tests
anything; it goes to a library suite and runs every night. A library test
that turns out to take seconds can move to the fast gate once its fixture
fits.

### 3.4 Scratch development

The test area is not a scratch area. When a test directory holds
"everything that is not production", it fills with programs that nobody runs,
owns or can delete, which is what the old `production/tests` became.
Unfinished work has three homes:

- **A local worktree.** Exploration, prototypes, experiments and throwaway
  diagnostics live in a git worktree on your own computer until they become
  something (see below).
- **A developer program.** A diagnostic or tool the team will run again, on
  real data, is a commander with a program in the UI at developer
  visibility (the default of `ui_program`). It gets parameter parsing,
  project I/O, logging and distributed execution for free, which a
  standalone program had to reinvent. It needs an owner and a one-line
  purpose, and it goes when it is no longer used.
- **A test.** Once the exploration settles, the claims that must stay true
  are written as assertions in a tester module (section 4).

What lands on master is production code, a developer program or a test.
Nothing else.

**Work on master; scratch in a worktree.** SIMPLE practises extreme
programming: everyone works on master and pushes to master, in small and
frequent steps, so that everyone's work is integrated, built and tested
every day. What makes this safe is the fast gate, and running it is up to
you: **before you push, build with `./compile_debug.sh` and
let the gate pass.** This is a strong recommendation, not a lock: nothing in
Git or CI blocks a push, because sometimes a push only moves code to another
machine (a cluster node, say). Then run the gate before the work counts as
done. The nightly GitHub Actions workflows build and gate master on Linux and
macOS using clean Release, Debug, NICE GUI, and coarray configurations, and
they run the broader platform, library, and high-level CTest labels. A separate
benchmark-dataset workflow runs short real-data subsets and publishes its FSC
and resolution gates. The README badges and linked workflow/report pages show
the next morning what slipped through. Do not keep work on a branch in the
online repository: nobody checks it, it drifts away from master, and the
merge gets harder every day it waits. A branch pushed to the online
repository is the exception and needs a strong reason (for example a
coordinated change that cannot build on master for a while); it is merged
back as soon as it can be.

A git worktree is the tool for scratch development. It is a branch with its
own physical copy of the source on your computer, next to your main
checkout, with its own build, so you can try several things side by side
without stashing or switching, and evaluate them against master. When you
are happy with the result, merge it into master and push; when it led
nowhere, remove it.

```sh
# a scratch branch with its own checkout beside the main one
git worktree add -b scratch/nu_probe ../SIMPLE-nu_probe master
cd ../SIMPLE-nu_probe && ./compile_debug.sh
# ... develop, build, run the tests ...
# what should stay: merge into master in the main checkout and push
cd ../SIMPLE && git pull && git merge scratch/nu_probe && git push
# clean up (git branch -D for a branch that is not merged)
git worktree remove ../SIMPLE-nu_probe && git branch -d scratch/nu_probe
```

`git worktree list` shows the worktrees you have. Scratch branches stay
local: they are not pushed.

## 4. How to implement a test

### 4.1 A unit test (fast or library tier)

A unit test is a subroutine in a tester module next to the code it tests.

- **The file** is `simple_<thing>_tester.f90` in the directory of the
  production module, named after what it tests (`simple_stack_io_tester.f90`
  beside `simple_stack_io.f90`). Files named `*_tester.f90` are compiled only
  with `BUILD_TESTS=ON`.
- **The header** is a single `!@descr:` line naming what the module tests and
  the production modules it covers; more description goes below it in plain
  `!` lines (`scripts/check_descr.py`, section 4.5).
- **The module** uses `simple_test_utils`, is `private`, and exports one
  entry, `run_all_<thing>_tests`, which prints
  `**** running all <thing> tests ****` and calls each test. Each
  `test_<what>` subroutine prints its own name first.
- **No process control in a tester**: no `stop`, `error stop`, `simple_end` or
  `report_summary`. The suite runner reports and sets the exit status.
  `THROW_HARD` is only for a broken fixture (a file that could not be
  written), never for a failed check.
- **Registration**: import `run_all_<thing>_tests` in
  `src/main/commanders/test/simple_commanders_test_class.f90` and add
  `call add_suite(s, n, '<sub-suite name>', run_all_<thing>_tests)` to the
  `suites_<area>` (or `suites_lib_<area>`) table. Names are at most 32
  characters. Add the selector to the `suite=` list of the area's program in
  `src/main/ui/simple_test/simple_test_ui_class.f90`, in table order; the
  registry check fails the build until the two agree.
- **A new area suite** (a new `unit_<area>` or `lib_<area>`) also needs its
  program in the test UI, a router case and a CTest registration, and it
  raises the process budget. That is an owner decision.

### 4.2 Fixtures

- Generate them in the test: images from `gauran`, `gauimg`, `square`;
  molecules from the embedded models (`simple_molecule_data`, `pdb2mrc`);
  particles and movies from the simulation commanders.
- Seed every draw with `set_fixed_seed(<n>)`. Call it **after**
  `parameters%new`, which reseeds.
- Keep them small and in memory. Write files only when file I/O is what the
  test is about, in a directory made by `enter_fixture(tag, ...)` and left by
  `leave_fixture(...)` (`simple_test_utils`): the directory is `<tag>_<pid>`
  under the working directory; it is removed when no check failed since the
  test entered it, and kept, with a message, when one did, so that a failing
  run can be inspected. Files outside a fixture directory are named
  `tmp_<tester>_*` and deleted at the end of the test, also when a check
  failed.
- A tester that loops over something printing through `logfhandle` diverts
  `logfhandle` to a scratch unit and restores it.

### 4.3 White-box self-tests

A self-test that must read private components of a type goes in a submodule
of the production module, in a file named `*_tester.f90` (so it is built only
with `BUILD_TESTS=ON`); the module declares its interface. It takes an
out-argument per check, and a tester asserts those by name
(`test_flex_pcg_operator`, declared in `simple_flex_pca_pcg` and implemented
in its submodule `simple_flex_pca_pcg_tester`, is the model). Every other
self-test belongs in a tester module: the review moved the self-tests of
`image`, `imghead`, `atoms`, `oris`, `ftiter`, `ftexp_shsrch`,
`bspline_smoother`, `online_var` and `aff_prop`, the flex PCA, UI hash,
cavg-quality and class-average registration self-tests out of their
production modules, and deleted those of unused code (`hclust`,
`srchspace_map2D_io`) and the dead ones (`CPlot2D`, `jpg`). No self-test is
left in a production module.

### 4.4 A workflow or platform test

A test that runs a pipeline is a commander: its type and `exec_test_<name>`
routine in the file of its family in `src/main/commanders/test/`
(`simple_commanders_test_highlevel.f90` and its `..._highlevel_<topic>.f90`
submodules, `..._single.f90`, `..._stream.f90`; the suites and platform tests
are in `..._class.f90`), a case in the router of the same family,
`src/main/exec/simple_test_exec_{highlevel,single,stream,class}.f90`, a program
in `src/main/ui/simple_test/simple_test_ui_highlevel.f90` (the suites and
platform tests are in `simple_test_ui_class.f90`), a registration in
`production/CMakeLists.txt` with `simple_add_test` (label, timeout, threads,
and `RUN_SERIAL` for a library or platform entry, `PROCESSORS <reservation>`
for a high-level one) and, for a high-level entry, its row in the test
registry (section 1.5) and its card (section 6.5). Section 6 says what it has to check. It still fails through
assertions or `THROW_HARD` with a message that names the metric, its value
and the floor.

### 4.5 Checks and records

Before asking for a build:

- `python3 scripts/check_descr.py .` — every source file starts with a
  one-line `!@descr:` and declares at most one module or submodule.
- `python3 scripts/check_test_registry.py . --verbose` — the CTest
  registrations, the test UI programs, the router cases and the `suite=`
  lists agree. The fast gate runs it on every build.
- `python3 scripts/check_flex_dag.py . --verbose` — each FLEX source declares
  one module or submodule, and the production graph follows its declared layers
  and remains acyclic. The fast gate runs it on every build.
- The default CMake build regenerates the test inventory and code-base map
  when their source inputs change.
- `./compile_debug.sh` — the build and the gate.
- `SIMPLE_UNIT_ORDER=reverse` on the area suite — the sub-suites do not
  depend on each other's state.

When a test is deleted, the commit says why. When a test finds a defect, the
defect is fixed in the same change and the test pins the fix.

Two kinds of code need a check that is easy to leave out:

- **Threaded paths.** The fast gate runs every suite on one OpenMP thread, so
  a routine with a threaded path gets a check that forces a team (the masks
  suite runs a team of three) and compares with the serial result. Nothing
  repeats the fast tier on more threads, so that check is the only one.
- **Readers and writers.** A file format is tested by round trips at the
  sizes where its layout changes: small, odd and non-square boxes and the
  sizes where a header or record length changes. The SPIDER header was wrong
  for every box below 43 and right for the 64-pixel boxes the tests used.

### 4.6 Fortran traps seen in the review

- A dummy argument or local named like a module constant hides it, silently
  (Fortran is case-insensitive: `nmics` hides `NMICS`).
- `set_fixed_seed` before `parameters%new` is overwritten by it.
- A constant integer division that truncates, `(BOX-1)/2`, draws a
  `-Winteger-division` warning; write the value.
- Comparing whole images after an interpolating rotation includes the
  circular wrap of `rtsq` at the corners; compare inside the inscribed disc.
- An `rmat` pointer (`get_rmat_ptr`) is the box, `ldim(1) x ldim(2) x ldim(3)`
  (since 2026-09-30): a whole-array expression on it is an expression on the
  image and conforms with a box-sized array, and `get_rmat()` returns a copy
  of the same box. Until then the pointer was the padded buffer of the
  in-place FFT, with extra rows in the first dimension; an unbounded
  whole-array expression on it did not conform with a box-sized array, which
  bounds checking stops and which passes by luck without it (the `cif2mrc`
  tester, 2026-09-25), and a whole-array reduction took the padding in. What
  is left of the trap: the box is not contiguous in memory, so pass the
  pointer to assumed-shape or pointer dummies (an explicit-shape dummy gets a
  copy), and an index past `ldim(1)` is now out of bounds. In real space the
  padding holds zeros; a test reads it through `max_abs_padding()`, never
  through `get_rmat_ptr_padded`, which is for the callers that hand Fourier
  data on in the FFT layout (`simple_image_slab`, `simple_stack_read_session`).

## 5. Library support for tests

### 5.1 `simple_test_utils`

The assertion and reporting layer. A tester uses the assertions and
`set_fixed_seed`; the suite runner uses the rest.

| routine | use |
|---|---|
| `assert_true(cond, msg)`, `assert_false(cond, msg)` | a logical guarantee |
| `assert_int(expected, actual, msg)` | an exact integer |
| `assert_real(expected, actual, tol, msg)` | a real within an absolute tolerance (`tol = 0.` for an exact value) |
| `assert_double(expected, actual, msg[, ulp_tol])` | a double within a number of ulps (10 by default) |
| `assert_char(expected, actual, msg)`, `assert_string_eq(char, string, msg)` | text, trimmed |
| `set_fixed_seed(base[, propagate])` | seed the intrinsic generator with a fixed state, base plus the seed offset; `propagate=.true.` also sets `SIMPLE_SEED`, from which the commanders run in process reseed |
| `fix_test_seed` | `SIMPLE_SEED` of the test process set to 20260923 plus the seed offset; called by `simple_test_exec` and the test dispatcher, not by testers |
| `set_test_seed_offset(n)`, `offset_seed(base)` | the seed offset of the process (`seed_offset=` of `simple_test_exec`), and a base seed moved by it (a 64-bit integer); for the front door and the dispatcher, and for a test that prints its effective seed |
| `enter_fixture(tag, cwd_saved, fixture_root)`, `leave_fixture(cwd_saved, fixture_root, nfail_before)` | the scratch directory of a test (section 4.2): made and entered, then left and removed unless a check failed since `nfail_before` |
| `begin_test_suite`, `end_test_suite`, `reset_test_report`, `report_summary`, `tests_run`, `tests_failed` | the runner's bookkeeping; not for testers |

A failed assertion records its message, with expected and actual values, and
the suite goes on: every check of a run is reported, not only the first
failure.

### 5.2 Workflow-gate support

Shared by the highlevel test commanders (`src/main/commanders/test/`, built
with the tests); a new gate uses these instead of its own copies.

| routine | use |
|---|---|
| `test_gate` (`simple_test_gate`): `new(fname)`, `check(name, ok)`, `metric(name, val, floor, ok[, upper])`, `report(name, val)`, `success(name, ok, floor)`, `passed()`, `kill` | a gate's checks and metrics, logged as PASS/FAIL and tabulated in a TSV file (name, value, floor, pass, sense); `upper=.true.` marks the floor as a ceiling; a success never fails the run, and `test=profile` holds the fraction of its runs that succeed to `floor` |
| `validate_reconstructed_volume` (`simple_test_truth_metrics`) | a workflow's final map against its simulation truth: box and sampling, docking in both hands, correlation and masked-FSC floors, each decision also a row of the gate it is given |
| `dock_both_hands`, `compare_to_truth` | the two steps of it, for a gate with its own floors (the masked FSC is the production `compare_volpair` of `simple_volpair_metrics`, which also gives a band-limited correlation) |
| `pair_pose_error` | poses against the simulation truth without docking (median error of relative rotation angles over seeded pairs) |
| `add_gaussian_blob` | a blob that breaks the symmetry of a model map, so that c1 poses are unique |
| `state_fixture` (`simple_test_state_fixture`) | a multi-state dataset with known states, poses, shifts and CTF: blob phantoms (box 64) or beta-galactosidase states (box 96 at 2.2 A; compositional: tetramer D2, tetramer less one chain C1, dimer C2; conformational: one chain rigidly rotated), one `simulate_particles` mixture with given populations, and a project per mode (search: no poses, no states; known: the true poses; views: exact repeated views, class = view); the truth is `truth.txt` beside it. Cached per build in `<SIMPLE_PATH>/test_runs/fixture_cache/<family>_<hash>`, where the hash (64-bit FNV-1a) is of a key text holding the commit of the build, the fixture version, every parameter of the dataset and its seed with the seed offset; the directory holds the key as `key.txt` and is built in a temporary directory renamed into place when complete |
| `point_group_check`, `betagal_frame` | a truth's self-correlation under its declared point group and under the extra operator of the next larger one; the D2 frame of 1JYX from its chains |
| `match_states`, `adjusted_rand_index`, `mapped_accuracy`, `accuracy_null`, `collapsed_accuracy`, `population_errors` | state labels against the truth: the one-to-one mapping that maximises agreement (at most eight states), the adjusted Rand index, the accuracy under the mapping, its null distribution from shuffled truth labels (counts kept) and the collapsed predictor |
| `match_maps`, `discrimination` | the map mapping from docked correlations (it must agree with the label mapping), and whether the maps discriminate the states inside the mask where the truths differ (margins and difference correlations; copies of one map fail) |
| `carry_pose`, `group_errors`, `group_operators` | a pose and shift carried into the truth's frame by a docking and hand (`dock_both_hands` returns them), and the rotation error under the operators of the truth's declared point group only (poses r and r g are equivalent) |

## 6. The nightly suite

### 6.1 What belongs in it

- **Library suites**: tests that meet every rule of the fast gate except its
  cost. Realistic box sizes and particle counts, minutes rather than
  seconds, still hermetic, deterministic and in-process, and they do not
  start `simple_exec` or workers. One suite per coherent part of the library.
- **Workflow gates**: simulated pipelines. They run production programs on
  data simulated from embedded atomic models (6VXX, 1JYX, the Pt
  nanoparticle), where the orientations, defocus, shifts and positions are
  known, and they compare the result with that truth.
- **Platform entries**, on a machine that has the capability.

What does not belong: anything that needs user data or a download, timing
benchmarks, and anything that cannot fail.

On the test machine the whole high-level tier runs as one command,
`simple_test_exec test=nightly` (section 1.5), which packs the entries on the
machine's cores and writes one report; `ctest -L highlevel -j <cores>` runs the
same entries packed by their `PROCESSORS`. Despite its name, `test=nightly` is
not a nightly of record with archives and notifications: it is the tier of the
tests too expensive for the build, today the simulated workflow gates (ruling
of 2026-10-18). The library, platform and fast tiers are not part of it. Real
datasets may join it as entries run from within SIMPLE once it is settled how
their data reaches the machine (section 6.6).

### 6.2 How a workflow gate is designed

- **It checks truth, not completion.** "The files exist and the counts
  match" is the floor, not the gate. The gate compares the result with the
  model that generated the data: the estimated defocus against the simulated
  one, the estimated frame shifts against the simulated drift, picked
  positions against the placed particles, the map against the truth map, the
  recovered orientations against the simulated ones.
- **The truth is independent of the code under test.** Compare the final map
  with the map `pdb2mrc` makes from the same coordinates on the same grid,
  not with a map from an earlier stage of the same run.
- **Maps are compared after docking and hand.** A `solve3D` map is
  neither docked to the truth map nor necessarily of the right hand. The gate
  docks it with `dock_vols` at a low-pass of 15 to 20 Å, keeps the hand (the
  map or its `mirror('x')`) that docks with the higher correlation, and takes
  the masked FSC at 0.143 against a declared floor. Orientations are compared
  after composing each one with the docking rotation (and the mirror), modulo
  the point group.
- **Floors come from a profile, by the floor rule**, never from a first run.
  `test=profile` runs the gate at seed offsets 0 to 19 at its CTest thread
  count (and at one thread where it is deterministic) and gives every
  metric's distribution. Let `s` be the robust spread: 1.4826 times the median
  absolute deviation, or the observed range over 4 if that is zero. For a
  metric where higher is better, a floor is admissible when it lies
  - at most at the worst of the 20 runs minus `s`;
  - at most at the median minus `3s`;
  - at least at the negative control's value plus `3s`.

  A ceiling is the mirror image. Every floor is written next to its check
  with its origin: the profile, its date, its median, spread and worst run,
  and the control. A floor moves only by this rule, with the new profile
  written beside it; it is never loosened to make a run pass. If no value is
  admissible, the metric is reported as a robustness profile (sense `report`)
  and the card says why. A contract the code must always keep (a value inside
  the interval it was asked for, a file written, a count) is a check, not a
  floor. There are no blessed baselines and no platform keys.
- **An outcome a single run may miss is a success count, not a gate** (ruling
  of 2026-10-18). Where the profile shows that a correct algorithm misses an
  outcome in some runs by the nature of the method (the de novo multi-state
  routes, which take the best-scoring state outright from random labels, so
  that a state which fits the noise slightly better grows), the outcome is a
  `success` row and `test=profile` holds the number of the 20 runs that
  succeed to a floor: the observed count less three binomial standard
  deviations (the square root of n·p·(1−p), with p the observed fraction),
  rounded down. The floor must be above the count of the shuffled-label
  control. The nightly run of such a route keeps only the gates that held in
  every run of its profile (path assertions, pose ceilings, finite values) and
  reports the rest.
- **Every gate has a negative control** that shows it can fail for the right
  reason: the same metric against a broken truth (shuffled labels, another
  molecule, another micrograph's defocus, copies of one map, shuffled poses)
  must fail the floor. The control's value is a row of the table.
- **Every workflow writes its metrics** to `metrics.tsv` in its run
  directory (the working directory it was started in, never a fixture
  directory it deletes), one `name value floor pass sense` line per check and
  metric, for the nightly report and the profiles. The sense says how the
  floor bounds the value: `min` (the value must reach the floor), `max` (the
  floor is a ceiling), `check` (1 or 0 against 1), `report` (no bound) or
  `success` (1 or 0, never failing the run; its floor is the smallest fraction
  of a profile's runs that must succeed).
  Values are written in exponent form, so a tolerance of 1e-6 keeps its
  digits. Every value a test decides on is a row, written before the decision,
  so a failing run still reports it. A name written more than once in a run (a
  value reported and then held to a floor, a check repeated) is one metric of
  that run, which passes only if every row of it passes.
- **It fails on a missed floor**, with a message naming the metric, its value
  and the floor. Its draws are fixed (`SIMPLE_SEED` and the seed offset), so a
  one-thread run is reproducible; a run on several threads is still a draw
  from a distribution, which is why its floors come from a profile.

### 6.3 The nightly GitHub Actions runs

Nightly integration runs on GitHub-hosted Linux and macOS runners through two
scheduled workflows. Both workflows can also be started manually with
`workflow_dispatch`.

The build workflow, `.github/workflows/ci_build.yml`, checks that the supported
build configurations remain viable:

1. It checks out the selected commit and installs the platform dependencies.
2. It builds SIMPLE in Debug mode with `./compile_debug.sh`.
3. It builds the NICE SIMPLE GUI with `./compile_gui.sh`.
4. In a separate matrix job, it builds the coarray-enabled configuration with
   `./compile_coarrays.sh` where the platform dependencies are available.
5. Each build script includes the test code and runs the `fast` CTest gate by
   default before installation. A compilation or fast-gate failure makes that
   matrix job and the build workflow fail.

The test workflow, `.github/workflows/ci_test.yml`, owns the broader CTest
execution:

1. It checks out the same branch on GitHub-hosted Linux and macOS runners and
   installs the platform dependencies.
2. It performs a clean release build with `./compile_clean.sh`. This build runs
   the `fast` gate before installation; a failure stops the job, so none of the
   longer labels run against a failed build.
3. From `build/`, it runs the capability-dependent platform entries with
   `ctest -L platform --no-tests=error --output-on-failure`.
4. It runs the longer in-process library suites with
   `ctest -L library --no-tests=error --output-on-failure`.
5. It explicitly runs the end-to-end workflow gates with
   `ctest -L highlevel --no-tests=error --output-on-failure`.

The build and test workflows are separate scheduled workflows and may execute
at the same time. The high-level gate is ordered after—and depends on—the clean
build and fast gate inside `ci_test.yml`; it does not wait for the independent
Debug, NICE GUI, or coarray jobs in `ci_build.yml`.

GitHub Actions retains the command output, CTest failure output, duration, job
status, commit, operating system, and dependency-installation log for each
matrix entry. A nonzero build-script or CTest status fails the corresponding
job and therefore the workflow. Tests must not use `continue-on-error` or shell
constructs that hide those exit codes.

The repository README reports the nightly state through two workflow badges:

- **Build SIMPLE** links to `ci_build.yml` and summarizes the Debug, NICE GUI,
  and coarray build workflow.
- **Test SIMPLE** links to `ci_test.yml` and summarizes the clean build, fast
  gate, platform, library, and high-level test workflow.

The badges are workflow-level summaries for the `master` branch. Detailed
results remain on the linked GitHub Actions pages; a badge does not replace the
per-job Linux/macOS results or the metric output produced by a high-level gate.
Every CTest entry retains its declared `TIMEOUT`, and changes that weaken a
timeout or scientific threshold require review.

### 6.4 Benchmark-dataset workflow and webpage

The repository-level build and CTest workflows are complemented by the
`SIMPLE_data_testing` workflow linked under **Benchmark Datasets** in the
README. This is a separate real-data integration layer: it exercises installed
SIMPLE workflows on representative cryo-EM datasets and publishes the results
through the linked GitHub Pages report.

These dataset runs are intentionally short enough for routine automation. Most
datasets are trimmed to 50 micrographs, with the exact subset and order fixed
in the dataset manifest. The report must identify any different limit. Trimming
keeps download, preprocessing, particle extraction, classification, and
reconstruction costs bounded while retaining real detector data and an
end-to-end workflow. It does not make the subset representative of full-scale
production throughput or final biological resolution.

For each dataset, the workflow should record and publish at least:

- dataset and immutable subset identity, including the number of micrographs;
- SIMPLE commit or release, workflow parameters, and execution environment;
- completion status and the stage that failed, if any;
- imported and accepted micrograph and particle counts;
- final-map dimensions, sampling distance, symmetry, and particle count;
- the final FSC curve, the declared FSC criterion, and the resolution obtained
  at that criterion;
- the dataset-specific resolution limit and an explicit FSC **PASS** or
  **FAIL** result;
- important runtime and resource measurements needed to recognize a major
  performance regression.

The FSC decision must be quantitative. The dataset manifest declares the FSC
criterion, commonly 0.143, and the maximum acceptable resolution in Angstrom.
The result passes only when the FSC data are finite and valid, the threshold
crossing can be determined, and the measured resolution is no worse than the
declared dataset-specific limit. A missing curve, invalid crossing, or missed
resolution limit is a failure. Changing that limit requires the same review
and justification as changing a high-level test floor.

Dataset status must not be report-only. A failed workflow stage, invalid
scientific output, or FSC failure must produce a nonzero status that fails the
`SIMPLE_data_testing` GitHub Actions workflow. Its README badge therefore
summarizes the actual benchmark-dataset gate, while the published webpage
provides the per-dataset FSC plot, resolution, threshold, PASS/FAIL decision,
and supporting metrics.

These real-data runs do not replace the hermetic fast, library, platform, or
high-level suites. CTest supplies reproducible generated-data gates within the
SIMPLE repository; the benchmark-dataset workflow supplies short real-data
evidence in its separately versioned data and reporting repository.

### 6.5 The test cards

One card per high-level entry, from the audit of the runs `test_truth` and `test_truth_sync`
(`doc/implementation_notes/planned/test_suite_truth_and_robustness.md`, Progress: Phase 1 for the
first profiles, Phase 5 for the repairs, Phase 7 for the entries whose outputs stayed byte-identical
through the merge with master, Phase 8 for the profiles of the merged code). Every floor named here
is written in the test next to its check with the profile it comes from (section 6.2). "Profile"
means 20 seed offsets at the entry's CTest thread count, and at one thread too where a refactoring
run checks the entry at one thread; the profiles are kept under the runs' `scratch/keep/profiles/`
on the test machine. A row of the nightly table whose sense is `report` is a robustness
measurement, not a gate; a `success` row is held over the seeds of a profile, not in one run.
Cost is the wall time in the `test=nightly` of 2026-10-19 on 24 cores (final code), with the
threads and the CPU reservation of the registry. No floor has been profiled on the maintainer's
Mac yet; the Mac's values of 2026-10-17 are recorded in the plan's Phase 8 entry.

**`mini_stream_6vxx`, `mini_stream_1jxy`**
- Claim: six micrographs simulated from the embedded model (8 particles each) go through CTF
  estimation, picking, extraction, `solve2D` and shape ranking, and the defocus and the picks agree
  with the simulation.
- Checks: every micrograph's defocus against its truth (ceiling 0.13 µm); picked over placed
  particles (floors 0.05 for 6vxx and 0.25 for 1jxy, ceiling 1.5); at least one class; the
  project's micrographs, particle counts, class assignments and populations, shape ranks and class
  averages consistent (checks).
- Path: the `solve2D` record and the `refine2D` iterations it ran (the strategies' benchmark
  records). Control: the fit against the truth of the micrograph whose defocus lies farthest must
  fail the ceiling.
- Profile (2026-10-18, still valid: byte-identical on the merged code): defocus error median
  0.002-0.003 µm, worst 0.125 (a rare wrong fit); pick ratio median 0.33 (6vxx, worst 0.19) and 0.49
  (1jxy, worst 0.40). The simulator saves no particle positions, so the picks are judged by count
  only and have no control.
- Cost: 6vxx 59 s at 8 threads (reserving 2 cores); 1jxy 66 s at 8 threads (reserving 3 cores); each commander one thread.

**`simulated_workflow_6vxx`, `simulated_workflow_1jxy`**
- Claim: simulated particles (6vxx: 160 in 16 movies; 1jxy: 100 in 10 movies) go through import,
  motion correction, CTF estimation, picking, extraction, `solve2D` and `solve3D`, and the map
  matches the map of the atomic model; 1jxy runs `solve3D` a second time distributed over two worker
  processes.
- Checks: at least 4 extracted particles and one class; the map written with the expected box and
  sampling. Reported (no admissible floor): the map's correlation with its truth and its FSC.
- Path: the distributed run wrote its second part's record. Control: the map docked onto the other
  suite's molecule must stay under the correlation every found map reached (0.80).
- Profile (2026-10-19, one thread): `solve3D` missed the map, at the control's correlation (median
  0.53-0.55), in 2 of 20 runs for 6vxx and in 1 of 20 for 1jxy in each of its two steps; the found
  maps correlate at 0.88 (s 0.02-0.03) with FSC 4.5-4.8 Å. A robustness finding left to the
  maintainer, with a proposal (hold "the map is found" as a success count).
- Identity: the declared outputs of a one-thread run (the entry is deterministic).
- Cost: 6vxx 1129 s at 1 thread (reserving 2 cores); 1jxy 999 s at 1 thread (reserving 3 cores) with two one-thread workers in its distributed step.

**`single_workflow_fcc`, `single_workflow_wurtzite`**
- Claim: a simulated Pt (face-centred cubic) or CdSe (wurtzite) nanoparticle goes through
  reprojection, a trajectory, denoising, `analysis2D_nano`, `autorefine3D_nano`, `refine3D_nano
  refine=cont` and the polish, and its map and poses match the truth.
- Checks: the map's box and sampling, its correlation with the simulated map (floor 0.90) and FSC
  (ceiling 5 Å); the median rotation error of the particles in the truth's frame, after one
  docking of the map to 2 Å and up to the lattice's group (O or D6) and the in-plane rotation by
  180° (both lattices are nearly centrosymmetric), reported for `autorefine3D_nano` and held for
  `refine=cont` and the polish to at most 0.5° worse than the stage before.
- Path: the stages' outputs. Identity: the declared outputs of a one-thread run (deterministic).
- Profile: fcc 2026-10-18 (byte-identical on the merged code): rotation error median 4.9°,
  stage-to-stage excess worst +0.19°. wurtzite 2026-10-19 (the merged code no longer clamps the
  matching low-pass): rotation error median 1.5°, excess worst −0.12°, map correlation worst 0.905,
  FSC worst 3.0 Å. Random poses err by tens of degrees.
- Cost: fcc 592 s at 1 thread (reserving 2 cores); wurtzite 586 s at 1 thread (reserving 2 cores).

**`simulate_particles`**
- Claim: `reproject` and `simulate_particles` make the images they are asked for.
- Checks: the images against the reprojections at their recorded poses (50 images without CTF:
  each closer to its own pose's reprojection than to another image's); with three blob volumes and
  given populations, exact state counts and every image closest to its own state's projection;
  every CTF and shift draw inside its requested interval (checks of the simulator's contract; the
  largest of many uniform draws sits just inside its bound by construction, so no floor applies,
  and the maxima are reported).
- Profile: passes 20 of 20 (2026-10-18, byte-identical on the merged code).
- Cost: 11 s at 8 threads (reserving 8 cores).

**`state_fixture`**
- Claim: the multi-state test fixture builds the datasets the multi-state tests rely on.
- Checks: the blob phantoms in the search, known-pose and exact-view modes (built twice, the second
  from the cache, the records equal) and the beta-galactosidase states, each truth holding its
  declared point group (correlation at least 0.99 over the group) and no larger one (at least 0.01
  above the next larger group).
- Profile: passes 20 of 20 (2026-10-18, its gate table identical on the merged code).
- Cost: 7 s at 8 threads (reserving 8 cores) (from the cache; about 25 s when it builds the datasets).

**`solve3D_addon`**
- Claim: `solve3D_addon` extends a frozen solution with new particles without disturbing it.
- Checks: the union map against the truth (correlation floor 0.90, at most 0.03 below the base map,
  FSC at most 1 Å worse than the base); coverage of the new cohort (floor 0.90); the cohort's pose
  error in excess of the frozen particles' own (ceiling 10°); the frozen project, map and sigma2
  unchanged, the manifests and bundles valid (checks).
- Control: the cohort's truth shuffled must exceed the ceiling.
- Profile (2026-10-19, eight threads and one thread, the speed week's one-thread case): excess
  median 0.49° / 0.11°, spread 0.80 / 1.10, worst 1.80° / 2.16°; the shuffled control at 28.6° /
  28.5° (worst 24.7°); the union map worst 0.938 / 0.937.
- Cost: 289 s at 8 threads (reserving 11 cores); about 25 minutes at one thread.

**`cont_refine3D_1jxy`**
- Claim: continuous Cartesian refinement recovers the poses of 1 000 simulated 1JYX particles from
  perturbed starts, in every route that uses it.
- Checks: median rotation error at most 1.5°, shift at most 0.2 px, at least 95 % inside the basin
  (a tenth of the injected error), on the shift-then-joint seam (stages a-c) and on the production
  route (stages i, j); the distributed run equal to the shared-memory one within the rounding of
  threaded evaluation; the polish and the continuation no worse than without them.
- Path: no polar in-plane route where none is expected; the distributed stage ran two parts.
- Control: every E25 stage's rotation error and basin fraction against the truth shuffled among
  the particles must fail the floors (the simulation has no shift, so the shift error has none).
- Profile (2026-10-19, eight threads and one thread): passes 40 of 40, every floor far from its
  bound (stage a rotation worst 1.00°, stage e worst 2.22° against 7.5°). Stages i, j and the
  controls (2026-10-20, eight threads): passes 20 of 20, stage i rotation worst 1.01°, stage j
  0.24°, basin 1.000 in every run; the shuffled controls at 128.8-135.2° and basin at most 0.001.
  With the continuous gather (same day): 20 of 20, every floor admissible (stage a worst 1.02°).
- Cost: a median of 390 s at 8 threads with stages i and j (the 2026-10-20 profile; reserving 10
  cores), up to 11 GiB; about 37 minutes at one thread (2026-10-19, before stages i and j).

**`multistate_assignment`**
- Claim: `refine3D nstates=3` (`refine=prob` and `prob_neigh`) recovers the known states of three
  blob phantoms (1 500 particles, populations 50/30/20 %, SNR 2) from the truths low-passed to 20 Å
  and random initial states, scored by correlation.
- Checks: adjusted Rand index (floors 0.89 / 0.80), accuracy (0.95 / 0.92), smallest per-state
  recall (0.93 / 0.85), largest population error (ceilings 0.03 / 0.04), no estimated state empty,
  the map mapping agreeing with the label mapping, the maps' discrimination margin (0.30), the
  poses and shifts of the correctly assigned particles (6.0°, 0.40 px).
- Controls: shuffled labels inside the accuracy's null; copies of one map and of the mean map not
  distinct.
- Profile (2026-10-19, eight threads and one thread): every state resolved in all 80 runs; the
  `prob` adjusted Rand index at one thread (median 0.941, s 0.015) moved its floor from 0.90 to 0.89.
- Cost: 78 s at 8 threads (reserving 11 cores); about 10 minutes at one thread.

**`multistate_workflows`**
- Claim: the multi-state workflows recover three compositional beta-galactosidase states (the
  tetramer, the tetramer less chain D, the A-B dimer; 50/30/20 %) in four routes: `cavgs`
  (`solve3D_cavgs` on 600 simulated class averages imported with `import_cavgs`), `shc` (`solve3D`
  from 3 000 particles), `states` (`refine3D_states` from the true poses and the consensus map) and
  `consensus` (single-state `solve3D` in C1, then `refine3D_states`: the route used on real data).
- Checks: the path, the controls and the pose ceilings of every route (class averages 4.5° /
  0.40 px, particles 2.5° / 0.15 px, `refine3D_states` 2.0° / 0.10 px, consensus 2.5° / 0.15 px);
  for `states` and `consensus` every true state resolved (recall above one half) with no estimated
  state empty, and their labels held to floors (`states`: adjusted Rand index, accuracy and recall
  0.99, population error 0.01, discrimination 0.40; `consensus`: adjusted Rand index 0.40, accuracy
  0.70, population error 0.30, discrimination 0.20, its recall reported).
- Successes (ruling of 2026-10-18): the de novo routes `shc` and `cavgs` start from random labels and
  take the best-scoring state outright, so whether every true state is resolved is decided early and
  a run may miss it; it is a `success` row, held by `test=profile` to 6 of 20 seeds (`shc`) and 2 of
  20 (`cavgs`), the shuffled-label control at 0; their labels are reported.
- Path: the refine mode, directions and states of every iteration from the refine3D strategy's
  record. Controls: shuffled labels, copies of one map.
- Profile (2026-10-19, eight threads; `consensus` at one thread too): `shc` resolves every state in
  13 of 20 runs, `cavgs` in 9, `states` and `consensus` in all (`consensus` 40 of 40, adjusted Rand
  index median 0.78, one state often only partly separated). The per-seed metrics of the de novo
  routes are the baseline of the planned annealed state draw.
- By hand: `suite=cavgs_single`, `suite=cavgs_as_particles`, `suite=<route>_identical`,
  `suite=hard_<route>`.
- Cost: 1809 s at 8 threads (reserving 11 cores) for the four routes; the consensus route alone about 70 minutes at one thread.

**`single_atoms_stats` and its CdSe variants (`_cdsew`, `_cdsez`, `_cdser`)**
- Claim: atom detection and statistics recover the atoms of simulated Pt and CdSe (wurtzite,
  zincblende, rocksalt) nanoparticles.
- Checks: recall and precision at least 0.90, position error at most 0.75 Å, density correlation
  at least 0.80, the coordination number of interior atoms and the fitted lattice constant within
  2 % (CdSe). Control: the atoms shifted by half a bond must fail.
- Profile: passes 20 of 20 (noise-free: the values repeat; byte-identical on the merged code).
- Cost: 19, 10, 12 and 14 s at 8 threads (reserving 3 or 4 cores).

**`species_discovery`**
- Claim: `detect_atoms` with a species list and with `discover_species=yes` recovers the elements of
  simulated Pt, Pt/Ni and Pt/Al particles.
- Checks: Pt recall 0.98, all species 0.95, the second element 0.89 (aluminium 0.88), at most one
  false atom, half-map label agreement 0.95. Control: `element=Pt` alone finds less than 0.90 of
  the aluminium.
- Profile (2026-10-18, its gate table identical on the merged code): the second element's recall
  median 0.977, worst 0.908.
- Cost: 46 s at 8 threads (reserving 8 cores).

**`stream_preproc`**
- Claim: five simulated movies go through the stream's preprocessing stage, its work done by worker
  jobs of the local queue, and the defocus, frame drift and image are recovered.
- Checks: defocus within 0.20 µm, drift within 0.5 px, the integrated micrograph's correlation with
  the optimal average at least 0.50 and 0.05 above another movie's.
- Path: the worker jobs finished, each with exit status 0. Control: the farthest other movie's
  defocus must fail the ceiling.
- Profile (byte-identical on the merged code): defocus worst 0.134 µm.
- Cost: 82 s at 8 threads (reserving 2 cores), its jobs one thread each.

**`flex_pca_blobs`**
- Claim: `flex_pca` separates two blob states at their true poses, in shared memory and with two
  worker processes; the per-state refinement holds; the class-average state route runs end to end.
- Checks: label accuracy at least 0.60 (one-to-one mapping) in both modes, the two maps matching
  the two truths, the per-state refined maps' correlation; the route's contract (projects, hand-back
  of the states, one map per state). Reported: the class-average route's state accuracy, the number
  of planted states its maps match, and its poses, on three beta-galactosidase states in exact views.
- Path: the distributed run's parts ran as worker jobs. Control: shuffled labels fail the floor; the
  route's accuracy beside its shuffled-label null.
- Profile (2026-10-19): label accuracy median 0.989, worst 0.907; the class-average route at the
  chance level (four runs deliver one state, one run's maps match a single planted state), a
  finding left to the maintainer.
- Cost: 404 s at 8 threads (reserving 9 cores), up to 8.7 GiB.

### 6.6 Real data: the dataset layer (planned)

Real datasets carry no truth of their own, so a real-data case never judges
poses or state labels. By the rulings of 2026-10-18 (recorded in
`doc/implementation_notes/planned/testing_and_benchmarking_plan.md`, section 5),
real data joins `test=nightly` as public datasets with previously solved maps:

- **the data**: subsets of public datasets (EMPIAR images whose maps are in
  EMDB, under a licence that allows a subset to be redistributed), kept as
  versioned Zenodo records with the subset, the reference map, the case's
  command lines and a checksum per file, fetched by name into a local cache
  and checked against the checksums; a missing or corrupt dataset skips its
  case with a visible status, never a pass;
- **the judgement**: the run's map docked in both hands onto the earlier map,
  their FSC and correlation in a band that stops short of the reference's own
  resolution, with floors by the floor rule from a 20-seed profile and another
  dataset's reference as the negative control. The reference is never an
  input of the run (that would measure reference bias), and it judges the
  map only: it came from the same images. A multi-state dataset is compared
  state by state only where per-state maps were deposited, and its labels are
  never judged;
- **robustness and performance as before**: does the workflow complete, does
  it hold its resolution limit (with cFAR, the conical FSC area ratio, beside
  the FSC=0.143 resolution; FSC=0.5 is never a gate on a dataset whose FSC runs
  along 0.5), how long does it take, how much do the results scatter across
  seeds;
- **the components**: the fetch and cache with checksums, a SIMPLE test program
  that runs each case through the same dispatcher as `test=nightly` and
  `test=profile`, and the same `metrics.tsv` rows (map-to-reference FSC and
  correlation as floors, resolution and cFAR as ceilings and floors, wall time
  and per-stage `BENCH` times as reports).

The full-size benchmark datasets that are not public (PfCRT, msp1, the
beta-galactosidase speed cases) stay as local copies behind a per-machine
manifest, for profiles and timings. **Open:** the design of each case (its
starting point, downsampling, number of images and budget, which fixes what it
can test in the time it is given), the choice of datasets and their licences,
how the fetch is done (a SIMPLE program or a one-time setup step per machine),
and whether `SIMPLE_data_testing`'s trimmed datasets and the Zenodo records are
the same data. No dataset has been chosen, downloaded or uploaded yet. Until the
layer exists, real-data runs stay by hand or in the `SIMPLE_data_testing`
workflow (section 6.4).

## 7. Why tests are no longer standalone programs

Until September 2026 a test could be a program unit in `production/tests`,
built by a glob into its own executable, installed, and called by name from
CI. That route is closed; `production/tests` and the glob are gone. The
reasons, all found in the review:

- **Two routes, drifting apart.** 83 standalone programs and 83
  `simple_test_exec` cases covered 130 identities; 36 existed on both routes,
  and the twins had diverged (`cmdline` duplicated the `command line` unit
  test; `imgfile` was a subset of `test_image`).
- **No failure path.** 63 identities could not fail. A program that runs to
  its end reports success whatever it computed, and several printed
  "PASSED" with no way to fail.
- **Process control in the test.** A `stop` ends the process: one failed
  check hid every later one, and `stop` with a message exits with status 0
  (`openmp_offload` had 65 of them).
- **No shared lifecycle.** Hand-parsed arguments instead of the UI and
  `parameters`, ad hoc working directories, files left behind, unseeded
  random numbers.
- **One process per test.** Nothing could be budgeted or run as a gate; each
  program needed its own executable, install rule and CI call, and many were
  never run at all.
- **Dead code kept alive.** Production routines survived because a program
  in the test area called them; the review deleted dozens of routines and
  modules that nothing else used.

A test procedure in a tester module has none of these problems: it shares
the runner, the assertions, the seeding and the report with every other test,
and a suite runs many of them in one process. Program units are for
production executables.

## 8. What the refactoring deleted, and why

The review gave every test identity in the inventory a verdict. The
retired-tests table of `doc/refactoring_notes/completed/test_review_record.md` has one
row per removed identity (136), with its reason and replacement. Most were
not lost: 79 were merged into tester modules with real assertions, 10 were
modified or moved, and 47 were deleted or retired outright. The reasons for
deletion fall into these groups.

| reason | tests |
|---|---|
| printed or plotted results for eyeballing; nothing to assert that other tests do not assert | `uniform_euler`, `uniform_rot`, `order_corr`, `phasecorr`, `ptcl_center` |
| benchmarks and timing loops | `rotate_ref`, `eval_polarftcc`, `io`, `io_parallel`, `star_export`, `openacc`, `openmp`, `simd` |
| duplicates of a route or of a gate sub-suite | `starfile`, `binoris_test`, `binoris_io_test`, `imgfile`, the standalone `mini_stream`, `gui_assembler`, `gui_metadata`, `project_merge`, `clustering`, `multinomal_test` |
| runners on user data or downloads that asserted nothing | `continuous_inplane_rotation2D`, `continuous_inplane_rotation2D_metadata`, `nu_filter`, `create_gain`, `search_gain_flips`, `atomfit`, `eo_diff`, `opt_lp`, `cif2mrc`, `cif2pdb` |
| not needed: a production program runs the same code, or replaces the test | `nu_envmask` (`nu_filt3D`), `phase_rand_fsc` (`fsc`), `angres` (now the program `measure_projspace_angres`) |
| drove code that had no production caller, deleted with it | `subproject_distr`, `ptcls_ppca_subproject_distr`, `socket_client`, `socket_comm_distr`, `socket_io`, `socket_server` |
| broken or empty | `install` (ran an executable that no longer existed), `nice` (posted to a server that does not exist), `stream_initial_analysis` (a commander smoke on a missing folder) |
| tested a hand-written command line instead of the production one | `abinitio2D_stream` |

## 9. Where is my test now?

Look your old test up below. "Deleted" means it was removed without a
replacement; the reason is in section 8 and in the inventory. A sub-suite is
written as `sub-suite` (entry).

| old test | now |
|---|---|
| `abinitio2D_stream` | deleted; `solve2D` runs nightly in `simulated_workflow_6vxx` and `simulated_workflow_1jxy` |
| `angres` | the program `simple_exec prg=measure_projspace_angres nspace=<n>`; the table of values is a comment above `find_angres` |
| `ansi_colors` | `string` (`unit_core`) |
| `assign_optics` | `optics assignment` (`lib_stream`) |
| `atomfit` | deleted, with `atoms%fit_bfactors`, which nothing called |
| `atoms_stats` | the high-level entries `single_atoms_stats` (Pt) and `single_atoms_stats_cdsew`, `_cdsez`, `_cdser` (`simple_test_exec test=single_atoms_stats`); the `nanoparticle atoms` sub-suite of `lib_single` no longer exists |
| `binoris`, `binoris_io`, `inside_write` | `binoris` (`unit_project`) |
| `binoris_test`, `binoris_io_test` | deleted (empty stubs); see `binoris` |
| `bounds_from_mask3D`, `bounds_from_mask3D_test`, `graphene_mask`, `mask`, `msk_routines` | `masks` (`unit_image`) |
| `cartesian_fourier` | `Cartesian Fourier` (`unit_cart_align3D`) |
| `cavg_quality_relations` | `cavg quality relations` (`unit_numerics`) |
| `cavg_registration` | `cavg registration` (`unit_pftc_align2D3D`) |
| `cc_connectivity`, `image_bin` | `binary image` (`unit_image`) |
| `cif2mrc`, `cif2pdb` | deleted; they ran the production programs of the same names |
| `class_sample`, `class_sample_test` | `class sample I/O` (`unit_core`) |
| `clustering` | deleted; it called `affinity propagation` (`unit_numerics`) |
| `cmdline` | `command line` (`unit_core`) |
| `coarrays` | the capability-gated platform suite `coarrays`; the coarray-linked private executable runs two images through initial and final synchronization barriers |
| `continuous_3D_pcg_reconstruction` | `observation noise` (`unit_reconstruction`) and `PCG half-set` (`lib_reconstruction`) |
| `continuous_inplane_cc_grad`, `continuous_inplane_hybrid_grad`, `continuous_inplane_rotation2D_stage1_validation`, `continuous_inplane_rotation2D_route_identity` | `continuous in-plane` (`unit_pftc_align2D3D`) |
| `continuous_inplane_refine3D` | `refine3D in-plane state` and `continuous in-plane` (`unit_pftc_align2D3D`) |
| `continuous_inplane_rotation2D`, `continuous_inplane_rotation2D_metadata` | deleted (a shell driver; a post-run scan of a user project) |
| `corrs2weights`, `corrs2weights_test`, `rank_weights` | `statistics` (`unit_numerics`) |
| `create_gain`, `search_gain_flips` | `motion gain` (`unit_motion`) |
| `ctf`, `ctf_test` | `CTF` (`unit_image`) |
| `detect_atoms`, `simulate_nanoparticle` | no test of their own; they run as stages of the high-level `single_atoms_stats` entries, `detect_atoms` also in `species_discovery`, and `simulate_nanoparticle` also in `single_workflow_fcc` and `single_workflow_wurtzite` |
| `detect_calpha` | `C-alpha finder` (`unit_single`) |
| `detect_calpha_molecules` | deleted; quantitative synthetic coverage remains in `C-alpha finder` (`unit_single`) |
| `diff_map_graphs` | `diffusion-map graphs` (`unit_numerics`) |
| `discrete_stack_io`, `stack_io` | `stack I/O` (`unit_core`) |
| `eigh`, `eigh_test` | `linear algebra` (`unit_numerics`) |
| `eo_diff`, `opt_lp` | deleted (needed refine3D volumes or a download; asserted nothing) |
| `eul_prob_tab2D_io` | `2D probability table I/O` (`unit_pftc_align2D3D`) |
| `extr_frac` | `decay schedules` (`unit_numerics`) |
| `flex_gpu` | deleted with the flex GPU implementation |
| `flex_pca` | `flex PCA` (`unit_heterogeneity`) |
| `flex_pcg` | `flex PCG operator` (`unit_heterogeneity`), `flex PCG operator 64` and `flex PCG solve sweep` (`lib_heterogeneity`) |
| `forked_process` | the platform entry `forked_process` |
| `ft_expanded` | `shift search` (`unit_numerics`) |
| `gencorrs_fft` | `polar correlation` (`unit_pftc_align2D3D`) |
| `gen_pickrefs` | `picking references` (`lib_stream`) |
| `gui_assembler`, `gui_metadata` | `GUI assembler`, `GUI metadata` (`unit_ui`) |
| `imgfile` | `image` (`unit_image`), the SPIDER and MRC round trips |
| `install` | deleted; the fast gate replaces it |
| `io`, `io_parallel`, `star_export` | deleted (benchmarks); the round trips are in `stack I/O` (`unit_core`) and `STAR file` (`unit_project`) |
| `kbinterpol_fast` | `Kaiser-Bessel kernel` (`unit_numerics`) |
| `lbfgsb`, `lbfgsb_cosine` | `optimisers` (`unit_numerics`) |
| `lplims`, `lpstages`, `lpstages_test` | `low-pass stages` (`unit_numerics`) |
| `master` | `stream heartbeat` (`forked_process`) |
| `maxnloc`, `maxnloc_test` | `search, sort, locate` (`unit_numerics`) |
| `mini_stream` | the high-level entries `mini_stream_6vxx` and `mini_stream_1jxy` (simulated micrographs) |
| `mrc2jpeg`, `mrc_validate` | `mrc2jpeg`, `mrc validate` (`unit_image`) |
| `multinomal_test`, `rnd_shuffle` | `random draws` (`unit_numerics`) |
| `nano_mask` | `nano mask` (`unit_image`) |
| `score_volume_shape` | `volume shape` (`unit_image`) |
| `neigh`, `sym`, `sym_test` | `symmetry` (`unit_ori`) |
| `nice` | deleted |
| `nu_envmask`, `nu_filter` | deleted; `simple_exec prg=nu_filt3D` |
| `openacc`, `openmp`, `simd` | deleted (demonstrations and timings) |
| `openmp_offload` | the platform entry `openmp_offload` |
| `order_corr` | deleted; its one check is asserted in `orientation collection` (`unit_ori`) |
| `ori`, `ori_test` | `orientation` (`unit_ori`) |
| `oris`, `oris_test` | `orientation collection` (`unit_ori`) |
| `otsu`, `otsu_test`, `peak_thres_fdr` | `segmentation` (`unit_image`) |
| `pca_all`, `pca_imgvar` | `PCA` (`unit_numerics`) |
| `pcg_frac_update` | the library entry `pcg_frac_update` |
| `rec3D_backends` | by hand, `simple_test_exec test=rec3D_backends` |
| `pcg_recon` | the library entry `pcg_recon` |
| `pdb2mrc` | `pdb2mrc` (`lib_single`) |
| `phase_rand_fsc` | deleted; `simple_exec prg=fsc` |
| `phasecorr`, `ptcl_center`, `rotate_ref`, `uniform_euler`, `uniform_rot`, `eval_polarftcc` | deleted; the image basics `ptcl_center` touched are in `image` (`unit_image`), `rotate_ref_8` is in `polar correlation` (`unit_pftc_align2D3D`) |
| `phshift_policy`, `ui_visibility` | `UI visibility` (`unit_ui`) |
| `phshift_star` | `STAR project` (`unit_project`) |
| `pick_extract` | `pick and extract` (`lib_stream`) |
| `pose_cont_refine3D_adapter` | the adapter is gone (pose_cont refactoring): `cart calculator`, `pose optimizer`, `pose strategy`, `pose statistics` (`unit_cart_align3D`), `pose 1JYX recovery` (`lib_cart_align3D`) and the workflow entry `cont_refine3D_1jxy`; the workflow sub-suite went with the removed program |
| `pose_cont_refinement` | `cart calculator` and `pose optimizer` (`unit_cart_align3D`); the pose refiner was split into them |
| `preproc` | the high-level workflow entry `stream_preproc` |
| `project_merge` | `project merge` (`unit_project`) |
| `projdir_accumulator` | retired with `projrec` (release 4) |
| `qsys_ctrl`, `qsys_env` | `qsys control`, `qsys environment` (`unit_parallel`) |
| `rec3D_backend` | `rec3D backend` (`unit_reconstruction`) |
| `reproject` | the workflow entry `simulate_particles` |
| `serialize` | `image serialisation` (`unit_image`) |
| `sieve_cavgs` | `particle sieve` (`unit_project`) |
| `sigma2_state` | `sigma2 state` (`unit_pftc_align2D3D`) |
| `simulate_particles` | the workflow entry `simulate_particles` |
| `simulated_workflow` | the workflow entries `simulated_workflow_6vxx`, `simulated_workflow_1jxy` |
| `single_workflow` | the workflow entries `single_workflow_fcc`, `single_workflow_wurtzite` |
| `socket_client`, `socket_comm_distr`, `socket_io`, `socket_server` | deleted with the socket modules; `IPC TCP socket` (`unit_ipc`) tests the live transport |
| `sp_project` | `project records` (`unit_project`) |
| `starfile`, `starfile_test` | `STAR file` (`unit_project`) |
| `stream_initial_analysis` | deleted |
| `stringmatch` | `string` (`unit_core`) |
| `subproject_distr`, `ptcls_ppca_subproject_distr` | deleted with the subproject scheduling code, which had no other caller |
| `trail_rec_blend` | `trailing-reconstruction blend` (`unit_image`) |
| `ui_hash_test` | `UI hash` (`unit_ui`) |
| `units` | `simple_test_exec test=units` (all fifteen area suites, not a CTest entry) |
| the `unit_<area>` suites | unchanged names, CTest label `fast` |
