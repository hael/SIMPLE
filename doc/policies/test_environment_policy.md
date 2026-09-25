# Test Environment Policy

## Purpose and scope

This policy says what a test is in SIMPLE, where each test goes, how to write
one, and what runs when. It applies to every Fortran test in the repository
and to everyone who adds, changes or removes one. It records the state after
the test-environment refactoring of September 2026.

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
| `fast` | 13 area suites `unit_<area>` | every `compile_*.sh` build (unless `--exclude-tests`), before installation | unit tests of the library: hermetic, in-process, one thread, seconds |
| `library` | 5 library suites `lib_<area>` | nightly | longer numerical tests on generated data: realistic sizes, minutes |
| `highlevel` | 9 high-level gates | explicit CTest command only | long simulated pipelines and commander integrations, including independent 6VXX/1JXY suites |
| `platform` | `forked_process`, plus `coarrays`, `flex_gpu`, `openmp_offload` when CMake finds the capability | by hand, nightly where the machine has the capability, and `coarrays` during `compile_coarrays.sh` | tests that need child processes, a launcher or a device; `coarrays` is a two-image integer-transfer smoke |

The fast tier is the build-time gate. Library and supported platform tests may
run overnight; high-level tests run only when explicitly selected with CTest.

**The fast gate is part of the build.** Every `compile_*.sh` build
runs `scripts/run_fast_gate.sh` between `make` and `make install`. It first
runs `scripts/check_test_registry.py` (section 4.5), then
`ctest -L fast` with half the cores, then `scripts/ctest_budget.py`, which
fails the build when an entry fails or the gate takes more than 30 s of real
time. A failed gate installs nothing. The per-entry timings are kept in
`build/test_runs/ctest_fast.log.timing.txt`; the gate takes about 5 s on the
reference Mac in Debug. Tests are on by default in every compile script;
`--exclude-tests` (`BUILD_TESTS=OFF`) builds the library and executables only
and skips the gate, for when only the executables are needed.

`compile_coarrays.sh` additionally runs the capability-gated
`coarrays` CTest entry after the fast gate and before installation. A failed
two-image smoke test therefore prevents a coarray build from being installed.

**The process budget.** The number of CTest entries is fixed in
`SIMPLE_CTEST_BUDGET` (`production/CMakeLists.txt`, currently 28: 13 fast,
5 library, 9 highlevel, 1 platform) and configuration fails when it does not
match. A CTest entry is an isolation unit, not a place for one more check:
checks are added inside existing suites. A new entry needs a stated reason and
the owner's agreement, and is recorded in the plan.

**Every entry** runs in its own directory `build/test_runs/<entry>`, with its
own `TIMEOUT`, an explicit `OMP_NUM_THREADS` (1 for the fast tier) and
`SIMPLE_SEED=20260923`. Each suite writes a report,
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
| `unit_numerics` | online variance, random draws, affinity propagation, statistics, linear algebra, Kaiser-Bessel kernel, search/sort/locate, decay schedules, PCA, cavg quality relations, diffusion-map graphs, optimisers, low-pass stages, shift search |
| `unit_project` | STAR file, STAR project, binoris, project records, project merge, class compatibility, particle sieve, motion gain, motion model |
| `unit_ui` | UI JSON, GUI metadata, GUI assembler, UI hash, UI visibility |
| `unit_ipc` | IPC TCP socket, HTTP POST, persistent worker message, persistent worker server (localhost only) |
| `unit_reconstruction` | rec3D backend, observation noise, class-average accumulator |
| `unit_pftc_align2D3D` | polar correlation, continuous in-plane, refine3D in-plane state, 2D probability table I/O, sigma2 state, cavg registration |
| `unit_cart_align3D` | Cartesian Fourier, pose refiner, pose adapter |
| `unit_heterogeneity` | flex PCA, flex PCG operator |
| `unit_parallel` | qsys control, qsys environment |
| `unit_single` | atoms, cif2mrc, C-alpha finder |

`simple_test_exec test=units` runs all thirteen in one process. It is a
convenience and deliberately not a CTest entry.

### 1.2 Long-running CTest entries

| entry | label | what it runs |
|---|---|---|
| `lib_reconstruction` | library | PCG half-set: independent half-set PCG solves against gridding |
| `lib_cart_align3D` | library | pose 1JYX recovery: 5 000 simulated 1JYX particles refined by the Cartesian pose refiner |
| `lib_heterogeneity` | library | flex PCA deconvolution of 20 000 particles, the PCG operator at box 64, the PCG solve sweep |
| `lib_single` | library | pdb2mrc coverage of the built-in molecular models |
| `lib_stream` | library | optics assignment, picking references, pick and extract |
| `mini_stream_6vxx`, `mini_stream_1jxy` | highlevel | independent embedded-model mini-stream validations |
| `simulated_workflow_6vxx`, `simulated_workflow_1jxy` | highlevel | simulated movies through import, motion correction, CTF, picking, extraction, `abinitio2D`, `abinitio3D` |
| `single_workflow` | highlevel | the SINGLE pipeline on a simulated Pt nanoparticle |
| `pcg_recon` | highlevel | gated stages of the PCG reconstruction operator |
| `simulate_particles` | highlevel | `reproject` and `simulate_particles` on the embedded 6VXX volume |
| `single_atoms_stats` | highlevel | simulated Pt nanoparticle atom detection and statistics |
| `stream_preproc` | highlevel | five simulated movies through the stream's preprocessing stage and its worker jobs |

### 1.3 Programs that are not registered

A few programs are reachable through `simple_test_exec` but are not CTest
entries because they need data a user supplies: `pcg_frac_update` and
`rec3D_backends`. No new program joins this list
without a reason: a diagnostic that runs on a user's data is a developer
program (section 3.4), not a test.

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
```

The `suite=` selector is the sub-suite name in lower case, with blanks and
hyphens written as underscores and commas and slashes dropped
(`stack I/O` is `stack_io`, `search, sort, locate` is `search_sort_locate`).
The selectors of each suite are listed in the description of its `suite=`
input in the test UI (`src/main/ui/simple_test/simple_test_ui_class.f90`). To
reproduce a CTest run by hand, export `SIMPLE_SEED=20260923` first.

`SIMPLE_UNIT_ORDER=reverse` is a diagnostic for the fused, in-process suites:
it reverses the sub-suite table so a result that differs from the normal run
exposes state leaked by one sub-suite into another. It does not reverse CTest
entries and is not part of the normal build gate; use it after changing shared
state, lifecycle or test-runner code.

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

## 3. Where a test goes

### 3.1 The decision

Answer these in order.

1. **Does it assert something (section 2.1)?** If not, it does not go in the
   test environment at all; see section 3.4.
2. **Which area does it test?** A test goes to the area of the code it
   tests, not the area of its author or of the workflow that uses the code.
   The areas are named after the machinery and kept short: `core`, `ori`,
   `image`, `numerics`, `project`, `ui`, `ipc`, `reconstruction`,
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
4. It is cheap. The whole gate must stay under 30 s of real time; keep every
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
done. The nightly CI builds and gates master on Linux and macOS, in Debug and
Release, and shows the next morning what slipped through. Do not keep work on a branch in the
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
  test is about; name them `tmp_<tester>_*` and delete them at the end of the
  test, also when a check failed.
- A tester that loops over something printing through `logfhandle` diverts
  `logfhandle` to a scratch unit and restores it.

### 4.3 White-box self-tests

A self-test that must read private components of a type may stay inside the
production module. It then takes an out-argument per check, and a tester
asserts those by name (`test_flex_pcg_operator` in `simple_flex_pca_pcg` is
the model). Every other self-test belongs in a tester module: the review
moved the self-tests of `image`, `imghead`, `atoms`, `oris`, `ftiter`,
`ftexp_shsrch`, `bspline_smoother`, `online_var` and `aff_prop`, the flex
PCA, UI hash, cavg-quality and class-average registration self-tests out of
their production modules, and deleted those of unused code (`hclust`,
`srchspace_map2D_io`) and the dead ones (`CPlot2D`, `jpg`). The only self-tests
left in production modules are this white-box one and the five `flex_gpu`
tests of the GPU platform entry.

### 4.4 A workflow or platform test

A test that runs a pipeline is a commander: its type and `exec_test_<name>`
routine in `src/main/commanders/test/simple_commanders_test_<area>.f90`, a
case in the router `src/main/exec/simple_test_exec_<area>.f90`, a program in
`src/main/ui/simple_test/simple_test_ui_<area>.f90`, and a registration in
`production/CMakeLists.txt` with `simple_add_test` (label, timeout, threads,
`RUN_SERIAL`). Section 6 says what it has to check. It still fails through
assertions or `THROW_HARD` with a message that names the metric, its value
and the floor.

### 4.5 Checks and records

Before asking for a build:

- `python3 scripts/check_descr.py .` — every source file starts with a
  one-line `!@descr:`.
- `python3 scripts/check_test_registry.py . --verbose` — the CTest
  registrations, the test UI programs, the router cases and the `suite=`
  lists agree. The fast gate runs it on every build.
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
  suite runs a team of three) and compares with the serial result. The
  nightly run repeats the fast tier on four threads.
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
- An `rmat` pointer (`get_rmat_ptr`) is the padded array, with extra rows in
  the first dimension for the in-place FFT; bound it by `ldim`,
  `rmat(:ldim(1),:ldim(2),:ldim(3))`, before a whole-array expression. The
  unbounded pointer does not conform with a box-sized array: bounds checking
  stops on it, and without bounds checking it can pass by luck (the `cif2mrc`
  tester, 2026-09-25). `get_rmat()` returns a copy of the box.

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
| `set_fixed_seed(base)` | seed the intrinsic generator with a fixed state |
| `begin_test_suite`, `end_test_suite`, `reset_test_report`, `report_summary`, `tests_run`, `tests_failed` | the runner's bookkeeping; not for testers |

A failed assertion records its message, with expected and actual values, and
the suite goes on: every check of a run is reported, not only the first
failure.

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
- **Maps are compared after docking and hand.** An `abinitio3D` map is
  neither docked to the truth map nor necessarily of the right hand. The gate
  docks it with `dock_vols` at a low-pass of 15 to 20 Å, keeps the hand (the
  map or its `mirror('x')`) that docks with the higher correlation, and takes
  the masked FSC at 0.143 against a declared floor. Orientations are compared
  after composing each one with the docking rotation (and the mirror), modulo
  the point group.
- **Floors are declared in the test**, next to where their value comes from:
  the first measured runs with a stated margin, or a physical argument. A
  floor is loosened only with a written justification in the commit. There
  are no blessed baselines and no platform keys.
- **Every workflow writes its metrics** to `metrics.tsv` in its run
  directory, one `name value floor pass` line per metric, for the nightly
  summary.
- **It fails on a missed floor**, with a message naming the metric, its value
  and the floor, and it is reproducible: `SIMPLE_SEED` fixes every draw.

### 6.3 The nightly run

The runner (Phase 5 of the plan, designed and written by Ruben:
`doc/refactoring_notes/phase5_workflow_gates_and_nightly_runner_handover.md`)
does the following on the dedicated machine:

1. Takes a lock, so two runs never overlap.
2. Builds a known commit with `./compile_clean.sh`, which
   also runs the fast gate. A failed gate stops the night.
3. Runs `ctest -L library` (the library suites may run side by side), then
   `ctest -L platform` where the machine has the capability. High-level tests
   are excluded and run only by an explicit `ctest -L highlevel` command.
4. Writes a dated summary outside `build/` with the commit, host, compiler,
   the status and time of every entry, every metric against its floor, and
   the tail of every failing log. It appends to a history file, so that a
   regression can be dated and a growing runtime seen.
5. Reports briefly where the team looks.

Every entry keeps its CTest `TIMEOUT` and the whole run must fit the night.
An entry that grows past its share is reported, and trimming it is a reviewed
change, as for the fast gate.

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
| duplicates of a route or of a gate sub-suite | `starfile`, `binoris_test`, `binoris_io_test`, `imgfile`, the standalone `mini_stream`, `gui_assembler`, `gui_metadata`, `project_merge` and `coarrays`, `clustering`, `multinomal_test` |
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
| `abinitio2D_stream` | deleted; `abinitio2D` runs nightly in `simulated_workflow_6vxx` and `simulated_workflow_1jxy` |
| `angres` | the program `simple_exec prg=measure_projspace_angres nspace=<n>`; the table of values is a comment above `find_angres` |
| `ansi_colors` | `string` (`unit_core`) |
| `assign_optics` | `optics assignment` (`lib_stream`) |
| `atomfit` | deleted, with `atoms%fit_bfactors`, which nothing called |
| `atoms_stats` | `nanoparticle atoms` (`lib_single`); high-level route `simple_test_exec test=single_atoms_stats` |
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
| `coarrays` | the capability-gated platform entry `coarrays`; two images transfer one known integer |
| `continuous_3D_pcg_reconstruction` | `observation noise` (`unit_reconstruction`) and `PCG half-set` (`lib_reconstruction`) |
| `continuous_inplane_cc_grad`, `continuous_inplane_hybrid_grad`, `continuous_inplane_rotation2D_stage1_validation`, `continuous_inplane_rotation2D_route_identity` | `continuous in-plane` (`unit_pftc_align2D3D`) |
| `continuous_inplane_refine3D` | `refine3D in-plane state` and `continuous in-plane` (`unit_pftc_align2D3D`) |
| `continuous_inplane_rotation2D`, `continuous_inplane_rotation2D_metadata` | deleted (a shell driver; a post-run scan of a user project) |
| `corrs2weights`, `corrs2weights_test`, `rank_weights` | `statistics` (`unit_numerics`) |
| `create_gain`, `search_gain_flips` | `motion gain` (`unit_project`) |
| `ctf`, `ctf_test` | `CTF` (`unit_image`) |
| `detect_atoms`, `simulate_nanoparticle` | `nanoparticle atoms` (`lib_single`) |
| `detect_calpha` | `C-alpha finder` (`unit_single`) |
| `detect_calpha_molecules` | deleted; quantitative synthetic coverage remains in `C-alpha finder` (`unit_single`) |
| `diff_map_graphs` | `diffusion-map graphs` (`unit_numerics`) |
| `discrete_stack_io`, `stack_io` | `stack I/O` (`unit_core`) |
| `eigh`, `eigh_test` | `linear algebra` (`unit_numerics`) |
| `eo_diff`, `opt_lp` | deleted (needed refine3D volumes or a download; asserted nothing) |
| `eul_prob_tab2D_io` | `2D probability table I/O` (`unit_pftc_align2D3D`) |
| `extr_frac` | `decay schedules` (`unit_numerics`) |
| `flex_gpu` | the platform entry `flex_gpu` |
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
| `mini_stream` | `simple_test_exec test=mini_stream`, by hand (needs a user's movies) |
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
| `pcg_frac_update`, `rec3D_backends` | by hand, `simple_test_exec test=pcg_frac_update`, `test=rec3D_backends` |
| `pcg_recon` | the workflow entry `pcg_recon` |
| `pdb2mrc` | `pdb2mrc` (`lib_single`) |
| `phase_rand_fsc` | deleted; `simple_exec prg=fsc` |
| `phasecorr`, `ptcl_center`, `rotate_ref`, `uniform_euler`, `uniform_rot`, `eval_polarftcc` | deleted; the image basics `ptcl_center` touched are in `image` (`unit_image`), `rotate_ref_8` is in `polar correlation` (`unit_pftc_align2D3D`) |
| `phshift_policy`, `ui_visibility` | `UI visibility` (`unit_ui`) |
| `phshift_star` | `STAR project` (`unit_project`) |
| `pick_extract` | `pick and extract` (`lib_stream`) |
| `pose_cont_refine3D_adapter` | `pose adapter` (`unit_cart_align3D`) and `pose 1JYX recovery` (`lib_cart_align3D`) |
| `pose_cont_refinement` | `pose refiner` (`unit_cart_align3D`) |
| `preproc` | the high-level workflow entry `stream_preproc` |
| `project_merge` | `project merge` (`unit_project`) |
| `projdir_accumulator` | `class-average accumulator` (`unit_reconstruction`) |
| `qsys_ctrl`, `qsys_env` | `qsys control`, `qsys environment` (`unit_parallel`) |
| `rec3D_backend` | `rec3D backend` (`unit_reconstruction`) |
| `reproject` | the workflow entry `simulate_particles` |
| `serialize` | `image serialisation` (`unit_image`) |
| `sieve_cavgs` | `particle sieve` (`unit_project`) |
| `sigma2_state` | `sigma2 state` (`unit_pftc_align2D3D`) |
| `simulate_particles` | the workflow entry `simulate_particles` |
| `simulated_workflow` | the workflow entries `simulated_workflow_6vxx`, `simulated_workflow_1jxy` |
| `single_workflow` | the workflow entry `single_workflow` |
| `socket_client`, `socket_comm_distr`, `socket_io`, `socket_server` | deleted with the socket modules; `IPC TCP socket` (`unit_ipc`) tests the live transport |
| `sp_project` | `project records` (`unit_project`) |
| `starfile`, `starfile_test` | `STAR file` (`unit_project`) |
| `stream_initial_analysis` | deleted |
| `stringmatch` | `string` (`unit_core`) |
| `subproject_distr`, `ptcls_ppca_subproject_distr` | deleted with the subproject scheduling code, which had no other caller |
| `trail_rec_blend` | `trailing-reconstruction blend` (`unit_image`) |
| `ui_hash_test` | `UI hash` (`unit_ui`) |
| `units` | `simple_test_exec test=units` (all thirteen area suites, not a CTest entry) |
| the `unit_<area>` suites | unchanged names, CTest label `fast` |
