# A uniform, two-tier test environment for SIMPLE

Date: 2026-09-22

Status: in progress. Phases 0, 1 and 2 are complete (2026-09-22): the fast
gate is seven `fast` area suites run by every `compile_*.sh --compile-tests`
under a 30 s budget, at 3.3 s real on the reference Mac in Debug, with the
process-count ratchet armed and every suite passing in both table orders.
Phase 3 (the review of everything else) is next. This is a large project
with four workstreams (section 1.1), delivered in slices that are each
useful on their own.

Validation level: static source inspection of SIMPLE, and of X's test system
(`production/CMakeLists.txt`, AGENTS.md "Test admission" and "Test tiering",
`docs/design_notes/completed/test_suite_refactoring.md`) as the reference
design. No SIMPLE source was compiled and no test executable was run while
preparing this plan; the timing inventory in Phase 0 is the first thing that
will change that.

This is the single living design record for the test-environment refactor.
Update it with inventory decisions, completed migrations and validation
evidence as work lands instead of creating parallel plans.

## 1. Objective

SIMPLE gets two test gates with one front door.

**The fast gate** runs as part of the build whenever `--compile-tests` is on:
`compile_*.sh --compile-tests` builds, installs, then runs `ctest` on the
`fast` label, and the whole label finishes within **30 s of wall time on the
reference Mac** (24 cores; CI runners get a timeout multiplier, not a bigger
budget). Every test in it can fail on its own assertions, runs in-process
inside a fused per-area suite, needs no network, no download and no external
data, and pins itself to one OpenMP thread. This is the everyday loop.

Its content is decided (owner decision, 2026-09-22): the fast gate is the
unit suites that `simple_test_units` runs today, 39 sub-suites (the union of
its two routes) over the `*_tester` modules and the `test_*` procedures of
core types, split into a handful of area suites. Nothing else in the tree is
a fast-tier candidate.
Another test joins only by passing the admission rules of section 5.1 through
a review verdict, never by default.

**The extensive gate** runs overnight on a dedicated machine and has two
kinds of entries. **Library suites** group the tests of one coherent part of
the library (Fourier transforms, geometry, masks, numerics, optimisation,
statistics, I/O, projects, search, ...) into one process each, built from
the standalone tests that survive the review; they can take minutes, use
real-sized fixtures, and must fail on their own assertions. **Workflow
gates** are Ruben's self-contained simulated workflows (`simulated_workflow`
on the embedded 6VXX/1JYX systems, `single_workflow`, `mini_stream`, the
stream suite, `pcg_recon`, `rec3D_backends`, ...), which
generate their own data from known atomic models. Because the truth is known,
these workflows can be gated on it (section 5.2) without the real-data
validation registry and blessed baselines that X's `validate` needs; that
machinery is too costly for SIMPLE and is a non-goal (section 17).

**One front door** serves both:

```text
simple_test_exec test=<suite-or-test> [arguments]
```

The standalone `program simple_test_*` sources under `production/tests` are
converted into callable procedures in grouped test modules and dispatched
through area test commanders. The test executable and every test-only source
stay out of production builds (`BUILD_TESTS=OFF`, the default of every
`compile_*.sh` since 2026-09-21).

Uniformity alone does not give a build-time gate. Three facts about the
current tree shape the plan: nothing today bounds what a test run costs;
every test is its own process, which is the direction X had to reverse after
reaching 172 CTest processes; and more than half of the existing tests cannot
fail, because they contain no assertion, no `error stop`, no `THROW_HARD` and
no `tests_failed` check (section 4.4). The tiering therefore comes first, the
fast tier has an admission rule, and quality and performance of the fast-tier
tests is a workstream of its own.

### 1.1 Workstreams

| | Workstream | Depends on | Value on its own |
|---|---|---|---|
| A | Fast gate scaffolding: `ctest` in the compile scripts, labels, working directories, thread pinning, timeouts, the budget check, the review dossier script, registration of the tests that already qualify | nothing | yes: a build-time gate over today's assertion-bearing tests |
| B | Review and unification: per-area review verdicts (section 9), then standalone programs into grouped modules and fused suites behind `simple_test_exec`; overlaps merged, dead tests deleted; the glob removed | A for the registration shape and the dossiers | yes, per area |
| C | Fast-tier performance and hermeticity: split `units` into area suites, reconcile its two routes, time each sub-suite, fit the 30 s budget, confine or move the sub-suites that use sockets, HTTP or child processes | Phase 0 timing of `units` | yes: the gate gets faster and more local with each step |
| D | Extensive tier: the library suites assembled from the review's survivors, simulation-truth gates on the simulated workflows, the overnight runner, result archiving | B for each area | yes, per library suite and per gated workflow |

Phase 0 and A come first. B, C and D proceed by test area and can
interleave.

## 2. What SIMPLE takes from X, and what it does not

X's `xlms_test` is, by its own comment, "SIMPLE's simple_test_exec pattern":
every unit test a commander, dispatched by name. On top of that X added, and
SIMPLE adopts:

- **Fused suites.** The fast checks run as about ten in-process
  bounded-context suites (`xlms_test unit_io`), each accumulating ordinary
  failures and stopping once at the process boundary. Focused selectors stay
  available for debugging. Only genuine process-death contracts, the binary
  smokes and resource-hungry cases keep their own process.
- **A ratcheted process budget.** `XLMS_CTEST_BUDGET` must equal the
  registered count or configuration fails; the number only goes down (172 to
  22). "Checks are free; processes are budgeted."
- **A time budget on the label.** The regression label stays under 15 s;
  the everyday loop is `ctest`, run by `compile_debug.sh` and `recompile.sh`
  right after the build. (X also asks any entry above 0.5 s in Release to
  say why; SIMPLE does not adopt that one, section 5.1 rule 4.)
- **Hermetic tests.** No outcome depends on free memory, core count, load,
  wall clock, network or another test's files; every entry owns its working
  directory; ordinary work runs on one OpenMP thread.
- **Abort tests only where dying is the claim**, accepted on the dedicated
  hard-error status, never on "any non-zero exit".

SIMPLE does not adopt X's `validate` (a registry of real datasets and arms run
through the production paths, compared against platform-keyed blessed
baselines with a report directory and gated blessing). The datasets, the
curation and the blessing discipline are more than SIMPLE can carry. The
extensive tier uses simulated data with known truth instead (section 5.2).

## 3. Preservation contract

The refactor must preserve:

1. Every supported test assertion, fixture, argument, deterministic seed,
   expected failure, output artifact and success/failure exit status. A test
   that has no failure path has no behaviour to preserve in this sense: it
   either gains assertions (workstream C) or is reclassified as extensive or
   manual; it is not admitted to the fast gate as it stands.
2. `simple_test_exec test=list` as the discoverable catalogue of suites and
   tests.
3. The ability to execute one named test independently, in its own process,
   for debugging.
4. Process isolation where an existing mother suite deliberately launches
   child cases, where a test's contract is that the process dies, or where a
   platform launcher (`cafrun`, `mpirun`) requires a separate process.
5. Conditional availability for coarray, MPI, OpenMP offload, OpenACC, CUDA,
   network and external-data tests.
6. Existing production-library behaviour. This is a test-architecture change,
   not a numerical, scientific or workflow refactor. Adding an assertion to a
   test is in scope; changing what the code under test does is not.
7. `BUILD_TESTS=OFF` excluding `simple_test_exec` and every test-only source
   from production builds (in place since 2026-09-21).

The migration must establish one authoritative implementation for each test.
Temporary wrappers may exist within a migration phase, but duplicate program
and commander bodies must not remain as the completed state.

## 4. Current-state audit

Four terms are used with fixed meanings throughout:

- **Canonical test identity**: one `test=<name>` value. There are 150. An
  identity may be implemented on one or both routes.
- **Route implementation**: one implementation of an identity, either a
  standalone `program simple_test_<name>` under `production/tests` or a
  commander procedure reached by `simple_test_exec test=<name>`. There are
  203 (107 standalone, 96 commander); 53 identities have both.
- **Unit sub-suite**: one `begin_test_suite` / `end_test_suite` group inside
  `units` (39 in the union of its two routes). Sub-suites are not identities
  and are not registered individually.
- **Registered CTest process**: one `add_test` entry. This is what the
  process budget counts; it is neither an identity nor a sub-suite.

All headline counts below are regenerated by `scripts/test_review_dossier.py`
from the inventory and are per identity unless stated per route.

At the 2026-09-22 baseline, static source inspection finds:

| Test representation | Count | Current ownership |
|---|---:|---|
| Standalone route implementations (`production/tests/simple_test_*.f90`) | 107 | Auto-globbed into separate executables and CTest entries by `production/CMakeLists.txt` |
| Commander route implementations (`simple_test_exec test=` dispatch cases) | 96 | Fourteen topic routers under `src/main/exec`; 95 `exec_test_*` procedures plus `exec_volume_shape_descriptors` |
| Canonical test identities | 150 | 107 + 96 − 53 on both routes |
| Identities with two route implementations | 53 | Two independently maintained implementations or wrappers |
| Reusable `*_tester.f90` modules | 24 | Scattered alongside the production domains they test |

The existing `simple_test_exec` route is:

```text
production/simple_test_exec.f90
  -> test UI registry
  -> topic execution router
  -> test commander type
  -> test body or reusable tester procedure
```

The standalone route is:

```text
production/tests/simple_test_<name>.f90
  -> program body
  -> optional suite-specific helper modules or *_tester modules
```

The standalone programs are discovered with a filename glob. Each is compiled,
linked, installed and registered separately with CTest, with no `TIMEOUT`, no
`WORKING_DIRECTORY`, no labels and no thread pinning. CI runs 30 of them by
name plus two `simple_test_exec test=...` calls, all sequentially.

### 4.1 Duplication and drift

Exact-name overlap does not guarantee equivalent behaviour. The standalone
`simple_test_units` and the `test=units` implementation in
`simple_commanders_test_class` maintain separate import lists and suite
schedules and have already diverged: the standalone program includes suites
the commander version does not, while the commander version has its own
additions and terminal message. Several other pairs contain copied program
bodies. A correction made in only one path leaves the other stale.

### 4.2 Registration overhead

The 96 `simple_test_exec` cases are represented in three parallel layers: UI
program construction under `src/main/ui/simple_test`, commander types and
implementations under `src/main/commanders/test`, and topic-specific
`select case` routing under `src/main/exec`. The UI and execution registration
remain useful because tests accept typed SIMPLE arguments and are discoverable
through `test=list`. The problematic part is treating each small test as its
own commander object when a coherent area commander could call a test-module
procedure directly.

### 4.3 Tests are not all the same kind

The current sources include small deterministic unit tests; numerical and
file-format regressions; multi-stage workflow tests; mother suites that launch
isolated child cases; performance and I/O benchmarks; coarray, MPI,
accelerator, SIMD and OpenMP platform tests; socket client/server roles; tests
requiring an external reference volume (11 sources read `vol1=`, `stk=` or
`projfile=`); and tests that download data or need other external services.
Uniform entry does not mean forcing these into one process or pretending they
have identical requirements. The tiers in section 5 make the requirements
explicit.

### 4.4 Failure paths

Only 36 sources use `simple_test_utils` (`assert_*`, `begin_test_suite`,
`report_summary`), and only 20 turn `tests_failed` into a non-zero exit. The
rest fail by `error stop` (31 standalone programs) or `THROW_HARD` (which ends
in `error stop 1`), or not at all. Per route, 61 of the 107 standalone
programs and 57 of the 96 commander procedures contain no assertion, no
`error stop`, no `THROW_HARD` and no `tests_failed` check; per identity, 84
of the 150 have no failure path on any route. They print, and they pass
whenever they do not crash. Examples: `angres`, `corrs2weights`, `eigh`, `ft_expanded`,
`gencorrs_fft`, `io`, `io_parallel`, `lbfgsb`, `mask`, `neigh`, `otsu`,
`pca_all`, `serialize`, `sym`, `uniform_euler`, `uniform_rot`. A build-time
gate made of these would be green without meaning anything. Workstream C
exists for them.

`simple_end` is not a terminator: it prints a banner and optionally touches a
file, then returns. The terminators that matter for a fused suite are
`error stop`, `THROW_HARD` and the bare `stop` on the `test=list` path of
`simple_test_exec`.

### 4.5 Timing

Unknown. No per-test run time has been recorded. Phase 0 records it.

### 4.6 `units`, the fast gate in embryo

`simple_test_units` (and its commander twin, `test=units` in
`simple_commanders_test_class`) already has the shape the fast gate needs:
one process, `seed_rnd`, its own working directory, and its sub-suites (38
on the standalone route, 36 on the commander route, 39 in the union) each
run inside `begin_test_suite` / `end_test_suite` over the `*_tester` modules
(`string`, `syslib`, `fileio`, the hashes, `linked list`, `cmdline`, `ori`,
`oris`, `rec_list`, `starfile`, `project merge`, `motion gain`, ...) and the
`test_*` procedures of core types (`image`, `imghead`, `ftiter`,
`ftexp_shsrch`, `online_var`, `bspline_smoother`, `aff_prop`, `hclust`,
`srchspace_map2D_io`), plus `validate_ui_json`, then `report_summary` and
`error stop 1` if anything failed. The two routes have diverged: the
standalone runs `class compatibility`, `particle sieve` and
`2D search-space map I/O`, which the commander does not; the commander runs
`atoms`, which the standalone does not. Four sub-suites (`IPC TCP socket`,
`HTTP POST`, `forked process`, `persistent worker server`) use sockets or
child processes and need checking against the hermeticity rule. It is one
CTest entry today (`units`, called by CI as `simple_test_units`), so a
failure anywhere fails the whole thing and `ctest --parallel` cannot spread
it.

Measured 2026-09-22 (Debug, one OpenMP thread, reference Mac, all sub-suites
passing): **21.6 s standalone, 21.5 s through `simple_test_exec`**, so the
front door costs nothing measurable and the two routes are equivalent apart
from the four sub-suites they do not share. (An earlier run the same day
gave 47 s and 71 s; the machine was loaded, and the difference between the
routes it suggested was noise.) `simple_test_utils` now times each sub-suite
(commit 045c6ba39); the per-sub-suite times are:

| sub-suite | s | note |
|---|---:|---|
| forked process | 12.85 | 28 assertions; the time is `c_usleep` polling at `FORK_POLL_TIME` (100 ms) around real child processes |
| Fourier shift search | 4.37 | `test_ftexp_shsrch`: 240-pixel box, 50 trials, 20 noisy images, at `-O0`; 2 hard checks |
| image | 0.87 | |
| orientation data | 0.73 | |
| HTTP POST | 0.56 | localhost |
| orientation collection | 0.34 | 155 checks |
| particle sieve | 0.32 | standalone route only |
| class compatibility | 0.20 | standalone route only |
| UI JSON | 0.18 | |
| every other sub-suite | < 0.15 | 30 sub-suites, 1.2 s together |

Two sub-suites were 80% of the time. Without `forked process` the whole of
`units` is 8.7 s single-threaded; without both it is 4.3 s. The 30 s budget
is therefore met already on one thread, and the area split is about failure
locality and parallel spread rather than about fitting the budget. Both
decisions were taken the same day (owner):

- **`forked process` is not part of the build.** It spawns real children and
  polls them on a clock, which is what the hermeticity rule excludes. It
  keeps its assertions and moves to the `platform` label in the nightly run
  (section 5.3); in Phase 1 it is still inside the provisional `units` entry,
  and Phase 2 takes it out.
- **`test_ftexp_shsrch` runs on a 128-pixel box** (was 240, `SQRAD` scaled
  60 to 32; `TRS`, `HP`, `LP` unchanged). The change also removed two things
  the test did without checking anything: 20 noisy images generated and never
  used, and `profile_corrs`, a benchmark that filled and FFT'd two 4096-pixel
  images to print CPU times. Committed 2026-09-22 and measured the same day:
  4.37 s before, **0.30 s after**, both routes, checks passing.

## 5. Tiers and admission rules

Every test identity in the inventory (section 8) gets exactly one tier. The
tier decides how the test is registered, where it runs and what it may cost.

### 5.1 Fast tier (`fast` label, part of the build)

Admission rules, all of them:

1. **It can fail.** The test makes at least one assertion through
   `simple_test_utils` (or an equivalent typed check) whose failure reaches the
   process exit status. "Completed without crashing" is admitted only for a
   closed set of binary smokes registered with a `PASS_REGULAR_EXPRESSION`
   (X has four; SIMPLE should need no more than a handful: one `simple_exec`
   program with a pass pattern, one `simple_private_exec`, one
   `single_exec`, one `simple_stream` invocation that starts and stops).
2. **It is hermetic.** No network, no download, no user-supplied file, no
   dependence on core count, load, wall clock or another test's files.
   Fixtures are generated in the test from a fixed seed, or committed to the
   repository and small.
3. **It runs in-process** inside its area suite (section 6) on one OpenMP
   thread, and returns on success. It does not `stop`, does not change the
   working directory without restoring it, and does not leave global state
   (random-number generator, module variables, open units) that the next test
   in the suite can see.
4. **It is cheap.** The whole `fast` label finishes in 30 s of real time on
   the reference Mac with `ctest --parallel`, in every build type. The 30 s
   is a ratchet: it does not go up, and additions that would exceed it are
   paid for by shrinking fixtures or moving something out. There is no
   per-entry time rule: X's "an entry over 0.5 s says why" guarded a process
   count that SIMPLE's `SIMPLE_CTEST_BUDGET` already guards, and it would
   measure the front door rather than the tests (a `simple_test_exec`
   process pays about 1 s at `-O0` for the test UI registry before its first
   check; 2026-09-22: `unit_numerics`, whose sub-suites take 0.1 s, runs in
   1.07 s). Instead `ctest_budget.py` writes the per-entry table beside the
   log on every build, so a suite that grows is seen when it grows.
5. **It has a `TIMEOUT`** (default 60 s) and its own `WORKING_DIRECTORY` under
   `build/test_runs/<suite>`.

Members: the sub-suites of `units` (section 4.6), regrouped into area suites
along these lines, to be settled by the Phase 0 timing of each sub-suite:

| suite | sub-suites of `units` | Debug, 1 thread (2026-09-22) |
|---|---|---:|
| `unit_core` | string, syslib, fileio, character hash, hash, value-reference hash, linked list, record list, command line | 0.2 s |
| `unit_ori` | orientation, orientation collection, orientation data, Euler shift | 1.2 s |
| `unit_image` | image, image header, Fourier iterator, Fourier shift search, B-spline smoother 2D and 3D | 1.3 s (Fourier shift search 0.30 s after the trim) |
| `unit_numerics` | online variance, multinomial random draw, straight-line fit, affinity propagation, hierarchical clustering | 0.1 s |
| `unit_project` | STAR file, project merge, class compatibility, particle sieve, 2D search-space map I/O, motion gain, atoms | 0.7 s |
| `unit_ui` | UI JSON, GUI metadata, GUI assembler | 0.2 s |
| `unit_ipc` | IPC TCP socket, HTTP POST, persistent worker server, persistent worker message — localhost only, bounded; `forked process` is excluded by decision and goes to `platform` | 0.6 s |

Measured through the gate on 2026-09-22 (Debug, `ctest -j12`, one thread
per entry): all seven pass, **3.3 s real**, 12.9 processor-seconds;
`unit_ori` 3.3 s, `unit_image` 2.6 s, `unit_ipc` 1.8 s, `unit_project`
1.7 s, `unit_ui` 1.25 s, `unit_core` 1.2 s, `unit_numerics` 1.1 s. About
1 s of every entry is the front door (the test UI registry at `-O0`), which
was paid once when `units` was one process and is paid seven times now;
the parallel spread more than covers it. The budget is comfortable, and the
process count for the ratchet is seven.

No other test identity in the tree is a fast-tier candidate by default. The
remaining 149 identities are reviewed for the extensive tier, manual use or
deletion (section 9); one of them joins the fast gate only through a `keep`
or `modify` verdict that states which admission rules it meets and what it
measured in Phase 0.

### 5.2 Extensive tier (`library` and `workflow` labels, overnight)

The extensive tier runs nightly on the dedicated machine as
`ctest -L "library|workflow"`. It has two kinds of entries.

#### 5.2.1 Library suites (`library` label)

A library suite is the tests of one coherent part of the library, run in one
process through `simple_test_exec test=lib_<area>`, one CTest entry per
suite. It is the extensive-tier counterpart of the fast area suites: the same
fused shape (each member inside `begin_test_suite` / `end_test_suite`,
`report_summary` at the end, non-zero exit on any failure, a focused selector
for one member), without the 30 s budget. Members are the standalone tests
the review keeps (`demote to lib_<area>`), so a numerical test that takes a
minute on a realistic box has a home and runs every night instead of never.

Admission rules for a library suite member:

1. it can fail: at least one assertion whose failure reaches the process
   exit status (a `demote` of a print-only test carries the assertion it
   must gain, section 9.4);
2. it runs unattended: no user-supplied file, no download, fixtures generated
   from a seed or committed; a member may take minutes and use full-sized
   boxes, but it must finish;
3. it is deterministic across runs on one machine (declared seed) and
   restores the working directory and any module state it changes, since it
   shares a process with its suite;
4. it is registered with a `TIMEOUT` and the suite's total is recorded in the
   nightly summary, so growth is visible.

Provisional suites, from the router areas and the standalone programs (to
be settled by the review):

| suite | drawn from |
|---|---|
| `lib_fft` | the `fft` area: FFT, Fourier-space operations, gencorrs, polar FTCC, rotations, expanded FT |
| `lib_geometry` | `geometry`: symmetry, Euler sampling, uniform rotations, angular resolution |
| `lib_masks` | `masks`: masking, envelopes, mask bounds, graphene, Otsu |
| `lib_numerics` | `numerics`: linear algebra, eigensolvers, PCA, KPCA, clustering, random draws |
| `lib_optimize` | `optimize`: L-BFGS-B, low-pass optimisation, shift search, LP limits |
| `lib_stats` | `stats`: correlations, weighting, rank statistics, sigma estimation |
| `lib_io` | `io`: image and stack I/O, MRC validation, STAR import/export, parallel I/O |
| `lib_project` | `class` and `utils`: project files, orientation documents, binoris, serialization |
| `lib_search` | the search-strategy and continuous-refinement suites (`continuous_inplane_*`, `pose_cont_refinement`, `continuous_3D_pcg_*`) once they run on generated data |
| `lib_parallel` | `parallel`: OpenMP correctness with an explicit small team |

Registration: one entry per suite, `LABELS library`, its own working
directory, `OMP_NUM_THREADS` set explicitly (a suite may use a team, since
the library suites run before the serial workflow gates and can be scheduled
against each other), a `TIMEOUT` sized to the suite. A library suite is not
a workflow: it does not start `simple_exec` or distributed workers; a test
that needs them is a workflow gate.

#### 5.2.2 Workflow gates (`workflow` label)

The workflow gates are the simulated workflows, gated on the truth they were
simulated from. Today `simulated_workflow`, `single_workflow`, `mini_stream`
and the stream suite check that the pipeline completes: files exist,
counts match, the heartbeat is well-formed, `abinitio2D` produced classes.
They do not compare the result with the model that generated the data. Since
the data comes from embedded atomic coordinates (6VXX, 1JYX) with known
orientations, defocus and B-factor, the comparison is available at no data
cost. Each workflow declares floors such as:

- resolution: FSC=0.143 between the reconstructed map and the ground-truth
  map simulated from the same coordinates, at or better than a declared
  value for that system, box and particle count;
- poses: fraction of particles whose recovered orientation is within a
  declared angular distance of the simulated one, after symmetry and
  hand alignment; shift error likewise;
- workflow structure: particle, class and state counts, heartbeat
  completeness, project-file consistency (the checks that exist today);
- stream: the same, per stage, over the simulated movie set.

Floors are declared in the test, versioned with it, and loosened only with a
written justification in the commit. There are no blessed baselines to
maintain, no platform keys and no report registry; the run archives its
summary (commit, host, compiler, per-workflow metrics and pass/fail) into a
directory on the dedicated machine so a regression can be dated.

Registration: one CTest entry per workflow, `LABELS workflow`,
`RUN_SERIAL TRUE` (each owns the machine's OpenMP team and may start
distributed workers), a long `TIMEOUT`, its own working directory. Expected
members: `simulated_workflow` (both systems, both pickers),
`single_workflow`, `mini_stream`, the stream suite, `pcg_recon`,
`pcg_frac_update`, `rec3D_backends`, `reproject`, and the nano workflows (`atoms_stats`,
`detect_atoms`, `detect_calpha*`, `simulate_nanoparticle`).

#### 5.2.3 The nightly run

`ctest -L "library|workflow"` after a clean `--compile-tests` build, started
by cron or a scheduler on the dedicated machine. The library suites run
first (in parallel, each pinned), then the workflow gates serially. The
summary (commit, host, compiler, per-suite and per-workflow times, metrics
against floors, pass/fail) is archived into a dated directory on that
machine and mailed or written where the team looks. The whole run must fit
the night; a suite or workflow that grows past its share is reported by the
summary, and trimming it is a reviewed change, as for the fast gate.

### 5.3 Isolated and platform tests (`platform` label)

Coarray (`cafrun -np 2`), MPI, OpenMP offload, OpenACC, CUDA, socket
client/server pairs, the `forked process` suite (real child processes,
clock-based polling) and expect-abort tests keep their own process. They are
registered only when CMake has confirmed the capability and launcher, carry
the `platform` label, and are excluded from the fast gate unless a specific
entry fits the budget on the reference Mac and is hermetic (the coarray test
currently runs in its own CI job and stays there). Socket tests must be
bounded and self-terminating or be driven by an orchestrating test that
starts both roles and checks both exit statuses. Expect-abort tests are
accepted on the SIMPLE hard-error status (`error stop 1` from
`simple_exception`), never on "any non-zero exit", and exist only where dying
is the contract.

### 5.4 Not registered

Benchmarks (`io_parallel`, `simd`, `openmp` as timing tools), tests that
download data, and tests that need a user-supplied volume or stack remain
runnable by name through `simple_test_exec` and are documented as manual.
They are not CTest entries.

## 6. Target architecture

```text
simple_test_exec
  -> test UI metadata and argument parsing
  -> topic execution router
  -> area test commander (suite: runs every test of its area in-process,
     accumulates failures, stops once at the process boundary;
     focused selector: runs one named test)
  -> grouped callable test module
       |-- focused test procedures
       |-- suite-owned fixtures, built once per process and handed over
       `-- reusable domain tester modules
  -> production SIMPLE APIs
```

| Layer | Owns | Does not own |
|---|---|---|
| `simple_test_exec` | Process front door, command parsing, timing, logging, memory monitoring, final process status | Individual test algorithms |
| Test UI | Suite and test names, arguments, help, required inputs, defaults | Test execution |
| Topic router | Routing a registered suite or test to one area commander | Test bodies |
| Area test commander | Running its area's tests in one process, shared fixtures and setup, failure accumulation, workflow orchestration in the extensive tier | Low-level assertions or copied test algorithms |
| Grouped test module | Callable test bodies and closely related suite helpers | CLI front-door behaviour or production orchestration policy |
| Existing `*_tester` or domain test APIs | Reusable focused checks close to the type or subsystem under test | Executable lifecycle |
| CMake/CTest | Build gating, launcher selection, labels, timeouts, working directories, thread pinning, the budgets | A second implementation of test behaviour |

### 6.1 Grouped callable modules

Related tests are grouped into modules named for their domain. A provisional
grouping, to be settled by the Phase 0 inventory:

```text
simple_test_cases_core
simple_test_cases_fft
simple_test_cases_geometry
simple_test_cases_io
simple_test_cases_masks
simple_test_cases_numerics
simple_test_cases_optimization
simple_test_cases_parallel
simple_test_cases_project
simple_test_cases_reconstruction
simple_test_cases_stream
simple_test_cases_workflows
simple_test_cases_utils
```

Large numerical suites that already have a mother module and several focused
helper modules keep that cohesion rather than being pasted into a monolithic
topic file. A grouped module is private by default and exports only its
callable entries:

```fortran
module simple_test_cases_io
use simple_cmdline, only: cmdline
implicit none
private

public :: test_imgfile
public :: test_mrc_validate
public :: test_stack_io

contains

subroutine test_stack_io(cline)
    class(cmdline), intent(inout) :: cline
    ! Migrated standalone test behaviour, returning on success.
end subroutine test_stack_io

end module simple_test_cases_io
```

A common dummy argument is not imposed where it adds nothing: a test with no
inputs stays a no-argument procedure; tests that need CLI values take
`cmdline` at the pre-parse boundary or, preferably, typed `parameters` after
their area commander has initialised them; tests that exercise one production
object may receive it directly.

### 6.2 Area test commanders and suites

One commander per coherent test area, not one per test:

```fortran
type, extends(commander_base) :: commander_test_io
contains
    procedure :: execute => exec_test_io
end type commander_test_io
```

`test=unit_core` runs its sub-suites in one process: it builds any shared
fixture once, calls each callable procedure inside a `begin_test_suite` /
`end_test_suite` pair, and after the last one calls `report_summary` and exits
non-zero if anything failed, exactly as `units` did for all of them at once.
(Implemented in `simple_commanders_test_class`: each area is a table of
`unit_suite` entries, a name and a no-argument procedure, built by
`suites_<area>` and run by `run_unit_suites`; procedures that take arguments
are wrapped.) The area suites registered with CTest are the authoritative
gate.
`test=units` remains as a developer convenience that runs every fast suite in
sequence in one process; it is not authoritative, and a failure it shows that
the CTest gate does not is a state leak between suites, to be fixed in the
suite that leaks. A focused selector runs one sub-suite for debugging:
`test=unit_core suite=hash`, where `suite` is an optional string input
registered on each area program in the test UI (its help lists that area's
sub-suite names, lowercase with underscores) and read from `cline` in the
area commander, in the way `simulated_workflow` reads `system`. Like
`system`, it needs a field in `simple_parameters` (`suite`), because the
command-line parser accepts only keys that the generated argument list
knows.
The area commander is what CTest registers; the focused names are what a
developer types.

For the extensive tier the area commanders are the library-suite commanders
(`test=lib_fft`, ..., section 5.2.1), one per suite and the same fused shape
as the fast ones, and the workflow commanders that exist
(`simulated_workflow`, `single_workflow`, `mini_stream`, the stream suite),
one CTest entry each.

Per-test commander types are removed as their bodies migrate. The fourteen
topic areas are a starting point; the inventory may split, merge or rename
them by real test ownership.

### 6.3 Test procedure lifecycle

Callable test procedures must:

1. return normally on success;
2. report checks through `simple_test_utils` (or an equivalent that reaches
   `tests_failed`);
3. signal failure through the accumulating path, reserving `THROW_HARD` for
   conditions that make continuing the suite meaningless (a missing fixture);
4. leave process-wide timing, Git-version output, log closure and memory
   monitoring to `simple_test_exec`;
5. not parse raw process arguments;
6. not call `stop` or any success-path terminator;
7. restore the working directory and release resources they own, including
   on early-failure paths, and reset any module-level state they changed;
8. read the fixture handed to them by the suite rather than rebuilding it, and
   build their own only when run alone.

Contained procedures of a standalone program move into the owning grouped
module or an existing suite-specific helper module, never into a commander.

### 6.4 Process isolation

Ordinary fast-tier tests share a process by design (section 5.1). Isolation
is kept for the cases in section 5.3 and for mother suites that already
launch child cases; those launch the same executable with another selector:

```text
simple_test_exec test=pose_cont_refinement case=<case-name>
```

An eventual `test=all` convenience runs the suites in-process one after
another and delegates only the isolated cases to subprocesses.

### 6.5 Specialized launchers

```text
cafrun -np 2 simple_test_exec test=coarrays
mpirun -np N simple_test_exec test=<mpi-case>
simple_test_exec test=<gpu-case>
```

CMake and CI select the launcher and register the test only when the required
capability is available.

## 7. Build and CTest design

1. **Gating.** `BUILD_TESTS` (SIMPLE's option, default ON in CMake, OFF in
   every `compile_*.sh` without `--compile-tests`) gates `simple_test_exec`,
   the test-only library sources and every CTest registration. Nothing test-
   related is built otherwise.
2. **The fast gate runs from the compile scripts.** With `--compile-tests`,
   each `compile_*.sh` runs, after `make install`:

   ```bash
   ctest --test-dir build -L fast --output-on-failure --parallel "$NJOBS" --timeout 120 \
       2>&1 | tee build/test_runs/ctest_fast.log
   scripts/ctest_budget.py build/test_runs/ctest_fast.log --budget 30
   ```

   Until Phase 2 declares the fast gate, the label expression is
   `"fast|provisional"` and `ctest_budget.py` runs with `--no-budget`, so the
   provisional `units` entry runs and is timed on every build without a
   budget it cannot yet meet. This is `scripts/run_fast_gate.sh`, called by
   every `compile_*.sh --compile-tests` after `make install` and by
   `make check`; `GATE_DECLARED` inside it is the Phase 2 switch.

   `NJOBS` is the core count divided by two (tests are pinned to one thread
   but do I/O). `ctest_budget.py` reads the `ctest` output, prints every
   entry sorted by time, fails if the run's real time exceeded 30 s or if
   any entry failed, and writes the table beside the log so the numbers are
   kept. The `check` target runs the same thing by hand.
3. **Registration by suite.** One `add_test` per fast area suite
   (`simple_test_exec test=unit_core`), per library suite
   (`test=lib_fft`), per workflow gate, per platform case and per binary
   smoke. Every registration sets `LABELS` (`fast`, `library`, `workflow`,
   `platform`), `TIMEOUT`, `WORKING_DIRECTORY` under `build/test_runs/<name>`
   (created at configure time) and an explicit `OMP_NUM_THREADS` (1 for the
   fast tier). Workflow entries set `RUN_SERIAL TRUE`.
4. **Process budget ratchet.** `SIMPLE_CTEST_BUDGET` in
   `production/CMakeLists.txt` must equal the registered count or
   configuration fails, as in X. It is armed at the end of Phase 2, once the
   area suites exist, at the count registered then; before that the count is
   allowed to change (Phase 1 registers one `units` process, Phase 2 replaces
   it with about seven). From then on it only goes down; raising it is an
   owner decision recorded here.
5. **No more standalone executables.** The `simple_test_*.f90` glob, its
   per-program `add_executable`, `install` and `add_test` are removed at the
   end of workstream B (Phase 7). Until then the glob keeps building the
   programs that have not migrated, but none of them is registered with
   CTest once Phase 1 lands (they are runnable by name and through
   `test_timing_run.sh`); CI keeps calling the not-yet-migrated ones by
   name until Phase 7 switches it to labels.
6. **Registry consistency.** A configure-time or `check`-time script compares
   the CTest registrations with `simple_test_exec test=list` and fails on a
   registered selector that is not dispatchable, or on a selector the test UI
   marks as registrable (area suites, library suites, workflow gates,
   platform cases) that is not registered. Manual tools, focused sub-suite
   selectors and platform cases whose capability is absent on this machine
   are dispatchable without being registered, by design. A generated common
   registry is a possible later improvement, not a prerequisite.
7. **Install.** A `--compile-tests` install contains `simple_test_exec` and
   no standalone test binaries.

NICE's Python tests remain under their Python runner. Python validators, shell
compatibility checks and external oracle packages are not converted into
Fortran modules. The uniformity goal applies to SIMPLE's Fortran tests.

## 8. Test inventory

Before moving code, create the migration inventory with one row per current
test identity (both routes), with at least:

| Field | Purpose |
|---|---|
| Canonical test ID | Final `test=<name>` value |
| Current sources | Standalone program, commander procedure, tester module, helpers |
| Authoritative behaviour | Which implementation or merged behaviour is retained |
| **Tier** | fast, extensive, platform, or not registered (section 5) |
| **Wall time / run state** | From Phase 0, Debug and Release, single-threaded: `measured <s>` (passed or failed on its own), `timed out`, `crashed`, `missing fixture` (refused for lack of arguments or data), `unsupported capability` (coarray, MPI, GPU or a launcher this machine lacks), `manual` (persistent server, interactive). A runtime is required only for `measured` rows; every other state is itself the Phase 0 result for that row |
| **Failure path** | assertion / error stop / THROW_HARD / none |
| **Performance action** | none, shrink fixture, share fixture, drop I/O, split, move tier |
| Target grouped module and area commander | Owners after migration |
| Arguments | Existing CLI keys, defaults, parsing path |
| Launcher/capability | Serial, child process, coarray, MPI, GPU, network |
| Fixtures | Generated, committed, downloaded, user-supplied |
| Working files | Products, cleanup, retained evidence |
| Baseline result | Exit status, success marker, assertions, tolerances |
| **Verdict** | keep, modify, merge into, demote, delete, retire, investigate (section 9.4), with reviewer and date |
| **Verdict note** | The required note for the verdict |
| Migration status | Unreviewed, baselined, callable, routed, old executable removed, validated |

The inventory is also the review record (section 9): every row carries a
verdict before its test is touched, and a **retired tests** table at the end
of the inventory lists every identity removed, with date, reason and
replacement.

The inventory lives at `doc/refactoring_notes/test_inventory.md`, one table
per area, generated by `scripts/test_review_dossier.py`; the verdict and note
columns are the only ones edited by hand and survive regeneration. It is
updated in the same commit as the verdicts it records, and
`scripts/test_review_dossier.py --coverage-after doc/refactoring_notes/test_inventory.md`
is the coverage accounting of section 9.5.

## 9. Test review

Unification and tiering move tests; they do not decide whether a test
deserves to exist. That decision is a review, and it is where most of the
overlap and dead weight will be found. The review is part of the plan, not a
side activity: no test is migrated (workstream B), given assertions
(workstream C) or gated (workstream D) before it has a recorded verdict.

### 9.1 What the review faces

Static inspection of the 150 identities (203 route implementations) gives
the starting picture:

- **Two route implementations: 53 identities.** Two implementations claiming
  the same test; the routes of `units` are the known divergent pair.
- **Overlap by what they exercise: 14 clusters covering 30 identities.**
  Comparing the production calls each identity makes (both routes pooled),
  17 pairs share at least half their footprint or one is a subset of the
  other: `ori`, `oris`, `eigh`, `maxnloc`, `starfile`, `class_sample` and
  `bounds_from_mask3D` against their `_test` twins; `io` inside
  `io_parallel`; `lbfgsb` and `lbfgsb_cosine`; `pca_all` and `pca_imgvar`;
  `cc_connectivity` and `image_bin`; `gen_pickrefs` and `preproc`;
  `reproject` and `simulate_particles`; and the `continuous_inplane_*`
  gradient tests with `eval_polarftcc`.
- **Too thin to compare: 57 identities** make fewer than three production
  calls. These are the print-only smokes and demos of section 4.4, and the
  prime candidates for deletion or merging.
- **No usable age signal.** Every test source was touched in 2026 (the test
  tree was reorganised this year), so "last modified" says nothing about
  whether a test is still meaningful. Staleness has to be judged from what
  the test exercises and whether that code path still exists.

- **The question is not "fast or extensive".** With the fast gate fixed as
  the `units` sub-suites, the review of the other 149 identities asks: which
  library suite does this belong to (it has, or can be given, a verifiable
  claim about a production path and runs unattended), is it a workflow gate,
  is it a manual tool worth keeping under its name, or does it go? The
  proposal column in the inventory is pre-filled on that basis.

The footprint comparison is a screen, not a verdict: it finds candidates,
and a person decides.

### 9.2 Unit, reviewer, order

The unit of review is the test identity. Review proceeds by area, in the
order the areas are migrated, so verdicts are fresh when the migration
happens. The reviewer is the owner of the production subsystem the area
tests (by authorship today: Ruben for `io`, `parallel`, `single`, `stream`,
`utils` and the workflow tests; Cyril for `masks`; Hans for the rest), with
Hans as the arbiter for disagreements and for every deletion. A reviewer may
review their own tests; a second person reads the batch summary.

### 9.3 The dossier

`scripts/test_review_dossier.py` generates one dossier per test identity
from the tree and the Phase 0 timing runs (`scripts/test_timing_run.sh`), so
the reviewer does not have to assemble it:

- name, route(s), source file(s), line count, author history;
- production modules imported and type-bound procedures called (the
  footprint), and which of those no other test touches (unique coverage);
- failure path (assertion, `error stop`, `THROW_HARD`, none) and what the
  test would report if it ran on wrong results;
- fixtures (generated, committed, downloaded, user-supplied), arguments,
  launcher, working files;
- measured wall time, Debug and Release;
- callers: CI, scripts, documentation, other tests;
- overlap candidates: same-name twin, footprint cluster members, subset
  relations;
- a proposed tier (section 5) from the rules, to be confirmed or overruled.

The dossiers are regenerated for each batch and are inputs, not records; the
record is the verdict.

### 9.4 Verdicts

Each test identity gets exactly one verdict, recorded in the inventory with
reviewer and date:

| Verdict | Meaning | Required note |
|---|---|---|
| `keep` | Migrate as is; tier assigned | tier |
| `modify` | Migrate with a stated change: add an assertion, shrink or share a fixture, split, pin threads, remove a workflow run | what changes, and why the test is worth it |
| `merge into <id>` | Its unique coverage is folded into the named test; this identity is deleted | which checks move |
| `demote` | Not fast-tier material; goes to a named library suite (`lib_<area>`), the workflow gates, platform or manual | destination, the assertion it must gain if it has none, and what it would take to promote it |
| `delete` | Removed without replacement | one of the reasons below |
| `retire` | Removed because the production path it tests is itself being removed | the production change |
| `investigate` | Parked: the verdict needs more than the dossier gives | owner, the open question, and the condition or date on which it is revisited; an `investigate` older than one migration phase is escalated to Hans |

A test is deleted when at least one of these holds and the note says which:
it duplicates another test's coverage entirely; it exercises code that no
longer exists or is deprecated for removal; it is a demo or benchmark with no
verifiable claim and no plausible assertion; it cannot be made hermetic and
has no extensive-tier value; or it has been broken with no caller for long
enough that nobody noticed. A test is kept when it has unique coverage of a
production path that matters, or when it is an isolation or launcher case
(section 5.3). "It might be useful someday" is not a reason to keep; the
history has it.

Deleted and retired identities go into a **retired tests** table in the
inventory (name, date, reason, replacement if any), so the same test is not
rewritten by accident and a reader of an old note or CI log can find out what
became of it.

### 9.5 Coverage accounting

Deleting tests must not silently drop coverage. For each batch the dossier
script reports the union of production modules and type-bound procedures
exercised by the area's tests before and after the verdicts are applied, and
lists every production module that loses all test coverage. That list is
part of the batch summary; each entry is either accepted with a reason
("was only exercised by a print-only smoke") or answered by a `modify` or
`merge`. This is a call-footprint proxy, not line coverage; if it proves too
coarse, an instrumented (`--coverage`) build of the extensive tier is the
stronger tool and can be added later.

The review also reads in the other direction: a test that only exercises
dead production code has found dead production code. Whether that code goes
too is an owner decision, taken separately from the test's verdict; when the
answer is yes, test and code go in one commit so neither outlives the other.
First instance, 2026-09-22: `subproject_distr` and
`ptcls_ppca_subproject_distr` were the only callers of the subproject
scheduling framework in `simple_qsys_ctrl`/`simple_qsys_env` (added with
them on 2026-04-06); both tests and the framework were removed together.

### 9.6 Mechanics and pace

The review of an area is one commit series: the verdicts (inventory rows
and the retired table) first, as reviewable text; then one commit per test
for the action taken (extract, assert, shrink, merge, delete). A verdict
takes about ten minutes with the dossier in hand; a test that needs longer is
marked `investigate` with a note and parked, not allowed to block the batch.
At that pace an area of fifteen tests is an afternoon, and the whole tree is
a few weeks of reviewer time spread across the migration.

Verdicts are revisited only through the same process: a `keep` that later
fails the budget goes back through review with its timing, not straight to
deletion.

## 10. Fast-tier performance

The 30 s budget will not be met by classification alone; the fast candidates
have to be made cheap. After the Phase 0 timing run, for every fast candidate
in descending order of wall time:

- **Shrink the fixture.** Box size, particle count, iteration count and
  number of orientations are chosen to exercise the code path, not to be
  realistic. A 64-pixel box and a few dozen particles test an FFT, a CTF or a
  correlation as well as a 256-pixel box with thousands.
- **Build once, share.** The suite builds one fixture per process and hands
  it to its tests; a test copies the immutable parts and writes only derived
  files. No fixture is cached in the build tree across runs.
- **Stay in memory.** Prefer in-memory images and stacks over writing and
  re-reading MRC files unless file I/O is what the test is about.
- **No workflow in a unit test.** A test that runs `cluster2D` or `refine3D`
  to check a utility is an extensive-tier test wearing a fast label; replace
  it with a direct check of the utility and let the workflow tests cover the
  integration.
- **No sleeps, polls or real timers.** Test scheduling and watchdog policies
  as state machines with injected time.
- **One thread.** OpenMP teams inside a `ctest --parallel` run oversubscribe
  the machine and make timings meaningless; parallel correctness is tested by
  the `parallel` suite with an explicit small team.
- **Measure again.** `ctest_budget.py` keeps the per-entry table beside the
  log on every build; a suite that grows is seen when it grows, and the 30 s
  label total is what fails the gate.

The same pass gives the print-only tests their assertions (section 4.4): the
person who shrinks a fixture is looking at what the test computes and can
state what the right answer is.

## 11. Command and naming contract

`test=` is the canonical selector:

```text
simple_test_exec test=list
simple_test_exec test=unit_core
simple_test_exec test=unit_core suite=hash
simple_test_exec test=lib_fft
simple_test_exec test=simulated_workflow system=6vxx
```

Suites are named `unit_<area>` for the fast tier; extensive workflows keep
their descriptive names. Documentation, CTest, CI and child-suite launchers
use this form. Historical `prg=` examples are normalized during migration;
whether a temporary `prg=` alias is retained is a compatibility decision, and
it must not remain the documented interface.

Test identifiers describe the behaviour under test and do not encode the
former executable shape. Existing stable IDs are retained unless misleading
or colliding. Renames require a documented alias or a coordinated update of
CI, scripts, implementation notes and user instructions.

## 12. Staged migration

| Phase | Workstream | Change | Exit gate |
|---:|---|---|---|
| 0 | all | **Timing and failure-path inventory.** Build with `--compile-tests` in Debug and Release; run `scripts/test_timing_run.sh` in each (every standalone binary and every `simple_test_exec` case, each in its own directory under a timeout, single-threaded). Run `scripts/test_review_dossier.py --timing ...` to generate the dossiers and the inventory (section 8): proposed tier, failure path, run state and time, overlap candidates, fixtures, launchers, callers. Propose the grouped-module and commander map and the area review order. | Every identity has a dossier, a run state (a measured time where it could run; otherwise timed out, crashed, missing fixture, unsupported capability or manual) and a proposed tier. |
| 1 | A | **Scaffolding and a provisional gate.** `ctest` after install in every `compile_*.sh --compile-tests`; `ctest_budget.py`; labels, timeouts, working directories, thread pinning. Register `units` as it is under the label `provisional`, not `fast`: it runs on every `--compile-tests` build and reports its time, but the budget is not enforced and nothing carries the `fast` label yet, because `units` still contains the socket, HTTP and child-process sub-suites that the fast admission rules exclude. Register the simulated workflows under `workflow` with their current checks. `SIMPLE_CTEST_BUDGET` is not yet set. | `compile_debug.sh --compile-tests` builds and runs `units` green; its per-sub-suite times are known; CI still passes. **Met 2026-09-22:** commit 8b7dfd4d7; `units` 17.5 s through the gate (Debug, 1 thread). |
| 2 | C | **Split `units` into hermetic area suites, then declare the fast gate.** Reconcile the two routes into one implementation (the union of their sub-suites), move the sub-suite lists into the grouped modules of section 6.1, move `forked process` out to its own `platform` entry (decided, section 4.6) and confirm the remaining `unit_ipc` sub-suites are localhost-only and bounded. Register one entry per area suite (section 5.1 table); when every registered suite meets the admission rules, relabel them `fast`, drop the `provisional` entry, set `SIMPLE_CTEST_BUDGET` to the registered count, and turn on the 30 s check in `ctest_budget.py`. Shrink what is over budget. Remove the standalone `simple_test_units` and its CI call. | Every area suite runs in one process and meets the admission rules; the `fast` label is under 30 s with `ctest --parallel`; the budget ratchet is armed; a failure names its suite. **Landed 2026-09-22** (`simple_commanders_test_class` rewritten as area tables over a `unit_suite` type, `suite=` input, `SIMPLE_UNIT_ORDER=reverse`, `forked_process` under `platform`, `SIMPLE_CTEST_BUDGET=19`, `GATE_DECLARED=yes`). **Met 2026-09-22:** the `--compile-tests` build passed 7/7 in 3.3 s real; every suite also passed with `SIMPLE_UNIT_ORDER=reverse`, so no sub-suite leaks state into its neighbours in either direction. |
| 3 | B | **Review everything else.** The other 149 identities, area by area (section 9): `demote` to a named library suite or the workflow gates, `keep` as manual, `merge`, `delete` or `retire`; the 53 two-route identities and 14 footprint clusters resolved to one implementation each; deletions applied with their retired-tests rows and coverage accounting. | Every identity has a verdict naming its destination; no pair or cluster retains two implementations of the same coverage. |
| 4 | B + D | **Build the library suites.** Area by area: the survivors move into the grouped module, gain the assertions their verdicts require, and are registered as one `lib_<area>` entry under `library`; standalone binaries removed as each suite completes. The first suite (`lib_fft` or `lib_geometry`) is the pilot for the fused extensive shape. | Each library suite runs in one process nightly with a recorded time; its members' binaries are gone. |
| 5 | D | **Simulation-truth gates and the nightly runner.** `simulated_workflow`, `single_workflow` and `mini_stream` compare against the generating model (FSC to the truth map, pose agreement) with declared floors; the nightly `ctest -L "library|workflow"` run and its archive on the dedicated machine. | The nightly run completes unattended and reports per-suite times and per-workflow metrics against floors. |
| 6 | B | **Mother suites, platform and socket cases** with explicit isolation and launcher policy. | Parent suites launch `simple_test_exec` children with full accounting; platform cases skip or register predictably and cannot hang the fast gate. |
| 7 | B | **Retire the glob.** Remove the standalone executable glob, switch CI to `ctest -L fast` plus the platform jobs, normalize documentation, delete stranded per-test commander types, routers and UI entries. | A clean `--compile-tests` build produces only `simple_test_exec`; registry consistency passes; CI uses no standalone Fortran test executable. |

Complete a vertical area slice (baseline, extract, assert, shrink, route,
switch callers, delete duplicate) before starting the next. Do not copy every
program into a module and leave both systems standing.

## 13. Migration rules for individual tests

For each test, in this order:

0. Confirm it has a review verdict (section 9). A `delete` or `retire` ends
   here: remove the sources, the UI entry, the router case and the callers in
   one commit, and add the retired-tests row. A `merge` moves its unique
   checks into the target before this identity is removed.
1. Record its unchanged baseline and its Phase 0 wall time before editing.
2. Identify all source and script callers by test executable name.
3. Compare any existing `simple_test_exec` implementation with the standalone
   program line by line and select or merge the authoritative behaviour.
4. If it is kept (fast or extensive) and has no failure path, write the
   assertion it was implicitly making (what would have been wrong if the
   printed numbers were wrong). A manual tool needs none.
5. Move the body into its grouped module; keep suite-specific helper modules
   that already express ownership.
6. Remove raw `get_command_argument` or `parse_oldschool` ownership from the
   body; register required keys in the test UI and consume the parsed command
   state through the established lifecycle.
7. Make the procedure return on success and accumulate failures; keep process
   lifecycle in `simple_test_exec`.
8. Apply the performance actions from the inventory row; re-measure.
9. Route it through its area suite without adding a one-test commander type.
10. Change CTest, CI, scripts and documentation to the suite or `test=<id>`.
11. Compare exit status, assertions, products, tolerances and expected-failure
    behaviour with the baseline.
12. Delete the standalone program and any duplicate commander body in the same
    completed slice.

Mechanical extraction, adding an assertion, and shrinking a fixture are three
separate commits per test, so each can be reviewed and reverted alone. None
of them changes what the production code does.

## 14. Risks and mitigations

| Risk | Level | Mitigation |
|---|---:|---|
| The fast gate is green but hollow because print-only tests were admitted. | High | Admission rule 1; the failure-path column in the inventory; Phase 1 registers only assertion-bearing tests. |
| The 30 s budget is missed and quietly raised. | High | The budget is a ratchet checked by `ctest_budget.py` on every `--compile-tests` build; raising it is an owner decision recorded here. |
| Fused suites leak state between tests (RNG, module variables, cwd, open units) and produce order-dependent results. | Medium-high | Lifecycle rule 7; during Phase 2 run each area suite in a few fixed alternate orders (its table order, reversed, and with any I/O or IPC sub-suites first) and require identical results, since the sub-suite procedures have different signatures and are not trivially shuffled; keep the focused selector so a failing test can be reproduced alone. |
| Duplicate implementations have diverged; selecting one loses coverage. | High | Review verdicts with the dossier's unique-coverage list; merge missing assertions before deleting either path (`units` is the known example); coverage accounting per batch. |
| The review deletes tests that were the only coverage of something that matters. | Medium-high | Section 9.5: every module that loses all coverage is listed in the batch summary and accepted with a reason or answered by a `modify`/`merge`; deletions go through Hans. |
| The review stalls on hard cases and blocks migration. | Medium | The `investigate` verdict parks a test without blocking its batch; ten-minute pace with the dossier. |
| Conversion from program to subroutine changes initialisation, working directory, logging or teardown. | Medium-high | Baselines; explicit lifecycle ownership; a file-producing test in the pilot. |
| A callable test calls `stop`, `error stop` or `THROW_HARD` on a path that should accumulate, ending the suite early. | Medium | Audit terminal calls in extraction; `THROW_HARD` only for missing fixtures; expect-abort tests isolated. |
| Mother suites lose crash isolation when child programs disappear. | High | Launch the same `simple_test_exec` binary as a subprocess for each isolated case. |
| Coarray/MPI/GPU tests compile but run under the wrong launcher or build. | High | Launcher and capability checks in CMake/CI; compile guards retained. |
| Simulation-truth floors are set from one run and become flaky. | Medium | Set floors with margin from several seeded runs; declare the seed; treat a floor change as a reviewed commit. |
| The overnight machine drifts (compiler, libraries) and the extensive tier fails for reasons unrelated to SIMPLE. | Medium | The archive records compiler and host; the runner does a clean build each night. |
| Grouped modules become dumping grounds. | Medium | Group by domain; keep cohesive suite helpers; split by responsibility when a file becomes hard to review. |
| CTest and `test=list` drift apart. | Medium | The registry-consistency check; stable canonical IDs. |
| `ctest --parallel` oversubscribes the machine and the budget is missed for scheduling reasons. | Medium | `OMP_NUM_THREADS=1` on fast entries; `NJOBS` at half the cores; explicit teams on the library suites and `RUN_SERIAL` on the workflow gates. |
| Removing many executables breaks scripts and historical validation packages. | Medium | Search all callers, update active ones atomically, document the command translation; historical evidence is preserved, not pretended runnable. |

## 15. Validation plan

Compilation and runtime validation are user-run unless separately authorised.

### 15.1 Static

1. Every inventory row carries a verdict with reviewer and date, and resolves
   to exactly one final callable implementation or a retired-tests row.
2. Every test UI entry has one reachable commander/module path.
3. Every CTest registration refers to a valid suite or `test=<id>`, and every
   registration has a label, a timeout and a working directory.
4. Every fast registration is assertion-bearing (grep for `assert_`,
   `tests_failed` or an equivalent typed check in its module).
5. No active invocation of a removed `simple_test_*` executable remains.
6. Test implementation modules contain no program units and no `stop`.
7. No migrated test has a new dedicated commander type whose only purpose is
   to call one procedure.
8. `BUILD_TESTS=OFF` source filtering excludes all test-only modules.
9. `git diff --check` and non-compiling syntax diagnostics pass for edited
   files.

### 15.2 Per-test behavioural

For every migrated test compare before and after: process exit status;
assertion and failure counts; expected normal-stop or failure marker; required
inputs and defaults; output files and retained artifacts; numerical values at
the pre-existing tolerance; cleanup and working-directory behaviour; launcher,
process count, threads and device requirements; skip behaviour when a
capability or fixture is absent; and wall time against the inventory.

### 15.3 Integration

1. `compile_debug.sh --compile-tests` and `compile_clean.sh --compile-tests`
   run the fast gate green in under 30 s on the reference Mac, in Debug and
   Release, and `ctest_budget.py` reports no violation.
2. `simple_test_exec test=list` contains every suite and canonical test ID once.
3. Each area suite passes in one process and each of its tests passes alone.
4. A deliberately failing assertion in one test fails its suite, fails the
   gate, and is reported by name.
5. Mother suites survive an intentionally failing child and return a failing
   aggregate status.
6. Coarray and MPI invocations use the required launcher and process count.
7. The nightly `library|workflow` run completes on the dedicated machine and
   its archive holds per-suite times and per-workflow metrics against
   declared floors.
8. A `BUILD_TESTS=OFF` build contains no test executable or test-only object;
   a `--compile-tests` install contains `simple_test_exec` and no standalone
   `simple_test_*` binaries.

## 16. Acceptance criteria

The project is complete when:

1. `compile_*.sh --compile-tests` builds and then runs a fast gate that
   finishes under 30 s on the reference Mac, in which every entry can fail.
2. The extensive gate runs nightly unattended: the library suites, one
   process per coherent part of the library, and the simulated workflows
   reporting resolution and pose-recovery metrics against declared floors.
3. `simple_test_exec` is the sole public Fortran test executable, and all
   supported former standalone programs are callable procedures in grouped
   modules run by area suites.
4. Every test has one authoritative implementation and a recorded review
   verdict; every removed test has a retired-tests row.
5. CTest, CI, documentation and child suites use suites or `test=<id>`.
6. Process isolation and special launchers are preserved where section 5.3
   requires them, and nowhere else.
7. The `production/tests/simple_test_*.f90` glob is gone.
8. Test-only code is absent when `BUILD_TESTS=OFF`.
9. The process budget and the time budget are enforced on every build.
10. Compilation and runtime checks not actually observed are listed as
    outstanding rather than claimed as passing.

## 17. Non-goals

- An X-style `validate` tier: a registry of real datasets and arms, blessed
  platform-keyed baselines, a report directory and gated blessing. Too costly
  for SIMPLE; the extensive tier gates on simulation truth instead.
- Adding test-specific commander types.
- Rewriting numerical algorithms while moving their tests.
- Forcing all tests to share one procedure signature when their dependencies
  differ.
- Running the extensive tier, or any subprocess-isolated case, inside the
  fast gate.
- Converting NICE/Python tests or external oracle analyzers into Fortran.
- Making network-dependent tests part of any registered tier.
- Replacing CTest, the SIMPLE command-line system or `simple_test_utils` as a
  prerequisite for consolidation.
- A generated universal test manifest before executable unification shows it
  is needed.

## 18. Effort and delivery shape

This is a large project, because two of the four workstreams are not
mechanical:

- Phase 0 and workstream A (the timing inventory and the fast-gate
  scaffolding over today's qualifying tests): days, and immediately useful.
- Workstream B (review of 149 identities, then unification of what
  survives; 53 two-route identities and 14 footprint clusters to reconcile): three
  to six weeks of area-by-area work, dominated by the review and by preserving
  the meaning and launch behaviour of what is kept. The review itself is a
  few weeks of owner and reviewer time, an afternoon per area, spread across
  the migration.
- Workstream C (splitting `units` into area suites, reconciling its routes,
  timing and trimming to 30 s, settling the socket and child-process
  sub-suites): about a week, and it delivers the whole fast gate. Assertions
  for print-only tests are needed only for the ones the review keeps for the
  extensive tier, and are written as those tests are migrated in Phase 4.
- Workstream D (the library suites as their areas are reviewed, the
  simulation-truth gates, the nightly runner): a few days per library suite
  once its review is done, one to two weeks for the first three workflow
  gates; more as workflows are added.

The order of value is A, then C (the fast gate is complete after it), then
the review and migration of the rest area by area as those areas are touched
in ordinary development, so the migration rides on work that is happening
anyway.
