# Nightly suite triage: a standing plan for an autonomous run

A plan that is run again whenever the nightly suite needs attention. It runs the
suite, reports the status of every entry, repairs what has an evident repair,
speeds up what has an evident speed-up, fuses entries that repeat each other's
work, brings the nightly test code to the style guides, describes with a
proposal every failure that it cannot repair, and proposes tests for the code
that changes most and is tested least. Written 2026-11-01 (base
`eefc7e1a1`) after `ci_test.yml` had failed every night since 2026-10-25.

The plan stays in `planned/`. Each run writes its report beside it
(`nightly_suite_triage_report.md`, replacing the previous one; git keeps the
earlier reports) and leaves its Progress section in place; the next run
clears Progress when it starts. For `arun`: `PLAN=PLAN_DONE=` this file,
`ALLOW_TABLE="### File table"`, `PHASES="0 1 2 3 4 5 6"`.

## The suite

The nightly is `.github/workflows/ci_test.yml`, on `ubuntu-26.04` and
`macos-latest`: `./compile_clean.sh` (Release, the fast gate), then
`ctest -L platform`, `ctest -L library` and `ctest -L highlevel`, in one step.
The entries and their contracts are in `doc/policies/test_environment_policy.md`
(sections 1.2 and 1.5) and `src/main/commanders/test/simple_test_registry.f90`;
`production/CMakeLists.txt` registers them. Ruben owns the runner's closeout
(`phase5_workflow_gates_and_nightly_runner_handover.md`); Hans reviews floors
and CTest entries.

GitHub's standard runners are small (Linux 4 cores and 16 GB, macOS 3 cores and
7 GB) and end a job at 6 hours. A local run on a large machine finds what is
broken; only the CI logs show what fails on those runners and on macOS.

## Standing decisions (Hans, 2026-11-01)

- The run reads the CI logs itself with `gh` (logged in on the machine; read
  access to `hael/SIMPLE` is enough).
- A production fix that makes a failing entry pass is applied in the working
  tree, verified, and listed apart from the test fixes in the report.
- `.github/workflows/ci_test.yml` may be changed in the working tree.
- Entries that should be fused are fused, and nightly test code that does not
  follow the style guides is brought to them.
- Hot areas (code that changes a lot) with poor coverage get proposed tests in
  the report; Hans decides which are written. The run writes none of them.

## Terms

- **Entry**: a CTest registration (`ctest -N` lists them). **Nightly entries**:
  the `platform`, `library` and `highlevel` labels.
- **Evident repair**: the cause is shown (a log line, a diff, a commit) and the
  repair restores what the test was written to check, without changing what it
  asserts.
- **Triage class** of a failure, one of:
  1. *stale test*: a renamed program or key, a moved or renamed output, a
     registry or CMake row out of step;
  2. *environment*: runner memory, cores or deadline, platform, the workflow;
  3. *deliberate divergence*: production changed on purpose and the test does
     not know;
  4. *regression*: production is broken;
  5. *unstable*: the entry passes and fails on the same code (seeds, threads).

## Rules

- Never lower a floor, widen a tolerance, weaken or remove an assertion, skip a
  check, or raise a deadline to make an entry pass. Such a change, when the
  evidence supports it, is a proposal in the report with that evidence. The
  test policy treats a widened tolerance as a finding.
- A failure at seed offset 0 is a failure. Runs at other seeds or thread counts
  (`test=profile`) characterise it and never overturn it.
- Shrinking a test (fewer particles, iterations or seeds, a smaller box)
  changes what it can detect: a proposal with the measured saving, not a
  change.
- A fusion or style refactoring keeps every check and the result: one-thread
  entries are compared with `test=identity`, the others by their gate rows at
  seed offset 0 against the baseline (equal, or inside the entry's profile).
- No new CTest entries. A fusion lowers `SIMPLE_CTEST_BUDGET` and updates, in
  one change, the registry row, the CTest registration, the test UI and router,
  `scripts/check_test_registry.py` if needed, and the policy's tables.
- A repair adds the files it changes to the file table, with the reason, before
  editing them.
- Git is read-only; nothing is committed. Older commits are built from
  `git archive` exports, never by checking them out.
- `gh` is read-only too: `gh run list`, `gh run view` and GET requests of
  `gh api`. The run never pushes, dispatches, re-runs, comments or changes a
  setting on GitHub, whatever the token would allow.
- Fortran: the `simple-modern-fortran` skill, the compile-time policy and the
  test policy (sections 2 and 5). New or changed OpenMP code gets a tester at
  24 threads compared bit for bit with one thread. THROW macros are never
  continued with `//&`; no `.and.` short-circuit on allocatable arrays.
- Run outputs go to a scratch directory on a disk with at least 60 GB free (the
  high-level entries alone grow to about 37 GB). `build/test_runs`, where CTest
  works, may be a link to it. Bulky outputs are deleted once their numbers are
  in the journal; reports, `metrics.tsv` and `run.log` files are kept.
- Nothing else heavy runs on the machine while timings are measured.
- Stop and report when the build breaks and the cause is not this run's, or a
  disk falls under 20 GB free.

## Phases

### Phase 0: evidence and baseline

1. Record the commit, compiler, CMake cache (build type, flags), machine and
   disks.
2. CI history: `gh run list -R hael/SIMPLE -w ci_test.yml -L 14`, and for each
   run `gh run view <id> -R hael/SIMPLE --log` (the job logs of both
   platforms). From the CTest summaries: per night and platform, each entry's
   verdict and time; *not run* for a tier a failed tier stopped, and
   *cancelled* for a job the matrix cancelled; the build or gate failure when
   the step before the tests failed.
3. Build (`make -C build -j`), the fast gate, install.
4. The baseline, one tier after the other on the otherwise idle machine:
   `ctest -L platform` and `ctest -L library` (with `--output-junit`), then the
   high-level entries through `simple_test_exec test=nightly ncunits=<cores>
   dir=<scratch>/nightly` (the same command lines and environment as
   `ctest -L highlevel`, plus each gate's value and floor, wall time, peak CPU
   and memory in `nightly_report.tsv`).

Exit: the status table, one row per nightly entry: local verdict, local wall
time and peak memory, CI verdict and time per platform over the last nights,
and for every failure its triage class with the evidence.

### Phase 1: repairs

For each failing entry, in the order of the status table (failures on every
platform first):

1. Read the entry's `run.log`, report and `metrics.tsv`, the test, and the
   history of the test and of the production code it drives since the entry
   last passed (`git log -p` on both); when the history does not decide it,
   build the older commit from a `git archive` export and run the entry there.
2. Classes 1 and 2: repair. Class 4: repair when the defect and its repair are
   evident, otherwise describe it. Class 3: describe what production does now,
   what the test expects, and the proposal (adapt the test, or restore the
   behaviour); adapt the test only when the commit that changed production
   says the change is intended and the test's guarantee survives it. Class 5:
   `test=profile target=<entry> nseeds=5` (more when the machine is free),
   report the failure rate, and propose (a profile entry, a stronger signal,
   or a fix of the instability's cause).
3. After each repair: build, the fast gate, the entry rerun alone in its own
   scratch directory. A production repair also reruns the fast suite of its
   area and every nightly entry that drives the changed code.

Exit: every failure is repaired and verified, or described with a proposal.

### Phase 2: speed

1. Rank the entries by core-hours (wall time times peak CPU) from the
   baseline, and estimate each one's time on the CI runners from the CI logs
   where they hold it.
2. Profile the slow entries' stages (their `BENCH` files, logs, `perf` where it
   helps). Evident speed-ups are made: the same fixture built several times,
   threads reserved but idle, work whose result nothing reads. Each is measured
   before and after and keeps the entry's results (as for refactorings).
3. Anything else that would shorten the suite is a proposal with its saving.

Exit: every change measured, every proposal priced.

### Phase 3: fusion and style

1. Fusion candidates: entries or sub-suites that build the same fixture or run
   the same pipeline prefix, or one entry that splits what one process could
   do without losing isolation that matters (a crash in one part hiding the
   others, memory, CI packing). A candidate is fused when it removes repeated
   work or setup, keeps every check, keeps each part runnable alone
   (`suite=`), and does not lengthen the CI wall time. Otherwise it is listed
   with the reason it stays.
2. Style: the nightly test code (the high-level test commanders, the library
   suites' testers, the registry and the runner) against the
   `simple-modern-fortran` skill and the test policy: `use ..., only:`, large
   objects `class(T), allocatable` in testers, `new`/`kill` symmetry,
   assertions through `simple_test_utils`, tolerances named and derived, draws
   from `set_fixed_seed`, gate rows through `simple_test_gate`
   (`metrics.tsv`), one `!@descr:` per file. Changes are verified as
   refactorings (Rules).

Exit: fusions done and verified, with the budget, registry, CMake, UI and
policy in step; style changes verified; what was left, with the reason.

### Phase 4: the CI workflow

Bring `ci_test.yml` to what the phases found:
- both platforms always finish (`fail-fast: false`);
- every tier runs whatever the tier before it did;
- the test reports, `metrics.tsv` files and CTest logs are uploaded as
  artifacts, also on failure;
- explicit job timeouts;
- when the baseline shows it pays, the tiers or groups of entries run as
  parallel jobs, each under its own 6-hour limit;
- entries the macOS runner cannot hold (registry `mem_gb` against 7 GB) are
  handled explicitly, never silently skipped.

The YAML is checked locally (parsed, and with `actionlint` if present). A
workflow is only proven by its first run on GitHub: the report says so.

Exit: the workflow changed and checked, or the reasons it was left.

### Phase 5: hot areas and their coverage

1. Churn: over the last eight weeks (`git log --since='8 weeks ago'
   --no-merges --numstat`, renames followed), per file and per subsystem under
   `src/`, testers excluded: commits, lines changed, and commits whose message
   names a fix. Bulk renames and moves are not churn. The hot areas are the top
   twenty files, and the subsystems they belong to.
2. Coverage: a scratch build directory configured with the coverage flags on
   the CMake command line (`-DCMAKE_Fortran_FLAGS=--coverage`, the same for C
   and C++, Release otherwise; nothing in the repository's CMake changes); the
   fast gate, `ctest -L library` and `test=nightly` run from it; `gcov
   --json-format` per object, summarised per file and per tier (fast alone,
   and all tiers). The run's timings are not baseline timings.
3. For each hot area: its line coverage per tier, which entries reach it, and
   whether a test asserts its behaviour or a workflow only executes it. A hot
   area that only the high-level workflows execute, or that nothing executes,
   is poorly covered.
4. For each poorly covered hot area, a proposed test: the guarantee it pins,
   its tier and suite (a sub-suite of an existing area suite before a new
   entry), how its expected values are derived independently of the code (test
   policy, section 2.1), its negative control, and its cost in run time.

Exit: the churn and coverage table of the hot areas, and the proposed tests.

### Phase 6: verification and report

1. Build from a clean configure in a scratch build directory, the fast gate,
   then the full nightly again (`platform`, `library`, `test=nightly`).
2. The report `nightly_suite_triage_report.md`: the status table before and
   after; each repair with its class, cause and evidence, test repairs and
   production repairs apart; each proposal; speed results and proposals;
   fusions and style changes; the CI changes; the hot areas, their coverage
   and the proposed tests; what is not verified (macOS, the CI runners) and
   how to verify it.
3. Add the run as evidence to `phase5_workflow_gates_and_nightly_runner_handover.md`.
4. End with `git status`.

Exit: the second nightly passes, or every failure in it is in the report with
its class and proposal.

### File table

| file | why | phases |
| --- | --- | --- |
| `doc/refactoring_notes/planned/nightly_suite_triage.md` | this plan and its progress | each |
| `doc/refactoring_notes/planned/nightly_suite_triage_report.md` | the run's report | 6 |
| `doc/refactoring_notes/planned/phase5_workflow_gates_and_nightly_runner_handover.md` | the run as nightly evidence | 6 |
| `doc/policies/test_environment_policy.md` | entries, fusions, budget, timings | 1 2 3 |
| `src/main/commanders/test/` | high-level tests, fixtures, registry, runner | 1 2 3 |
| `src/main/exec/simple_test_exec_highlevel.f90`, `src/main/apis/simple_test_exec_api.f90` | routers of the test programs | 1 3 |
| `src/main/ui/simple_test/` | test UI programs and `suite=` help | 1 3 |
| `production/CMakeLists.txt` | CTest registrations, deadlines, budget | 1 2 3 |
| `scripts/check_test_registry.py`, `scripts/ctest_budget.py`, `scripts/run_fast_gate.sh` | the registry lint and the gate | 1 3 |
| `scripts/simple_args_generator.pl` | keys of the test front door | 1 3 |
| `.github/workflows/ci_test.yml` | the nightly workflow | 4 |

A repair, speed-up or fusion that touches other files (a tester beside the code
it tests, production code) adds them here with the reason first.

Added by the run of 2026-11-01:

| file | why | phases |
| --- | --- | --- |
| `src/main/**/*_tester.f90` reached by the `lib_<area>` suites (flex, nano, volume, stream, strategies, solve, commanders/simple) | Phase 3 style review: imports, assertions that could not fail, lifecycle, seeds that did not move with `seed_offset` | 3 |
| `.github/workflows/ci_push.yml` | a build and fast gate on every push to `master` (ruling of 2026-11-02) | 4 |
| `src/main/volume/simple_reconstructor_pcg.f90`, `src/main/volume/simple_wave_projector.f90` | production repair: measured FFTW plans only without a fixed seed (one-thread runs were not reproducible) | 3 |

## Progress

Run of 2026-11-01 on the Dell (interactive session, not arun), base
`eefc7e1a1`; scratch `/data2/agent_runs/nightly_2026-11-01/` (journal
`progress.md`, `build/test_runs` linked to its `ctest_runs/`).

### Phase 0 (started 2026-11-01)

- Build, fast gate (15 of 15, 8 s) and install at `eefc7e1a1`.
- CI: every `ci_test.yml` run from 2026-10-25 to 2026-11-01 failed: CTest exit 8
  in the test step (Ubuntu), the build step's exit 2 on 10-28 (macOS) and 10-31
  (Ubuntu); macOS was cancelled by the matrix whenever Ubuntu failed first; the
  test step took 3 h 20 min to 5 h 10 min.
- Baseline started 18:45 (platform, library, then `test=nightly ncunits=24`).
- CI history from the logs of all 45 runs since 2026-09-29 (`gh`; scratch
  `ci_entries.tsv`, `ci_times.tsv`). Last green night 2026-10-01. Of the
  current failures:
  - `multistate_assignment`, both platforms every night since 10-26, 4 s in:
    `refine=prob|prob_neigh needs objfun=euclid` (`simple_eul_prob_tab.f90`),
    added by `811934ca5` (10-24); the test runs `objfun=cc`. Class 3.
  - `species_discovery`, Ubuntu since 10-29: `density_gradient: largest error
    of the measured transfer = 0.0510`, floor 0.050. Suspect `62ffd5924`
    (10-28, nano density table bounds). Class 3 or 4, open.
  - `multistate_workflows` timed out (7200 s) on about half the Ubuntu nights;
    when it passed, median 5 498 s, maximum 6 312 s. Class 2 (speed).
  - Build breaks before any test: 10-28 (macOS) `simple_nano_density.f90`
    missing from `8e0ea79a9` (added by `bae73b8ec`, "forgotten files");
    10-31 (Ubuntu) link error, `shiftmap_peak`/`shiftmap_objfun_grid`
    undefined at `653c5cd21`; 10-28 (Ubuntu) the bundled OpenBLAS failed to
    build (`99175130e`, `f24d336b7` era). All gone at HEAD.
  - Earlier failures (`simulated_workflow_6vxx`, `single_workflow_fcc`,
    `solve3D_addon`, `flex_pca_blobs`, `lib_stream`) pass on every night since
    10-25.
  - Ubuntu job time: the CI medians sum to 4.3 h of the 6 h limit;
    `multistate_workflows` 1.5 h, `lib_stream` 0.7 h.
  - macOS has not finished a test step since 10-01: cancelled by the matrix
    whenever Ubuntu failed first.
- Baseline done (54 min): platform 1/1, library 6/7 (`lib_cart_align3D`, local
  only), high-level 18/19 (`multistate_assignment`); `species_discovery` passes
  here at 0.0486.

### Phase 1 (done 2026-11-01)

Every failure classified with evidence (report, section 3):
`multistate_assignment` class 3, `species_discovery` and `lib_cart_align3D`
class 5, `multistate_workflows` class 2 (speed), the build breaks class 2. No
evident repair among them; one repair of the environment class (the fast
gate's warning about untracked sources).

### Phase 2 (done 2026-11-01)

`multistate_workflows` timed per route on 4 pinned cores (81 min, no
test-side waste); no evident test-side speed-up elsewhere; the CI layout
carries the gain (Phase 4).

### Phase 3 (done 2026-11-01)

Two read-only reviews (high-level commanders; library testers); the
conformance items that change no result made, the rest listed as proposals;
no fusion pays. Verification found `single_workflow_wurtzite` not reproducible
at one thread (measured FFTW plans under a fixed seed, a production repair) and
two stale identity rules (repaired).

### Phase 4 (done 2026-11-01)

`ci_test.yml` rewritten: four parallel groups per platform, `fail-fast` off,
reports uploaded, `ctest --parallel` by `PROCESSORS`.

### Phase 5 (done 2026-11-02)

Churn of eight weeks; coverage of the fast tier with every file instrumented,
of the library and high-level tiers with the 22 hot files only (with every
file instrumented the threaded entries ran 13-39 times slower). Nine proposed
tests (report, section 9).

### Phase 6 (done 2026-11-02)

Final run on the working tree: fast 15/15, platform and library 7/8, high-level
18/19, the same two failures as the baseline; the one-thread workflows
byte-identical (section 7 of the report). Report written; the Phase 5 handover
links it.

Rulings of 2026-11-02 on the report's decisions, carried out the same day:
`multistate_assignment` retired, the `species_discovery` floors re-derived, the
`lib_cart_align3D` perturbation moved inside the rotation bound,
`multistate_workflows` 4 h on Linux only, `ci_push.yml` added; then the stale
checks removed, the quick tests' seeds fixed (and `mini_stream`'s defocus
ceiling re-derived), `class(T), allocatable`, the fixture's `new`/`kill`, merged
helpers, rows before stops, the floors' origins, and two new tests (PCG density
support, volume programs). `classify3D_refs` and the probability-table fill stay
proposals: a workflow and a library-tier fixture, not quick tests.

Lessons for the next run: a fully instrumented coverage build is unusable for
the threaded tiers (instrument the target files only); journal times from
`date`, never estimated; `test=identity` between the baseline and the final run
is the cheapest proof that test edits changed no result, and it found the one
production defect of this run. Price a proposed test from the code it must run (a whole
workflow, a fixture with noise spectra) before offering it as a quick one.
