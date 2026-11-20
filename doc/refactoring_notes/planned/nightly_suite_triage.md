# Nightly suite triage: a standing plan for an autonomous run

Status, 2026-11-19: standing plan, revised for its second run, whose target is
real-data `solve3D` (decision 13). Written 2026-11-01 (base `eefc7e1a1`); first
run 2026-11-01/02, ruled and committed as `4f3fed7e5`.

A plan that is run again whenever the nightly suite needs attention. It runs the
suite and reports the status of every entry, looks for regressions in entries
that still pass, repairs what has an evident repair, speeds up what has an
evident speed-up, fuses entries that repeat each other's work, brings the
nightly test code to the style guides, describes with a proposal every failure
it cannot repair, and proposes tests for the code that changes most and is
tested least.

The plan stays in `planned/`. Before a run, the launcher fills section 1; the
standing body (sections 2 to 10) changes only by a ruling. The run writes its
report beside the plan (`nightly_suite_triage_report.md`; git keeps the earlier
ones) and its Progress at the end; when it starts, it moves the previous run's
Progress into one line of "Previous runs". For `arun`: `PLAN=` this file
(`PLAN_DONE` unset: the plan stays in `planned/`), `ALLOW_TABLE="## 10. File
table"`, `PHASES=` the phases of section 1.

## 1. This run

| item | value |
| --- | --- |
| Base | the run's HEAD (`state/base_head`; `7e88a551e` when this was written). The run works on its own copy of the dell's checkout in `<run root>/repo`; this plan travels as a launch difference unless it is committed first (decision 15) |
| Reference | `603c44490` (2026-11-15), a default (section 3): the commit at which Hans read the three real sets as perfect (Known 7) |
| Machine | dell2: Ubuntu 26.04 (the CI runner's release), 16 cores (4 P, 8 E, 4 low-power E), 61 GB, idle |
| Run root | `/data2/agent_runs/nightly_triage` behind `~/agent_runs/nightly_triage` on dell2; the scratch is its `scratch/` |
| Phases | 0 1 2 3 5 7 (4 and 6 left out: both ran on 2026-11-01/02) |
| Target | real-data `solve3D`: about 1 run in 10 fails subtly (Hans, 2026-11-19) |
| Real-data cases | the three real sets of the class and pose student runs (decision 18), below |
| Real-data stacks | on dell2 under `/data2/dell_data/`, in the layout of the dell's `/data`, with the dell's checksums in `/data2/dell_data/MD5SUMS` (copied 2026-11-19, Known 8). Phase 0 remaps frozen copies of the input projects to them |
| Phase 1 budget | 72 h of wall time on dell2 (guessed; Phase 0 prices it again from one run per set) |

Real-data cases (decision 18). Each run is a plain `solve3D` on a frozen copy of
the set's selection project, single state, `nthr=16`, with `SIMPLE_SEED=<seed>`
(seeds 1 to 10, the same on both builds), each run in its own directory. The
keys are the student runs' (`pose_student_physical_raw_picks.md`, the teachers'
`tools/teacher.sh`):

| set | input project | keys | particles selected | final box, pixel size |
| --- | --- | --- | --- | --- |
| bgal | `/data2/bgal_addon/3_selection/bgal_addon.simple` | `nsample=1500 mskdiam=180 pgrp=d2 pgrp_start=c1` | 2 756 | read in Phase 0 |
| trpv1 | `/data2/trpv1_c2/3_selection/trpv1_c2.simple` | `nsample=5000 mskdiam=175 pgrp=c2 pgrp_start=c2` | 35 504 | read in Phase 0 |
| strep | `/data2/dell_data/10335-streptav/tst/3_selection/tst.simple` | `nsample=5000 mskdiam=85 pgrp=d2 pgrp_start=c1` | 6 821 | read in Phase 0 |

The failure criterion, frozen here before any run (the method of the
small-membrane note, Known 7):

- A run *fails* when it does not complete, or when its final FSC=0.143
  resolution (the `POSTPROCESS` low-pass value) is coarser than its set's
  ceiling. The ceiling is the median of the reference's ten runs on the set,
  plus the larger of `3s` and one Fourier shell at the final box and pixel size
  (`s` is 1.4826 times the median absolute deviation of the ten). It is
  computed when the reference's ten runs of the set are done, and frozen before
  any base run of the set is read. The worst reference run is not the anchor,
  because the reference may fail one run in ten itself. The origin is measured
  (the reference runs), the rule guessed.
- The consolidation readout, reported for every run and not a gate: overlap
  at least 0.30 and angular change at most 30° at the end of stage 3 or 4, read
  from the iteration log as the note did. A run that passes its ceiling without
  consolidating, or consolidates and then fails, is flagged.
- Reported, not gated: the final map, the halves and the log of every failed
  or flagged run, and of seed 1 on both builds, kept under
  `<scratch>/keep/real/` with their paths in the report, for Hans to read. He
  judges a "subtle" failure the ceiling misses.

Known at writing (2026-11-19, against CI and the tree at `7e88a551e`; Phase 0
checks each again):

1. `ci_test.yml` did not run from 2026-11-09 to 11-16. The last night every
   job passed is 11-08, at `913d50765`: 91 commits before `7e88a551e`, 64 of
   them in `solve/`, `strategies/`, `sigma2/`, `volume/` or `nu_filt/`.
2. 11-17, 11-18 and 11-19 are red on Ubuntu only, all nine failed jobs (3, 4
   and 2 of the 4 Ubuntu groups) for one reason: in the build step's fast gate,
   `unit_numerics` fails HHCNN's `test_thread_identity` ("the chunked reduction
   gives the same weights at one thread and at a team of 4, bit for bit",
   `src/main/learn/simple_hhcnn_tester.f90`). On 11-17 and 11-19 the other
   Ubuntu jobs of the same commit passed it, macOS passed every time, and
   `ec990429f` ("fixed flaky hhcnn test", 11-16) is in all three commits. A failed build step runs no
   test: Ubuntu has not run the `workflows` group since 11-08.
3. Green is not right. On 11-08 `simulated_workflow_6vxx`'s map was wrong on
   both platforms (truth correlation 0.52, below the wrong-truth control's
   0.60; FSC 31 and 50 Å) and passed because its map rows were reports then.
   `3c3e0ee8e` (11-12) raised the fixture's SNR and restored the hard map and
   pose gates (`MIN_VOL_CORR` 0.80, FSC 6 Å, `MAX_PAIR_POSE_ERR` 15°,
   `MIN_PAIR_FRAC5` 0.45). The floors carry no profile; the test card asks for
   a fresh 20-seed profile (the last, 2026-10-19: 2 of 20 missed the map).
4. `simulated_workflow_1jxy` was removed on 11-17 (`fd9af0add`): its `solve3D`
   failed the pose gate in 13 of 40 solves, the pair fraction within 5° at
   0.48-0.50 while the map passed. That is a subtle failure of the kind the
   target names. No registered entry runs `solve3D` over parts since.
5. CI keeps the test reports 14 days. The 11-08 reports are saved on dell2 in
   `/data2/agent_runs/ci_artifacts/42449321961_913d50765_2026-11-08/` (37
   `metrics.tsv`); those of 11-17 to 11-19 expire on 12-01 to 12-03.
6. A run records its 64-bit seed in the project's info segment
   (`adopt_run_seed`): `rnd_seed_w1` to `rnd_seed_w4`, 16 bits each, low word
   first, and `rnd_version` (`print_project_field oritype=projinfo`).
   `SIMPLE_SEED=<seed>` replays it, except draws made before the seed was set
   in a run started without `SIMPLE_SEED`. `restarted_exec` (`nrestarts=`)
   gives every restart the same environment and does not vary the seed, so
   each triage run is its own process with its own `SIMPLE_SEED`. Whether a
   real-data run repeats at 16 threads is not known; Phase 0 measures it.
7. The [small-membrane note](../../implementation_notes/completed/solve3D_small_membrane_regression.md)
   (archived 2026-11-17) found the same kind of failure on PorF at
   `603c44490`: one restart in five or six at 7.43 Å against 4.11-4.22 Å,
   separated from the successes by a stage-3/4 consolidation readout. It read
   streptavidin, beta-galactosidase and TRPV1 as perfect at `603c44490`. The
   `solve3D` sampling changes since then came from it: 31 first-parent commits
   from `603c44490` to `7e88a551e`, among them `360fbeb4f`, `8363af240`,
   `e94681175` and `11f832a8e`. The 3D in-plane posterior insertion, off by
   default at `603c44490`, was on by default from `1ea0f8340` until
   `fd9af0add` removed it. `nsample_ramp` is not part of the library: added
   by `76110d3f5` as a switch off by default, which the real-data cases never
   set, it was withdrawn by `6668fb842`. On 2026-11-19 it had no trace in
   `src/`, the skills, the generated indexes or the build; only the solve3D
   policy names it, as withdrawn. The note's seeded restarts and per-stage
   summary line (its decisions 6 and 7) were never implemented.
8. dell2 has the three sets' project folders (`/data2/bgal_addon`,
   `/data2/trpv1_c2`), but not their particle stacks. The projects point at the
   dell's `/data/bgal/DAT/sumstack.mrc` (1.46 GB), `/data/trpv1/flipnew.mrc`
   (9.37 GB) and `/data/10335-streptav/` (2.5 GB, with streptavidin's project).
   dell2 has no `/data`. On 2026-11-19 they were copied to dell2 under
   `/data2/dell_data/`, keeping the layout below `/data`. The projects still name
   `/data`, so Phase 0 writes remapped copies with `update_project
   old_root=/data new_root=/data2/dell_data`. The remap matches whole path
   components, so the `/data2/...` class-average paths stay as they are.

## 2. The suite

The nightly is `.github/workflows/ci_test.yml`: Ubuntu and macOS jobs in four
groups (library and platform, workflows, rest, multi-state workflows; the last
Linux only, decision 8). Each builds Release with the fast gate, then runs its
CTest selection with explicit timeouts and report uploads. Re-read the workflow
and the registry at the start of every run; group membership and the process
budget change. The contracts are in `doc/policies/test_environment_policy.md`,
the registry rows in `src/main/commanders/test/simple_test_registry.f90`, the
CTest registrations in `production/CMakeLists.txt`.

The entries that run `solve3D`, for decision 13: `simulated_workflow_6vxx`
(single state, shared memory, one thread, deterministic; 1 129 s, reserving 2
cores), `solve3D_addon` (eight threads, 289 s, reserving 11) and
`multistate_workflows` (eight threads, 1 913 s for all routes, reserving 11;
`suite=` runs one route). The times are the test cards' and skills'; Phase 0
measures them on the run's machine.

Ruben owns runner closeout in
[the consolidated testing plan](../../implementation_notes/planned/testing_and_benchmarking_plan.md#6-consolidated-test-and-ci-work);
Hans reviews floors and CTest entries.

GitHub's standard runners are small (Linux 4 cores and 16 GB, macOS 3 cores and
7 GB) and end a job at 6 hours. A local run finds what is broken; only the CI
logs show what fails on those runners and on macOS.

## 3. Decisions and rulings

Standing (Hans, 2026-11-01):

1. The run reads the CI logs and test reports itself with `gh` (logged in on
   the machine; read access to `hael/SIMPLE` is enough).
2. A production fix that makes a failing entry pass is applied in the working
   tree, verified, and listed apart from the test fixes in the report.
3. `.github/workflows/ci_test.yml` and `ci_push.yml` may be changed in the
   working tree.
4. Entries that should be fused are fused, and nightly test code that does not
   follow the style guides is brought to them.
5. Hot areas (code that changes a lot) with poor coverage get proposed tests in
   the report; Hans decides which are written. The run writes none of them.

Rulings of 2026-11-02 on the first run's report (carried out in `4f3fed7e5`):

6. `multistate_assignment` is retired and not restored.
7. The `species_discovery` floors were re-derived by the floor rule; the
   `lib_cart_align3D` perturbation lies inside the rotation bound.
8. `multistate_workflows` has a 4 h deadline and runs on Linux only in CI.
9. `ci_push.yml` builds and runs the fast gate on every push to `master` that
   touches code.
10. `classify3D_refs` and the probability-table fill stay proposals: a workflow
    and a library-tier fixture, not quick tests.
11. The stream chain tests (`lib_stream`'s chain sub-suite, `stream_chains`)
    belong to Joe and Ruben: a failure there is described with its evidence,
    not repaired.

Rulings of 2026-11-19 (this revision):

12. Fast entries that fail in CI are in scope and triaged as nightly entries: a
    failed fast gate stops its CI group before any nightly entry runs.
13. Regressions in entries that still pass are in scope (Phase 1). The target
    of the 2026-11 run is real-data `solve3D`, which fails subtly in about 1
    run in 10.
14. Real data stays outside the suite (test policy, section 6.6): real-data
    runs are made by hand, from the cases of section 1, on the export builds
    Phase 1 names; reported, not gated.
15. The run stays on its base: `master` is not pulled during the run.
16. The real-data runs follow section 1: its cases, its seeds and its
    failure criterion, frozen before any run.
17. A comparison of one entry across commits holds its test fixed (Phase 1,
    step 2): the test commanders run production in process, so a test change
    between two commits is a change in what is measured.
18. (Hans) The real-data cases are the real sets of the class and pose
    student runs: beta-galactosidase, TRPV1 and streptavidin.

Defaults taken without a question (2026-11-19; vetoed by editing section 1
before launch):

- Reference `603c44490`, where Hans read the three real sets as perfect (Known
  7). The last green CI night, `913d50765`, is older and its own
  `simulated_workflow_6vxx` seed 0 missed the map (Known 3). Phase 1 measures
  the reference before it attributes anything to the range. If the real-data
  runs already failed at `603c44490`, Hans names an earlier commit.
- The real-data protocol of section 1: plain `solve3D` on each set's selection
  project with the student runs' keys, `nthr=16`, seeds 1 to 10 per build, and
  the failure ceiling of the reference's median plus the larger of `3s` and one
  shell.
- The real-data stacks live on dell2 under `/data2/dell_data/`, and the inputs
  are remapped copies (Known 8), so no `/data` is needed on dell2.
- Phases 4 (fusion, style) and 6 (coverage) are left out.
- Simulated entries: 20 seeds per arm, 100 for a lead on a target entry, 50 for
  any other lead (Phase 1, step 4).
- `simulated_workflow_1jxy` is profiled as a lead: at the reference, where it
  still exists, and at the base only as a scratch patch on an export that
  restores its suite (removed in `fd9af0add`), never in the working tree.
- The Phase 1 budget of section 1.

## 4. Not in this plan

- Writing the proposed tests (decision 5); new CTest entries; restoring
  `simulated_workflow_1jxy` in the working tree.
- The real-data dataset layer (test policy, section 6.6; testing plan,
  section 5) and real-data CTest entries.
- Stream chain repairs (decision 11).
- Any write on GitHub (section 7).
- Files outside the file table.

## 5. Terms

- **Entry**: a CTest registration (`ctest -N` lists them). **Nightly entries**:
  the `platform`, `library` and `highlevel` labels; fast entries that fail in
  CI count as nightly entries (decision 12).
- **Base**, **reference**: the two commits of section 1. **Export build**: a
  Release build of a `git archive` export in its own directory under the
  scratch, configured with the checkout's CMake options (`-DFORTLS=OFF`); never
  the checkout's `build/`.
- **Regression**: a gate metric of an entry, or the failure count of a
  real-data case, worse at the base than at the reference by the rule of
  Phase 1, step 4. **Lead**: a difference the screen or the rule points at but
  does not decide.
- **Evident repair**: the cause is shown (a log line, a diff, a commit) and the
  repair restores what the test was written to check, without changing what it
  asserts.
- **Triage class** of a failure or a regression, one of:
  1. *stale test*: a renamed program or key, a moved or renamed output, a
     registry or CMake row out of step;
  2. *environment*: runner memory, cores or deadline, platform, the workflow;
  3. *deliberate divergence*: production changed on purpose and the test does
     not know;
  4. *regression*: production is broken;
  5. *unstable*: the entry passes and fails on the same code (seeds, threads).

## 6. Change types and gates

| Phase | Kind | Gate |
| --- | --- | --- |
| 0, 1, 6 | measurement; no edit but this plan's section 1 and Progress | none |
| 2 | test repairs: tests only. Production repairs: identical or result-changing | identical: `test=identity` at one thread, base export against the repaired build. Result-changing: the rule of Phase 1, step 4: the repaired build no worse than the reference, and its repaired metric better than the base |
| 3, 4 | identical | `test=identity` at one thread; the other entries' gate rows at seed offset 0 equal, or inside the base profile where Phase 1 measured one |
| 5 | the CI workflow | parsed, and `actionlint` if present; proven only by its first run on GitHub |
| 7 | verification | the full nightly; the Phase 1 rule for every repaired regression |

## 7. Rules

- Never lower a floor, widen a tolerance, weaken or remove an assertion, skip a
  check, or raise a deadline to make an entry pass. Such a change, when the
  evidence supports it, is a proposal in the report with that evidence. The
  test policy treats a widened tolerance as a finding.
- A failure at seed offset 0 is a failure. Runs at other seeds or thread counts
  (`test=profile`) characterise it and never overturn it.
- Shrinking a test (fewer particles, iterations or seeds, a smaller box)
  changes what it can detect: a proposal with the measured saving, not a
  change.
- A fusion or style refactoring keeps every check and the result (section 6).
  `test=identity` between the baseline and the final run is the cheapest proof
  that test edits changed no result; in the first run it found the one
  production defect.
- No new CTest entries. A fusion lowers `SIMPLE_CTEST_BUDGET` and updates, in
  one change, the registry row, the CTest registration, the test UI and router,
  `scripts/check_test_registry.py` if needed, and the policy's tables.
- A repair adds the files it changes to the file table, with the reason, before
  editing them.
- Git is read-only; nothing is committed. Older commits are built from
  `git archive` exports, never by checking them out. A bisection is done by
  hand over `git rev-list`, since `git bisect` writes repository state.
- `gh` is read-only too: `gh run list`, `gh run view`, `gh run download` and
  GET requests of `gh api`. The run never pushes, dispatches, re-runs, comments
  or changes a setting on GitHub, whatever the token would allow.
- Fortran: the `simple-modern-fortran` skill, the compile-time policy and the
  test policy (sections 2 and 5). New or changed OpenMP code gets a tester that
  forces a team and compares it with one thread bit for bit, and a stress run
  by hand at the machine's full thread count (16 on dell2, 24 on the dell).
  THROW macros are never continued with `//&`; no `.and.` short-circuit on
  allocatable arrays.
- Run outputs go to the scratch, on a disk with at least 60 GB free (the
  high-level entries alone grow to about 37 GB). `build/test_runs`, where CTest
  works, may be a link to it. Bulky outputs are deleted once their numbers are
  in the journal; reports, `metrics.tsv`, `executions.tsv` and `run.log` files
  are kept. On dell2 both drives are QLC and slow down under long sustained
  writes; delete bulky outputs as soon as their numbers are read.
- Timings: nothing else heavy runs on the machine while timings are measured.
  On dell2 the four low-power cores make thread timings noisy: a timing records
  the cores it ran on, and a before/after pair runs pinned to the same cores
  (`taskset`). Timings are compared only within one machine.
- Journal times come from `date`, never estimated.

## 8. Stop conditions

Stop and report when:

- the base build or the reference export build breaks (an export at a bisection
  step that does not build is skipped and recorded, not a stop);
- the run's repository differs at launch from its HEAD by more than this plan
  and what the spec's rulings name (`state/launch_status.txt`);
- a real-data stack is missing on dell2, or its checksum differs from the
  dell's (Phase 0, step 5);
- a disk falls under 20 GB free;
- a repair would need a floor, tolerance, assertion or deadline to move.

Not a stop: Phase 1 reaching its budget (what is left is reported as leads).

Gates no ruling may waive: the fast gate after every edit; `test=identity` for
every change claimed identical; Phase 7's full nightly; the Phase 1 rule for
every result-changing production repair.

## 9. Phases

### Phase 0: evidence and baseline

1. Record the base (`state/base_head`) and the launch differences
   (`state/launch_status.txt`) in section 1, then the reference, compiler,
   CMake cache (build type, flags), machine, load and disks.
2. CI history of every `ci_test.yml` run since 2026-11-02: `gh run list -R
   hael/SIMPLE -w ci_test.yml`, for each run `gh run view <id> --log` and,
   while it exists, `gh run download <id>` into the scratch (the 11-08 reports
   are already saved, Known 5). Per night, platform and group: the build step's
   result and, for a failed fast gate, the failing check; each entry's verdict
   and time; each gate metric's value from `metrics.tsv`; *not run* for a group
   whose build step failed, *cancelled* for a job the matrix cancelled.
3. Build the run's repository in its `build/` (Release), the fast gate,
   install.
   Build the base and the reference again as export builds
   (`<scratch>/keep/build_base`, `<scratch>/keep/build_ref`), with each one's
   fast gate.
   They stay untouched until Phase 7.
4. The baseline at the base, one tier after the other on the otherwise idle
   machine: `ctest -L platform` and `ctest -L library` (with
   `--output-junit`), then `simple_test_exec test=nightly ncunits=<cores>
   dir=<scratch>/nightly_base` (the same command lines and environment as
   `ctest -L highlevel`, plus each gate's value and floor, wall time, peak CPU
   and memory in `nightly_report.tsv`). Then `test=nightly` of the reference
   export build (`dir=<scratch>/nightly_ref`).
5. Real data (section 1, decision 16):
   1. Check the stacks of Known 8 with `md5sum -c /data2/dell_data/MD5SUMS`.
   2. Freeze a remapped copy of each input project under
      `<scratch>/keep/real/inputs/`, with its checksum: `update_project
      projfile=<copy> old_root=/data new_root=/data2/dell_data
      projfile_out=<frozen>` with the base build. bgal and trpv1 are copied
      from `/data2/bgal_addon` and `/data2/trpv1_c2`, and their checksums are
      compared with the dell's in `MD5SUMS`. Every run copies its input from
      there.
   3. Read each set's final box and pixel size, and so the width of one
      Fourier shell at the reference's median resolution, once that is known.
   4. Check each set's command line against the UI's required keys on both
      builds.
   5. Run one base run per set at seed 1 and record its wall time and peak
      memory. If 10 seeds × 3 sets × 2 builds does not fit the budget beside
      the simulated profiles, trpv1 drops to 5 seeds per build first; the
      report labels it a lead.
   6. Run bgal at seed 1 a second time. If the two runs' final maps are
      identical, real-data runs are replayable at 16 threads and Phase 1's
      bisection by seed applies to them; otherwise real data gets rates only.

Exit: the status table, one row per nightly entry: local verdict at base and
reference, wall time and peak memory, the CI verdict and time per platform
since 2026-11-02, and for every failure its triage class with the evidence;
the entries removed since the last run, and what nothing gates since; both
export builds and both nightly reports; the real-data stacks on dell2 with
matching checksums, the frozen inputs, the run times, and the replay check.

### Phase 1: regressions in entries that pass

The question is whether a change between the reference and the base made an
entry, or a real-data case, worse while it still passes. One run cannot show a
rate of 1 in 10. So verdicts come from paired profiles, and attribution comes
from deterministic runs.

1. Screen every entry from Phase 0: per gate metric, its seed-0 value at
   reference and base, its CI series since 11-02, and the profile median and
   spread written beside its floor, where there is one. A metric worse at the
   base by more than its spread `s`, or, without a recorded spread, any change
   in a one-thread entry, is flagged. The screen only chooses what is
   profiled; its threshold is guessed.
2. Hold the test fixed. Before comparing one entry across commits, list the
   commits that change its test: `git log <ref>..<base>` on its test commander
   and fixtures. If there are any, the older export takes the base's versions
   of those files, as a patch. If they do not compile there, the comparison's
   older end moves forward to the first commit where they do, and the report
   says so. Since `603c44490` the workflow test changed only by `fd9af0add`
   (the 1jxy suite removed), so 6vxx comparisons from the default reference
   already hold the test fixed. An older reference crosses `3c3e0ee8e`, which
   changed the 6vxx fixture (Known 3).
3. Profile on the reference and base export builds, at the entry's CTest
   thread count, with `test=profile target=<entry> nseeds=20 keep_outputs=no
   dir=<scratch>/profiles/<build>_<entry>`. Compare with `test=profile dir1=
   dir2=`, and read the values per seed from each `executions.tsv`. Profiled
   whatever the screen says: the target entries `simulated_workflow_6vxx` and
   `solve3D_addon`; `multistate_workflows`, one route, when a real-data case is
   multi-state; and the lead `simulated_workflow_1jxy` (section 3, defaults).
   Then every flagged entry, while the budget lasts. The base profile of
   `simulated_workflow_6vxx` is also the fresh profile its test card asks for.
4. The verdict rule, frozen here before any profile is run:
   - Executions are paired by seed offset. A seed *fails* a gate metric as
     `test=profile` counts it: a failing row, no row, or a failed or timed-out
     execution.
   - Failure count: `b` seeds fail at the base and pass at the reference, `c`
     the reverse. A regression is a one-sided exact sign test on `b` against
     `c` at p ≤ 0.05 (for example 5 against 0, 7 against 1, 9 against 2).
   - Value: the per-seed difference in the metric's good direction, ties
     dropped. A regression is a median difference worse than one `s` (the
     larger robust spread of the two profiles) together with a one-sided sign
     test at p ≤ 0.05 (for example 15 of 20 seeds worse).
   - A lead (`b > c`, or a median worse by more than half an `s`, without the
     test) is profiled once more: 100 seeds per arm on a target entry, 50 on
     any other. The first 20 seeds of a one-thread entry must repeat their
     gate rows exactly, a check of the export builds. The same rule then decides; what it
     does not decide stays a lead, not a finding.
   - Power, stated before the run: a rise from 0 to 1 failure in 10 reaches
     `b` ≥ 5 at `c` = 0 with 20 seeds about 4 % of the time, with 50 about
     57 %, with 100 about 98 %. A rate like 1jxy's (13 in 40) shows at 20.
   - Control: a threaded entry gets a second base profile at the same seeds.
     The rule applied base against base must find nothing; otherwise that
     entry's verdicts are leads only. A one-thread entry's control is the
     repetition above.
   - The report lists every entry and metric tested beside the findings: at
     p ≤ 0.05 about one test in twenty finds a difference by chance.
5. Attribute each regression of a one-thread entry. The commits are `git
   rev-list --first-parent --reverse <ref>..<base>`.
   - For each of up to three seeds that fail at the base and pass at the
     reference, search by bisection for the first commit where the seed fails:
     one run per step, from an export build with the test held fixed. A commit
     that does not build is replaced by its neighbour and recorded. The range
     is not monotone where a behaviour was switched on and later removed (the
     in-plane insertion, `1ea0f8340` to `fd9af0add`, Known 7). A bisection that
     lands inside that span runs the seed at `fd9af0add` as well; when it
     passes there, the search continues from `fd9af0add` to the base.
   - When the seeds agree on a commit C, profile C and its parent with the
     regression's number of seeds; the rule of step 4 confirms C or not.
   - When the seeds point to different commits and none is confirmed, any
     result-changing commit can flip a seed's trajectory. The report gives the
     commits and the profiles as a lead.
   - A regression of a threaded entry is not bisected by seed. Its commit range
     and the `solve3D` commits in it are reported; a profile per bisection step
     needs a ruling.
6. Real data (section 1, decision 16). These runs come before the simulated
   profiles, one run at a time on the whole machine.
   1. The reference's runs on each set, seeds 1 to 10. Then each set's ceiling
      goes into Progress, before any base run of that set is read.
   2. The base's runs on each set, at the same seeds.
   3. The rule of step 4 on the failure counts, paired by seed, per set and
      pooled over the three sets. Power, stated before the run: with one
      failure in ten at the base and none at the reference, 10 pairs reach
      `b` ≥ 5 almost never, and 30 pooled pairs about 18 % of the time. The
      rates are therefore a lead unless the failures are many; replaying
      single seeds is what decides.
   4. Each seed that fails at the base and passes at the reference is
      bisected, when Phase 0 found real-data runs replayable. The search runs
      over `git rev-list --first-parent --reverse <ref>..<base>` as in step 5,
      one real-data run per step: 31 commits take five or six runs. Two or
      more seeds, of any set, that agree on a commit C attribute the failure to
      C; one seed alone is a lead. C is then checked on the simulated target
      entries' profiles at C and its parent (steps 3 and 4). When runs are not
      replayable, the report gives the base's failing seeds and the commit
      range.
   5. A seed that fails on both builds predates the range: it is reported, and
      nothing is attributed.
7. Classify each regression (class 3 or 4, from its commit's message and the
   policy notes of the area) and hand it to Phase 2 with its evidence.

Exit: for every target entry and every flagged metric, the reference and base
profiles and the verdict of step 4 (regression, lead or none) with `b`, `c`,
the median shift, the seeds and the control. For every regression, its
confirmed commit, its commit range, or its spread over several commits. For
every real-data set: its ceiling; a table of every seed on both builds
(completed, FSC=0.143 and 0.5, the consolidation readout, the flags); the
verdict; and, for each seed that fails only at the base, its commit or the
reason there is none. The test-fixed comparisons and their moved older ends.
The budget used.

### Phase 2: repairs

For each failure and each Phase 1 regression, in this order: failures on every
platform, fast-gate failures that stop CI groups, failures on one platform,
then regressions:

1. Read the entry's `run.log`, report and `metrics.tsv`, the test, and the
   history of the test and of the production code it drives since the entry
   last passed (`git log -p` on both). When the history does not decide it,
   build the older commit from a `git archive` export and run the entry there.
2. Classes 1 and 2: repair. Class 4: repair when the defect and its repair are
   evident, otherwise describe it. Class 3: describe what production does now,
   what the test expects, and the proposal (adapt the test, or restore the
   behaviour). Adapt the test only when the commit that changed production
   says the change is intended and the test's guarantee survives it. Class 5:
   `test=profile target=<entry> nseeds=20` (5 at least when time is short,
   reported as a lead), the failure rate, and a proposal: a profile entry, a
   stronger signal, or a fix of the instability's cause. A failure seen only in
   CI is reproduced under the runner's conditions: the fast label on four
   pinned cores, its entries in parallel as `compile_clean.sh` runs them, at
   least 20 times. Its local failure rate is reported beside the CI one.
3. After each repair: build, the fast gate, the entry rerun alone in its own
   scratch directory. A production repair also reruns the fast suite of its
   area and every nightly entry that drives the changed code, and a
   result-changing one passes its gate (section 6).
4. The floors of `simulated_workflow_6vxx` get their origin from the Phase 1
   base profile, in the code and the test card, when the floor rule admits
   them. Otherwise they are a proposal.

Exit: every failure and regression repaired and verified, or described with a
proposal.

### Phase 3: speed

1. Rank the entries by core-hours (wall time times peak CPU) from the
   baseline, and estimate each one's time on the CI runners from the CI logs
   where they hold it.
2. Profile the slow entries' stages (their `BENCH` files, logs, `perf` where it
   helps). Evident speed-ups are made: the same fixture built several times,
   threads reserved but idle, work whose result nothing reads. Each is measured
   before and after under the timing rule (section 7) and keeps the entry's
   results (section 6).
3. Anything else that would shorten the suite is a proposal with its saving.

Exit: every change measured, every proposal priced.

### Phase 4: fusion and style

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
   (`metrics.tsv`), one `!@descr:` per file. Changes are verified as identical
   (section 6).

Exit: fusions done and verified, with the budget, registry, CMake, UI and
policy in step; style changes verified; what was left, with the reason.

### Phase 5: the CI workflow

Bring `ci_test.yml` to what the phases found:
- both platforms always finish (`fail-fast: false`);
- the test reports, `metrics.tsv` files and CTest logs are uploaded as
  artifacts, also on failure;
- explicit job timeouts;
- groups of entries run as parallel jobs, each under its own 6-hour limit,
  where the baseline shows it pays;
- entries the macOS runner cannot hold (registry `mem_gb` against 7 GB) are
  handled explicitly, never silently skipped;
- whether a failed fast gate must keep a group's nightly entries from running
  (Known 2). A change keeps the job red when the fast gate fails.

The report's retention (14 days) against the history a triage needs is
Ruben's (testing plan, section 6): a proposal, not a change.

The YAML is checked locally (parsed, and with `actionlint` if present). A
workflow is only proven by its first run on GitHub; the report says so.

Exit: the workflow changed and checked, or the reasons it was left.

### Phase 6: hot areas and their coverage

1. Churn: over the last eight weeks (`git log --since='8 weeks ago'
   --no-merges --numstat`, renames followed), per file and per subsystem under
   `src/`, testers excluded: commits, lines changed, and commits whose message
   names a fix. Bulk renames and moves are not churn. The hot areas are the top
   twenty files, and the subsystems they belong to.
2. Coverage: a scratch build directory configured with the coverage flags on
   the CMake command line (`--coverage` for Fortran, C and C++, Release
   otherwise; nothing in the repository's CMake changes). Instrument the hot
   files only: with every file instrumented, the threaded entries of the first
   run were 13 to 39 times slower. Run the fast gate, `ctest -L library` and
   `test=nightly` from it; summarise `gcov --json-format` per file and per
   tier (fast alone, and all tiers). The run's timings are not baseline
   timings.
3. For each hot area: its line coverage per tier, which entries reach it, and
   whether a test asserts its behaviour or a workflow only executes it. A hot
   area that only the high-level workflows execute, or that nothing executes,
   is poorly covered.
4. For each poorly covered hot area, a proposed test: the guarantee it pins,
   its tier and suite (a sub-suite of an existing area suite before a new
   entry), how its expected values are derived independently of the code (test
   policy, section 2.1), its negative control, and its cost in run time,
   priced from the code it must run (a whole workflow, a fixture with noise
   spectra) before it is offered as a quick test.

Exit: the churn and coverage table of the hot areas, and the proposed tests.

### Phase 7: verification and report

1. Build from a clean configure in a scratch build directory, the fast gate,
   then the full nightly again (`platform`, `library`, `test=nightly`);
   `test=identity` of every one-thread entry against the Phase 0 baseline,
   except the entries a result-changing repair changes, which pass the Phase 1
   rule instead.
2. The report `nightly_suite_triage_report.md`, in the order of the
   [plan guide](../../how2s/how2write_agent_plans.md) (section 8): what was
   learned, Phase 1's verdicts first; what changed, phase by phase; decisions
   taken under delegation; deviations, and every file changed outside the
   table; the evidence (the status table before and after, profiles, identity
   comparisons); what is left open, as reported-not-gated items,
   data-dependent items (real data) and follow-up runs; a commit message.
   It also says what is not verified (macOS, the CI runners) and how to verify
   it, and resolves each item of Known at writing.
3. Add the run as evidence to section 6 of
   [the consolidated testing plan](../../implementation_notes/planned/testing_and_benchmarking_plan.md#6-consolidated-test-and-ci-work).
4. End with `git status`.

Exit: the second nightly passes, or every failure in it is in the report with
its class and proposal; every Phase 1 verdict is in the report.

## 10. File table

| file | why | phases |
| --- | --- | --- |
| `doc/refactoring_notes/planned/nightly_suite_triage.md` | this plan: section 1 and Progress | each |
| `doc/refactoring_notes/planned/nightly_suite_triage_report.md` | the run's report | 7 |
| `doc/implementation_notes/planned/testing_and_benchmarking_plan.md` | the run as nightly evidence | 7 |
| `doc/policies/test_environment_policy.md` | entries, fusions, budget, timings, the 6vxx card's profile | 2 3 4 |
| `src/main/commanders/test/` | high-level tests, fixtures, registry, runner | 2 3 4 |
| `src/main/exec/simple_test_exec_highlevel.f90`, `src/main/apis/simple_test_exec_api.f90` | routers of the test programs | 2 4 |
| `src/main/ui/simple_test/` | test UI programs and `suite=` help | 2 4 |
| `production/CMakeLists.txt` | CTest registrations, deadlines, budget | 2 3 4 |
| `scripts/check_test_registry.py`, `scripts/ctest_budget.py`, `scripts/run_fast_gate.sh` | the registry lint and the gate | 2 4 |
| `scripts/simple_args_generator.pl` | keys of the test front door | 2 4 |
| `.github/workflows/ci_test.yml`, `.github/workflows/ci_push.yml` | the nightly and push workflows | 5 |
| `src/main/learn/simple_hhcnn_tester.f90` | the CI-only fast-gate failure (Known 2) | 2 |
| `src/main/learn/simple_hhcnn.f90` (if needed) | a production repair of the chunked reduction's thread identity | 2 |

A repair, speed-up or fusion that touches other files adds them below, with the
reason, before editing them.

Added by this run:

| file | why | phases |
| --- | --- | --- |

## Previous runs

- 2026-11-01/02 on the dell, interactive, base `eefc7e1a1`: every red entry
  classified; CI split into four parallel groups and `ci_push.yml` added; one
  production repair (measured FFTW plans only without a fixed seed, for
  one-thread reproducibility); the test-code style review; nine proposed tests.
  Ruled and committed as `4f3fed7e5`. Its Progress is in this file at
  `cfbad5a26` and its report at `4f3fed7e5`
  (`doc/refactoring_notes/planned/nightly_suite_triage_report.md`). Its
  rulings are decisions 6 to 11; its lessons are in sections 7 and 9.

## Progress

(written by the run, one entry per phase)
