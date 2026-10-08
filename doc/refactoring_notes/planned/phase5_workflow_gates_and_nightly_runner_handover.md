# Phase 5: workflow gates and nightly runner handover

**Status, 2026-10-08:** core workflow gates and scheduled CI are implemented;
closeout is pending the remaining truth checks and a recorded complete nightly
run with its artifacts. Keep this note in `planned/` until the checklist below
is closed or its remaining requirements are explicitly deferred by the owners.
This is the single living Phase 5 record, replacing the September handover's
implementation to-do list. Source audit baseline: `43808e946`.

The original contract is in
[the completed test-environment plan](../completed/uniform_test_environment_refactoring.md),
sections 5.2.2, 5.2.3 and 12. Day-to-day rules live in
[the test policy](../../policies/test_environment_policy.md). Part A makes
simulated workflows fail against known truth; Part B runs and records the
extensive tier unattended. Implementation evidence and observed passing runs
are distinguished below: an implemented assertion is not proof that its latest
workflow has passed on every platform.

## Current test manifest

`simple_test_exec` remains the only test executable. The authoritative registry
is `production/CMakeLists.txt`: **34 baseline entries**, with additional
capability-gated coarray and OpenMP-offload entries. Working directories are
`build/test_runs/<entry>`; registered environments include
`SIMPLE_SEED=20260923` and explicit thread counts.

| Tier | Registered entries | Execution |
|---|---|---|
| `fast` | 15 `unit_<area>` suites | Every normal compile script, before installation |
| `library` | `lib_reconstruction`, `lib_cart_align3D`, `lib_heterogeneity`, `lib_single`, `lib_stream` | Scheduled CI and explicit runs |
| `highlevel` | 13 entries listed below | Scheduled CI explicitly selects this label; excluded from the build gate |
| `platform` | `forked_process`; capability-gated `coarrays`, `openmp_offload` | Scheduled CI selects registered platform tests |

The high-level entries are `mini_stream_6vxx`, `mini_stream_1jxy`,
`simulated_workflow_6vxx`, `simulated_workflow_1jxy`, `single_workflow_fcc`,
`single_workflow_wurtzite`, `pcg_recon`, `simulate_particles`, `solve3D_addon`,
`cont_refine3D_1jxy`, `single_atoms_stats`, `stream_preproc`, and `flex_pca_blobs`.
The old `workflow` label and single combined `single_workflow` registration
are superseded. New entries still require a justified process-budget change;
new checks should normally extend existing suites.

The executable fast-gate budget is now **60 seconds for the whole run**, not
60 seconds per test or a required duration. Failed assertions still fail the
build. The change is limited to `scripts/run_fast_gate.sh` and the default in
`scripts/ctest_budget.py` (`43808e946`); older comments elsewhere may still say
30 seconds. Individual CTest timeouts are separate limits. `lib_stream` now
has 3600 seconds and `flex_pca_blobs` 7200 seconds.

## Part A: implementation and remaining coverage

| Workflow or support | Implemented evidence in current source | Remaining qualification |
|---|---|---|
| Shared truth comparison | `simple_test_truth_metrics.f90` validates dimensions and sampling, docks direct and mirrored maps, selects by docking correlation, and gates finite correlation and masked FSC against declared limits | Use the actual fixture grid and declared floors; retain the distinction between implementing a gate and observing it pass |
| Molecular simulated workflows | `simple_commanders_test_highlevel_workflow.f90` runs 6VXX/1JXY pipelines through final-map truth validation; fixed seed is propagated; `segdiam` and `new` picker routes are accepted | Registered entries select only the default picker. Per-movie CTF/motion and matched-pick precision/recall are not yet equivalent to the full original contract. The routine still uses `NTHR=4` despite an eight-thread CTest allocation |
| SINGLE workflows | `simple_commanders_test_single.f90` has independent FCC Pt and wurtzite CdSe workflows, stage checks, final-map truth comparison, and pairwise pose-error checks; CTest explicitly passes `smpd=0.358` | Record full-run evidence and disposition of the requested final atomic-model comparison; do not treat a separate atom-statistics test as that end-to-end check |
| Stream preprocessing | `simple_commanders_test_stream.f90` retains simulated parameters and optimal averages, matches movies to truth, checks defocus within 0.10 micrometres, frame motion within 0.50 pixels, and image correlation against truth and a wrong-movie control | Confirm the current motion-model changes with a complete workflow run |
| Mini-stream | `simple_commanders_test_highlevel_stream.f90` now provides generated-data 6VXX/1JXY suites with quantitative checks | The original statement that mini-stream is only a manual user-data program no longer describes these registered suites |
| Particle simulation | The same high-level stream submodule checks stack dimensions, counts, sampling, finite nonconstant density, orientation records, CTF metadata and shift bounds | Bounds and nonzero variance are not an independent measurement of requested SNR or proof of exact applied CTF/orientation identity |
| PCG reconstruction | The staged fixed-seed operator gate remains a `highlevel` entry in `simple_commanders_test_highlevel_pcg.f90` | Retain it there until measured runtime and any proposed reclassification are reviewed |
| Add-on, continuous refinement and FLEX | Dedicated high-level gates exercise truth comparisons and application paths; their implementations use `simple_test_gate` and write `metrics.tsv` | Full current runs and archived metrics remain part of nightly sign-off |

`simple_test_gate.f90` defines the machine-readable columns
`name`, `value`, `floor`, `pass`, with finite-value checking and accumulated
failure reporting. The add-on, continuous-refinement and FLEX gates already
use this writer. It is **not yet a uniform output contract across every
workflow**; older workflows still report quantitative results primarily to
logs. Numerical floors remain in source and changes need a scientific reason,
not a platform-specific exemption.

## Work completed and defects addressed this week

These are source/history milestones, not a claim that every affected extensive
suite has passed after every commit. The earlier shared-truth work
(`0c9271c89`, `e443948cd`, `74273134f`) supplies the foundation.

| Change | Evidence | Effect on Phase 5 |
|---|---|---|
| Reproducible molecular and SINGLE fixtures | `30f12c4c5` | Fixed-seed propagation through workflow stages |
| Deterministic atomic density convolution | `931264408` | Removes an OpenMP source of simulation variability |
| Command vocabulary migration | `5b8ad9f06` | Current workflows use `solve2D`, `solve3D`, `refine2D`; historical names are not supported aliases |
| Sparse 2D assignment correction | `9145b617e` | Selects the best frontier for top-1 assignment |
| Stream chain fixture stabilization | `c03b029fc`, `75cba3c19` | Repairs multistate sieve-to-3D and movies-to-3D fixtures |
| Sieve cleanup correction | `97984a77b` | Preserves registered class averages needed downstream |
| Empty-state reconstruction lifecycle | `ed3cb0e41`, `e1577a335`, `bd78813fd` | Retains required volumes, removes stale final artifacts and avoids postprocessing empty states |
| Add-on gate and stream cohort rules | `adfbd36cf`, `85afad3f5`, `90bc27f32` | Aligns class balancing, first ingestion and add-on cohort limits with production behavior |
| Stream initial-analysis selection | `787cc3594` | Prevents a cycle from leaving an empty class selection |
| Motion model migration | `3666eb7fc` | Consolidates motion metadata and adds motion-model unit coverage; extensive workflow compatibility still needs run evidence |
| Continuous-pose aftermath fixes | `bb67eb88d` | Repairs refinement/convergence behavior and updates quantitative tester coverage; details are in the completed pose-cont aftermath note |
| FLEX integration and Linux compatibility | `b744c4a9c`, `457e33c3c` | Integrates fractional state weights/per-state estimation and fixes Linux compilation issues in the PCG fit path |
| Extensive-tier timeouts | `1aaa42ba4`, `92d684c6f` | Gives progressing stream/FLEX workflows appropriate time without changing their assertions |
| Fast-gate portability | `43808e946` | Raises the aggregate time ceiling to 60 seconds for cross-system runtime variation |

The old blanket claims that `lib_single` and stream tests cannot fail are
obsolete: current testers contain assertions, and long nanoparticle work has
its own `single_atoms_stats` entry. Some narrower review findings remain:
`simple_stream_tester.f90` checks picking-reference counts/dimensions but does
not thereby prove every rotation/mirror/normalisation, and its pick/extract
fixture still caps picks at the asserted count. Track those as coverage gaps,
not as unimplemented entire suites. The companion
[single](single_area_tests_handover.md) and
[stream](stream_area_tests_handover.md) handovers retain historical requests;
current source is the authority for their implementation status.

## Part B: nightly execution and handover

### Nightly workflow inventory

All three workflows below also support `workflow_dispatch`. The two SIMPLE
workflows declare cron `0 0 * * *` with `timezone: America/New_York`; the
external dataset workflow declares the same cron without a timezone field.
The SIMPLE YAML comments still describe UTC; the schedule fields above are
what is actually checked in. Scheduled executions can start later than their
nominal time. Workflow definitions establish intended execution; run links
below establish observed outcomes.

| Workflow | Repository / definition | Runners and jobs | Outputs |
|---|---|---|---|
| [Build SIMPLE (Linux/MacOS)](https://github.com/hael/SIMPLE/actions/workflows/ci_build.yml) | `.github/workflows/ci_build.yml` | `Build` and `Build_coarray`, each on `ubuntu-26.04` and `macos-latest` | Build and fast-gate logs; GUI and coarray build results |
| [Build SIMPLE & Run Tests (Linux/MacOS)](https://github.com/hael/SIMPLE/actions/workflows/ci_test.yml) | `.github/workflows/ci_test.yml` | `Build_and_Test` on both platforms | Build, platform, library and high-level CTest logs |
| [Run & Display SIMPLE Data Sets](https://github.com/rmeanapa/SIMPLE_data_testing/actions/workflows/data_testing.yml) | External repository, `.github/workflows/data_testing.yml` | `build_self_hosted` on `[self-hosted, linux]`; dependent `deploy_pages` on `ubuntu-latest` | Downloadable report artifact and published GitHub Pages |

**Build workflow.** Each job checks out the repository with `actions/checkout@v4`
and installs its dependencies. `Build` runs `./compile_debug.sh`, including
the fast gate, followed by `./compile_gui.sh`. `Build_coarray` separately runs
`./compile_coarrays.sh`, which includes the fast gate and the two-image
coarray synchronization check when configured. Ubuntu's coarray branch uses
OpenMPI/OpenCoarrays; macOS installs `opencoarrays` through Homebrew. Both
matrices also declare `node-version: [24]`, but do not invoke `setup-node`.
There is no explicit artifact upload, workflow concurrency group or job timeout
in this definition.

**Build-and-test workflow.** After checkout/dependency setup, it runs
`./compile_clean.sh`, including the fast gate. It then runs, in order:

```bash
ctest -L platform --no-tests=error --output-on-failure
ctest -L library --no-tests=error --output-on-failure
ctest -L highlevel --no-tests=error --output-on-failure
```

The nightly invocation of `highlevel` is explicit; those tests remain outside
the ordinary compile gate. The workflow respects each test's registered
`RUN_SERIAL` and timeout settings. Failed steps stop later work under the
normal Actions shell behavior; default matrix fail-fast can cancel the other
platform, so a cancelled job is not a passing platform result.

For the regular Linux jobs, apt dependencies include GCC, CMake, image codecs,
FFTW, gnuplot, SQLite, CUDA toolkit, MPICH, curl, OpenBLAS, LAPACK and ARPACK.
The macOS jobs install Homebrew dependencies including jbigkit, Python 3.10,
FFTW, curl, OpenBLAS, LAPACK and ARPACK. These are workflow dependencies, not
proof that a hardware offload test was registered or executed. Neither SIMPLE
nightly workflow has a custom JUnit/metrics collector or explicit report
artifact upload.

**Real-data workflow.** The self-hosted job has a **1440-minute** timeout.
After checking out `SIMPLE_data_testing`, it updates the separate SIMPLE
checkout at `/home/meanapanedar2/SIMPLE` using `git pull --ff-only` and runs
`./compile_clean.sh`. Runtime settings are `SIMPLE_QSYS=local`,
`SIMPLE_PATH=/home/meanapanedar2/SIMPLE/build`, and the build's `scripts` and
`bin` directories prepended to `PATH`; `SIMPLE_EMAIL` is set in the YAML.
The seven dataset scripts run **sequentially** with `bash -e` in this Actions
workflow. The independent Slurm launcher described below is a separate route.

After the real datasets, the workflow runs four generated-data gates:

```bash
simple_test_exec test=simulated_workflow suite=6vxx > LOG_6vxx
simple_test_exec test=simulated_workflow suite=1jxy > LOG_1jxy
simple_test_exec test=single_workflow suite=fcc > LOG_single_fcc
simple_test_exec test=single_workflow suite=wurtzite > LOG_single_wurtzite
```

`make_report.sh` receives the seven dataset names and four generated-workflow
result directories. It prepares `build/`, which is uploaded as
`simple-data-sets-report-${{ github.run_number }}` by `upload-artifact@v4`,
with **31-day retention** and `if-no-files-found: error`. The same directory
is uploaded with `upload-pages-artifact@v4`. The dependent deployment job uses
`deploy-pages@v4`, the `github-pages` environment, and Pages/id-token write
permissions to publish [the dataset report](https://rmeanapa.github.io/SIMPLE_data_testing/).
These report/upload steps do not use `if: always()`, so a failed test can
prevent the report and artifact from being produced. The data-repository
commit alone does not identify the SIMPLE source: archive the separate
pulled SIMPLE commit as well.

**Documentation deployment is not nightly.**
`.github/workflows/docs_build_and_deploy.yml` runs manually or on pushes to
`master` affecting `doc/**` or `mkdocs.yml`. It installs Python/MkDocs Material,
uses a weekly-keyed pip/cache directory, and runs `mkdocs gh-deploy --force`.
Its `contents: write` permission and resulting Pages deployment are for the
documentation site, not test validation. Do not count a green documentation
or `pages-build-deployment` run as a green build-and-test night.

### Observed GitHub Actions evidence (2026-10-08)

The following statuses were read from GitHub Actions during this update;
in-progress runs need their final verdict added before sign-off.

| Run | Commit | Observed result |
|---|---|---|
| [Build SIMPLE, 37801608397](https://github.com/hael/SIMPLE/actions/runs/37801608397) | `457e33c3c` | **Success**: Linux and macOS Debug/GUI jobs and both coarray jobs completed successfully |
| [Build and test, 37822274500](https://github.com/hael/SIMPLE/actions/runs/37822274500) | `43808e946` | **In progress**; started 18:10:49 UTC |
| [Previous build and test, 37726150380](https://github.com/hael/SIMPLE/actions/runs/37726150380) | `92d684c6f` | **Failure** in Ubuntu's test step; macOS matrix job cancelled |
| [Data sets, 37824252590](https://github.com/rmeanapa/SIMPLE_data_testing/actions/runs/37824252590) | data repo `04c8cc1d8` | **In progress**; started 18:26:10 UTC |
| [Previous data sets, 37724747677](https://github.com/rmeanapa/SIMPLE_data_testing/actions/runs/37724747677) | data repo `3d4d56b21` | **Failure**; does not provide full benchmark sign-off |

This closes the absence of any recorded cross-platform build evidence in the
old handover. It does not close the full extensive-tier run requirement. The
source contains scheduled automation rather than a dedicated
`scripts/nightly_run.sh`; acceptance of Actions as the replacement operational
contract should be explicit.

### Original runner contract: remaining differences

| Original requirement | Current disposition |
|---|---|
| Scheduled build and extensive-tier execution | Implemented in GitHub Actions; latest complete success still needs a run link |
| Known source and environment provenance | Checkout and CI logs exist; no consolidated compiler/CMake/host manifest is written; the data job also pulls a separate SIMPLE checkout |
| Prevent overlapping nights and bound the whole run | No explicit concurrency group in the inspected definitions; SIMPLE jobs have no explicit timeout, while the self-hosted data job declares 1440 minutes |
| Library/platform/high-level status and times | CTest reports them in logs; no checked-in runner collects a complete per-entry history |
| Four-thread reverse-order fast tier | Not invoked by the current workflow; must run outside CTest's one-thread environment |
| Dated archive outside the disposable build tree | SIMPLE workflows retain Actions logs but have no explicit test-artifact upload/history writer; the data workflow uploads its report for 31 days and deploys Pages after success |
| Uniform metric collection and failed-log tails | Shared metrics writer exists, but workflow-wide collection/summary is outstanding |
| Notification and install/runbook ownership | Agree the accepted Actions notification/archive destination and whether this replaces the dedicated-machine contract |
| Weekly coverage | Optional; not a prerequisite unless retained by owner decision |

A green fast gate is not a green nightly run. The original requirement to
record every library/high-level entry's runtime remains open until a complete
run is attached. Do not invent timings from CTest timeout values.

## Observed validation and cluster handover

On 2026-10-08, x3's Release build initially passed all 15 fast suites but took
37.19 seconds, exceeding the former 30-second ceiling. `unit_stream` took
37.12 seconds, including 17.90 seconds in file-heavy watcher fixtures.
A temporary local-fixture experiment reduced the gate to 6.49 seconds, but
that implementation was explicitly **reverted** in favour of the 60-second
budget. There is no retained `SIMPLE_TEST_FIXTURE_ROOT` feature.

After removal of that experiment, requested incremental Release rebuilds and
installation on x3 succeeded with normal fixtures: all 15 fast suites passed
in **27.21 seconds** and **28.61 seconds**. These are observed Linux Release
results, not macOS Debug, BOX, offload or full nightly results. The last build
log was `/tmp/simple-gate-doublecheck.sXQK3QBv.log` on x3; temporary logs are
not a durable nightly archive. The build-tree report is
`build/test_runs/ctest_fast.log.timing.txt` and is overwritten on reruns.

The separate `SIMPLE_data_testing` repository has a Slurm launcher for seven
real datasets: betagal, apof, clc, motab, proteasome, trpm4 and nanox. Each is
submitted independently, with local SIMPLE workers inside its allocation,
80 CPUs, 128 GB and a 24-hour limit. The observed cluster partition is `norm`
(the initially requested `normal` is not its name). Each invocation creates
a unique results directory; dataset program output is mostly redirected to
`<dataset>/LOG`, not the Slurm stdout/stderr files.

The first observed seven-dataset run is **not a successful benchmark gate**:
six jobs returned exit code 1 during `solve3D_cavgs` within a 15-second window;
NanoX was still running at inspection. The home-directory allocation reported
256 GB used and zero available, while the run occupied about 42 GB. This is
strong evidence for a storage-related interruption, not proof of a numerical
regression. Results should use an adequately provisioned shared data directory
such as the user's BeeGFS workspace, with capacity checked before submission.
No successful rerun or final FSC verdict was observed. This operational work
complements, but does not replace, generated-data nightly validation.

## Closeout checklist before moving to completed

- [x] Replace the stale 25-entry/13-fast manifest with the current 34-entry registry.
- [x] Record implemented map-truth, stream-preprocessing, SINGLE pose and newer application gates.
- [x] Record scheduled CI, this week's defect fixes and verified x3 build/fast-gate results.
- [ ] Attach a complete current Linux/macOS nightly run, commit, per-entry verdicts and runtimes; explicitly identify unsupported platform cases.
- [ ] Finish or explicitly defer the molecular workflow's per-stage truth checks, second-picker execution and thread-count mismatch, particle-simulation SNR/identity checks, and final SINGLE atomic-model comparison.
- [ ] Finish or explicitly revise the nightly archival, uniform metrics, overlap control, total timeout, threaded reverse-order and notification requirements.
- [x] Link successful Linux/macOS Debug, GUI and coarray build jobs (`37801608397`).
- [ ] Record the carried-over `BUILD_TESTS=OFF` link check and hardware-offload evidence; capture compiler/bounds flags from the successful macOS Debug log if exact gfortran-runtime sign-off is required.
- [ ] Have the owners accept the remaining scope decisions; then move this same note to `completed/` and update references.

Hans reviews numerical floors and new CTest entries; Ruben owns runner closeout;
Joseph owns the stream production behavior under test. The outstanding items
are a finite handover list. Moving the file is the final bookkeeping step,
not a substitute for those decisions or their evidence.
