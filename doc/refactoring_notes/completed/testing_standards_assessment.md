# SIMPLE's testing against modern practice: an assessment

Date: 2026-09-24

Status (2026-09-25): **closed.** Every recommendation was either carried out,
handed to the document that already plans it, or dropped with a reason (table
below), and the note moved to `completed/`. Sections 1 to 4 are kept as written
on 2026-09-24, the morning before the test review finished; several of their
facts are out of date, and the current state is in the policy
(`doc/policies/test_environment_policy.md`), the archived plan
([`uniform_test_environment_refactoring.md`](uniform_test_environment_refactoring.md))
and the review record
([`test_review_record.md`](../../code_overview/test_review_record.md)).

| Recommendation | Outcome |
|---|---|
| 1. Fast gate on every push | Not enforced, by decision (Hans): a push sometimes only moves code to a cluster machine. Running the gate before a push is a strong recommendation in the policy (section 3.4) and `AGENTS.md`; the nightly CI gates master on Linux and macOS in Debug (`compile_debug.sh`, `compile_gui.sh`) and Release (`compile_clean.sh`), all with `--compile-tests`. The hand-kept CI test list went with Phase 7: CI runs `ctest` by label and name. |
| 2. Simulation-truth floors | Phase 5, Ruben: `doc/refactoring_notes/phase5_workflow_gates_and_nightly_runner_handover.md`. |
| 3. Bounds checking on macOS Debug, FP traps | Linux Debug builds, the developers' and the nightly CI's, run `-fcheck=all`; the SPIDER header defect was caught by a bounds-checked build. macOS keeps the reduced checks until someone isolates, on a Mac, which gfortran 16 descriptor check crashes on the FFTW-backed pointers. FP traps and signalling-NaN initialisation are an opt-in: `./compile_debug.sh --traps` (`SIMPLE_DEBUG_TRAPS`), for a diagnostic run. |
| 4. Warnings visible | `-w` removed from the Release Fortran flags (vendored sources keep their own). Warnings-as-errors: not before the warnings in a Release build have been seen and cleared. |
| 5. Test rules where contributors read them | The policy; a Tests section in `AGENTS.md`; a pointer in `.github/copilot-instructions.md`; the architecture skill points to the policy. |
| 6. Coverage once | In the Phase 5 handover: an optional weekly instrumented build of the fast and library tiers with a `gcovr` summary in the nightly archive, as a report. |
| 7. Generative tests | Round trips at the sizes where a file layout changes: a rule in the policy (section 4.5) and `test_file_roundtrip_sizes` in `simple_image_tester` (MRC and SPIDER, odd, small and non-square images, stacks and volumes). Corruption fuzzing of the readers: dropped for now; the readers' inputs are SIMPLE's own and the common external formats, and the defect that prompted the rule was found by a round trip at a small box. Property tests: done where the invariants are natural (symmetry groups, rotations, address maps, transfer functions). |
| Parallel correctness (section 2) | A rule in the policy (section 4.5): a threaded path gets a check that forces a team; the Phase 5 runner repeats the fast tier on four threads each night. |

Validation level: static inspection only. Read for this assessment:
`.github/workflows/ci_build_and_test.yml`, `CMakeLists.txt`,
`production/CMakeLists.txt`, `cmake/CompilerConfig.cmake`, the
`compile_*.sh` scripts, `scripts/run_fast_gate.sh`, `AGENTS.md`,
`.github/copilot-instructions.md`, the refactor record and the inventory.
Nothing was compiled or run. Counts and timings are the refactor record's and
the inventory's own, as of the dates they carry.

## 1. Summary

SIMPLE is in the middle of a test overhaul modelled on X's, and its
*structure* is now close to modern practice: a fast gate that is part of the
build, a process budget, hermeticity rules, fixed seeds, and a nightly CI on
two operating systems that runs that gate. Its weaknesses are elsewhere. A
large share of the suite still cannot fail — it prints and passes. The
simulated workflows are checked for completion, not against the truth they
were simulated from. The reference development machine runs Debug without
bounds checking, and Release builds silence every warning. The review now
under way is closing the first gap and finding real defects as it does; the
others are cheap to close once named.

## 2. The standards, and where SIMPLE stands

| Practice | What current practice expects | SIMPLE, 2026-09-24 | Verdict |
|---|---|---|---|
| Fast automated suite | A suite of seconds to minutes, run on every change | Thirteen `unit_<area>` suites under the `fast` label, each one process, pinned to one OpenMP thread, with a 30 s budget checked by `scripts/ctest_budget.py`. Every `compile_*.sh --compile-tests` runs the gate between build and install, and a failed gate installs nothing and fails the script (`scripts/run_fast_gate.sh`). The gate passed in 3.3 s real when it was declared with seven suites (2026-09-22) | Strong |
| Tests that can fail | Every test ends in an assertion whose failure reaches the exit status | The 2026-09-22 audit found 84 of 150 test identities with no failure path on any route: they printed and passed unless they crashed. Only 36 sources used `simple_test_utils` assertions. The inventory now reads 63 of 130 identities without a failure path, and 36 still implemented on two routes (a standalone program and a commander case), down from 53 | The central weakness; closing |
| Process budget | Isolation units counted and bounded | `SIMPLE_CTEST_BUDGET` = 25 (13 fast, 5 library, 6 workflow, 1 platform); configure fails on a mismatch | Strong |
| Hermetic and deterministic | No dependence on clock, network, machine state or other tests; seeded randomness | Written rules (refactor record section 5); `SIMPLE_SEED` under CTest and a reseed before every sub-suite; the forked-process and socket cases moved out of the fast gate; sleeps removed from the IPC gate; every fast suite also passes in reverse sub-suite order | Good |
| Continuous integration | Every push builds and tests automatically; failures block or at least flag the change | GitHub Actions on `ubuntu-26.04` and `macos-latest`: Debug, GUI, coarray and clean builds, each with `--compile-tests`, so the fast gate runs in each and fails the job. Triggers are the nightly schedule and manual dispatch only. The clean-build job then runs a hand-kept list of standalone programs and `simple_test_exec` calls, some of which are identities the inventory lists as having no failure path, so they cannot turn the job red; the list has drifted as tests migrated (the record notes CI losing `neigh`, `ctf` and `lplims`, and its `test=units` call stalling until 5f3c705ce) | Present and gating nightly; not per push |
| Runtime checking | Bounds and pointer checks, FP traps and uninitialised-value detection in Debug; sanitizers in CI | Linux Debug: `-fcheck=all -fbounds-check`. macOS Debug: `-fcheck=do,mem` only, because gfortran 16's debug runtime segfaults on the C-interoperable FFTW pointers — so the reference Mac, where the fast gate is budgeted, runs without bounds checking. No FP traps, no `-finit-real=snan`, no sanitizers | Weak on the reference platform |
| Compiler warnings | Warnings visible in every build; warnings as errors in CI | `-Wall` in the base flags, `-Wuninitialized -Wunused` in Debug; Release adds `-w`, which suppresses every warning. No warnings-as-errors option | Weak |
| Coverage | Line or branch coverage measured, used to find unexercised code | A call-footprint proxy (`scripts/test_review_dossier.py`, record section 9.5): before any test is deleted, the modules and type-bound procedures that would lose all tests are listed and answered. Not line coverage; an instrumented `--coverage` build is noted as a later option | Partial, deliberately |
| Oracles | Expected values independent of the code under test | Varies by test; the review adds assertions area by area. For the workflows the natural oracle exists and is unused (next row) | Mixed |
| Scientific validation | Results checked against a known truth | Six `workflow` entries run simulated pipelines; the `simulated_workflow` pair generates its data from embedded atomic models (6VXX, 1JYX) with known orientations, defocus and B-factor. The workflows check completion — files exist, counts match, the heartbeat is well formed. Gating on the truth (FSC to the truth map, pose recovery, declared floors) is phase 5 of the refactor and not yet built. A real-data validation registry of X's kind is an explicit non-goal | Planned; the right design |
| Parallel correctness | Parallel results agree with serial | Parallel-equals-serial checks in some batches (six routines in one, a relative-L2 tolerance in another); no general rule that every parallel path ships with one | Partial |
| Generative testing | Fuzzing of file readers; property-based tests of invariants | None found | Gap |
| Performance tracking | Timings recorded and compared across builds | The fast gate writes a per-entry timing table beside its log on every build (`ctest_fast.log.timing.txt`); the nightly archive of suite times is part of phase 5 | Emerging |
| Rules where contributors read them | Test rules in the contributor instructions | The admission rules live in the refactor record; `AGENTS.md` and `.github/copilot-instructions.md` carry no test guidance beyond "do not compile", and there is no testing skill under `.github/skills/` | Gap |
| Merge gate | Changes reach the main line through a gate | Focused commits pushed to `origin/master` (the XP workflow in `AGENTS.md`); the nightly CI flags after the fact | By choice |

## 3. What the review has already found

The value of making tests able to fail is not hypothetical. By the refactor
record's count, the tests written in the review have found and fixed
thirteen production defects, among them: the mask mirror asymmetry and the
disc padding count; the Otsu bin edge and two-valued input; `norm_2` and
`vabs` returning 0 on macOS through Accelerate's `snrm2`; the reversal of
even-length double arrays; the Nystroem kPCA feature and projection
scaling; the non-convergent cosine pre-image; and the residual BIC of the
PPCA rank scan. The review also removed thirty-one dead routines and three
unused optimisers, and the subproject scheduling framework that only two
tests still called; a queue-environment test turned out never to have run.

Every one of those defects sat under tests that printed plausible output and
passed. That is the concrete cost of the no-failure-path problem, and the
strongest reason to finish the review before anything else on this list.

## 4. Recommendations, in order of value for effort

1. **Run the fast gate on every push.** Add a `push` trigger for
   `master` to `.github/workflows/ci_build_and_test.yml`, for one job at
   first (Linux Debug, `--compile-tests`), and keep the full matrix
   nightly. The gate itself is seconds; the build with tests dominates.
   In the same edit, drop from the clean-build test step every entry the
   inventory marks as having no failure path — a step that cannot fail
   makes a green run mean less than it appears to — and let the list
   shrink with the review until phase 7 replaces it with `ctest` labels.
2. **Build the simulation-truth floors (phase 5).** FSC against the map
   simulated from the same coordinates, pose and shift recovery after
   symmetry and hand alignment, each with a declared floor versioned with
   the test. This is SIMPLE's counterpart to real-data validation at no
   data cost, since the truth is generated with the data; it is the
   single highest-value scientific item here.
3. **Restore bounds checking on macOS Debug.** Compile only the
   FFTW-binding modules with the reduced checks (per-source
   `COMPILE_OPTIONS`) and give the rest of the library `-fcheck=all`; add
   `-ffpe-trap=invalid,zero,overflow` and `-finit-real=snan` to Debug on
   both platforms, and run the fast gate once to see what they flag.
4. **Make warnings visible.** Remove `-w` from
   `CMAKE_Fortran_FLAGS_RELEASE` and add an opt-in warnings-as-errors
   option for CI, so a warning introduced on one platform is seen on both.
5. **Put the test rules where contributors read them.** A short testing
   section in `AGENTS.md` (admission: can fail, hermetic, one process per
   suite, the budget ratchet, verify in Debug and Release) and a
   `simple-tests` skill under `.github/skills/`, pointing to the refactor
   record for the design rather than repeating it.
6. **Check the coverage proxy once.** An instrumented (`--coverage`)
   build of the fast and library tiers with `gcovr`, read against the
   dossier's call-footprint report, to learn whether the proxy misses
   whole branches.
7. **Later: generative tests.** A seeded corruption generator (truncate,
   flip bytes, drop fields) over the MRC, STAR and `.simple` project
   readers — formats other programs write — asserting a clean error or a
   valid parse, never a crash. Property tests on geometry, where the
   invariants are natural: rotations compose and invert, symmetry groups
   close, the reprojection of a rotated volume equals the rotated
   reprojection.

