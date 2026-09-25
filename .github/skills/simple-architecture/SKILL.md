---
name: simple-architecture
description: Use when working anywhere in the SIMPLE repository and you need a high-level architectural map of the cryo-EM platform, including the relationship between production executables, the src core library, src/main subsystems, tests, scripts, docs, and the NICE layer.
---

# SIMPLE Architecture

Use this skill first when the request spans multiple parts of SIMPLE.

## Mental Model

SIMPLE is a scientific application platform, not a narrow library.

- `production/` contains thin executable entrypoints such as `simple_exec`, `single_exec`, `simple_stream`, `simple_test_exec`, and `simple_private_exec`.
- `src/` builds one large static core library.
- `src/main/` contains most application/domain logic.
- `src/fileio/`, `src/utils/`, and `src/defs/` provide shared infrastructure.
- `doc/` contains architecture and policy notes that are often more useful than comments in code.
- `nice/` is the optional web/UI layer.

## Runtime Layering

Read these first for end-to-end flow:

1. `production/simple_exec.f90`
2. `src/main/apis/simple_exec_api.f90`
3. `src/main/simple_cmdline.f90`
4. `src/main/ui/simple_ui.f90`
5. `src/main/params/simple_parameters.f90`
6. `src/main/simple_builder.f90`

The common pattern is:

`executable -> API aggregator -> UI/CLI lookup -> exec router -> commander -> strategies/domain objects -> file/project/image/reconstruction objects`

## Build And Docs

- Build shape: `CMakeLists.txt`, `src/CMakeLists.txt`, `production/CMakeLists.txt`
- Repo overview: `README.md`, `doc/code_overview/code_base_map.md`
- Policy docs worth checking before refactors: `doc/policies/*.md`, `doc/refactoring_notes/*.md`

## Tests

- Read `doc/policies/test_environment_policy.md` before adding, moving or
  deleting a test. The standalone `production/tests` programs, their glob and
  `VOL1_FIXTURE_TESTS` are gone (September 2026).
- `simple_test_exec` is the only test executable (built by every `compile_*.sh` unless `--exclude-tests`).
  A unit test is a subroutine in a `simple_<thing>_tester.f90` module next to the
  code it tests, asserting through `simple_test_utils`, registered with
  `add_suite` in the `suites_<area>` table of
  `src/main/commanders/test/simple_commanders_test_class.f90` and listed in the
  area's `suite=` help in `src/main/ui/simple_test/simple_test_ui_class.f90`.
- CTest registers one entry per suite (labels `fast`, `library`, `workflow`,
  `platform`); `SIMPLE_CTEST_BUDGET` fixes their number. The 13 `fast` suites
  run on every `compile_*.sh` build (unless `--exclude-tests`) through
  `scripts/run_fast_gate.sh` (registry check, 30 s budget); the rest run nightly.
- A program that only runs commanders or prints results is not a test; use a
  branch or a developer-visibility program.

## Working Style

- Treat `src/main/exec` as routing/orchestration.
- Treat `src/main/commanders` as command objects and high-level workflow owners.
- Treat `src/main/strategies` as algorithm/policy execution layers.
- Treat `src/main/image`, `project`, `ori`, `pftc`, `volume`, `nu_filt`,
  `params`, `ctf`, `motion`, and `opt` as core scientific subsystems.
- Expect modern Fortran patterns: modules, abstract types, type-bound procedures, generics, allocatables, submodules, OpenMP, and C interop.

## Cross-Cutting Hotspots

- `refine3D` spans `ui`, `exec`, `commanders`, `strategies`, `volume`, `project`, `pftc`, and `ori`.
- `cluster2D` spans `ui`, `exec`, `commanders`, `strategies`, `class`, `image`, `project`, and `stream`.
- Streaming work often crosses `stream`, `project`, `motion`, `pick`, `ctf`, and `nice`.
