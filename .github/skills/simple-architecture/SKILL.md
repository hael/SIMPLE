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

- Every `production/tests/simple_test_*.f90` is auto-globbed into its own
  executable; adding a test needs no CMake edit to build.
- Acceptance is direct execution of the `simple_test_*` binaries (usually with
  `key=value` args parsed via `parse_oldschool`), not CTest. CTest entries are
  an automatic side effect of the same glob loop.
- A test whose fixture requires an external reference volume (`vol1=`) must be
  added to `VOL1_FIXTURE_TESTS` in `production/CMakeLists.txt`; that list
  guards the auto-registration so the test is only registered with CTest when
  `-DSIMPLE_TEST_VOL1=<path.mrc>` is supplied instead of appearing as a
  guaranteed failure.
- New validation fixtures should be modeled on
  `simple_test_continuous_inplane_rotation2D_stage1_validation.f90`
  (builder + pftc setup from a volume, FD gradient checks, PASS/FAIL lines,
  `error stop` on failure).

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
