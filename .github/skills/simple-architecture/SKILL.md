---
name: simple-architecture
description: Use when working anywhere in the SIMPLE repository and you need a high-level architectural map of the cryo-EM platform, including the relationship between production executables, the src core library, src/main subsystems, the strategy-owned distributed masters, the simple_test_exec test environment (fast gate, test registry, nightly/profile/identity), scripts, docs, and the NICE layer.
---

# SIMPLE Architecture

Use this skill first when the request spans multiple parts of SIMPLE.

## Mental Model

SIMPLE is a scientific application platform, not a narrow library.

- `production/` contains thin executable entrypoints: `simple_exec`, `single_exec`, `simple_stream`, `simple_test_exec`, `simple_private_exec` (worker steps, dispatched in `simple_private_exec_driver.f90`) and `simple_persistent_worker`.
- `src/` builds one large static core library (no CMake targets for subsets of the sources).
- `src/main/` contains most application/domain logic; `src/main/solve/` holds the `solve2D`/`solve3D` controllers.
- `src/fileio/`, `src/utils/`, and `src/defs/` provide shared infrastructure.
- `doc/` contains architecture and policy notes that are often more useful than comments in code.
- `nice/` is the optional web/UI layer.

## Runtime Layering

Read these first for end-to-end flow:

1. `production/simple_exec.f90`
2. `src/main/apis/simple_exec_api.f90`
3. `src/main/simple_cmdline.f90`
4. `src/main/ui/simple_ui.f90`
5. `src/main/exec/simple_exec_<topic>.f90`
6. `src/main/params/simple_parameters.f90`
7. `src/main/simple_builder.f90`

The common pattern is:

`executable -> API aggregator -> UI/CLI lookup -> exec router -> commander -> strategy/domain objects -> file/project/image/reconstruction objects`

Rules for the layers, for where each kind of step lives (program, master step, worker
step) and for project reads and writes: `doc/policies/layering_and_project_ownership_policy.md`.

Distributed `refine3D`/`refine2D`: the strategy's master owns its project. Volume
assembly, `prob_align`(`_neigh`), `prob_align2D`, class-average assembly and the
sigma2 commit are strategy routines in `src/main/strategies/parallelization/`
run on the master's `(params, build)`; their commanders are thin `execute(cline)`
wrappers for command-line and worker callers. See
`doc/policies/3D/refine3D_policy.md` (section 4) and `doc/policies/2D/solve2D_policy.md`.

## Build And Docs

- Build shape: `CMakeLists.txt`, `src/CMakeLists.txt`, `production/CMakeLists.txt`, `cmake/CompilerConfig.cmake`
- Generated at build time (never edit): `simple_args.f90` (accepted keys, from `simple_parameters.f90`), `simple_ui_default_values.f90` (`scripts/default_audit.py`), the git hash, and the untracked `doc/code_overview/code_base_map.md` / `test_inventory.md`.
- Repo overview: `README.md`, `doc/code_overview/code_base_map.md`, `doc/code_overview/fortran-indexes/module_index.md`
- Policy docs worth checking before refactors: `doc/policies/**/*.md` (including `compile_time_policy.md`, `test_environment_policy.md`, `ui_layer_policy.md`), `doc/refactoring_notes/{planned,completed}/`

## Tests

- Read `doc/policies/test_environment_policy.md` before adding, moving or
  deleting a test.
- `simple_test_exec` is the only test executable (`BUILD_TESTS=ON`, the
  `compile_*.sh` default; `--exclude-tests` turns it off). Every test is a
  procedure in the library reached through `simple_test_exec test=<name>`.
- A unit test is a subroutine in a `simple_<thing>_tester.f90` module next to the
  code it tests, asserting through `simple_test_utils`, registered with
  `add_suite` in the `suites_<area>` table of
  `src/main/commanders/test/simple_commanders_test_class.f90` and listed in the
  area's `suite=` help in `src/main/ui/simple_test/simple_test_ui_class.f90`.
- A high-level test is a commander (`commanders/test/`), with a router case
  (`exec/simple_test_exec_<area>.f90`), a test-UI program, a CTest registration
  and a row in the test registry (`commanders/test/simple_test_registry.f90`).
- CTest registers one entry per suite, labels `fast`, `library`, `highlevel`,
  `platform`. `SIMPLE_CTEST_BUDGET` in `production/CMakeLists.txt` fixes their
  number and configure fails on a mismatch; checks go into existing suites.
- The 15 `fast` suites run on every `compile_*.sh` build (unless
  `--exclude-tests`) through `scripts/run_fast_gate.sh`, after
  `scripts/check_test_registry.py`, `check_descr.py` and `check_flex_dag.py`;
  a failed gate installs nothing. `library` runs nightly; `highlevel` only when
  selected, or as a whole with `simple_test_exec test=nightly` (also
  `test=profile` over seeds and `test=identity` between run directories).
- A program that only runs commanders or prints results is not a test; use a
  branch or a developer-visibility program.

## Working Style

- Treat `src/main/exec` as routing/orchestration.
- Treat `src/main/commanders` as command objects and high-level workflow owners.
- Treat `src/main/strategies` as algorithm/policy execution layers.
- Treat `src/main/image`, `project`, `ori`, `pftc`, `volume`, `nu_filt`,
  `params`, `ctf`, `motion`, `opt`, `sigma2` and `flex` as core scientific subsystems.
- Expect modern Fortran patterns: modules, abstract types, type-bound procedures, generics, allocatables, submodules, OpenMP, and C interop.

## Cross-Cutting Hotspots

- `refine3D` spans `ui`, `exec`, `commanders`, `strategies`, `volume`, `project`, `pftc`, `sigma2`, and `ori`.
- `refine2D` (a public `simple_exec` program) spans `ui`, `exec`, `commanders`, `strategies`, `class`, `image`, `project`, and `stream` (the 2D pool runs distributed `refine2D`).
- `solve2D`/`solve3D` span `solve`, `commanders`, `strategies` and the refine workflows they stage.
- Streaming work often crosses `stream`, `project`, `motion`, `pick`, `ctf`, and `nice`.
