# Compile-Time Policy

## Purpose and scope

SIMPLE is a large modern-Fortran code base built with GNU Fortran (gfortran)
through CMake's Make generator. A clean build compiles about 840 sources. This
policy says how to keep build time down in day-to-day development and how to
reduce it when it has grown. It applies to every SIMPLE-owned source, to the
CMake files and to the build scripts. Every rule below comes from a measured
experiment; the evidence is summarized at the end.

Two quantities matter and they are not the same:

- **Wall time**: how long a clean build takes. This is what developers wait for.
- **Summed compile time**: the per-file compile times added up. It measures the
  total work, not the waiting.

With few cores (a laptop, `-j11`), wall time follows the total work. With many
cores (the 80-thread Linux build host), wall time follows the **critical
chain**: the longest sequence of compiles in which each must wait for the one
before, because it imports that file's module. A change can cut summed time
substantially and leave wall time unchanged, or the reverse.

## Where the time goes

Three mechanisms explain nearly all of SIMPLE's compile time.

1. **Inline derived-type graphs.** When a type holds another large type as a
   plain (non-allocatable, non-pointer) component, the compiler expands the
   initialization, copy and clean-up of the whole nested tree as straight-line
   code: in the module that defines the type, at every intrinsic assignment and
   at the end of every procedure that declares a local of the type. Types such
   as `parameters` (about 100 string components and 195 fixed-length character
   components) and `cmdline` (a fixed array of 100 key/value records) make every
   enclosing type expensive. A test module that declares 14 local stage objects
   pays that cost 14 times.
2. **Module dependency chains.** A source cannot compile until every module it
   imports at module level has been compiled. One slow file in a chain delays
   everything after it. Umbrella modules that re-export many others
   (`simple_core_module_api`, `simple_commanders_api`) turn a few real
   dependencies into dependencies on everything they re-export.
3. **Target boundaries.** With the Make generator, a target that depends on
   another target (for example by linking it) does not start compiling until
   the whole other target is built, even if only a few of its modules are
   needed. Splitting sources into extra targets therefore creates waiting, not
   parallelism.

## Rules for day-to-day development

These rules prevent compile-time growth. Follow them when writing or reviewing
code; they cost nothing at run time.

### Types and state

- **Do not embed heavy types inline.** A component of type `parameters`,
  `cmdline`, a stage, a sieve or a similar aggregate is `allocatable`, allocated
  where it is first needed (right before its constructor) and deallocated in
  `kill`, on both the full and the early-return path. Keep `new`/`kill`
  symmetric and idempotent.
- **Do not snapshot `parameters`.** An object that needs a few settings copies
  those fields into its own components. Copying a whole `parameters` object
  (`self%params = params`) costs compile time at the copy and run time at every
  call, and usually hides that most of it is never read.
- **Remove state that is never read.** An assignment with no reader is dead
  state, not a safety copy.
- **Keep the type small that many places declare.** The cost of a type grows
  with the number of procedures that declare it as a local. A type declared in
  many tests or helpers must stay light.
- **Put state where it is used.** Module variables used only by one submodule
  are declared in that submodule. State shared by several submodules and not
  by the parent goes into a small companion module (for example
  `simple_nu_filter_vars`) that the parent imports. Do not declare it in the
  parent: gfortran reports it as unused, and it widens the parent's interface.

### Tests

- **Declare large test objects as polymorphic allocatables.** In testers, a
  local stage, sieve or similar object is `class(T), allocatable :: x`
  followed by `allocate(x)`; helpers take `class(T)` arguments. Its clean-up
  then goes through one routine the compiler generates once per type, instead
  of being expanded in every test procedure.
- **Keep one test-only file from ending the build.** A large test commander
  that nothing else imports is split into a parent module (the commander types
  and the interfaces of their `execute` procedures) and submodules grouped by
  test family, so the families compile in parallel. This applies only to such
  leaf files (see "What does not work").

### Imports

- **Import with `use ..., only:` and only what the module itself needs.** Every
  module-level import is a build-order edge for every file that imports this
  module. An import needed by one procedure goes inside that procedure, or into
  the submodule that holds it.
- **Do not add new umbrella modules or re-exports.** New code imports the
  modules it uses, not an umbrella. Existing umbrellas stay as they are (see
  "What does not work").
- **Import only the names a third-party module provides that you use.** A bare
  `use FoX_dom`, for example, brought in its own `len`, which caused a warning,
  and its `//` for joining text and integers, which code had come to rely on
  instead of `int2str`. Name the procedures and types you call.

### Build configuration

- **One library target for SIMPLE.** Do not create additional targets for
  subsets of the sources (tests, vendored code, a subsystem) to give them
  different flags or to organize them; with the Make generator each new
  dependent target waits for the whole target it depends on.
- **No per-file compiler options for SIMPLE-owned sources.** Warnings are fixed
  in the code. The one exception is a single `-w` on the bundled third-party
  sources (`SIMPLE_VENDOR_SOURCES` in `src/CMakeLists.txt`), which stay in the
  library target.
- **Release flags are `-O3 ${ARCH_FLAG} -fPIC`.** `-Wmaybe-uninitialized` is off
  for all build types (it produces false positives on the hidden fields of
  allocatable arrays). A new global optimization flag needs a run-time
  benchmark that shows a gain worth its compile-time cost; `-funroll-loops`
  did not (see the evidence).
- **Builds use a bounded job count.** The `compile_*.sh` scripts use the number
  of available processors or `CMAKE_BUILD_PARALLEL_LEVEL` (`compile_clean.sh`
  also takes `--jobs N`); never an unbounded `make -j`.
- **Day-to-day builds may leave out the tests** (`./compile_clean.sh
  --exclude-tests`). Test code is about a fifth of the summed compile time.
  Builds that precede a commit include the tests and the fast gate.

## Measuring

Use `scripts/profile_build.sh`. It configures and builds from scratch and records
each compile's start, end and duration.

```bash
./scripts/profile_build.sh clean --label <name> -j <N>       # clean Release build, tests on
./scripts/profile_build.sh clean --label <name> --no-tests -j <N>
./scripts/profile_build.sh compare build_profile/<before> build_profile/<after>
./scripts/profile_build.sh touch <source.f90> --label <name> -j <N>   # rebuild after one edit
```

Profiles are written to `build_profile/` (not tracked). Then run the fast gate:
`ctest --test-dir build -L fast --output-on-failure`.

Rules for comparisons:

- **Like for like.** Same machine, compiler, build type, test inclusion and job
  count; one change between the two profiles.
- **Check the machine.** Compare the times of the files the change did not
  touch. If they moved by more than a few percent (a busy or warm machine), the
  wall time of the pair is not comparable; repeat the profile.
- **Read the shape, not only the total.** From `compile.tsv`: how many compiles
  run in each 10 s interval (concurrency), which compiles end last (the tail),
  and the chain behind the last compile (for each file, the imported module that
  finished last before it started). A gap of low concurrency names the problem;
  the total does not.
- **Summed compile time is wall time per command added up**, not CPU
  accounting. Use it to compare work, never to predict wall time.
- **Single samples.** A profile is one measurement. Do not claim a wall-time
  gain under about 5 s from one pair.

## Reducing compile time when it has grown

1. Profile a clean build and the fast gate on the current source.
2. Classify the problem from the profile:
   - **One file far slower than its peers** (tens of seconds): look for a large
     inline type graph in it or in the types it declares (rules under "Types
     and state" and "Tests").
   - **Low concurrency at the end of the build:** find the file that ends the
     build and the chain behind it. If the last file is a leaf that nothing
     imports, split it into submodules. If the chain runs through modules that
     others import, the fix is fewer module-level imports in those modules,
     not splitting them.
   - **Low concurrency at the start or in the middle:** look for a target
     boundary or a module that everything imports and that compiles slowly.
   - **High concurrency throughout, but slow:** the total work is the problem;
     reduce the slowest files one by one, each with a stated cause.
3. Make one change, profile again on the same machine, run the fast gate.
4. Keep the change if wall time drops by at least 5 s, or if it removes a clear
   local cost (a file several times faster) without adding structure. Otherwise
   revert it.
5. Record what was tried and the numbers, including what was reverted, in the
   change's notes.

## What does not work

Each of these was built, measured and reverted.

| Approach | Result |
| --- | --- |
| Test-only sources in their own library target, compiled at `-O1` | Test compile time fell from 270 s to 165 s, wall time only from 179 s to 171 s: every test source waited for the whole main library. |
| Bundled third-party sources in their own object-library target | No SIMPLE source started during the first 14.7 s of the build; concurrency in the first 20 s fell from 3.9 to 1.7 jobs. |
| Splitting chain-link commanders (`refine3D`, `solve3D`) into parent and submodules | No gain: the parents still waited for the modules they themselves import. |
| Contract/implementation split of the whole commander layer (September 2026) | 16 % faster clean build, no gain for edits or interface changes, 56 extra files; rejected (`doc/refactoring_notes/rejected/contract_submodule_architecture.md`). |
| `-funroll-loops` in the Release flags | Kernel suites unchanged, one end-to-end workflow about 3 % faster: not worth the compile time. |

## Evidence

Clean Release builds with tests; GNU Fortran 16.2.0 unless stated.

| Change | Machine | Measurement |
| --- | --- | --- |
| Initial-analysis stage: sieve and parameters allocatable | Linux, 80 jobs | Stage source 91.6 s to 16.0 s; build 225.3 s to 138.4 s |
| Stream testers: `class(T), allocatable` stage locals | Mac, 11 jobs | Eight testers 74.9 s to 30.4 s summed; initial-analysis tester 18.9 s to 4.2 s |
| Vendor sources back in the library target, single `-w` | Mac, 11 jobs | First SIMPLE compile at 0 s instead of 14.7 s; about 9 s of wall time |
| High-level test commander split into seven submodules | Mac, 11 jobs | Its span 18.0 s to 7.2 s; it no longer ends the build |
| Whole follow-up batch | Mac, 11 jobs | 179.3 s to 159.9 s wall, 1375 s to 1256 s summed, no warnings, fast gate passing |
| Both batches (commit 8aedefec5) | Linux, 80 jobs, GNU Fortran 15.2.1 | 118.7 s wall, 1358 s summed; the same host took 225.3 s before the first batch (GNU Fortran 16.2.0) |

The 2026-10-04 profiles are kept in `build_profile/` on the machine that made
them. The first change is reported in full in
`doc/refactoring_notes/completed/stream_owned_state_compile_time_report.md`.
