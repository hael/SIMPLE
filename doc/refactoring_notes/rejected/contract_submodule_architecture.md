# Contract/Implementation Module Architecture for Fast Builds

Date: 2026-09-14

**Status:** Rejected 2026-09-21 after measurement (written 2026-09-14). Kept
for the record; not a plan. The commander layer was converted on a branch and
built against master: clean builds were 16% faster, body edits and core
interface changes not at all, which does not justify the added structure.
Section 11 has the numbers. The text above it is left as written; read it
with this header.

Purpose: one set of rules, explainable in five minutes, that gives SIMPLE a
clean layered architecture and fast clean and incremental builds, without
repeated `use` lists or long `public` lists in every module.

## 1. Decision summary

Every module above the utility layer is split into a **contract** (the module
file: types, constants, interfaces, nothing executable) and one or more
**implementations** (submodules holding all bodies). Umbrella modules such as
`simple_core_module_api` and `simple_commanders_api` are kept, but they are
consumed only by implementations and programs, never by contracts. Layers
depend downward only. Tests leave the default build target and the build
moves to Ninja.

The effect is structural, not cosmetic: a body edit recompiles exactly one
file; a contract edit stops cascading after one hop because contracts import
almost nothing; `.mod` files shrink by an order of magnitude because nothing
re-exports an umbrella; and on a clean build all bodies become leaves of the
dependency graph and compile fully in parallel behind a thin chain of small
contracts.

## 2. Why the build is slow today

Measured on the 2026-09-14 incremental profile (369 files, 436 s CPU) and the
current `build/modules` directory:

- `simple_core_module_api` re-exports 35 modules and is `use`d without
  `only:` by ~220 files; `simple_commanders_api` re-exports parameters,
  builder, image, sp_project and more to 68 commander files.
- Modules default to public and `use` umbrellas without `only:`, so every
  module's `.mod` re-exports its whole closure: `simple_commanders_api.mod`
  is 345 KB and `simple_commanders_refine3d.mod` 371 KB gzip-compressed; the
  module directory is 251 MB across 1217 files.
- Two consequences follow. An interface change inside the umbrella closure
  (a new `oris` method) recompiles most of the tree — the profile above was
  such a cascade: 369 of 621 sources, with the core modules themselves not
  even among them. And every file that `use`s an umbrella decompresses and
  parses megabytes of symbols first, which is why the 533-line
  `single_commanders_tseries.f90` takes 13 s.
- The critical path (module graph weighted by compile time) is ~37 s and
  ends in a single 2.9k-line file, `simple_commanders_test_highlevel.f90`
  (20 s); `stream_p03_initial_analysis` -> `ptcl_sieve` is 24 s. Fat late
  files, not the sum, bound the wall time once cores are plentiful.
- `BUILD_TESTS=ON` puts ~90 test programs (152 s CPU, 35% of the profile,
  plus one full static link each) into the default target.

## 3. The rules

There are five. A developer needs to know all five and nothing else.

**Rule 1 — Layers depend downward only.**

```
programs      production/*.f90                       (simple_exec, simple_stream, tests)
exec          main/exec, main/apis/*_exec_api        (command dispatch)
workflow      main/commanders, main/stream, main/strategies, main/sieve, ...
domain        main/pftc, main/nu_filt, main/flex, main/motion, main/ctf, main/pick, ...
core          main/image, main/ori, main/project, main/params, simple_builder, simple_cmdline
utils         utils/, fileio/, defs/, extlibs/
```

A file in a layer may depend on its own layer and on layers below it. Never
upward. Fortran already forbids cycles; this rule forbids the long way round.

**Rule 2 — A module file is a contract.**

Above the utils layer a module file contains only: derived-type definitions,
named constants, module variables that are part of the API, generic
interfaces, and `interface` blocks declaring `module subroutine` /
`module function` signatures. It contains no executable code.

A contract imports with `use ..., only:` and imports only what its
signatures need — a handful of types. That list is short by construction, not
by discipline.

**Rule 3 — A submodule is an implementation.**

All procedure bodies live in submodules. A submodule sees everything its
parent contract imported (host association), so it repeats no `use`
statements for those. It may `use` the umbrella of the layer below freely,
without `only:`, because nothing a submodule imports is ever re-exported and
no other file depends on a submodule for compilation.

Bodies are written as `module procedure name` ... `end procedure`, so the
signature exists in exactly one place, the contract. Nothing is repeated.

**Rule 4 — Umbrellas are for implementations and programs.**

`simple_core_module_api`, `simple_commanders_api` and the other `*_api`
modules exist so that implementations do not carry repeated `use` lists.
They are `use`d only from submodules and programs. A contract never `use`s
an umbrella. One umbrella per layer, re-exporting that layer and below.

**Rule 5 — One module per file, named alike.**

File `simple_foo.f90` holds `module simple_foo`; its implementations are
`simple_foo_<part>.f90` holding `submodule (simple_foo) simple_foo_<part>`.
(`simple_core_api.f90` holding `module simple_core_module_api` is the kind of
mismatch this rule ends.)

That is the whole system. No `private`/`public` bookkeeping is required:
because a contract holds nothing but API, default public is correct, and
because it imports only the few types in its signatures, what leaks through
re-export is negligible. Where a contract imports a type purely to name it in
a signature and must not re-export it, one line suffices:
`private :: image, oris`.

## 4. What each rule buys

| Rule | Incremental build | Clean build | Architecture |
|---|---|---|---|
| 2+3 contract/implementation | Body edit: 1 file recompiles. Contract edit: its submodules + direct users, no further, because users import only types. | All bodies are graph leaves; the critical path is the chain of small contracts plus one body. | The API of a subsystem is readable in one file; the implementation can be split by topic without changing the API. |
| 4 umbrellas only in implementations | No `.mod` carries an umbrella closure; cascades from utility changes stop at the contracts that name the changed type. | `.mod` files shrink ~10x; every compile reads far less. | Contracts document their real dependencies; umbrellas stay the convenience they were meant to be. |
| 1 layers | A change in `pftc` cannot recompile `image`. | Layer order is the build order. | Dependency direction is explicit and reviewable. |
| 5 naming | Tooling (module graph, indexes, dependency scanner) is exact. | — | — |

## 5. What it looks like

Contract, `main/nu_filt/simple_nu_filter.f90` (already close to this form):

```fortran
module simple_nu_filter
use simple_image, only: image
implicit none
private :: image
integer, parameter :: NU_DMAT_CANDIDATE_CAP = 24
type :: nu_highres_extension_stats
    ...
end type
interface
    module subroutine setup_nu_dmats( vol_even, vol_odd, n_highres_steps, evidence_source, fsc_res )
        class(image),               intent(in) :: vol_even, vol_odd
        integer,          optional, intent(in) :: n_highres_steps
        character(len=*), optional, intent(in) :: evidence_source
        real,             optional, intent(in) :: fsc_res
    end subroutine
    module real function get_nu_filter_bank_finest_lp()
    end function
end interface
end module simple_nu_filter
```

Implementation, `main/nu_filt/simple_nu_filter_bank.f90`:

```fortran
submodule (simple_nu_filter) simple_nu_filter_bank
use simple_core_module_api        ! umbrella: fine here, never in the contract
implicit none
contains
    module procedure setup_nu_dmats
        ! vol_even, vol_odd, n_highres_steps ... are declared once, in the contract
        ...
    end procedure
    module procedure get_nu_filter_bank_finest_lp
        ...
    end procedure
end submodule
```

A commander, `main/commanders/simple/simple_commanders_refine3D.f90`, becomes
a contract of ~40 lines (the commander types and their `exec` interfaces) plus
`simple_commanders_refine3D_auto.f90`, `simple_commanders_refine3D_states.f90`,
... as submodules that `use simple_commanders_api`. Editing `exec_refine3D_auto`
then recompiles one file instead of the 68 users of the commander API.

## 6. Migration

The transformation is mechanical and can be scripted (`scripts/`, alongside
`clean_simple_uses.pl`): for each module, lift every procedure signature
(header through the last dummy-argument declaration) into an `interface`
block in the contract, move the bodies into a submodule as
`module procedure`, move the module's `use` statements that only bodies need
into the submodule, and leave in the contract the `use ..., only:` of the
types named in signatures. The script reports what it could not classify; a
human resolves those.

Order, by payoff:

1. **Build system first** (one commit, no source change): tests under
   `EXCLUDE_FROM_ALL` behind a `tests` target; Ninja as the documented
   generator; `-j` bounded to the core count in `compile_*.sh`. This is
   immediate and independent.
2. **Workflow layer**: `main/commanders`, `main/stream`, `main/strategies`,
   `main/sieve`. These are the 68 umbrella users with the largest `.mod`
   files and the most frequent edits; they are also where bodies are long
   and contracts are trivially small (commander types + `exec`).
3. **Domain layer**: `flex`, `pftc`, `motion`, `ctf`, `pick`, `nano`.
   `nu_filt` is nearly done; convert its `module subroutine` bodies to
   `module procedure` to drop the duplicated declarations.
4. **Core layer**: `image`, `ori`, `parameters`, `project` already use
   submodules; audit them for Rule 4 and Rule 2 (no umbrella `use`, no bodies
   in the contract). `simple_builder` and `simple_cmdline` are the two to
   split.
5. **Utils**: leave as plain modules. They are small, stable and cheap.

Delivery follows the existing convention in this directory: small commits on
`master`, one subsystem at a time, each leaving the tree compiling and the
unit tests passing. A converted subsystem must not depend on an unconverted
one having been converted; the rules are per file and compose.

Two files are worth splitting regardless of the rules because they sit alone
on the critical path: `simple_commanders_test_highlevel.f90` (2.9k lines,
20 s) into per-topic submodules, and `simple_flex_pca_em_iter.f90` (2.9k
lines, 20 s, only 9 `use` lines, so its cost is intrinsic — profile it with
`gfortran -ftime-report` before deciding how).

## 7. Build system changes

- `production/CMakeLists.txt`: `add_executable(${test_exe} EXCLUDE_FROM_ALL ...)`
  and `add_custom_target(tests DEPENDS ${all_test_exes})`; keep `add_test`
  as is. The everyday `make`/`ninja` builds the library and six executables.
- Generator: `cmake -G Ninja` in `compile_*.sh` (CMake >= 3.20 handles Fortran
  module dependencies natively). Ninja schedules the critical path better
  than Make's recursive Fortran scanning, no-op builds are instant, and it
  does not oversubscribe the machine the way a bare `make -j` does.
- Flags: debug stays `-O0 -g -fbacktrace -fcheck=do,mem`; consider `-g1` for
  the everyday build. Add a `-O2` dev-release configuration for performance
  work; keep `-O3 -funroll-loops` for installs. `-fPIC` on a static library
  is inert and can go.
- `Fortran_MODULE_DIRECTORY` stays single; `.smod` files land there too.

## 8. Targets and how to check them

Add `scripts/build_profile.sh` (wraps the compiler with per-file timing, as
the 2026-09-14 profile did) and record before/after for:

- Clean build wall time on the reference Mac (`ninja` after `rm -rf build`).
- Touch-a-body rebuild: `touch` one commander submodule -> exactly one
  compile job.
- Touch-a-contract rebuild: `touch simple_oris.f90` -> its submodules and the
  contracts that name `oris`, nothing in `commanders/` beyond submodules.
- `du -sh build/modules`: target below 40 MB (from 251 MB).
- Largest `.mod`: target below 100 KB (from 371 KB).
- Critical path from the weighted module graph
  (`doc/code_overview/fortran-indexes/module_graph.dot` + per-file times):
  target below 15 s.

## 9. Gotchas

- `module procedure` bodies must match the contract exactly; the compiler
  enforces this, which is the point.
- A submodule inherits the contract's `use ..., only:` imports and may add
  its own. It does not inherit imports of sibling submodules; put shared
  implementation-only helpers in a parent submodule
  (`submodule (simple_foo) simple_foo_common`, then
  `submodule (simple_foo:simple_foo_common) simple_foo_part`).
- Module state that is implementation detail can live in a submodule's
  specification part and is shared with its descendants; keep it out of the
  contract unless it is API.
- Type-bound procedures: the type stays in the contract with
  `procedure :: foo => foo_impl`; `foo_impl` is declared in the contract's
  interface block and implemented in a submodule. Nothing changes for
  callers.
- Generic interfaces and operators stay in the contract.
- CMake's Makefile generator already avoids recompiling users when a `.mod`
  is byte-identical; Ninja does too. The cascades today are real interface
  changes amplified by re-export, which is exactly what Rule 4 removes.
- ccache cannot cache Fortran module side outputs; do not plan on it.

## 10. Out of scope

Behavioural changes of any kind, renaming public procedures, moving
directories between layers (only the dependency direction is enforced, the
present tree already nearly satisfies it), and the domain-model migration in
`staged_domain_driven_design_refactor.md`, which this proposal makes cheaper
but does not depend on.

## 11. Evaluation and rejection (2026-09-21)

All 56 modules in `src/main/commanders` were converted by
`scripts/split_contract.py` on the branch `contract-submodules` (worktree
`~/src/SIMPLE-contract`, head 72bbc0ac6). Each contract imported only
`cmdline` and `commander_base`; bodies moved to `<module>_impl.f90`. Master
and branch were built with `scripts/profile_build.sh` on the reference Mac
(arm64, GCC 16.2.0, Release, `BUILD_TESTS=ON`, Make generator, `make -j24`).
The generator was kept as Make by decision; section 7 was not evaluated.

| scenario | master | branch |
|---|---|---|
| clean build, wall | 142.0 / 142.9 / 142.3 s | 120.5 / 118.8 s (-16%) |
| clean build, compile CPU | 1338 s | 1396 s (+4%) |
| edit a commander body (`touch`) | 1 compile, 24.9 s | 1 compile (same) |
| interface change in `commanders_refine3D` (`probe`) | 7 files, 34.7 s | same 7 + contract, 24.8 / 25.6 s |
| interface change in `oris` (`probe`) | 643 files, 124.5 s | same 643 + 56 contracts, 106.3 / 101.4 s |
| `build/modules` | 265 MB | 237 MB (target was < 40 MB) |
| largest `.mod` | 685 KB | 556 KB (target was < 100 KB) |

Why the claims in sections 1, 4 and 5 did not hold:

- Body edits were already one compile. CMake's Makefile generator does not
  recompile users of a byte-identical `.mod`, and a body edit does not change
  it (section 9 notes this, but section 5's "one file instead of the 68
  users" contradicts it). About 10 s of the 25 s for such an edit is linking
  and regenerating `simple_ui_default_values`, not compiling.
- Rule 3 lets every implementation `use` the umbrella, so a core interface
  change recompiles every implementation: the same 643 files on both trees.
  506 of the 699 files recompiled on the branch never mention `ori`/`oris`;
  they recompile only because `simple_core_module_api` re-exports it.
- The measured gains are purely scheduling: users compile against a 0.2-0.5 s
  contract instead of waiting behind a 15 s body. Mean clean-build
  parallelism rose from 9.5 to 11.9 on 24 cores. A critical-path simulation
  (module graph weighted by measured compile times) over-predicted the
  measured gain by about 20% and puts a conversion of all of `src/main`
  (375 modules) at roughly 95-105 s clean, with body edits and core changes
  unchanged.
- The cost was visible immediately: 56 extra files, each commander's API and
  body in separate files, commander bodies 5% more compile CPU, and
  `scripts/default_audit.py` silently losing every commander's defaults until
  it was taught `module procedure`.

Cheaper levers found along the way, none requiring a source restructuring:
building without tests day to day (test programs compile in the last ~12 s
of the master build; `profile_build.sh clean --no-tests`, not yet measured);
the ~10 s link and UI-defaults step that follows every edit; and, if core
cascades matter, the umbrella re-export itself.
