---
name: simple-modern-fortran
description: Use when editing SIMPLE's Fortran code and you need repo-specific guidance on the modern Fortran style used here, including modules, submodules, type-bound procedures, abstract command interfaces, generated sources, THROW macro and no-short-circuit traps, OpenMP parallel-region rules, compile-time rules (allocatable heavy components, lean imports, no per-file flags), comment rules, and common patterns for extending the code safely.
---

# Modern Fortran In SIMPLE

SIMPLE already uses modern Fortran heavily. Follow the local style instead of introducing foreign patterns.

## Patterns To Expect

- Module-based APIs with `use ..., only: ...`
- Abstract base types with deferred procedures, e.g. command interfaces
- Type-bound procedures and generic bindings on domain objects
- Submodules for large types such as `image`, `oris`, `sp_project`, and `polarft_calc`
- Allocatable state and explicit lifecycle methods such as `new`, `copy`, `kill`
- OpenMP-aware loops and shared-memory execution paths
- C interoperability for FFTW and selected low-level helpers

## Best Entry Files

- `src/main/simple_builder.f90`
- `src/main/simple_cmdline.f90`
- `src/main/image/simple_image.f90`
- `src/main/project/simple_sp_project.f90`
- `src/main/commanders/simple/simple_commander_base.f90`

## Local Conventions

- Prefer adding behavior through existing types and submodules rather than creating detached helper files.
- Keep orchestration in commanders/strategies and numerical work in domain modules.
- Reuse `parameters`, `cmdline`, and `builder` instead of inventing parallel configuration plumbing.
- In commanders, normalize and validate `cmdline` before the single `params%new(cline)` call; never construct a `parameters` object twice in one commander execution path.
- After `params%new(cline)` or a builder parameter initializer, treat the typed
  `parameters` object as authoritative. Use `params%...` fields in workflow
  logic instead of re-reading parsed values with `cmdline%get_iarg`,
  `get_rarg`, or `get_carg`.
- New SIMPLE command-line arguments must be represented in
  `src/main/params/simple_parameters.f90` and registered through the normal
  parsing/checking path before downstream code consumes them.
- For project-file workflows, do not pass native `box`/`smpd` through child
  command lines. Project metadata is the native sampling authority; pass
  `box_crop` only when a deliberate crop is requested and let parsing derive
  `smpd_crop`.
- Keep child command lines sparse. Set only values that change behavior or
  preserve execution mode rather than restating defaults already supplied by the
  child commander or `parameters`.
- Do not declare a dummy argument `optional` when every valid caller supplies
  it; make the dependency explicit and update all callers together.
- Preserve `new`/`kill` lifecycle symmetry for stateful types.
- A new stateful module is an encapsulated class: `private` components,
  `new`/`kill`, behavior as type-bound procedures, private helpers, getters
  only for what production code reads, and no unrelated utilities exported
  beside the type, which is named after its module (`simple_<x>` holds `<x>`).
  Do not model it as a public record passed to free procedures, even where an
  older module does (`simple_sigma2_state`); free-procedure modules are for
  stateless utilities.
- Generated at build time, never edited by hand: `simple_args.f90` (the accepted
  keys, from the declarations in `simple_parameters.f90`, by
  `scripts/simple_args_generator.pl`), `simple_ui_default_values.f90`
  (`scripts/default_audit.py`) and the git hash (`scripts/insert_git_commit_hash.pl`).
- Layout: in a run of `use mod, only:` lines align `only:`; continue a long list
  on the next line at statement indent + 4 with `&`; align a procedure's dummy
  declarations so `::` lines up.
- Fortran does not short-circuit `.and.`/`.or.`. Never put a guard and the
  access it protects in one expression (`flag .and. arr(i)`,
  `allocated(a) .and. a(i) > 0`); nest the `if`s. The bounds-checked Debug
  build (`./compile_debug.sh`) catches many of these; Release does not.
- No `goto` or numbered-label control flow; restructure instead of copying
  legacy routines that use it (`simple_opt_lbfgsb.f90`).
- Do not rely on integer overflow (wrap-around hashes or random generators): it
  is undefined and the Release build is `-O3`. Use 64-bit integers and explicit
  `modulo`.
- Delete dead code and state nothing reads; do not keep it for safety.
- Release 4 keeps no backwards compatibility: no aliases, shims or old-format
  readers. The one exception is reading `.simple` project files of earlier
  releases.
- Do not pass the pointer result of a polymorphic type-bound function straight
  as an actual argument (`call json%add(p, item%jsonise())` with `item` a
  `class(...)` dummy): assign it to a local pointer first. gfortran 15.2
  corrupted the json-fortran tree when the result was passed straight on.
- Keep new function-like preprocessor macro invocations on one physical line.
  `THROW_HARD(...)` and `THROW_WARN(...)` are expanded before Fortran parses the
  source, so a Fortran continuation inside their argument is flattened. The
  common `//&` form then becomes an invalid expression:

  ```fortran
  ! Wrong: fails after preprocessing.
  THROW_HARD('first part '//&
      &'second part')

  ! Right: shorten the message.
  THROW_HARD('concise error message')
  ```

  If the full text is important, build it with ordinary Fortran continuation in
  a named character variable or parameter, then pass only that identifier to
  the macro: `THROW_HARD(error_message)`. Do not copy legacy multiline macro
  calls: some compile only because preprocessing flattens them, which can leave
  continuation markers in the emitted diagnostic text.
- Before handing off Fortran edits, scan fatal forms with
  `rg -n "THROW_(HARD|WARN)\\(.*//&" src` and inspect newly added macro lines
  with `git diff -U0 -- '*.f90' | rg '^\\+.*THROW_(HARD|WARN)\\(.*&[[:space:]]*$'`.

## Compile Time

The rule set is `doc/policies/compile_time_policy.md`. gfortran expands the
initialization, copy and clean-up of every plain (non-allocatable) component tree
inline, at each declaration and assignment, so large inline types cost compile time
everywhere they are declared.

- A component of a large type (`parameters`, `cmdline`, a stream stage, `ptcl_sieve`,
  any aggregate with many allocatable or string fields) is `allocatable`. Allocate
  it right before its constructor (`if( .not. allocated(self%x) ) allocate(self%x)`)
  and deallocate it in `kill` on every path, including the early return of a
  `kill` that checks `l_exists`.
- Never store a copy of a whole `parameters` (`self%params = params`). Copy the
  fields the object reads; delete a copy nothing reads.
- In testers, a local stage, sieve or similar object is
  `class(T), allocatable :: x` followed by `allocate(x)`; helper dummies are
  `class(T)`.
- Module variables used only by one submodule are declared in that submodule;
  state shared by several submodules goes into a small companion module the parent
  imports (`simple_nu_filter_vars`), not into the parent.
- Module-level `use` lines are build-order edges for every importer: import with
  `only:` and only what the module itself needs; put procedure-only imports in the
  procedure or its submodule. Do not add umbrella modules or re-exports.
- With a third-party module, name the entities you use (`use FoX_dom, only: ...`).
- Fix warnings in code; SIMPLE-owned sources take no per-file compile options and
  no new CMake targets.

## Comments

Comments are short, snappy and accurate; long ones have been found to rot.
See `doc/refactoring_notes/completed/comment_novel_inventory_2026-09-30.md`.

- File header: one `!@descr:` line first (`scripts/check_descr.py`), then at
  most 5 lines.
- Procedure header: at most 3 lines. In-body comments: 1-2 lines that say why
  or state a non-obvious invariant, not what the next line does.
- No dates, names, commit hashes, run measurements, "used to" narration,
  rejected alternatives, `file.f90:NNN` references or out-of-tree paths in
  code. Git history keeps them; decision logs live in `doc/`.
- Do not restate a policy note or list what the compiler knows (dependencies,
  routine lists, field lists); point to the doc instead.
- When code changes, fix or delete the comments it falsifies. A comment that
  names a switch or path must name one that exists.

## OpenMP

One thread runs a parallel loop serially, so one-thread identity never shows
thread safety. For every parallel region you add or change:

- No string or allocatable-character temporaries and no calls into routines
  with module-level buffers inside the region; build file names and messages
  before it.
- Size per-thread work arrays by the team the region actually runs
  (`omp_get_max_threads()` read just before it), not by `nthr_glob` or
  `params%nthr`.
- An error inside the region sets a per-task status; raise it after the region.
- Add a tester check that forces a team and compares with the serial result
  (the fast gate runs at one thread, teams capped at 4). Stress new parallel
  code by hand at the machine's full thread count, repeated, bit for bit where
  the result is deterministic.

## Debugging And Build Notes

- Debug build: `./compile_debug.sh` (runtime checks including bounds; `--traps`
  for floating-point traps); `./compile_clean.sh` for Release. Both run the
  fast gate unless `--exclude-tests`.
- Helpful notes: `doc/how2s/how2debug_fortran.md`
- Many objects expose derived-type state through type-bound methods; prefer following those paths before adding temporary print code.

## When Extending Code

- Match the existing folder ownership before adding files.
- If a large module already uses submodules, extend the submodule structure instead of inflating the parent interface file.
- If a workflow already has `ui -> exec -> commander -> strategy`, add the feature in the narrowest layer that owns it.
- For parameter parsing changes, read `.github/skills/simple-main-params/SKILL.md`.
