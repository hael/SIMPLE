---
name: simple-modern-fortran
description: Use when editing SIMPLE's Fortran code and you need repo-specific guidance on the modern Fortran style used here, including modules, submodules, type-bound procedures, abstract command interfaces, generated sources, OpenMP-aware code, and common patterns for extending the code safely.
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
- Preserve `new`/`kill` lifecycle symmetry for stateful types.
- Watch for generated sources from `scripts/simple_args_generator.pl` and git-hash insertion during builds.

## Debugging And Build Notes

- Debug build: `cmake .. -DCMAKE_BUILD_TYPE=Debug`
- Helpful notes: `doc/how2s/how2debug_fortran.txt`
- Many objects expose derived-type state through type-bound methods; prefer following those paths before adding temporary print code.

## When Extending Code

- Match the existing folder ownership before adding files.
- If a large module already uses submodules, extend the submodule structure instead of inflating the parent interface file.
- If a workflow already has `ui -> exec -> commander -> strategy`, add the feature in the narrowest layer that owns it.
- For parameter parsing changes, read `.github/skills/simple-main-params/SKILL.md`.
