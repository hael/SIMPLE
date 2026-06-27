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
