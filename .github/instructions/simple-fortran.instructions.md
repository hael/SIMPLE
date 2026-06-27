---
applyTo: "src/**/*.f90,production/**/*.f90"
---

# SIMPLE Fortran Editing

Use the `.github/skills/simple-modern-fortran/SKILL.md` and the relevant
`.github/skills/simple-main-*/SKILL.md` file before making nontrivial changes.

Prefer small, owner-aligned edits. Follow existing module imports, type-bound methods,
allocatable state patterns, OpenMP-aware loops, and `new`/`kill` lifecycle conventions.

For command-line or UI-visible behavior, update the owning parameter, parser, UI metadata,
commander, execution-router, and project/reporting paths consistently. Check for generated
argument sources before assuming a handwritten file is authoritative.

