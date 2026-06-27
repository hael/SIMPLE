# SIMPLE Copilot Instructions

SIMPLE is a scientific cryo-EM application platform. Treat this repository as a large
Fortran application with established ownership boundaries, generated metadata, and
workflow-specific orchestration layers.

## Use The Repo Skills

When working in Copilot agent mode, prefer the focused skills in `.github/skills/`.
Select the narrowest applicable skill before editing code:

- `simple-architecture`: broad repository orientation and cross-cutting workflows.
- `simple-modern-fortran`: Fortran style, lifecycle, generated-source, and module guidance.
- `simple-cluster2d`, `simple-refine3d`, `simple-cluster-cavgs-quality`: workflow-specific behavior.
- `simple-main-*`: subsystem guidance for `src/main` areas such as `ui`, `root`,
  `commanders`, `strategies`, `project`, `ori`, `pftc`, `image`, `volume`, and related modules.

If a task spans multiple areas, read `simple-architecture` first, then the most specific
subsystem skill. Do not guess ownership from filenames alone; follow the established
`ui -> exec -> commander -> strategy/domain object` flow.

## Engineering Defaults

- Prefer existing SIMPLE patterns over new abstractions.
- Keep changes scoped to the owning subsystem.
- Preserve dirty worktree changes that are not part of the task.
- Use `rg` or `rg --files` first when searching.
- Use direct, local Fortran APIs and structured project/parameter/orientation helpers
  instead of ad hoc parsing or parallel configuration paths.
- Avoid broad refactors while fixing a narrow behavior.
- Do not run expensive builds or broad test suites unless requested. Prefer targeted
  checks and tell the user what was not run.

## SIMPLE Structure

- `production/`: thin executable entrypoints.
- `src/`: core library.
- `src/main/`: application and domain logic.
- `src/defs`, `src/fileio`, `src/utils`: shared infrastructure.
- `src/main/ui`: command and parameter metadata exposed to CLI/NICE.
- `src/main/exec`: execution routers.
- `src/main/commanders`: high-level workflow command objects.
- `src/main/strategies`: algorithm and execution-policy layers.
- `doc/`: architecture, policy, and refactoring notes that may be more current than comments.

## Fortran Conventions

- Use `use ..., only: ...` imports.
- Match local module, submodule, type-bound procedure, and lifecycle conventions.
- Preserve `new`/`kill` symmetry for stateful types.
- Keep orchestration in commanders/strategies and numerical work in domain modules.
- Reuse `parameters`, `cmdline`, `builder`, project, and orientation APIs.
- In commanders, normalize and validate `cmdline` before `params%new(cline)`.
- Be aware that argument metadata and git-hash sources may be generated during builds.

