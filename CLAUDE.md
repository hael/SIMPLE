# SIMPLE

SIMPLE is a scientific cryo-EM application platform: a large modern-Fortran
application with established ownership boundaries, generated metadata, and
workflow-specific orchestration layers. Treat it as an application platform,
not a narrow library.

## Use the repo skills

Skills live in `.claude/skills/` (a symlink to `.github/skills/`, shared with
Copilot). **Select the narrowest applicable skill before editing code.** Only
name + description stay in context; the skill body loads on demand.

- `simple-architecture`: read first when a task spans multiple subsystems.
- `simple-modern-fortran`: Fortran style, lifecycle, generated sources, modules.
- Workflow skills: `simple-abinitio2d`, `simple-refine3d`,
  `simple-abinitio3d-importance-sampling`, `simple-cluster-cavgs-quality`,
  `simple-microchunk-rejection`, `simple-frac-update-trailing`,
  `simple-nonuniform-regularization`.
- Subsystem skills: `simple-main-*` for `ui`, `root`, `commanders`,
  `strategies`, `project`, `ori`, `pftc`, `image`, `params`, `nu-filt`,
  `volume`, `ctf`, `motion`, `opt`, `pca`, `pick`, `star`, `stream`, `nano`,
  and related modules.

Routing:
- Multi-area task → `simple-architecture` first, then the most specific skill.
- refine3D / abinitio3D sampling or reconstruction → `simple-refine3d`, then the
  narrower sampling / fractional-update / nonuniform skill if the task touches
  those contracts.
- 2D workflow or class-average restoration → `simple-abinitio2d`.
- Streaming microchunk rejection / `model_cavgs_rejection` → read
  `simple-microchunk-rejection` before changing stream lifecycle or
  particle-state cleanup.

Do not guess ownership from filenames. Follow the flow:
`ui -> exec -> commander -> strategy/domain object`.

## Structure

- `production/`: thin executable entrypoints (`simple_exec`, `single_exec`,
  `simple_stream`, `simple_test_exec`, `simple_private_exec`).
- `src/`: core static library.
- `src/main/`: application and domain logic.
- `src/defs`, `src/fileio`, `src/utils`: shared infrastructure.
- `src/main/ui`: command/parameter metadata exposed to CLI/NICE.
- `src/main/params`: typed `parameters` object, parsing, derived settings, validation.
- `src/main/exec`: execution routers.
- `src/main/commanders`: high-level workflow command objects.
- `src/main/strategies`: algorithm and execution-policy layers.
- `src/main/nu_filt`: nonuniform filtering used by volume assembly.
- `doc/`: architecture, policy, and refactoring notes — often more current than
  code comments. Check `doc/policies/*.md` and `doc/refactoring_notes/*.md`
  before refactors.
- `nice/`: optional web/UI layer.

## Engineering defaults

- Prefer existing SIMPLE patterns over new abstractions.
- Keep changes scoped to the owning subsystem; avoid broad refactors while
  fixing a narrow behavior.
- Preserve dirty worktree changes that are not part of the task.
- Use `rg` / `rg --files` first when searching.
- Use direct, local Fortran APIs and structured project/parameter/orientation
  helpers instead of ad hoc parsing or parallel configuration paths.
- Do not run expensive builds or broad test suites unless requested. Prefer
  targeted checks and tell the user what was not run.

## Fortran conventions

- Use `use ..., only: ...` imports.
- Match local module, submodule, type-bound procedure, and lifecycle conventions.
- Preserve `new`/`kill` symmetry for stateful types.
- Keep orchestration in commanders/strategies and numerical work in domain modules.
- Reuse `parameters`, `cmdline`, `builder`, project, and orientation APIs.
- In commanders, normalize and validate `cmdline` before `params%new(cline)`.
- For CLI/UI-visible behavior, update the owning parameter, parser, UI metadata,
  commander, exec router, and project/reporting paths consistently.
- Argument metadata and git-hash sources may be generated during builds — check
  for generated sources before assuming a handwritten file is authoritative.

## Maintaining this config

Propose small, incremental updates to `.github/skills/` (and this file) after
significant tasks when a pattern or guardrail is worth capturing. Apply skill
edits only when explicitly requested or approved. Keep updates concise and
aligned with observed repository practice.
