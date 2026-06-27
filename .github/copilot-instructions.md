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
- `simple-2d-classification-restoration`: 2D probabilistic sampling, online class-average restoration, and even/odd conventions.
- `simple-bayesian-3d-refinement`: probabilistic/Bayesian-style 3D refinement, matcher I/O, assembly handoffs, automasking, and multi-state 3D behavior.
- `simple-abinitio3d-importance-sampling`: abinitio3D `update_frac`/`nsample*`, `sampled`/`updatecnt`, `prob_align` reuse, and trailing-reconstruction coupling.
- `simple-cartesian-frac-update-trailing`: reference contract for Cartesian fractional updates, online reconstruction I/O, previous halfmap/rho handoffs, and obsfield mirrors.
- `simple-nonuniform-regularization`: `filt_mode=nonuniform|nonuniform_lpset`, `nu_refine`, `_nu_filt`/`_nu_locres` products, automask/mask precedence, and `simple_nu_filter` lifecycle.
- `simple-main-*`: subsystem guidance for `src/main` areas such as `ui`, `root`,
  `commanders`, `strategies`, `project`, `ori`, `pftc`, `image`, `params`,
  `nu_filt`, `volume`, and related modules.

If a task spans multiple areas, read `simple-architecture` first, then the most specific
subsystem skill. Do not guess ownership from filenames alone; follow the established
`ui -> exec -> commander -> strategy/domain object` flow.

For refine3D or abinitio3D sampling/reconstruction questions, prefer
`simple-bayesian-3d-refinement` and then the narrower sampling,
fractional-update, or nonuniform skill when the task touches those contracts.
For 2D class-average restoration, read `simple-2d-classification-restoration`
after `simple-cluster2d`.

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

## Instruction And Skill Maintenance

- Propose small, incremental updates to `.github/instructions/` and `.github/skills/`
  after significant tasks when patterns or guardrails should be captured.
- Only apply edits to instruction or skill files when explicitly requested or approved.
- Keep instruction and skill updates concise, specific, and aligned with observed
  repository practice.

## SIMPLE Structure

- `production/`: thin executable entrypoints.
- `src/`: core library.
- `src/main/`: application and domain logic.
- `src/defs`, `src/fileio`, `src/utils`: shared infrastructure.
- `src/main/ui`: command and parameter metadata exposed to CLI/NICE.
- `src/main/params`: typed `parameters` object, parsing, derived settings, and validation.
- `src/main/exec`: execution routers.
- `src/main/commanders`: high-level workflow command objects.
- `src/main/strategies`: algorithm and execution-policy layers.
- `src/main/nu_filt`: nonuniform filtering implementation used by volume assembly.
- `doc/`: architecture, policy, and refactoring notes that may be more current than comments.

## Fortran Conventions

- Use `use ..., only: ...` imports.
- Match local module, submodule, type-bound procedure, and lifecycle conventions.
- Preserve `new`/`kill` symmetry for stateful types.
- Keep orchestration in commanders/strategies and numerical work in domain modules.
- Reuse `parameters`, `cmdline`, `builder`, project, and orientation APIs.
- In commanders, normalize and validate `cmdline` before `params%new(cline)`.
- Be aware that argument metadata and git-hash sources may be generated during builds.
