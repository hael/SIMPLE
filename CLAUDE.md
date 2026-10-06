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
- Workflow skills: `simple-solve2d`, `simple-refine3d`,
  `simple-solve3d-importance-sampling`, `simple-cluster-cavgs-quality`,
  `simple-frac-update-trailing`, `simple-nonuniform-regularization`.
- Subsystem skills: `simple-main-*` for `ui`, `root`, `commanders`,
  `strategies`, `project`, `ori`, `pftc`, `image`, `params`, `nu-filt`,
  `volume`, `ctf`, `motion`, `opt`, `pca`, `pick`, `star`, `stream`, `nano`,
  and related modules.

Routing:
- Multi-area task → `simple-architecture` first, then the most specific skill.
- refine3D / solve3D sampling or reconstruction → `simple-refine3d`, then the
  narrower sampling / fractional-update / nonuniform skill if the task touches
  those contracts.
- 2D workflow or class-average restoration → `simple-solve2d`.
- Streaming particle sieving / sieve rejection → read `simple-main-stream`
  and `doc/policies/sieving_and_rejection/ptcl_sieve_policy.md` before
  changing sieve lifecycle or particle-state cleanup.
- Other stream pipeline changes → read `simple-main-stream` and the matching
  policy in `doc/policies/stream/` (reference generation, 3D ingestion, pool 2D,
  IPC, restart) before changing that contract.
- `model_cavgs_rejection` / class-average quality backend →
  `simple-cluster-cavgs-quality`.

Do not guess ownership from filenames. Follow the flow:
`ui -> exec -> commander -> strategy/domain object`.

`solve3D` is de novo map determination (ab initio 3D reconstruction coupled
with initial 3D refinement); `solve2D` is its 2D equivalent. Before 2026-10-03
they and their variants were named `abinitio*` (`src/main/abinitio/` is now
`src/main/solve/`), and `refine2D`, the 2D counterpart of `refine3D` whose
stages `solve2D` runs, was named `cluster2D`. The old names are not aliased
(release 4 keeps no backwards compatibility); use the new names in code and
living docs and leave the old ones in dated history docs.

## Structure

- `production/`: thin executable entrypoints (`simple_exec`, `single_exec`,
  `simple_stream`, `simple_test_exec`, `simple_private_exec`). `simple_test_exec`
  and the `*_tester` modules are built only with `BUILD_TESTS=ON`, the
  `compile_*.sh` default (`--exclude-tests` turns them off).
- `src/`: core static library.
- `src/main/`: application and domain logic.
- `src/defs`, `src/fileio`, `src/utils`: shared infrastructure.
- `src/main/ui`: command/parameter metadata exposed to CLI/NICE.
- `src/main/params`: typed `parameters` object, parsing, derived settings, validation.
- `src/main/exec`: execution routers.
- `src/main/commanders`: high-level workflow command objects; `commanders/stream` holds the
  stream pipeline's (the p00 master and p01-p07), driving the stage types of
  `src/main/stream/stages` (the master's parts are in `src/main/stream/master`, the modules
  the stages share in `src/main/stream/shared`, the 2D pool in
  `src/main/stream/pool2D`).
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
- Do not compile, link, run CMake builds, or execute test binaries unless the user
  explicitly requests it. The user performs compilation to avoid spending agent
  credits on builds.
- Validate edits with editor or language-server syntax diagnostics and other
  lightweight checks that do not compile code. Tell the user that compilation and
  runtime tests were left for them.

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

## Compile time

`doc/policies/compile_time_policy.md` is the rule set; follow it in every change.
- No large type (`parameters`, `cmdline`, a stage, a sieve) as a plain component:
  make it `allocatable`, allocated before its constructor and released in `kill`.
  Never copy a whole `parameters`; copy the fields you need.
- In testers, declare large objects `class(T), allocatable` and `allocate` them.
- Module-level imports only for what the module itself needs, with `only:`; no
  new umbrella modules or re-exports.
- No new CMake targets for subsets of sources, no per-file compile options for
  SIMPLE-owned code (fix warnings in code), no new global optimization flags
  without a run-time benchmark.
- Claim a compile-time gain only from a like-for-like `scripts/profile_build.sh`
  comparison.

## Maintaining this config

Propose small, incremental updates to `.github/skills/` (and this file) after
significant tasks when a pattern or guardrail is worth capturing. Apply skill
edits only when explicitly requested or approved. Keep updates concise and
aligned with observed repository practice.
