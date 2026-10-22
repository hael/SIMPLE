# SIMPLE Copilot Instructions

SIMPLE is a scientific cryo-EM application platform. Treat this repository as a large
Fortran application with established ownership boundaries, generated metadata, and
workflow-specific orchestration layers.

## Use The Repo Skills

When working in Copilot agent mode, prefer the focused skills in `.github/skills/`.
Select the narrowest applicable skill before editing code:

- `simple-architecture`: broad repository orientation and cross-cutting workflows.
- `simple-modern-fortran`: Fortran style, lifecycle, generated-source, and module guidance.
- `simple-solve2d`: solve2D orchestration, refine2D stages (the public program and its private worker), 2D probabilistic tables, class sums in the 2D matcher, class-average restoration and convergence in the 2D strategy, distributed class-average assembly, and even/odd conventions.
- `simple-refine3d`: probabilistic 3D refinement, the thin commanders and the strategy routines they wrap (volume assembly, `prob_align`, the sigma2 commit), the distributed master as project owner, matcher I/O and the particle buffer, automasking, nonuniform filtering, and multi-state 3D behavior.
- `simple-solve3d-importance-sampling`: solve3D `nsample`/`update_frac` and the full-sampling switch, stage-wise sampling, `prob_align` sampling reuse and the master's draw, `sampled`/`updatecnt`, and the realized update fraction trailing reconstruction consumes.
- `simple-frac-update-trailing`: reference contract for fractional updates, the particle buffer and cropped records feeding partial reconstruction, and the per-state trailing accumulator chain with its population blend.
- `simple-nonuniform-regularization`: `filt_mode=nonuniform|nonuniform_lpset`, `_nu_filt`/`_nu_locres` products, spherical NU support, the NU matching low-pass handoff, ML-regularized auxiliary coupling, and the `simple_nu_filter` lifecycle.
- `simple-cluster-cavgs-quality`: `model_cavgs_rejection` and the class-average quality backend.
- `simple-main-*`: subsystem guidance for `src/main` areas such as `ui`, `root`,
  `commanders`, `strategies`, `project`, `ori`, `pftc`, `image`, `params`,
  `nu_filt`, `volume`, and related modules.

If a task spans multiple areas, read `simple-architecture` first, then the most specific
subsystem skill. Do not guess ownership from filenames alone; follow the established
`ui -> exec -> commander -> strategy/domain object` flow.

For refine3D or solve3D sampling/reconstruction questions, prefer
`simple-refine3d` and then the narrower sampling, fractional-update, or
nonuniform skill when the task touches those contracts.
For 2D workflow or class-average restoration work, read `simple-solve2d`.
For streaming particle sieving or sieve rejection, read `simple-main-stream` and
`doc/policies/sieving_and_rejection/ptcl_sieve_policy.md` before changing sieve
lifecycle or particle-state cleanup behavior. For other stream pipeline changes,
read `simple-main-stream` and the matching policy in `doc/policies/stream/`
(reference generation, 3D ingestion, pool 2D, IPC, restart) before changing that
contract. For `model_cavgs_rejection` and the
class-average quality backend, read `simple-cluster-cavgs-quality`.

`solve3D` is de novo map determination (ab initio 3D reconstruction coupled
with initial 3D refinement); `solve2D` is its 2D equivalent. Before 2026-10-03
they and their variants were named `abinitio*` (`src/main/abinitio/` is now
`src/main/solve/`), and `refine2D`, the 2D counterpart of `refine3D` whose
stages `solve2D` runs, was named `cluster2D`. The old names are not aliased
(release 4 keeps no backwards compatibility); use the new names in code and
living docs and leave the old ones in dated history docs.

## Engineering Defaults

- Prefer existing SIMPLE patterns over new abstractions.
- Keep changes scoped to the owning subsystem.
- Preserve dirty worktree changes that are not part of the task.
- Use `rg` or `rg --files` first when searching.
- Use direct, local Fortran APIs and structured project/parameter/orientation helpers
  instead of ad hoc parsing or parallel configuration paths.
- Avoid broad refactors while fixing a narrow behavior.
- Compile and run tests as part of the development cycle. Build incrementally in the
  existing `build/` (`make -C build -j<n>`, then `scripts/run_fast_gate.sh build` and
  `make -C build install`); `./compile_clean.sh` and `./compile_debug.sh` rebuild from
  scratch. Run the touched area's `build/bin/simple_test_exec test=unit_<area>` and the
  high-level entries a change affects. Report what was built and run, and the results.
- Never commit or push. The maintainer decides what is committed and pushed: leave changes
  in the working tree and say which are ready.

- Tests follow `doc/policies/test_environment_policy.md`: a test can fail, lives in a
  `simple_<thing>_tester.f90` next to the code it tests, and runs through `simple_test_exec`;
  see the Tests section of `AGENTS.md`.

## Layering And Project I/O

`doc/policies/layering_and_project_ownership_policy.md` is the rule set; read it
before adding or changing a commander, a strategy routine, a worker step or any
code that reads or writes an `sp_project`. Its section 8 lists where the code does
not follow it yet; do not add to that list.

- Search for the step before writing it (`rg` the verb and the data it touches);
  extend an existing routine, never copy its frame.
- A part of a stage that another workflow or the master runs is a strategy routine on
  the caller's `(params, build)` with a `run_<x>(cline)` entry; its commander is a thin
  `execute(cline)` wrapper. A strategy never calls a commander for a part of a stage.
- One project owner per run: workers read the file and write part files, never the
  project; the master merges in memory and writes the changed segments once per
  iteration. Do not read what `build` already holds.
- Shared project state (even/odd halves, the sample, sigma2) has one maker and verifying
  consumers (`ensure_field_halves`).

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
- `src/main/commanders`: high-level workflow command objects; `commanders/stream` holds the
  stream pipeline's (the p00 master and p01-p07), driving the stage types of
  `src/main/stream/stages` (the master's parts are in `src/main/stream/master`, the modules
  the stages share in `src/main/stream/shared`, the 2D pool in
  `src/main/stream/pool2D`).
- `src/main/strategies`: algorithm and execution-policy layers. In `parallelization/` the
  refine2D/refine3D strategies and the distributed master's steps (volume and class-average
  assembly, probability alignment, the sigma2 commit) are routines on the master's
  `(params, build)`; their commanders are thin `execute(cline)` wrappers.
- `src/main/solve`: the solve2D/solve3D stage controllers.
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

## Compile Time

`doc/policies/compile_time_policy.md` is the rule set; follow it in every change.

- No large type (`parameters`, `cmdline`, a stage, a sieve) as a plain component: make it
  `allocatable`, allocated before its constructor and released in `kill`. Never copy a whole
  `parameters`; copy the fields you need.
- In testers, declare large objects `class(T), allocatable` and `allocate` them.
- Module-level imports only for what the module itself needs, with `only:`; no new umbrella
  modules or re-exports.
- No new CMake targets for subsets of sources, no per-file compile options for SIMPLE-owned
  code (fix warnings in code), no new global optimization flags without a run-time benchmark.
- Claim a compile-time gain only from a like-for-like `scripts/profile_build.sh` comparison.
