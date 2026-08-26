# SIMPLE Agent Instructions

These instructions apply to the entire repository.

## Guidance Locations

- `.github/copilot-instructions.md` - repository instructions and path-specific guidance.
- `.codex/AGENTS.md` - generated maps, workflow notes, and agent-facing references when deliberately force-added to Git.

## Initialize Before Editing

1. Inspect the current branch and worktree with:
   `git branch --show-current` and `git status --short --branch`.
2. Read `.github/copilot-instructions.md`.
3. Read every applicable file under `.github/instructions/`.
4. Read the smallest relevant set of repository skills under `.github/skills/`. For Fortran changes, always read `.github/skills/simple-modern-fortran/SKILL.md`.
5. Normalize the current branch name by replacing `/` with `-`, then read `.codex/branch-context/<normalized-branch>.md` when that file exists.
6. Read relevant workflow notes under `.codex/` before building or testing.

Branch-context files are local development memory and may describe work that is not present in another checkout. Verify important claims against current source and Git history.

## Repository Boundaries

- Keep orchestration in controllers, commanders, and strategies.
- Keep numerical algorithms in their owning domain modules.
- Reuse `parameters`, `cmdline`, `builder`, and existing lifecycle methods.
- Keep workflow changes aligned with the owning UI, execution, commander, strategy, and domain layers.

## Development Notes

- When work needs a pre-implementation design record, create one concise,
  living note that contains the contract, implementation plan, review findings,
  and validation criteria.
- Do not create parallel specification and plan documents for the same
  development unless the user explicitly requests separate documents.
- Update the single note as the design is reviewed and implemented; rely on Git
  history instead of retaining redundant companion notes.

## Git and Validation

- Preserve unrelated user changes.
- Do not compile, link, run CMake builds, or execute test binaries unless the user explicitly requests it. The user performs compilation to avoid spending agent credits on builds.
- Validate edits with editor or language-server syntax diagnostics and other lightweight checks that do not compile code. Report that compilation and runtime tests were left to the user.
- For this project, pushing a focused, validated commit to `origin/master` is the normal XP workflow unless the user gives a different direction. Never force-push.
- Create commits only when the user asks, using focused messages.
- Do not claim Linux/BOX tests passed unless their output was actually observed.
- Record verified branch milestones and outstanding tests in the matching `.codex/branch-context/` note.
