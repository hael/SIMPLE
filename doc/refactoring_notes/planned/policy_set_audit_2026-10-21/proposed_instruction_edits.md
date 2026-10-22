# Edits to the instruction files (applied 2026-10-21 with the maintainer's approval; text as proposed, minor wording changes)

## CLAUDE.md and .github/copilot-instructions.md: new section after "Engineering defaults"

## Layering and project I/O

`doc/policies/layering_and_project_ownership_policy.md` is the rule set; read it
before adding or changing a commander, a strategy routine, a worker step or any
code that reads or writes an `sp_project`.
- Search for the step before writing it (`rg` the verb and the data it touches);
  extend an existing routine, never copy its frame.
- A part of a stage that another workflow or the master runs is a strategy routine
  on the caller's `(params, build)` with a `run_<x>(cline)` entry; its commander is a
  thin `execute(cline)` wrapper. A strategy never calls a commander for a part of a stage.
- One project owner per run: workers read the file and write part files, never the
  project; the master merges in memory and writes the changed segments once per
  iteration. Do not read what `build` already holds.
- Shared project state (even/odd halves, the sample, sigma2) has one maker and
  verifying consumers (`ensure_field_halves`).

## .github/skills/simple-architecture/SKILL.md: add to "Runtime Layering"

Rules for layers, step kinds and project I/O: `doc/policies/layering_and_project_ownership_policy.md`.

## .github/skills/simple-main-commanders/SKILL.md: add to "Working Rule"

A worker step's frame (params/builder set-up, reproduction of the master's sample,
empty-partition contribution, batching, completion) is written once; a new
variant supplies only what differs. See the layering policy, section 3.

## .github/skills/simple-main-project/SKILL.md: add to "Persistence"

No file locking exists: one process writes a project file at a time. Workers do not
write the project. Shared state that must agree across processes has one maker (the
master) and verifiers (`ensure_field_halves`). See the layering policy, section 5.
