---
name: simple-main-commanders
description: Use when working in SIMPLE's src/main/commanders subsystem, including abstract command interfaces plus the SIMPLE, SINGLE, and test commander families that own high-level workflow behavior for most CLI programs.
---

# SIMPLE `src/main/commanders`

This is the command object layer.

## Read First

- `simple/simple_commander_base.f90`
- `simple/simple_commanders_refine3D.f90`
- `simple/simple_commanders_cluster2D.f90`
- `simple/simple_commanders_preprocess.f90`
- `simple/simple_commanders_project_core.f90`

Then branch into `simple/`, `single/`, or `test/` as needed.

## Role

- Owns high-level command behavior for named programs
- Bridges CLI/UI requests to strategies and domain objects
- Is usually the right place for workflow sequencing, defaults, and validation

## Working Rule

Keep commanders relatively thin: orchestration here, heavy algorithmic work in `strategies/` or domain folders.

## Parameter Construction Invariant

- Never create a `parameters` object twice in a commander execution path.
- Normalize and validate `cmdline` defaults before the single `params%new(cline)` call.
- After `params%new(cline)` creates/chdirs into an execution directory, set child command lines to `mkdir=no` before invoking nested commanders.
