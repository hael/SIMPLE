---
name: simple-main-commanders
description: Use when working in SIMPLE's src/main/commanders subsystem, including the abstract commander_base interface (execute(cline) only) plus the SIMPLE, SINGLE, stream and test commander families that own high-level workflow behavior for most CLI programs, thin wrappers of strategy routines, and the high-level test commanders run by simple_test_exec.
---

# SIMPLE `src/main/commanders`

This is the command object layer.

## Read First

- `simple/simple_commander_base.f90`
- `simple/simple_commanders_refine3D.f90`
- `simple/simple_commanders_refine2D.f90`
- `simple/simple_commanders_preprocess.f90`
- `simple/simple_commanders_project_core.f90`

Then branch into `simple/`, `single/`, `stream/` (see `simple-main-stream`) or `test/` as needed.

## Role

- Owns high-level command behavior for named programs
- Bridges CLI/UI requests to strategies and domain objects
- Is usually the right place for workflow sequencing, defaults, and validation

## Working Rule

Keep commanders relatively thin: orchestration here, heavy algorithmic work in `strategies/` or domain folders.

- `commander_base` defers one procedure, `execute(self, cline)`. Do not add
  typed entry points beside it; a step another workflow runs on its own
  `(params, build)` is a strategy routine (`<x>_exec`) with a `run_<x>(cline)`
  entry, and the commander calls `run_<x>` (`commander_volassemble`,
  `commander_prob_align`(`_neigh`), `commander_prob_align2D`,
  `commander_cavgassemble`, `commander_commit_sigma2`).
- A worker step's frame (params and builder set-up, reproduction of the master's
  sample with `sample4update_reprod`, the empty-partition contribution, batching,
  completion) is written once; a new variant supplies only what differs. The
  three probability-table workers in `simple_commanders_prob.f90` each carry their
  own copy today; do not add a fourth. See the layering policy, section 3.
- `refine2D` and `refine3D` commanders drive a strategy object
  (`create_refine2D_strategy`, `create_refine3D_strategy`); the `*_distr_worker`
  types are the per-partition workers (`simple_private_exec prg=refine2D|refine3D`).
- Shared defaults that several commanders apply live in
  `src/defs/simple_default_clines.f90` (`set_refine2D_defaults`); the UI's
  display defaults are traced from them by `scripts/default_audit.py`.

## Parameter Construction Invariant

- Never create a `parameters` object twice in a commander execution path.
- Normalize and validate `cmdline` defaults before the single `params%new(cline)` call.
- After `params%new(cline)` creates/chdirs into an execution directory, set child command lines to `mkdir=no` before invoking nested commanders.

## Test Commanders

- `test/simple_commanders_test_class.f90` holds the unit-suite tables
  (`suites_<area>`, `add_suite`) behind `test=unit_<area>`, `lib_<area>` and `units`.
- High-level tests are commanders: `simple_commanders_test_highlevel.f90` (types
  and interfaces) with one submodule per family, plus `_single` and `_stream`.
  `_suite` runs `test=nightly`, `test=profile` and `test=identity` through the
  dispatcher (`simple_test_dispatch.f90`).
- `simple_test_registry.f90` lists every high-level CTest entry (command line,
  threads, CPU reservation, deadline, identity outputs); a new entry gets its
  registry row and its `production/CMakeLists.txt` registration together, and
  `scripts/check_test_registry.py` must pass. Rules:
  `doc/policies/test_environment_policy.md`.
