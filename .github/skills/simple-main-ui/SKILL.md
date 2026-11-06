---
name: simple-main-ui
description: Use when working in SIMPLE's src/main/ui subsystem, including program registration, categories and visibility levels, parameter metadata, generated display defaults (simple_ui_default_values.f90 from scripts/default_audit.py), private-program required keys, CLI-visible program definitions, requirement groups, test UI groups, SINGLE/simple/stream UI groups, JSON metadata for NICE (print_ui_json, print_ui_stream) and other interfaces, and the chain for adding a program or exposing a key.
---

# `src/main/ui`: the program registry and its descriptors

Owns the `ui_program` descriptors of every program (name, executable, project requirement, inputs,
visibility, requirement groups), the shared input metadata, the private programs' required keys and
the JSON written from them. Not here: the accepted keys and their defaults (`params/`; the UI
publishes defaults), routing (`exec/`), behaviour. Rule set: `doc/policies/ui_layer_policy.md`.

## Rules

Silent ones first: nothing reports these.

- A program's inputs are its help, GUI and required keys, not the keys it accepts: `cmdline` accepts
  every `parameters` field for every program. So no input, and no parameter, for an internal
  mechanism or a test switch (every field is a key of every program): decide from existing state,
  force alternative paths in A/B identity runs with an existing key, ask Hans before adding a
  parameter. An input whose key the program does not read is accepted and ignored, and nothing checks
  that (`trajectory_graphene_subtr` offered `stk_traj`/`stk_backgr` and read `stk`/`stk2`, fixed
  2026-11-05); `scripts/check_cli_chain.py` fails the gate only on an input that is no accepted key.
- Exception: `solve3D_addon` and `solve3D_chunks` refuse keys their descriptor lacks
  (`ui_program%accepts`), and `solve3D_chunks` forwards to the add-on only what `solve3D_addon`
  accepts, dropping the rest without a word.
- `ui_program%new` arguments are behaviour: `executable` decides which binary may run the program;
  `sp_required` makes `params%new` find the project in the cwd and, under `mkdir=yes`, copy it into
  the execution directory, so the run writes the copy.
- Required inputs and requirement groups are enforced by `cmdline`, also on workers that share a
  public program's name (`refine2D`, `refine3D`, `reconstruct3D`, ...): a new required input stops
  the command lines their masters write. Every other worker of `dispatch_private_prg` has an entry in
  `simple_private_prgs.f90` (`add_private_prg`, then `push_req`/`push_opt`; the count follows), which
  `scripts/check_cli_chain.py` enforces; a required key there stops the masters' command lines that
  lack it.
- Defaults are generated (`simple_ui_default_values.f90` in the build tree, by
  `scripts/default_audit.py`): the declaration, `init_dynamic_defaults`, and a literal
  `if( .not. cline%defined(k) ) call cline%set(k, v)` in the routed commander's `execute` before
  `params%new`, the front door, the router case or a helper named in `TRACED_DEFAULT_HELPERS` or
  `DELEGATED_WORKFLOWS`. Anywhere else (after `new`, `strategies/`, `solve/`) the UI shows the
  baseline instead. The rule: a default the UI must show is set in the commander before
  `params%new`. No literal UI default that repeats a baseline.
- `prg=print_ui_stream` (`print_stream_ui_json`) is written by hand, keys, defaults and program names
  included: a stream key or default change is made there too. Every build checks it against the
  descriptors (`validate_complete_ui_json.py`: its processes are simple_stream programs, each user input
  it gives a process is defined, is an input of that program and has its default), the gate checks its
  keys are accepted, and `validate_stream_ui_contract` compares it with
  `production/stream_launcher_contract.json`.
- Visibility is presentation only. Required inputs are Standard (anything else stops); optional
  inputs default to Advanced, programs to Developer.
- Loud: `make_ui` runs in every front door and stops every program on a bad descriptor (summary
  30-100 characters, display name, placeholder form, a default outside its choices, a duplicate
  name, incomplete category). Testers pin some input sets (`simple_ui_visibility_tester.f90`;
  `inpl_cont` in `strategies/search/simple_strategy3D_inplane_tester.f90`). A test program's
  `suite=` help lists its sub-suites in table order (`scripts/check_test_registry.py`).

## Adding a program

Descriptor in its category module (`simple/simple_ui_<cat>.f90`, `single/single_ui_<cat>.f90`,
`simple_test/`) with the module's `UI_CATEGORY`, shared inputs from `simple_ui_params_common.f90`;
then a route (skill `simple-main-exec`) and a commander. Renames get no alias. Every step and what
fails silently when it is missed: `references/cli_change_chain.md`.

## Where to start

`simple_ui.f90` (`make_ui`, lookups, JSON writers), `simple_ui_program.f90` (descriptor,
`add_input`, required keys, `accepts`), `simple_ui_param.f90` (one input, choices, placeholders).

## Tests

- `build/bin/simple_test_exec test=unit_ui` (`suite=` `ui_json`, `gui_metadata`, `gui_assembler`,
  `ui_hash`, `ui_visibility`); under a second.
- Every build regenerates and validates the whole UI JSON (`validate_complete_ui_json`).
- simple_stream programs and the launcher: `cmake --build build --target validate_stream_ui_contract`;
  update the references (`update_stream_ui_contract`) only after review. The programs' reference
  (`production/stream_ui_contract.json`) differed from the descriptors already at `6565f5ba5`: its
  review is the maintainer's.
- The chain: `python3 scripts/check_cli_chain.py` (in the fast gate).
- Defaults by hand: `cmake --build build --target audit_ui_defaults`, then read
  `build/production/simple_ui_default_audit.md` (unresolved and unexposed defaults).

## See also

`doc/policies/ui_layer_policy.md`, `doc/for_developers/onboarding/gui_onboarding.md`; skills
`simple-main-params`, `simple-main-exec`.
