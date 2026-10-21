---
name: simple-main-ui
description: Use when working in SIMPLE's src/main/ui subsystem, including program registration, categories and visibility levels, parameter metadata, generated display defaults, private-program required keys, CLI-visible program definitions, test UI groups, SINGLE/simple/stream UI groups, and JSON metadata generation for NICE and other interfaces.
---

# SIMPLE `src/main/ui`

This folder defines the program registry and parameter metadata. The rules are
in `doc/policies/ui_layer_policy.md`.

## Read First

- `simple_ui.f90`
- `simple_ui_program.f90`
- `simple_ui_param.f90`
- `simple_ui_params_common.f90`
- `simple_ui_simple_group.f90`

Then read the specific group file under `simple/`, `single/`, or `simple_test/` that matches the program family.

## Role

- Defines available programs and their parameters
- Supports listing, help/describe flows, and JSON metadata emission
- Acts as the contract between CLI/UI and execution routing

## Registration

- Each program module (`simple/simple_ui_<cat>.f90`, `single/single_ui_<cat>.f90`)
  owns one `UI_CATEGORY` descriptor and a `construct_<cat>_programs` routine that
  the group module (`simple_ui_simple_group.f90` etc.) calls; `add_ui_program`
  rejects duplicate names and incomplete categories.
- Every program and program input has a visibility (`UI_VIS_STANDARD`,
  `UI_VIS_ADVANCED`, `UI_VIS_DEVELOPER`); required inputs are standard.
  Visibility never changes what the CLI accepts.
- Stacks and volumes are `UI_IMG`, every other path `UI_FILE`, scalars `UI_PARM`.
  Shared input metadata is defined once in `simple_ui_params_common.f90`.
- Programs are renamed without aliases.
- Private worker programs: required keys in `simple_private_prgs.f90`; a private
  program with a public program's name (`refine2D`, `refine3D`) is checked
  against the public descriptor (`cmdline%parse_private`).
- Test programs live in `simple_test/simple_test_ui_class|highlevel.f90`; the
  `suite=` help lists the sub-suites in table order and must match
  `scripts/check_test_registry.py`.

## Defaults

- The UI publishes defaults; it does not own them. Baselines come from
  `simple_parameters.f90` declarations and `init_dynamic_defaults`; exact
  commander overrides (`if( .not. cline%defined(k) ) call cline%set(k, v)`)
  are traced by `scripts/default_audit.py` into the generated
  `simple_ui_default_values.f90`. A cline helper that sets defaults must be
  listed there (`TRACED_DEFAULT_HELPERS` and the other tables) or its keys stay
  unresolved.
- Do not add literal UI defaults that duplicate a baseline.

## Checks

- `cmake --build build --target validate_stream_ui_contract` guards the
  `simple_stream` screens (`production/stream_ui_contract.json`); update the
  reference only after review (`update_stream_ui_contract`).
- `unit_ui` covers the UI JSON, hash and visibility contracts.

## Working Rule

If a parameter or program is “missing from the surface,” the fix is often here rather than in commanders.
