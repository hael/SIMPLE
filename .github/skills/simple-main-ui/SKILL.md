---
name: simple-main-ui
description: Use when working in SIMPLE's src/main/ui subsystem, including program registration, parameter metadata, CLI-visible program definitions, test UI groups, SINGLE/simple/stream UI groups, and JSON metadata generation for NICE and other interfaces.
---

# SIMPLE `src/main/ui`

This folder defines the program registry and parameter metadata.

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

## Working Rule

If a parameter or program is “missing from the surface,” the fix is often here rather than in commanders.
