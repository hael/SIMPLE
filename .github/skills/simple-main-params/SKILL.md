---
name: simple-main-params
description: Use when working in SIMPLE's src/main/params subsystem, including the parameters type, command-line parameter binding, typed parsing, derived parameter phases, dynamic defaults, filt_mode flags, and validation of parsed workflow state.
---

# SIMPLE `src/main/params`

This folder owns the typed `parameters` object and the parse/derive/validate
pipeline behind `params%new(cline)` and builder parameter initialization.

## Read First

- `simple_parameters.f90`
- `simple_parameters_parse.f90`
- `simple_parameters_phases.f90`
- `simple_parameters_core.f90`
- `simple_parameters_registry.f90`

## Parameter Lifecycle

To add or change a command-line parameter:

1. Declare the field and declaration-time default in `simple_parameters.f90`.
2. Bind command-line parsing in `simple_parameters_parse.f90` with the matching
   registry call, such as `add_char`, `add_int`, `add_real`, `add_file`, or `add_dir`.
3. Derive dependent values in the right phase of `simple_parameters_phases.f90`.
4. Add mode restrictions and sanity checks in `validate_parameter_consistency`.
5. Touch `simple_parameters_core.f90` only for dynamic string defaults or shared
   lifecycle utilities.

## Working Rules

- After parsing, consume typed `params%...` fields in workflow logic.
- Use `cmdline` mainly before parsing, for command-shape validation, and for
  constructing sparse child command lines.
- Keep parser registration and validation close to the parameter's ownership.
- Check existing derived logicals before inventing new downstream tests. For
  example, `filt_mode` currently derives `l_lpauto`, `l_nonuniform`, and
  `l_nonuniform_lpset`.
