---
name: simple-main-params
description: Use when working in SIMPLE's src/main/params subsystem, including the parameters type, command-line parameter binding, the generated list of accepted keys, typed parsing, derived parameter phases, dynamic defaults, filt_mode flags, and validation of parsed workflow state.
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

`params%new(cline)` runs `init_dynamic_defaults`, reseeds the random generator,
then `parse_inputs`, `setup_execution_context`, `resolve_parameter_sources`,
`derive_sampling_settings`, `derive_parallel_settings`, `derive_image_settings`
and `validate_parameter_consistency`. A test's `set_fixed_seed` goes after it.

To add or change a command-line parameter:

1. Declare the field and declaration-time default in `simple_parameters.f90`.
   Its scalar declarations also become the accepted keys: the build generates
   `simple_args.f90` from them (`scripts/simple_args_generator.pl`), and
   `cmdline` rejects unknown keys.
2. Bind command-line parsing in `simple_parameters_parse.f90` with the matching
   registry call: `add_char`, `add_int`, `add_real`, `add_file`, or `add_dir`.
3. Derive dependent values in the right phase of `simple_parameters_phases.f90`.
4. Add mode restrictions and sanity checks in `validate_parameter_consistency`.
5. Touch `simple_parameters_core.f90` only for dynamic string defaults
   (`init_dynamic_defaults`) or shared lifecycle utilities.
6. Expose it through the owning UI program (`.github/skills/simple-main-ui/SKILL.md`).

## Working Rules

- After parsing, consume typed `params%...` fields in workflow logic.
- Use `cmdline` mainly before parsing, for command-shape validation, and for
  constructing sparse child command lines.
- Keep parser registration and validation close to the parameter's ownership.
- Check existing derived logicals before inventing new downstream tests. For
  example, `filt_mode` derives `l_lpauto`, `l_nonuniform`, and
  `l_nonuniform_lpset`.
- `parameters` is a large type (hundreds of string and character components):
  never a plain component of another type and never copied whole; copy the
  fields an object reads (`doc/policies/compile_time_policy.md`).
- No aliases for renamed or retired keys (release 4 keeps no backwards
  compatibility).
