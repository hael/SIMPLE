# UI Descriptor Layer Policy

## 1. Scope

This document describes the current command-descriptor system in
`src/main/ui/` and defines the rules that changes to that system must follow.
The descriptor system is shared by the command-line interface and clients
that consume `simple_ui.json`.

The planned structural and schema changes are recorded separately in
[`../refactoring_notes/ui_descriptor_layer_refactoring.md`](../refactoring_notes/ui_descriptor_layer_refactoring.md).

This document does not cover:

- scientific parameter parsing in `src/main/params/`;
- commander dispatch or workflow execution;
- the stream runtime-status metadata in `src/utils/gui/metadata/`;
- `gui_assembler`, which sends live stream results and status to NICE.

The runtime-status classes and the command-descriptor classes both serve GUI
clients, but they are separate systems with different data and lifetimes.

## 2. What the UI descriptor layer does

The descriptor layer supplies five services:

1. It registers the programs supported by `simple_exec`, `single_exec`,
   `simple_stream`, and `simple_test_exec`.
2. It records the input keys associated with each program.
3. It tells command-line parsing which program/executable combinations are
   valid and which input keys are required.
4. It prints program help and command-line input summaries.
5. It exports program and input metadata as JSON for external UI clients.

It describes commands; it does not execute them. The commander's parameter
and validation paths remain authoritative for scientific behavior.

## 3. Current object model

The core model uses composition. `ui_program` does not extend `ui_param`, and
`ui_param` does not extend another UI class. A `ui_program` contains linked
lists whose elements are copies of `ui_param` values.

`ui_hash` is the only inheritance relationship in the descriptor layer:
it extends the generic `vrefhash` and adds typed accessors for `ui_program`
and `ui_param`.

```text
vrefhash
   |
   +-- ui_hash
         |
         +-- references module-owned ui_program objects by program name

ui_program
   |
   +-- img_ios      -- linked_list of ui_param values
   +-- parm_ios     -- linked_list of ui_param values
   +-- alt_ios      -- linked_list of ui_param values
   +-- srch_ctrls   -- linked_list of ui_param values
   +-- filt_ctrls   -- linked_list of ui_param values
   +-- mask_ctrls   -- linked_list of ui_param values
   +-- comp_ctrls   -- linked_list of ui_param values
```

### 3.1 `ui_param`

`ui_param` is declared in `src/main/ui/simple_ui_param.f90`. It describes one
CLI input.

| Field | Current meaning |
| --- | --- |
| `key` | Stable CLI key, for example `nthr` or `mskdiam`. |
| `keytype` | String type tag used by UI consumers, such as `num`, `file`, `dir`, `str`, `binary`, or `multi`. |
| `descr_short` | Short user-facing description. Usually rendered as a field label. |
| `descr_long` | Longer user-facing help text. |
| `descr_placeholder` | Placeholder text. It also currently encodes choices and defaults for some inputs. |
| `gui_submenu` | Optional GUI group name within a program. |
| `active_flags` | Optional string condition controlling when the input is active. |
| `exclusive_group` | Optional identifier for mutually exclusive inputs. |
| `cval_default` | Character default for non-numeric inputs. |
| `rval_default` | Real default for numeric inputs. |
| `required` | Whether the CLI requires the key. |
| `advanced` | Current two-level GUI visibility flag. |
| `online` | Marks an input as applicable to online operation. |

`set_param` has numeric-default and character-default overloads.
`apply_gui_overrides` applies program-specific GUI grouping, activation,
visibility, and online overrides after a common parameter has been copied.

The current setters only retain a supplied default when `required` is false.
Callers still pass a dummy default for required inputs because of the
constructor signature.

### 3.2 `ui_program`

`ui_program` is declared in `src/main/ui/simple_ui_program.f90`. It describes
one value accepted by the `prg=` or `test=` argument.

| Field | Current meaning |
| --- | --- |
| `name` | Stable program key, for example `refine3D`. |
| `descr_short` | Short program description used in listings and JSON. |
| `descr_long` | Long description printed by `describe=yes` and exported to JSON. |
| `executable` | Owning executable, or `all` where supported. |
| `gui_submenu_list` | Comma-separated GUI group names. |
| `advanced` | Current two-level program visibility flag. |
| seven linked lists | Inputs divided into image, parameter, alternative, search, filter, mask, and computer sections. |
| `sp_required` | Whether normal execution requires a SIMPLE project file. |
| `exists` | Internal lifecycle guard used by `new` and `kill`. |

`new` initializes the program. `add_input` has three overloads:

- create a numeric input from individual arguments;
- create a character-valued input from individual arguments;
- copy a common `ui_param` and apply program-specific overrides.

The section selector passed to `add_input` is one of `UI_IMG`, `UI_PARM`,
`UI_ALT`, `UI_SRCH`, `UI_FILT`, `UI_MASK`, or `UI_COMP`.

`ui_program` also:

- prints detailed UI metadata;
- prints CLI help grouped by section;
- prints its long description;
- writes a per-program JSON document;
- returns the program name and executable;
- counts and returns required CLI keys;
- reports whether a project file is required.

### 3.3 Common parameter descriptors

`src/main/ui/simple_ui_params_common.f90` owns reusable `ui_param` objects.
`set_ui_params` initializes them before programs are constructed. Program
modules copy these objects with `add_input` and may override descriptions,
requiredness, grouping, activation, visibility, and online behavior.

Common descriptors avoid repeating the base metadata for widely used keys.
An override belongs in a program module only when that program presents or
uses the input differently.

### 3.4 Program modules and groups

Program objects are module-owned `type(ui_program), target` variables. Domain
modules such as `simple_ui_preproc` and `simple_ui_refine3D` initialize those
objects and add them to a registry.

The group modules aggregate construction:

| Module | Registry contents |
| --- | --- |
| `simple_ui_simple_group` | Normal SIMPLE programs. |
| `simple_ui_stream_group` | SIMPLE stream programs. |
| `simple_ui_single_group` | SINGLE programs. |
| `simple_ui_test_group` | Test programs in a separate test registry. |

`simple_ui_utils::add_ui_program` rejects duplicate program keys and stores a
reference to the module-owned program object in `ui_hash`.

### 3.5 Registry and construction flow

`simple_ui` owns the production registry `prgtab`, the test registry `tsttab`,
and their sorted key lists.

```text
executable
   |
   +-- make_ui / make_test_ui
          |
          +-- set_ui_params
          |
          +-- group constructors
                 |
                 +-- domain program constructors
                        |
                        +-- ui_program%new
                        +-- ui_program%add_input
                        +-- add_ui_program -> ui_hash
```

`make_ui` registers SIMPLE, stream, and SINGLE programs in the same production
registry. `make_test_ui` builds the separate test registry.

`get_prg_ptr` and `get_test_prg_ptr` return pointers to registered
`ui_program` objects.

### 3.6 CLI consumers

The executables call `make_ui` before command parsing. `simple_cmdline` uses
the selected `ui_program` to:

- reject unknown program names;
- check that the program belongs to the selected executable;
- print the long description for `describe=yes`;
- obtain the statically required keys;
- print grouped command help when required inputs are missing.

The parameters layer also retains the selected `ui_program` pointer and reads
`sp_required` while preparing the execution context.

The UI descriptor is therefore part of the CLI contract. Removing a program,
changing its executable, changing `required`, or moving a key out of a program
can change command-line behavior.

### 3.7 JSON consumers

`simple_ui` currently provides:

- `print_ui_json`, which prints the full registry;
- `write_ui_json`, which writes `simple_ui.json`;
- `validate_ui_json`, which constructs, writes, and parses that file;
- `print_stream_ui_json`, a separate hand-built stream workflow description.

`ui_program%write2json` writes a single program to its own JSON file.

The full-registry print and write routines contain separate copies of the
program and parameter serialization logic. The per-program writer contains a
third closely related implementation. For `binary` and `multi` inputs these
serializers currently derive the option list by parsing text between
parentheses in `descr_placeholder`.

## 4. Current descriptor declaration pattern

New code normally reuses a common input:

```fortran
call refine3D%add_input(UI_MASK, mskdiam, &
    gui_submenu='mask', gui_advanced=.false.)
```

It declares an input inline when the key is program-specific:

```fortran
call automask%add_input(UI_IMG, 'vol1', 'file', &
    'Odd volume', 'Odd volume', 'e.g. vol1.mrc', .true., '')
```

The first form copies the common parameter value. The second constructs a new
temporary value. In both cases `add_input` stores a copy in the selected
program list.

## 5. Policy for the existing contract

### 5.1 One CLI/GUI descriptor

- CLI keys, accepted values, defaults, requiredness, dependencies, and
  execution ownership must not be redefined in a GUI client.
- A GUI command must pass through the normal CLI/parameter validation path.
- A supported CLI command must remain representable by the structured UI
  model, including Developer controls.
- User-facing text may change without changing stable CLI keys.

### 5.2 Visibility

Every program and parameter must be classified as:

| Level | Definition |
| --- | --- |
| `standard` | Needed or normally inspected for a typical successful run. |
| `advanced` | A scientific, quality, or resource choice with a suitable normal default. |
| `developer` | Diagnostic, compatibility, experimental, testing, or specialist recovery control. |

Required inputs must never remain hidden. A client must reveal a conditionally
required input when its condition is active.

Visibility changes presentation only. It must not make a CLI input invalid or
unavailable.

### 5.3 User-facing text

The target presentation fields and limits are:

| Field | Use | Maximum |
| --- | --- | --- |
| Program display name | Human-readable program title. | 45 characters |
| Program summary | One-line program list text. | 90 characters |
| Program help | Explanation of purpose and expected result. | 350 characters |
| Input label | Form label. | 45 characters |
| Input help | Meaning, effect, and units where needed. | 350 characters |
| Placeholder | Example value or expected format only. | 40 characters |

Text must follow these rules:

- Use sentence case and plain English.
- Preserve necessary scientific terms, but explain an abbreviation or
  specialist term on first use.
- Put units in structured metadata and repeat them in help when clarity
  requires it.
- Do not expose internal CLI keys as labels unless the key is itself the
  clearest established name.
- Do not encode choices, defaults, units, or explanatory sentences in a
  placeholder.
- Do not put implementation notes or developer instructions in user-facing
  help.

### 5.4 Change ownership

- Program modules own program membership, scientific wording, grouping,
  visibility, and program-specific overrides.
- The common-parameter module owns reusable input definitions.
- The UI core types own representation and lifecycle.
- A single serializer must own the exported schema.
- Commander and parameter modules remain authoritative for execution and
  scientific validation.

### 5.5 Review and testing

A descriptor change must be checked for:

- CLI key and required-input compatibility;
- executable and project-file ownership;
- valid visibility;
- length limits and plain-language wording;
- structured choices and defaults;
- JSON serialization;
- dependent input visibility.

Scientific wording and Standard/Advanced classification require review by
someone familiar with the affected workflow.
