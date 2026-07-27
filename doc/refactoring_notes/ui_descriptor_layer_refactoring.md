# UI Descriptor Layer Refactoring

## 1. Objective

Refactor the command descriptor layer so that:

- programs and inputs have `standard`, `advanced`, or `developer` visibility;
- user-facing names, summaries, help, placeholders, units, choices, and
  defaults are separate fields;
- JSON is produced by one versioned serializer;
- program and input descriptions can be reviewed in one Markdown catalog and
  compiled back into the UI layer;
- CLI program selection, required-key checks, help, and project-file behavior
  remain unchanged during the metadata migration.

The current implementation is documented in
[`../policies/ui_layer_policy.md`](../policies/ui_layer_policy.md).

## 2. Non-goals

This refactor will not:

- move scientific defaults or validation out of the parameters layer;
- change commander dispatch;
- change program names or CLI keys as part of wording cleanup;
- merge the command descriptor layer with runtime stream GUI metadata;
- implement the batch scheduler or job monitor.

## 3. Target object model

### 3.1 New support types

Add `src/main/ui/simple_ui_descriptor_types.f90` containing:

```fortran
integer, parameter :: UI_VIS_STANDARD  = 1
integer, parameter :: UI_VIS_ADVANCED  = 2
integer, parameter :: UI_VIS_DEVELOPER = 3

type :: ui_choice
    type(string) :: value
    type(string) :: label
    type(string) :: help
end type ui_choice

type :: ui_condition
    type(string)              :: key
    type(string), allocatable :: values(:)
end type ui_condition
```

The module must provide validation and conversion between visibility constants
and the JSON strings `standard`, `advanced`, and `developer`. Free-form
visibility strings must not be stored in Fortran objects.

`ui_choice%value` is the exact CLI value. `label` and `help` are presentation
text. This permits a GUI to display plain English while preserving CLI values
such as `prob_neigh`.

`ui_condition` replaces encoded expressions such as
`quality_mode=apply|analyze|evaluate`.

### 3.2 Target `ui_param`

Retain these execution fields:

- `key`
- `keytype`
- `cval_default`
- `rval_default`
- `required`
- `online`
- `exclusive_group`

Add:

| Field | Type | Purpose |
| --- | --- | --- |
| `label` | `type(string)` | Human-readable form label. |
| `help` | `type(string)` | Full plain-English explanation. |
| `placeholder` | `type(string)` | Short example or format hint. |
| `units` | `type(string)` | Machine-readable/displayable units. Empty when dimensionless. |
| `visibility` | integer visibility constant | Standard, Advanced, or Developer. |
| `section` | integer section constant | Existing image/parameter/alternative/search/filter/mask/compute classification. |
| `group` | `type(string)` | Optional within-program GUI group, replacing `gui_submenu`. |
| `choices(:)` | allocatable `type(ui_choice)` | Structured choices for `binary` and `multi` inputs. |
| `visible_when(:)` | allocatable `type(ui_condition)` | Conditions that activate/display the input. |
| `required_when(:)` | allocatable `type(ui_condition)` | Conditions that make an otherwise optional input required. |
| `has_default` | logical | Distinguishes no default from an empty-string or zero default. |

Remove after migration:

| Current field | Replacement |
| --- | --- |
| `descr_short` | `label` |
| `descr_long` | `help` |
| `descr_placeholder` | `placeholder` |
| `advanced` | `visibility` |
| `gui_submenu` | `group` |
| `active_flags` | `visible_when(:)` |

Stop encoding allowed values and defaults in any description or placeholder.

### 3.3 Target `ui_program`

Retain:

- `name`
- `executable`
- `sp_required`
- `exists` unless lifecycle work later makes it unnecessary

Add:

| Field | Type | Purpose |
| --- | --- | --- |
| `display_name` | `type(string)` | Plain-English program title, separate from the CLI `name`. |
| `summary` | `type(string)` | Short program-browser description. |
| `help` | `type(string)` | Full description printed by `describe=yes`. |
| `visibility` | integer visibility constant | Standard, Advanced, or Developer. |
| `inputs` | `type(linked_list)` | One ordered list of `ui_param` values. |

Remove after migration:

| Current field | Replacement |
| --- | --- |
| `descr_short` | `summary` |
| `descr_long` | `help` |
| `advanced` | `visibility` |
| `gui_submenu_list` | Derived from distinct non-empty `ui_param%group` values. |
| `img_ios`, `parm_ios`, `alt_ios`, `srch_ctrls`, `filt_ctrls`, `mask_ctrls`, `comp_ctrls` | `inputs`, with `ui_param%section` identifying the section. |

Using one list removes seven repeated traversal paths from printing, required
key collection, validation, and JSON export. Existing section order must be
preserved when rendering CLI help and JSON.

`UI_IMG`, `UI_PARM`, `UI_ALT`, `UI_SRCH`, `UI_FILT`, `UI_MASK`, and `UI_COMP`
may remain as compatibility aliases during migration. Rename them to explicit
`UI_SECTION_*` constants in the target API.

## 4. Target JSON schema

Write one root object with an integer schema version:

```json
{
  "schema_version": 2,
  "programs": [
    {
      "name": "refine3D",
      "display_name": "3D refinement",
      "summary": "Refine a 3D map from aligned particle images.",
      "help": "...",
      "executable": "simple_exec",
      "requires_project": true,
      "visibility": "standard",
      "inputs": []
    }
  ]
}
```

Each input object must contain:

```json
{
  "key": "pgrp",
  "type": "multi",
  "label": "Point-group symmetry",
  "help": "...",
  "placeholder": "e.g. C1",
  "units": "",
  "section": "search",
  "group": "model",
  "visibility": "standard",
  "required": false,
  "has_default": true,
  "default": "C1",
  "choices": [
    {"value": "c1", "label": "C1", "help": "No imposed symmetry."}
  ],
  "visible_when": [],
  "required_when": [],
  "exclusive_group": "",
  "online": false
}
```

Rules:

- `default` is emitted only when `has_default` is true.
- `choices` contains exact CLI values plus display text.
- Conditions are arrays of objects, never parsed expressions.
- `section`, `visibility`, and `type` are validated enumerations.
- Program and input order is deterministic.

## 5. Serializer consolidation

Add `src/main/ui/simple_ui_serializer.f90`. It owns conversion of a registry,
program, or input to a `json_value`.

Change:

- `print_ui_json` to serialize the registry and print it;
- `write_ui_json` to serialize the same registry and write it;
- `ui_program%write2json` to call the same program serializer.

Remove:

- the internal `create_program_entry` and `create_section_list` copies from
  `print_ui_json`;
- the duplicate copies from `write_ui_json`;
- `create_section_from_list` from `simple_ui_program`;
- all parsing of `descr_placeholder` to obtain choices.

Keep `print_stream_ui_json` separate until its hand-built workflow schema is
either retired or given its own refactoring. It is not the generic command
descriptor serializer.

## 6. Description catalog workflow

### 6.1 Decision

Use one reviewable Markdown catalog as the source of user-facing text and
visibility. Generate Fortran from that catalog; do not patch arbitrary
Fortran call sites with text substitutions.

The catalog will be:

`doc/ui/ui_descriptions.md`

Each program has a machine-readable fenced block. A restricted TOML subset is
recommended because it remains readable and supports repeated input records:

````markdown
## `simple_exec/refine3D`

```toml
name = "refine3D"
display_name = "3D refinement"
summary = "Refine a 3D map from aligned particle images."
help = "..."
visibility = "standard"

[[input]]
key = "pgrp"
label = "Point-group symmetry"
help = "Symmetry applied while building the map. Use C1 when no symmetry is known."
placeholder = "e.g. C1"
units = ""
visibility = "standard"

[[input.choice]]
value = "c1"
label = "C1"
help = "Do not impose rotational symmetry."
```
````

Each program input appears explicitly, even when its execution metadata was
copied from `simple_ui_params_common`. The stable identity is `(executable,
program name, input key)`. Display text must never be used as an identity.
Keeping exact program/input records avoids an implicit text fallback that
could apply generic wording to a program with different semantics.

### 6.2 Tools

Add two maintenance tools under `scripts/ui/`:

1. `export_ui_descriptions.py`
   - reads the current versioned UI JSON;
   - creates or refreshes the Markdown catalog;
   - preserves reviewed text when the stable identity already exists;
   - marks new and removed programs/inputs explicitly;
   - writes entries in deterministic order.

2. `generate_ui_description_catalog.py`
   - parses and validates the fenced catalog blocks;
   - checks visibility values, duplicate identities, required fields, length
     limits, placeholder rules, and choice/default consistency;
   - generates
     `src/main/ui/generated/simple_ui_description_catalog.f90`;
   - supports `--check`, which fails if generated code is stale.

The generated module contains lookup routines keyed by program and input key.
It applies:

- program display name, summary, help, and visibility;
- input label, help, placeholder, units, visibility, and choice labels.

The generated file must contain a header saying not to edit it manually.
Runtime code must not parse Markdown.

### 6.3 Injection point

During program construction:

1. `ui_program%new` establishes execution metadata.
2. The generated catalog applies program presentation metadata by program key.
3. `add_input` copies or creates the execution descriptor.
4. The catalog applies the exact program/input presentation record.
5. The completed `ui_param` is appended to `ui_program%inputs`.

This keeps CLI and GUI consumers on the same assembled `ui_program` object
while removing prose from hundreds of constructor calls. A missing
program/input catalog record is a construction error once migration is
complete; it must not fall back silently to the CLI key.

### 6.4 Catalog editing process

Descriptions are updated by domain group:

1. Generate the catalog from the current UI registry.
2. Select one group, beginning with project and preprocessing programs.
3. Classify each program and input as Standard, Advanced, or Developer.
4. Rewrite display name/label, summary/help, placeholder, and units.
5. Have a workflow expert check scientific meaning and defaults.
6. Run the catalog validator and regenerate the Fortran catalog.
7. Generate JSON and review the rendered Standard-only form.
8. Commit one domain group at a time.

Do not combine wording changes with CLI key, default, requiredness, or
commander behavior changes.

## 7. Implementation sequence

### Step 1: Add characterization tests

Before changing types:

- snapshot registered program names, executables, `sp_required`, input keys,
  sections, requiredness, defaults, and input order;
- test `get_nrequired_keys` and `get_required_keys`, including alternative
  inputs;
- test the current JSON output and all supported key types;
- record every placeholder from which an option list is currently parsed.

Acceptance criterion: the baseline tests pass without changing production
behavior.

### Step 2: Add visibility and structured choices

- Add `simple_ui_descriptor_types`.
- Add `visibility`, `choices`, `units`, and `has_default` to the existing
  types without removing old fields.
- Accept temporary `gui_visibility` arguments in `new`, `add_input`, and
  `apply_gui_overrides`.
- Map old `gui_advanced=.false.` to Standard and `.true.` to Advanced.
- Classify unspecified legacy descriptors as Developer so that an omission
  cannot accidentally expose a specialist control.
- Populate structured choices alongside the legacy placeholder encoding.

Acceptance criterion: old CLI behavior and old JSON remain unchanged; the new
fields can be validated independently.

### Step 3: Add the catalog round trip

- Implement the Markdown export and generator tools.
- Generate the initial catalog without rewriting text.
- Add the generated catalog lookup during `new` and `add_input`.
- Add a CI/test target that runs the generator with `--check`.

Acceptance criterion: exporting, generating, rebuilding, and exporting again
does not change stable identities or execution metadata.

### Step 4: Consolidate serialization

- Add `simple_ui_serializer`.
- Switch all generic JSON entry points to it.
- Emit schema version 2 and the new structured fields.
- Optionally provide a temporary schema version 1 compatibility writer if a
  current external client still consumes the old shape.

Acceptance criterion: all entry points serialize identical program/input
content for the same registry.

### Step 5: Replace the seven lists

- Add `section` to `ui_param`.
- Add the single `inputs` list to `ui_program`.
- Make `add_input` translate the current `UI_*` selector to `section`.
- Update required-key collection, CLI help, validation, and serialization to
  traverse `inputs`.
- Compare output with the characterization baseline.
- Remove the seven old lists after all consumers use `inputs`.

Acceptance criterion: program input order, CLI help grouping, and required-key
behavior match the baseline.

### Step 6: Migrate conditions

- Convert every `active_flags` expression to `visible_when`.
- Add `required_when` where program semantics require a value conditionally.
- Teach the GUI client to reveal any active required input regardless of
  visibility.
- Remove string condition parsing after all declarations are migrated.

Acceptance criterion: condition behavior has tests for every controlling key.

### Step 7: Rewrite descriptions by domain

Use the catalog process in section 6.4. Suggested order:

1. project and preprocessing;
2. 2D classification and class-average processing;
3. ab initio 3D and refinement;
4. image, volume, mask, filtering, and utility programs;
5. stream programs;
6. SINGLE programs;
7. test programs.

Acceptance criterion: every shipped descriptor passes policy validation and
has domain review.

### Step 8: Remove compatibility fields and arguments

Remove only after schema-2 clients and all declarations have migrated:

- `descr_short`, `descr_long`, and `descr_placeholder`;
- description arguments and description-override arguments from `set_param`,
  `ui_program%new`, and `ui_program%add_input`;
- `advanced` and `gui_advanced`;
- `gui_submenu` and `gui_submenu_list`;
- `active_flags`;
- placeholder option/default encodings;
- the seven program input lists;
- the old serializer implementations;
- schema-1 output, after an agreed compatibility window.

Acceptance criterion: repository search finds no legacy field, argument, or
encoded-placeholder use, and all CLI/UI tests pass.

## 8. Test plan

Add a focused descriptor-policy test suite covering:

- valid and invalid visibility constants;
- unique program keys and unique input keys per program;
- supported key types and sections;
- explicit `has_default` behavior for zero and empty-string defaults;
- binary and multi choice counts;
- choice values matching accepted CLI values;
- label, summary, help, and placeholder limits;
- placeholders containing no `{default}` or `(a|b)` encodings;
- condition references to keys present in the same program;
- required inputs visible under active conditions;
- deterministic JSON and catalog generation;
- CLI/GUI round-trip for representative Standard, Advanced, and Developer
  programs.

Keep the existing GUI metadata and assembler tests. They test the separate
runtime-status GUI subsystem and do not replace descriptor tests.

## 9. Review checkpoints

Pause for review after:

1. agreeing the target field list and schema-2 JSON example;
2. generating the first Markdown catalog;
3. completing one representative program end to end;
4. migrating each domain group;
5. removing schema-1 compatibility.

The first representative program should contain numeric, file, binary, and
multi inputs, a common-parameter override, a conditional input, and inputs at
all three visibility levels. `refine3D` is a suitable candidate once its
scientific wording owner is available.
