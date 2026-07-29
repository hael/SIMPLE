# UI Descriptor Layer Refactoring

## 1. Objective

Refactor the command descriptor layer so that:

- programs and inputs have `standard`, `advanced`, or `developer` visibility;
- user-facing names, summaries, help, placeholders, units, choices, and
  defaults are separate fields;
- JSON is produced by one versioned serializer;
- the suite and group structure currently shown by `prg=list` becomes explicit
  metadata shared by CLI and GUI consumers;
- one group-owned TOML catalog is the single editable source for its
  group, programs, inputs, navigation, and presentation metadata;
- generated Fortran constructs the descriptor registry from that catalog;
- CLI program selection, required-key checks, help, and project-file behavior
  remain unchanged during the metadata migration.

The current implementation is documented in
[`../policies/ui_layer_policy.md`](../policies/ui_layer_policy.md).

## 2. Non-goals

This refactor will not:

- move scientific defaults or validation out of the parameters layer;
- generate commander, execution-router, or scientific implementation code;
- change commander dispatch;
- change program names or CLI keys as part of wording cleanup;
- merge the command descriptor layer with runtime stream GUI metadata;
- implement the batch scheduler or job monitor.

## 3. Target object model

### 3.1 New support types

Keep `src/main/ui/simple_ui_visibility.f90` as the small, independent owner of
the visibility enumeration. Add `src/main/ui/simple_ui_descriptor_types.f90`
for the remaining shared presentation types, beginning with:

```fortran
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

The visibility module must provide validation and conversion between visibility
constants and the JSON strings `standard`, `advanced`, and `developer`.
Free-form visibility strings must not be stored in Fortran objects.

`ui_choice%value` is the exact CLI value. `label` and `help` are presentation
text. This permits a GUI to display plain English while preserving CLI values
such as `prob_neigh`.

`ui_condition` replaces encoded expressions such as
`quality_mode=apply|analyze|evaluate`.

### 3.2 Suites and program groups

The existing `prg=list` headings are part of the UI information architecture,
not incidental formatting. Today this structure is implicit in module
boundaries and repeated `print_*_programs` routines:

- `simple_exec` has grouped project, preprocessing, 2D, 3D, filtering,
  image, volume, validation, utility, and related program families;
- `single_exec` has time-series, trajectory, Nano 2D, Nano 3D, atom, map, and
  validation families;
- `simple_test_exec` has class, FFT, geometry, high-level, I/O, mask, network,
  numerics, optimization, parallel, SINGLE, statistics, and utility families;
- `simple_stream` has one flat Stream Workflows family.

Represent this explicitly with typed suite and group descriptors. A target
shape is:

```fortran
integer, parameter :: UI_LAYOUT_GROUPED = 1
integer, parameter :: UI_LAYOUT_FLAT    = 2

type :: ui_program_group
    type(string)              :: id
    type(string)              :: title
    integer                   :: display_order
    type(string), allocatable :: program_names(:)
end type ui_program_group

type :: ui_suite
    type(string)                         :: id
    type(string)                         :: executable
    type(string)                         :: display_name
    integer                              :: layout
    type(ui_program_group), allocatable  :: groups(:)
end type ui_suite
```

Initially the stable suite ID is the executable name: `simple_exec`,
`single_exec`, `simple_test_exec`, or `simple_stream`. Group IDs are stable,
lowercase machine identifiers such as `project`, `preprocessing`, and
`refinement3d`. Group titles such as "Project management" are presentation
text supplied by the description catalog.

The suite `executable` identifies the entry point that presents the suite. A
program record retains its own `executable` field. Most programs match their
suite entry point, but shared programs legitimately use `all` and may appear
in more than one executable's command surface under the current CLI routing
rule. Suite placement and program execution scope must therefore remain
separate fields.

The catalog directory is the single editable source of truth for the UI
descriptor layer:

- the suite catalog owns suite ID, executable, display name, and layout;
- each group TOML file owns its group ID, title, order, programs, inputs, and
  all descriptor fields; program declaration order is its display order;
- generated Fortran owns no independent metadata and must never be edited;
- `prg=list`, CLI help, required-key checks, JSON, and GUI clients consume the
  generated registry rather than separate declarations.

All `ui_suite`, `ui_program_group`, `ui_program`, and `ui_param` construction
code is derived from the catalog. Handwritten code outside the descriptor
layer continues to own executable entry points, dispatch, typed runtime
parameters, scientific validation, and workflow behavior. The catalog does
not generate commanders or algorithms. Program names, input keys, types,
defaults, and requiredness are boundary contracts: they are declared once for
the descriptor layer in the catalog and checked against the handwritten
execution/parameters layer.

The generic Fortran type definitions, generated-registry facade interface,
serializer, renderers, and validators remain handwritten infrastructure
because they implement behavior rather than program-specific data. They must
contain no program names, group names, input keys, descriptions, or other
catalog facts.

The hierarchy is open data, not a closed Fortran enumeration. Adding a group
means adding a group TOML file with a new stable ID and order. Moving a
program means moving its complete catalog block between group files. Changing
a title is a presentation-only edit. None of these operations should require
serializer, registry-construction, CLI-rendering, or GUI implementation
changes.

The catalog must support these navigation-only operations:

- add a group by adding one TOML file with a unique ID, title, and sparse
  integer order;
- rename or reorder a group without changing its stable ID;
- move or reorder a program by moving or editing its complete program block;
- remove an empty group after all of its programs have been reassigned.

Use order values with gaps, such as 10, 20, and 30, so most insertions do not
renumber unrelated records. Adding an executable program requires a catalog
record plus its handwritten execution routing and implementation, but no
handwritten UI descriptor constructor. Schema version 2 supports one group
level below each suite. If a real use case later requires nested groups, add an
explicit, versioned parent/child model rather than encoding hierarchy in group
IDs or titles.

The generator compiles the suite and group catalogs into the complete
descriptor registry. Static validation must reject duplicate program
definitions, duplicate group IDs or order values, invalid descriptor fields,
and incompatible suite/executable combinations. Post-build contract tests
must reject catalog programs without execution routing and descriptor keys
that are not recognized by the parameters layer. Thus the catalog can define
and reorganize descriptors but cannot silently invent executable behavior.

Group IDs should remain stable because clients may retain expanded/collapsed
state or saved batch views by ID. A deliberate ID replacement is a migration,
whereas changing the displayed title is not. `simple_stream` still has one
group in the model, but clients may suppress its heading because its suite
layout is `flat`.

The compiled grouped registry becomes the source for both `prg=list` and JSON.
Once it is established, remove the handwritten `print_*_programs` traversals
rather than maintaining a second navigation hierarchy.

### 3.3 Target `ui_param`

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

### 3.4 Target `ui_program`

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
| `suite_id` | `type(string)` | Stable owning suite ID, normally equal to `executable`. |
| `group_id` | `type(string)` | Stable owning program-group ID. |
| `display_order` | integer | Order within the owning program group. |
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
  "suites": [
    {
      "id": "simple_exec",
      "executable": "simple_exec",
      "display_name": "SIMPLE",
      "layout": "grouped",
      "groups": [
        {
          "id": "refinement3d",
          "title": "3D refinement",
          "order": 60,
          "programs": [
            {
              "name": "refine3D",
              "executable": "simple_exec",
              "display_name": "3D refinement",
              "summary": "Refine a 3D map from aligned particle images.",
              "help": "...",
              "order": 50,
              "requires_project": true,
              "visibility": "standard",
              "inputs": []
            }
          ]
        }
      ]
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
- Every program belongs to exactly one group in its executable suite.
- Suite, group, program, and input order is deterministic.
- A `flat` suite still serializes its single group so clients receive one
  consistent schema.

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

Do not use one repository-wide file or one file per program. Use one TOML file
per program group, plus one small suite TOML file. TOML is source code for the
generator, not prose documentation: it gives strict parsing and validation
without Markdown fences or a second syntax to maintain. This follows the
ownership boundaries already visible in `prg=list` while keeping the number and
size of editable files manageable:

```text
src/main/ui/catalog/
├── simple_exec/
│   ├── suite.toml
│   ├── project.toml
│   ├── preprocessing.toml
│   ├── cluster2d.toml
│   └── refinement3d.toml
├── single_exec/
│   ├── suite.toml
│   ├── time_series.toml
│   └── trajectory.toml
├── simple_test_exec/
│   ├── suite.toml
│   ├── io.toml
│   └── fft.toml
└── simple_stream/
    ├── suite.toml
    └── workflows.toml
```

`suite.toml` contains only suite metadata:

```toml
id = "simple_exec"
executable = "simple_exec"
display_name = "SIMPLE"
layout = "grouped"
order = 10
```

Each group TOML file contains its navigation metadata and complete descriptor
records for every program and input in that group. The routine form is compact:

```toml
suite_id = "simple_exec"
group_id = "refinement3d"
group_title = "3D refinement"
group_order = 60

[[program]]
name = "refine3D"
executable = "simple_exec"
display_name = "3D refinement"
summary = "Refine a 3D map from aligned particle images."
help = "..."

[[program.input]]
key = "pgrp"
label = "Point-group symmetry"
help = "Symmetry applied while building the map. Use C1 when no symmetry is known."
placeholder = "e.g. C1"
choices = ["cn", "dn", "t", "o", "i"]
```

The generator applies these defaults when a field is omitted:

| Scope | Default |
| --- | --- |
| group | `review_status = "legacy"` |
| program | `help = ""`, `visibility = "standard"`; declaration order is display order |
| input | `help = ""`, `placeholder = ""`, `units = ""`, `visibility = "standard"`, `choices = []` |
| string choice | identical CLI value and display label, with no choice-specific help |

Do not write empty strings or `visibility = "standard"` merely to restate a
default. Labels, display names, summaries, stable IDs, executable scope, and
group order remain explicit because they do not have a safe general default.

Use a string list for ordinary choices. Use an expanded inline table only when
the visible label or choice-specific help differs from the CLI value:

```toml
choices = [
  { value = "prob_snhc", label = "Probabilistic search" },
  { value = "prob", label = "Probabilistic alignment", help = "Use the previous alignment model." },
]
```

The validator rejects the old `[[program.input.choice]]` table form. This keeps
the common case short and makes a detailed choice an intentional exception.

`review_status` is `legacy` or `reviewed`. Group files marked `legacy` may
preserve current wording and encoded placeholders during structural migration.
Group files marked `reviewed` must pass all wording, length, placeholder, and
choice-help rules. CI must always apply the strict rules to every file marked
`reviewed`; strict validation must not depend on a developer remembering an
optional flag.

Each program input appears explicitly, even when its execution metadata was
copied from `simple_ui_params_common`. Stable identities are:

- suite: `suite_id`;
- group: `(suite_id, group_id)`;
- program: `(suite_id, program name)`;
- input: `(suite_id, program name, input key)`.

Display text and file paths must never be used as identities. Moving a program
between groups means moving its complete block, but does not change its stable
program identity.

The catalog is authoritative after migration:

- developers edit catalog TOML, never generated Fortran;
- generated Fortran may be deleted and reproduced without information loss;
- JSON, `prg=list`, CLI help, and GUI metadata are downstream products;
- the legacy importer is a one-time migration tool, not a reverse
  synchronization path;
- missing catalog data is an error, never a fallback to handwritten descriptor
  text.

### 6.2 Tools

Add two maintenance tools under `scripts/ui/`:

1. `import_legacy_ui_catalog.py`
   - reads the characterized legacy UI registry and `prg=list` grouping;
   - creates the initial suite and group TOML files;
   - places programs with unresolved membership in an explicit `unassigned`
     group and report rather than guessing;
   - writes only to an empty target catalog unless an explicit migration
     override is given;
   - is retired from the normal workflow after catalog cutover.

2. `generate_ui_description_catalog.py`
   - discovers suite and group TOML files;
   - parses each file as TOML;
   - performs the static validation in section 6.3;
   - generates complete suite, group, program, and input registry-construction
     code under `src/main/ui/generated/`;
   - supports `--check`, which fails if generated code is stale.

Generate one Fortran module per catalog group so source and generated compiler
units have the same manageable ownership boundary. A small generated facade
constructs the complete registry. Generated code supplies:

- suite layout and display name;
- group title, order, and program membership;
- every `ui_program` identity, execution-contract field, ordering field, and
  presentation field;
- every `ui_param` identity, type, section, default, requiredness, condition,
  choice, ordering, and presentation field.

Generated files must contain a header saying not to edit them manually.
Runtime code must not parse TOML. No tool may regenerate or update
the authoritative catalog from generated Fortran or emitted JSON after
cutover.

### 6.3 Validation and code-generation safety

Validation has two distinct stages.

Static validation runs before generation and must reject:

- invalid TOML or unrecognized fields;
- missing or wrongly typed required fields;
- invalid suite layout, visibility, input type, or section values;
- invalid stable IDs and CLI keys;
- duplicate suite, group, program, input, or choice identities;
- duplicate or negative suite or group order values;
- empty group files;
- a catalog program assigned to multiple groups;
- control characters, unsupported multiline values, and text exceeding the
  applicable character limit;
- legacy placeholder choice/default encodings in `reviewed` files;
- empty labels or required help text in `reviewed` files;
- a default whose value does not match the declared input type or choices;
- invalid required/default combinations;
- conditions referring to unknown controlling keys.

The post-build contract test compares the generated descriptor registry with
the handwritten execution and parameters layers and must reject:

- catalog programs that have no route in the named executable;
- programs intended for CLI/GUI exposure that are missing from the catalog;
- catalog input keys or types unknown to the typed parameters layer;
- defaults or accepted choice values that disagree with runtime parameter
  contracts;
- suite/executable mismatches;
- differences between `prg=list`, JSON grouping, and the compiled grouped
  registry.

After cutover, a repository ownership check must also reject production
handwritten `ui_program` or `ui_param` descriptor construction outside the
generated directory. This prevents a second editable registry from
reappearing.

The generator must never copy arbitrary source fragments into Fortran. It
accepts scalar TOML data only, validates identifier fields against restricted
patterns, escapes every string through one tested Fortran-literal encoder, and
writes deterministic output. The validation pipeline is:

1. parse;
2. schema and policy validation;
3. generate to a temporary location;
4. compare or atomically replace generated files;
5. compile the generated modules;
6. run the registry contract test.

The CMake target `simple_ui_catalog_check` performs static validation,
generated-source `--check`, and the source-ownership check before the SIMPLE
library is built. A separate test target performs the post-build execution and
parameters contract comparison. A successful generator run alone is not
sufficient evidence that the catalog is correct.

### 6.4 Assembly and rendering

The generated facade constructs each `ui_suite`, `ui_program_group`,
`ui_program`, and `ui_param` directly from catalog data and inserts the
completed objects into the registry. There is no handwritten base descriptor
followed by a generated presentation overlay.

This keeps CLI and GUI consumers on the same assembled `ui_program` object and
the same grouped registry while removing the handwritten constructor calls.
Missing catalog data is always a construction or validation error, including
during the `legacy` wording phase; `review_status` relaxes wording policy only,
not structural completeness.

`prg=list` renders the grouped registry. JSON serializes the same suites,
groups, and ordering. The GUI may render a `grouped` suite as sections and a
`flat` suite without its sole group heading.

### 6.5 Catalog editing process

Descriptions are updated by domain group:

1. Select one authoritative group file, beginning with project and
   preprocessing programs.
2. Edit only that catalog file; do not edit generated Fortran.
3. Adjust group membership or ordering deliberately when the current
   `prg=list` organization needs improvement.
4. Classify each program and input as Standard, Advanced, or Developer.
5. Rewrite display name/label, summary/help, placeholder, units, and choice
   help.
6. Have a workflow expert check scientific meaning, choices, and defaults.
7. Change the completed group file to `review_status = "reviewed"`.
8. Run static validation, regenerate, compile, and run the registry contract
   test.
9. Generate JSON and review `prg=list` plus the rendered Standard-only form.
10. Commit one group at a time.

Do not combine wording changes with CLI key, default, requiredness, or
commander behavior changes.

### 6.6 Prototype disposition

The single-file catalog and generated lookup module created during the initial
experiment were prototypes. The legacy catalog has been split into group-owned
TOML source files under `src/main/ui/catalog/` and the prototype has been
removed; do not reintroduce it. The current lookup generator remains an
interim compatibility layer until complete generated registry construction
replaces the presentation overlay. The importer is a one-time migration tool,
not a normal catalog-editing path.

### 6.7 Current stopping point (2026-07-29)

The completed groundwork is deliberately limited to presentation and
navigation support:

- `ui_visibility`, `ui_choice`, and generic suite/group descriptors exist,
  with focused visibility, catalog, and navigation tests;
- `src/main/ui/catalog/` contains 155 legacy presentation records in 27 group
  TOML files for `simple_exec`, `single_exec`, and `simple_stream`;
- the catalog generator validates plain TOML, compact and expanded choices,
  explicit fields, and generated-output freshness;
- `simple_ui_catalog_check` is a build dependency; and
- generated presentation data overlays the existing handwritten UI
  constructors without changing the CLI contract.

This is **not yet** the complete descriptor source of truth. Handwritten
constructors still own program construction, input type and section,
requiredness, defaults, conditions, `sp_required`, and execution ownership.
The current `prg=list` traversal is also still handwritten. There is no
`simple_test_exec` catalog yet because the seed presentation source did not
contain test-program records.

Start the next session by completing the Step 1 characterization baseline for
the full registry: program membership, executable scope, input order and
section, requiredness, defaults, conditions, `sp_required`, and JSON output.
Then extend the TOML schema with those semantic fields and generate complete
registry construction for one suite at a time. Keep the existing handwritten
constructors until the generated suite matches that baseline; do not begin the
plain-English wording review or remove constructors before this equivalence
test exists.

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

### Step 3: Formalize extensible suite and group navigation

- Add generic `ui_suite` and `ui_program_group` descriptors; do not use a
  closed enum of group names.
- Characterize the current `simple_exec`, `single_exec`, `simple_test_exec`,
  and `simple_stream` headings, membership, and order.
- Define the suite/group catalog schema and exercise generic rendering with
  characterization fixtures.
- Represent `simple_stream` as a `flat` suite with one explicit group.
- Do not create a second permanent handwritten group registry.

Acceptance criterion: the generic CLI/JSON/GUI rendering APIs accept arbitrary
validated group fixtures without group-specific code.

### Step 4: Add the grouped catalog round trip

- Adapt the exporter into a one-time importer that creates `suite.toml` and
  one TOML catalog per group.
- Add the complete suite, group, program, and input schema, including ordering
  and `review_status`.
- Implement static validation and safe deterministic Fortran generation.
- Generate the initial `legacy` catalog without rewriting text.
- Generate complete registry construction and switch `make_ui`/`make_test_ui`
  to the generated facade.
- Compare the generated registry with the characterization baseline.
- Remove the handwritten UI program constructors, group print routines, and
  `simple_ui_params_common` after no production consumers remain.
- Add `simple_ui_catalog_check` and a post-build registry contract test.
- Remove the monolithic prototype catalog and generated module after the
  grouped outputs replace them.

Acceptance criterion: the group TOML catalog is the only editable
descriptor definition; generated sources are reproducible, the generated
registry matches the legacy characterization baseline, and malformed fixtures
are rejected at each validation layer.

### Step 5: Consolidate serialization

- Add `simple_ui_serializer`.
- Switch all generic JSON entry points to it.
- Emit schema version 2 and the new structured fields.
- Optionally provide a temporary schema version 1 compatibility writer if a
  current external client still consumes the old shape.

Acceptance criterion: all entry points serialize identical program/input
content for the same registry.

### Step 6: Replace the seven lists

- Add `section` to `ui_param`.
- Add the single `inputs` list to `ui_program`.
- Make `add_input` translate the current `UI_*` selector to `section`.
- Update required-key collection, CLI help, validation, and serialization to
  traverse `inputs`.
- Compare output with the characterization baseline.
- Remove the seven old lists after all consumers use `inputs`.

Acceptance criterion: program input order, CLI help grouping, and required-key
behavior match the baseline.

### Step 7: Migrate conditions

- Convert every `active_flags` expression to `visible_when`.
- Add `required_when` where program semantics require a value conditionally.
- Teach the GUI client to reveal any active required input regardless of
  visibility.
- Remove string condition parsing after all declarations are migrated.

Acceptance criterion: condition behavior has tests for every controlling key.

### Step 8: Rewrite descriptions by domain

Use the catalog process in section 6.5. Suggested order:

1. project and preprocessing;
2. 2D classification and class-average processing;
3. ab initio 3D and refinement;
4. image, volume, mask, filtering, and utility programs;
5. stream programs;
6. SINGLE programs;
7. test programs.

Acceptance criterion: every shipped descriptor passes policy validation and
has domain review.

### Step 9: Remove compatibility fields and arguments

Remove only after schema-2 clients and all declarations have migrated:

- `descr_short`, `descr_long`, and `descr_placeholder`;
- description arguments and description-override arguments from `set_param`,
  `ui_program%new`, and `ui_program%add_input`;
- `advanced` and `gui_advanced`;
- `gui_submenu` and `gui_submenu_list`;
- `active_flags`;
- placeholder option/default encodings;
- the seven program input lists;
- handwritten UI program/group constructor and print modules;
- `simple_ui_params_common`, after its descriptor declarations are generated;
- the legacy import tool from normal maintenance workflows;
- the old serializer implementations;
- schema-1 output, after an agreed compatibility window.

Acceptance criterion: repository search finds no legacy field, argument, or
encoded-placeholder use, and all CLI/UI tests pass.

## 8. Test plan

Add a focused descriptor-policy test suite covering:

- valid and invalid visibility constants;
- unique suite and group IDs plus valid `grouped`/`flat` layouts;
- exactly one group definition for every catalog program;
- execution routes intended for CLI/GUI exposure covered by the catalog;
- deterministic suite, group, and program ordering;
- extensibility fixtures that add, reorder, rename, and remove a group and
  move a program without changing renderer code;
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
- identical grouping and ordering in `prg=list`, schema-2 JSON, and the
  constructed UI registry;
- generator fixtures for invalid TOML, unknown fields, duplicate identities,
  invalid ordering, unknown program/input references, control characters, and
  difficult Fortran string escaping;
- compilation of all generated Fortran modules;
- deletion and exact regeneration of every generated descriptor module;
- stale-check failure after a deliberate generated-source edit;
- source-ownership failure for a handwritten production descriptor fixture;
- CLI/GUI round-trip for representative Standard, Advanced, and Developer
  programs.

Keep the existing GUI metadata and assembler tests. They test the separate
runtime-status GUI subsystem and do not replace descriptor tests.

## 9. Review checkpoints

Pause for review after:

1. agreeing the target field list and schema-2 JSON example;
2. agreeing the suite/group ownership and extensibility model;
3. generating the first split suite/group catalog;
4. demonstrating add-group, move-program, and rename-title operations;
5. cutting over one suite to catalog-only generated descriptor construction;
6. completing one representative program end to end;
7. migrating each domain group;
8. removing schema-1 compatibility.

The first representative program should contain numeric, file, binary, and
multi inputs, a common-parameter override, a conditional input, and inputs at
all three visibility levels. `refine3D` is a suitable candidate once its
scientific wording owner is available.
