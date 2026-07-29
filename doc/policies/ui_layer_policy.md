# UI Descriptor Layer Policy

## Purpose and scope

The UI descriptor layer is the shared description of SIMPLE commands. It is
used by both the command-line interface and graphical clients. This policy
applies to the Fortran code under `src/main/ui/` and to the JSON produced from
that code.

The descriptor layer describes programs and their inputs. It does not execute
programs, independently own scientific defaults, or replace validation
performed by the parameters, commander, strategy, or domain layers. It does
expose the effective user-visible defaults supplied by the parameter layer.

The implementation work is tracked separately in
[`../refactoring_notes/ui_descriptor_layer_refactoring.md`](../refactoring_notes/ui_descriptor_layer_refactoring.md).

## One source of truth

The Fortran backend is the only source of truth for UI descriptors.

- Program modules define program names, categories, descriptions, visibility,
  executable ownership, project requirements, and input membership.
- Declaration-time initialization in `type(parameters)`, together with
  `init_dynamic_defaults`, defines the baseline parameter defaults.
- `ui_param` objects define input names, types, descriptions, choices, units,
  visibility, dependencies, and selection rules. They record a read-only snapshot of those
  baseline parameter defaults; they do not maintain an independent default.
- JSON is a serialization of the in-memory Fortran objects. JSON must not add,
  omit, reinterpret, or override descriptor data.
- TOML, Markdown, spreadsheets, and other documents must not contain an
  independently maintained copy of UI metadata.
- Python or other code generators must not generate program-specific UI
  descriptor code.

The descriptor is not extended by GUI-specific overlays, overrides, or client
metadata. A GUI consumes the same completed Fortran object as the CLI. If a
field is needed to render, validate, or submit an input, it belongs in that
shared descriptor with a client-neutral name and meaning.

Documentation may explain the model and the editing procedure, but the values
rendered by the CLI and GUI must come from Fortran.

## Default ownership and precedence

An input default is execution behavior, not display text. The parameter layer
therefore owns defaults, and the UI descriptor layer publishes them.

SIMPLE currently obtains missing input values from several places:

- declaration-time initialization of `type(parameters)`;
- initialization of dynamic `type(string)` components;
- conditional program or workflow assignments of the form
  `if (.not. cline%defined(key)) call cline%set(key, value)`; and
- values derived at runtime from other inputs, project state, or input data.

The initializations remain where they are. They are already the source used by
the parameter lifecycle and are available before any program runs. The UI
must obtain them through a read-only default snapshot: construct only the
baseline `parameters` state and apply `init_dynamic_defaults`; do not call
`parameters%new`, parse a command line, inspect a project, or run a commander.

Conditional `if (.not. cline%defined(key)) call cline%set(key, value)` calls
remain program execution behavior. They can legitimately override a baseline
value for one workflow, and the GUI cannot treat them as its universal default
without reimplementing commander logic.

This gives two distinct meanings that must not be confused:

- the UI and JSON publish the baseline value from `parameters` when one is
  available;
- the executed program may choose a different value when an optional key is
  omitted and its own command setup supplies a value.

Consequently, a GUI must omit an untouched optional input from the submitted
command line. It sends a value only when the user has deliberately changed it
(or when the input is required). This preserves the CLI's program-specific
default behavior. A GUI must not claim that a displayed baseline is always the
final value used by every program.

Runtime-derived values remain automatic; the UI omits a default when it cannot
show a sensible baseline value. Defaults must not be encoded in `placeholder`.

During migration, an existing literal default passed to `add_input` may be
retained only as an assertion against the baseline snapshot. Once a category
is validated, that literal must be removed.

## Current class design

The descriptor model uses composition:

```text
ui_hash
  +-- references module-owned ui_program objects by program name

ui_program
  +-- program identity and presentation fields
  +-- image input/output       -- linked_list of ui_param values
  +-- parameter input/output   -- linked_list of ui_param values
  +-- alternative inputs       -- linked_list of ui_param values
  +-- search controls          -- linked_list of ui_param values
  +-- filter controls          -- linked_list of ui_param values
  +-- mask controls            -- linked_list of ui_param values
  +-- computer controls        -- linked_list of ui_param values
```

`ui_hash` extends the generic `vrefhash` and supplies typed access to
`ui_program` and `ui_param`. `ui_program` and `ui_param` do not inherit from
one another.

### `ui_param`

`simple_ui_param.f90` defines one command input. Its current fields provide:

| Field | Meaning |
| --- | --- |
| `key` | Exact CLI key. This is a stable interface identifier. |
| `keytype` | Input type used by CLI and GUI consumers. |
| `label` | Short field label. |
| `help` | Full help text. |
| `placeholder` | Short example or entry hint. Legacy choice syntax may still be present during migration. |
| `cval_default`, `rval_default`, `has_default` | Current internal storage and presence flag for a display default. |
| `units` | Display units, empty when not applicable. |
| `choices` | Structured values accepted by binary and multiple-choice inputs. |
| `visibility` | Standard, Advanced, or Developer. |
| `gui_submenu`, `exclusive_group`, `active_flags`, `online` | Legacy presentation and interaction fields; scheduled for replacement by shared descriptor semantics. |

The two `set_param` overloads construct numeric and character inputs.
`refresh_legacy_choices` converts the existing placeholder choice syntax into
structured Fortran choices. `apply_gui_overrides` is legacy compatibility
machinery, not the target construction model.

Program modules must declare the final descriptor for each program input.
Common parameter definitions may provide reusable starting values, but they
must not be supplemented by renderer-only metadata or a GUI-specific overlay.
Any program-specific difference belongs in the shared descriptor before it is
registered.

The GUI needs one thing only: a sensible value to display before the user
edits an optional input. It does not need to know how SIMPLE obtained that
value. The default export therefore has no `default_kind` classification.

### `ui_program`

`simple_ui_program.f90` defines one program accepted through `prg=` or
`test=`. Its current fields provide:

| Field | Meaning |
| --- | --- |
| `name` | Exact CLI program name. |
| `category` | Stable category identifier inherited from the owning program module. |
| `category_display_name` | Plain-English category heading shared by CLI listings and JSON. |
| `category_order` | Display order within the owning executable. |
| `display_name` | Plain-English GUI title, separate from the CLI name. |
| `summary` | Short program summary. |
| `help` | Full program help. |
| `executable` | Executable that accepts the program. |
| `visibility` | Standard, Advanced, or Developer. |
| seven input lists | Inputs grouped into the existing CLI sections. |
| `sp_required` | Whether a SIMPLE project is required. |

`ui_program%new` creates the program. Its `add_input` overloads create or copy
`ui_param` values into one of the seven input lists. The same object supplies
required-key checks, CLI help, program descriptions, and JSON.

### Construction and registration

Common inputs are initialized by `set_ui_params`. The program modules then
construct programs and register them:

```text
make_ui / make_test_ui
  +-- set_ui_params
  +-- executable group constructor
      +-- simple_ui_* or single_ui_* constructor
          +-- ui_program%new
          +-- ui_program%add_input
          +-- add_ui_program
              +-- assign module category
              +-- register in ui_hash
```

Each program-defining module declares one module-local `UI_CATEGORY` descriptor
with an identifier, a plain-English heading, and an order. The identifier is
the lowercase suffix of its module name. For example, `simple_ui_denoise` owns
`denoise` with the heading `Denoising`, and `single_ui_atom` owns `atom` with
the heading `Atom Analysis`. Every call to `add_ui_program` passes that one
descriptor. Registration rejects incomplete category metadata and duplicate
program names.

This convention is deliberately open-ended. A new category is created by
adding a new program module and including its constructor in the appropriate
group module. Moving a program to a different category means moving its
Fortran construction to the new owning module.

The public `prg=list` paths traverse registered programs and group them by
this metadata. JSON serializes the same identifier, heading, and order for
each program. There are no handwritten program-list headings; changing a
category descriptor changes CLI listings and JSON together.

## Visibility

Every program and every input must have exactly one of these visibility
levels:

- **Standard**: needed for the usual workflow and suitable for most users.
- **Advanced**: useful for experienced users or less common workflows.
- **Developer**: diagnostic, experimental, implementation-specific, or unsafe
  without detailed knowledge of SIMPLE.

Fortran stores these as the constants `UI_VIS_STANDARD`,
`UI_VIS_ADVANCED`, and `UI_VIS_DEVELOPER`. JSON stores the corresponding
lowercase names.

Visibility is a shared descriptor field, not a GUI add-on. During migration,
the legacy `gui_visibility` constructor keyword assigns that field. New APIs
must use `visibility` instead. Until every descriptor is reviewed, an input or program without an explicit visibility defaults to
Developer. This preserves the conservative migration: unreviewed controls
are not exposed as ordinary Standard or Advanced controls.

Visibility changes presentation only. They must never change whether a CLI
key is accepted or how a program executes.

## User-facing text

User-facing text must explain the task in plain English. Keep scientific terms
when they are necessary for accuracy, but do not expose internal class names,
variable names, abbreviations, or implementation details as labels. CLI
program names, input keys, file formats, and established scientific terms may
remain exact.

### Program summary

For production programs, `summary` is the one-line answer to “what does this
program do?” It is shown where a user chooses a program; it is not a title and
must therefore carry useful meaning without the reader knowing the CLI name.

- Use one active phrase or sentence fragment of 30–100 characters. The
  100-character ceiling is enforced by `ui_program%new`.
- Start with the action and name the main object, result, or workflow. For
  example: `Estimate CTF parameters from micrographs`.
- State the user-visible outcome, not only an algorithm name. Include scope
  such as `streaming`, `2D`, or `nanoparticle` when it changes what the
  program is for.
- Do not repeat the program name, list controls, or describe implementation
  details that belong in `help`.
- `help` remains the complete explanation of purpose, workflow, constraints,
  and scientific behavior. The summary must agree with it and with the
  current implementation.

Developer-only test programs may retain their technical test names until their
own UI is reviewed; they are not part of the production program chooser.

### Program display name

`display_name` is the program's plain-English title in a GUI. It identifies a
program in a chooser, heading, breadcrumb, or batch-job card. `name` remains
the exact CLI identifier and must never be changed to improve presentation.

- Use a short, scannable title: normally two to seven words and preferably no
  more than 60 characters. The hard storage limit is 100 characters so an
  established technical title is not silently truncated.
- Use title case and describe the user task or result: `Estimate CTF`,
  `Create 2D Class Averages`, or `Import Particle Data`.
- Expand implementation-style CLI spellings into ordinary words where that is
  clear. Retain established scientific terms and acronyms such as CTF, FSC,
  2D, and 3D when expansion would reduce clarity.
- Do not include `prg=`, executable names, underscores, parameter keys, or a
  trailing full stop. Do not make the title a sentence-length explanation;
  that belongs in `summary` and `help`.
- The title and summary serve different reading situations: the title answers
  “which tool?”, while the summary answers “what will it do?”. They may share
  wording during migration, but category reviews must replace the fallback
  with an explicit title when a shorter or clearer one is available.

Every `ui_program` has a populated `display_name`. To preserve current UI
coverage while titles are reviewed category by category, `ui_program%new`
copies the existing plain-English `summary` when no explicit `display_name`
is supplied. Callers can supply `display_name=` now; the fallback is a
compatibility path, not another source of metadata.

### Input placeholder

`placeholder` is a compact entry hint, not a second label or help paragraph.
It must contain at most 40 characters; the descriptor constructor enforces a
single standard representation for every rendered input.

| Input kind | Rendered placeholder |
| --- | --- |
| number (`num`, `int`, `float`) | `e.g. 10` |
| file | `e.g. input.mrc` |
| directory | `e.g. /path/to/folder` |
| free text | `e.g. value` |
| choice, binary, or hidden input | empty |

Choice widgets already render their accepted values, so their placeholder
must be empty. Units, ranges, choice lists, default values, and explanatory
prose belong respectively in `units`, `help`, `choices`, `default`, and
`help`—never in a placeholder.

The current Fortran declarations still pass legacy choice syntax into the
constructor so it can initialize `ui_param%choices`. The constructor consumes
that syntax and replaces it with the standard display placeholder before JSON
or CLI help is rendered. The next choice-constructor refactoring will remove
that legacy input syntax entirely.

## JSON contract

All JSON writers must serialize the same Fortran fields with the same meaning.
At minimum, a program record includes its exact CLI `name`, `category`,
`category_display_name`, `category_order`, `display_name`, descriptions, executable, and visibility. An input record includes its key, type,
descriptions, required state, an optional display `default`, visibility, and any
applicable units, choices, and GUI behavior.

When present, `default` is a display value and may be serialized as a string
regardless of the input's eventual CLI type. The normal CLI parser remains
responsible for interpreting the submitted value. When no sensible baseline
value is available, the JSON record simply omits `default`.

Choice values in JSON must come from `ui_param%choices`, not from parsing the
placeholder again. Choice values remain exact CLI values.

Visibility is the sole presentation classification. JSON must not serialize a
second Boolean visibility flag.

## Change rules

- Preserve the one-to-one relationship between registered programs and CLI
  programs, and between registered inputs and accepted CLI keys.
- Do not rename a program or input as part of a wording cleanup.
- Define shared input metadata once in `simple_ui_params_common.f90`; override
  it in a program module only when that program genuinely differs.
- Do not add GUI-only descriptor fields, client overlays, JSON patches, or
  renderer-side defaults. Add a shared semantic field only when both the CLI
  and GUI can interpret the same meaning.
- Keep declaration-time and dynamic `parameters` initialization authoritative.
  Do not add a literal UI default that duplicates it.
- Keep conditional commander defaults, runtime-derived values, and fixed
  internal execution assignments distinct from baseline UI defaults.
- Keep category ownership in the program module. Do not add a second category
  table in a renderer, JSON writer, or GUI.
- Update all JSON paths when a descriptor field changes until serialization is
  consolidated.
- Validate duplicate keys, visibility values, categories, choices, and JSON
  before merging descriptor changes.
