# UI Descriptor Layer Policy

## Purpose and scope

The UI descriptor layer is the shared description of SIMPLE commands. It is
used by both the command-line interface and graphical clients. This policy
applies to the Fortran code under `src/main/ui/` and to the JSON produced from
that code.

The descriptor layer describes programs and their inputs. It does not execute
programs, independently own scientific defaults, or replace validation
performed by the parameters, commander, strategy, or domain layers. CMake
generates a read-only Fortran lookup module from the existing parameter and
commander declarations; the UI uses that module only to expose display
defaults.

## One source of truth

The Fortran backend is the only source of truth for UI descriptors.

- Program modules define program names, categories, descriptions, visibility,
  executable ownership, project requirements, and input membership.
- Declaration-time initialization in `type(parameters)`, together with
  `init_dynamic_defaults`, defines the baseline parameter defaults.
- `ui_param` objects define reusable input names, types, descriptions, choices,
  units, and display defaults. Program-specific semantics are stored in
  `ui_program_input` bindings.
- JSON is a serialization of the in-memory Fortran objects. JSON must not add,
  omit, reinterpret, or override descriptor data.
- TOML, Markdown, spreadsheets, and other documents must not contain an
  independently maintained copy of UI metadata.
- Python or other code generators must not generate program metadata or UI
  constructors. The CMake-generated Fortran default lookup is permitted
  because it is derived only from the existing parameter and commander source
  and contains no UI structure, descriptions, categories, or choices.

JSON and graphical clients do not apply a separate runtime overlay: they
consume the Fortran descriptor registered for the CLI.

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
the parameter lifecycle and are available before any program runs. CMake runs
`scripts/default_audit.py` to generate `simple_ui_default_values.f90` in the
build tree. It extracts baseline values and only exact, statically verified
missing-key commander overrides. `ui_program` consults the generated
`get_ui_default` routine while constructing each input. Commanders and the
parameter lifecycle never read this UI-only module.

The generator preserves the parameter parser's scalar meaning before it emits
a display value. Integer and real values are represented in the UI's
real-valued numeric storage (for example, an integer `3` becomes `3.0`). The
legacy numeric token `no` represents the command line's initialized numeric
value and is exported as `0.0`. Any other nonnumeric value for a known numeric
parameter is a generation error; an invalid value must never reach a UI
descriptor.

A global parameter baseline can be valid for the CLI but outside the narrower
choice set declared by one program. Before applying a generated binary or
multiple-choice value, `ui_param` checks the program's declared choices. An
incompatible baseline leaves that program's already-validated local default in
place; without such a default, descriptor construction fails. This is shared
Fortran descriptor validation, not a GUI override.

Conditional `if (.not. cline%defined(key)) call cline%set(key, value)` calls
remain program execution behavior. They can legitimately override a baseline
value for one workflow, and the GUI cannot treat them as its universal default
without reimplementing commander logic.

The current implementation has two distinct meanings that must not be
confused:

- the UI and JSON publish the generated baseline or verified program value
  when one is available; and
- the executed program may choose a different value when an optional key is
  omitted and its own command setup supplies a value.

## Current class design

The descriptor model uses composition:

```text
ui_hash
  +-- references module-owned ui_program objects by program name

ui_program
  +-- program identity and presentation fields
  +-- program-local groups     -- derived from input bindings
  +-- image input/output       -- linked_list of ui_program_input values
  +-- file input/output        -- linked_list of ui_program_input values
  +-- parameter input/output   -- linked_list of ui_program_input values
  +-- search controls          -- linked_list of ui_program_input values
  +-- filter controls          -- linked_list of ui_program_input values
  +-- mask controls            -- linked_list of ui_program_input values
  +-- computer controls        -- linked_list of ui_program_input values

ui_program_input
  +-- reusable ui_param definition
  +-- section, group, visibility, and activation semantics

ui_requirement_group
  +-- named key set and minimum/maximum selected counts
```

`UI_ALT` and `alt_ios` have been removed. Inputs remain in their meaningful
CLI sections—image input/output, file input/output, parameter input/output,
search, filter, mask, or computer controls. A relationship between inputs is
represented separately by a program-owned requirement group; it never creates
another input section.

`UI_IMG` is reserved for image stacks and volumes, including output stacks and
output volumes. `UI_FILE` contains every other path-like input or output:
project, STAR, orientation, CTF, coordinate, and table files, together with
directories. A list *of* stacks or volumes is still a table file and therefore
uses `UI_FILE`. `UI_PARM` is for non-path scalar, choice, and Boolean settings.

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
| `placeholder` | Short example or entry hint; empty for choice inputs. |
| `cval_default`, `rval_default`, `has_default` | Current internal storage and presence flag for a display default. |
| `units` | Display units, empty when not applicable. |
| `choices` | Structured values accepted by binary and multiple-choice inputs. |

The two `set_param` overloads construct numeric and character inputs. Binary
and multiple-choice inputs must pass `choices=ui_choices([...])`; this creates
an explicit `ui_choice` array with exact CLI values and matching display labels.
The descriptor rejects missing, empty, duplicate, or incorrectly sized choice
lists, as well as an optional default that is not one of the declared values.
`ui_param` deliberately contains no visibility, group, activation, or renderer
field.

### `ui_program_input` and groups

`ui_program_input` is the registered use of one `ui_param` in one program. It
owns the input's CLI-help section, Standard/Advanced/Developer visibility,
optional input group, and optional activation predicate. The three
`ui_program%add_input` overloads construct this binding directly; their
client-neutral optional arguments are `group=`, `visibility=`, and
`activation=`.

`group` creates or reuses a program-local `ui_input_group` with a stable id,
plain-English label, and first-use order. `ui_program%groups` is derived while
bindings are added, so no separate menu-list field can drift from actual input
membership. JSON emits the ordered program group list and an input's group
object.

The current structured activation form is
`ui_activation_equals_any(key, values)`. It records that a binding applies
when the controlling CLI key equals one of explicit values. All former
`quality_mode=...` pipe-delimited strings now use this form and JSON emits an
`activation` object with `key` and `equals_any`. The predicate is descriptor
data, not a renderer expression.

The GUI needs one thing only: a sensible value to display before the user
edits an optional input. It does not need to know how SIMPLE obtained that
value. The default export therefore has no `default_kind` classification.

### Requirement groups

`ui_requirement_group` records an input condition in shared descriptor data.
It has a stable id, a plain-English label and explanation, a set of registered
CLI keys, and inclusive `min_selected` and `max_selected` cardinalities. For
example, an image operation may require exactly one of `stk` and `vol1`.

Requirement members remain in their normal sections and retain their own
labels, help, visibility, activation, and defaults. A group is not a GUI-only
radio widget: CLI, JSON, and GUI clients all receive the same cardinality
rule. The registry rejects an empty group, duplicate member, duplicate group
id, invalid cardinality, or member key that is not an input of that program.

The command-line parser evaluates requirement groups after it has parsed all
provided keys. It prints the program's command guidance only when a group is
unsatisfied, including the plain-English rule, accepted keys, supplied count,
and required cardinality. When the group is satisfied, normal parsing proceeds
without printing usage. Requirement groups describe unconditional key presence
only; dependent or value-specific rules remain in activation predicates and
commander validation until a richer shared rule is defined.

Requirement guidance is intentionally compact. `Accepted keys`, `Supplied`,
and `Required` are each one trimmed line; fixed-length internal buffers must
be trimmed before they are concatenated for CLI output. Do not print empty
sections or padded lines.

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
| `groups` | Program-local groups derived from input bindings. |
| seven input lists | Inputs grouped into the existing CLI sections. |
| `requirements` | Program-owned input cardinality rules. |
| `sp_required` | Whether a SIMPLE project is required. |

`ui_program%new` creates the program. Its `add_input` overloads create or copy
a `ui_param` into a `ui_program_input` binding in one of the seven input
lists. The same object supplies required-key checks, CLI help, program
descriptions, group metadata, and JSON.

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
              +-- create ui_program_input binding
              +-- derive program group metadata
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

Visibility is a shared descriptor field, not a GUI add-on. Program constructors
and program-input bindings use the client-neutral `visibility` argument. Until
every descriptor is reviewed, an input or program without an explicit value
defaults to Developer. This preserves the conservative migration: unreviewed
controls are not exposed as ordinary Standard or Advanced controls.

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

Choice values are declared explicitly with `ui_choices([...])`. They are never
parsed from placeholders. The placeholder for a choice is empty before JSON or
CLI help is rendered.

## JSON contract

All JSON writers must serialize the same Fortran fields with the same meaning.
At minimum, a program record includes its exact CLI `name`, `category`,
`category_display_name`, `category_order`, `display_name`, descriptions,
executable, visibility, and its derived groups. An input record includes its
key, type, descriptions, required state, an optional display `default`,
visibility, and any applicable units, choices, group, and activation object.
Each program record also includes its requirement groups with keys and
minimum/maximum selection counts.

Every serialized program includes all seven section arrays, including `image
input/output` and `file input/output`, even when a section is empty. The CLI
uses the same lists for its headings, so a client cannot reclassify a file
independently of the command description.

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
- Put every stack or volume path in `UI_IMG`; put every other file, table, or
  directory path in `UI_FILE`. Do not use `UI_PARM` for a path-like value.
- Define shared input metadata once in `simple_ui_params_common.f90`; override
  it in a program module only when that program genuinely differs.
- Do not add GUI-only descriptor fields, client overlays, JSON patches, or
  renderer-side defaults. Add a shared semantic field only when both the CLI
  and GUI can interpret the same meaning.
- Keep declaration-time and dynamic `parameters` initialization authoritative.
  Do not add a literal UI default that duplicates a global baseline. An
  explicit program default is permitted only when that program's declared
  choice set narrows the global values and no exact routed commander default
  has been generated; it must itself be one of the declared choices.
- Keep conditional commander defaults, runtime-derived values, and fixed
  internal execution assignments distinct from baseline UI defaults.
- Keep category ownership in the program module. Do not add a second category
  table in a renderer, JSON writer, or GUI.
- Update all JSON paths when a descriptor field changes. They are duplicated
  today and must remain equivalent until serialization is consolidated.
- Validate duplicate keys, visibility values, categories, choices, and JSON
  before merging descriptor changes.

## Planned changes

1. **Consolidate JSON serialization.** Make `print_ui_json`,
   `write_ui_json`, and `ui_program%write2json` call one descriptor-to-JSON
   implementation, preserving the current JSON contract.
2. **Complete generated-default coverage and auditing.** Extend the static
   route analysis to the remaining stream-backed and helper-mediated routes.
   The current UI safely retains a program's validated choice default whenever
   a global baseline lies outside its choice set; replace that fallback with an
   exact routed default where the execution source makes one statically
   knowable. Add a review report for unmatched routes, ambiguous routes, and
   UI defaults that differ from their generated source value. Keep commander
   and parameter code out of this work.
3. **Extend requirement groups after audit.** The current groups cover
   unconditional key-presence cardinality and have replaced `alt_ios`. Add
   richer alternatives only when execution logic establishes a shared rule:
   for example, nested all-of/one-of conditions or value-dependent
   requirements. Do not restore the removed singleton `gui_exclusive_group`,
   and do not add an execution-context field unless a later audit finds an
   actual shared rule.
4. **Use activation predicates for CLI applicability validation.** The current
   binding and JSON carry `equals_any` activation data. Add validation only
   after defining how a supplied but inactive key should be reported, then
   extend the predicate form deliberately with explicit all/any/not
   composition rather than reintroducing expression strings.
5. **Review visibility by owning module.** Give every program and input a
   deliberate Standard, Advanced, or Developer value. Do not mechanically
   promote the current conservative Developer defaults.
6. **Remove obsolete listing routines.** The public listings already use
   registered category descriptors. Remove the now-unused module-local list
   printers after confirming no non-public caller remains.
7. **Review titles and descriptions by category.** Replace the temporary
   `display_name`-from-`summary` fallback with explicit concise titles, then
   review summaries, labels, help, units, and placeholders against this
   policy. A temporary exported Markdown sheet may help a domain review, but
   it is review evidence only: accepted text must be applied back to Fortran
   and the sheet must never become a source of metadata or build input.
8. **Extend descriptor validation.** Validate bindings, groups, activation
   rules, requirement-group member keys and cardinalities, choice integrity,
   user-facing text limits, compact CLI output, and JSON equivalence to the
   completed Fortran descriptors.
