# UI Descriptor Layer Refactoring

## Decision

The UI descriptor layer will remain entirely in Fortran. The Fortran objects
used by the CLI are also the objects serialized for graphical clients.

There will be no TOML catalog, Markdown catalog, Python generator, generated
program catalog, or runtime overlay. JSON is output, not an editable source.

There will also be no GUI-specific descriptor extension layer. The Fortran
descriptor is complete before registration: CLI help and validation, JSON, and
GUI rendering consume the same fields. A client must not patch, reinterpret,
or attach additional program or input metadata.

The governing rules and current class design are documented in
[`../policies/ui_layer_policy.md`](../policies/ui_layer_policy.md).

## Outcome

The refactoring will provide:

- one Fortran source of truth for CLI and GUI descriptors;
- a stable category on every program, derived from the owning
  `simple_ui_*` or `single_ui_*` module;
- Standard, Advanced, or Developer visibility on every program and input;
- JSON that directly mirrors the in-memory Fortran descriptor;
- baseline defaults exported from the existing `parameters` initialization to
  UI descriptors;
- structured defaults, units, and choices instead of duplicating or
  extracting those values from display text;

Program names, CLI keys, executable routing, required-key behavior, project
requirements, and scientific execution are not being changed by this work.

## Text standardization on `ui_refac`

Placeholders are standardized by the Fortran descriptor constructor, so every
rendered input has a short type-specific hint and choice controls have no
placeholder. The legacy placeholder argument is still consumed internally to
extract existing choice values into `ui_param%choices`; it is not rendered.

Program summaries now have a 100-character enforced upper bound. Production
summaries shorter than 30 characters, and the prior overlong summary, were
rewritten to state the program action and primary outcome using the current
program help as the implementation reference. Developer-only test summaries
remain technical until their separate UI review. The policy defines the
required summary and placeholder formats in detail.

## Fields

### Program fields

The target program metadata is:

| Field | Status | Purpose |
| --- | --- | --- |
| `name` | Existing | Exact CLI program name. |
| `category` | Added | Stable owning-module category identifier. |
| `category_display_name`, `category_order` | Added | Shared category heading and display order. |
| `display_name` | Added | Plain-English GUI title, separate from the CLI name. |
| `summary` | Renamed | Short plain-English program summary. |
| `help` | Renamed | Full program help. |
| `executable` | Existing | Owning executable. |
| `visibility` | Added | Standard, Advanced, or Developer. |
| `sp_required` | Existing | Project requirement used by the CLI. |
| seven input lists | Existing | Existing help and input sections. |

`display_name` is stored directly on `ui_program` and serialized from that
object. It is never supplied by a catalog. The constructor currently copies a
program's reviewed `summary` when an explicit `display_name=` is not yet
present, ensuring every existing GUI record has a plain-English title without
changing a CLI identifier. Replace that migration fallback with concise,
explicit titles category by category.

### Input fields

The target input metadata is:

| Field | Status | Purpose |
| --- | --- | --- |
| `key`, `keytype` | Existing | Exact CLI identity and type. |
| `label` | Renamed | Short plain-English field label. |
| `help` | Renamed | Full help. |
| `placeholder` | Renamed, to simplify | Short entry example only. |
| `required` | Existing | CLI requiredness. |
| `default` | To add to the UI export | One optional display value, sourced from `parameters`. |
| `has_default` | Existing internal field | Controls whether the current descriptor emits a default. |
| `cval_default`, `rval_default` | Existing internal storage | Compatibility storage until the default export is simplified. |
| `units` | Added | Units separate from descriptions. |
| `choices` | Added | Structured CLI values and display text. |
| `visibility` | Added | Standard, Advanced, or Developer. |
| legacy `gui_*` fields | To retire | GUI-specific submenu, activation, exclusion, online, and visibility override machinery. |

Choice values are currently initialized from the legacy placeholder once,
inside the Fortran descriptor construction path. JSON now consumes `choices`;
it must not parse the placeholder independently.

### Descriptor terminology

The current API and JSON use these presentation names:

| Scope | Current fields |
| --- | --- |
| Input | `label`, `help`, `placeholder` |
| Program | `summary`, `help` |
| Program help method | `print_help` |

These names are semantic only: they do not change program membership, CLI
keys, descriptions, or execution behavior.

### Shared descriptor semantics

The target model does not preserve a separate layer of GUI behavior. The
following legacy mechanisms require migration rather than expansion:

| Legacy mechanism | Shared replacement |
| --- | --- |
| `gui_visibility` | `visibility` on the program or input. |
| `gui_submenu` | A neutral input group, also used by CLI help and JSON. |
| `gui_active_flags` | A declarative activation condition that CLI validation and GUI rendering evaluate consistently. |
| `gui_exclusive_group` and `alt_ios` | A program-owned input-requirement group with explicit minimum/maximum selection counts. |
| `gui_online` | A declared execution-context capability, if it remains needed after audit. |
| `apply_gui_overrides` | Final shared descriptor data declared in the owning program module. |

Do not add another GUI-prefixed field or a renderer-specific override. Reuse is
still desirable, but common parameter descriptors must resolve to complete
program-specific descriptors before registration.

## Category design

Categories use the organization that already exists in the source tree.
Every program-defining module declares one constant:

```fortran
module simple_ui_denoise
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = &
    category_descriptor('denoise', 'Denoising', 70)
```

Every program registration passes it:

```fortran
call add_ui_program('icm2D', icm2D, prgtab, UI_CATEGORY)
```

`add_ui_program` requires complete category metadata, assigns it to the
`ui_program`, and then registers the program. This prevents new registrations
from silently omitting an identifier, heading, or order.

The convention supports expansion without a closed enumeration. Add a new
category by adding a new `simple_ui_*` or `single_ui_*` module and wiring its
constructor into the correct executable group. Do not maintain a parallel
category catalog.

The public CLI listing paths now traverse the registry and group programs by
their registered category descriptor. JSON serializes the same category
identifier, heading, and order. The obsolete handwritten heading routines are
no longer used by any public listing path.

## Default-value design

### Current duplication

Defaults are currently represented independently in several places:

- declaration-time component initialization in `type(parameters)`;
- `init_dynamic_defaults` for dynamic string components;
- conditional `cline%set` calls in commanders, strategies, and executable
  entry points; and
- literal default arguments and `{default}` placeholder fragments in UI
  constructors.

These values can differ by design. For example, `autoscale` initializes to
`no` in `parameters`, while `cluster2D` conditionally supplies `yes` when the
user omits the key. The GUI must not duplicate either value in a second source
of truth. It publishes the initialization as its baseline and leaves an
untouched optional key out of the submitted command line, allowing the
program-specific behavior to remain intact.

### Target ownership

Keep the existing declaration-time values in `type(parameters)` and the
dynamic string values in `init_dynamic_defaults`. They are the only source of
baseline defaults. Do not move them to a catalog, generator, or second Fortran
table.

Add a read-only default-snapshot API in `src/main/params/`. It must construct
only the base `parameters` state and dynamic defaults; it must not parse a
command line, inspect files or a project, derive settings, or execute a
commander. Extend the existing typed parameter registry as needed to look up a
snapshot value by CLI key without creating another key-to-field table. The
snapshot API returns the canonical CLI text representation needed for display;
the GUI does not need typed default metadata.

After UI programs are constructed, the GUI/JSON bootstrap path will copy that
snapshot into each optional `ui_param`. This is a one-way export from
`parameters` to the UI. Existing literal UI defaults become temporary
assertions during migration and are removed after the snapshot agrees with
them.

Conditional commander and strategy defaults remain in their owning code. They
are program-specific execution behavior, not part of the baseline UI snapshot.
For example, `cluster2D` may select `autoscale=yes` when the key is absent even
though the baseline `parameters` value is `no`.

The GUI must submit only user-changed optional values. An untouched displayed
baseline value is omitted, allowing the existing commander default to take
effect. This avoids changing CLI behavior merely because a GUI rendered a
field before program setup.

### UI and JSON representation

An input either emits one optional JSON `default` display value or omits it.
The exported value may be text for every key type; normal CLI parsing remains
responsible for converting a submitted value. Do not add `default_kind`,
runtime/default provenance, or a second GUI-specific default table.

Placeholders become examples or format hints only. They must not contain
`{default}` syntax.

### Migration sequence

1. Add the read-only `parameters` default-snapshot and typed lookup API without
   changing parameter initialization or execution behavior.
2. Verify the snapshot against a default-constructed `parameters` object plus
   `init_dynamic_defaults`.
3. Add one UI/JSON bootstrap synchronization point after `make_ui` and before
   UI consumption. Do not create a UI-to-parameters module dependency.
4. Make `ui_program%add_input` use snapshot values. Treat an old literal UI
   default as a temporary assertion, then remove it.
5. Omit a display default when no sensible baseline is available, and keep
   forced internal values outside the UI default contract.
6. Ensure GUI submission omits unchanged optional inputs so commander defaults
   continue to work.
7. Remove default fragments from placeholders and remove the old literal
   default constructor arguments after validation.

This stays entirely in Fortran and does not introduce a catalog or generator.

## Implementation status on `ui_refac`

Completed in the current working change:

- added `ui_program%category`, category accessors, and registration
  validation;
- added one `UI_CATEGORY` constant to all 40 program-defining
  `simple_ui_*`, `single_ui_*`, and test UI modules;
- changed those constants into module-owned descriptors carrying the category
  identifier, heading, and display order;
- updated all 241 `add_ui_program` calls to pass their module category;
- replaced public CLI listing calls with registry traversal grouped by the
  registered category descriptor;
- added category heading and order to all program JSON serialization paths;
- retained the Fortran visibility type and descriptor support types;
- preserved Developer as the conservative default for unspecified
  descriptors during migration;
- removed the legacy Boolean `advanced` flag; `visibility` is now the sole
  presentation classification;
- added program categories to all three JSON serialization paths;
- changed JSON option export to consume `ui_param%choices`;
- added `has_default` and units to JSON and stopped emitting invented defaults.
- renamed legacy `descr_*` fields and JSON members to input `label`, `help`,
  and `placeholder`, and program `summary` and `help`.
- added `display_name` to `ui_program`, all JSON paths, program information,
  descriptor validation, and the NICE Lite batch-program chooser. Existing
  programs receive their reviewed summary as an initial title until an owning
  module supplies a shorter explicit `display_name=`.

Validation completed on macOS with GNU Fortran 15.2.0:

- CMake reconfigured successfully after the catalog sources were removed;
- the SIMPLE library and `simple_exec`, `single_exec`, `simple_stream`,
  `simple_test_exec`, and `simple_private_exec` compiled;
- `simple_test_ui_visibility` passed all 40 assertions, including
  representative SIMPLE, SINGLE, stream, and test categories;
- a full `simple_ui.json` was generated and parsed with `jq`;
- all serialized production programs had non-empty categories and valid
  visibility names;
- all serialized production inputs had valid visibility names;
- representative choice/default serialization was checked;
- all four executable listing paths completed successfully.
- `simple_exec` and `simple_private_exec` rebuilt successfully after the
  Boolean visibility cleanup; the JSON emitted by `print_ui_json` contains no
  `advanced` member.

The current conservative production baseline is 53 Standard and 103
Developer programs. No program is explicitly Advanced yet. Inputs currently
contain 100 Standard, 33 Advanced, and 1166 Developer descriptors. The
category-by-category visibility review must decide which Developer programs
and inputs should move to Advanced; this must not be guessed mechanically.

Linux/BOX tests have not been run.

## Next refactoring steps

Resume in this order:

1. **Export baseline parameter defaults.** Implement the read-only snapshot
   and UI synchronization described above before further serializer or
   description work. Preserve `parameters` initialization and commander
   defaults exactly as they are.
2. **Consolidate JSON serialization.** `print_ui_json`, `write_ui_json`, and
   `ui_program%write2json` still contain closely related implementations.
   Move descriptor-to-JSON behavior into one serializer and make all three
   entry points call it.
3. **Replace GUI-specific metadata with shared semantics.** Rename
   `gui_visibility` to `visibility`; introduce neutral grouping, activation,
   execution-context, and input-requirement semantics only where their CLI
   meaning is defined. Migrate one category at a time, then remove
   `apply_gui_overrides`, `gui_*` fields, and `alt_ios`.
4. **Review visibility module by module.** Replace inherited defaults with
   deliberate `visibility=` values for each program and input. Work one
   category per change so domain owners can review the choices.
5. **Make choices fully explicit in Fortran.** Add a concise constructor API
   for choices, migrate binary and multiple-choice inputs, and then stop
   decoding choice values from `placeholder`.
6. **Remove obsolete listing routines.** The public listing paths now use the
   registry. Remove unused module-local `print_*_programs` routines after
   confirming no non-public caller depends on them.
7. **Review display names and descriptions by category.** Replace the
   migration title fallback with concise explicit `display_name=` values, then
   apply the workflow below to summaries, labels, and help. Keep wording-only
   changes separate from structural changes.
8. **Add descriptor validation.** Report missing categories, invalid
   visibility, malformed or duplicate choices, empty user-facing text, and
   length-limit violations.

## Description update workflow

Descriptions remain in their owning Fortran modules. To make large wording
reviews manageable without creating another source of truth:

1. Select one category/module, for example `simple_ui_denoise.f90`.
2. Export a temporary Markdown review sheet containing program/input identity,
   current summary or label, help, placeholder, units, choices, and visibility.
3. Edit and review that temporary sheet with the domain owner.
4. Apply the accepted text back to the Fortran constructors.
5. Regenerate the sheet and compare it with the reviewed version.
6. Delete or archive the sheet as review evidence; never compile or load it at
   runtime.
7. Compile the module, run descriptor validation, and inspect its JSON.

The exporter/importer, if created, is an editing aid only. It must validate
program names and input keys against the current Fortran source, refuse
unknown or duplicate records, show an explicit diff, and never become part of
the build. The reviewed Fortran diff is the authoritative change.

For each description and title:

- preserve exact program names, CLI keys, and accepted values;
- give every production program a concise `display_name` distinct from the
  CLI identifier; use the policy's title rules and retain established
  scientific acronyms where useful;
- keep labels and summaries near 60 characters or fewer;
- keep placeholders near 40 characters or fewer;
- use placeholders only for examples or formats;
- move units into `units`;
- move choice values into `choices`;
- put explanations and constraints in `help`;
- avoid internal variable names and implementation terminology.

## Validation strategy

Structural validation should check:

- every registration has a non-empty category;
- every visibility value is valid;
- program names remain unique;
- input keys remain unique within their intended scope;
- required/default state is internally consistent;
- each emitted UI default equals the textual `parameters` baseline snapshot
  for its key;
- a key without a sensible baseline value omits `default`;
- GUI submission omits unchanged optional inputs, allowing commander defaults
  to apply;
- forced internal assignments are not mislabeled as baseline UI defaults;
- binary inputs have two choices and multiple-choice inputs have at least two;
- JSON categories, visibility, defaults, units, and choices match Fortran.

Migration-time wording validation should report, rather than immediately fail,
labels or summaries over 60 characters and placeholders over 40 characters.
These reports provide the worklist for category-by-category cleanup.
