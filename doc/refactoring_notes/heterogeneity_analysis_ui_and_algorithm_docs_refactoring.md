# Heterogeneity Analysis UI and Algorithm Documentation Refactoring

## Status

Deferred proposal for developer review, 2026-09-09. Implementation has not
started. This note is the single design, implementation, review, and validation
record for the change when it enters the schedule.

## Objective

Introduce a user-facing **Heterogeneity Analysis** category containing exactly:

- `flex_pca`;
- `refine3D_states`;
- `classify3D_refs`.

Align the algorithm documentation with that taxonomy by giving each of these
three major programs its own document under a heterogeneity-analysis folder.
This is a presentation and documentation refactor. Program names, command-line
contracts, execution routing, commanders, algorithms, and project artifacts
must remain unchanged.

## Current State

The UI category is assigned when a `ui_program` is registered. The same
metadata drives `simple_exec prg=list` and the complete UI JSON contract.
There is no separate handwritten program list.

The three programs are currently split across presentation categories:

| Program | Current UI module and category | Execution router |
| --- | --- | --- |
| `flex_pca` | `simple_ui_denoise`, Denoising | `simple_exec_denoise` |
| `refine3D_states` | `simple_ui_refine3D`, Refine 3D Workflows | `simple_exec_refine3D` |
| `classify3D_refs` | `simple_ui_refine3D`, Refine 3D Workflows | `simple_exec_refine3D` |

The existing UI policy gives a category its own program-construction module.
Moving only the category argument on three registrations would work
technically, but it would spread ownership of one category across unrelated UI
modules and weaken that policy.

Algorithm documentation is currently flat under `doc/algorithms`. The
continuous model is documented in `flex_pca.md`, while `refine3D_states` and
`classify3D_refs` share `heterogeneous_refinement.md`.

## Approved Target Design

### Category contract

Use one new descriptor:

```fortran
type(category_descriptor), parameter :: UI_CATEGORY = &
    category_descriptor('heterogeneity', 'Heterogeneity Analysis', 65)
```

Order `65` places the category between Refine 3D Workflows (`60`) and
Denoising (`70`) without renumbering any existing category. The identifier
`heterogeneity` is the stable machine-facing value; the heading
`Heterogeneity Analysis` is the user-facing value.

### UI ownership

Add `src/main/ui/simple/simple_ui_heterogeneity.f90`. It owns:

- the category descriptor;
- the `flex_pca`, `refine3D_states`, and `classify3D_refs` `ui_program`
  objects;
- their existing `new_*` constructors;
- `construct_heterogeneity_programs`.

Move the three constructors mechanically from `simple_ui_denoise.f90` and
`simple_ui_refine3D.f90`. Their program names, display names, summaries, help,
input bindings, defaults, requirements, ordering, and visibility must not be
edited in the relocation pass.

Import and call `construct_heterogeneity_programs` from
`simple_ui_simple_group.f90`, between the refine3D and denoise constructors.
The registry ultimately sorts categories by `category_order`, but matching the
call order to the visible order keeps the source legible.

The recursive CMake source glob already includes a new Fortran file below
`src/main`; no handwritten source list should be introduced.

### Execution ownership

Do not move or duplicate execution cases:

- `flex_pca` remains routed by `simple_exec_denoise`;
- `refine3D_states` and `classify3D_refs` remain routed by
  `simple_exec_refine3D`.

UI modules describe how programs are presented. Execution modules describe
which commander handles a program. Forcing these two taxonomies to match would
create code churn without improving the user contract.

No commander, strategy, parameter parser, generated argument source, or
scientific module belongs in this refactor.

### Algorithm documentation

Create the following structure in a separate, reviewable phase:

```text
doc/algorithms/heterogeneity_analysis/
|-- README.md
|-- flex_pca.md
|-- refine3d_states.md
`-- classify3d_refs.md
```

The folder README explains the decision boundary:

- `flex_pca` estimates continuous variability and derives discrete state maps
  from a fixed-pose latent model;
- `refine3D_states` refines same-lineage conformational states from an
  existing particle/reference scaffold;
- `classify3D_refs` classifies particles against a complete external reference
  set and then reconstructs data-derived state maps.

Move the existing `flex_pca.md` content into its program document. Split
`heterogeneous_refinement.md` into the two workflow-specific documents rather
than copying the common text. Put genuinely shared explanation in the folder
README and link to the general `refine3d` and reconstruction chapters.

Update `doc/algorithms/README.md` and all internal links. This documentation
move must not change the UI category contract and may be scheduled separately
from the Fortran relocation.

## Compatibility

The public command names and accepted arguments do not change, so scripts and
project workflows are unaffected. The complete UI JSON changes only these
fields for the three programs:

```text
category              = heterogeneity
category_display_name = Heterogeneity Analysis
category_order        = 65
```

Repository-side list and JSON generation already consume category metadata
dynamically. Before release, audit any separately distributed UI client for a
hard-coded category allowlist or persisted category identifier. Saved items
keyed by program name require no migration; state keyed by the old category
identifier may need to be refreshed or mapped.

## Implementation Sequence

1. Add `simple_ui_heterogeneity.f90` and mechanically relocate the three UI
   constructors and objects.
2. Remove their old constructor calls and register the new category constructor
   in `simple_ui_simple_group.f90`.
3. Add category regression assertions before making any wording cleanup.
4. Inspect the CLI listing and generated JSON for exact membership and order.
5. In a separate commit, reorganize and split the algorithm documentation.
6. Record completed validation and any external-client migration in this note,
   then move it to `doc/refactoring_notes/completed/`.

## Validation Criteria

Static and registry-level acceptance criteria:

- `flex_pca`, `refine3D_states`, and `classify3D_refs` are registered exactly
  once and have category `heterogeneity`, heading `Heterogeneity Analysis`, and
  order `65`;
- no other program has category `heterogeneity`;
- `refine3D` remains in Refine 3D Workflows and `icm2D` remains in Denoising;
- category identifiers and category orders remain unique within
  `simple_exec`;
- the three relocated program descriptors are otherwise unchanged;
- `simple_exec prg=list` shows the new category between Refine 3D Workflows and
  Denoising with the three programs alphabetized;
- complete UI JSON generation and validation succeed;
- command dispatch for all three names resolves through the existing routers;
- every internal algorithm-documentation link resolves after the folder move.

Add the category assertions to
`production/tests/simple_test_ui_visibility.f90`, using its existing
`assert_registered_category` helper. Compilation and runtime execution remain
user-owned under repository policy.

## Explicit Non-Goals

- Renaming any command, module outside the UI layer, or project artifact.
- Combining the three scientific algorithms or their commanders.
- Moving `flex_pca` execution into the refine3D router.
- Introducing a second category table in a renderer, JSON writer, or GUI.
- Rewording program descriptors during the mechanical UI move.
- Reorganizing unrelated algorithm-documentation categories in the same
  change.

