# Heterogeneity Analysis UI and Algorithm Documentation Refactoring

## Status

Implemented 2026-09-13. This note is the single design, implementation,
review, and validation record for the change.

## Objective

Introduce a user-facing **Heterogeneity Analysis** category containing exactly:

- `flex_pca`;
- `refine3D_states`;
- `classify3D_refs`.

Align the algorithm documentation with that taxonomy by giving each of these
three major programs its own document under a heterogeneity-analysis folder.

Remove the `ppca_volvar` program. It is a volume-variance utility that sits
next to `flex_pca` in the Denoising module, has no algorithm documentation, no
test, and no user-guide coverage, and its role is superseded by `flex_pca`.
Leaving it in Denoising after the move would invite the question of why a
variability tool is not under Heterogeneity Analysis; moving it would give the
new category a program without a documented contract.

This is a presentation, documentation, and program-removal refactor. Apart
from the removal of `ppca_volvar`, program names, command-line contracts,
execution routing, commanders, algorithms, and project artifacts must remain
unchanged.

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

The existing UI policy (`doc/policies/ui_layer_policy.md`) gives a category
its own program-construction module, with the category identifier equal to the
lowercase suffix of the module name. Moving only the category argument on
three registrations would work technically, but it would spread ownership of
one category across unrelated UI modules and violate that policy.

Both `simple_ui_denoise` and `simple_ui_refine3D` depend only on
`simple_ui_modules`; the three constructors reference no module-private
helpers, so the relocation is a mechanical cut-and-paste.

`ppca_volvar` is defined in these places and nowhere else:

| Item | Location |
| --- | --- |
| UI object, constructor `new_ppca_volvar`, registration | `src/main/ui/simple/simple_ui_denoise.f90` |
| Execution case and `xppca_volvar` commander instance | `src/main/exec/simple_exec_denoise.f90` |
| `commander_ppca_volvar` type and `exec_ppca_volvar` | `src/main/commanders/simple/simple_commanders_volops.f90` |
| `make_pcavol` (used only by `exec_ppca_volvar`) | `src/main/image_processing/simple_imgproc.f90` |
| Generated symbol indexes | `doc/code_overview/fortran-indexes/` |

Algorithm documentation is currently flat under `doc/algorithms`. The
continuous model is documented in `flex_pca.md`, while `refine3D_states` and
`classify3D_refs` share `heterogeneous_refinement.md`. The documentation site
is built by MkDocs (`mkdocs.yml`, `docs_dir: doc`, no explicit `nav`,
`use_directory_urls: false`), so file paths under `doc/` are the published
URLs.

## Proposed Target Design

### Category contract

Use one new descriptor:

```fortran
type(category_descriptor), parameter :: UI_CATEGORY = &
    category_descriptor('heterogeneity', 'Heterogeneity Analysis', 65)
```

Order `65` places the category between Refine 3D Workflows (`60`) and
Denoising (`70`) without renumbering any existing category. `simple_exec`
currently uses orders 10 to 190 in steps of 10, so 65 is free. The identifier
`heterogeneity` is the stable machine-facing value and matches the module
suffix required by the UI policy; the heading `Heterogeneity Analysis` is the
user-facing value.

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
edited in the relocation pass. The new module needs only `use simple_ui_modules`,
as the source modules do.

Import and call `construct_heterogeneity_programs` from
`src/main/ui/simple_ui_simple_group.f90`, between the refine3D and denoise
constructors. The registry ultimately sorts categories by `category_order`,
but matching the call order to the visible order keeps the source legible.

The recursive CMake source glob (`file(GLOB_RECURSE SIMPLE_src "main/*")` in
`src/CMakeLists.txt`) already includes a new Fortran file below `src/main`; no
handwritten source list should be introduced.

### Execution ownership

Do not move or duplicate execution cases:

- `flex_pca` remains routed by `simple_exec_denoise`;
- `refine3D_states` and `classify3D_refs` remain routed by
  `simple_exec_refine3D`.

UI modules describe how programs are presented. Execution modules describe
which commander handles a program. Forcing these two taxonomies to match would
create code churn without improving the user contract.

No commander, strategy, parameter parser, generated argument source, or
scientific module belongs in this refactor, except the `ppca_volvar` removal
below.

### `ppca_volvar` removal

Delete, in one commit separate from the UI relocation:

- `ppca_volvar` object, `new_ppca_volvar`, and its call in
  `construct_denoise_programs` (`simple_ui_denoise.f90`);
- the `case( 'ppca_volvar' )` branch, the `xppca_volvar` instance, and the
  `commander_ppca_volvar` import (`simple_exec_denoise.f90`);
- `commander_ppca_volvar` and `exec_ppca_volvar`
  (`simple_commanders_volops.f90`);
- `make_pcavol` in `simple_imgproc.f90`, which has no other caller.

Do not remove parameters that `ppca_volvar` shared with surviving programs
(`vol1`, `outstk`, `smpd`, `neigs`, `mskdiam`, `nthr`, `kpca_ker`,
`kpca_backend`); they are owned by other programs. Regenerate the Fortran
indexes under `doc/code_overview/fortran-indexes/` with
`scripts/gen_fortran_indexes.pl` after the deletion rather than editing them
by hand.

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

Internal links to update are exactly:

- `doc/algorithms/README.md` lines 52, 151, and 154;
- `doc/algorithms/refine3d.md` line 12;
- `doc/algorithms/flex_pca.md` line 8 (becomes a sibling link after the move).

No file outside `doc/algorithms` links to either document, including
`doc/policies/heterogeneity/refine3D_states_policy.md`.

In `doc/algorithms/README.md`, keep the numbered narrative and replace items 12
and 13 under **Heterogeneity** with three numbered items pointing into the
subfolder (`refine3d_states.md`, `classify3d_refs.md`, `flex_pca.md`), each
one sentence, plus a lead-in sentence linking the subfolder README. Do not
renumber unrelated chapters beyond the shift this causes.

Because MkDocs has no explicit `nav`, no `mkdocs.yml` change is needed for the
new folder to appear. The move does, however, change the published URLs
`algorithms/flex_pca.html` and `algorithms/heterogeneous_refinement.html`.
Leave a one-paragraph stub at each old path that links to the new location; do
not add a redirects plugin for two pages.

This documentation move must not change the UI category contract and may be
scheduled separately from the Fortran relocation.

## Compatibility

The public command names and accepted arguments of the three relocated
programs do not change, so scripts and project workflows are unaffected. The
complete UI JSON changes only these fields for the three programs:

```text
category              = heterogeneity
category_display_name = Heterogeneity Analysis
category_order        = 65
```

and drops the `ppca_volvar` descriptor. `simple_exec prg=ppca_volvar` will
fail with the standard unknown-program error; there is no deprecation alias
because the program has no documented workflow or test coverage.

Repository-side list and JSON generation already consume category metadata
dynamically. The bundled GUI client under `nice/` does not key on the category
identifier outside its tests (grep of `nice/` for `category` finds no
non-test use), and `production/stream_ui_contract.json` covers only the stream
executable, so no external-client migration is required. Saved items keyed by
program name require no migration.

## Implementation Sequence

1. Add `simple_ui_heterogeneity.f90` and mechanically relocate the three UI
   constructors and objects.
2. Remove their old constructor calls and register the new category constructor
   in `simple_ui_simple_group.f90`.
3. Add category regression assertions (see Validation Criteria).
4. Regenerate the complete UI JSON and diff it against the previous build; the
   diff must touch only the three category fields of the three programs.
   Inspect `simple_exec prg=list` for membership and order.
5. In a separate commit, remove `ppca_volvar` as listed above and regenerate
   the Fortran indexes. Regenerate and diff the UI JSON again; the diff must
   be exactly the dropped descriptor.
6. In a separate commit, reorganize and split the algorithm documentation and
   leave stubs at the old paths.
7. Record completed validation in this note, then move it to
   `doc/refactoring_notes/completed/`.

## Validation Criteria

Enforced automatically at registration by `validate_category_metadata` in
`src/main/ui/simple_ui.f90` (a violation aborts startup, so no test is needed):

- category identifiers and category orders are unique within `simple_exec`;
- heading and order are consistent for every program in a category.

Static and registry-level acceptance criteria to assert in
`production/tests/simple_test_ui_visibility.f90`:

- `flex_pca`, `refine3D_states`, and `classify3D_refs` have category
  `heterogeneity`, heading `Heterogeneity Analysis`, and order `65`
  (`assert_registered_category`, three calls);
- `refine3D` remains `refine3d` / Refine 3D Workflows / 60 (new call; the
  existing `icm2D` / `denoise` / Denoising / 70 assertion already covers the
  Denoising side);
- exactly three programs have category `heterogeneity` (new helper that walks
  the program table and counts programs per category; the existing helper
  asserts one program at a time and cannot express membership);
- `ppca_volvar` is not registered (new negative helper, or the same
  table-walk helper).

Inspection criteria, checked in step 4 and step 5 and recorded here:

- `simple_exec prg=list` shows the new category between Refine 3D Workflows
  and Denoising with the three programs in registry order
  (`classify3D_refs`, `flex_pca`, `refine3D_states`);
- complete UI JSON generation and validation succeed, and the JSON diff is
  limited to the fields named under Compatibility;
- command dispatch for all three names resolves through the existing routers;
- `simple_exec prg=ppca_volvar` reports an unknown program;
- every internal algorithm-documentation link resolves after the folder move,
  and the two stubs render.

Compilation and runtime execution remain user-owned under repository policy.

## Implementation and Validation Record

Implemented on `master` on 2026-09-13:

- added `simple_ui_heterogeneity` with the `heterogeneity` / Heterogeneity
  Analysis / 65 category contract and mechanically relocated the three
  program constructors;
- kept the three execution cases in their existing denoise and refine3D
  routers;
- removed the `ppca_volvar` UI, dispatch, commander, and caller-less
  `make_pcavol` helper;
- added registry assertions for all category boundaries, exact membership,
  and `ppca_volvar` absence;
- split the algorithm documentation into the heterogeneity-analysis folder
  and retained stubs at both old URLs;
- regenerated the code overview and Fortran indexes.

Validation observed:

- the three relocated constructor bodies compare byte-for-byte with their
  pre-move versions;
- source scans confirm the three execution cases remain in their original
  routers and no removed `ppca_volvar` implementation symbol remains;
- all local Markdown links under `doc/algorithms` resolve;
- regenerated indexes contain `simple_ui_heterogeneity` and no removed
  implementation symbol;
- `fprettify -d` parsed the edited Fortran units; warnings in the relocated
  constructors are the pre-existing overlength descriptor lines preserved by
  the mechanical move;
- `git diff --check` passes.

Per repository policy, no compilation or executable was run. The maintainer
still needs to rebuild, run `simple_test_ui_visibility`, inspect
`simple_exec prg=list`, generate and validate the complete UI JSON, confirm
the expected JSON-only category changes and descriptor removal, and verify
the unknown-program response for `ppca_volvar`.

## Explicit Non-Goals

- Renaming any command, module outside the UI layer, or project artifact.
- Combining the three scientific algorithms or their commanders.
- Moving `flex_pca` execution into the refine3D router.
- Moving `cls_split` or any other Denoising program; only `ppca_volvar` is
  affected, and it is removed, not moved.
- Introducing a second category table in a renderer, JSON writer, or GUI.
- Rewording program descriptors; any wording review is a separate note.
- Removing shared parameters or the PCA/PPCA/KPCA modules that `ppca_volvar`
  used; they have other callers.
- Reorganizing unrelated algorithm-documentation categories in the same
  change.
