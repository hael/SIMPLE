# UI Category Membership Tidying

## Status

Implemented 2026-09-13 (proposed, extended with the Reconstruct 3D and
Post-processing categories and the `export_manifoldem_starproject` removal,
and implemented the same day). This note is the single design,
implementation, review, and validation record for the change. See
"Implementation and Validation Record" at the end.

## Objective

Four things, each in its own commit:

1. Introduce a user-facing **Reconstruct 3D Workflows** category containing
   exactly `reconstruct3D` and `bootstrap_rec3D`, which are currently listed
   under Refine 3D Workflows.
2. Introduce a user-facing **Post-processing** category containing exactly
   `postprocess` and `postprocess_nu`, currently split between Refine 3D
   Workflows and Filtering.
3. Relocate six `simple_exec` programs whose UI category disagrees with what
   they do, so that every program is found where a user would look for it.
4. Remove the `export_manifoldem_starproject` program.

The categories and the relocations follow the pattern established by
`../completed/heterogeneity_analysis_ui_and_algorithm_docs_refactoring.md`:
move the Fortran construction to the module that owns the target category,
leave execution routing and commanders untouched, and assert the result in
`simple_test_ui_visibility`. The removal follows the `ppca_volvar` removal in
that same note.

Apart from the removal, program names, command-line contracts, execution
routing, commanders, algorithms, and project artifacts remain unchanged.
Descriptor wording is not touched except for one summary that is factually
wrong (`automask`); the broader wording pass is a separate note.

## Current State

The category survey behind this note (all 19 `simple_exec` categories, from
the generated complete UI JSON) found the taxonomy sound overall, with these
misplacements:

| Program | Current UI module (category) | Target UI module (category) | Execution router (unchanged) | Commander module (unchanged) |
| --- | --- | --- | --- | --- |
| `reconstruct3D` | `simple_ui_refine3D` (Refine 3D Workflows) | `simple_ui_reconstruct3D` (Reconstruct 3D Workflows, new) | `simple_exec_refine3D` | `simple_commanders_refine3D` |
| `bootstrap_rec3D` | `simple_ui_refine3D` (Refine 3D Workflows) | `simple_ui_reconstruct3D` (Reconstruct 3D Workflows, new) | `simple_exec_refine3D` | `simple_commanders_refine3D` |
| `postprocess` | `simple_ui_refine3D` (Refine 3D Workflows) | `simple_ui_postprocess` (Post-processing, new) | `simple_exec_refine3D` | `simple_commanders_volops` |
| `postprocess_nu` | `simple_ui_filter` (Filtering) | `simple_ui_postprocess` (Post-processing, new) | `simple_exec_filter` | `simple_commanders_postprocess_nu` |
| `automask` | `simple_ui_refine3D` (Refine 3D Workflows) | `simple_ui_mask` (Masking) | `simple_exec_refine3D` | `simple_commanders_mask` |
| `ptcl3D_state_consensus` | `simple_ui_project` (Project Management) | `simple_ui_heterogeneity` (Heterogeneity Analysis) | `simple_exec_project` | `simple_commanders_project_core` |
| `cls_split` | `simple_ui_denoise` (Denoising) | `simple_ui_cluster2D` (Cluster2D Workflows) | `simple_exec_denoise` | `simple_commanders_denoise` |
| `fractionate_movies` | `simple_ui_other` (Other Utilities) | `simple_ui_preproc` (Pre-processing) | `simple_exec_other` | `simple_commanders_misc` |
| `split` | `simple_ui_other` (Other Utilities) | `simple_ui_image` (General Image Processing) | `simple_exec_other` | `simple_commanders_distr` |
| `split_stack` | `simple_ui_other` (Other Utilities) | `simple_ui_image` (General Image Processing) | `simple_exec_other` | `simple_commanders_project_ptcl` |

Rationale per program:

- `reconstruct3D` and `bootstrap_rec3D` are reconstruction workflows, not
  refinement: `reconstruct3D` builds volumes from already-oriented particles
  and is invoked internally by `abinitio3D`, `refine3D`, and
  `refine3D_states`; `bootstrap_rec3D` is the complete final-reconstruction
  sequence (sigma2 seeding plus ML-regularized reconstruction). Listing them
  under Refine 3D Workflows hides the step users run after refinement or
  heterogeneity analysis. `sigma2_convert` and `center` were considered and
  deliberately left where they are.
- `postprocess` ("Filter and sharpen a reconstructed density map for
  interpretation") and `postprocess_nu` ("Nonuniform evidence-bounded
  postprocessing of even/odd half-maps") are the uniform and nonuniform
  variants of the same final step, yet a user finds one under Refine 3D and
  the other under Filtering. `nu_filt3D` was considered and left in
  Filtering: it is a local low-pass filter applied by name, not a
  sharpening/interpretation step.
- `automask` performs envelope masking (`exec_automask` calls
  `mskvol%automask3D` and writes `automask3D_masked_vol.mrc`); its commander
  already lives in the masking module and its help text says "automated
  envelope masking". Only its UI registration and its summary ("Create a
  spherical mask from the estimated particle diameter", which describes
  `auto_spher_mask`, not this program) place it elsewhere. Masking rather
  than Post-processing because it is also used ahead of refinement.
- `ptcl3D_state_consensus` builds a consensus particle-state assignment from a
  file table of projects, i.e. it post-processes `refine3D_states` /
  `classify3D_refs` runs. It is the fourth heterogeneity-analysis program.
- `cls_split` ("Split classes with latent clustering") splits 2D/3D particle
  classes into subclasses by diffusion-map or kPCA embedding and k-medoids,
  then regenerates class averages through `make_cavgs`. It operates on class
  memberships, like `sample_classes` and `bootstrap_cavgs`, and has nothing
  to do with denoising beyond sharing the PCA machinery; its commander's
  module name is historical. `reimport_particles` ("Re-import denoised
  particle stack") was briefly considered for Denoising during review and
  rejected: it replaces the project particle stack while preserving
  particle/class metadata, and a denoised stack is only its typical input.
  It stays in Project Management.
- `fractionate_movies` re-generates micrographs from selected movie frames;
  it is a pre-processing step with a distributed workflow, like
  `motion_correct`.
- `split` and `split_stack` are stack utilities; their counterpart `stack`
  ("Combine image files or stacks into one stack") is in General Image
  Processing.

Note that the `split` UI object is named `split_` in `simple_ui_other.f90`
(`type(ui_program), target :: split_`, constructor `new_split_`); the
relocation must keep that spelling.

`export_manifoldem_starproject` exports a `particles3D_manifoldem.star` file
for ManifoldEM. It is defined in these places and nowhere else:

| Item | Location |
| --- | --- |
| UI object, constructor `new_export_manifoldem_starproject`, registration | `src/main/ui/simple/simple_ui_project.f90` |
| Execution case, `xexport_manifoldem_starproject` instance, and the commander import | `src/main/exec/simple_exec_project.f90` |
| `commander_export_manifoldem_starproject` and `exec_export_manifoldem_starproject` | `src/main/commanders/simple/simple_commanders_starproject.f90` |
| `export_manifoldem_ptcls3D` (type-bound procedure, no other caller) | `src/main/star/simple_starproject.f90` |
| Generated symbol indexes | `doc/code_overview/fortran-indexes/` |

No test, script, algorithm chapter, or user guide references it; the only
prose mentions are in a rejected refactoring note and a completed
implementation note, which are historical records and stay as they are.

Two further findings from the survey are recorded here as decisions, not as
part of this change:

- `extract_subproj` ("extraction of a subproject of time-series of metallic
  nanoparticles") is a nanoparticle time-series program registered under
  `simple_exec` Project Management. The natural home is the `single_exec`
  Time-series Pre-processing category, but that changes the executable and is
  therefore a command-line contract change, out of scope for a mechanical
  move. Decide separately whether to move it (with a deprecation alias) or
  leave it.
- `mini_stream` ("standalone mini_stream for a quick look") sits in
  Validation. No other `simple_exec` category fits better, and the stream
  programs live in the separate `simple_stream` executable. Leave it; fix its
  summary in the wording pass.

Near-duplicate pairs were checked and are not duplicates: `automask`
(envelope) and `auto_spher_mask` (spherical) differ; `noisevol` (volumes,
`nstates`; called internally by `simple_commanders_abinitio2D`) and
`simulate_noise` (2D images, `nptcls`; its summary "images or volumes" is
wrong) differ. Both pairs stay.

## Proposed Target Design

### New category contracts

Two new descriptors, one per new module:

```fortran
! src/main/ui/simple/simple_ui_reconstruct3D.f90
type(category_descriptor), parameter :: UI_CATEGORY = &
    category_descriptor('reconstruct3d', 'Reconstruct 3D Workflows', 68)

! src/main/ui/simple/simple_ui_postprocess.f90
type(category_descriptor), parameter :: UI_CATEGORY = &
    category_descriptor('postprocess', 'Post-processing', 69)
```

Orders `68` and `69` sit after Heterogeneity Analysis (`65`) and before
Denoising (`70`), so the listing reads Refine 3D → Heterogeneity Analysis →
Reconstruct 3D → Post-processing → Denoising: final maps, then their
interpretation, after either refinement path. `simple_exec` orders 10 to 190
in steps of 10 plus 65 are taken; 68 and 69 are free. Both identifiers are
the lowercase suffix of their module name, as the UI policy requires
(compare `simple_ui_refine3D` → `refine3d`). The category identifier
`postprocess` coincides with the program name `postprocess`; the two are
separate namespaces (`ui_program%category` versus the `ui_hash` key) and
nothing in the registry, JSON writer, or validator conflates them. The
heading "Reconstruct 3D Workflows" parallels "Refine 3D Workflows"; "3D
Reconstruction" was avoided because `single_exec` already uses that heading
for `nano3d`. "Post-processing" parallels "Pre-processing" (`preproc`, 20).

Add `src/main/ui/simple/simple_ui_reconstruct3D.f90` owning its descriptor,
the `reconstruct3D` and `bootstrap_rec3D` `ui_program` objects, their
existing `new_*` constructors moved mechanically from
`simple_ui_refine3D.f90`, and `construct_reconstruct3D_programs`. Add
`src/main/ui/simple/simple_ui_postprocess.f90` owning its descriptor, the
`postprocess` and `postprocess_nu` objects, their constructors moved
mechanically from `simple_ui_refine3D.f90` and `simple_ui_filter.f90`, and
`construct_postprocess_programs`. Each module needs only
`use simple_ui_modules`.

Import and call both constructors from
`src/main/ui/simple_ui_simple_group.f90` immediately after
`construct_heterogeneity_programs`, in the order reconstruct3D then
postprocess, matching the visible order. The recursive CMake glob picks up
the new files; no source list is edited.

### UI ownership for the relocations

For each remaining row of the table, cut the
`type(ui_program), target :: <name>` declaration, the `new_<name>`
constructor, and the `call new_<name>(prgtab)` line from the source module
and paste them into the target module. Each constructor's final
`add_ui_program(...)` call already passes the module-local `UI_CATEGORY`, so
the category follows the module automatically. Program names, display names,
summaries, help, input bindings, defaults, requirements, ordering, and
visibility must not be edited in the relocation pass, with the single
exception below.

All source and target modules depend only on `simple_ui_modules`; no
module-private helper is referenced by any of the moved constructors, so no
`use` statement changes.

Place each pasted constructor call at the end of the target module's
`construct_*_programs` list. Registry order is not display order (programs
are listed alphabetically within a category), so this is a legibility choice
only.

Apart from the two new lines, `simple_ui_simple_group.f90` needs no change:
every other target module is already imported and constructed there.

### The one wording change

Replace the `automask` summary

```text
Create a spherical mask from the estimated particle diameter
```

with

```text
Create an automatic envelope mask for a 3D volume
```

in the same commit as its relocation, because the old summary is false and
would become more visible next to `auto_spher_mask` in the Masking listing.
`display_name` falls back to `summary`, so this also corrects the GUI title.
No other descriptor text changes.

### Execution ownership

Nothing moves. `simple_exec_refine3D`, `simple_exec_filter`,
`simple_exec_denoise`, `simple_exec_project`, and `simple_exec_other` keep
their `case` branches and commander instances. After the change `simple_exec_refine3D` routes
programs presented under four headings (Refine 3D, Reconstruct 3D,
Post-processing, Masking); that asymmetry is expected under the
UI-versus-execution split and is not a reason to touch the routers. The
internal `reconstruct3D` invocations from `abinitio3D`, `refine3D`, and
`refine3D_states` go through `cline%set('prg', ...)` and are unaffected by
presentation.

### `export_manifoldem_starproject` removal

Delete, in one commit separate from the moves:

- the `export_manifoldem_starproject` object, `new_export_manifoldem_starproject`,
  and its call in `construct_project_programs` (`simple_ui_project.f90`);
- the `case( 'export_manifoldem_starproject' )` branch and the
  `xexport_manifoldem_starproject` instance (`simple_exec_project.f90`). The
  router currently imports `simple_commanders_starproject` twice, once with
  the ManifoldEM commander and once without; delete the two-line import that
  names it and keep the one-line import, which then covers everything the
  router still uses;
- `commander_export_manifoldem_starproject` and
  `exec_export_manifoldem_starproject` (`simple_commanders_starproject.f90`);
- the `export_manifoldem_ptcls3D` type-bound procedure and its binding in the
  `starproject` type (`simple_starproject.f90`). The helpers it calls
  (`initialise`, `propagate_optics`, `propagate_optics_box`,
  `get_stkname_and_ind`) have other callers and stay.

Regenerate the Fortran indexes under `doc/code_overview/fortran-indexes/`
with `scripts/gen_fortran_indexes.pl` rather than editing them by hand.
`simple_exec prg=export_manifoldem_starproject` will then fail with the
standard unknown-program error; there is no deprecation alias because the
program has no documented workflow or test coverage. `export_starproject`
and `export_relion` remain the supported STAR exports.

### Heterogeneity Analysis membership

The completed heterogeneity note fixed the category at exactly three
programs. This note supersedes that count: after the move the category holds
`classify3D_refs`, `flex_pca`, `ptcl3D_state_consensus`, `refine3D_states`.
Update the `count_prgs_in_category('heterogeneity')` assertion from 3 to 4
and add one line to `doc/algorithms/heterogeneity_analysis/README.md` noting
that `ptcl3D_state_consensus` combines state assignments across runs (no
algorithm chapter is needed for a metadata utility).

## Compatibility

Command names and accepted arguments of the surviving programs are
unchanged. The complete UI JSON changes `category`,
`category_display_name`, and `category_order` for the ten moved programs,
`summary` and `display_name` for `automask`, and drops the
`export_manifoldem_starproject` descriptor. The external-client audit from
the heterogeneity note still holds: the bundled `nice/` client does not key
on category identifiers.

## Implementation Sequence

1. Add `simple_ui_reconstruct3D.f90`, relocate `reconstruct3D` and
   `bootstrap_rec3D` into it, and register the constructor in
   `simple_ui_simple_group.f90`.
2. Add `simple_ui_postprocess.f90`, relocate `postprocess` and
   `postprocess_nu` into it, and register the constructor.
3. Relocate the other six constructors and objects; fix the `automask`
   summary.
4. Extend `simple_test_ui_visibility` (see Validation Criteria).
5. Regenerate the complete UI JSON and diff it against the previous build;
   the diff must be limited to the fields named under Compatibility, minus
   the dropped descriptor.
6. Inspect `simple_exec prg=list` for the two new headings and the ten
   programs under their new headings.
7. In a separate commit, remove `export_manifoldem_starproject` as listed
   above. Regenerate and diff the UI JSON again; the diff must be exactly the
   dropped descriptor.
8. Update the heterogeneity README line and regenerate the Fortran indexes.
9. Record completed validation here, then move this note to
   `doc/refactoring_notes/completed/`.

## Validation Criteria

Registry assertions in `production/tests/simple_test_ui_visibility.f90`,
using the existing helpers:

- `assert_registered_category` for each moved program with its new
  identifier, heading, and order: `reconstruct3D` and `bootstrap_rec3D` →
  `reconstruct3d` / Reconstruct 3D Workflows / 68; `postprocess` and
  `postprocess_nu` → `postprocess` / Post-processing / 69; `automask` →
  `mask` / Masking / 100; `ptcl3D_state_consensus` → `heterogeneity` /
  Heterogeneity Analysis / 65; `fractionate_movies` → `preproc` / Pre-processing / 20; `split` and
  `split_stack` → `image` / General Image Processing / 90; `cls_split` →
  `cluster2d` / Cluster2D Workflows / 30;
- `count_prgs_in_category('reconstruct3d')` equals 2,
  `count_prgs_in_category('postprocess')` equals 2, and
  `count_prgs_in_category('heterogeneity')` equals 4;
- `assert_program_not_registered('export_manifoldem_starproject')`;
- the anchors `refine3D` (60), `icm2D` (70), `filter` (80), `new_project`
  (10), `reimport_particles` (10), `stack` (90), `motion_correct` (20),
  `abinitio2D` (30), and `export_starproject` (10) keep their categories, so a mistaken move or deletion in the other
  direction is caught.

Inspection criteria, checked in steps 5 to 7 and recorded here:

- the JSON diffs are limited to the fields named under Compatibility;
- the two new headings appear between Heterogeneity Analysis and Denoising
  in that order, the ten programs appear under their new headings in
  `prg=list`, and no heading has gone empty (`Other Utilities` retains
  `cif2pdb` and `sigma2_convert`; `Refine 3D Workflows` retains `refine3D`
  and `refine3D_auto`; `Filtering` retains `filter`, `nu_filt3D`,
  `uniform_filter2D`, `uniform_filter3D`; `Denoising` retains `icm2D`,
  `icm3D`, `ppca_denoise`, `ppca_denoise_classes`, `denoise_project`,
  `map_params_from_den`);
- command dispatch for the ten moved names resolves through the existing
  routers;
- `simple_exec prg=export_manifoldem_starproject` reports an unknown
  program.

Category-order uniqueness and per-category metadata consistency are enforced
at registration by `validate_category_metadata`; no test is needed.
Compilation and runtime execution remain user-owned under repository policy.

## Explicit Non-Goals

- Adding `sigma2_convert` or `center` to Reconstruct 3D Workflows, or
  `nu_filt3D` or `automask` to Post-processing; each new category holds
  exactly its two named programs.
- Moving `extract_subproj` between executables, or moving `mini_stream`.
- Removing or merging `noisevol` / `simulate_noise` or `automask` /
  `auto_spher_mask`.
- Touching execution routers, commanders, or the `split_` object name.
- Removing `export_starproject`, `export_relion`, or any STAR helper that has
  a surviving caller.
- Any descriptor wording other than the `automask` summary. The wording pass
  (copy-pasted `gen_pspecs_and_thumbs` summary; identical summaries on
  `cluster_stack` / `match_cavgs` / `match_stacks`; `volanalyze` "emsemble";
  `simulate_noise` "images or volumes"; lowercase summaries on
  `extract_subproj`, `prune_project`, `replace_project_field`,
  `bootstrap_rec3D`, `auto_spher_mask`, `check_refpick`, `mini_stream`) is a
  separate note under `doc/policies/ui_layer_policy.md`'s category-review
  rule.
- Renumbering existing categories.

## Implementation and Validation Record

Implemented on `master` on 2026-09-13, in the working tree alongside the
heterogeneity-analysis refactor (not yet committed; the sequence above still
applies when committing).

Source changes:

- Added `src/main/ui/simple/simple_ui_reconstruct3D.f90`
  (`reconstruct3d` / Reconstruct 3D Workflows / 68) holding `reconstruct3D`
  and `bootstrap_rec3D`, and `src/main/ui/simple/simple_ui_postprocess.f90`
  (`postprocess` / Post-processing / 69) holding `postprocess` and
  `postprocess_nu`; both constructors registered in
  `simple_ui_simple_group.f90` directly after
  `construct_heterogeneity_programs`.
- Relocated `automask` → `simple_ui_mask`, `ptcl3D_state_consensus` →
  `simple_ui_heterogeneity`, `cls_split` → `simple_ui_cluster2D`,
  `fractionate_movies` → `simple_ui_preproc`, `split_` and `split_stack` →
  `simple_ui_image`. (`reimport_particles` was moved to `simple_ui_denoise`
  in a first pass and moved back to `simple_ui_project` on review; its
  category is unchanged, but its declaration, constructor call, and
  constructor now sit last in that module rather than at their original
  positions.) Each pasted constructor call sits at the end of its
  target `construct_*_programs` list.
- `automask` summary changed to "Create an automatic envelope mask for a 3D
  volume".
- Removed `export_manifoldem_starproject` from `simple_ui_project.f90`,
  `simple_exec_project.f90` (the duplicated `simple_commanders_starproject`
  import collapsed to the one-line form), `simple_commanders_starproject.f90`,
  and the `export_manifoldem_ptcls3D` binding and body in
  `simple_starproject.f90`. No other symbol lost its last caller.
- `simple_test_ui_visibility.f90`: category assertions for the ten moved
  programs, counts for `reconstruct3d` (2), `postprocess` (2), and
  `heterogeneity` (4), negative assertion for
  `export_manifoldem_starproject`, and anchors `filter`, `new_project`,
  `export_starproject`, `stack`, `motion_correct` alongside the existing
  `refine3D` and `icm2D`.
- `doc/algorithms/heterogeneity_analysis/README.md`: one bullet for
  `ptcl3D_state_consensus`. `doc/code_overview/code_base_map.md`: the two new
  modules. Fortran indexes regenerated with
  `perl scripts/gen_fortran_indexes.pl --root src --out doc/code_overview/fortran-indexes`.

Static validation performed (no compilation or executable run, per
repository policy):

- Nine of the ten relocated constructor blocks are byte-identical to their
  pre-move definitions. `automask` differs in exactly two lines: the summary
  replacement described above, and its `subroutine new_automask( prgtab )`
  header, which was unindented in `simple_ui_refine3D.f90` and now carries
  the standard four-space indent. One trailing-whitespace comment line inside
  the moved `split_` block was stripped so that `git diff --check` passes.
- No moved constructor references a module-level object of its source
  module; every touched UI module still depends only on
  `simple_ui_modules`; no duplicate object or constructor names in any UI
  module; subroutine/end-subroutine and module/end-module pairs balance in
  every edited unit.
- Every moved program still has its `case` branch in its original router
  (`simple_exec_refine3D`, `simple_exec_filter`, `simple_exec_denoise`,
  `simple_exec_project`, `simple_exec_other`).
- Repository-wide scan finds no remaining reference to
  `export_manifoldem_starproject`, `commander_export_manifoldem_starproject`,
  or `export_manifoldem_ptcls3D` outside historical notes.
- `git diff --check` passes.

Outstanding for the maintainer: rebuild; run `simple_test_ui_visibility`;
inspect `simple_exec prg=list` for the two new headings between Heterogeneity
Analysis and Denoising; regenerate the complete UI JSON and diff it against
the previous build (expected: category fields for ten programs, `summary` and
`display_name` for `automask`, one dropped descriptor); confirm
`prg=export_manifoldem_starproject` reports an unknown program. Record the
outcome here.
