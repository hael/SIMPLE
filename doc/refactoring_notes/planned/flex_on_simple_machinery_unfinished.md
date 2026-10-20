# FLEX on SIMPLE machinery: open work

Status (2026-10-20): the refactoring itself is implemented (phases 0 to 6). What is left is
verification of the last two phases and seven follow-ups. Background, rulings, contracts and the
day-by-day record are in
[`../completed/flex_on_simple_machinery_completed.md`](../completed/flex_on_simple_machinery_completed.md);
the explanatory report for the FLEX developer is
[`../completed/flex_on_simple_machinery_report.md`](../completed/flex_on_simple_machinery_report.md).

FLEX is the `flex_pca` heterogeneity workflow. "Phase 5" is the reuse work (FLEX now uses shared
SIMPLE code for sigma2, particle preparation, insertion and gather, the PCG lattice, clustering,
statistics). "Phase 6" is the document and comment updates.

The learning so far: Phases 0 to 4 ran on the test machine with measured exits. Phases 5 and 6 were
written by hand, one diff per item, and the agent could not build them. The maintainer built them
first on 2026-10-08 and fixed what the compiler and the fast gate found. The tier runs and the
before/after FLEX comparison on those diffs were not recorded as done.

## 1. Verification still to do

Run in this order; each step guards the next. Tolerances are the ones of section 8.2 of the
completed note.

1. **Fast gate on a Debug build** (`scripts/run_fast_gate.sh build`, or `ctest -L fast`). Covers the
   clustering testers, the weight-set withdrawal, the weighted insertion and FLEX's testers.
2. **Hard paths byte for byte.** Phase 5 touched two production paths that must not change:
   particle preparation for gridding reconstruction (`prep_imgs4rec`) and the PCG lattice layer
   (`pcg_lattice`). Run `simple_test_exec test=pcg_recon`, and the four reference reconstructions
   of the report (gridding and PCG, each shared memory and `nparts=2`) from the beta-galactosidase
   project; every `.mrc` must be identical to the kept references.
3. **Library tier:** `test=lib_reconstruction` (the `fractional_reconstruction` sub-suite, with 0/1
   weights byte-identical to hard labels).
4. **FLEX application phantom:** `test=lib_heterogeneity`. Per arm: recoded label accuracy at
   least the Phase 0 value minus 0.02 and at least 0.95; matched-map correlation within 0.02 of
   the earlier value; each matched map prefers its own truth. State numbers may permute (the
   clustering diffs relabel clusters by decreasing population).
5. **FLEX workflow gate:** `test=flex_pca_blobs nthr=4`. Its distributed mode exercises the worker
   commander, the plane cache and the embedding artifact.
6. **Refinement entries** that passed at Phase 0: `simulated_workflow suite=1jxy` and `suite=6vxx`,
   `solve3D_addon nthr=8`, `cont_refine3D_1jxy nthr=8`.

For any Phase 5 diff that does not touch clustering, the FLEX A/B demands identical recoded hard
labels and an embedding identical bit for bit (or within 1e-4 relative after matching component
order and sign where the plan says so). For a clustering diff the FLEX margins apply instead.
Record the outcome per diff; this note moves to `completed/` when steps 1 to 6 pass.

## 2. Open points

Each is a follow-up, none blocks the verification above.

1. **Polar bank and symmetry (a correctness fix to make).** The E-step's polar bank snaps each
   particle to the bank direction nearest its raw plane normal, without first mapping the
   orientation into the asymmetric unit. SIMPLE-refined particles already lie in the bank's
   coverage (on a d2 project all 5,513 particles snapped within 3.05 degrees). Imported or
   symmetry-expanded orientations do not: over random d2 orientations the mean snap error was 31.7
   degrees without the mapping and 1.25 degrees with it. Only polar rings above the hybrid radius
   use the snapped direction. Fix: apply the symmetry mapping first. Files:
   `src/main/flex/fit/simple_flex_pca_polar.f90`, `simple_flex_probe_fit_estep.f90`. Acceptance:
   a test over random orientations of a symmetric group shows the snap error near 1.25 degrees,
   and the d2 application phantom is unchanged.
2. **UMAP default.** The maintainer ruled UMAP off by default; the code still defaults to
   `umap=yes` (`src/main/params/simple_parameters.f90`, input in
   `src/main/ui/simple/simple_ui_heterogeneity.f90`). Change the default and its UI text together.
3. **Reader for the registered embedding.** The embedding artifact is registered in the project
   (`out` segment) but nothing reads the entry. A states-only rerun still takes `infile=`. Add the
   reader, or drop the registration; the maintainer has not ruled which.
4. **Weights between 0 and 0.01 on PCG.** Ruled and implemented: `PCG_WEIGHT_THRESHOLD` stays 0.01,
   so in a PCG per-state refinement rows with a weight in (0, 0.01) are aligned by
   `refine3D_auto state=X m_estimator=flex` but not reconstructed. Nothing to change; keep the
   behaviour in mind when reading per-state results, and document it where users look if it
   surprises anyone.
5. **Compile chain of the clustering modules.** `simple_kmeans`, `simple_kcenter` and
   `simple_xd_gmm` import `simple_clustering_utils` only for the seeding, which pulls the
   distance-matrix algorithms into their compile chain. Split the seeding into a small module so
   the three modules stop importing the heavy one (see `doc/policies/compile_time_policy.md`).
6. **Emptying a project segment.** Writing a segment inside a project file does nothing for an
   empty table. Emptying on disk is now the explicit `sp_project%clear_segment_inside`
   (`binoris%empty_segment_inside`); the weight-set withdrawal uses it. Other code that removes the
   last entry of a segment must call it too. Existing callers were not audited: audit them.
7. **Tracked Fortran indexes and one skill.** The tracked indexes under
   `doc/code_overview/fortran-indexes/` still list deleted modules until
   `scripts/gen_fortran_indexes.pl` is run. The skill `simple-main-nu-filt` still lists `flex_pca`
   as a caller of the nonuniform filter; it no longer is. Skill edits need the maintainer's
   approval (a proposed update is with him).

## 3. Rulings that bind this work

- Sampling, the refinement and FLEX keep the contracts of section 3 of the completed note; the
  goal, tolerances and exits change only by a maintainer ruling.
- `PCG_WEIGHT_THRESHOLD` stays 0.01 (ruled 2026-10-08).
- `refine3D_states` refuses `m_estimator=flex`; `refine3D` refuses it with more than one state.
- The plane cache is deleted by the run that built it once the master has its embedding.
- No new `goto`, no `//&` continuation in `THROW_HARD`/`THROW_WARN` arguments, no
  `flag .and. array(i)` on one line where the array may be unallocated.
- Checks to keep green: `scripts/check_test_registry.py`, `scripts/check_descr.py`,
  `scripts/check_flex_dag.py`; no file-mode changes.
