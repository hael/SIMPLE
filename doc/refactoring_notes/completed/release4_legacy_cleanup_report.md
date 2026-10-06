# Release 4 legacy cleanup: refactoring report

Date: 2026-10-03.
Author of the changes: Claude (agent), working on the instructions of Hans
Elmlund, the SIMPLE maintainer, who made every decision recorded here.

## 1. Summary

SIMPLE release 4 is a deliberate clean break with earlier releases. Before it,
the code base carried a large amount of support for the past:

- old program and parameter names, kept as aliases;
- readers for old file layouts;
- automatic repairs for files written by earlier versions;
- diagnostics and options left over from finished migrations;
- code that nothing could reach any more.

This refactor removed that baggage.

In numbers: about 210 files changed, 17 files were deleted, and roughly 5,800
lines were removed against 1,800 added. The regenerated, machine-written code
indexes under `doc/code_overview/fortran-indexes/` are not counted. About 3,800
of the removed lines are Fortran source.

The change is complete in the working tree but not yet committed. It sits on
top of commit `5b8ad9f06`, which renamed the de novo reconstruction programs
`abinitio2D`/`abinitio3D` to `solve2D`/`solve3D` and the 2D refinement program
`cluster2D` to `refine2D`. The maintainer has compiled the change and run the
fast test gate (the set of quick unit tests run on every build) twice. Each
failure found was fixed, as described in section 9. A final build and test run
on the maintainer's machine is the remaining step before commit.

One related change was planned but deliberately not made: how trailing
reconstruction starts. Section 12 explains what it is and where the plan is.

## 2. Terms used in this report

**Project file.** A SIMPLE project is stored in one binary `.simple` file. It
holds segments of records:

- micrographs;
- particle stacks, the image files that hold the extracted particle images;
- particles, in two segments: `ptcl2D` for 2D analysis and `ptcl3D` for 3D
  analysis;
- class averages and output volumes.

Each particle is a fixed-width record of floating-point values, one value per
named parameter (orientation, shift, defocus, state and so on).

**Stack indices.** Each particle row points at its image in two steps.
`stkind` names the stack row (the image file). `indstk` names the image's
position inside that file. A stack row describes:

- the project rows that belong to it, as the range `fromp` to `top` with the
  count `nptcls`;
- the number of images physically in the file, `nptcls_stk`.

After particles are removed from a project without rewriting the image files
("pruning"), the project rows and the physical images no longer line up. Only
`indstk` then says which image a particle is.

**sigma2.** The per-particle noise power spectra used by the Euclidean
(likelihood-based) alignment objective. They are stored in a versioned binary
state file.

**Trailing reconstruction.** An option of 3D refinement (`trail_rec=yes`).
When each iteration updates only a fraction of the particles, the new partial
reconstruction is blended with what came before instead of replacing it. The
"before" is kept as an accumulator chain: running sums of Fourier data and
sampling weights, stored on disk between iterations.

**Reconstruction backends.** `gridding` is the default 3D reconstruction.
`pcg` is an iterative alternative that solves the reconstruction as a
least-squares problem with the preconditioned conjugate-gradient method.

**NICE.** SIMPLE's web interface, a Django application under `nice/`. It
stores the jobs users launch in a database.

**Program descriptor.** Every SIMPLE program is described once in the UI layer
(`src/main/ui/`). The descriptor holds its name, summary, help, inputs and the
title shown in NICE. The command-line parser and NICE are both driven by it.

## 3. Principles decided for this refactor

The maintainer set these rules before and during the work:

1. **No backwards compatibility.** The maintainer: "No backwards
   compatibility. This is the current standard."
   - There are no aliases for old program or parameter names, no readers for
     old file formats, and no way to resume runs started by an earlier
     release.
   - When a file format changes, its version number goes up and older files
     are refused with a clear message.
2. **One exception: project files.** Projects are users' data and must keep
   working across releases. Project files whose particle records are narrower
   than today's (written before newer parameters existed) are still read, with
   zeros in the missing values.

   The exception covers the record width only, not old bugs. Projects that
   lack the physical stack indices (section 7) are not supported at run time,
   because silently guessing the indices was a bug; a separate repair program
   converts them where the answer can be proven.
3. **Areas left alone.** Another developer is working on the graphical
   interface and the streaming (on-the-fly) pipeline, so those parts were not
   changed except where unavoidable (section 11):
   - `src/main/stream/`;
   - `production/simple_stream.f90`;
   - `src/utils/gui/`;
   - `nice/`.

   GPU offloading is work in progress and was also left alone. MPI support
   stays.
4. **Kept on purpose:**
   - several search modes that are reachable only through undocumented
     parameter values;
   - unused GPU benchmark kernels and the GPU reconstruction routine;
   - all compile scripts.

## 4. What users will notice

### 4.1 Programs

- **Old program names no longer run.** `abinitio2D`, `abinitio3D` and their
  variants, and `cluster2D` and its variants, now fail with "not recognized".
  Scripts must use `solve2D`, `solve3D` and `refine2D`. The NICE database
  migration added with the rename still converts stored jobs.
- **`validate_projfile` is now `fix_projfile`.** The program changed purpose,
  from reporting problems to converting old projects (section 7), and was
  renamed without an alias. A new NICE migration (`0007`) renames stored jobs.
- **Every program has an explicit title in NICE.** Previously, 131 programs
  had no title of their own and NICE showed their one-line summary instead.
  Each now has a short title, reviewed and signed off by the maintainer during
  the release-4 cleanup. The six
  streaming programs keep their previous wording, because the streaming
  interface fixes it (`production/stream_ui_contract.json`). A title is now a
  required part of every program descriptor.

### 4.2 Command-line parameters

These parameters were removed. A command line that still passes one is
rejected with "argument is not allowed".

- **`inpl_cont`.** It switched the continuous optimisation of the in-plane
  angle during discrete orientation search. That optimisation is now always
  on. The continuous Cartesian refinement mode (`refine=cont`) still turns it
  off internally, and the convergence report still shows how often it improved
  a pose.
- **`pcgop`.** It selected the operator of the PCG backend. Production only
  ever accepted the fast kernel operator; the slow reference operator remains
  in the tests.
- **`projrec`.** It switched on an experimental reconstruction that first
  summed 2D Fourier slices per projection direction. Neither the PCG backend
  nor continuous refinement supported it, and it was removed completely.
- **`euclid_diag`.** It printed diagnostics written for a finished change to
  the amplitude scaling of reference maps.
- **`nu_msk_beta`, `nu_msk_dens` and `nu_msk_rel`.** Tuning knobs of the mask
  built from nonuniform-filter evidence in `nu_filt3D`. The filtering policy
  allows only `nu_msk_sig` and `amsklp`; the others now always take the fixed
  production values.
- **`nchunks` on `particle_sieving`.** It was an alias of `nparts`.
- **32 parameters that were parsed but never read:**

  `nptcls_per_subcls`, `pdfile`, `system`, `detector`, `automatic`,
  `thres_low`, `thres_up`, `cn_type`, `column_sampling`, `corr_thres`,
  `cs_thres`, `dfsdev`, `fraczero`, `gauref_last_stage`, `icm_stage`,
  `linethres`, `lp_discrete`, `nang_nbrs`, `ncls_sub`, `ninit`,
  `nparts_per_part`, `nsearch`, `nspace_max`, `pcrot`, `phranlp`,
  `protocol`, `randomise`, `refine_type`, `shift_stage`, `smpd_pickrefs`,
  `use_thres`, `wiener_const`.

These choices were removed:

- **`refine=snhc` and `refine=shc_neigh` for `refine3D` and `nano3D`.** They
  were offered but always stopped with "unsupported".
- **`picker=old`.** It was offered but always stopped with "no longer
  supported". An unknown picker name now stops the run, where before it was
  silently ignored.
- **`sigma2_convert sigma_action=parts_import`.** It converted sigma2 files
  written before the current state-file format.
- **`extract_substk state<0`.** This meant "extract every particle regardless
  of state" and is now refused.

Error messages that said "no longer supported" about parameters that are
still valid in other programs now say "does not take" or "is not available".

### 4.3 Files

- **Project files: particle records.** Each particle carries 53 named values.
  On disk each record now has 64 slots, of which slots 54 to 64 are written as
  zeros. These spare slots let later releases add parameters without changing
  the record width again.
  - Slot 42, which held an unused value (`cc_nonpeak`), is now free.
  - In memory a particle still holds only the 53 named values, so the reserve
    costs disk space only: 44 bytes per particle in each particle segment.
  - Older, narrower records are read with zeros in the missing slots; wider
    records are refused.
- **Project files: stack indices.** Every stack row that has particles must
  record `nptcls_stk`, and every particle must have an `indstk` between 1 and
  `nptcls_stk`.

  Earlier releases filled in a missing value with a guess (the particle's
  position within its stack's project range). That guess is wrong as soon as a
  stack has been pruned, so it could read the wrong image. Any operation that
  needs a missing index now stops with an error:
  - reading a particle image;
  - pruning a project;
  - merging projects;
  - writing a RELION STAR file.

  `fix_projfile` (section 7) converts old projects where the indices can be
  proven.
- **sigma2 state files: format version 2.** An always-zero checksum field and
  8 reserved bytes per particle, both kept only "for on-disk compatibility",
  are gone. Version-1 files are refused.
- **Class-average quality model files: version 11.** These are the learned
  models used to reject poor 2D class averages.
  - Ten threshold fields that the classifier no longer read are gone, from the
    file format and from the three built-in models.
  - Files of any other version are refused, as are two retired field names.
  - The learning program's training tables must now declare their version and
    context, and must contain every feature column.
- **Trailing-reconstruction chain files.** Recognition of the version-1 index
  file is removed. An old chain was already discarded and rebuilt; it now
  simply counts as unreadable.
- **solve3D add-on manifests.** These are the run records of the
  `solve3D_addon` workflow. Fields belonging to removed features are no longer
  silently accepted.
- **Flexibility-analysis caches.** A run of the flexible-heterogeneity
  analysis can no longer be resumed from the old separate deconvolution cache
  file, which nothing has written for some time.
- **Reference volumes.** Maps written under an old amplitude convention used
  to be detected by a heuristic and silently rescaled; that could misfire on
  legitimate low-contrast maps. The rescaling is gone, and the warning about
  suspicious scaling stays.

### 4.4 Build and scripts

- **Version.** The project version is 4.0.0, so the compiled library is now
  named `SIMPLE4.0.0`.
- **librt.** The build no longer requires or links the POSIX real-time library
  librt. Nothing called the message-queue functions it provided, and current C
  libraries contain the clock function directly. This still needs to be
  confirmed on the Oracle Linux test machine.
- **`BUILD_DOCS`.** The CMake option did nothing usable and is gone.
- **Benchmark files.** Refinement now writes one benchmark file per parallel
  partition, `REFINE3D_BENCH_ITERnnn_PARTppp.txt`. It no longer also writes a
  copy for partition 1 without the suffix.
  - The two benchmark tools (`plot_refine3d_bench.py`, `parse_bench.pl`)
    represent each iteration by partition 1.
  - `parse_bench.pl` keeps the benchmark families apart (the new `benchmark`
    column) and has a `--self-test` option.
- **Scripts deleted:**
  - unused developer scripts;
  - personal research scripts with hard-coded paths;
  - three 2020-era user tools (`avg_sym_stats.py`, `relion2emanbox.pl`,
    `relion2simple.pl`); `relion2simplegui.pl` remains;
  - a single-file metrics parser that duplicated
    `parse_solve3D_metrics_all.pl`.
- **Other build fixes:**
  - the persistent worker prints the real git version instead of a hard-coded
    one;
  - the tcsh setup template uses valid syntax;
  - a placeholder sentence in the comments of the compile scripts was
    removed.

## 5. Bugs fixed along the way

- **Programs that swap a particle stack read the wrong images.**
  `reimport_particles` and `trajectory_swap_stack` replace a project's images
  with one new stack. They kept each particle's old in-stack index. After a
  prune those indices pointed past or into the wrong part of the new stack.
  Both now set every particle to stack 1 with the index of its row.
- **RELION STAR import stored no in-stack index.** Imported projects relied
  entirely on the range guess. Import now stores the index given in the STAR
  file's image name.
- **Splitting a stack.** Splitting into parts wrote in-stack indices only if
  the project already had them; it now always writes them.
- **Offered but unsupported choices** (`refine=snhc`, `refine=shc_neigh`,
  `picker=old`) are gone, and an unknown picker name is no longer silently
  ignored.

## 6. Internal simplifications (no change in behaviour)

- **Class-average quality.** An unreachable fallback classifier (k-medoids
  clustering with Otsu thresholding, about 600 lines) is removed, with its
  cache and with diagnostics that always reported zero.
- **Old command-line paths.**
  - The pre-UI command-line parser and its printing routine are removed.
  - The table of internal ("private") programs drops ten entries nothing
    could reach.
  - Its parameter dictionary drops 266 of its 289 keys.
- **Dead programs and tests.**
  - An edge-detection program that no router dispatched is removed.
  - So is an unbound test of the streaming mini-workflow.
  - Unreachable aliases and two unreachable input checks are removed.
- **solve3D.** A disabled "gold-standard stage" and the branches behind it are
  removed. Gold-standard refinement belongs to `refine3D_auto`.
- **Projection-direction reconstruction** (`projrec`, section 4.2) is removed
  from the matcher, the reconstructor, the class-average accumulator code, the
  manifest and the memory benchmark.
- **Euclidean diagnostics.**
  - The arrays and report behind `euclid_diag` are removed, and so is a test
    column for an old de-apodisation variant.
  - The Cartesian calculation of the sigma2 contribution now computes
    reference power, particle power and normalised loss only when a test asks
    for them. Production asks only for the residual, as the polar calculation
    already did.
- **Distributed execution.**
  - The routine that marks a partition as finished had two names; one is
    gone.
  - A command-line flag passed between strategy and matcher
    (`force_volassemble`) is removed, with the routine that read it. That
    routine renamed output volumes only inside distributed workers, where
    nothing used the names.
- **Housekeeping.**
  - The deletion of file names from older runs (`polar_refs*.bin`, unpadded
    alignment documents) is removed.
  - Stripping an obsolete path from old projects' computing environment is
    removed; the runtime value is always used.
  - The Windows build loses eight stubs for removed message-queue functions.
  - `.gitignore` is cleaned, and a tracked CMake cache file is removed from
    git.
- **Wording.** Comments, messages and names that called the current default
  path "legacy" are corrected. The optimiser constructor `new_legacy` is now
  `new_alternating`, because it alternates discrete in-plane and shift
  optimisation.

## 7. Converting old projects: `fix_projfile`

`fix_projfile projfile=<name>.simple` reads a project from an earlier release
and writes `<name>_fixed.simple`. It writes the result only when every problem
was either repaired or never existed; otherwise it lists the problems and
writes nothing.

It proceeds in this order:

1. **Stack ranges.** Each stack row's range of project rows is rebuilt from
   the particles: the 2D segment, or the 3D segment when there is no 2D one.
   A range pointing outside the project is ignored, with a warning, and
   rebuilt from the stack's particle count. Ranges whose total cannot match
   the number of particles are errors.
2. **Stack assignment.** In both particle segments, every particle is given
   the stack whose range contains its row. A stored stack number that names a
   different stack is replaced and reported. A particle that no stack's range
   contains is an error.
3. **Physical image counts.** A stack row without `nptcls_stk` gets the image
   count from its file header. Such a stack predates the fix, so its stored
   in-stack indices are not trusted:
   - when the file holds exactly one image per project row, every particle
     gets index = its position within the stack's range, checked to lie
     inside the file;
   - otherwise the images cannot be matched to particles, and the stack is
     reported (the particles must be re-extracted or re-imported).
4. **Existing indices.** A stack that already records `nptcls_stk` keeps its
   particles' indices; a missing or out-of-range one is reported.
5. **Missing stacks.** A project with particles but no stack rows is an error.

Writing the result also brings the particle records to the current width.
`doc/policies/project_stack_indexing_policy.md` was rewritten to state these
rules and the required fields. It had still described the run-time guesses.

## 8. Independent review

An independent review of the change set found seven issues. All are resolved.

1. **The repair could write an inconsistent project.** If a 3D-segment
   particle named a valid but wrong stack, the derived index came from the
   wrong stack and was not checked. Fixed by step 2 of section 7 and by the
   range check in step 3, with a new test.
2. **The 64-slot record also cost memory.** At first the spare slots existed
   in memory as well, about 420 MiB per particle segment for ten million
   particles. Fixed by the split described in section 4.3: 53 values in
   memory, 64 slots on disk. The maintainer still has to confirm this split,
   because it changes his original decision of 64 slots everywhere.
3. **A project with particles but no stacks was "fixed".** It is now an error
   (step 5 of section 7), with a new test.
4. **The benchmark parser mixed partitions and families.** Fixed as described
   in section 4.4, with a self-test.
5. **Jobs stored in NICE kept the old program name.** Fixed by migration 0007
   and an updated NICE test fixture.
6. **Diagnostic arithmetic ran in production.** Fixed as described in
   section 6. A first version allocated the diagnostic sums only when they
   were requested, which drew a "may be used uninitialised" warning from
   gfortran. The few hundred values are now always allocated, and only the
   arithmetic is conditional.
7. **The checked-in code indexes described removed routines.** The maintainer
   regenerated them.

## 9. Problems found by the maintainer's builds

| Symptom | Cause | Fix |
| --- | --- | --- |
| Syntax error in `simple_binoris.f90` | An error message passed to the `THROW_HARD` macro was continued over two lines; the preprocessor cannot handle that. | The message is now built in a variable first. |
| `./compile_debug.sh`: permission denied | Editing the compile scripts in place removed their executable permission (12 files). | Permissions restored; file modes are now checked before each hand-over. |
| Test-registry check failed | A deleted test suite was still listed in the help of the reconstruction tests. | The listing and a stale command in a policy document were removed. |
| Project tests failed in `fix_projfile` | The tool counted a range it had just repaired as an error, and so refused to write. | A repaired range is now a warning. |
| Cartesian alignment tests failed with "nptcls_stk not present" | Four test projects were built without the now-required `nptcls_stk`. | The test projects were completed. |
| Six streaming programs had no title | The script that added titles recognised programs by executable name and missed `simple_stream`. | Titles added; every program definition is now checked by type. |

## 10. Tests

**New tests:**

- `fix_projfile` on:
  - a three-image old stack with a missing index, an out-of-range index and a
    wrong range;
  - a two-stack project in which a 3D particle names the wrong stack;
  - a project with particles but no stacks, which must produce no output.

  The program gained an optional error count for the last case, so a test can
  observe the refusal without the program stopping.
- Project records:
  - reading particle records of 53 and 40 values;
  - the layout, spare slots and both widths of the current record.
- The layout of version-2 sigma2 files.
- `perl scripts/parse_bench.pl --self-test`.

**Changed tests:**

- The project-merge test uses a valid second project.
- The solve3D manifest test checks that fields of removed features are
  refused.
- The in-plane test checks that no refinement program exposes `inpl_cont`.
- The computing-environment test checks that the runtime SIMPLE path wins.
- Four test projects record `nptcls_stk`.

**Removed tests:**

- the class-average accumulator tests, with `projrec`;
- the test of the old program-name table;
- the old `validate_projfile` test.

## 11. Changes in areas otherwise left alone

- **Streaming:**
  - one import line in `src/main/apis/simple_stream_api.f90`, for the renamed
    partition-completion routine;
  - explicit titles in the streaming program descriptors, with unchanged
    wording.
- **GPU:** the GPU reconstruction routine lost one call and one import when
  the `force_volassemble` flag was removed. The maintainer approved this.
- **NICE:** the new migration and one test fixture.

## 12. What remains

- **Before committing:**
  - the final build and test run;
  - the maintainer's confirmation of the in-memory/on-disk record split;
  - a build on the Oracle Linux test machine without librt;
  - when v4.0.0 is tagged, the version strings in `README.md` and
    `doc/installation.md`, which still name v3.0.0.
- **Trailing reconstruction: done (2026-10-04).** When `trail_rec=yes`
  starts without an accumulator chain, both backends used to blend finished
  half maps (restored, regularised volumes) with the previous ones for one
  iteration. The maintainer ruled that this must go: "No trailing should ever
  happen on finished halfmaps." That iteration now seeds the chain from the
  current sample at full mass and ships the current sample's map, with no
  extra reconstruction and without reading the previous maps. Validated on
  beta-galactosidase against a baseline of the old code: final resolutions
  are unchanged within run-to-run spread. Design, record and report:
  `doc/refactoring_notes/completed/trailing_reconstruction_without_halfmap_blend.md`
  and `trailing_reconstruction_without_halfmap_blend_report.md` beside it.
- **`solve2D_chunks`: kept and moved (2026-10-06).** The experimental program
  that runs `solve2D` on stack-bound subsets of a project stays (formerly item
  C15 of the inventory). It has left the stream: its commander,
  `commander_solve2D_chunks`, sits beside `commander_solve2D` in
  `simple_commanders_solve2D`, and its chunk type is `solve2D_chunk` in
  `src/main/solve/simple_solve2D_chunk.f90`. The program name, its UI and its
  behaviour are unchanged. Plan:
  `doc/refactoring_notes/planned/solve2D_chunks_move_plan_2026-10-06.md`.
- **`simple_nice` and `simple_guistats`: retired (2026-10-06)**, formerly items
  A1 and C13 of the inventory, about 2,300 lines.
  - `simple_nice_comm` sent messages without a `version`, which NICE's API
    answers with HTTP 400. Nothing it sent arrived, and no answer ever set
    its `exit` or `stop`. The nine commanders that created it (`refine2D`,
    `project_core`, `project_cls`, `project_ptcl`, `starproject`, `imgops` and
    `single`'s `nano2D`) no longer do. Batch commanders report through
    `simple_gui_communicator`.
  - `guistats` was left only as the pool's fallback for a missing
    `cavgs_iterNNN.jpg`, which refine2D writes after every iteration. Its
    single-column strip did not match the GUI's tile grid. The pool's own
    sprite sheet (`write_jpeg`) covers the case.
  - `GUISTATS_FILE` went with it.
- **Streaming and NICE clean-ups.** Dead routines, stale constants, unused
  NICE fields, a squash of the NICE migrations, old statistics channels,
  compatibility parameters and module renames. These wait for the developer
  working in that area, and are listed in
  `doc/refactoring_notes/planned/release4_legacy_cleanup_inventory.md`.
- **GPU.** An unreferenced debug compile script waits for the GPU work.
