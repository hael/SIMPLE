# SIMPLE `scripts/` inventory and cleanup candidates

Status: initial cleanup implemented

Inventory date: 2026-09-22

Scope: the retained contents of `scripts/` after the owner-approved removal
pass. Removed files remain recoverable from Git history.

## 1. Executive summary

The directory now contains **62 files**:

| Kind | Count |
|---|---:|
| Perl | 23 |
| Python | 27 |
| Shell | 5 |
| Shell setup templates | 2 |
| TSV and text policy/data files | 2 |
| CMake, Markdown, JSON, and directory metadata | 3 |

Most of the disorder is not caused by the number of files. It comes from five
different classes of artifact sharing one installable directory:

1. build-time generators and CI gates;
2. repository review and source-rewriting tools;
3. supported user-facing utilities;
4. scientific analysis and benchmarking scripts;
5. one-off, machine-specific experiments and historical migration scripts.

The owner decision is to continue installing the directory wholesale for now.
That is reasonable while the script set is being sorted out, but the install
rule should at least exclude `__pycache__/`, `*.pyc`, editor artifacts, and
other ignored/generated files. Otherwise the installation depends on what
happens to be present in the source checkout at install time.

The cleanup should begin by removing scripts whose obsolete external
dependencies are known, then tidy the retained flat script set without
changing the current broad installation policy.

## 2. Status vocabulary

| Status | Meaning |
|---|---|
| **Core** | Referenced by the build, compile wrappers, or an active repository contract. Keep and protect. |
| **Supported** | Plausible user-facing utility. Keep and document dependencies. |
| **Internal** | Maintainer/review utility. Keep in the repository; it remains installed while the directory-wide policy is active. |
| **Consolidate** | Useful behavior overlaps another script or belongs in a clearer tool family. |
| **Review** | Purpose is plausible, but support, ownership, portability, or callers are unclear. |
| **Archive candidate** | Strong evidence of a personal, obsolete, broken, or superseded artifact. All candidates identified in this inventory were approved and removed in the initial cleanup. |

Repository reference counts are only evidence, not proof of use. Many scripts
are intentionally invoked by people and will have no in-tree caller.

## 3. Immediate findings

### 3.1 Installation boundary

The owner has chosen to retain the broad `install(DIRECTORY ...)` rule for now.
Harden it with exclusions for generated and local artifacts. The eventual
conceptual groups remain:

- installed user utilities;
- build-only generators;
- maintainer tools;
- analysis/benchmark tools;
- examples or archived scripts.

Installing the complete tracked folder is acceptable. An install should still
be reproducible from tracked content and must not depend on ignored or
untracked files in a developer checkout.

### 3.2 Static health

Lightweight parsing, without running any script, found:

- all 23 retained Perl scripts pass `perl -c`;
- all tracked shell scripts pass `bash -n`;
- 26 of 27 retained Python files parse as Python 3;
- `avg_sym_stats.py` is the remaining Python file that does not parse as
  Python 3;
- ignored `scripts/__pycache__/` and `scripts/ui/__pycache__/` directories were
  present locally; `scripts/ui/` contained bytecode but no tracked source.

The parsing result establishes syntax compatibility only. It does not establish
that optional packages, external programs, fixtures, or scientific behavior
are available or correct.

### 3.3 Clear duplication or succession

- `relion2simple.pl` and `relion2simplegui.pl` share almost all parsing logic;
  the latter mainly adds `_rlnImageName` handling.
- `nanokpcadn.py` and `nanokpcadn2.py` are two hard-coded variants of the same
  experiment. Nanoparticle tooling is part of SIMPLE and should be retained,
  but these variants should become one parameterized tool.
- the retained nanoparticle scripts form a coherent family. They will remain
  in the flat top-level script directory rather than being moved into a new
  hierarchy.

### 3.4 Remaining cleanup signals

- `delete_unused_variables.pl` was removed; its dry-run-first replacement is
  `fix_warnings.py`.
- `nanokpcadn.py`, `nanokpcadn2.py`, `figs_core_nptcls.py`, and
  `figs_nptcls.py`: retained nanoparticle tools that need parameters in place
  of hard-coded workstation or cluster datasets;
- `tab2space.pl` is a minimally validated script that mutates every local
  Fortran file and remains under review;
- `avg_sym_stats.py` needs a Python 3 decision;
- `relion2simple.pl` and `relion2simplegui.pl` still need a consolidation or
  retirement decision.

## 4. Full inventory

### 4.1 Build, registry generation, test gates, and packaging

| File | Purpose | Status | Recommendation |
|---|---|---|---|
| `CMakeLists.txt` | Installs the scripts tree. | **Core** | Retain directory-wide installation for now, but exclude caches and untracked/generated artifacts. |
| `default_audit.py` | Generates/audits command defaults used by CMake. | **Core** | Keep at a stable path; do not install. |
| `insert_git_commit_hash.pl` | Updates generated Git-version content during the build. | **Core** | Keep at a stable path; document its environment and generated outputs. |
| `simple_args_generator.pl` | Generates the argument list/module from typed parameters. | **Core** | Keep at a stable path; do not merge casually with maintenance utilities. |
| `ctest_budget.py` | Enforces fast-test total and per-test budgets. | **Core** | Keep with test infrastructure; standard-library-only and documented. |
| `run_fast_gate.sh` | Runs the provisional/fast CTest gate after test builds. | **Core** | Keep with test infrastructure; referenced by compile wrappers and CMake. |
| `test_args.tsv` | Per-test arguments for timing characterization. | **Core** | Keep beside the timing runner; treat as test policy data. |
| `test_review_dossier.py` | Generates canonical test dossiers and inventory. | **Core** | Keep with test infrastructure; output remains generated review evidence. |
| `test_timing_run.sh` | Runs test routes in bounded isolated directories. | **Core** | Keep with test infrastructure; it is not a general installed command. |

### 4.2 Repository audit and navigation tools

| File | Purpose | Status | Recommendation |
|---|---|---|---|
| `capture_cli_review.py` | Captures UI registry JSON and readable CLI snapshots. | **Internal** | Keep in the flat directory as a review tool. |
| `check_add_ui.py` | Finds mismatched UI registration names. | **Internal** | Keep; add a short usage contract or fold into a unified UI audit command. |
| `check_descr.py` | Checks description headers under source trees. | **Internal** | Keep if the description-header policy is active; otherwise retire with that policy. |
| `gen_fortran_indexes.pl` | Generates Fortran navigation indexes. | **Internal** | Keep at a documented stable path. |
| `generate_codeoverview.pl` | Generates the repository code overview. | **Internal** | Keep at a documented stable path. |
| `list_procedures.pl` | Lists/parses Fortran procedures. | **Internal** | Retain as a maintainer tool; document its parser limitations. |
| `placeholder_audit.py` | Generates UI placeholder review evidence. | **Internal** | Keep with the UI audit family. |
| `scripts.inf` | Directory description consumed by code-overview generation. | **Internal** | Keep at the directory root while the generator requires this convention. |

### 4.3 Source maintenance and mechanical refactoring

| File | Purpose | Status | Recommendation |
|---|---|---|---|
| `clean_simple_uses.pl` | Rewrites redundant SIMPLE umbrella-module imports. | **Internal** | Keep only as an explicit source migration tool; dry-run support would reduce risk. |
| `fix_warnings.py` | Dry-run-first cleanup of mechanical gfortran warnings. | **Internal** | Keep as the supported warning-fix tool. |
| `move_constructors.pl` | Mechanical constructor relocation helper. | **Review** | Add a shebang/usage contract or archive after its refactor is complete. |
| `program2mod_sub.pl` | Converts Fortran programs into callable subroutines. | **Internal** | Keep during the test refactor. |
| `remove_git_version_conflict.pl` | Repairs generated version-call conflicts in production sources. | **Internal** | Keep only if this conflict remains a recurring workflow; otherwise archive. |
| `split_commanders.sh` | Historical mechanical commander splitter. | **Review** | Treat as a migration tool, not a supported general parser. |
| `tab2space.pl` | Rewrites tabs in all local `*.f90` files. | **Review** | Replace with formatter/editor policy or add dry-run/path arguments. |

### 4.4 Build and performance profiling

| File | Purpose | Status | Recommendation |
|---|---|---|---|
| `profile_build.sh` | Profiles clean, touched, and interface-changing builds. | **Internal** | Keep adjacent to its launcher in the flat directory. |
| `profile_build_launcher.sh` | CMake compiler launcher recording per-source timings. | **Internal** | Keep adjacent to `profile_build.sh`; not independently installed. |
| `calc_avg_exec_time.pl` | Summarizes elapsed times and iterations from output logs. | **Review** | Group with benchmark parsers and correct the stale usage name in its header. |
| `parse_bench.pl` | Converts SIMPLE benchmark report sections into tables. | **Internal** | Keep with benchmark/reporting tools. |
| `plot_refine3d_bench.py` | Produces CSV/SVG summaries from refine3D benchmark reports. | **Internal** | Keep with benchmark/reporting tools. |

### 4.5 Memory-estimator subsystem

This is already the best-organized script family. Its README distinguishes the
installed estimator from calibration and reporting tools.

| File | Purpose | Status | Recommendation |
|---|---|---|---|
| `memory_estimator.py` | User-facing resource estimator. | **Supported** | Keep as an intentional installed utility with its model JSON. |
| `memory_estimator_models.json` | Active fitted estimator coefficients. | **Supported** | Install beside the estimator; validate schema/version compatibility. |
| `memory/README.md` | Maintenance contract for benchmarks and model fitting. | **Internal** | Keep. |
| `memory/benchmark_motion_correct.py` | Generates/benchmarks motion-correction cases. | **Internal** | Keep under `memory/`; never install by directory accident. |
| `memory/benchmark_abinitio2d.py` | Generates/benchmarks abinitio2D cases. | **Internal** | Keep under `memory/`. |
| `memory/benchmark_abinitio3d.py` | Generates/benchmarks abinitio3D cases. | **Internal** | Keep under `memory/`. |
| `memory/fit_models.py` | Fits and optionally updates estimator coefficients. | **Internal** | Keep under `memory/`; preserve explicit write opt-in. |
| `memory/report_abinitio3d.py` | Builds abinitio3D memory reports. | **Internal** | Keep under `memory/`; document its substantial Python dependencies. |

### 4.6 User-facing data preparation and inspection

| File | Purpose | Status | Recommendation |
|---|---|---|---|
| `add2.bashrc.template` | Generates the installed Bash environment snippet. | **Core** | Retain: CMake, the installation guide, onboarding notes, and the post-install message rely on the generated `add2.bashrc`. |
| `add2.tcshrc.template` | Generates the installed tcsh environment snippet. | **Core** | Retain with the Bash template as part of the current post-install environment contract. |
| `check_dims_movs.pl` | Checks MRC/EER/TIFF movie dimensions. | **Supported** | Keep; document supported formats and external assumptions. |
| `defocus_stats_by_cls2d.py` | Reports class-wise defocus statistics from a project. | **Supported** | Keep; it is documented and uses SIMPLE executables explicitly. |
| `display_nu_locres.py` | Opens a density/local-resolution pair in ChimeraX. | **Supported** | Keep; clearly declare ChimeraX and Python-mode requirements. |
| `filetab_movs.pl` | Generates movie file tables. | **Supported** | Keep; already referenced by the user guide. |
| `filetab_mrc.pl` | Generates MRC/MRCS file tables. | **Supported** | Keep; already referenced by the user guide. |
| `generate_jpg_report.pl` | Produces an HTML image report for test/artifact directories. | **Supported** | Keep in the flat directory and clarify whether testing is its sole consumer. |
| `manual_particle_picker.py` | Interactive MRC particle picker. | **Supported** | Keep if this is a supported utility; document `mrcfile`, NumPy, and Matplotlib. |
| `maxdiam.pl` | Derives a maximum diameter from image/volume information. | **Review** | Clarify input contract and whether a SIMPLE command supersedes it. |
| `select_mrc_classes.py` | Interactive class-average selector. | **Supported** | Keep; document optional viewer/conversion dependencies. |
| `simple_term_stream.pl` | Creates the stream termination sentinel. | **Supported** | Keep if this remains the documented stream stop interface; a shell-independent implementation may be clearer. |
| `update_cls2d_os_out_paths.py` | Updates class-average output paths in projects. | **Supported** | Keep, but document backup/dry-run behavior because it mutates projects. |

### 4.7 Format conversion and interoperability

| File | Purpose | Status | Recommendation |
|---|---|---|---|
| `chimera_fitmap.py` | Runs a fixed Chimera fit/resample workflow on fixed filenames. | **Supported** | Retain because Chimera tooling is active; parameterize filenames without breaking the Chimera interpreter contract. |
| `relion2emanbox.pl` | Converts RELION coordinate STAR data to EMAN box files. | **Supported** | Retain as executable documentation of the coordinate/box convention, even if the conversion workflow itself is rarely used. |
| `relion2simple.pl` | Converts RELION STAR fields to SIMPLE key/value records. | **Consolidate** | Share one parser with `relion2simplegui.pl` or retire both in favor of current import commands. |
| `relion2simplegui.pl` | RELION-to-SIMPLE conversion including image/frame paths. | **Consolidate** | Merge its additional fields into one maintained converter with explicit modes. |

### 4.8 Generic filesystem helpers

No generic filesystem helper remains in this category after the approved
cleanup pass.

### 4.9 Scientific analysis and plotting

| File | Purpose | Status | Recommendation |
|---|---|---|---|
| `avg_sym_stats.py` | Aggregates symmetry-ranking statistics. | **Review** | Does not parse as Python 3; modernize and document input schema or archive. |
| `parse_abinitio_metrics.pl` | Extracts metrics from an abinitio run. | **Internal** | Keep with abinitio analysis tools. |
| `parse_abinitio_metrics_all.pl` | Aggregates metrics across abinitio runs. | **Internal** | Keep with its smaller parser; remove personal absolute paths from usage examples. |
| `plot_fsc_area_score.py` | Plots FSC-area score data under the active policy. | **Internal** | Keep; it has an explicit policy/documentation owner. |

### 4.10 Site-specific research and experiment scripts

| File | Purpose | Status | Recommendation |
|---|---|---|---|
| `figs_core_nptcls.py` | Generates nanoparticle-core figures and Chimera movies. | **Supported** | Retain as active Chimera/nanoparticle tooling; replace hard-coded paths and broad destructive globs with arguments and bounded work directories. |
| `figs_nptcls.py` | Generates nanoparticle CN/DOI/ANI/RADS figures. | **Supported** | Retain with the related Chimera/nanoparticle family and parameterize its paths. |
| `nanokpcadn.py` | Hard-coded kernel-PCA denoising experiment. | **Consolidate** | Retain the capability and merge it with its variant into one parameterized Python tool. |
| `nanokpcadn2.py` | Second hard-coded kernel-PCA denoising experiment. | **Consolidate** | Merge into the maintained nanoparticle denoising tool. |

## 5. Flat-directory policy

Do not reorganize the scripts into new category folders. Retained top-level
scripts remain directly under `scripts/`, preserving short paths and the
current installation behavior. The existing cohesive `memory/` subtree may
remain, but this cleanup will not introduce `build/`, `test/`, `audit/`,
`converters/`, `research/`, or similar subdirectories.

Purpose and support status will be expressed through this inventory, script
headers, and naming rather than directory hierarchy.

### 5.1 Language direction

Do **not** convert the directory wholesale to Python. Use this policy instead:

- New maintained data-processing, parsing, reporting, and file-manipulation
  tools should default to Python 3.
- Convert retained Perl utilities when they need substantive maintenance, when
  two tools are being consolidated, or when conversion removes a real external
  dependency. Preserve inputs and outputs with fixtures before rewriting them.
- Do not rewrite stable build generators merely for language uniformity.
  `simple_args_generator.pl` and `insert_git_commit_hash.pl` are short,
  build-critical, and already integrated with CMake. A rewrite is worthwhile
  only if SIMPLE deliberately removes Perl as a build dependency.
- Keep shell for thin process orchestration such as `run_fast_gate.sh`,
  `test_timing_run.sh`, and the compiler launcher. Python would add ceremony
  without improving their ownership or behavior.
- Respect embedded interpreter constraints. Active classic Chimera scripts may
  be tied to Chimera's Python version and API; do not apply an automatic Python
  3 conversion until the supported Chimera/ChimeraX runtime is named.
- Remove obsolete FREALIGN, EMX, and EMAN2/SPARX paths instead of porting them.

The practical target is therefore **Python 3 for maintained application-like
scripts, shell for orchestration, and unchanged Perl for stable build plumbing**.
This yields convergence without turning cleanup into a large behavior-preservation
project.

## 6. Implemented removal pass

The owner approved every archive candidate and confirmed that the hook samples
are unused. The following 23 files were removed:

```text
align_from_chimera.py
convert_frealign2simple.pl
convert_simple2relion.py
cp_early_stream_UT_results.pl
cp_filetab.pl
delete_unused_variables.pl
dm42mrc.pl
emx2simple.pl
kill_descendents.sh
mv_dir_content.sh
parse_NP_seg_inf.pl
parse_cn_dependent_stats.pl
parse_det_atms.pl
parse_nanoparticle_stats.pl
pre-commit.sample
pre-push.sample
puts_gitignore.pl
relion2simplectf.pl
rename.pl
s2h_min.pl
simple_args_varlist.pl
simple_star_test.sh
split_filetab.pl
```

The shell setup templates were initially included in the removal pass, then
restored after verifying that CMake generates and installs `add2.bashrc` and
`add2.tcshrc` and that the post-install documentation relies on them.
`fix_warnings.py` no longer names its removed predecessor.

## 7. Remaining work and decisions

- Continue installing the complete retained scripts folder for now, while
  excluding caches, bytecode, editor files, and other generated/local
  artifacts from installation.
- Retain the Bash and tcsh setup templates until the documented post-install
  environment mechanism is deliberately replaced across CMake, the README,
  installation documentation, onboarding material, and post-install message.
- Keep the directory flat; the existing `memory/` subtree remains the only
  deliberate functional subtree.
- Retain nanoparticle and Chimera tooling. Document the exact Chimera/ChimeraX
  interpreter contract and parameterize hard-coded nanoparticle paths.
- Retain `relion2emanbox.pl` as documentation of the coordinate/box convention.
  Decide whether to consolidate or retire `relion2simple.pl` and
  `relion2simplegui.pl`.
- Decide whether to port or retire the remaining Python-2-era
  `avg_sym_stats.py`.
- Add cheap validation for Python 3 parsing, `perl -c`, `bash -n`, shebangs,
  and the installed file list.

Future removals still require an ownership decision; absence of an in-tree
caller alone is not sufficient because scripts can be human-invoked
interfaces.
