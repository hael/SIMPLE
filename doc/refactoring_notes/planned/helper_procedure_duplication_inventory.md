# Helper procedure duplication inventory

Reviewed 2026-11-14 on `master`, HEAD `85213734e9eac15a23ef634e4a8bf6cfe3712652`, **including the current uncommitted changes**. This complements [the constant inventory](module_constant_inventory.md). It inventories source and proposes reuse; it does not change Fortran or infer who wrote a procedure.

## Findings and scope

The inventory identifies **23 families worth addressing**: nine with existing library operations to reuse, and fourteen needing a shared implementation or an API extension. The best first changes are ascending sort indices, POSIX quoting, median/standard-deviation helpers, STAR path trimming, and image-array destruction. Some replacements need small adapters; several are deliberately not classified as drop-in replacements.

Unlike constants, these procedures belong in their numerical, text, I/O, image, or workflow owners, not in definitions modules. Keep a contained procedure when it expresses local policy or adapts a shared API. The problem is duplicated implementation, not `contains` itself.

| Lexical census | Count |
| --- | ---: |
| SIMPLE-owned `.f90`/`.f08` files scanned under `src`, excluding `extlibs` | 872 |
| Procedure implementation blocks indexed, including tests | 11,766 |
| Production procedure blocks | 8,695 |
| Production internal procedures nested inside another procedure | 898 |
| Production commander procedures, all scopes | 491 |
| Of those, internal procedures | 225 |

These are discovery counts, **not duplicate counts**. The census excludes interface declarations and separates `_tester.f90`, `commanders/test`, and the test-exec sources from production. It recognizes explicit subroutine/function implementations, including `module subroutine/function`; it is not a compiler call graph. Candidate discovery used names, normalized bodies, renamed-body comparisons, and searches of existing utility APIs, followed by reading the candidate bodies and their proposed replacements. This is a reviewed inventory of concrete candidates, not a proof that every duplication has been found. Vendor code and independent test oracles are outside the cleanup scope.

Source links below refer to this worktree snapshot. The SHA-256 over sorted owned Fortran paths and contents, separated by NUL bytes, is `32456f4c583c6918c49ed711597b8b7beaf9a70988c46bc211ff30dc8ed41a19`.

## A. Reuse existing library operations

### A01. Ascending sort permutation — direct replacement

- **Local implementation:** contained `sort_order` in [single_commanders_trajectory_solve3D.f90:628](../../../src/main/commanders/single/single_commanders_trajectory_solve3D.f90#L628).
- **Existing API:** [simple_srch_sort_loc::dists2order:767](../../../src/utils/math/simple_srch_sort_loc.f90#L767).
- **Evidence/action:** both copy the real vector, construct indices `1:n`, and call the same `hpsort(copy, indices)`. Replace the local helper with `dists2order`. The caller at line 599 needs ascending order; `scores2order` reverses the result and is the wrong replacement. Both implementations preserve the input and inherit the same heap-sort tie behavior.

### A02. POSIX shell quoting — two implementations

- **Local implementations:** [simple_refine_motion_model_strategy::shell_quote:1304](../../../src/main/strategies/parallelization/simple_refine_motion_model_strategy.f90#L1304) and [simple_qsys_ctrl::shell_quote:691](../../../src/utils/qsys/simple_qsys_ctrl.f90#L691).
- **Existing API:** [simple_string_utils::posix_quote:630](../../../src/utils/text/simple_string_utils.f90#L630).
- **Evidence/action:** all surround a trimmed argument with single quotes and escape embedded apostrophes. The queue version uses the same escape sequence as the shared routine; the motion version uses an equivalent shell sequence. Use `posix_quote`, retaining a `type(string)` conversion where the queue API needs it. Preserve the existing trailing-blank trimming. Check apostrophes, empty strings, spaces, dollar signs, and backticks; compare shell argument values rather than requiring the two equivalent encodings to have identical text. Gnuplot quoting is a different contract (see exclusions).

### A03. Median — three implementations, different adaptation requirements

| Local helper | Existing implementation and action |
| --- | --- |
| Contained [simple_atomic_species::median_of:152](../../../src/main/nano/simple_atomic_species.f90#L152) | Use `simple_stat::median` for nonempty real arrays. Retain the local empty-input result of zero as a guard/adapter. It currently sorts a copy and averages the central pair for even lengths. |
| Private [simple_nano_species::median_dp:499](../../../src/main/nano/simple_nano_species.f90#L499) | Use the existing DP overload. This helper currently converts DP input to default real, sorts, and converts the median back to DP. Reuse removes that precision loss, so this is a numerical change to validate, not a bitwise-preserving substitution. |
| Public, pure [simple_state_assignment::median_of:296](../../../src/main/solve/simple_state_assignment.f90#L296) | Replace the insertion-sort implementation with shared median logic only after resolving the purity/API contract. The existing shared functions are not declared `pure`; a direct call inside this pure function will not compile. Preserve or deliberately revise the public interface and define empty-input handling. |

The existing [SP median:1352](../../../src/utils/math/simple_stat.f90#L1352) and [DP median:1376](../../../src/utils/math/simple_stat.f90#L1376) already implement the odd/even convention using selection on a copy. Neither supplies the atomic helper's empty-array policy. For the state-assignment API, first check public and test callers, then either remove an unnecessary purity promise or provide a genuinely pure shared implementation; do not simply add `pure` without checking its callees.

### A04. Sample standard deviation — two implementations

- **Local implementations:** contained [simple_atomic_species::sdev:145](../../../src/main/nano/simple_atomic_species.f90#L145) and [simple_nanoparticle::sdev_of:1829](../../../src/main/nano/simple_nanoparticle.f90#L1829).
- **Existing API:** [simple_stat::avg_sdev, SP implementation:79](../../../src/utils/math/simple_stat.f90#L79).
- **Evidence/action:** each computes the mean and the sample deviation with denominator `n-1`. A small value-returning adapter can call `avg_sdev` if expression-style call sites are useful. The shared routine already returns zero deviation for zero or one sample; `sdev_of` does too, while `sdev` currently divides by `size(vv)` before its denominator guard. Validate the intended empty-input behavior and preserve sample rather than population normalization.

### A05. Cached root-relative path trimming — five contained copies

- **Local implementations:** `get_relative_path_here` in [simple_starfile:334](../../../src/main/star/simple_starfile.f90#L334), [430](../../../src/main/star/simple_starfile.f90#L430), [585](../../../src/main/star/simple_starfile.f90#L585), and [simple_starproject_stream:192](../../../src/main/star/simple_starproject_stream.f90#L192), [261](../../../src/main/star/simple_starproject_stream.f90#L261).
- **Existing API:** [simple_fileio::get_relative_path:1181](../../../src/fileio/simple_fileio.f90#L1181), including its optional `trimlength` cache.
- **Evidence/action:** all cache the first substring position of a root and return the suffix beginning at that position, or the original path when no position is found. Use the shared API with the correct root and cache, with a character/string adapter if necessary. These routines keep the matched root in the output; they are not general filesystem-relative-path calculations.
- **Preserve:** fixed-length `XLONGSTRLEN` return behavior where required, cache lifetime/reset, and the parallel writer's thread-private `pathtrim` and string-operation critical sections. An overload taking characters may be preferable to repeated temporary string construction in a hot loop. Do not turn the cache into shared module state.

### A06. Destroying an allocatable image array — two loops

- **Local implementations:** contained [simple_extract_strategy::killimgbatch:497](../../../src/main/strategies/parallelization/simple_extract_strategy.f90#L497) and [simple_matcher_ptcl_io::killimgbatch:173](../../../src/main/strategies/search/simple_matcher_ptcl_io.f90#L173).
- **Existing API:** [simple_imgarr_utils::dealloc_imgarr:87](../../../src/main/image_processing/simple_imgarr_utils.f90#L87).
- **Evidence/action:** all check allocation, kill each image, and deallocate the array. Delegate the bodies to `dealloc_imgarr(imgs)` or `dealloc_imgarr(build%imgbatch)`. The public matcher wrapper may remain to preserve its builder-oriented API. The shared implementation runs the kill loop under OpenMP, whereas these copies are serial; check calling contexts and nested parallelism before switching.

### A07. Integer-array reversal inside a sorting helper — partial replacement

- **Local implementation:** final swap loop in [simple_cls_expansion_strategy::sort_order_by_weight_desc:1479](../../../src/main/strategies/parallelization/simple_cls_expansion_strategy.f90#L1479).
- **Existing API:** [simple_srch_sort_loc::reverse_iarr:509](../../../src/utils/math/simple_srch_sort_loc.f90#L509), exposed through `reverse`.
- **Evidence/action:** retain the subset/weight mapping and the existing integer heap sort; replace only the final manual reversal with `call reverse(order)`. The real-only `scores2order` is not an equivalent substitute for this integer-weight operation.

### A08. Paired value/weight sorting inside weighted percentile — partial replacement

- **Local implementation:** nested exchange-sort loops in contained [simple_class_compatibility::percentile_weighted_real:660](../../../src/main/class/simple_class_compatibility.f90#L660).
- **Existing API:** [simple_srch_sort_loc::hpsort_real_real:243](../../../src/utils/math/simple_srch_sort_loc.f90#L243), through `hpsort(sorted, wsorted)`.
- **Evidence/action:** the helper implements its own ascending paired sort before accumulating nonnegative weights. Reuse the paired sort; retain its percentile definition, mismatch fallback, zero-weight fallback, and clamped quantile. Equal-value ordering may change the order of floating-point weight accumulation, so validate boundary quantiles and ties. This finding concerns the sort, not a claim that the complete weighted-percentile API already exists.

### A09. Directory removal through an ad hoc shell command — contract check required

- **Local implementation:** contained [single_commanders_trajectory_solve3D::simple_rmdir_checked:340](../../../src/main/commanders/single/single_commanders_trajectory_solve3D.f90#L340).
- **Existing API:** [simple_syslib::simple_rmdir:948](../../../src/fileio/simple_syslib.f90#L948), bound to [remove_directory:1654](../../../src/fileio/simple_posix.c#L1654).
- **Evidence/action:** the helper checks for a directory and constructs `rm -rf ` plus an unquoted path through `exec_cmdline`; directory-tree removal already has a library owner. Use that owner after checking the required behavior, and preserve any desired workflow-specific error reporting. Despite its name, the local helper does not explicitly examine a removal status.
- **Not a blind substitution:** the current [recursive C implementation:1610](../../../src/fileio/simple_posix.c#L1610) uses `stat` on children and can traverse directory symlinks; shell `rm -rf` has different link behavior. Establish safe symlink semantics in the common API before adopting it for arbitrary trees. Paths with spaces/metacharacters and failed removals also deserve coverage. This is duplicated orchestration of a library capability, not a duplicated recursive algorithm.

## B. Consolidate copies or extend the owning API first

These findings have an existing counterpart or a closely related library operation, but not an appropriate drop-in public API for every caller. Suggested new entry points below are proposals, not existing functions.

### B01. Three-point parabolic peak offset — three exact helpers and one extended helper

- **Copies:** contained `interp_peak` in [simple_image_calc:1955](../../../src/main/image/simple_image_calc.f90#L1955), [simple_motion_align_hybrid:569](../../../src/main/motion/simple_motion_align_hybrid.f90#L569), and [simple_motion_align_nano:226](../../../src/main/motion/simple_motion_align_nano.f90#L226). They implement the same zero-on-small-denominator rule and `0.5*(alpha-gamma)/(alpha+gamma-2*beta)` offset.
- **Related copy:** [simple_tseries_graphene_subtr::interp_peak:209](../../../src/main/nano/simple_tseries_graphene_subtr.f90#L209) repeats the offset calculation and also updates angle and peak height. The image correlation code additionally repeats the formula inline, e.g. [simple_image_calc:1519](../../../src/main/image/simple_image_calc.f90#L1519).
- **Owner/action:** a small scalar interpolation routine in the math/interpolation layer, taking the three samples explicitly. Keep peak selection, periodic indexing, angular conventions, and any requirement that the centre exceed its neighbors in the caller. Preserve each caller's strict/non-strict `TINY` comparison when reviewing the inline copies.

### B02. RGB packing — two exact contained helpers

- **Copies:** [simple_oris_io::rgb_code:354](../../../src/main/ori/simple_oris_io.f90#L354) and [simple_ptcl_sieve::rgb_code:2134](../../../src/main/sieve/simple_ptcl_sieve.f90#L2134).
- **Evidence:** both clamp channels to 0–255, pack them with shifts 16/8/0, then divide by `16777215.0` into a real value.
- **Owner/action:** put one pure packing routine beside the image/color-output utilities. Do not make orientation I/O depend on sieve code or vice versa. Preserve clamping and the normalized packed representation, which is not luminance or three independent normalized channels. Low-cost, low-priority consolidation.

### B03. Extraction-box containment — two exact copies

- **Copies:** contained [simple_extract_strategy::box_inside:507](../../../src/main/strategies/parallelization/simple_extract_strategy.f90#L507) and private [simple_reextract_strategy::box_inside:719](../../../src/main/strategies/parallelization/simple_reextract_strategy.f90#L719).
- **Owner/action:** expose one predicate with the extraction geometry owner used by both strategies. It translates a zero-based corner with `coord+1`, uses inclusive `from+box-1`, and checks only the first two image dimensions. Test boundary boxes and preserve that convention; do not substitute a centre-coordinate or three-dimensional check.

### B04. Filename-based file size — two exact helpers

- **Copies:** contained [simple_volassemble_strategy::trail_chain_file_size:283](../../../src/main/strategies/parallelization/simple_volassemble_strategy.f90#L283) and private [simple_frozen_accum::file_size:386](../../../src/main/volume/simple_frozen_accum.f90#L386).
- **Owner/action:** one file-I/O routine returning `integer(kind=8)` bytes and `-1` when absent. [simple_fileio::funit_size:346](../../../src/fileio/simple_fileio.f90#L346) is related but accepts an open unit and is not a direct replacement. Preserve the sentinel and large-file range; decide how an inquiry error is represented when defining the common API.

### B05. Corrupt-stack header validation — two exact routines

- **Copies:** [simple_procimgstk::raise_exception_imgfile:12](../../../src/main/image_processing/simple_procimgstk.f90#L12) and private [simple_stackops::raise_exception:15](../../../src/main/image_processing/simple_stackops.f90#L15).
- **Evidence:** the same `n < 1 .or. any(ldim == 0)` predicate, diagnostics, and `THROW_HARD('procimgfile exception')`.
- **Owner/action:** the first routine is already public, but choose a common stack/header-I/O owner instead of introducing a large image-processing-module dependency just for this check. Preserve caller diagnostics and the existing predicate; changing it to reject every negative dimension is a separate behavior change.

### B06. Low-pass stage scale information — two exact contained routines in one module

- **Copies:** `calc_scaleinfo` in [simple_estimate_ssnr::lpstages:425](../../../src/utils/filter/simple_estimate_ssnr.f90#L425) and [lpstages_fast:487](../../../src/utils/filter/simple_estimate_ssnr.f90#L487).
- **Owner/action:** one private procedure in `simple_estimate_ssnr`, with box, sampling, and the stage record passed explicitly. Both compute target sampling, call `autoscale` with `minbox=88`, set the same translation bound, and record whether cropping occurred. This needs no new module or public API.

### B07. Particle-update coverage — two near copies plus one policy variant

- **Copies:** contained [simple_commanders_refine3D::read_update_coverage:1627](../../../src/main/commanders/simple/simple_commanders_refine3D.f90#L1627) and [simple_refine3D_states_controller::read_update_coverage:854](../../../src/main/simple_refine3D_states_controller.f90#L854).
- **Variant:** contained [simple_commanders_solve3D::read_update_coverage:2208](../../../src/main/commanders/simple/simple_commanders_solve3D.f90#L2208).
- **Evidence/action:** all obtain `state` and `updatecnt`, count active/updated rows, and derive missing rows. The first two also duplicate project-segment loading, allocation cleanup, missing-field rejection, and reporting, with different label sources. Factor the common calculation into the orientation/project layer or a focused shared workflow service. Keep project ownership and loading explicit.
- **Do not merge policy blindly:** solve3D permits a single-state path without updates and changes `nmissing` to zero when every active particle will be reconstructed. Its treatment of an absent `updatecnt` differs from the other two. Reporting and workflow-specific rejection should remain with each controller/commander.

### B08. PCG raw-part cleanup — two copies with renamed locals

- **Copies:** [simple_rec3D_service::remove_pcg_raw_files:357](../../../src/main/strategies/parallelization/simple_rec3D_service.f90#L357) and [simple_refine3D_strategy::remove_pcg_raw_files:660](../../../src/main/strategies/parallelization/simple_refine3D_strategy.f90#L660).
- **Evidence/action:** both loop over states, even/odd halves, and parts; generate the same raw accumulator name; delete both it and its `.tmp` companion. Share this in reconstruction/strategy support, using explicit counts and filename-width inputs where practical. Both already reuse `refine3D_pcg_raw_accum_fname`; the duplicated part is the cleanup policy around it. Do not move file deletion into a definitions module.

### B09. Stateless PCG timing — two copies of a two-procedure pair

- **Copies:** `pcg_tic`/`pcg_toc` in [simple_pcg_solver:330](../../../src/main/opt/simple_pcg_solver.f90#L330) and [simple_reconstructor_pcg:239](../../../src/main/volume/simple_reconstructor_pcg.f90#L239).
- **Related API:** [simple_timer::tic:80](../../../src/utils/simple_timer.f90#L80), [tdiff:97](../../../src/utils/simple_timer.f90#L97), [toc:114](../../../src/utils/simple_timer.f90#L114).
- **Owner/action:** provide a stateless clock-read/elapsed-time API in `simple_timer`, then use it from both PCG owners. The reconstructor explicitly documents concurrent even/odd solves. Existing `tic`/`toc` modify `last_time_point`, and `tdiff` may initialize the global rate cache, so replacing the helpers with those APIs without examining concurrency would discard a deliberate property. Retain `real(dp)` elapsed values, including under `TIMER_I4` configurations.

### B10. Double-precision cross product — missing overload

- **Local implementation:** pure [simple_ori_utils::cross3:204](../../../src/main/ori/simple_ori_utils.f90#L204).
- **Related API:** [simple_math::cross:1121](../../../src/utils/math/simple_math.f90#L1121), currently default-real and not declared pure.
- **Owner/action:** add a DP/pure-capable cross-product interface to the math owner, then reuse it. Both implement the same three-component determinant. Converting the orientation input to default real to use the present routine would lose precision; simply replacing the call would also encounter the purity mismatch where a pure caller requires it.

### B11. Double-precision descending ordering — shared sorting capability gap

- **Local implementations:** [simple_nano_species::sort_desc:512](../../../src/main/nano/simple_nano_species.f90#L512), an insertion sort of indices, and [simple_flex_pca_targets::sort_block_desc:458](../../../src/main/flex/states/simple_flex_pca_targets.f90#L458), a selection sort that swaps eigenvector columns with eigenvalues.
- **Related API:** [simple_srch_sort_loc::scores2order:755](../../../src/utils/math/simple_srch_sort_loc.f90#L755) and its `hpsort` overloads are default-real/integer, without the needed DP ordering overload.
- **Owner/action:** extend sorting at the owning library layer, then keep only the column/permutation application in FLEX. These helpers have related goals, not identical contracts: nano's strict comparisons preserve input order among ties, whereas the FLEX selection swaps do not promise stability. Preserve or explicitly specify the tie rule and DP precision. Do not round values to default real to fit the current API.

### B12. Filename extension stripping — overlapping capability, stronger local semantics

- **Local implementation:** contained [simple_commanders_resolest::strip_extension:208](../../../src/main/commanders/simple/simple_commanders_resolest.f90#L208).
- **Related APIs:** [simple_fileio::get_fbody:484](../../../src/fileio/simple_fileio.f90#L484) and [fname2ext:541](../../../src/fileio/simple_fileio.f90#L541).
- **Owner/action:** extend the file-I/O owner with the required last-extension operation and use it here. The local helper preserves extensionless names, dotfiles, and dots occurring only in the directory portion. `get_fbody` searches for a supplied suffix, while `fname2ext` scans the full path for its last dot; composing them is not a proven replacement for these cases. This is an ownership/API-gap candidate, not an exact duplicate.

### B13. Motion-shift RMSD — two implementations differing only in receiver type

- **Copies:** [simple_motion_align_hybrid::calc_rmsd:635](../../../src/main/motion/simple_motion_align_hybrid.f90#L635) and [simple_motion_align_nano::calc_rmsd:283](../../../src/main/motion/simple_motion_align_nano.f90#L283).
- **Owner/action:** one motion utility taking the two shift arrays. Both sum squared two-coordinate shifts and divide by frame count before taking the square root. Preserve per-frame normalization; a generic elementwise RMS over the flattened array would introduce a factor of `sqrt(2)` difference. Small, secondary cleanup.

### B14. Motion correlation weights — two larger supporting routines

- **Copies:** [simple_motion_align_hybrid::gen_weights:580](../../../src/main/motion/simple_motion_align_hybrid.f90#L580) and [simple_motion_align_nano::gen_weights:236](../../../src/main/motion/simple_motion_align_nano.f90#L236).
- **Evidence/action:** the 46-line bodies differ in their receiver type, with the same Fourier limits, B-factor calculation, cosine tapers, OpenMP loop, and squared combined weights. Extract the numerical operation into motion support with explicit geometry/filter inputs. Keep each type's lifecycle method as an adapter. This is larger than the initial small-helper examples but a clear nearby numerical duplicate; preserve rectangular-image frequency scaling and squared weights rather than substituting a superficially similar image bandpass operation.

## Commander-specific disposition

This table collects the commander findings so the initial cleanup need not start by traversing the whole library.

| Commander/helper | Disposition |
| --- | --- |
| `single_commanders_trajectory_solve3D::sort_order` | A01: use `dists2order`. |
| `simple_commanders_refine3D::read_update_coverage` | B07: share calculation/loading where appropriate; retain workflow failure policy. |
| `simple_commanders_solve3D::read_update_coverage` | B07: reuse the common calculation and preserve its single-state exception. |
| `simple_commanders_resolest::strip_extension` | B12: extend filename utilities before replacing. |
| `single_commanders_trajectory_solve3D::simple_rmdir_checked` | A09: use the directory-removal owner after checking link/error semantics. |
| `simple_commanders_mkcavgs::rank_of` | Optional intrinsic simplification below; retain its failure check. |
| `shuffle_ints`, checked filesystem wrappers, `rank_cavgs`, formatting/report helpers | Keep when they adapt existing APIs or express local behavior; examples below. |

## C. Similar names and small helpers that should not be counted as duplicate algorithms

| Candidate inspected | Why it is not a direct replacement finding |
| --- | --- |
| [simple_commanders_mic_subset::csv_field:158](../../../src/main/commanders/simple/simple_commanders_mic_subset.f90#L158) vs [simple_string_utils::csv_field:713](../../../src/utils/text/simple_string_utils.f90#L713) | The commander **encodes/quotes** one field and doubles embedded quotes. The library routine **extracts** a comma-separated field by number and does not implement that encoder. If sharing is wanted, introduce a separately named CSV encoder in text utilities; do not call the current parser instead. |
| [simple_commanders_mkcavgs::shuffle_ints:355](../../../src/main/commanders/simple/simple_commanders_mkcavgs.f90#L355) | Already delegates to `ran_tabu%shuffle` and manages its lifecycle plus the small-array guard. It does not duplicate a shuffle algorithm. |
| [simple_commanders_resolest::gq:511](../../../src/main/commanders/simple/simple_commanders_resolest.f90#L511) | Already calls shared `fortran_quote`; it is a trimmed-return adapter. |
| [simple_refine_motion_model_strategy::gnu_quote:1282](../../../src/main/strategies/parallelization/simple_refine_motion_model_strategy.f90#L1282) | Escapes double quotes, backslashes, backticks, and control characters for its plotting syntax. Neither POSIX argument quoting nor `fortran_quote` has the same contract. Keep separate from A02. |
| [simple_commanders_class_student::r2s:558](../../../src/main/commanders/simple/simple_commanders_class_student.f90#L558) vs [simple_string_utils::realsp2str:424](../../../src/utils/text/simple_string_utils.f90#L424) | The former specifies fixed decimal places and trims a report field; the latter uses the string type's general numeric conversion. Formatting is observable output, so this is not a drop-in replacement. |
| [trajectory `chdir_checked`:332](../../../src/main/commanders/single/single_commanders_trajectory_solve3D.f90#L332) and `require_file` nearby | Workflow adapters around shared filesystem operations, with explicit status or existence checks and workflow-specific diagnostics. Reusing the underlying library is already the pattern. The separate removal helper is A09. |
| [class-student `enter_subdir`:362](../../../src/main/commanders/simple/simple_commanders_class_student.f90#L362) and [pose-student `enter_subdir`:312](../../../src/main/commanders/simple/simple_commanders_pose_student.f90#L312) | Exact six-line orchestration copies, already calling shared mkdir/chdir operations. Low-value extraction on its own; a shared student workflow owner is reasonable if these workflows are refactored together. |
| [simple_cavg_quality_stats::median_by_state:113](../../../src/main/cavg_quality/simple_cavg_quality_stats.f90#L113) | Already packs the mask and calls library `median`, adding shape and empty-selection handling. This is the adapter pattern to retain. |
| [simple_cavg_quality_feats::normalize_robust:465](../../../src/main/cavg_quality/simple_cavg_quality_feats.f90#L465) and [simple_cavg_quality_relations::normalize_relation_feature:160](../../../src/main/cavg_quality/simple_cavg_quality_relations.f90#L160) | Both already use shared median/MAD; their array rank, eligibility masks, clipping and rejected-row handling express feature policy. Do not replace them with a generic normalization that omits those rules. |
| [extract `prepimgbatch`:476](../../../src/main/strategies/parallelization/simple_extract_strategy.f90#L476) vs [matcher `prepimgbatch`:135](../../../src/main/strategies/search/simple_matcher_ptcl_io.f90#L135) | Related allocation patterns, but matcher reuse also checks box and pixel size and supports overrides. A06 is the shared destruction operation; allocation-policy unification needs a separate contract. |
| [solve3D `rank_cavgs`:567](../../../src/main/commanders/simple/simple_commanders_solve3D.f90#L567) vs [stream_pool2D `rank_cavgs`:1092](../../../src/main/stream/pool2D/simple_stream_pool2D.f90#L1092) | Different ranking criteria, project segments and outputs; both delegate the substantive operation. The same helper name is not evidence of duplicate ranking code. |
| [CTF `poly2val`:1379](../../../src/main/ctf/simple_ctf_estimate_fit.f90#L1379) vs [motion `poly2val`:785](../../../src/main/motion/simple_motion_align_hybrid.f90#L785) | One evaluates a ten-term bivariate polynomial in DP, the other a four-term univariate cubic returning default real. Both already use `dot_product`; merging them solely by name would obscure the basis contract. |
| [`percentile_real`:637](../../../src/main/class/simple_class_compatibility.f90#L637) vs [`student_boot_quantile`:97](../../../src/main/strategies/parallelization/simple_strategy_student_common.f90#L97) | Different order-statistic conventions, and the student procedure includes deterministic bootstrap resampling. A generic percentile substitution requires agreeing on the rank definition first. |

An optional, small intrinsic simplification is [simple_commanders_mkcavgs::rank_of:364](../../../src/main/commanders/simple/simple_commanders_mkcavgs.f90#L364): `findloc(ranked, pind, dim=1)` supplies the first exact-match index and zero when absent. Keep the existing `THROW_HARD` check. `simple_srch_sort_loc::find`/`locate` are ordered-table/nearest or bracket searches and are not equivalents. This is separate from the 23 library-reuse/consolidation families.

Do not collapse independent test calculations into the production routine they validate. Likewise, similar type-bound constructors, destructors, execution routers, and stream-stage lifecycle methods are not automatically utility duplicates; their receiver types and ownership are part of the behavior.

## Suggested implementation order and validation

1. **Existing API reuse:** A01 and A02, followed by A03/A04 with their explicit empty-input, precision and purity decisions. Keep focused imports and adapters where they preserve an interface.
2. **Shared I/O/lifecycle operations:** A05/A06/A09 and B03/B04/B05/B08. Preserve caches, project/resource ownership, threading, filename conventions, symlink semantics and failure behavior.
3. **Numerical consolidation:** B01, B06, B10/B11 and B13/B14. Use closed-form or independent expectations for numerical tests, cover ties/boundaries, and compare affected workflows when claiming unchanged behavior.
4. **Policy-sensitive work:** B07, B09 and B12 need the common API contract written down before substitution. A07/A08 can be small separate cleanup changes with the stated ordering checks.

For an implemented change, follow the existing test policy: use the owning area's registered tester/sub-suite, build incrementally, run the fast gate and relevant area/workflow entries. Do not create standalone scratch test programs or a new CMake target to test these helpers. Tests should cover observable contracts, not merely reproduce the helper body.

**Inventory validation:** source bodies and proposed counterpart implementations were inspected; relative links and source anchors were checked. No runtime equivalence claim is made. No Fortran was edited and no build or tests were run for this documentation-only inventory. The previous constant-review build/gate results do not validate any of the proposed helper refactors.
