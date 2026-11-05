# Unreferenced code inventory

Status: completed review; archived on 2026-11-04. The two approved changes,
native EPU XML reading and loading the test molecules from PDB files, were
delivered in `ef5fe84d5`. The other deletion proposals were withdrawn by the
owner's ruling below. This record authorizes no further code removal.
Its measurements and reachability claims describe the dated audit, not a new
validation or deletion list.

Historical inventory, 2026-11-02, against master at `703a7e825`. **Identification only: nothing was deleted by the audit.** The
question was which code we compile and test over and over although nothing uses it. Appendices A-C list
what the analysis found; appendix D says how to rerun it.

> **Owner ruling 2026-11-02: the lists of sections 3, 4 and 6 and appendices A-C are not deletion lists.**
> They contain mission-critical code. That no procedure in master calls a routine does not make it dead,
> and this inventory was wrong to call it that and to propose deleting it. Nothing on these lists is to be
> removed on the strength of this analysis. What stands: section 1 (the XML library, replaced) and
> section 2 (the molecule data, moved to PDB files read at run time).

## Summary

| category | procedures | lines | share of compile time |
| --- | ---: | ---: | ---: |
| SIMPLE production code, live | 6 595 | 214 000 | |
| SIMPLE code nothing reaches | 442 | 8 600 | about 2% |
| SIMPLE production code only the tests reach (without the molecule data) | 247 | 4 000 | about 1% |
| embedded molecule data, only the tests reach | 2 | 56 500 | 0.4% |
| vendored libraries, unreachable parts (XML 25 000, json-fortran 4 900, others 600) | | 30 500 | |
| the whole XML library (one use, see 1) | | 53 600 | 1.7% |

The finding that matters most: **dead code is a small part of the compile time.** Of the 1 001 s a clean
build spends compiling (sum over 935 jobs), SIMPLE's production code takes 69%, the test code 28% and all
vendored libraries together 3%. SIMPLE's own dead code is about 2%, and the slowest files are slow because
of what they import, not because of their size (`simple_commanders_ori.f90` compiles for 11.9 s; its two
dead commanders are a quarter of its lines). Removing everything in this inventory makes the tree smaller
and the reviews and tests easier; it does not make the build much faster. The levers on build time are
the test code (`--exclude-tests` already drops it from a developer build) and the umbrella imports of the
commanders, which `doc/policies/compile_time_policy.md` already targets.

What I would do, in order of payoff:

1. **Replace the XML library** (67 files, 53 600 lines, 1.7% of the compile time) by one routine of about
   40 lines: SIMPLE uses it only to read two numbers from EPU metadata (section 1). *Done 2026-11-02.*
2. ~~Delete the code nothing reaches~~ *Withdrawn: see the ruling above.*
3. ~~Rule on the owner decisions below~~ *Withdrawn with it, except decision 2.*

Decisions for you:

1. The GPU and OpenMP-offload code paths are unreachable (`simple_reconstructor_openmpoffload.f90`, the
   cuFFT/GPU distance routines of `simple_polarft_corr.f90`, about 600 lines). Retire them, or are they
   the start of a GPU port you want to keep?
2. The 56 500-line `simple_molecule_data.f90` (the atoms of 1JYX and 6VXX) is test data compiled into the
   core library. Move it to PDB files the tests read, or leave it? *Ruled 2026-11-02: move it to PDB
   files in the repository, read when needed.*
3. Production code that only the tests reach (section 4, appendix B): retire the features among it with
   their tests (Canny edge detection, the batch polar projector, ...)? The test oracles and helpers stay.
4. 65 registered programs have no test, no internal caller and no NICE reference (appendix C). They are
   live code, so only you can say which are obsolete.

## Method

A reachability analysis of every Fortran source the build compiles (`src/`, `production/`, the two
generated sources `simple_args.f90` and `simple_ui_default_values.f90`, the `#include`d macro files):

- **Nodes**: every procedure, derived type, generic interface and module variable.
- **Roots**: the production programs (`simple_exec`, `single_exec`, `simple_stream`, `simple_private_exec`,
  `simple_persistent_worker`), separately `simple_test_exec`, and every `bind(c)` procedure (444).
- **Edges**: an identifier in a reached node reaches every node of that name in a module visible from it
  (the node's own module and, transitively, what it and its procedure `use`, so re-exports count);
  `%name` reaches the bindings of that name of every reached type, generic bindings their specifics;
  generic interfaces their specifics; a reached type its final procedures and defined operators and
  assignments, a reached module its operator and assignment interfaces; use-renames are followed.
  Both branches of every `#ifdef` count, so code only a GPU or Windows build uses is not called dead.
- **Classes**: *live* (reached from a production program), *test-only* (reached only from
  `simple_test_exec`), *unreachable* (reached from neither).

The analysis can only err towards calling dead code live (a dead procedure sharing a name with a visible
live one is kept), not the reverse, except through a parsing error. Three checks against that:

- **The linker.** A build with `-ffunction-sections -Wl,--gc-sections`, then the symbol tables of all six
  executables. Of the 442 unreachable SIMPLE procedures, 178 are in no executable. 255 are type-bound
  procedures of live types: a type's vtable references every binding, so the linker cannot drop an uncalled
  method; the analysis found no `%name` call of them in reachable code. The remaining 9 the linker keeps are
  called only from those uncalled methods (checked one by one). No procedure the analysis calls unreachable
  is kept by the linker for any other reason.
- **Spot checks** of the largest candidates by hand (`rg` of every name): no caller outside dead code, or
  only mentions in strings, logs and umbrella import lists.
- **No implicit calls**: SIMPLE has no defined derived-type I/O and no `declare reduction`; C code calls
  Fortran only through `bind(c)`; procedure pointers are assigned by name, which the analysis follows.

What it does not see: programs a user starts by hand are live by definition (section 6 lists the ones
nothing else touches), and a binding called only through a name the analysis cannot resolve would be
missed. Neither weakens the "unreachable" verdicts above.

## Where the compile time goes

One profiled clean build of `703a7e825` (`scripts/profile_build.sh clean`, `-j16`, Release, tests on;
the machine was also running the nightly suite, so the shares are reliable, the absolute times are not
a benchmark):

| part | compile s | share |
| --- | ---: | ---: |
| SIMPLE production code | 691 | 69.0% |
| test code (testers, test commanders, test routers and UI) | 278 | 27.8% |
| XML library | 17 | 1.7% |
| starfile, gnufor, json-fortran, jpeg, CPlot2D, unix, curl, tiff | 16 | 1.6% |
| total | 1 001 | wall time 101 s |

The dead code's share, by lines within each file: unreachable code about 18 s, test-only production code
about 10 s, the molecule data 3.8 s. A compile-time gain may be claimed only from a like-for-like
`scripts/profile_build.sh compare` after the change (compile-time policy).

## 1. The XML library

**Done 2026-11-02.** `src/fileio/simple_epu_xml.f90` (`read_epu_beamshift`, with the `unit_core` sub-suite
`EPU XML`) replaces the library at the three call sites, which now log a file without a BeamShift and go
on (the DOM lookup crashed on one); `src/extlibs/xml` and its CMake glob are gone.

`src/extlibs/xml` (FoX: 67 files, 53 600 lines, 17 s) has three users, and all three do the same thing:
parse an EPU metadata file and read `BeamShift/a:_x` and `BeamShift/a:_y`
(`simple_relion.f90:588`, `simple_starproject.f90:918`, `simple_preprocess_strategy.f90:305`, the same
seven lines each). 25 000 of the library's lines are unreachable even so, 24 of its files entirely.

Proposal: one routine in `src/fileio` that reads the file, finds `<BeamShift>` and returns the two numbers
between `<a:_x>`/`</a:_x>` and `<a:_y>`/`</a:_y>`, called from the three places; then delete
`src/extlibs/xml` and its CMake glob. Needs one or two real EPU XML files as a unit-test fixture.

## 2. The embedded molecule data

**Ruled 2026-11-02: move the two molecules to PDB files in the repository, read when needed. Done the
same day:** `testdata/pdb/1jyx.pdb` and `6vxx.pdb` (4.8 MB, written from the old array constructors and
read back bit-identical, every field of all 60 210 atoms); `simple_molecule_data.f90` is now a reader of
about 90 lines with the same interface.

`src/main/nano/simple_molecule_data.f90`, 56 457 lines: `betagal_1jyx` and `sars_cov2_spkgp_6vxx` return
the atoms of two PDB entries as array constructors. Only tests call them (the stream, workflow, addon and
continuous-refinement test commanders, the state fixture, `simple_pdb2mrc_tester`,
`simple_strategy3D_cont_1jyx_tester`); no production executable contains them. Production uses only the
`molecule_data` type, through `atoms%new_from_molecule`. It costs 3.8 s per clean build and a recompile
whenever the file is touched. Options: keep the type in `simple_atoms`, move the two functions to the
test code, or replace them by two PDB files the tests read.

## 3. SIMPLE code nothing reaches

442 procedures, 8 600 lines, 128 files; appendix A lists every one with its linker verdict.

**Whole files** (nothing in them is reached):

| file | lines | what it is |
| --- | ---: | --- |
| `utils/math/simple_testfuns.f90` | 861 | 32 optimiser test functions; not even the optimiser tests use them |
| `utils/comm/simple_stream_communicator.f90` | 461 | `stream_http_communicator`, a heartbeat/checksum client no one instantiates |
| `main/volume/simple_reconstructor_openmpoffload.f90` | 388 | OpenMP-offload reconstruction; no file uses the module, with or without `USE_OPENMP_OFFLOAD` (decision 1) |
| `extlibs/unix/src/unix_utsname.F90`, `unix_ipc.F90` | 85 | wrappers nothing uses |

**Commanders no router reaches** (12, 472 lines; none is a registered program, so no user-visible program
is broken): `check_states`, `oriconsensus` (`simple_commanders_ori.f90`), `map_cavgs_states`
(`simple_commanders_cavgs.f90`), `trajectory_backgr_subtr` (`single_commanders_trajectory.f90`),
`map_validate`, `model_validate_eo`, `map2model_fsc` (`simple_commanders_atoms.f90`),
`simulate_subtomogram` (`simple_commanders_sim.f90`), `concatenate_projects`
(`simple_commanders_project_core.f90`), `random_rec_commander` (`simple_commanders_rec.f90`),
`commander_mkdir` (`simple_commanders_misc.f90`), and the internal `execute_solve2d_staged`
(`simple_commanders_solve2D.f90`). With them goes `sp_project%append_project`, which only
`concatenate_projects` calls.

**GPU paths** (decision 1): `gen_many2many_euclids_cufft`, its kernel, `gen_many_euclids_gpu` and the
bindings that name them in `simple_polarft_calc.f90`, with the offload reconstructor above.

**The rest**, largest first: the 1D hierarchical-clustering and filter helpers of `simple_math.f90`
(`hac_1d_fast`, `hac_med_thres`, `hac_1d`, `savitzkygolay_filter`; 572 lines, half the file's
procedures), unused correlations of `simple_image_calc.f90` (376), 88% of `simple_timer.f90`'s profiling
API (339), unused filters of `simple_image_filt.f90` (`gabor_filter2d`, the even/odd non-local means,
`whiten_noise_power`, `glcm`; 274), statistics of `simple_stat.f90` (268), five methods of
`simple_ui_program` including `print_ui` and `write2json` (118), and 255 uncalled methods of live types
across image, ori, project, reconstructor and structs.

## 4. Production code only the tests reach

247 procedures, 4 000 lines in 95 files, plus the molecule data (appendix B). Two kinds:

- **Test oracles and helpers**, which belong with the tests and should stay (possibly moved into the
  testers): the reference terms of `simple_cartft_calc.f90` (`*_full_disk`, `objective_gradient`), the
  kernel and accumulator checks of `simple_reconstructor_pcg.f90` (`measure_kernel_scale`,
  `compare_raw_accum`, ...), `files_identical` and `filelength` in `simple_fileio.f90`, the
  `simple_forked_process` queries, `simple_ui%validate_ui_json` and `make_test_ui`.
- **Features with tests but no production caller**, which go with their tests if retired (decision 3):
  Canny and Sobel edge detection and the FDR peak thresholds of `simple_segmentation.f90` (320 lines), the
  batch polar projector `simple_projector_pft_batch.f90` (103), unused getters and setters of the
  starfile wrapper (177), `write_ptcl2d_table_parallel` (147), parts of `simple_particle_extractor.f90`
  (135), the Gaussian images of `simple_image_geom.f90` (110), the shift-map routines of the polar code.

Sorting the 95 files into these two kinds needs a look at each; appendix B is the list to do it from.

## 5. Vendored libraries

Apart from XML: json-fortran has 4 900 unreachable lines of 9 800, gnufor 270, curl 196 (curl itself is
used, by `simple_http_post` for the GUI). A vendored library is normally used in part; I would leave them
whole, except the XML library, which can go entirely.

## 6. Programs nothing else touches

165 registered programs. 65 have no test (by name in the test code), no internal caller (`prg` set by
another commander) and no NICE reference (appendix C, with the number of mentions in `doc/` and
`scripts/`). They are live code, and a user may run any of them; this is the list on which only you can
mark the obsolete ones. Several are in daily use (`vizoris`, `oriops` and `oristats`, for example, per
the owner); a heuristic like this one cannot tell. A retired program typically takes a commander, a
router case, a UI entry and their keys, as the PCA retirement did.

## 7. Small items

- 156 module-level constants and variables nothing reads, most in `simple_fftw3` (33, the C interface),
  `simple_defs` (28), `simple_defs_fname` (26) and `simple_defs_stream` (16); negligible compile cost.
- 17 unused derived types: the 11 commander types of section 3, `stream_http_communicator`, `sigma_array`
  and four FFTW interface types.

## 8. Found on the way

The final verification run of the nightly-suite triage commit (on `build/`, started 07:59) hit a hang in
`lib_stream`: the stream-chain sub-suite stopped after `MOVIES TO PROCESS: 0` at 08:07:56 and waited until
the 7 200 s timeout. The suite passed in the three earlier library runs (606-622 s). No machine activity of
mine overlapped the hang. The stream chain tests are Joe's and Ruben's; this is worth passing to them, as
tonight's CI may show the same timeout. The rest of that run: the other seven library entries passed, and
all 18 high-level entries passed (10 394 s). `build/` was rebuilt as a debug build of the retirement tree
at 09:22, mid-run, so the high-level tier verified that debug build, and its identity checks against the
earlier Release runs differ at rounding level from the first stage on, as debug and Release numerics do.

## 9. Dead by three signs (2026-11-02, after the ruling)

A second pass asked for three signs at once instead of "no caller": nothing imports the module, or the
code is no one's interface (not public, not a method, not behind a public generic); git history shows
when and in which commit the last caller went; and the code is not maintained (no edit in 120 days
apart from sweeping commits). Public procedures never count as dead here, only as "ask". Evidence for
all 150 candidates: `/data2/agent_runs/dead_code_2026-11-02/history_evidence.json` on the build box.

**Dead** (1 284 lines):

| code | lines | evidence |
| --- | ---: | --- |
| `utils/math/simple_testfuns.f90` | 861 | no file imports it; its last user went in "standalone test programs and the utils area retired" (2026-09-24) |
| commander types `commander_check_states`, `commander_map_cavgs_states`, `commander_simulate_subtomogram`, `commander_map2model_fsc` | 232 | their programs were removed on purpose in 8126d7911 "Removed empty applications UIs and unused commanders" (2025-12-16); the types stayed behind |
| `random_rec_commander` | 18 | its program left the UI in 2021 |
| `commander_map_validate`, `commander_concatenate_projects`, `commander_trajectory_backgr_subtr` | 112 | never registered or routed under any name |
| `commander_mkdir` | 7 | router case removed 2025-12-11, no program |
| private `bessi0_dp`, `bessi1`, `bessi0f` (`simple_kbinterpol.f90`) | 54 | superseded by the faster KB routines (2026-01-31; `bessi0f` 2021) |

**Ask** (129 procedures, 2 345 lines): every public procedure without a caller (the statistics routines,
`has_ldim_nptcls`, the starfile wrapper's setters, the timers, the memory estimator, ...) and everything
that lost its last caller in the last six months. The one that looks like an accident:
`simple_stream_communicator.f90` lost its only user (nice-lite, added 2026-09-09) in the stream refactor of
2026-10-05; a question for Joe and Ruben. **Keep** (9): edited in the last four months. The 255 methods of
live types were not assessed: they are interface by this rule.

## Proposed next steps

1. The XML replacement, as its own diff (a fixture from a real EPU file, a unit test, the three call sites,
   the library removed).
2. One diff deleting appendix A minus the GPU paths: build, fast gate, library tier, and the high-level
   entries that touch the edited areas; compile time compared like for like with `profile_build.sh`.
3. Your rulings on decisions 1-4, each a follow-up diff.

## Appendix A: SIMPLE code nothing reaches

Tags: none = absent from every executable; **M** = method of a live type, kept only by its vtable;
**K** = kept only because an M calls it.

| file | lines | procedures |
|---|---:|---|
| `utils/math/simple_testfuns.f90` | 668 | `get_testfun`, `testfun17`, `testfun3`, `testfun18`, `testfun11`, `testfun10`, `testfun26`, `testfun14`, `testfun2`, `testfun13`, `testfun1`, `testfun6`, `testfun5`, `testfun19`, `testfun12`, `testfun9`, `testfun8`, `testfun7`, `testfun4`, `testfun32`, `testfun16`, `testfun15`, `testfun30`, `testfun23`, `testfun22`, `testfun20`, `testfun31`, `testfun29`, `testfun28`, `testfun27`, `testfun25`, `testfun24`, `testfun21` |
| `utils/math/simple_math.f90` | 572 | `hac_1d_fast`, `hac_med_thres`, `hac_1d`, `savitzkygolay_filter`, `calc_score_thres`, `round2odd`, `great_circle_samples`, `nextpow2`, `sinc`, `vis_3dreal_mat`, `vis_3dinteger_mat`, `vis_2dreal_mat`, `vis_2dinteger_mat`, `shft`, `rotate_vec`, `gaussian3d` K, `lcm`, `gaussian1d` |
| `main/image/simple_image_calc.f90` | 376 | `frame_ref_phase_correlation` M, `corr_shifted` M, `fcorr_shift` M, `avg_loc_sdev` M, `loc_var3d` M, `loc_var` M, `gen_argtransf_comp` M, `dead_hot_positions` M, `masked_subtr_corr` M, `checkimg4nans` M, `masked_subtr_sqsum` M, `kurt` M, `masked_sqsum` M |
| `utils/comm/simple_stream_communicator.f90` | 369 | `curl_request` M, `background_heartbeat` M, `evaluate_checksum` M, `create` M, `term` M, `send_jobstats` M, `join_background_heartbeat` M, `send_heartbeat` M, `get_json_arg_5` M, `update_json_4` M, `update_json_3` M, `update_json_2` M, `get_json_arg_4` M, `get_json_arg_3` M, `get_json_arg_2` M, `get_json_arg_1` M, `update_json_1` M, `add_to_json_5` M, `add_to_json_2` M, `add_to_json_1` M, `remove_from_json_if_present` M, `destroy_arg` M, `add_to_json_4` M, `add_to_json_3` M, `stop_status` M, `exit_status` M, `arg_associated` M, `active_status` M |
| `main/volume/simple_reconstructor_openmpoffload.f90` | 354 | `insert_all_slices`, `insert_slices_batch`, `calc_3drec_gpu`, `update_ptcls_arrays` |
| `utils/simple_timer.f90` | 339 | `timer_profile_setup`, `timer_profile_report`, `timer_loop_end`, `timer_profile_break`, `timer_profile_start`, `tocprint`, `timer_loop_start`, `in_timer_loop`, `toc_i4`, `tickrate`, `tdiff_i4`, `tickrate_i4`, `tic_i4`, `reset_timer` |
| `main/image/simple_image_filt.f90` | 274 | `gabor_filter2d` M, `nlmean3d_eo` M, `whiten_noise_power` M, `nlmean2d_eo` M, `glcm` M |
| `utils/math/simple_stat.f90` | 268 | `prepare_tree_sub_distmat`, `moment_serial`, `kurtosis_2` K, `kurtosis_1`, `merge_smats`, `robust_scaling`, `medoid_ranking_from_smat`, `norm_corr_8`, `norm_corr`, `robust_sigma_thres`, `robust_normalize_minmax`, `robust_z_scores`, `robust_normalization`, `corr2distweight`, `spear` |
| `main/pftc/simple_polarft_corr.f90` | 246 | `gen_many2many_euclids_cufft` M, `gen_many_euclids_gpu` M, `gen_many2many_euclids_cufft_kernel`, `gen_many_objfun_vals` M |
| `main/commanders/simple/simple_commanders_ori.f90` | 183 | `check_states` M, `exec_oriconsensus` M |
| `main/image/simple_image_fft.f90` | 150 | `norm_noise_fft` M, `mask_fft` M, `img2ft` M |
| `main/project/simple_sp_project_core.f90` | 149 | `append_project` M, `set_sp_oris` M, `get_sp_oris` M |
| `main/ori/simple_oris_getters.f90` | 142 | `get_arr` M, `gen_ptcl_mask` M, `get_tseries_neighs` M, `gen_mask` M, `projs_exist` M, `get_all_rmats` M |
| `main/image/simple_image_freq_anal.f90` | 130 | `fcomps_below_noise_power_stats` M, `fsvar` M, `img2spec` M, `frc_pspec` M, `ccf` M |
| `main/ui/simple_ui_program.f90` | 118 | `print_ui` M, `print_param_list` K, `write2json` M, `get_display_name` M, `get_category` M |
| `main/interp/simple_edges_sqwins.f90` | 114 | `cosedge_2`, `cosedge_1`, `cosedge_inner_2`, `cosedge_inner_1`, `sqwin_3d_2`, `sqwin_3d_1`, `sqwin_2d_2`, `sqwin_2d_1`, `hardedge_inner_4`, `hardedge_inner_2`, `hardedge_inner_3`, `hardedge_inner_1` |
| `utils/math/simple_rnd.f90` | 103 | `r8po_fa`, `mnorm_smp`, `irnd_gasdev` K, `rnd_4dim_sphere_pnt`, `randn_2`, `randn_1`, `ran3arr_2`, `ran3arr_1` |
| `fileio/simple_imghead.f90` | 92 | `setdim` M, `getdim` M, `has_ldim_nptcls` |
| `main/commanders/simple/simple_commanders_atoms.f90` | 92 | `exec_map_validate` M, `exec_model_validate_eo` M, `exec_map2model_fsc` M |
| `main/image/simple_image_msk.f90` | 91 | `env_rproject` M, `binary_imat_to_pix`, `automask2d_support_pix` |
| `main/ori/simple_oris_setters.f90` | 91 | `rnd_proj_space` M, `gen_balanced_partitions` M, `zero_projs` M, `e3swapsgn` M |
| `main/strategies/search/simple_corrmat.f90` | 88 | `calc_inpl_invariant_cc_nomirr` |
| `main/image/simple_image_geom.f90` | 87 | `shift2d_demoivre` M, `gauimg3d` M, `corners` M |
| `utils/structs/simple_rec_list.f90` | 86 | `get_projfiles` M, `get_processed_flags` M, `get_busy_flags` M, `push2project_list` M |
| `utils/structs/simple_linked_list.f90` | 82 | `at_real` M, `at_logical` M, `at_int` M, `at_complex` M, `push_front` M, `push_back_real` M, `push_back_logical` M, `push_back_complex` M |
| `main/volume/simple_reconstructor.f90` | 76 | `write_rho_as_mrc` M, `write_absfc_as_mrc` M, `read_rho` M, `set_sh_lim` M |
| `utils/qsys/simple_qsys_funs.f90` | 76 | `qsys_watcher_diag`, `qsys_watcher_2`, `qsys_watcher_1` K |
| `utils/math/simple_ran_tabu.f90` | 74 | `ne_mnomal_iarr` M, `mnomal` M, `irnd_gau` M, `irnd_pair` M |
| `main/image/simple_image_ops.f90` | 73 | `iter_gridcorr` M, `salt_n_pepper` M, `init_gridcorr_mats` M, `remove_neg` M, `zero_neg` M, `inv` M |
| `main/opt/simple_opt_spec.f90` | 73 | `eval_fdf_4` M, `eval_f_8` M, `set_fdfcostfun` M, `eval_f_4` M, `set_inipop` M |
| `main/image/simple_image_access.f90` | 71 | `set_within` M, `get_subimg` M, `set_cmat`, `set_slice` M, `get_slice` M |
| `utils/simple_progress.f90` | 71 | `progressfile_update_part`, `lastfoundfile_update`, `progressfile_complete_parts`, `progressfile_init_part`, `progressfile_init_parts`, `progress_estimate_preprocess_stream` |
| `utils/text/simple_string_utils.f90` | 71 | `removepunct`, `removesp`, `split_str`, `fortran_logical` |
| `main/nano/simple_nanoparticle.f90` | 67 | `center_on_atom` M, `write_individual_atoms` M, `get_valid_corrs` M, `set_ncc` M, `get_img_raw` M, `get_img` M |
| `main/pftc/simple_polarft_geom.f90` | 67 | `mirror_ref_pft_1` M, `rotate_ptcl` M, `rotate_ctf` M, `mirror_ref_pft_2` M |
| `utils/qsys/simple_mem_estimator.f90` | 65 | `estimate_mem_usage_preprocess`, `estimate_mem_usage_pick`, `estimate_mem_usage_motion_corr`, `estimate_mem_usage_extract`, `estimate_mem_usage_ctf_estimation` |
| `main/image/simple_ft_expanded.f90` | 64 | `shift_and_add` M, `corr_unnorm` M, `normalize_mat` M, `get_hp_lp` M, `add_uncond` M |
| `main/image/simple_image_ctf.f90` | 64 | `calc_ice_frac` M |
| `main/image_processing/simple_segmentation.f90` | 64 | `automatic_thresh_sobel` |
| `main/class/simple_classaverager_core.f90` | 61 | `write_cmat_debug` M, `read_cmat` M, `write_ctfsq_debug` M, `write_ctfsq` M, `read_ctfsq` M |
| `main/star/simple_starfile_wrapper.f90` | 61 | `starfile_table__setvalue_bool`, `starfile_table__getvalue_float`, `starfile_table__getvalue_bool`, `starfile_table__setvalue_float`, `starfile_table__append` |
| `main/commanders/simple/simple_commanders_cavgs.f90` | 60 | `exec_map_cavgs_states` M |
| `main/ori/simple_ori.f90` | 59 | `rnd_shift` M, `get_3dshift` M, `write2bild` M, `rnd_inpl` M, `isstate` M, `reset_pparms` M, `kill_hash` M, `kill_chash` M |
| `main/commanders/single/single_commanders_trajectory.f90` | 58 | `exec_trajectory_backgr_subtr` M |
| `utils/simple_is_check_assert.f90` | 58 | `check4nans_2` K, `check4nans_1` K, `check4nans3d_2` K, `check4nans3d_1` K, `check4nans2d_2`, `check4nans2d_1` |
| `main/image/simple_image_seg.f90` | 57 | `taper_edges_hann` M, `zero2one` M, `remove_edge` M, `logical2bin` M |
| `main/image_processing/simple_stackops.f90` | 57 | `make_pcavec_stack` |
| `utils/filter/simple_estimate_ssnr.f90` | 57 | `lowpass_from_klim`, `gaussian_filter`, `ssnr2optlp`, `ssnr2fsc`, `mskdiam2streamresthreshold` |
| `utils/math/simple_histogram.f90` | 56 | `find_hill` M, `kld` M, `find_next_valley` M, `topdf` M, `get_nbins` M |
| `main/ctf/simple_ctf_estimate_fit.f90` | 55 | `write_star` M, `get_parms` M |
| `main/ori/simple_oris_stats.f90` | 55 | `set_extremal_vars` M, `extremal_bound` M, `median_1` M |
| `fileio/simple_fileio.f90` | 54 | `get_open_funits`, `file2imat`, `imat2file` |
| `main/interp/simple_kbinterpol.f90` | 54 | `bessi1`, `bessi0_dp`, `bessi0f` |
| `main/image/simple_image_arith.f90` | 53 | `add_cmats_to_cmat` M, `add_mat2cmat_2` M, `add_mat2cmat_1` M, `div_rmat_at` M, `conjugate` M, `add_cmat_at` M |
| `main/ori/simple_sym.f90` | 49 | `nearest_sym_neighbors` M, `apply2all` M |
| `utils/simple_jiffys.f90` | 48 | `haloween_end` |
| `main/image_processing/simple_micproc.f90` | 45 | `read_mic_subtr_backgr`, `sample_filetab` |
| `main/image_processing/simple_procimgstk.f90` | 45 | `random_cls_from_imgfile` |
| `main/simple_builder.f90` | 45 | `init_params_spproj_tbox2d` M |
| `utils/filter/simple_fsc.f90` | 43 | `plot_fsc2` |
| `main/exec/simple_exec_helpers.f90` | 42 | `async_exec`, `gen_exec_cmd`, `copy_project_file_to_root_dir` |
| `main/image/simple_image.f90` | 39 | `ellipse` M |
| `main/flex/simple_flex_pca_util.f90` | 38 | `project_onto_target_polyline` |
| `main/project/simple_sp_project_mic.f90` | 38 | `get_movies_table` M, `get_mic_kind` M |
| `utils/clustering/simple_kmedoids.f90` | 38 | `merge_ranked` M |
| `main/motion/simple_motion_align_hybrid.f90` | 36 | `align_discrete` M, `set_weights` M, `get_coords` M, `get_corrs` M, `is_fitshifts` M |
| `main/ori/simple_oris_io.f90` | 36 | `write2bild` M |
| `main/ori/simple_oris_reshape.f90` | 36 | `split_state` M, `merge_classes` M |
| `main/strategies/search/simple_srchspace_map.f90` | 36 | `get_inds_in_full` M, `sub2full` M, `full2sub` M, `has_data` M, `get_nsub` M |
| `utils/qsys/simple_qsys_env.f90` | 36 | `get_qsys` M, `exec_simple_prg_in_queue` M |
| `utils/simple_map_reduce.f90` | 36 | `split_pairs_in_parts` |
| `main/project/simple_sp_project_out.f90` | 35 | `get_imgdims_from_osout` M |
| `main/image_processing/simple_imgarr_utils.f90` | 34 | `join_imgarrs`, `write_junk_cavgs` |
| `main/nano/simple_atoms.f90` | 32 | `model_validate_eo` M |
| `main/pftc/simple_polarft_ops_io.f90` | 32 | `write_ptcl_pft_range` M |
| `main/star/simple_starproject_stream.f90` | 32 | `copy_optics` M |
| `main/pftc/simple_polarft_access.f90` | 30 | `get_ptcl_line` M, `ptcl_iseven` M, `get_roind_fast` M, `is_with_ctf` M |
| `main/pick/simple_picksegdiam.f90` | 29 | `write_diameters` M |
| `main/ori/simple_oris_sampling.f90` | 27 | `select_particles_set` M, `set_nonzero_updatecnt` M |
| `main/commanders/simple/simple_commanders_sim.f90` | 26 | `exec_simulate_subtomogram` M |
| `utils/math/simple_math_ft.f90` | 26 | `fplane_get_ctfsq`, `fplane_get_cmplx` |
| `utils/simple_timer_omp.f90` | 23 | `tocprint_omp`, `now_omp`, `tickrate_omp`, `reset_timer_omp` |
| `main/pftc/simple_polarft_core.f90` | 22 | `set_ref_fcomp` M, `set_ptcl_fcomp` M, `cp_odd2even_ref` M, `cp_refs` M |
| `utils/structs/simple_kd_tree.f90` | 21 | `node_height`, `get_kd_tree_height` M, `get_kd_tree_nnodes` M |
| `main/pftc/simple_polarft_shiftmap.f90` | 20 | `shiftmap_peak` M |
| `main/class/simple_class_frcs.f90` | 19 | `print_frcs` M |
| `main/commanders/simple/simple_commanders_project_core.f90` | 19 | `exec_concatenate_projects` M |
| `main/image/simple_image_core.f90` | 19 | `cylinder` M |
| `main/pick/simple_pickref.f90` | 19 | `write_boxfile` M, `get_scores` M |
| `main/commanders/simple/simple_commanders_rec.f90` | 18 | `exec_random_rec` M |
| `main/motion/simple_motion_patched.f90` | 18 | `get_poly_coeffs` M, `set_poly_coeffs` M, `set_nframes` M |
| `main/cavg_quality/simple_cavg_quality_feats.f90` | 16 | `cavg_quality_feature_direction`, `cavg_quality_feature_description` |
| `main/flex/fit/simple_flex_pca_pcg.f90` | 16 | `is_ready` M, `get_npk` M, `get_npairs` M, `get_nexp` M |
| `main/ori/simple_oris_dists.f90` | 16 | `min_euldist` M |
| `main/project/simple_sp_project_cls.f90` | 16 | `get_selected_clsinds` M |
| `utils/math/simple_linalg.f90` | 16 | `projz`, `pythag_sp`, `pythag_dp` |
| `utils/persistent_worker/simple_persistent_worker_server.f90` | 15 | `get_queue_pressure_workers_required` M |
| `fileio/simple_nrtxtfile.f90` | 14 | `writecommentline` M |
| `main/image/simple_image_norm.f90` | 12 | `scale_pixels` M |
| `main/flex/fit/simple_flex_pca_basis.f90` | 11 | `particle_contrast` |
| `main/interp/simple_winfuns.f90` | 11 | `eval_instr` M, `get_whalf` M |
| `main/nano/single_tseries_tracker.f90` | 11 | `tracker_get_nnfname`, `tracker_get_stkname`, `tracker_get_nnn` |
| `main/strategies/search/simple_strategy3D_matcher.f90` | 11 | `print_flags` M |
| `main/image/simple_projector.f90` | 10 | `reset_expanded` M, `is_expanded` M |
| `main/ori/simple_oris_neig.f90` | 10 | `replace_with_closest` M |
| `main/commanders/simple/simple_commanders_solve2D.f90` | 9 | `execute_solve2d_staged` |
| `main/ori/simple_oris_life.f90` | 9 | `kill_chash` M |
| `main/strategies/search/simple_matcher_refvol_utils.f90` | 9 | `reprojection_model_available` |
| `utils/simple_magic_boxes.f90` | 9 | `find_magic_boxes4scale` |
| `main/image/simple_image_vis.f90` | 8 | `print_rmat` M, `print_cmat` M |
| `main/strategies/search/probabilistic/simple_eul_prob_tab2D.f90` | 8 | `fill_tab_prob_snhc` M, `fill_tab` M |
| `main/commanders/simple/simple_commanders_misc.f90` | 7 | `exec_mkdir` M |
| `fileio/simple_imgfile.f90` | 6 | `getdim` M |
| `main/image/simple_ftexp_shsrch.f90` | 6 | `set_factr_pgtol` M |
| `main/volume/simple_dock_vols.f90` | 6 | `set_dock_info` M |
| `main/image/simple_image_checks.f90` | 5 | `same_dims` M |
| `main/motion/simple_motion_align_nano.f90` | 5 | `get_reference` M |
| `main/params/simple_parameters_core.f90` | 5 | `is_final_planned_iter` M |
| `main/pick/simple_segdiam_bin_picker.f90` | 5 | `get_accepted_bins` M |
| `main/preprocess/simple_pspecs.f90` | 5 | `get_nspecs` M |
| `main/project/simple_state_weight_set.f90` | 5 | `get_manifest_path` M |
| `main/volume/simple_reconstructor_pcg.f90` | 5 | `get_solvent_lambda` M |
| `main/solve/simple_solve3D_chunks_plan.f90` | 4 | `get_nsteps` M |
| `main/solve/simple_solve3D_controller.f90` | 4 | `solve3d_cavgs_early_nstages` |
| `main/strategies/search/probabilistic/simple_strategy3D_prob.f90` | 4 | `has_valid_assignment` M |
| `main/volume/simple_wave_projector.f90` | 4 | `get_model` M |
| `main/flex/fit/simple_flex_probe_fit_estep.f90` | 3 | `cov_polar_enabled` |
| `main/image/simple_memoize_ft_maps.f90` | 3 | `ft_map_get_lfny` |

## Appendix B: production code only the tests reach

Without `simple_molecule_data.f90` (section 2) and the vendored libraries.

| file | lines | procedures |
|---|---:|---|
| `main/image_processing/simple_segmentation.f90` | 320 | `detect_peak_thres_fdr`, `canny_edge`, `non_max_supp`, `canny`, `double_thresh`, `sobel`, `detect_peak_thres_for_npeaks` |
| `main/volume/simple_reconstructor_pcg.f90` | 293 | `accumulate_abst2`, `measure_kernel_scale`, `fourier_dot`, `compare_raw_accum`, `build_kernel`, `build_operators`, `set_lambda_relative`, `get_raw_accum`, `build_precond`, `forward_plane`, `adjoint_plane_add`, `set_deapod`, `get_rhs`, `get_ml_prior`, `extract_native_plane`, `get_lims3`, `get_invenv`, `get_env` |
| `main/star/simple_starfile_wrapper.f90` | 177 | `starfile_table__getvalue_string`, `starfile_table__getnames`, `starfile_table__write_omem`, `starfile_table__getcomment`, `starfile_table__read`, `starfile_table__getvalue_int`, `starfile_table__getvalue_double`, `starfile_table__haslabel`, `starfile_table__setcomment`, `starfile_table__hascomment`, `starfile_table__numberofobjects`, `starfile_table__nextobject`, `starfile_table__firstobject` |
| `main/cftc/simple_cartft_calc.f90` | 176 | `pose_normal_terms_full_disk`, `shift_normal_terms_full_disk`, `objective_gradient` |
| `fileio/simple_fileio.f90` | 168 | `files_identical`, `filelength`, `fname2iter`, `file2drarr`, `get_relative_path`, `file2lmat`, `rm_from_fbody`, `lmat2file`, `is_funit_open` |
| `main/interp/simple_kbinterpol.f90` | 154 | `apod_mat_3d_fast_grad`, `apod_kb15_a2`, `apod_fast_value_deriv`, `apod_mat_2d`, `kb_windim` |
| `main/star/simple_starfile.f90` | 147 | `write_ptcl2d_table_parallel` |
| `main/preprocess/simple_particle_extractor.f90` | 135 | `test_movie_outlier_curation`, `display`, `get_particle_mask`, `has_weights`, `has_model_frameweights`, `has_isoshifts`, `has_hotpix_coords`, `has_frames`, `get_l_neg`, `from_model`, `does_exist` |
| `main/image/simple_image_geom.f90` | 110 | `gauimg_1`, `gauimg_2`, `add_window`, `square` |
| `main/nano/simple_atoms.f90` | 110 | `map_validate`, `get_beta`, `get_radius`, `get_occupancy`, `get_name`, `get_element`, `get_atomicnumber` |
| `main/image/simple_projector_pft_batch.f90` | 103 | `fproject_polar_batch_mirr`, `fproject_polar_batch` |
| `main/nano/simple_nano_species.f90` | 88 | `fit_gauss_width`, `enclosed_fraction` |
| `utils/structs/simple_linked_list.f90` | 81 | `pop_front`, `replace_with`, `at_char`, `push_back_char`, `front`, `push_back_int`, `end_iter` |
| `main/pftc/simple_polarft_shiftmap.f90` | 69 | `shiftmap_objfun_grid`, `shiftmap_cross`, `shiftmap_dims`, `shiftmap_cell_area` |
| `utils/math/simple_linalg.f90` | 65 | `plane_from_points`, `eigh_sp`, `vox2ang` |
| `main/motion/simple_motion_model.f90` | 62 | `display` |
| `utils/simple_forked_process.f90` | 62 | `has_exited`, `await_final_status`, `get_exit_code`, `signal_group`, `is_sig_ign` |
| `main/image/simple_image_bin.f90` | 59 | `elim_ccs`, `max_dist` |
| `main/pftc/simple_polarft_core.f90` | 56 | `vol_pad2ref_pfts`, `set_ptcl_pft`, `set_ptcl_ctfmat`, `cp_even_ref2ptcl`, `set_with_ctf` |
| `utils/math/simple_srch_sort_loc.f90` | 52 | `reorder_2`, `reorder_1`, `peakfinder`, `dists2order` |
| `utils/simple_gpu_utils.f90` | 52 | `print_gpu_specs` |
| `main/interp/simple_cartesian_fourier.f90` | 51 | `gather_packed_window_grad`, `center_crop_real3d` |
| `main/pftc/simple_polar_shiftmap.f90` | 50 | `cross_spectrum`, `eval_rotcoeffs`, `peak` |
| `main/project/simple_sp_project_mic.f90` | 49 | `add_single_movie` |
| `utils/qsys/simple_qsys_ctrl.f90` | 46 | `set_jobs_status`, `get_jobs_status`, `failed_parts`, `clear_stack`, `print_jobs_status`, `free_all_cunits` |
| `main/image/simple_image_calc.f90` | 45 | `calc_gradient` |
| `utils/math/simple_math.f90` | 44 | `bounds_from_mask3d`, `gauwfun` |
| `main/pftc/simple_polarft_geom.f90` | 43 | `gen_shmat_8`, `shift_ptcl`, `gen_shmat` |
| `main/volume/simple_wave_projector.f90` | 41 | `get_exit_wave`, `linearize_planes`, `adjoint_planes` |
| `utils/structs/simple_vrefhash.f90` | 40 | `del_str`, `del_char` |
| `main/ui/simple_ui.f90` | 39 | `validate_ui_json`, `count_prgs_in_category`, `make_test_ui` |
| `main/volume/simple_pcg_lattice.f90` | 39 | `scatter_window_nowrap`, `scatter_window_real` |
| `main/class/simple_class_frcs.f90` | 37 | `crop` |
| `utils/structs/simple_chash.f90` | 36 | `move_key_to_front_1`, `reverselookup`, `move_key_to_front_2` |
| `utils/structs/simple_rec_list.f90` | 36 | `get_nptcls_sel_tot`, `subset_to`, `replace_with`, `end_iter` |
| `main/ori/simple_oris_setters.f90` | 34 | `revshsgn`, `revorisgn`, `rnd_lps`, `rnd_corrs` |
| `main/ori/simple_sym.f90` | 32 | `get_subgrp`, `get_pgrp`, `is_platonic`, `is_icosahedral`, `is_dihedral`, `is_circular` |
| `main/sigma2/simple_sigma2_state.f90` | 32 | `sigma2_state_reduce` |
| `main/pftc/simple_polarft_corr.f90` | 31 | `gen_euclid_dists_bank`, `gen_euclid_dists`, `gen_raw_euclid_vals` |
| `main/image/simple_image_ctf.f90` | 30 | `ctf2img` |
| `main/nu_filt/simple_nu_filter_stats.f90` | 30 | `pack_filtmap_lowpass_limits` |
| `utils/text/simple_string.f90` | 28 | `readfile`, `to_dble` |
| `main/ori/simple_oris_transform.f90` | 27 | `rot_1`, `rot_2`, `round_shifts` |
| `main/flex/run/simple_flex_pca_state_service.f90` | 26 | `auto_min_neff`, `auto_box_crop`, `auto_state_count` |
| `main/pftc/simple_pftc_shsrch_grad.f90` | 26 | `shift_prior_penalty`, `shift_prior_penalty_ptcl`, `clear_shift_prior`, `uses_joint_inplane`, `has_shift_prior` |
| `main/pftc/simple_polarft_access.f90` | 25 | `get_ptcl_pft`, `get_ref_pft`, `allocate_ptcl_pft` |
| `main/star/simple_starproject.f90` | 25 | `export_mics` |
| `main/flex/simple_flex_cls_kernel.f90` | 23 | `preimage` |
| `main/ori/simple_oris_sampling.f90` | 23 | `sample4update_updated` |
| `main/strategies/search/probabilistic/simple_eul_prob_tab_utils.f90` | 23 | `insert_table_get_pind`, `insert_table_get_nextra`, `insert_table_get_committed_weight`, `neff_readout_nkept`, `insert_table_get_kmax` |
| `main/image/simple_projector.f90` | 21 | `interp_fcomp_oversamp` |
| `main/ori/simple_oris_dists.f90` | 19 | `find_angres` |
| `main/project/simple_state_weight_set.f90` | 19 | `get_selection`, `get_parent_layout_digest`, `get_parent_generation` |
| `utils/math/simple_rnd.f90` | 19 | `multinomal` |
| `main/opt/simple_opt_spec.f90` | 18 | `set_gcostfun`, `change_opt` |
| `main/image/simple_image_access.f90` | 17 | `get_spat_freq`, `get_lhp`, `is_wthreads` |
| `main/ori/simple_oris_getters.f90` | 17 | `get_update_frac` |
| `main/ori/simple_oris_stats.f90` | 16 | `order_cls` |
| `main/sieve/simple_ptcl_sieve.f90` | 16 | `get_n_pass_2_non_rejected_ptcls`, `get_n_accepted_micrographs` |
| `main/volume/simple_reconstructor.f90` | 16 | `apply_weight`, `get_rho_copy` |
| `main/simple_refine3D_startup.f90` | 13 | `write_cline` |
| `main/image/simple_ftiter.f90` | 12 | `get_spat_freq`, `get_lhp` |
| `main/image/simple_wave_field.f90` | 12 | `crop_t`, `crop` |
| `main/image/simple_image_checks.f90` | 11 | `max_abs_padding` |
| `fileio/simple_imghead.f90` | 9 | `getiform` |
| `main/flex/fit/simple_flex_probe_fit_engine.f90` | 9 | `write_probe_part` |
| `main/sigma2/simple_euclid_sigma2.f90` | 9 | `get_sigma2_part` |
| `main/stream/shared/simple_stream_job_sets.f90` | 9 | `get_job_dir`, `get_counter` |
| `main/image/simple_wave_imaging.f90` | 8 | `get_exit_wave` |
| `main/nano/simple_nano_density.f90` | 8 | `ntab`, `halo_kernel` |
| `utils/gui/metadata/simple_gui_metadata_vol3D.f90` | 8 | `get_oridist` |
| `utils/qsys/simple_qsys_env.f90` | 8 | `n_warmup_envs`, `is_warmup_owner` |
| `main/solve/simple_solve3D_addon_report.f90` | 7 | `get_dock_angle` |
| `main/cftc/simple_cartft_pose_opt.f90` | 6 | `get_motion` |
| `main/flex/fit/simple_flex_pca_pcg.f90` | 6 | `set_lambda_relative` |
| `main/flex/fit/simple_flex_pca_plane_cache.f90` | 6 | `cache_pristine` |
| `main/stream/simple_stream_watcher.f90` | 6 | `clear_history` |
| `main/volume/simple_wave_multislice.f90` | 6 | `get_exit_wave` |
| `utils/gui/metadata/simple_gui_metadata_base.f90` | 6 | `set_assigned` |
| `fileio/simple_binoris.f90` | 5 | `get_first_data_byte` |
| `main/ori/simple_ori.f90` | 5 | `round_shifts` |
| `main/ui/simple_ui_program.f90` | 5 | `get_name` |
| `utils/comm/simple_ipc_tcp_socket_client.f90` | 5 | `set_retry_backoff_ms` |
| `utils/comm/simple_ipc_tcp_socket_helpers.f90` | 5 | `accept_connection` |
| `utils/gui/metadata/stream/simple_gui_metadata_stream_snapshot.f90` | 5 | `get_states` |
| `utils/gui/simple_gui_assembler.f90` | 5 | `is_associated` |
| `defs/simple_refine3D_fnames.f90` | 4 | `refine3d_frozen_manifest_fname` |
| `fileio/simple_imgfile.f90` | 4 | `getiform` |
| `main/class/simple_cavg_sums.f90` | 4 | `get_part` |
| `main/simple_cmdline.f90` | 4 | `get_argcnt` |
| `main/stream/shared/simple_stream_pipe.f90` | 4 | `is_broken` |
| `main/volume/simple_trail_chain_manifest.f90` | 4 | `get_ncontrib` |
| `utils/clustering/simple_hac.f90` | 4 | `get_nclust` |
| `utils/comm/simple_ipc_tcp_socket_server.f90` | 4 | `is_listening` |
| `utils/simple_memory_monitor.f90` | 3 | `mem_monitor_is_enabled` |

## Appendix C: programs with no test, no internal caller and no NICE reference

A heuristic: "test" means the program name in the test code, so a program a test drives through its
commander object is listed although it is exercised.

| program | executable | visibility | doc/script mentions | UI file |
|---|---|---|---:|---|
| `automask2D` | simple_exec | advanced | 6 | `simple_ui_mask.f90` |
| `cluster_cavgs_selection` | simple_exec | advanced | 0 | `simple_ui_cavgproc.f90` |
| `cluster_stack` | simple_exec | advanced | 1 | `simple_ui_cavgproc.f90` |
| `ctf_correct` | simple_exec | advanced | 6 | `simple_ui_image.f90` |
| `ctfops` | simple_exec | advanced | 1 | `simple_ui_image.f90` |
| `detect_calpha` | simple_exec | advanced | 6 | `simple_ui_volume.f90` |
| `fsc_area_score` | simple_exec | advanced | 5 | `simple_ui_res.f90` |
| `map_cavgs_selection` | simple_exec | advanced | 10 | `simple_ui_refine2D.f90` |
| `merge_projects` | simple_exec | advanced | 12 | `simple_ui_project.f90` |
| `noisevol` | simple_exec | advanced | 3 | `simple_ui_solve3D.f90` |
| `nu_filt3D` | simple_exec | advanced | 19 | `simple_ui_filter.f90` |
| `print_dose_weights` | simple_exec | advanced | 0 | `simple_ui_print.f90` |
| `print_fsc` | simple_exec | advanced | 0 | `simple_ui_print.f90` |
| `sample_classes` | simple_exec | advanced | 3 | `simple_ui_refine2D.f90` |
| `select_micrographs` | simple_exec | advanced | 2 | `simple_ui_project.f90` |
| `select_particles` | simple_exec | advanced | 3 | `simple_ui_project.f90` |
| `simulate_noise` | simple_exec | advanced | 1 | `simple_ui_sim.f90` |
| `uniform_filter2D` | simple_exec | advanced | 1 | `simple_ui_filter.f90` |
| `uniform_filter3D` | simple_exec | advanced | 1 | `simple_ui_filter.f90` |
| `vizoris` | simple_exec | advanced | 0 | `simple_ui_ori.f90` |
| `volanalyze` | simple_exec | advanced | 2 | `simple_ui_dock.f90` |
| `volcluster` | simple_exec | advanced | 3 | `simple_ui_dock.f90` |
| `assign_optics_groups` | simple_exec | developer | 6 | `simple_ui_preproc.f90` |
| `consolidate_particles` | simple_exec | developer | 4 | `simple_ui_project.f90` |
| `export_relion` | simple_exec | developer | 9 | `simple_ui_project.f90` |
| `select_clusters` | simple_exec | developer | 0 | `simple_ui_cavgproc.f90` |
| `solve3D_pyramid` | simple_exec | developer | 8 | `simple_ui_solve3D.f90` |
| `auto_spher_mask` | simple_exec | standard | 1 | `simple_ui_mask.f90` |
| `cif2pdb` | simple_exec | standard | 4 | `simple_ui_other.f90` |
| `estimate_lpstages` | simple_exec | standard | 0 | `simple_ui_solve3D.f90` |
| `import_starproject` | simple_exec | standard | 3 | `simple_ui_project.f90` |
| `info_image` | simple_exec | standard | 0 | `simple_ui_print.f90` |
| `info_stktab` | simple_exec | standard | 0 | `simple_ui_print.f90` |
| `make_oris` | simple_exec | standard | 1 | `simple_ui_ori.f90` |
| `match_cavgs` | simple_exec | standard | 4 | `simple_ui_cavgproc.f90` |
| `match_stacks` | simple_exec | standard | 1 | `simple_ui_cavgproc.f90` |
| `model_validate` | simple_exec | standard | 0 | `simple_ui_validation.f90` |
| `oriops` | simple_exec | standard | 0 | `simple_ui_ori.f90` |
| `oristats` | simple_exec | standard | 0 | `simple_ui_ori.f90` |
| `print_magic_boxes` | simple_exec | standard | 0 | `simple_ui_print.f90` |
| `replace_project_field` | simple_exec | standard | 1 | `simple_ui_project.f90` |
| `symmetry_test` | simple_exec | standard | 1 | `simple_ui_sym.f90` |
| `unbootstrap_cavgs` | simple_exec | standard | 2 | `simple_ui_refine2D.f90` |
| `write_mic_filetab` | simple_exec | standard | 0 | `simple_ui_project.f90` |
| `write_classes` | simple_exec | default | 1 | `simple_ui_refine2D.f90` |
| `zero_project_shifts` | simple_exec | default | 0 | `simple_ui_project.f90` |
| `cavgseoproc_nano` | single_exec | developer | 0 | `single_ui_validate.f90` |
| `cavgsproc_nano` | single_exec | developer | 1 | `single_ui_validate.f90` |
| `conv_atom_denoise` | single_exec | developer | 3 | `single_ui_map.f90` |
| `core_atoms_analysis` | single_exec | developer | 0 | `single_ui_atom.f90` |
| `ptclsproc_nano` | single_exec | developer | 1 | `single_ui_validate.f90` |
| `trajectory_make_projavgs` | single_exec | developer | 0 | `single_ui_trajectory.f90` |
| `trajectory_reconstruct3D` | single_exec | developer | 0 | `single_ui_trajectory.f90` |
| `tsegmaps_core_finder` | single_exec | developer | 0 | `single_ui_map.f90` |
| `validate_cavgs_vs_model` | single_exec | developer | 1 | `single_ui_validate.f90` |
| `atoms_register` | single_exec | standard | 0 | `single_ui_atom.f90` |
| `atoms_rmsd` | single_exec | standard | 0 | `single_ui_atom.f90` |
| `crys_score` | single_exec | standard | 0 | `single_ui_atom.f90` |
| `import_tilt_series` | single_exec | standard | 4 | `single_ui_nano3D.f90` |
| `refine_tomo` | single_exec | standard | 9 | `single_ui_nano3D.f90` |
| `trajectory_graphene_subtr` | single_exec | standard | 0 | `single_ui_trajectory.f90` |
| `tseries_import` | single_exec | standard | 0 | `single_ui_tseries.f90` |
| `tseries_make_pickavg` | single_exec | standard | 0 | `single_ui_tseries.f90` |
| `tseries_motion_correct` | single_exec | standard | 0 | `single_ui_tseries.f90` |
| `tseries_track_particles` | single_exec | standard | 1 | `single_ui_tseries.f90` |

## Appendix D: rerunning the analysis

The tools are outside the repository, in `/data2/agent_runs/dead_code_2026-11-02/tools/` on the build box
(`deadcode.py`, about 500 lines of Python 3.12, read-only). `deadcode.py <generated sources>` writes every
node with its class as JSON (77 s for the whole tree); `why.py <name> [P|T]` prints the chain of
references by which the analysis reaches a node; `summarise.py` and `files.py` tabulate. The generated
sources are `lib/simple/simple_args.f90` and `lib/simple/simple_ui_default_values.f90` of any build tree.
The linker check is a Release build with `-ffunction-sections -fdata-sections` and
`-Wl,--gc-sections`, then `nm --defined-only` of the six executables. If the analysis is worth keeping,
it can move to `scripts/` with a short usage note.
