# Release 4: Legacy-Support Inventory

**Status:** Inventory for developer review, 2026-10-03. No code was changed for
this document; it is the input to a cleanup plan once items are ruled on.

**Base:** HEAD `3a93d357a` plus the uncommitted solve/refine2D rename
(`abinitio*` to `solve*`, `cluster2D` to `refine2D`).

Every item has an ID so it can be ruled on. Line numbers are approximate, and
sizes count the lines that would go. Items fall into five groups:

- **A. Dead or unreachable code.** Removing it changes no behaviour, so it can go without a ruling.
- **B. Clean-break removals that change behaviour for old data or old command lines.** Each needs a ruling.
- **C. Decisions about current features that carry legacy baggage.**
- **D. Wording only.** The code uses "legacy" for what is now the default path.
- **E. Checked and kept.** These are current, a required test oracle, or already strict.

Totals: about 4,500 lines in A, about 400 lines in B, and up to about 1,000 lines in C.
The largest single items are the dead NICE channel (A1, about 1,750 lines), unused
stream routines (A2, about 900 lines), the unreachable class-average-quality path
(A3, about 600 lines), and search modes reachable only through undocumented
values (C3/C4, about 900 lines).

---

## A. Dead or unreachable code (remove, no behaviour change)

| ID | What | Where | Size | Notes |
|---|---|---|---|---|
| A1 | Old NICE message channel `simple_nice`. It sends `{"jobid":..,"job":..}` without `version`. `nice_lite/api.py:189` answers every such message with HTTP 400, so it does nothing. | `src/utils/gui/simple_nice.f90`. Still instantiated in the `simple_commanders_` modules `refine2D`, `project_core`, `project_cls`, `project_ptcl`, `starproject` and `imgops`, in `single_commanders_nano2D`, and in the `commanders`/`stream` APIs | ~1,750 + 9 call sites | `simple_gui_communicator` is the live channel. Before deleting, check each commander for reporting that only goes through `simple_nice` (that reporting fails today anyway). |
| A2 | Stream subroutines with no callers. | `simple_stream_utils.f90`: `process_selected_refs`, `process_selected_refs_2`, `wait_for_folder`, `stream_datestr`. `simple_stream_chunk2D_utils.f90`: `analyze2D_new_chunks`, `memoize_chunks`, `update_chunks`. `simple_stream_refine2D_utils.f90`: `update_user_params2D`, which makes `test_repick`, `write_repick_refs` and the repick state dead too, plus `set_dimensions` and `set_resolution_limits`. `p03_initial_analysis`: `micimporter`, `run_cavg_quality_selection`, `run_cavg_size_selection`. Also `progress_estimate_preprocess_stream`. | ~900 | Several of these modules can then come off the `-Wno-unused-function` list (`src/CMakeLists.txt:140-153`). |
| A3 | Retired linear-model k-medoids/Otsu "compatibility cache" in class-average quality. The source says "no longer reachable"; `read_model` already refuses linear models. | `cavg_quality/simple_cavg_quality_model.f90:27-60, 1416-1915`. `simple_cavg_quality_learn.f90`: the `caches` argument (898-916), `classify_training_dataset_cached_detail` (dead after a `return`), and `write_otsu_ablation_diagnostics`, which always reports zero. `normalize_quality_dmat`. | ~600 | Doc fix: `model_cavgs_rejection.md:96-106` wrongly says this path is kept. See C10 for the retired fields still stored in model files. |
| A4 | `cmdline%parse_oldschool` and `print_cmdline_oldschool`. No callers; the help text names `simple_distr_exec`, which no longer exists. | `simple_cmdline.f90:405-454`, `simple_private_prgs.f90:176-217` | ~90 | Enables A5. |
| A5 | Unused entries in the private-program dictionary and table. Of `cmd_dict`'s 289 keys, only 68 are used. Table entries for programs that are now in the main UI are never consulted, and `check_2Dconv`, `check_3Dconv`, `edge_detect` and `track_trajectory` have no dispatch case. | `simple_private_prgs.f90:219-766` | ~320 | The table must be renumbered. |
| A6 | `commander_edge_detect`: never dispatched. | `simple_commanders_imgops.f90:124-238` | ~115 | |
| A7 | `exec_test_mini_stream_legacy`: not bound (`execute => exec_test_mini_stream_quantitative`). | `simple_commanders_test_highlevel.f90:84-225` | ~140 | |
| A8 | Aliases that can never be reached. `simple_exec` and `cmdline%parse` reject non-UI names first. | `report_selection` FIX4NOW (`simple_cmdline.f90:104`, `simple_exec_project.f90:89`); `nonu_filt3D` (`simple_exec_filter.f90:35`); `initial_analysis`/`opening2D` (`simple_stream.f90:71`, "to maintain GUI support"; NICE never launches them) | ~5 | |
| A9 | `mskfile` traps that can never be reached: `mskfile` is not a parameter, so the parser rejects it first. | `simple_parameters_phases.f90:755-757`, `simple_commanders_mask.f90:40` | 4 | |
| A10 | Parameters that are parsed but never read. Certain: `nptcls_per_subcls` (marked "legacy"), `pdfile`, `system`. Probably dead (C9): `cn_type, column_sampling, corr_thres, cs_thres, dfsdev, fraczero, gauref_last_stage, icm_stage, linethres, lp_discrete, nang_nbrs, ncls_sub, ninit, nparts_per_part, nsearch, nspace_max, pcrot, phranlp, protocol, randomise, refine_type, shift_stage, smpd_pickrefs, use_thres, wiener_const`, and the derived flags `l_corrw, l_frac_worst, l_graphene, l_rec_states`. | `params/simple_parameters*.f90` | ~3 lines per key | Removing a key makes old command lines that pass it fail with "argument is not allowed". |
| A11 | GPU euclid benchmark kernels with no callers. | `pftc/simple_polarft_corr.f90:1445-1684`, bindings in `simple_polarft_calc.f90` | ~240 | |
| A12 | `simple_reconstructor_openmpoffload.f90`: `calc_3Drec_gpu` is not `use`d anywhere. | file + `src/CMakeLists.txt:123,133` | ~385 | Move to C if the GPU roadmap still wants it. |
| A13 | Unused constants. | `simple_defs_stream.f90`: `CHUNK_CLS_REJECTED, REJECTED_CLS_STACK, SIEVING_REFS_FNAME, FLUSH_TIMELIMIT, NMICS_DELTA, POOL_FREQ_REJECTION, SIEVING_MATCH_CAVGS_MAX, SIEVING_REF_CAVGS_MAX, STREAM_NMOVS_SET_TIFF, STREAM_SRCHLIM, FRAC_SKIP_REJECTION`. `simple_defs_fname.f90`: `PAUSE_STREAM`, `STREAM_DESELECTED_REFS`, `POLAR_REFS_FBODY` (after B7) | ~15 | |
| A14 | Commented-out blocks and the dead routine behind one of them. | `p03_initial_analysis.f90:484-495`, which also prints the stale `SIMPLE_GEN_PICKREFS NORMAL STOP` at 497 and 856. `p04_refpick_extract_new.f90:347-349, 448-456`, which makes `validate_ptcl2D_star_inputs` (465-505) dead. `p06_pool2D_new.f90:224-232`, `simple_stream_utils.f90:566-593`, `pool2D_utils.f90:1004,1007`, `refine2D_utils.f90:92`. | ~120 | Keep the intent note in `p00_master.f90:1443-1447` if you like it. |
| A15 | Win32 `mq_*` stubs left after message-queue IPC was removed (bc5ed208f). | `src/fileio/posix_stubs_win32.c:30-38` | 9 | |
| A16 | CMake options that do nothing. `USE_MPI` defines `USING_MPI`, which no source uses. `BUILD_DOCS` adds `doc/`, which has no `CMakeLists.txt`, so it would break if turned on. | `CMakeLists.txt:30,33,115-117`, `cmake/Dependencies.cmake:85-94,616-620` | ~20 | |
| A17 | Repo clutter. `CMakeFiles/cmake.check_cache` is tracked in git; `.gitignore` has entries for retired build systems (`Makefile_macros`, `simple_user_input.pm`, `obj/`, `doc/html/`, `docs/`, `/node_modules`, `/dist`, duplicate `*.hed`/`*.spi`); `scripts/CMakeLists.txt:97` excludes `backups/*`, which no longer exists. | | trivial | Local only, not in git: `SIMPLE_TEST_lib_cart_align3D_20261002/`, `fsc{n,t,u}_state01.bin`, `Testing/`. |
| A18 | Dead scripts: `nanokpcadn.py` and `nanokpcadn2.py` (hard-coded personal paths), `remove_git_version_conflict.pl` (it targets calls that no longer exist), `tab2space.pl`, `avg_sym_stats.py`, `relion2emanbox.pl`, `relion2simple.pl` (superseded by `relion2simplegui.pl`), `check_add_ui.py`, `clean_simple_uses.pl`. | `scripts/` | ~600 | |
| A19 | NICE leftovers. `WorkspaceModel.nstr` ("backward-compatibility payload"; never read or written). `Project.new` has a two-signature compatibility form that only tests use. `index.html:138-146` migrates an old session-storage key. `batch_views.py:708` falls back to `latest_micrographs`, a stream key. Unused images: `simple_stream.png`, `simple_stream_min.png`, `single_logo.png`. | `nice/nice_lite/` | ~40 + a migration | `nstr` is best removed together with N2/C12. |

## B. Clean-break removals (behaviour changes for old data or old command lines)

| ID | What it keeps working | Where | Size | Effect of removal |
|---|---|---|---|---|
| B1 | **Retired program names** (`abinitio*`, `cluster2D*`), added today on your instruction. | `simple_ui_legacy_names.f90` and calls in the 5 entry points, `simple_cmdline.f90`, `simple_ui.f90` and the private driver; test `test_legacy_program_names` | ~135 | Old scripts and logged command lines fail with "not recognized". NICE migration 0006 still rewrites stored jobs, independent of this table. This is the most direct conflict with "clean break". |
| B2 | **Project files written before ori slots 41-53 existed.** `read_particle_record` accepts any width from 1 to 53 reals. | `fileio/simple_binoris.f90:666-682`; tests `simple_binoris_tester.f90:177-260` | ~5 + tests | Every `.simple` file with narrower records becomes unreadable. Related: `I_CC_NONPEAK = 42 ! unused` is a dead slot kept to hold the layout; compacting it shifts every later slot index (ask separately). |
| B3 | **sigma2 `parts_import`**: converts pre-canonical `infile<NN>.dat` per-part files into `sigma2_state.bin`. | `simple_commanders_euclid.f90:176-290`, `SIGMA2_PROV_LEGACY_PARTS`, the `sigma_action` choice, `simple_ui_other.f90:50-55` | ~40 | State files with provenance 3 (only ever made by this converter) would fail validation. |
| B4 | **sigma2 state/range layout**: an always-zero particle-checksum header field and 8 bytes per particle of reserved integrity space, "kept for on-disk compatibility". | `fileio/simple_sigma2_state_file.f90` (header 79-80, layout 118-119, writers, range files 604-690); tester `:195-263` | ~35 | Needs `SIMPLE_SIGMA2_V2`; existing state files are refused. Dropping the dead `checksums` argument costs nothing even without the version bump. |
| B5 | **Trail-chain manifest version 1 detection.** It is only used to print "older-format manifest"; the chain is re-seeded either way. | `volume/simple_trail_chain_manifest.f90:7-8,102-116`, `rec_distr.f90:318-321`; tester `:44-49` | ~20 | Version 1 falls into "unreadable", which also re-seeds. |
| B6 | **Reference maps in the pre-2026-08 1/box amplitude convention**, auto-rescaled by a heuristic on every reference read. | `simple_matcher_refvol_utils.f90:347-351, 476-495` (`autorescale_old_convention`); warning text in `simple_euclid_sigma2.f90:442-449` | ~20 | No silent rescale any more. Today it can also misfire on a legitimate low-variance external map. Keep the mis-scaling diagnostic, but reword it. |
| B7 | Cleanup of file names left by older runs: deletes `polar_refs*.bin` and old unpadded `ALGN_FBODY<N>` part documents. | `simple_matcher_refvol_utils.f90:87-90`, `simple_commanders_project_ptcl.f90:519-523` | ~5 | Stale files from old runs stay on disk. |
| B8 | **Flex PCA resume from the old separate deconvolution cache** (`flex_pca_embedding_deconv.bin`, which nothing writes now). | `flex/run/simple_flex_pca_state_service.f90:101-143` | ~45 | A pre-release-4 flex run cannot be resumed. |
| B9 | **`simple_path` stripped from older projects' compenv.** It is harmless: the runtime overrides it anyway. | `simple_sp_project_core.f90:239`, `simple_sp_project_io.f90:1185,1216`; `simple_qsys_env_tester.f90:38-48` | ~15 | An old project keeps its stale path field, which is then ignored. |
| B10 | **`picker=old`**: offered in the UI, then rejected at run time with "Old picker no longer supported". | `simple_ui_params_common.f90:649-651`, `simple_pick_strategy.f90:538`, comment in `simple_parameters.f90:308` | ~5 | The bad choice disappears from the UI. |
| B11 | **`nchunks` alias for `nparts` on `particle_sieving`** ("Backward-compatible alias"). | `simple_ui_preproc.f90:373-375`, `simple_commanders_sieve.f90:29-43` | ~5 | `nchunks=` on the command line stops working. The internal `params%nchunks` stays, because stream `sieve_cavgs` and the engine use it. |
| B12 | **Bench file written twice**: partition 1 also writes the plain per-iteration file "the existing parsers read". | `simple_strategy3D_matcher.f90:256-260`, `simple_refine3D_fnames.f90:268-279` | ~5 | `plot_refine3d_bench.py` already globs `_PART` files (it probably double-counts today). `parse_bench.pl` and `plot_refine3d_bench.py` also still look for `REFINE3D_STAGE_BENCH_ITER*`, which is no longer written. Fix the scripts in the same change. |
| B13 | **File-based stream parameter updates** (`stream_user_params.txt`). Nothing in SIMPLE or NICE writes it any more; NICE pushes updates over HTTP. | `simple_stream_utils.f90:117-210` | ~95 | Only someone hand-writing that file is affected. |

## C. Decisions about current features with legacy baggage

| ID | Question | Where | Size if removed |
|---|---|---|---|
| C1 | **Keep `inpl_cont=no` as a user option?** It is the preserved "legacy in-plane route" and the documented regression reference (`doc/for_developers/preserving_legacy_inplane_path.md`, `solve2D_policy.md:230-255`). `refine=cont` forces it internally, so that internal use stays either way. Either way, rename `new_legacy` in `simple_pftc_shsrch_grad.f90`: it is the production candidate optimizer, not legacy. | ~15 branch sites in strategy2D/3D, `corrmat` and the in-plane tester | 60-100 |
| C2 | **Search modes reachable only through undocumented values, 3D:** `shc_smpl`, `snhc_smpl`, `greedy_inpl`, plus `greedy_smpl` (only the first iteration of the two `*_smpl` modes). No controller sets them, the UI doesn't offer them and no test covers them. Does nano or anyone else use them? | `simple_strategy3D_{shc_smpl,snhc_smpl,greedy_inpl,greedy_smpl}.f90`, matcher dispatch | ~430 |
| C3 | **Search modes reachable only through undocumented values, 2D:** `greedy_smpl`, `inpl_smpl`, `snhc_smpl_many`, and with them the many-reference pftc machinery (`nmany_refs`, `crmat_many`, `plan_bwd_many_refs`, `gen_many_*`). `solve2D_controller.f90:51` still accepts `snhc_smpl_many`. | `simple_strategy2D_*`, `pftc/simple_polarft_{core,calc,corr,access,memo}.f90` | ~520 |
| C4 | **UI bug, which you need to resolve either way:** refine3D (and nano3D) offer `refine=snhc` and `shc_neigh`, but the 3D matcher has no case for them and stops with "unsupported". Remove them from the UI, or implement them. The `shc_neigh` branch in `simple_strategy3D_utils.f90:119` is dead, and `simple_private_prgs.f90:446` advertises stale modes. | `simple_ui_refine3D.f90:85-88`, `single_ui_nano3D.f90:160` | ~5 |
| C5 | **Gold-standard solve3D stage** is switched off with `GOLD_STD_STAGE = TURNED_OFF`: two dead branches, and the demotion to `nonuniform_lpset` always fires. Has gold-standard solve3D been given up? | `simple_solve3D_controller.f90:22,35,411-414,627-630`; `solve3D_policy.md:441` | ~10 |
| C6 | **Trailing bootstrap "legacy previous-halfmap blend"**: the volume-domain blend still seeds the accumulator chain on its first iteration. Replacing it changes numbers. Keep it for release 4? | `simple_commanders_rec_distr.f90:485-525`, `simple_rec3D_pcg_strategy.f90:498` | ~45 (numerical change) |
| C7 | **`euclid_diag=yes`**: diagnostics built for the finished box-division migration, exposed in four UIs. Also the `rec3D_backends` "legacy padded-period instrument function" column (`test_highlevel.f90:2939-2960`). | `simple_euclid_sigma2.f90`, `simple_polarft_corr.f90:1282-1290`, UIs | ~150 |
| C8 | **`pcgop` user key**: production accepts only `kernel`. The matrix-free operator stays as the required test oracle; only the key would go. | `simple_parameters.f90:298`, `simple_rec3D_pcg_strategy.f90:1459` | ~10 |
| C9 | **The probably-dead parameter list in A10.** Some may be experiment knobs someone still wants. | `params/` | ~80 |
| C10 | **Class-average quality model format v11**: drop the retired thresholding fields that `apply_pairwise_logistic` never reads, and regenerate the three built-in presets. The training reader also tolerates a missing `quality_context` and missing overfit columns from older training files, and never checks `cavg_quality_training_version`. | `simple_cavg_quality_model.f90:102-117,~1080-1100`, `simple_cavg_quality_types.f90`, `learn.f90:746,839` | ~80 + presets |
| C11 | **Project repair tolerances**: `map_ptcl_ind2stk_ind` falls back when `nptcls_stk`/`indstk` are missing, and `prune_particles` has a "backwards compatibility" re-read. Every current import writes these fields. `validate_projfile` (~300 lines) is also a general repair tool. Keep the tool and drop only the tolerances? | `simple_sp_project_ptcl.f90:58-110,529-539`, `simple_projfile_utils.f90:1307-1600`; `simple_project_merge_tester.f90:117-161` | 30-330 |
| C12 | **NICE migrations 0001-0006**: squash them into a single fresh `0001_initial`. Existing installs then need a fresh database or `migrate --fake-initial`. Also remove the build-time `makemigrations` (`nice/CMakeLists.txt:104-105`), which can silently create migrations on user installs. | `nice/nice_lite/migrations/` | 6 files to 1 |
| C13 | **Old GUI stats channel** (`.guistats`/`.poolstats`): nothing in SIMPLE or NICE reads `.poolstats`. Keep `generate_pool_jpeg` (`pool_jpeg_map` is still used). | `src/utils/gui/simple_guistats.f90`, `pool2D_utils.f90:974-1028` | ~680 |
| C14 | **Stream parameter shims marked "backwards compatibility"**: `ncls_start`, `nparts_pool`, `nthr2D`, `nparts_chunk`. These are a refactor rather than a deletion; the chunk path still reads them. | `chunk2D_utils.f90:43,47,126`, `pool2D_utils.f90:83,142`, `solve2D_chunks.f90:37,85` | small, but cascades into params and UI |
| C15 | **`solve2D_chunks` together with `simple_stream_chunk.f90`**: an experimental workflow with no tests ("for experimentation" in its first commit). Retire it, or keep it? | `simple_stream_solve2D_chunks.f90`, `simple_stream_chunk.f90` | ~700 |
| C16 | **`*_new` suffixes on stream modules p01/p02/p04/p05/p06.** The old versions are gone, so this is a rename only. | 5 files, ~14 `use` lines, CMake | rename |
| C17 | **`qsys_job_finished`**, a "backward-compatible spelling" of `qsys_declare_part_finished` (33 call sites). Pick one name. | `simple_qsys_funs.f90:69-74` | rename |
| C18 | **UI display-name fallback** ("temporary migration fallback"): 135 of 189 programs have no `display_name`. Make it required (135 descriptor edits), or keep the fallback. | `simple_ui_program.f90:143-145` | ~3 + 135 edits |
| C19 | Smaller items: `extract_substk state<0` "legacy include-all" (`single_commanders_trajectory.f90:435`); `moldiam_max` "kept for API compatibility" (`simple_mini_stream_utils.f90:257,624`); the NU envelope knobs `nu_msk_rel/beta/dens`, exposed although `nonuniform_filtering_policy.md:231-236` allows only `nu_msk_sig` and `amsklp`; `force_volassemble` "legacy handshake" (could become an explicit argument). | | small |
| C20 | **Build and scripts**: `compile_windows.sh` and `compile_gpu_debug.sh` (referenced nowhere); `compile_csbclust.sh` (site-specific); the `compile_gui.sh` TMPDIR workaround; required `LibRt`; personal research scripts (`figs_*`, `chimera_fitmap.py`, ...); the two overlapping `parse_solve3D_metrics*.pl`; `updatesubmissiontemplate.py`; `LD_LIBRARY_PATH` in `wsgi.py.in`. | | varies |

## D. Wording only ("legacy" now names the current default)

Rename these, with no behaviour change:

- solve3D "legacy path / legacy rule" for the ordinary stage ladder (`simple_solve3D_utils.f90:43-47,1035-1038`, `simple_commanders_solve3D.f90:1488`).
- solve2D `fillin=yes` "Legacy" handshake (`simple_solve2D_controller.f90:378`).
- "Legacy FSC/cFAR representation" in the gridding reconstructor (`simple_reconstructor.f90:14-16,316,326`). The comment is stale; also check whether anything still produces the "legacy fplanes" at `:881`. If nothing does, remove that fallback.
- `SHSRCH_LEGACY` / `new_legacy` (C1), and the "legacy score exp(-L)" in `simple_polarft_corr.f90`.
- The "single-fit (legacy)" flex comment (`simple_flex_pca_fit_types.f90:31`); single fits are live.
- `merge_projects` "mixture of canonical and legacy sigma2 chunk projects" (`simple_projfile_utils.f90:125,636`), which really means "has no sigma2 state".
- The `OLD DIRECTORIES` label above `STDERROUT_DIR` (`simple_defs_fname.f90:95`).
- Old-command-line traps on keys that are still valid ("no longer supported" for `nsample_start/stop` in solve3D, 2D nonuniform, `vol1` for pick). Keep the checks and drop "no longer".
- Stale docs: `model_cavgs_rejection.md:7,22,27,96-106` and `distance_transform_shape_rejection_plan.md:10,396-401` cite `simple_cluster2D_rejector`/`simple_microchunked2D`, which no longer exist.

## E. Checked and kept

- **Strict single-version formats with nothing to remove:** the solve3D run manifest, the frozen context and set (schema 2), the PCG raw accumulator, flex weights v2, Cartesian reference volumes, `cavg_sums`, the motion model, the flex caches, and flex probe parts (two current layouts).
- **Required test oracles:** the PCG matrix-free operator and monolithic accumulators, and the `legacy_*` oracles in `simple_cartesian_fourier_tester.f90`.
- **Current features:** `rec_backend=gridding` (the default), `objfun=cc`, `projrec`, `refine=cont`/`pose_cont`, `expand_ft`, `old_root` project remapping, the `msk=` trap (still reachable), `solve3D_cavgs nrestarts_collapse` (used by the stream), and `TERM_STREAM`.
- **The `opening2D` naming in GUI keys and NICE stats** is a coordinated SIMPLE plus NICE rename. Leave it out of release 4 unless you want it.
- **Test asserts that retired test programs stay gone** (`simple_ui_visibility_tester.f90:249-261`).

## Found in passing (bugs, not legacy)

- `production/simple_persistent_worker.f90:289` hard-codes git hash `'342d1484'` instead of `SIMPLE_GIT_VERSION`.
- `scripts/add2.tcshrc.template:7-8`: `setenv SIMPLE_EMAIL="..."` is invalid tcsh syntax.
- `CMakeLists.txt:13` still says `project(SIMPLE VERSION 3.0.0)`.
- All `compile_*.sh` files and `run_fast_gate.sh:11` contain placeholder text ("as in X") in a comment.
