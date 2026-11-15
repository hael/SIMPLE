# Module constant ownership inventory

Status, 2026-11-15: the commander constants were moved in `656a50e55` to `src/defs/simple_defs_center.f90`,
`simple_defs_flex.f90`, `simple_defs_solve3D.f90` and `simple_defs_stream_commanders.f90` (some with a workflow
prefix, e.g. `MASTER_`, `P03_`, `P07_`, `CENTER_`, `CHUNKS_`), and the solve3D_pyramid commander was deleted; the homes
proposed below differ where the code now says otherwise. The rest of the inventory is still a proposal.

Reviewed 2026-11-14 on `master`, source revision `85213734e9eac15a23ef634e4a8bf6cfe3712652`.
The worktree was clean when the inventory began. This is an inventory and proposed ownership map; no Fortran has been changed. Authorship is not inferred from the declarations.

## Findings and scope

The strongest initial cleanup is **48 constants in eight production commander modules**, followed by stream resource/stage policy, solve schedules, FLEX settings and convergence policy. The complete commander module census contains **51** constants; three are private structural details I would retain locally.

The source census found:

| Scope | Named constants | Files |
| --- | ---: | ---: |
| Production commanders, module/submodule scope | 51 | 8 |
| Production commanders, procedure/block scope, counted separately | 161 declaration occurrences | 26 |
| Production commanders, both scopes | 212 declaration occurrences | 31 |
| All SIMPLE-owned production modules/submodules outside `src/defs` | 1,525 distinct file/name pairs | 262 |
| Test modules/submodules, excluded from production candidates | 975 declaration occurrences | 151 |

The broad count is **not a count of misplaced constants**. It includes API enums, binary-format layouts, mathematical coefficients and private implementation constants. Appendix A lists every production module/submodule constant found outside `src/defs`; Appendix B lists all procedure-local commander declarations, including those that should remain local.

Method: scan the 917 Fortran source files under `src/` and `production/`, join continuations, ignore comments, split multiple entities in one declaration, and distinguish module/submodule specifications from procedures, interfaces and derived types. Include intrinsic and derived-type `parameter` declarations, including `double precision`. Exclude `src/extlibs`, existing `src/defs`, test commanders, `simple_test_*` modules and `*_tester` modules from the production census. Arrays count as one named constant. Four conditional-compilation alternatives account for 1,529 source declarations versus 1,525 distinct file/name pairs. Both branches are recorded. This is a lexical source census, not compiler reachability analysis; preprocessor/include-file constants, mutable initialized variables, type-component defaults and unnamed literals are outside the count. Selected repeated literals are discussed where they expose a shared default.

## 1. Production commanders: complete module-scope inventory

Values and spelling below are the current source's. Proposed moves should preserve kinds, expressions and values.

| Source | Constants | Proposed home / decision |
| --- | --- | --- |
| [simple_commanders_center_ptcls.f90:18](../../../src/main/commanders/simple/simple_commanders_center_ptcls.f90#L18) | `CENTER_NSAMPLE=10000`, `CENTER_NCLS=10`, `CENTER_REF_MAXITS=15`, `CENTER_LP=20`, `CENTER_LPSTOP=12`, `CENTER_SHIFT_SIGMA=10` | New `simple_defs_center`: six particle-centering policy settings. Keep these distinct from unrelated equal-valued defaults. |
| [simple_commanders_center_ptcls.f90:24](../../../src/main/commanders/simple/simple_commanders_center_ptcls.f90#L24) | `SAMPLE_PROJ='center_sample.simple'`, `PASS_PROJ='center_pass.simple'`, `PASS_DIR='center_pass'` | Existing `simple_defs_fname`, with centering-specific names. |
| [simple_commanders_flex_pca.f90:10](../../../src/main/commanders/simple/simple_commanders_flex_pca.f90#L10) | `COV_LP_OVER_NYQUIST=2.5`, `COV_LP_DEFAULT=16`, `COV_MINBOX=64`, `FLEX_PCG_MAXITS_DEFAULT=4`, `FLEX_PCG_RTOL_DEFAULT=0` | Existing `simple_defs_flex`. Let the command's default-setting routine and UI consume the same definitions. |
| [simple_commanders_solve3D.f90:38](../../../src/main/commanders/simple/simple_commanders_solve3D.f90#L38) | `CAVGS_SIGMA2_STATE_FNAME='solve3D_cavgs_sigma2_state.bin'` | Existing `simple_defs_fname`. |
| [simple_commanders_solve3D_chunks.f90:21](../../../src/main/commanders/simple/simple_commanders_solve3D_chunks.f90#L21) | `SOURCE_PROJ_NAME='source.simple'`, `BASE_PROJ_NAME='base.simple'` | Existing `simple_defs_fname`, with chunk-workflow prefixes. |
| [simple_commanders_solve3D_pyramid.f90:70](../../../src/main/commanders/simple/simple_commanders_solve3D_pyramid.f90#L70) | `NPTCLS_REFS=100000`, `NCAVGS_STATE_MIN=20`, `STATE_VOL_CC_MIN=0.80` | New `simple_defs_solve3D`: reference-subset size and state acceptance policy. |
| [simple_commanders_solve3D_pyramid.f90:71](../../../src/main/commanders/simple/simple_commanders_solve3D_pyramid.f90#L71) | `REFS_SUBSET_DIR='refs_subset'`, `ASSIGN_DIR='assign'`, `REPORT_FNAME='solve3D_pyramid_report.txt'` | Existing `simple_defs_fname`, with pyramid-specific names. |
| [simple_commanders_solve3D_pyramid.f90:74](../../../src/main/commanders/simple/simple_commanders_solve3D_pyramid.f90#L74) | `NSTEPS_MAX=32` | Retain locally: capacity of this controller's timed-job report. |
| [simple_commanders_stream_p00_master.f90:67](../../../src/main/commanders/stream/simple_commanders_stream_p00_master.f90#L67) | `POOL2D_NCLS=150`, `MASTER_NCUNITS=16`, `MASTER_QSYS_NTHR=16`, `HEARTBEAT_S=5`, `LAST_HEARTBEAT_TRIES=3`, `LAST_HEARTBEAT_WAIT_S=1`, `MEMLOG_EVERY=12`, `OPTICS_STOP_TIMEOUT_S=60`, `STOP_TIMEOUT_S=600`, `STARTUP_STOP_TIMEOUT_S=60`, `LISTENER_SLEEP_US=10000` | Existing `simple_defs_stream`: 11 settings, with master-specific prefixes. Resource resolution stays in `simple_stream_master_resources`. |
| [simple_commanders_stream_p00_master.f90:70](../../../src/main/commanders/stream/simple_commanders_stream_p00_master.f90#L70) | `PICKREFS_3D_KEYS(7)`, `STOP_ONE_BY_ONE(3)` (line 82) | Retain locally: the command's forwarding list and shutdown order. Moving them does not establish a shared policy default and can pull stage-layer dependencies into defs. |
| [simple_commanders_stream_p03_initial_analysis.f90:44](../../../src/main/commanders/stream/simple_commanders_stream_p03_initial_analysis.f90#L44) | `NPARTS2D=8`, `NTHR2D_JOBS=16`, `NCHUNKS=4`, `NSTATES_PICKREFS=3`, `NSTAGES_PICKREFS=3`, `LPSTOP_PICKREFS=8`, `NSPACE_PICKREFS=50`, `NTHR3D_PICKREFS=16`, `NRESTARTS_COLLAPSE=3`, `LPSTART_INI3D=100`, `LPSTOP_INI3D=20` | Existing `simple_defs_stream`: all 11. Rename `NPARTS2D` to describe its actual standalone persistent-worker count. Do not merge it with the master's job-part default of 1. |
| [simple_commanders_stream_p07_solve3D_multistate.f90:45](../../../src/main/commanders/stream/simple_commanders_stream_p07_solve3D_multistate.f90#L45) | `NSTATES3D=3`, `NPARTS3D=8`, `NTHR3D=8` | Existing `simple_defs_stream`; share the appropriate 3D job resource definitions with the master's resource table. |

Total: 48 relocation candidates, three retain-local declarations. This does not include the 161 procedure-local declarations below.

## 2. Commander constants hidden inside procedures

These would be missed by an audit that stops at the module's `contains`.

| Source | Inventory / recommendation |
| --- | --- |
| [simple_commanders_refine3D.f90:512](../../../src/main/commanders/simple/simple_commanders_refine3D.f90#L512) | **32 local declarations overall.** Extract workflow budgets and defaults: `MAXITS_REFINE3D_CAP`, `TARGET_UPDATES_PER_PARTICLE_REFINE3D`, `REGPASS_NSPACE`, the `*_REFINE3D_STATES` overlap/band values, cavgs lineage `STATE_OVERLAP_TARGET`, `LPSTART_DEFAULT`, `LPSTOP_DEFAULT`, `MIN_STATE_FRAC`, `MIN_SUBAVGS`, `FLEX_MIN_NEFF`, `FLEX_MIN_FIT_ROWS`, and the `CLASSIFY3D_REFS` sampling/search/iteration/band block at line 1305. New `simple_defs_refine3D`; namespaced constants for genuinely different workflows. Move `WORK_PROJFILE`, `WORK_SIGMA2_STATE`, `MEAN_VOL_FNAME` to `simple_defs_fname`. Retain diagnostic labels/debug flags; review `SMPD_TOL`, `MINBOX` and `STEP_IT` with their individual owner. |
| [single_commanders_trajectory_solve3D.f90:40](../../../src/main/commanders/single/single_commanders_trajectory_solve3D.f90#L40) | **25 declarations.** Eight output/report names (`RECORDS`, `KEPT`, `DENOISED`, `STARTVOL`, `GATE_TABLE`, `GATE_REJ`, `GATE_KEPT`, `GATE_REPROJ`) belong in `simple_defs_fname`. `NCHUNKS`, `NKEPT_PAIRS`, `LP_EST_DIAM`, `MSK_FAC` belong in new `simple_defs_nano`. Retain the private `NSTAGES`, 11 `S_*` IDs and `STAGE_PRGS` beside the stage sequencing. |
| [single_commanders_nano2D.f90:162](../../../src/main/commanders/single/single_commanders_nano2D.f90#L162) | `STARTVOL='startvol.mrc'` and `LP_EST_DIAM=3` repeat the trajectory workflow's settings. Share the filename and nano policy definitions; include `NCLS_CEN_NANO` at line 34 in the nano default inventory. |
| [single_commanders_tomo.f90:130](../../../src/main/commanders/single/single_commanders_tomo.f90#L130) | Two output names and nine policy/default constants (`LPSTART_DEFAULT`, `LPSTOP_DEFAULT`, `AMSKLP_DEFAULT`, `TRS_DEFAULT`, `ATHRES_DEFAULT`, `MASK_GROW_A`, `NSTAGES_DEFAULT`, `MAXITS_DEFAULT`, `MAXITS_PCG_DEFAULT`) merit `simple_defs_fname` and new `simple_defs_tomo`. The separate `TXT_READ=1` at line 77 copies `simple_nrtxtfile`'s private read-mode value; address that API boundary separately, rather than adding another defs copy. |
| [simple_commanders_validate.f90:32](../../../src/main/commanders/simple/simple_commanders_validate.f90#L32) | `NCLS_MIN=10`, `NCLS_MAX=100`, `LPSTOP=8` are duplicated in two routines (again at line 136), with equivalent initial-analysis-stage defaults. Use explicitly scoped shared refpick-validation settings in `simple_defs_stream`; keep unrelated solve2D defaults distinct. Move its three directory/project names to filename definitions. |
| [simple_commanders_resolest.f90:791](../../../src/main/commanders/simple/simple_commanders_resolest.f90#L791) | `LPSTART_LB`, `LPSTART_DEFAULT`, `LPSTOP_BOUNDS`, `LPSTART_BOUNDS`: resolution-planning policy; check against `simple_defs_solve3D` without assuming the two planners use identical ranges. `TAB_CHAR` and `PLOT_RES` can stay local to formatting/plotting. |
| [simple_commanders_cavgs.f90:1305](../../../src/main/commanders/simple/simple_commanders_cavgs.f90#L1305) | `LP_DEFAULT`, `TRS_DEFAULT`, `NSPACE_DEFAULT`: class-state matching defaults for a scoped definitions module. `NCLUST_MAX=65` repeats at lines 103/1040 and `LP_BIN=20` at lines 859/1368; centralize by owning operation after checking meaning. `STATE_SCORE_FNAME` (line 626) belongs in filename definitions. Five private `QUALITY_MODE_*` dispatch IDs and debug flags can stay. |
| [simple_commanders_solve2D.f90:889](../../../src/main/commanders/simple/simple_commanders_solve2D.f90#L889) | `DIR_PROJS`, `WAITTIME=5`, `TARGET_NCLS_PER_CHUNK=100`: solve2D streaming policy. Use `simple_defs_stream`, but preserve the distinction from its existing `WAITTIME=10`; matching names do not establish equivalent behavior. |
| [simple_commanders_solve3D.f90:924](../../../src/main/commanders/simple/simple_commanders_solve3D.f90#L924) | `ADDON_DIAG_DIR`, `FROZEN_COPY_DIR`, plus `INI3D_DIR` at line 2307: filename definitions. |
| [simple_commanders_flex_pca.f90:145](../../../src/main/commanders/simple/simple_commanders_flex_pca.f90#L145) | `NATIVE_MEAN_FNAME`: filename definitions or the FLEX-specific filename section. |
| [simple_commanders_project_core.f90:413](../../../src/main/commanders/simple/simple_commanders_project_core.f90#L413) | `CONSENSUS_REPORT` and `DEFTAB` (line 773): filename definitions. Keep the private remapping key tables and exact-assignment size bound with their implementation unless reused. |
| [simple_commanders_project_ptcl.f90:805](../../../src/main/commanders/simple/simple_commanders_project_ptcl.f90#L805) | `STKFBODY`: filename definitions. `BUFBYTES`, `REPORT_SECS`, `RATE_WEIGHT`, `MAX_NCUNITS` are secondary I/O/progress/resource-policy candidates; establish an owner before creating a new general-purpose bag of constants. |
| [simple_commanders_volops.f90:831](../../../src/main/commanders/simple/simple_commanders_volops.f90#L831) | `SYMSHTAB` and `SYMAXISTAB`: filename definitions. `MAXBOX=128` appears twice (lines 916/1005); review whether the two symmetry operations should share one cap. `DOCK_HP` and `EDGE` need operation-specific ownership. |
| [simple_commanders_stream_p00_master.f90:312](../../../src/main/commanders/stream/simple_commanders_stream_p00_master.f90#L312) | The two identical local `POLL_US=200000` declarations (also line 424) should use one master polling definition in `simple_defs_stream`. |

Other low-priority local declarations, including simulation batch sizes, student report tables, plotting labels and small debug settings, are listed in Appendix B. This inventory does not recommend moving all of them.

## 3. Larger blocks outside commanders

Counts in this table describe each source block, not an unconditional instruction to move every name.

| Priority | Source / block | Proposed destination and boundary |
| --- | --- | --- |
| High | [simple_stream_master_resources.f90:40](../../../src/main/stream/master/simple_stream_master_resources.f90#L40) — 17 resource defaults | Existing `simple_defs_stream`. Keep the resource type, environment overrides, precedence and logging here. Reconcile the distinct master/standalone names and values before replacing declarations. |
| High | [simple_solve3D_controller.f90:7](../../../src/main/solve/simple_solve3D_controller.f90#L7) — 43 constants | New `simple_defs_solve3D`: stage counts/transitions, iteration and search schedules, filtering and sampling settings. `REC_FBODY` goes to `simple_defs_fname`. The controller keeps its schedule-building logic. |
| High | [simple_solve2D_controller.f90:21](../../../src/main/solve/simple_solve2D_controller.f90#L21) — 16 constants | New `simple_defs_solve2D`: dimensions, stage transitions, sampling and stage-omission policy. Update callers currently importing constants from the controller directly. |
| High | [simple_refine3D_states_controller.f90:23](../../../src/main/simple_refine3D_states_controller.f90#L23) — 10 constants | New `simple_defs_refine3D`, together with the related commander-local defaults; retain distinct per-lineage policies. |
| High | [simple_convergence.f90:13](../../../src/main/simple_convergence.f90#L13) — 2 constants | Existing `simple_defs_conv`: `STAGE_PLATEAU_PATIENCE`, `STAGE_PLATEAU_MI_GAIN`. This is an especially clear existing-home mismatch. |
| High | [simple_state_assignment.f90:79](../../../src/main/solve/simple_state_assignment.f90#L79) — 26 constants | New `simple_defs_solve3D` for state seeding, matching, acceptance and rejection policy; filename definitions for its stage directories and output names. |
| High | [simple_stream_stage_initial_analysis.f90:146](../../../src/main/stream/stages/simple_stream_stage_initial_analysis.f90#L146) — 26 constants | Existing `simple_defs_stream` for 12 settings plus `PICKREFS_SELECTION`; retain 13 private `INIT_*`/`ALL_*` state-machine IDs. |
| High | [simple_stream_stage_pool2D.f90:95](../../../src/main/stream/stages/simple_stream_stage_pool2D.f90#L95) — 11; [simple_stream_stage_preprocess.f90:96](../../../src/main/stream/stages/simple_stream_stage_preprocess.f90#L96) — 8; [simple_stream_stage_refpick.f90:92](../../../src/main/stream/stages/simple_stream_stage_refpick.f90#L92) — 5 | Existing `simple_defs_stream`: publication/iteration schedules, import/presentation limits, gain batches, quiet periods and artifact names. Keep GUI schema tables with the GUI contract where appropriate. |
| High | [simple_stream_stage_solve3D.f90:120](../../../src/main/stream/stages/simple_stream_stage_solve3D.f90#L120) — 20 | Existing `simple_defs_stream` for retry/retention policy and directory/report names; retain the 12 private `PHASE_*`/`JOB_*` IDs. `L_MODEL_LOGGED` is a policy/debug review item, not a required move. |
| High | [simple_ptcl_sieve.f90:65](../../../src/main/sieve/simple_ptcl_sieve.f90#L65) — 25 | New `simple_defs_sieve` for the 11 `DEFAULT_*` settings and sieve marker/artifact names; the sieve also serves non-stream callers. Avoid making the generic sieve depend on the stream implementation. Keep its debug flag and internal labels if private. |
| High | [simple_flex_pca_basis.f90:38](../../../src/main/flex/fit/simple_flex_pca_basis.f90#L38) — 9; [simple_flex_pca_pcg.f90:29](../../../src/main/flex/fit/simple_flex_pca_pcg.f90#L29) — 4; [simple_flex_probe_fit.f90:24](../../../src/main/flex/fit/simple_flex_probe_fit.f90#L24) — 5 | Existing `simple_defs_flex`: rank/calibration budgets and shared fit/PCG policy. Put artifact names in a filename section; retain serialization versions with the owning format unless shared. The basis currently exports policy constants alongside numerical routines. |
| High | [simple_cls_expansion_strategy.f90:30](../../../src/main/strategies/parallelization/simple_cls_expansion_strategy.f90#L30) — 11 | Existing `simple_defs_flex` for three numerical/default settings; `simple_defs_fname` for eight artifact names. Consumers should not import a strategy module solely for `FLEX_CLS_CAVGS_FILE`. |
| Medium | [simple_flex_cls_expansion.f90:16](../../../src/main/flex/simple_flex_cls_expansion.f90#L16) — 9; [simple_flex_cls_kernel.f90:18](../../../src/main/flex/simple_flex_cls_kernel.f90#L18) — 7; [simple_flex_traj_denoise.f90:24](../../../src/main/flex/simple_flex_traj_denoise.f90#L24) — 18 distinct names; [simple_flex_pca_deconv.f90:16](../../../src/main/flex/states/simple_flex_pca_deconv.f90#L16) — 5; [simple_flex_pca_merge.f90:24](../../../src/main/flex/states/simple_flex_pca_merge.f90#L24) — 6; [simple_umap.f90:11](../../../src/main/flex/simple_umap.f90#L11) — 11 | Review cohesive parameter groups for `simple_defs_flex`: CV grids, iteration budgets, regularization and state-merge thresholds. Keep loop tiling, build-dependent BLAS flags and purely private numerical details local. Equal `RIDGE`/`WIENER`/`LWMIN` values need semantic review before unification. |
| High | [simple_nu_filter.f90:44](../../../src/main/nu_filt/simple_nu_filter.f90#L44) — 38 | New lightweight `simple_defs_nu_filter` for bank, objective, label-smoothing, envelope and evidence policy shared by its submodules. The unary cap is also duplicated in the image layer: centralizing that value avoids importing the NU algorithm into `simple_image`. Keep private sentinels/algorithm tags if they have only one owner. |
| Medium | [simple_cavg_quality_feats.f90:53](../../../src/main/cavg_quality/simple_cavg_quality_feats.f90#L53) — 33; [simple_cavg_quality_learn.f90:28](../../../src/main/cavg_quality/simple_cavg_quality_learn.f90#L28) — 14; [simple_cavg_quality_lattice.f90:27](../../../src/main/cavg_quality/simple_cavg_quality_lattice.f90#L27) — 20 | New `simple_defs_cavg_quality` for rejection thresholds, feature scales and learning policy. The existing `simple_cavg_quality_types` already owns shared context/rejection/schema constants: do not create competing definitions. Feature tables depend on `cavg_quality_feature_def`; moving them blindly into a defs module imported by types would create a cycle. |
| Medium | [simple_quality_logistic.f90:22](../../../src/main/quality/simple_quality_logistic.f90#L22) — 13 | If centralizing its learning thresholds/grids, use separate new `simple_defs_quality`, since this fitter is generic. Keep model-size/layout constants with its API if that remains clearer. |
| High, narrow | [simple_picker_utils.f90:14](../../../src/main/pick/simple_picker_utils.f90#L14), [simple_picksegdiam.f90:15](../../../src/main/pick/simple_picksegdiam.f90#L15), [simple_segdiam_bin_picker.f90:33](../../../src/main/pick/simple_segdiam_bin_picker.f90#L33) | Reuse existing `simple_defs` picker settings where the units/meaning match: `PICKER_SHRINK`, `PICKER_SHRINK_REFINE`, `PICKER_OFFSET`. Also review the two procedure-local `SMPD_SHRINK1` copies in `simple_mini_stream_utils`. |
| High, narrow | [simple_nanoparticle.f90:29](../../../src/main/nano/simple_nanoparticle.f90#L29), [simple_nanoparticle_utils.f90:23](../../../src/main/nano/simple_nanoparticle_utils.f90#L23), [simple_calpha_finder.f90:13](../../../src/main/nano/simple_calpha_finder.f90#L13) | Share `NSTRAIN_COMPS=7` in `simple_defs_atoms`; reuse `PI` from `simple_defs` for `PI_LOCAL`, after confirming the required kind/rounding. |
| Medium | [simple_nano_density.f90:30](../../../src/main/nano/simple_nano_density.f90#L30) — 32; [simple_nano_species.f90:20](../../../src/main/nano/simple_nano_species.f90#L20) — 16; [simple_atomic_species.f90:22](../../../src/main/nano/simple_atomic_species.f90#L22) — 12; [simple_nanoparticle.f90:22](../../../src/main/nano/simple_nanoparticle.f90#L22) — 20 | New `simple_defs_nano` for named detection/species/fit policies that need coordinated tuning; filename definitions for stats products. Keep quadrature constants, private table dimensions, schemas and algorithm details with their owners unless genuinely shared. |
| High, narrow | [simple_qsys_local.f90:10](../../../src/utils/qsys/simple_qsys_local.f90#L10) and six other qsys backends | Seven copies of `MAXENVITEMS=100`: new `simple_defs_qsys` or an existing suitably low-level qsys definition owner. Shared capacity is a better target than broad extraction of all backend internals. |
| Medium | [simple_stream_pipe.f90:46](../../../src/main/stream/shared/simple_stream_pipe.f90#L46), [simple_persistent_worker_server.f90:30](../../../src/utils/persistent_worker/simple_persistent_worker_server.f90#L30), [simple_qsys_ctrl.f90:34](../../../src/utils/qsys/simple_qsys_ctrl.f90#L34), [simple_qsys_funs.f90:12](../../../src/utils/qsys/simple_qsys_funs.f90#L12) | Retry, polling, timeout and scaling policy: consider new `simple_defs_ipc`/`simple_defs_qsys`. Keep stream-specific defaults in `simple_defs_stream`; preserve `SHORTTIME=1` versus `SHORTTIME=3` as distinct operations. Keep OS bit flags beside their bindings. |
| Medium | [simple_motion_gain_analysis.f90:12](../../../src/main/motion/simple_motion_gain_analysis.f90#L12) — 8; [simple_motion_align_hybrid.f90:15](../../../src/main/motion/simple_motion_align_hybrid.f90#L15) — 11 | New `simple_defs_motion` if gain convergence and alignment budgets are to be maintained together. Retain purely numerical implementation tolerances where they are not shared. |
| Medium | [simple_halfmap_enhancer.f90:17](../../../src/main/volume/simple_halfmap_enhancer.f90#L17) — 8; [simple_halfmap_diagnostics.f90:18](../../../src/main/volume/simple_halfmap_diagnostics.f90#L18) — 5; [simple_pcg_order_prior.f90:35](../../../src/main/volume/simple_pcg_order_prior.f90#L35) — 8 | Consider separate volume/postprocessing and PCG definitions for policy grids/thresholds. These have correct domain owners already, so they rank below commanders and duplicated shared definitions. |
| Medium | [simple_eul_prob_tab_utils.f90:54](../../../src/main/strategies/search/probabilistic/simple_eul_prob_tab_utils.f90#L54) — 12; [simple_strategy2D_srch.f90:20](../../../src/main/strategies/search/simple_strategy2D_srch.f90#L20) and [simple_strategy3D_srch.f90:15](../../../src/main/strategies/search/simple_strategy3D_srch.f90#L15) — 4 each | Search/prior/readout settings could use new `simple_defs_search`. The four `CONT_ROUTE_*` status values are exactly duplicated between 2D/3D; a shared definition is reasonable if they are intended as one vocabulary. Do not move file magic/version or private numeric floors merely because they are parameters. |
| Medium | [simple_sigma2_state_file.f90:25](../../../src/fileio/simple_sigma2_state_file.f90#L25) — 49; [simple_state_weights_file.f90:22](../../../src/fileio/simple_state_weights_file.f90#L22) — 9; FLEX artifact/cache modules in Appendix A | Move stable workflow filenames where useful; keep most of these large blocks local because they define binary layouts, header slots and versions. A high constant count by itself is misleading here. |

Related small policy blocks worth including when their owner is extracted: `simple_solve3D_chunks_plan` (6), `simple_solve3D_addon_report` (11), `simple_solve3D_utils` (3), `simple_calc_pspec_strategy` (2), `simple_view_partition_sampling` (6), and `simple_project_superset` (2). Appendix A provides the exact symbols and declaration lines.

## 4. Verified duplication and distinctions

| Evidence | Consequence for extraction |
| --- | --- |
| FLEX commander `COV_LP_DEFAULT=16` versus the explicit UI value `16.0` in [simple_ui_heterogeneity.f90:138](../../../src/main/ui/simple/simple_ui_heterogeneity.f90#L138); the commander also documents duplicated PCG defaults | Make UI metadata and execution use the same low-level definitions; check generated UI through the normal generator. |
| `simple_defs_stream::REFPICK_NTHR=4` versus `simple_stream_master_resources::REFPICK_NTHR=8` | They currently have different values. Consolidation must preserve standalone versus master policy or make a separately reviewed behavior decision. Similar resource settings already exist in both modules. |
| p03 `NPARTS2D=8` versus resource-table `INITIAL_ANALYSIS_NPARTS=1` | **Different controls**, despite the misleading p03 name: [simple_commanders_stream_p03_initial_analysis.f90:114](../../../src/main/commanders/stream/simple_commanders_stream_p03_initial_analysis.f90#L114) uses 8 for standalone persistent workers; the master passes 1 as job parts. Keep separate named constants. |
| `NU_OBJECTIVE_UNARY_CAP` at [simple_nu_filter.f90:60](../../../src/main/nu_filt/simple_nu_filter.f90#L60) versus local `NU_LOSS_CAP` at [simple_image_calc.f90:1803](../../../src/main/image/simple_image_calc.f90#L1803) | Same `1./epsilon(1.)` bound, explicitly intended to match. A dependency-light shared definition is a concrete cross-domain extraction. |
| `NSTRAIN_COMPS=7` in nanoparticle and nanoparticle-utils; `MAXENVITEMS=100` in seven qsys backends; `CONT_ROUTE_*` values 0/1/2/3 in both search strategies | Actual duplicate declarations, with shared meaning or a strong reason to consolidate. |
| `LP_EST_DIAM=3` and `STARTVOL='startvol.mrc'` in nano2D and trajectory_solve3D; the trajectory source explicitly says “as analysis2D_nano” | Shared nano workflow settings and filename, with direct source evidence of intended agreement. |
| `SHORTTIME=1` in qsys control versus `SHORTTIME=3` in qsys helpers; `BUFFSZ_DEFAULT=5` in chash versus `10` in hash | Same spelling does not mean same policy. Namespace any extracted definitions; preserve the values. |

## 5. What should stay near its implementation

- Test fixtures, independently chosen tolerances and golden data. The 221 module-scope declarations under commander test/tester files include 42 in `simple_commanders_test_highlevel_msworkflows`, 28 in `simple_polar_matching_bench`, 22 in `simple_test_flex_phantom_fixture`, and 18 golden arrays in `simple_test_wave_projector_golden_fixture`. Move shared test data into test fixture modules when useful, never into production defs merely to reduce the count. Keep test oracle expectations independent of the production constant being tested.
- Private state-machine IDs, local phase/report tables, diagnostic labels and debug flags, unless another component consumes that vocabulary.
- Binary-format magic, versions, header offsets, serialization lengths and algorithm coefficients where the format/algorithm is the sole owner. Examples: `simple_sigma2_state_file`, `simple_eer_factory`, `simple_winfuns`, `simple_rnd` mixing constants and the private L-BFGS-B coefficients.
- Constants that form a small module's own API can remain with that API. Existing shared enum/type ownership, such as `simple_cavg_quality_types`, is not the same problem as policy accumulating in commanders.
- `simple_defs_hhcnn` explicitly states that network `CNN_*` constants and trainer step settings stay private to `simple_hhcnn` and `simple_polar_student`. Follow that existing boundary rather than treating their parameter counts as cleanup defects.
- Repeated UI category counts, plotting palettes, table headers and buffer shapes are low-priority implementation details. They are included in the census for completeness.

## 6. Suggested extraction order and validation

1. Move the 48 commander module candidates, the convergence pair, and the clearly shared duplicate definitions. Use existing `simple_defs_stream`, `simple_defs_flex`, `simple_defs_fname`, `simple_defs_conv`, `simple_defs_atoms` and the existing general picker/math definitions first.
2. Add small `simple_defs_solve2D`, `simple_defs_solve3D` and `simple_defs_refine3D` modules for the related controller and commander-local policy blocks. Add `simple_defs_center`, `simple_defs_nano` and `simple_defs_tomo` only for the identified cohesive families. Keep command-line default application in `simple_default_clines` or the owning normalization routine; constant storage must not introduce dependencies on commanders, parameters or builders.
3. Extract stream-stage/sieve policy and the selected FLEX/NU/quality groups. Review remaining medium-priority domain defaults in their owning subsystem. Avoid creating every proposed module preemptively.
4. Import definitions directly with `only:`; do not re-export them through commanders or create umbrella modules. Prefix generic names (`MAXITS`, `NSPACE`, `LPSTOP`, `WAITTIME`, `ASSIGN_DIR`) so semantically different settings remain distinct. A constants module must not import the algorithm whose constants it supplies. Preserve values, kinds, filename bytes and conditional compilation exactly; do not fold resource-policy changes into mechanical relocation.

For a future implementation, validation is incremental compilation in the existing build, the fast gate, affected area suites and affected high-level entries, plus checks that UI/default generation still agrees. Source/diff checks should show no changed numerical expressions or filenames. Any intended default-value reconciliation requires its own explicit behavioral validation. No build-speed claim is warranted without the policy's like-for-like build measurements.

**Validation of this inventory:** checked declaration coverage against a second source-line search; checked the eight commander blocks and the cited duplicate values/usages; validated local links and source line anchors. No build or runtime test was run because this change is documentation only. No Linux/BOX result is claimed. No code relocation, commit or push was performed.

## Appendix A. Full production module/submodule census outside `src/defs`

This appendix is deliberately broader than the recommended moves. Each name has its declaration line (`@line`). Counts deduplicate conditional alternatives within a file; both source locations are retained. Test and vendor code are excluded. Submodule-local specifications are included.

### src/fileio

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_bench_io.f90:12](../../../src/fileio/simple_bench_io.f90#L12) | 1 | `BENCH_DIR` @12 |
| [simple_binoris.f90:21](../../../src/fileio/simple_binoris.f90#L21) | 5 | `MAX_N_SEGMENTS` @21, `N_VARS_HEAD_SEG` @22, `N_BYTES_HEADER` @23, `THREAD_NSTRINGS` @24, `PTCL_BYTES_PER_REC` @25 |
| [simple_class_sample_io.f90:14](../../../src/fileio/simple_class_sample_io.f90#L14) | 1 | `NHEAD` @14 |
| [simple_eer_factory.f90:14](../../../src/fileio/simple_eer_factory.f90#L14) | 22 | `EER_IMAGE_WIDTH` @14, `EER_IMAGE_HEIGHT` @15, `EER_IMAGE_PIXELS` @16, `TIFF_COMPRESSION_EER8bit` @17, `TIFF_COMPRESSION_EER7bit` @18, `EER_BYTE_MASK` @19, `EER_XMASK` @20, `EER_YSHIFT` @21, `EER_7BIT_SKIP` @22, `EER_SYMBOL_XOR` @23, `EER_8K_XMASK` @24, `EER_8K_YMASK` @25, `EER_8K_SCALE` @26, `EER_8K_XSHIFT` @27, `EER_8K_YSHIFT` @28, `EER_16K_XMASK` @29, `EER_16K_YMASK` @30, `EER_16K_SCALE` @31, `EER_16K_YSHIFT` @32, `one2` @33, `EER_7BIT_MASK8` @34, `EER_SYMBOL_MASK8` @35 |
| [simple_image_slab.f90:17](../../../src/fileio/simple_image_slab.f90#L17) | 3 | `SLAB_OK` @17, `SLAB_OVER_BUDGET` @17, `SLAB_NO_MEMORY` @17 |
| [simple_imgfile.f90:26](../../../src/fileio/simple_imgfile.f90#L26) | 7 | `TIFF_SAMPLEFORMAT_UINT` @26, `TIFF_SAMPLEFORMAT_INT` @27, `TIFF_SAMPLEFORMAT_IEEEFP` @28, `TIFF_SAMPLEFORMAT_VOID` @29, `TIFF_SAMPLEFORMAT_COMPLEXINT` @30, `TIFF_SAMPLEFORMAT_COMPLEXIEEEFP` @31, `FLOAT16_WRITE_BUFFER_ELEMS` @32 |
| [simple_imghead.f90:36](../../../src/fileio/simple_imghead.f90#L36) | 9 | `NREALS` @36, `NREMAINS` @37, `CLOSE2THEANSWER` @38, `MRCHEADSZ` @39, `NLABL` @40, `MRC_MODE_FLOAT32` @41, `MRC_MODE_COMPLEX_FLOAT32` @42, `MRC_MODE_FLOAT16` @43, `MRC_NVERSION_20141` @44 |
| [simple_nrtxtfile.f90:15](../../../src/fileio/simple_nrtxtfile.f90#L15) | 2 | `OPEN_TO_READ` @15, `OPEN_TO_WRITE` @16 |
| [simple_sigma2_state_file.f90:25](../../../src/fileio/simple_sigma2_state_file.f90#L25) | 49 | `SIGMA2_STATE_FNAME` @25, `SIGMA2_STATE_NEXT_FNAME` @26, `STATE_MAGIC` @27, `RANGE_MAGIC` @28, `SIGMA2_STATE_VERSION` @34, `SIGMA2_GROUP_GLOBAL` @35, `SIGMA2_PROV_PSPEC` @36, `SIGMA2_PROV_RESIDUAL` @37, `SIGMA2_RANGE_LOCAL` @38, `SIGMA2_STATE_CANDIDATE` @39, `SIGMA2_STATE_COMMITTED` @40, `STATE_HEADER_BYTES` @42, `RANGE_HEADER_BYTES` @43, `STATE_NWORDS` @44, `RANGE_NWORDS` @45, `FNV_OFFSET` @46, `FNV_PRIME` @47, `W_VERSION` @49, `W_HEADER_BYTES` @49, `W_REAL_BYTES` @49, `W_KFROM` @49, `W_KTO` @49, `W_NPTCLS` @50, `W_BOX` @50, `W_NGROUPS` @50, `W_GROUPING` @50, `W_GENERATION` @50, `W_PROVENANCE` @51, `W_STATE` @51, `W_RANGE_IO` @51, `W_LAYOUT_DIGEST` @51, `W_GLOBAL_OFFSET` @52, `W_PARTICLE_OFFSET` @52, `W_FILE_BYTES` @52, `W_RETIRED_CHECKSUM` @55, `W_SMPD_BITS` @55, `W_RETIRED_HEADER_CHECKSUM` @55, `RW_VERSION` @57, `RW_HEADER_BYTES` @57, `RW_REAL_BYTES` @57, `RW_KFROM` @57, `RW_KTO` @57, `RW_FIRST` @58, `RW_LAST` @58, `RW_GENERATION` @58, `RW_LAYOUT_DIGEST` @58, `RW_DATA_OFFSET` @59, `RW_FILE_BYTES` @59, `RW_RETIRED_CHECKSUM` @60 |
| [simple_stack_read_session.f90:25](../../../src/fileio/simple_stack_read_session.f90#L25) | 11 | `SPAN_GAP_BUDGET` @25, `SPAN_SCRATCH_CAP` @26, `PREFETCH_WINDOW_BYTES` @27, `PREFETCH_MAX_STACKS` @28, `ST_OK` @30, `ST_OPEN` @30, `ST_HEADER` @30, `ST_DIMS` @30, `ST_INDEX` @30, `ST_SHORT` @30, `ST_READ` @30 |
| [simple_state_weights_file.f90:22](../../../src/fileio/simple_state_weights_file.f90#L22) | 9 | `STATE_WEIGHTS_MANIFEST_FNAME` @22, `STATE_WEIGHTS_FBODY` @23, `FILE_MAGIC` @24, `MANIFEST_FORMAT` @25, `MANIFEST_VERSION` @26, `NHEADER_WORDS` @27, `HEADER_BYTES` @28, `STATE_WEIGHTS_KIND_PARTITION` @30, `STATE_WEIGHTS_KIND_KERNEL` @32 |
| [simple_syslib.f90:354](../../../src/fileio/simple_syslib.f90#L354) | 1 | `MAX_IO_DEVS` @354 |

### src/main/cavg_quality

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_cavg_quality_feats.f90:53](../../../src/main/cavg_quality/simple_cavg_quality_feats.f90#L53) | 33 | `LOG_EPS` @53, `FOREGROUND_SEG_LP` @56, `SIGNAL_METRIC_LP` @57, `OVERFIT_SIGNAL_BP_HP` @58, `OVERFIT_SIGNAL_BP_LP` @59, `CAVG_RES_HARD_REJECT_A` @60, `POP_FRACTION_HARD_REJECT` @61, `BP_CENTER_EDGE_VAR_HARD_REJECT_MIN` @62, `FUZZY_BALL_SIGNAL_HARD_REJECT_MIN` @63, `SIEVE_BP_CENTER_EDGE_VAR_HARD_REJECT_MIN` @64, `CHUNK_LOCVAR_FG_HARD_REJECT_MAX` @65, `POOL_RES_HARD_REJECT_A` @66, `POOL_POP_FRACTION_HARD_REJECT` @67, `POOL_BP_CENTER_EDGE_VAR_HARD_REJECT_MIN` @68, `LOCVAR_WINDOW` @69, `MASK_HARD_OUTSIDE_PIXELS` @70, `ANALYSIS_BOXSIZE` @71, `I_LOG_POP` @73, `I_NEG_LOG_RES` @74, `I_CENTERED` @75, `I_LOCVAR_FG` @76, `I_LOCVAR_BG` @77, `I_CORR_FRC` @78, `I_CENTER_EDGE_SNR` @79, `I_CC_AREA_FRAC` @80, `I_PRESENCE` @81, `I_NEG_LOCVAR_FG` @82, `I_NEG_LOCVAR_BG` @83, `I_LOG_LOCVAR_FG_BG_RATIO` @84, `I_BP40_100_CENTER_EDGE_VAR` @85, `I_FUZZY_BALL_SIGNAL` @86, `FEATURE_DEFS` @89, `FEATURE_DEFS_NANO` @120 |
| [simple_cavg_quality_lattice.f90:27](../../../src/main/cavg_quality/simple_cavg_quality_lattice.f90#L27) | 20 | `LAT_NFEATS` @27, `LAT_LATTICE_PEAK` @28, `LAT_SPOT_COUNT` @28, `LAT_SPOT_ARC` @28, `LAT_SPOT_DIAM` @28, `LAT_SPOT_SNR` @28, `LAT_BAND_PRESENCE` @29, `NLOC` @30, `SPOT_HALFWIN` @31, `NWIN` @32, `BAND_LOW_RES` @33, `LAT_PEAK_MIN` @40, `LAT_PEAK_MIN_RAW` @46, `K_MODEL` @54, `ND` @55, `NOFF` @56, `D0` @57, `DSTEP` @57, `OFF0` @57, `OFFSTEP` @57 |
| [simple_cavg_quality_learn.f90:28](../../../src/main/cavg_quality/simple_cavg_quality_learn.f90#L28) | 14 | `CAVG_QUALITY_LEARN_N_STANDARD_POLICIES` @28, `CAVG_QUALITY_LEARN_N_NANO_POLICIES` @29, `LEARN_RECALL_ONLY_FLOOR` @30, `LEARN_RECALL_ONLY_PENALTY` @31, `LEARN_BALANCED_FN_TOLERANCE_FRAC` @32, `LEARN_BALANCED_FN_RATE_PENALTY` @33, `LEARN_BAD_OVERFIT_LOSS_MULT` @34, `LEARN_OVERFIT_FOCUS_BAD_LOSS_SCALE` @35, `LEARN_BAD_OVERFIT_LOWVAR_FG_MIN` @36, `LEARN_BAD_OVERFIT_SUPPORT_MAX` @37, `LEARN_OVERFIT_FP_FOCUS_MIN_TRAINABLE_FRAC` @38, `LEARN_OVERFIT_FP_ACCEPT_RATE_TARGET` @39, `LEARN_OVERFIT_FP_EXCESS_RATE_PENALTY` @40, `LEARN_OVERFIT_FP_BP_Z_MAX` @41 |
| [simple_cavg_quality_model.f90:32](../../../src/main/cavg_quality/simple_cavg_quality_model.f90#L32) | 9 | `CAVG_QUALITY_MODEL_CHUNK_DEFAULT` @32, `CAVG_QUALITY_MODEL_SIEVE_DEFAULT` @33, `CAVG_QUALITY_MODEL_POOL_DEFAULT` @34, `CAVG_QUALITY_MODEL_NANO_DEFAULT` @35, `CAVG_QUALITY_BUILTIN_MODELS` @36, `BUILTIN_MODEL_NAMES` @39, `NANO_LATTICE_PEAK_SCALE` @47, `CHUNK100MICS_FEATURE_POLICY` @53, `CAVG_QUALITY_LOGISTIC_WEIGHTS` @54 |
| [simple_cavg_quality_relations.f90:23](../../../src/main/cavg_quality/simple_cavg_quality_relations.f90#L23) | 1 | `CAVG_RELATIONAL_FEATURE_NAME` @23 |
| [simple_cavg_quality_types.f90:41](../../../src/main/cavg_quality/simple_cavg_quality_types.f90#L41) | 25 | `CAVG_QUALITY_NFEATS` @41, `CAVG_QUALITY_MAX_INTERACTIONS` @42, `CAVG_QUALITY_TRAINING_VERSION` @43, `EPS` @44, `CLIP_Z` @45, `CAVG_RELATIONAL_SCHEMA_NONE` @46, `CAVG_RELATIONAL_SCHEMA_CORR_KNN_SIGNAL_V1` @47, `CAVG_RELATIONAL_DEFAULT_KNN` @48, `CAVG_RELATIONAL_DEFAULT_CORR_HP` @49, `CAVG_RELATIONAL_DEFAULT_CORR_LP` @50, `CAVG_RELATIONAL_DEFAULT_CORR_TRS` @51, `CAVG_RELATIONAL_NANO_CORR_HP` @53, `CAVG_RELATIONAL_NANO_CORR_LP` @54, `CAVG_RELATIONAL_NANO_CORR_TRS` @55, `CAVG_QUALITY_CONTEXT_CHUNK` @63, `CAVG_QUALITY_CONTEXT_SIEVE` @64, `CAVG_QUALITY_CONTEXT_POOL` @65, `CAVG_QUALITY_CONTEXT_NANO` @66, `CAVG_QUALITY_CONTEXT_HARD_GATES` @67, `CAVG_REJECT_REASON_NONE` @69, `CAVG_REJECT_REASON_POP` @70, `CAVG_REJECT_REASON_BAD_PIXELS` @71, `CAVG_REJECT_REASON_NO_COMPONENT` @72, `CAVG_REJECT_REASON_MASK_GEOMETRY` @73, `CAVG_REJECT_REASON_MODEL` @74 |

### src/main/cftc

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_cartft_calc.f90:24](../../../src/main/cftc/simple_cartft_calc.f90#L24) | 5 | `NUMERIC_FLOOR` @24, `CURVATURE_SCALE_MIN` @26, `CURVATURE_SCALE_MAX` @26, `CART_REFVOLS_FORMAT_VERSION` @27, `CART_REFVOLS_NHEADER` @28 |
| [simple_cartft_pose_opt.f90:24](../../../src/main/cftc/simple_cartft_pose_opt.f90#L24) | 24 | `CARTFT_NOT_ATTEMPTED` @24, `CARTFT_INVALID_PREPARATION` @25, `CARTFT_ACCEPTED` @26, `CARTFT_NO_IMPROVEMENT` @27, `CARTFT_NO_RELIABLE_UPDATE` @28, `CARTFT_BOUND_REJECTED` @29, `CARTFT_INVALID_NUMERICS` @30, `CARTFT_ITERATION_LIMIT` @31, `CARTFT_STAGE_SHIFT` @32, `CARTFT_STAGE_JOINT` @32, `NUMERIC_FLOOR` @35, `SHIFT_RESIDUAL_FLOOR` @36, `LM_INITIAL_DAMPING` @37, `LM_INITIAL_REJECTION_MULT` @38, `LM_MAX_DAMPING` @39, `LM_ACCEPTED_RELATIVE_TOL` @40, `LM_PREDICTED_RELATIVE_TOL` @41, `LM_MAX_CONSECUTIVE_REJECTIONS` @42, `DEFAULT_MAXITS` @44, `DEFAULT_ROTATION_SCALE` @45, `DEFAULT_SHIFT_STEP_NATIVE` @46, `NRATIO_BINS` @48, `RATIO_EDGES` @49, `RATIO_LABELS` @50 |

### src/main/class

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_cavg_sums.f90:18](../../../src/main/class/simple_cavg_sums.f90#L18) | 9 | `CAVG_SUMS_STATE` @18, `CAVG_SUMS_CONTRIB` @19, `CAVG_SUMS_OK` @20, `CAVG_SUMS_MISSING` @21, `CAVG_SUMS_CORRUPT` @22, `CAVG_SUMS_MISMATCH` @23, `FORMAT_VERSION` @24, `MAGIC` @25, `TMP_SUFFIX` @26 |
| [simple_cavgs_state_match.f90:29](../../../src/main/class/simple_cavgs_state_match.f90#L29) | 3 | `LP_SAMPLING` @29, `HP_MATCH` @30, `BOX_MATCH_MIN` @31 |
| [simple_class_compatibility.f90:28](../../../src/main/class/simple_class_compatibility.f90#L28) | 14 | `PREPROCESS_BOXSIZE` @28, `PREPROCESS_MORPH_SIZE` @29, `NRELAX` @32, `NQ` @33, `DEFAULT_RELAX` @34, `DEFAULT_QLO` @35, `DEFAULT_QHI` @36, `CONV_ABS_EPS` @37, `CONV_REL_EPS` @38, `RESCUE_EDGE_FRAC` @39, `SUPPORT_EDGE_SOFT_FRAC_MAX` @40, `RELAX_GRID` @41, `QLOW_GRID` @42, `QHIGH_GRID` @43 |
| [simple_classaverager.f90:107](../../../src/main/class/simple_classaverager.f90#L107) | 2 | `READBUFFSZ` @107, `CAVG_KB_WDIM` @108 |

### src/main/commanders

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_commanders_center_ptcls.f90:18](../../../src/main/commanders/simple/simple_commanders_center_ptcls.f90#L18) | 9 | `CENTER_NSAMPLE` @18, `CENTER_NCLS` @19, `CENTER_REF_MAXITS` @20, `CENTER_LP` @21, `CENTER_LPSTOP` @22, `CENTER_SHIFT_SIGMA` @23, `SAMPLE_PROJ` @24, `PASS_PROJ` @25, `PASS_DIR` @26 |
| [simple_commanders_flex_pca.f90:10](../../../src/main/commanders/simple/simple_commanders_flex_pca.f90#L10) | 5 | `COV_LP_OVER_NYQUIST` @10, `COV_LP_DEFAULT` @11, `COV_MINBOX` @13, `FLEX_PCG_MAXITS_DEFAULT` @18, `FLEX_PCG_RTOL_DEFAULT` @19 |
| [simple_commanders_solve3D.f90:38](../../../src/main/commanders/simple/simple_commanders_solve3D.f90#L38) | 1 | `CAVGS_SIGMA2_STATE_FNAME` @38 |
| [simple_commanders_solve3D_chunks.f90:21](../../../src/main/commanders/simple/simple_commanders_solve3D_chunks.f90#L21) | 2 | `SOURCE_PROJ_NAME` @21, `BASE_PROJ_NAME` @22 |
| [simple_commanders_solve3D_pyramid.f90:70](../../../src/main/commanders/simple/simple_commanders_solve3D_pyramid.f90#L70) | 7 | `NPTCLS_REFS` @70, `REFS_SUBSET_DIR` @71, `ASSIGN_DIR` @72, `REPORT_FNAME` @73, `NSTEPS_MAX` @74, `NCAVGS_STATE_MIN` @76, `STATE_VOL_CC_MIN` @77 |
| [simple_commanders_stream_p00_master.f90:67](../../../src/main/commanders/stream/simple_commanders_stream_p00_master.f90#L67) | 13 | `POOL2D_NCLS` @67, `PICKREFS_3D_KEYS` @70, `MASTER_NCUNITS` @73, `MASTER_QSYS_NTHR` @74, `HEARTBEAT_S` @75, `LAST_HEARTBEAT_TRIES` @76, `LAST_HEARTBEAT_WAIT_S` @77, `MEMLOG_EVERY` @78, `OPTICS_STOP_TIMEOUT_S` @79, `STOP_TIMEOUT_S` @80, `STOP_ONE_BY_ONE` @82, `STARTUP_STOP_TIMEOUT_S` @83, `LISTENER_SLEEP_US` @84 |
| [simple_commanders_stream_p03_initial_analysis.f90:44](../../../src/main/commanders/stream/simple_commanders_stream_p03_initial_analysis.f90#L44) | 11 | `NPARTS2D` @44, `NTHR2D_JOBS` @46, `NCHUNKS` @47, `NSTATES_PICKREFS` @49, `NSTAGES_PICKREFS` @50, `LPSTOP_PICKREFS` @51, `NSPACE_PICKREFS` @52, `NTHR3D_PICKREFS` @53, `NRESTARTS_COLLAPSE` @54, `LPSTART_INI3D` @55, `LPSTOP_INI3D` @56 |
| [simple_commanders_stream_p07_solve3D_multistate.f90:45](../../../src/main/commanders/stream/simple_commanders_stream_p07_solve3D_multistate.f90#L45) | 3 | `NSTATES3D` @45, `NPARTS3D` @46, `NTHR3D` @47 |

### src/main/ctf

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_ctf_estimate_fit.f90:17](../../../src/main/ctf/simple_ctf_estimate_fit.f90#L17) | 6 | `TOL` @17, `NSTEPS` @18, `NBINS` @18, `POLYDIM` @18, `DEBUG_HERE` @19, `BENCH` @20 |

### src/main/flex

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_flex_pca_basis.f90:38](../../../src/main/flex/fit/simple_flex_pca_basis.f90#L38) | 9 | `COV_EIG_REL_FLOOR` @38, `COV_CALIB_MAX_PTCLS` @40, `COV_MAX_DTILDE` @42, `COV_DEFAULT_DTILDE` @45, `COV_SAMPLES_PER_PARAM` @46, `COV_UNIT_CONTRAST` @47, `COV_MEAN_FROM_DATA` @48, `COV_PROBE_META` @49, `MEAN_SCALE_FNAME` @51 |
| [simple_flex_pca_crossfsc.f90:17](../../../src/main/flex/fit/simple_flex_pca_crossfsc.f90#L17) | 4 | `COV_XFSC_MAGIC` @17, `COV_XFSC_VERSION` @18, `COV_XFSC_FNAME` @19, `COV_XFSC_H_FLOOR` @21 |
| [simple_flex_pca_embed.f90:29](../../../src/main/flex/fit/simple_flex_pca_embed.f90#L29) | 3 | `COV_EMBED_CONTRAST_GRID` @29, `GRAM_DIAG_STRIDE` @30, `EMBED_STATS_VERSION` @31 |
| [simple_flex_pca_pcg.f90:29](../../../src/main/flex/fit/simple_flex_pca_pcg.f90#L29) | 4 | `FLEX_PCG_KSCALE` @29, `FLEX_PCG_RIDGE_REL` @31, `FLEX_PCG_LAMBDA_REL_DEFAULT` @33, `FLEX_PCG_SUPPORT_FLOOR` @36 |
| [simple_flex_pca_plane_cache.f90:27](../../../src/main/flex/fit/simple_flex_pca_plane_cache.f90#L27) | 13 | `PLANE_CACHE_MAGIC` @27, `PLANE_CACHE_VERSION` @28, `NHEAD` @29, `IH_NROWS` @30, `IH_SMPD` @30, `IH_NSEL` @30, `IH_HIGHROW` @30, `CACHE_DIR_ENV` @31, `FREE_SPACE_DENOM` @32, `HASH_M1` @34, `HASH_A1` @34, `HASH_M2` @35, `HASH_A2` @35 |
| [simple_flex_pca_polar.f90:18](../../../src/main/flex/fit/simple_flex_pca_polar.f90#L18) | 1 | `POLAR_NANG_MIN` @18 |
| [simple_flex_pca_posterior.f90:23](../../../src/main/flex/fit/simple_flex_pca_posterior.f90#L23) | 1 | `COV_PINV_RCOND` @23 |
| [simple_flex_probe_fit.f90:24](../../../src/main/flex/fit/simple_flex_probe_fit.f90#L24) | 5 | `COV_PROBE_CONV` @24, `COV_EM_MIX` @27, `COV_EM_DEFLATE` @30, `PROBE_PART_VERSION` @31, `MIX_ZSUB_MAX` @32 |
| [simple_flex_reconstructor_latent_ops.f90:26](../../../src/main/flex/fit/simple_flex_reconstructor_latent_ops.f90#L26) | 3 | `LATENT_WDIM` @26, `COUPLED_MSTEP_RIDGE_REL` @28, `COUPLED_DENSITY_FLOOR` @29 |
| [simple_flex_pca_application.f90:45](../../../src/main/flex/run/simple_flex_pca_application.f90#L45) | 2 | `MIN_NSTATES` @45, `TRIAL_PREFIX` @47 |
| [simple_flex_pca_artifacts.f90:20](../../../src/main/flex/run/simple_flex_pca_artifacts.f90#L20) | 5 | `FLEX_PCA_PART_MAGIC` @20, `MERGED_PC_FBODY` @23, `MERGED_META` @24, `MERGED_EIG_FNAME` @25, `PAIRED_MANIFEST` @26 |
| [simple_flex_pca_embedding_io.f90:20](../../../src/main/flex/run/simple_flex_pca_embedding_io.f90#L20) | 4 | `COV_CACHE_MAGIC` @20, `COV_CACHE_VERSION` @22, `COV_TRAIL_MAGIC` @25, `PART_SUFFIX` @26 |
| [simple_flex_pca_project_gateway.f90:27](../../../src/main/flex/run/simple_flex_pca_project_gateway.f90#L27) | 1 | `SIGMA_STATE_FNAME` @27 |
| [simple_flex_pca_stages.f90:20](../../../src/main/flex/run/simple_flex_pca_stages.f90#L20) | 7 | `PCA_STAGE_EMBED` @20, `PCA_STAGE_PROBE` @23, `PCA_STAGE_FULLSET` @27, `FLEX_FIT_ALL` @31, `FLEX_FIT_A` @32, `FLEX_FIT_B` @33, `FLEX_MOD4_PAIRING` @36 |
| [simple_flex_pca_state_service.f90:32](../../../src/main/flex/run/simple_flex_pca_state_service.f90#L32) | 3 | `PARTITION_ROW_TOL` @32, `AUTO_NSTATES` @38, `POP_FLOOR_MAX_NSTATES` @40 |
| [simple_flex_cls_expansion.f90:16](../../../src/main/flex/simple_flex_cls_expansion.f90#L16) | 9 | `FLEX_CLS_N_ALS` @16, `FLEX_CLS_N_PPCA` @17, `FLEX_CLS_N_RESTARTS` @18, `FLEX_CLS_RIDGE_REL` @19, `FLEX_CLS_CONV_TOL` @20, `FLEX_CLS_WIENER_EPS` @21, `FLEX_CLS_NUIS_RIDGE` @22, `FLEX_CLS_INV_RIDGE` @23, `FLEX_CLS_GROUP_ALPHA` @24 |
| [simple_flex_cls_kernel.f90:18](../../../src/main/flex/simple_flex_cls_kernel.f90#L18) | 7 | `TAUS` @18, `RANKS` @19, `BLOCK` @20, `NSAMPLE` @21, `RIDGE` @22, `WIENER` @23, `LWMIN` @24 |
| [simple_flex_pca_plot.f90:14](../../../src/main/flex/simple_flex_pca_plot.f90#L14) | 8 | `PANEL` @14, `GUTTER` @15, `NBINS` @16, `DOT` @17, `LEGEND` @18, `PCT` @19, `TAB10` @21, `GREY_NONE` @24 |
| [simple_flex_pca_util.f90:16](../../../src/main/flex/simple_flex_pca_util.f90#L16) | 1 | `COV_SIGNAL_FACTOR` @16 |
| [simple_flex_traj_denoise.f90:24](../../../src/main/flex/simple_flex_traj_denoise.f90#L24) | 18 | `TAUS` @24, `RANKS` @25, `HS` @26, `BLOCKSZ` @27, `NVMAX` @28, `NQMAX` @29, `NMIN` @30, `RIDGE` @31, `WIENER` @32, `LWMIN` @33, `SNR_BAND` @34, `KBAND_MIN` @35, `CHUNK` @36, `NVB_MIN` @37, `NVB_MAX` @38, `LIBCOLS` @39, `LIBCHUNK` @40, `L_BLAS_TEAM_SAFE` @51,53 |
| [simple_flex_tree.f90:14](../../../src/main/flex/simple_flex_tree.f90#L14) | 5 | `TREE_UNASSESSED` @14, `TREE_SPLIT` @14, `TREE_MERGE` @14, `TREE_LEAF` @14, `TREE_JUNK` @14 |
| [simple_umap.f90:11](../../../src/main/flex/simple_umap.f90#L11) | 11 | `NDIM` @11, `NNN` @12, `SMALL_N` @13, `NEPOCH_SMALL` @14, `NEPOCH_LARGE` @14, `NBISECT` @15, `NEG_RATE` @16, `CLIP` @17, `INIT_SCALE` @18, `KERN_A` @21, `KERN_B` @21 |
| [simple_flex_pca_deconv.f90:16](../../../src/main/flex/states/simple_flex_pca_deconv.f90#L16) | 5 | `XD_MAXIT` @16, `XD_TOL` @17, `XD_CV_MAX` @18, `XD_PI_MIN` @19, `XD_KMAX` @20 |
| [simple_flex_pca_merge.f90:24](../../../src/main/flex/states/simple_flex_pca_merge.f90#L24) | 6 | `VIEW_DOF` @24, `VIEW_CHI2_CRIT` @28, `VIEW_MAD_K` @34, `MAD2SIGMA` @36, `VIEW_FOLD_MIN_R` @41, `MERGE_R_DEFAULT` @44 |

### src/main/image

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_ft_expanded.f90:11](../../../src/main/image/simple_ft_expanded.f90#L11) | 2 | `num` @11, `JJ` @12 |
| [simple_ftexp_shsrch.f90:13](../../../src/main/image/simple_ftexp_shsrch.f90#L13) | 2 | `TOL` @13, `MAXITS` @14 |
| [simple_image.f90:20](../../../src/main/image/simple_image.f90#L20) | 2 | `FFTW_MEASURE_TIME_LIMIT` @20, `FFTW_NO_TIME_LIMIT` @21 |
| [simple_image_msk.f90:15](../../../src/main/image/simple_image_msk.f90#L15) | 2 | `DEBUG` @15, `L_WRITE` @16 |
| [simple_image_polar.f90:6](../../../src/main/image/simple_image_polar.f90#L6) | 2 | `POLMEM_NORMAL` @6, `POLMEM_OVERSAMP` @6 |
| [simple_image_seg.f90:6](../../../src/main/image/simple_image_seg.f90#L6) | 1 | `DT_INF` @6 |
| [simple_wave_imaging.f90:19](../../../src/main/image/simple_wave_imaging.f90#L19) | 4 | `WAVE_PROJECTION` @19, `WAVE_FIRST_BORN` @20, `WAVE_RYTOV` @21, `WAVE_MSLICE` @22 |

### src/main/interp

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_winfuns.f90:64](../../../src/main/interp/simple_winfuns.f90#L64) | 3 | `ps` @64, `qs` @65, `thresh` @67 |

### src/main/learn

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_class_student.f90:21](../../../src/main/learn/simple_class_student.f90#L21) | 3 | `CLASS_STUDENT_VERSION` @21, `CAND_MASS` @22, `MODEL_TITLE` @23 |
| [simple_hhcnn.f90:42](../../../src/main/learn/simple_hhcnn.f90#L42) | 22 | `CNN_NBLK` @42, `CNN_NREAD` @43, `CNN_KREAD` @44, `CNN_MREAD` @45, `CNN_LT` @46, `CNN_NQ` @47, `CNN_NCHUNK` @48, `CNN_SLICE` @49, `CNN_LOOP_ROWS` @50, `CNN_MATMUL_ROWS` @51, `CNN_TILE` @52, `CNN_TCOLS` @53, `CNN_WEIGHTS_MAGIC` @54, `CNN_WEIGHTS_VERSION` @55, `SCALE_INIT_LOGIT` @56, `GATE_DELTA2` @57, `NORM_EPS` @58, `ADAM_BETA1` @59, `ADAM_BETA2` @60, `ADAM_EPS` @61, `L_BLAS_TEAM_SAFE` @67,69, `L_MATMUL_PRODUCTS` @72,74 |
| [simple_hhcnn_bench.f90:7](../../../src/main/learn/simple_hhcnn_bench.f90#L7) | 3 | `BENCH_BATCH` @7, `BENCH_NONE` @8, `BENCH_NPRED` @9 |
| [simple_polar_feats.f90:23](../../../src/main/learn/simple_polar_feats.f90#L23) | 25 | `POLAR_FEATS_NHEAD` @23, `FEATS_CACHE_MAGIC` @24, `FEATS_CACHE_VERSION` @25, `IH_MAGIC` @26, `IH_VERSION` @26, `IH_PATH` @26, `IH_MTIME` @26, `IH_ROWS` @26, `IH_BOX` @26, `IH_BOXCROP` @26, `IH_SMPD` @26, `IH_NSEL` @26, `IH_SELHASH` @26, `IH_HIGHROW` @26, `IH_PFTSZ` @26, `IH_KFROM` @26, `IH_KTO` @26, `IH_WHASH` @26, `IH_CTF` @26, `FREE_SPACE_DENOM` @29, `NET_ANGLE_OVERSAMPLING` @30, `HASH_M1` @32, `HASH_A1` @32, `HASH_M2` @33, `HASH_A2` @33 |
| [simple_polar_student.f90:24](../../../src/main/learn/simple_polar_student.f90#L24) | 7 | `RECORD_CHUNK` @24, `BATCH_MAX` @25, `BATCH_MIN` @26, `BATCH_FRAC` @27, `MAX_EXTRA` @28, `LEARNING_RATE` @29, `WEIGHT_DECAY` @30 |
| [simple_pose_student.f90:23](../../../src/main/learn/simple_pose_student.f90#L23) | 3 | `POSE_STUDENT_VERSION` @23, `CONE_MASS` @24, `MODEL_TITLE` @25 |

### src/main/motion

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_dose_series_manifest.f90:13](../../../src/main/motion/simple_dose_series_manifest.f90#L13) | 8 | `VERSION_LINE` @13, `LINE_CHUNK` @14, `MAX_LINE` @14, `MAX_PATH` @14, `SMPD_ATOL` @15, `VALUE_RTOL` @15, `KV_ATOL` @15, `ANGLE_ATOL` @16 |
| [simple_motion_align_hybrid.f90:15](../../../src/main/motion/simple_motion_align_hybrid.f90#L15) | 11 | `SMALLSHIFT` @15, `SRCH_TOL` @16, `CFTOL` @17, `CGTOL` @17, `MINITS_DCORR` @18, `MAXITS_DCORR` @18, `MINITS_CORR` @19, `MAXITS_CORR` @19, `POLYDIM` @20, `NRESUPDATES` @21, `LBFGSB_MAXITS` @22 |
| [simple_motion_align_nano.f90:12](../../../src/main/motion/simple_motion_align_nano.f90#L12) | 2 | `MINITS_DCORR` @12, `MAXITS_DCORR` @12 |
| [simple_motion_correct.f90:63](../../../src/main/motion/simple_motion_correct.f90#L63) | 2 | `NSIGMAS` @63, `FITSHIFTS` @64 |
| [simple_motion_correct_iter.f90:15](../../../src/main/motion/simple_motion_correct_iter.f90#L15) | 3 | `PATCH_FIT_THRESHOLD` @15, `SMPD4VIZ_NANO` @16, `speckind` @17 |
| [simple_motion_gain_analysis.f90:12](../../../src/main/motion/simple_motion_gain_analysis.f90#L12) | 8 | `GAIN_LP_RES_A` @12, `GAIN_FIRST_ANALYSIS_AT` @13, `GAIN_ANALYSIS_STEP` @14, `GAIN_CONV_STREAK_REQ` @15, `GAIN_CONV_EPS_DELTA` @16, `GAIN_CONV_EPS_SMOOTH_CORR` @17, `GAIN_CONV_EPS_SMOOTH_DELTA` @18, `GAIN_CONV_SMOOTH_ALPHA` @19 |
| [simple_motion_gain_helpers.f90:11](../../../src/main/motion/simple_motion_gain_helpers.f90#L11) | 1 | `EER_THUMB_UPSAMPLING` @11 |
| [simple_motion_model.f90:14](../../../src/main/motion/simple_motion_model.f90#L14) | 3 | `MODEL_VERSION` @14, `FILE_VERSION` @15, `MODELSZ` @16 |
| [simple_motion_patched.f90:20](../../../src/main/motion/simple_motion_patched.f90#L20) | 3 | `TOL` @20, `TRS_DEFAULT` @21, `PATCH_PDIM` @22 |

### src/main/nano

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_atomic_species.f90:22](../../../src/main/nano/simple_atomic_species.f90#L22) | 12 | `CN_CLOSE_PACKED` @22, `V_UNDEFINED` @23, `INT_FLOOR_FRAC` @24, `APER_NN_FRAC` @25, `NSHELL_MAX` @26, `NSHELL_ATOMS` @27, `LOW_SNR` @28, `INNER_FRAC` @29, `SURF_CONFINED` @30, `NTAU_ITER` @31, `TAU_TOL` @32, `NTAU_MIN` @33 |
| [simple_atoms.f90:14](../../../src/main/nano/simple_atoms.f90#L14) | 5 | `pdbfmt` @14, `xpdbfmt` @15, `pdbfmt_write` @16, `xpdbfmt_write` @17, `pdbfmt_anisou` @18 |
| [simple_calpha_finder.f90:13](../../../src/main/nano/simple_calpha_finder.f90#L13) | 5 | `PI_LOCAL` @13, `N_CA_BOND` @15, `CA_C_BOND` @16, `N_CA_C_ANGLE` @17, `MIN_TARGET_SIGMA` @18 |
| [simple_nano_density.f90:30](../../../src/main/nano/simple_nano_density.f90#L30) | 32 | `CN_FULL` @30, `CONTACT_NN_FRAC` @31, `FIT_NN_FRAC` @32, `MOVE_NN_FRAC` @33, `BSTART_NN_FRAC` @34, `BFAC_MIN` @35, `BFAC_MAX` @36, `NGRID_B` @37, `NGOLD_B` @38, `NSWEEPS` @39, `NROUNDS_MAX` @44, `CURVE_TOL` @45, `ROUND_TOL` @46, `NROUND_FIELDS` @47, `TRANSFER_DS` @48, `POWER_FRAC` @49, `TAB_DR_VOX` @50, `TAB_R_MAX` @51, `NTAB_B` @52, `TAB_DS_INT` @53, `REACH_FRAC` @54, `FFT_NTHR_MAX` @55, `GAUGE_RES_LO` @60, `GAUGE_RES_HI` @61, `PRIOR_SD_START` @62, `TAU_B_MIN` @63, `STRONG_SNR` @64, `NBIN_MIN` @65, `NPHANTOM` @66, `RENDER_NSD` @67, `GOLD` @69, `L_CHECK_RESIDUAL` @70 |
| [simple_nano_species.f90:20](../../../src/main/nano/simple_nano_species.f90#L20) | 16 | `MAX_NSPECIES` @20, `MIN_SEPARATION` @21, `MIN_CLASS_ATOMS` @22, `MIN_CLASS_FRAC` @23, `MAXITS_EM` @24, `TOL_EM` @25, `PI_MIN_EM` @26, `U_MAX` @27, `BFAC_MIN` @28, `BFAC_MAX` @29, `NGRID_B` @30, `NGOLD_B` @31, `NITER_SLOPE` @32, `SLOPE_TOL_REL` @33, `TOL_SLOPE_EM` @34, `TOL_SLOPE_FIN` @35 |
| [simple_nanoparticle.f90:22](../../../src/main/nano/simple_nanoparticle.f90#L22) | 20 | `NBIN_THRESH` @22, `NVOX_THRESH` @23, `WRITE_OUTPUT` @24, `ATOMS_STATS_OMIT` @25, `SOFT_EDGE` @26, `CNMIN` @27, `CNMAX` @28, `NSTRAIN_COMPS` @29, `ATOMS_STATS_FILE` @30, `NP_STATS_FILE` @31, `CN_STATS_FILE` @32, `MSKDIAM_OVER_DIAM` @33, `B_REF_PSEUDO` @36, `RADIUS_NN_FRAC` @37, `SPLIT_EXCL_NN_FRAC` @38, `CSCORE_INNER_FRAC` @39, `ATOM_STATS_HEAD` @41, `ATOM_STATS_HEAD_OMIT` @49, `NP_STATS_HEAD` @56, `CN_STATS_HEAD` @78 |
| [simple_nanoparticle_utils.f90:22](../../../src/main/nano/simple_nanoparticle_utils.f90#L22) | 4 | `DEBUG` @22, `NSTRAIN_COMPS` @23, `RMAX_NN_FRAC` @24, `U_WZ` @25 |
| [simple_tseries_graphene_subtr.f90:11](../../../src/main/nano/simple_tseries_graphene_subtr.f90#L11) | 4 | `angular_threshold` @11, `REMOVAL_HWIDTH1` @12, `REMOVAL_HWIDTH2` @13, `REMOVAL_HWIDTH3` @14 |
| [single_tseries_tracker.f90:14](../../../src/main/nano/single_tseries_tracker.f90#L14) | 5 | `SPECW` @14, `BSLAMBDA` @15, `BFACTOR` @16, `NNN` @17, `TESTER_WRITE_MULT` @18 |

### src/main/nu_filt

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_nu_filter.f90:44](../../../src/main/nu_filt/simple_nu_filter.f90#L44) | 38 | `NU_BANK_FSC_HEADROOM` @44, `NU_DMAT_CANDIDATE_CAP` @46, `NU_OBJECTIVE_SMOOTH_AWF` @53, `NU_OBJECTIVE_SMOOTH_RADIUS_FRAC` @54, `NU_OBJECTIVE_SMOOTH_MAX_RADIUS_A` @55, `NU_OBJECTIVE_UNARY_CAP` @60, `DISCONT_STEP_THRESH` @63, `NU_LABEL_SMOOTH_MAXITS` @64, `NU_LABEL_SMOOTH_STEP_TOL` @67, `NU_LABEL_SMOOTH_NNEIGH` @68, `NU_LABEL_SMOOTH_NCOLORS` @69, `NU_LABEL_SMOOTH_BETA_FRAC` @70, `NU_LABEL_SMOOTH_QUAD_FRAC` @71, `NU_LABEL_SMOOTH_TIE_EPS` @72, `NU_ENVMASK_BETA` @77, `NU_ENVMASK_DENS_WEIGHT` @79, `NU_ENVMASK_RELATIVE` @81, `NU_ENVMASK_MINVOL_FRAC` @83, `NU_ENVMASK_GROW_A` @89, `NU_ENVMASK_EDGE_A` @91, `NU_ENVMASK_EXCLUDED_SCORE` @94, `NU_ENVMASK_MIN_NULL_VOX` @97, `NU_ENVMASK_MIN_NULL_FRAC` @98, `NU_EVIDENCE_INVALID` @102, `NU_EVIDENCE_NBANDS` @105, `NU_EVIDENCE_BAND_LIMITS` @106, `NU_EVIDENCE_MAX_NBANDS` @110, `NU_EVIDENCE_BAND_RATIO` @111, `NU_EVIDENCE_MIN_BAND_SUPPORT` @113, `NU_ALIGN_LP_MIN_ASSIGNED_PCT` @116, `NU_ALIGN_LP_MIN_SIGNAL_PCT` @120, `NU_EVIDENCE_UNCERTAIN_ENTROPY` @121, `NU_SHARP_BFAC_FINEST_A` @126, `NU_EVIDENCE_MIN_NULL_FRAC` @130, `NU_EVIDENCE_MAX_NULL_FRAC` @131, `NU_EVIDENCE_SOURCE_BASE` @132, `NU_EVIDENCE_SOURCE_PREV` @135, `NU_EVIDENCE_ALGORITHM` @137 |
| [simple_nu_filter_vars.f90:8](../../../src/main/nu_filt/simple_nu_filter_vars.f90#L8) | 2 | `NU_LABEL_KIND` @8, `lowpass_limits` @9 |

### src/main/opt

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_opt_de.f90:12](../../../src/main/opt/simple_opt_de.f90#L12) | 15 | `N_COMMON` @12, `N_GENERAL` @12, `N_VALLEY` @12, `N_MULTIMODAL` @12, `N_FLAT` @12, `F_COMMON` @13, `F_GENERAL` @13, `F_VALLEY` @13, `F_MULTIMODAL` @13, `F_FLAT` @13, `X_COMMON` @14, `X_GENERAL` @14, `X_VALLEY` @14, `X_MULTIMODAL` @14, `X_FLAT` @14 |
| [simple_opt_lbfgsb.f90:43](../../../src/main/opt/simple_opt_lbfgsb.f90#L43) | 10 | `zero` @43, `one` @43, `p5` @43, `p66` @43, `xtrapl` @44, `xtrapu` @44, `big` @45, `xtol_lbfgsb` @46, `two` @47, `three` @47 |
| [simple_pcg_solver.f90:12](../../../src/main/opt/simple_pcg_solver.f90#L12) | 4 | `PCG_STOP_INDEFINITE` @12, `PCG_XTOL` @13, `PCG_RESID_REPLACE` @14, `PCG_RHO_FLOOR_FRAC` @15 |

### src/main/ori

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_ori.f90:127](../../../src/main/ori/simple_ori.f90#L127) | 1 | `zvec` @127 |
| [simple_ori_utils.f90:15](../../../src/main/ori/simple_ori_utils.f90#L15) | 1 | `zvec` @15 |
| [simple_sym.f90:79](../../../src/main/ori/simple_sym.f90#L79) | 11 | `ntet` @79, `noct` @80, `nico` @81, `delta2` @82, `delta3` @83, `delta5` @84, `deltan` @85, `gamma` @86, `tet_theta_lim` @87, `oct_theta_lim` @88, `ico_theta_lim` @89 |

### src/main/pftc

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_pftc_shsrch_grad.f90:18](../../../src/main/pftc/simple_pftc_shsrch_grad.f90#L18) | 10 | `coarse_num_steps` @18, `SHSRCH_SHIFT` @19, `SHSRCH_JOINT` @19, `FIT_EDGE` @21, `FIT_CURVATURE` @21, `FIT_STEP` @21, `FIT_LIMITS` @21, `FIT_VALUE` @21, `FIT_NREASONS` @21, `FIT_MIN_CURV_RATIO` @22 |
| [simple_polar_shiftmap.f90:17](../../../src/main/pftc/simple_polar_shiftmap.f90#L17) | 1 | `SHIFTMAP_JTOL` @17 |
| [simple_polarft_calc.f90:15](../../../src/main/pftc/simple_polarft_calc.f90#L15) | 9 | `REF_EVEN` @15, `REF_ODD` @16, `N_REF_PARITIES` @17, `SHERRSQ` @18, `SHIFTMAP_FULL` @20, `SHIFTMAP_SCREEN` @20, `SHIFTMAP_WINDOW` @20, `NSHIFTMAPS` @20, `zero` @1137 |

### src/main/pick

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_picker_utils.f90:14](../../../src/main/pick/simple_picker_utils.f90#L14) | 6 | `SMPD_SHRINK1` @14, `SMPD_SHRINK2` @14, `GAUSIG` @15, `OFFSET` @16, `L_WRITE` @17, `L_DEBUG` @18 |
| [simple_pickref.f90:14](../../../src/main/pick/simple_pickref.f90#L14) | 2 | `LP_PICK` @14, `OFFSET_DEFAULT` @15 |
| [simple_pickref_corr_batch.f90:12](../../../src/main/pick/simple_pickref_corr_batch.f90#L12) | 1 | `DEFAULT_BATCH_SIZE` @12 |
| [simple_picksegdiam.f90:15](../../../src/main/pick/simple_picksegdiam.f90#L15) | 2 | `SMPD_SHRINK1` @15, `FRAC_FG` @16 |
| [simple_segdiam_bin_picker.f90:31](../../../src/main/pick/simple_segdiam_bin_picker.f90#L31) | 7 | `MOLDIAMS_PICK` @31, `NBINS` @32, `SMPD_SHRINK1` @33, `SIGMA_CRIT` @34, `BOXFAC_MAX` @35, `BOXFAC_MIN` @36, `BOXFAC_DECAY_PX` @37 |

### src/main/preprocess

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_particle_extractor.f90:16](../../../src/main/preprocess/simple_particle_extractor.f90#L16) | 3 | `POLYDIM` @16, `NSIGMAS` @17, `DEBUG_HERE` @18 |

### src/main/project

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_class_candidates.f90:20](../../../src/main/project/simple_class_candidates.f90#L20) | 8 | `IMGKIND` @20, `MAGIC` @21, `MAGIC_PART` @22, `VERSION` @23, `CANDS_OK` @24, `CANDS_NONE` @25, `CANDS_UNREADABLE` @26, `CANDS_STALE` @27 |
| [simple_eo_partition.f90:18](../../../src/main/project/simple_eo_partition.f90#L18) | 11 | `EO_COMPLETE` @18, `EO_MISSING` @18, `EO_PARTIAL` @18, `EO_INCONSISTENT` @18, `EO_SEED_DEFAULT` @19, `MINSTD_A` @20, `MINSTD_M` @20, `MIN_GROUPS_WARN` @21, `IMBALANCE_WARN` @22, `KEY_GROUP` @23, `KEY_INDEX` @23 |
| [simple_project_superset.f90:18](../../../src/main/project/simple_project_superset.f90#L18) | 2 | `MIN_COHORT_STATE_POP` @18, `COHORT_WARN_FRAC` @19 |
| [simple_sp_project.f90:24](../../../src/main/project/simple_sp_project.f90#L24) | 1 | `MAXN_OS_SEG` @24 |
| [simple_state_weight_set.f90:29](../../../src/main/project/simple_state_weight_set.f90#L29) | 2 | `ROW_SUM_TOL` @29, `POPULATION_TOL` @31 |

### src/main/quality

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_quality_logistic.f90:22](../../../src/main/quality/simple_quality_logistic.f90#L22) | 13 | `QUALITY_MAX_FEATS` @22, `QUALITY_MAX_INTERACTIONS` @23, `LEARN_ROLE_SKIP` @24, `LEARN_ROLE_BALANCED` @25, `LEARN_ROLE_RECALL_ONLY` @26, `LEARN_ROLE_SPECIFICITY_ONLY` @27, `LEARN_ROBUST_TAIL_FRAC` @28, `LEARN_ROBUST_TAIL_WEIGHT` @29, `LEARN_BALANCED_GOOD_LOSS_WEIGHT` @30, `LEARN_BALANCED_BAD_LOSS_WEIGHT` @31, `LOGISTIC_LAMBDAS` @32, `LOGISTIC_THRESHOLDS` @34, `LOGISTIC_COEFF_ABS_BOUND` @38 |

### src/main/sieve

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_ptcl_sieve.f90:65](../../../src/main/sieve/simple_ptcl_sieve.f90#L65) | 25 | `DEBUG` @65, `DEFAULT_COARSE_POP_THRESHOLD` @66, `DEFAULT_FINE_POP_THRESHOLD` @67, `DEFAULT_NCLS` @68, `DEFAULT_WALLTIME` @69, `DEFAULT_COARSE_BOX` @70, `DEFAULT_FINE_BOX` @71, `DEFAULT_COARSE_NSAMPLE` @72, `DEFAULT_FINE_NSAMPLE` @73, `DEFAULT_LPSTART` @74, `DEFAULT_COARSE_LP` @75, `DEFAULT_FINE_LP` @76, `LABEL_COARSE` @79, `LABEL_FINE` @80, `REJECTION_FAILED` @81, `FINAL_INGESTION_MARKER` @84, `CHUNK_JOB_LABEL` @87, `CHUNK_EXIT_CODE` @88, `CHUNK_INPUT_PROJFILE` @89, `CHUNK_ATTEMPT1` @90, `CHUNKED_MICS` @93, `CHUNK_MICS` @96, `CONSUMED_COARSE` @99, `COARSE_PREFIX` @101, `FINE_PREFIX` @102 |

### src/main/sigma2

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_sigma2_anchor.f90:12](../../../src/main/sigma2/simple_sigma2_anchor.f90#L12) | 2 | `MAGIC` @12, `HEADER_BYTES` @13 |
| [simple_sigma2_bootstrap.f90:20](../../../src/main/sigma2/simple_sigma2_bootstrap.f90#L20) | 1 | `PSPEC_TEMPLATE_KEYS` @20 |
| [simple_sigma2_pspec_sums.f90:21](../../../src/main/sigma2/simple_sigma2_pspec_sums.f90#L21) | 1 | `PSPEC_SUMS_MAGIC` @21 |
| [simple_sigma2_state.f90:31](../../../src/main/sigma2/simple_sigma2_state.f90#L31) | 2 | `BLOCK_ROWS` @31, `SIGMA2_RANGE_FBODY` @32 |

### src/main/simple_cmdline.f90

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_cmdline.f90:13](../../../src/main/simple_cmdline.f90#L13) | 3 | `MAX_CMDARGS` @13, `DEBUG_HERE` @14, `INPROCESS_KEYS` @19 |

### src/main/simple_convergence.f90

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_convergence.f90:13](../../../src/main/simple_convergence.f90#L13) | 2 | `STAGE_PLATEAU_PATIENCE` @13, `STAGE_PLATEAU_MI_GAIN` @14 |

### src/main/simple_refine3D_startup.f90

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_refine3D_startup.f90:8](../../../src/main/simple_refine3D_startup.f90#L8) | 3 | `STARTUP_ORDINARY` @8, `STARTUP_REUSE_CHAIN` @9, `STARTUP_SEED_CHAIN` @10 |

### src/main/simple_refine3D_states_controller.f90

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_refine3D_states_controller.f90:23](../../../src/main/simple_refine3D_states_controller.f90#L23) | 10 | `NSAMPLE_PER_STATE_REFINE3D_STATES` @23, `NSAMPLE_REFINE3D_STATES_CAP` @24, `STAGE2_NSPACE` @25, `STAGE2_NSPACE_SUB` @26, `FREQUENCY_BLOCK_NITS` @27, `STAGE2_MINITS` @28, `MINITS_REFINE3D_STATES` @29, `MAXITS_REFINE3D_STATES_CAP` @30, `TARGET_UPDATES_PER_PARTICLE_REFINE3D_STATES` @31, `MIN_STATE_FRAC_FLEX` @32 |

### src/main/solve

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_solve2D_chunk.f90:51](../../../src/main/solve/simple_solve2D_chunk.f90#L51) | 1 | `DEBUG_HERE` @51 |
| [simple_solve2D_controller.f90:21](../../../src/main/solve/simple_solve2D_controller.f90#L21) | 16 | `SMPD_TARGET` @21, `MINBOXSZ` @22, `NSTAGES_CLS` @24, `ITS_INCR` @25, `PHASES` @26, `EXTR_LIM_LOCAL` @27, `PROBREFINE_STAGE` @28, `STICKY_SAMPL_STAGE` @29, `FRAC_UPDATE_STAGE` @30, `FIRST_SHIFT_STAGE` @31, `GOLD_STD_STAGE` @32, `NSAMPLE_DEFAULT_2D` @33, `COVERAGE_NITS_CAP` @34, `EO_STAGE` @35, `SKIP_MIN_STAGE` @37, `FRC05_SKIP_PTCL_FRAC` @38 |
| [simple_solve3D_addon_report.f90:17](../../../src/main/solve/simple_solve3D_addon_report.f90#L17) | 11 | `ADDON_REPORT_FNAME` @17, `SHELL_TOL` @18, `FSC_CRIT` @19, `DOCK_CORR_FLOOR` @20, `DOCK_HP` @21, `DOCK_LP_MIN` @22, `VERDICT_LEN` @23, `V_NONE` @24, `V_UP` @24, `V_SAME` @24, `V_DOWN` @24 |
| [simple_solve3D_chunks_plan.f90:17](../../../src/main/solve/simple_solve3D_chunks_plan.f90#L17) | 6 | `CHUNKS_REPORT_FNAME` @17, `MAX_STATE_ARTIFACTS` @18, `VERDICT_LEN` @19, `NSAMPLE_TOTAL_CHUNK` @20, `NPTCLS_CHUNK_MIN` @21, `NPTCLS_CHUNK_MAX` @22 |
| [simple_solve3D_controller.f90:7](../../../src/main/solve/simple_solve3D_controller.f90#L7) | 43 | `REC_FBODY` @7, `NSTAGES` @10, `NSTAGES_INI3D_MAX` @11, `MAXITS` @12, `NSPACE` @13, `NSPACE_CAVGS_EARLY` @14, `MAXITS_CAVGS_EARLY` @15, `CAVGS_EARLY_NSTAGES` @16, `NSPACE_REGISTER` @17, `NSPACE_SUB` @18, `NSPACE_SUB_BASE` @19, `NSPACE_MULTI_MAX` @20, `GAUREF_LAST_STAGE` @23, `ML_REG_START_STAGE` @24, `SYMSRCH_STAGE` @25, `NOSYMSRCH_STAGE` @26, `VOL_SUPPLIED_STAGE` @27, `PROB_REFINE_STAGE` @28, `TRAILREC_STAGE_SINGLE` @29, `NU_FILTER_STAGE` @30, `GOLD_STD_STAGE` @31, `PCG_START_STAGE` @32, `PROB_NEIGH_REFINE_STAGE` @33, `NSTAGES_INDEPENDENT` @34, `FRACSRCH_SCAN_EXIT` @35, `TRAILREC_STAGE_MULTI` @36, `PROB_NEIGH_MODE_EARLY` @37, `PROB_NEIGH_MODE_LATE` @38, `PROB_NEIGH_MODE_MULTI` @39, `LPSTOP_BOUNDS` @42, `FSC05_PROMOTE_MIN_STAGE` @43, `STAGE_MINITS` @44, `LPSTART_MSKDIAM_DIV` @45, `LPSTART_SIZE_BOUNDS` @46, `LPSTOP_MSKDIAM_DIV` @47, `CENLP_DEFAULT` @48, `LPSYMSRCH_LB` @49, `LPSTART_INI3D` @50, `LPSTOP_INI3D` @51, `LPSTOP_INDEPENDENT` @52, `UPDATE_FRAC_MAX` @55, `FULL_SAMPLE_SWITCH_FRAC` @56, `NSAMPLE_SOLVE3D_DEFAULT` @57 |
| [simple_solve3D_manifest.f90:23](../../../src/main/solve/simple_solve3D_manifest.f90#L23) | 9 | `MANIFEST_SCHEMA` @23, `MANIFEST_VERSION` @24, `MANIFEST_FNAME` @25, `MANIFEST_PROJINFO_KEY` @26, `RUN_ID_PROJINFO_KEY` @27, `KLEN` @28, `MANIFEST_INPUT_KEYS` @33, `MANIFEST_REPLAY_KEYS` @48, `MANIFEST_OVERRIDABLE_KEYS` @59 |
| [simple_solve3D_utils.f90:36](../../../src/main/solve/simple_solve3D_utils.f90#L36) | 3 | `FINAL_PCG_MAXITS_FLOOR` @36, `RESEED_MIN_POP` @38, `RESEED_MIN_FRAC` @39 |
| [simple_state_assignment.f90:79](../../../src/main/solve/simple_state_assignment.f90#L79) | 26 | `MIN_PTCLS_PER_STATE` @79, `NCLS2D_FIRST` @83, `NCLS_EXPANDED` @84, `LPSTART_CAVGS3D` @86, `LPSTOP_CAVGS3D` @87, `EO_ANGDIST_MAX` @89, `EO_OVERLAP_MIN` @90, `STATE_FRAC_MIN` @91, `NSTATES_MIN` @92, `REFS_PAIR_CC_MAX` @93, `DOCK_HP_REFS` @94, `LP_ASSIGN` @96, `NSPACE_ASSIGN` @97, `FIT_NMAD` @98, `FIT_REJECT_MAX` @99, `MARGIN_AMBIGUOUS` @100, `SOLVE2D_DIR` @102, `CLSEXP_DIR` @103, `CAVGS3D_DIR` @104, `ASSIGN_LABEL` @105, `REFINE_DIR` @106, `REFS_DIR` @107, `REFS_RECORD` @108, `QUALITY_DIR` @109, `SELECTED_CAVGS` @110, `REJECTED_CAVGS` @111 |

### src/main/star

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_relion.f90:33](../../../src/main/star/simple_relion.f90#L33) | 8 | `MIC_EXT_LEN` @33, `FORCTF_EXT_LEN` @34, `STK_EXT_LEN` @35, `STK_PREFIX_LEN` @36, `MICRON2ANGS` @37, `MICPATH` @38, `PTCLPATH` @39, `MOVIEPATH` @40 |
| [simple_starproject_utils.f90:13](../../../src/main/star/simple_starproject_utils.f90#L13) | 2 | `LEN_FLAG` @13, `RELION_PHASE_DEG2RAD` @14 |

### src/main/strategies

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_calc_pspec_strategy.f90:18](../../../src/main/strategies/parallelization/simple_calc_pspec_strategy.f90#L18) | 2 | `SIGMA2_BOOTSTRAP_MAX_PARTICLES` @18, `SIGMA2_BOOTSTRAP_SELECTION_FNAME` @19 |
| [simple_cls_expansion_strategy.f90:30](../../../src/main/strategies/parallelization/simple_cls_expansion_strategy.f90#L30) | 11 | `FLEX_CLS_NEIGS_DEFAULT` @30, `FLEX_CLS_SEED_BASE` @31, `FLEX_CLS_BAND_FIRST_LP` @32, `FLEX_CLS_MAP_FILE` @33, `FLEX_CLS_ALLOC_FILE` @34, `FLEX_CLS_GROUPS_FBODY` @35, `FLEX_CLS_CAVGS_FILE` @36, `FLEX_CLS_CAVGS_EVEN` @37, `FLEX_CLS_CAVGS_ODD` @38, `FLEX_CLS_FRCS_TXT` @39, `FLEX_CLS_FRCS` @40 |
| [simple_rec3D_pcg_strategy.f90:46](../../../src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90#L46) | 4 | `PCG_MASTER_NTHR_CAP` @46, `DEBUG` @51, `ATM_AMP_MAXITS` @53, `ATM_AMP_RTOL` @54 |
| [simple_rec3D_service.f90:28](../../../src/main/strategies/parallelization/simple_rec3D_service.f90#L28) | 11 | `REC3D_BACKEND_INVALID` @28, `REC3D_BACKEND_GRIDDING` @29, `REC3D_BACKEND_PCG` @30, `REC3D_WEIGHTS_HARD` @32, `REC3D_WEIGHTS_SET` @34, `REC3D_WEIGHTS_TABLE` @36, `REC3D_OUTPUT_RAW` @38, `REC3D_OUTPUT_DELIVERED` @39, `REC3D_DISPATCH_INPROC` @40, `REC3D_DISPATCH_QUEUE` @41, `PINDS_FNAME` @42 |
| [simple_rec3D_wave_strategy.f90:40](../../../src/main/strategies/parallelization/simple_rec3D_wave_strategy.f90#L40) | 3 | `STEP_SUBSET_FRAC` @40, `STEP_SUBSET_MIN` @41, `STEP_SUBSET_MAX` @42 |
| [simple_refine3D_strategy.f90:36](../../../src/main/strategies/parallelization/simple_refine3D_strategy.f90#L36) | 1 | `PROB_MULTISTATE_MIN_POP` @36 |
| [simple_refine_motion_model_strategy.f90:24](../../../src/main/strategies/parallelization/simple_refine_motion_model_strategy.f90#L24) | 5 | `DOSE_FRACTION_RELTOL` @24, `DOSE_FRACTION_ABSTOL` @25, `DOSE_VOLTAGE_TOL` @26, `DOSE_PARTITION_ROW_WIDTH` @27, `DOSE_MANIFEST_FNAME` @28 |
| [simple_volassemble_strategy.f90:33](../../../src/main/strategies/parallelization/simple_volassemble_strategy.f90#L33) | 7 | `TC_NREP` @33, `TC_NSMP` @34, `TC_NNEW` @35, `TC_NSEED` @36, `TC_NCREP` @37, `TC_NCSEED` @38, `TC_N` @39 |
| [simple_eul_prob_tab2D.f90:27](../../../src/main/strategies/search/probabilistic/simple_eul_prob_tab2D.f90#L27) | 1 | `NHOOD_FRAC` @27 |
| [simple_eul_prob_tab_neigh.f90:33](../../../src/main/strategies/search/probabilistic/simple_eul_prob_tab_neigh.f90#L33) | 1 | `MIN_STOCH_FINE_REFS` @33 |
| [simple_eul_prob_tab_utils.f90:54](../../../src/main/strategies/search/probabilistic/simple_eul_prob_tab_utils.f90#L54) | 12 | `SHIFT_PRIOR_SIGMA_MIN` @54, `SHIFT_PRIOR_LP_FRAC` @55, `SCREEN_QUANTILE` @58, `PROB_NEFF_ROWS` @81, `PROB_NEFF_NMIN` @82, `PROB_NEFF_MIN_ROWS` @83, `PROB_READOUT_3D_DENSE` @132, `PROB_READOUT_3D_NEIGH` @132, `PROB_READOUT_2D_DENSE` @132, `PROB_READOUT_2D_SPARSE` @132, `PROB_ROW_FILE_MAGIC` @300, `PROB_ROW_FILE_VERSION` @301 |
| [simple_matcher_ptcl_io.f90:29](../../../src/main/strategies/search/simple_matcher_ptcl_io.f90#L29) | 5 | `BUFFER_FILLED` @29, `BUFFER_OFF` @29, `BUFFER_OVER_BUDGET` @29, `BUFFER_NO_MEMORY` @29, `BUFFER_RECORDS` @30 |
| [simple_strategy2D_inpl_smpl.f90:14](../../../src/main/strategies/search/simple_strategy2D_inpl_smpl.f90#L14) | 1 | `DEBUG` @14 |
| [simple_strategy2D_snhc.f90:13](../../../src/main/strategies/search/simple_strategy2D_snhc.f90#L13) | 1 | `DEBUG` @13 |
| [simple_strategy2D_snhc_smpl.f90:14](../../../src/main/strategies/search/simple_strategy2D_snhc_smpl.f90#L14) | 1 | `DEBUG` @14 |
| [simple_strategy2D_snhc_smpl_many.f90:15](../../../src/main/strategies/search/simple_strategy2D_snhc_smpl_many.f90#L15) | 1 | `DEBUG` @15 |
| [simple_strategy2D_srch.f90:20](../../../src/main/strategies/search/simple_strategy2D_srch.f90#L20) | 4 | `CONT_ROUTE_NOT_ATTEMPTED` @20, `CONT_ROUTE_IMPROVED` @21, `CONT_ROUTE_NO_IMPROVEMENT` @22, `CONT_ROUTE_INVALID` @23 |
| [simple_strategy2D_tseries.f90:15](../../../src/main/strategies/search/simple_strategy2D_tseries.f90#L15) | 1 | `TRSSTEP` @15 |
| [simple_strategy2D_utils.f90:24](../../../src/main/strategies/search/simple_strategy2D_utils.f90#L24) | 11 | `W_HIST_SIG` @24, `W_POW_SIG` @24, `W_SIG` @25, `W_CC` @25, `W_RES` @25, `W_SCORE_HOMO` @26, `W_SCORE_RES` @26, `W_SCORE_CLUST` @26, `NHISTBINS` @28, `HP_SPEC` @29, `LP_SPEC` @29 |
| [simple_strategy3D_srch.f90:15](../../../src/main/strategies/search/simple_strategy3D_srch.f90#L15) | 4 | `CONT_ROUTE_NOT_ATTEMPTED` @15, `CONT_ROUTE_IMPROVED` @16, `CONT_ROUTE_NO_IMPROVEMENT` @17, `CONT_ROUTE_INVALID` @18 |
| [simple_view_partition_sampling.f90:19](../../../src/main/strategies/search/simple_view_partition_sampling.f90#L19) | 6 | `VIEW_PARTITION_FBODY` @19, `VIEW_PARTITION_CRIT` @20, `VIEW_PARTITION_HP` @21, `VIEW_PARTITION_LP` @22, `VIEW_PARTITION_MIN_FRAC` @23, `VIEW_PARTITION_TRS` @24 |

### src/main/stream

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_stream_master_resources.f90:40](../../../src/main/stream/master/simple_stream_master_resources.f90#L40) | 17 | `PREPROCESS_NTHR` @40, `PREPROCESS_NPARTS` @41, `OPTICS_NTHR` @42, `INITIAL_ANALYSIS_NTHR` @43, `INITIAL_ANALYSIS_NTHR2D` @44, `INITIAL_ANALYSIS_NTHR3D` @45, `INITIAL_ANALYSIS_NPARTS` @46, `INITIAL_ANALYSIS_NCHUNKS` @47, `REFPICK_NTHR` @48, `REFPICK_NPARTS` @49, `SIEVE_NTHR` @50, `SIEVE_NCHUNKS` @51, `POOL2D_NTHR` @52, `POOL2D_NPARTS` @53, `SOLVE3D_NTHR` @54, `SOLVE3D_NTHR3D` @55, `SOLVE3D_NPARTS3D` @56 |
| [simple_stream_master_stage.f90:45](../../../src/main/stream/master/simple_stream_master_stage.f90#L45) | 1 | `UPDATE_PARTIAL_RETRIES` @45 |
| [simple_stream_master_stage_ids.f90:34](../../../src/main/stream/master/simple_stream_master_stage_ids.f90#L34) | 11 | `STAGE_PREPROCESS` @34, `STAGE_ASSIGN_OPTICS` @35, `STAGE_INITIAL_ANALYSIS` @36, `STAGE_REFERENCE_PICKING` @37, `STAGE_PARTICLE_SIEVING` @38, `STAGE_POOL2D` @39, `STAGE_SOLVE3D` @40, `NSTAGES` @41, `OP_OPEN` @44, `OP_CLOSE` @45, `OP_CLOSE_OTHER` @46 |
| [simple_stream_pool2D.f90:50](../../../src/main/stream/pool2D/simple_stream_pool2D.f90#L50) | 3 | `POOL_NHISTORY` @50, `POOL_JOB_LABEL` @51, `DEBUG_HERE` @52 |
| [simple_stream_job_sets.f90:45](../../../src/main/stream/shared/simple_stream_job_sets.f90#L45) | 1 | `SOURCE_KEY` @45 |
| [simple_stream_pipe.f90:46](../../../src/main/stream/shared/simple_stream_pipe.f90#L46) | 3 | `MAX_SEND_RETRIES` @46, `RETRY_SLEEP_US` @47, `READ_CHUNK_BYTES` @48 |
| [simple_stream_watcher.f90:12](../../../src/main/stream/simple_stream_watcher.f90#L12) | 2 | `RATE_INTERVAL` @12, `HISTORY_CAPACITY0` @13 |
| [simple_stream_stage_initial_analysis.f90:146](../../../src/main/stream/stages/simple_stream_stage_initial_analysis.f90#L146) | 26 | `NMICS_PLAN` @146, `MAX_PROJECTS_IMPORT` @147, `NTHUMB_MAX` @148, `NCLS_MIN` @149, `NCLS_MAX` @149, `NSAMPLE2D` @150, `LPSTOP2D` @151, `EXTRACT_NTHR` @152, `STATE_DOMINANT_FRAC` @157, `STATE_POP_FLOOR` @158, `STATE_RES_FACTOR` @159, `TARGET_NCLS` @160, `PICKREFS_SELECTION` @161, `INIT_SETUP` @164, `INIT_PICK` @164, `INIT_EXTRACT` @164, `INIT_CLASSIFY` @164, `INIT_SELECT` @164, `ALL_COLLECT` @166, `ALL_SIEVE` @166, `ALL_CLASSIFY` @166, `ALL_SELECT` @166, `ALL_BALANCE` @166, `ALL_SOLVE3D` @166, `ALL_REPROJECT` @166, `ALL_DONE` @166 |
| [simple_stream_stage_optics.f90:62](../../../src/main/stream/stages/simple_stream_stage_optics.f90#L62) | 2 | `MAX_PROJECTS_IMPORT` @62, `NMAPS_KEPT` @63 |
| [simple_stream_stage_pool2D.f90:95](../../../src/main/stream/stages/simple_stream_stage_pool2D.f90#L95) | 11 | `NPTCLS_PER_CLS_MIN` @95, `EARLY_RATE_FACTOR` @96, `LATE_RATE_FACTOR` @97, `LATE_ITER` @98, `FINAL_ITER` @99, `MSKDIAM_SWITCH_ITER` @100, `FIRST_EXPORT_ITER` @101, `EXPORT_START_ITER` @104, `NMICS_FIRST3D` @105, `SIEVE_CAVGS_DIR` @108, `NPUBLICATIONS_KEPT` @110 |
| [simple_stream_stage_preprocess.f90:96](../../../src/main/stream/stages/simple_stream_stage_preprocess.f90#L96) | 8 | `GUI_THRESHOLDS` @96, `PLOT_WINDOW` @98, `NTHUMBS` @99, `STAR_EVERY_NMICS` @100, `STAR_STEP_NMICS` @101, `GAIN_BATCH_NMOVIES` @102, `PARTIAL_SET_QUIET_S` @103, `GENERATED_GAINREF` @105 |
| [simple_stream_stage_refpick.f90:92](../../../src/main/stream/stages/simple_stream_stage_refpick.f90#L92) | 5 | `NTHUMB_MAX` @92, `DEFAULT_EXTRACT_BOX` @93, `STAR_EVERY_NMICS` @94, `STAR_STEP_NMICS` @95, `PICKREFS_JOB_DIR` @96 |
| [simple_stream_stage_sieve.f90:77](../../../src/main/stream/stages/simple_stream_stage_sieve.f90#L77) | 1 | `MAX_PROJECTS_IMPORT` @77 |
| [simple_stream_stage_solve3D.f90:120](../../../src/main/stream/stages/simple_stream_stage_solve3D.f90#L120) | 20 | `PHASE_IMPORTING` @120, `PHASE_SOLVE2D` @121, `PHASE_CLSEXP` @122, `PHASE_CAVGS3D` @123, `PHASE_IDLE` @124, `PHASE_EXPAND` @125, `PHASE_ASSIGN` @126, `PHASE_REFINE` @127, `JOB_NONE` @129, `JOB_ROUTE` @130, `JOB_FALLBACK` @131, `JOB_ASSIGN` @132, `ROUTE_TRIES` @135, `L_MODEL_LOGGED` @136, `ASSIGN_DIR` @138, `EXPAND_DIR` @139, `RESULT_PROJ` @140, `NQUALITY_KEPT` @143, `QUALITY_RUNS` @144, `REJECTED_PUBLICATIONS` @146 |

### src/main/ui

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_ui_cavgproc.f90:6](../../../src/main/ui/simple/simple_ui_cavgproc.f90#L6) | 1 | `UI_CATEGORY` @6 |
| [simple_ui_dock.f90:6](../../../src/main/ui/simple/simple_ui_dock.f90#L6) | 1 | `UI_CATEGORY` @6 |
| [simple_ui_filter.f90:6](../../../src/main/ui/simple/simple_ui_filter.f90#L6) | 1 | `UI_CATEGORY` @6 |
| [simple_ui_heterogeneity.f90:6](../../../src/main/ui/simple/simple_ui_heterogeneity.f90#L6) | 1 | `UI_CATEGORY` @6 |
| [simple_ui_image.f90:6](../../../src/main/ui/simple/simple_ui_image.f90#L6) | 1 | `UI_CATEGORY` @6 |
| [simple_ui_mask.f90:6](../../../src/main/ui/simple/simple_ui_mask.f90#L6) | 1 | `UI_CATEGORY` @6 |
| [simple_ui_ori.f90:6](../../../src/main/ui/simple/simple_ui_ori.f90#L6) | 1 | `UI_CATEGORY` @6 |
| [simple_ui_other.f90:6](../../../src/main/ui/simple/simple_ui_other.f90#L6) | 1 | `UI_CATEGORY` @6 |
| [simple_ui_postprocess.f90:6](../../../src/main/ui/simple/simple_ui_postprocess.f90#L6) | 1 | `UI_CATEGORY` @6 |
| [simple_ui_preproc.f90:6](../../../src/main/ui/simple/simple_ui_preproc.f90#L6) | 1 | `UI_CATEGORY` @6 |
| [simple_ui_print.f90:6](../../../src/main/ui/simple/simple_ui_print.f90#L6) | 1 | `UI_CATEGORY` @6 |
| [simple_ui_project.f90:6](../../../src/main/ui/simple/simple_ui_project.f90#L6) | 1 | `UI_CATEGORY` @6 |
| [simple_ui_reconstruct3D.f90:6](../../../src/main/ui/simple/simple_ui_reconstruct3D.f90#L6) | 1 | `UI_CATEGORY` @6 |
| [simple_ui_refine2D.f90:6](../../../src/main/ui/simple/simple_ui_refine2D.f90#L6) | 1 | `UI_CATEGORY` @6 |
| [simple_ui_refine3D.f90:6](../../../src/main/ui/simple/simple_ui_refine3D.f90#L6) | 1 | `UI_CATEGORY` @6 |
| [simple_ui_res.f90:6](../../../src/main/ui/simple/simple_ui_res.f90#L6) | 1 | `UI_CATEGORY` @6 |
| [simple_ui_sim.f90:6](../../../src/main/ui/simple/simple_ui_sim.f90#L6) | 1 | `UI_CATEGORY` @6 |
| [simple_ui_solve3D.f90:6](../../../src/main/ui/simple/simple_ui_solve3D.f90#L6) | 1 | `UI_CATEGORY` @6 |
| [simple_ui_stream.f90:6](../../../src/main/ui/simple/simple_ui_stream.f90#L6) | 1 | `UI_CATEGORY` @6 |
| [simple_ui_sym.f90:6](../../../src/main/ui/simple/simple_ui_sym.f90#L6) | 1 | `UI_CATEGORY` @6 |
| [simple_ui_validation.f90:6](../../../src/main/ui/simple/simple_ui_validation.f90#L6) | 1 | `UI_CATEGORY` @6 |
| [simple_ui_volume.f90:6](../../../src/main/ui/simple/simple_ui_volume.f90#L6) | 1 | `UI_CATEGORY` @6 |
| [simple_private_prgs.f90:24](../../../src/main/ui/simple_private_prgs.f90#L24) | 2 | `NMAX_PRIVATE_PRGS` @24, `NMAX_CMD_DICT` @29 |
| [simple_ui.f90:25](../../../src/main/ui/simple_ui.f90#L25) | 3 | `UI_FNAME` @25, `UI_JSON_REAL_FORMAT` @26, `DEBUG` @27 |
| [simple_ui_param.f90:10](../../../src/main/ui/simple_ui_param.f90#L10) | 1 | `UI_PLACEHOLDER_MAX_LEN` @10 |
| [simple_ui_program.f90:14](../../../src/main/ui/simple_ui_program.f90#L14) | 12 | `UI_IMG` @14, `UI_PARM` @14, `UI_FILE` @14, `UI_SRCH` @14, `UI_FILT` @14, `UI_MASK` @14, `UI_COMP` @14, `UI_SUMMARY_MIN_LEN` @15, `UI_SUMMARY_MAX_LEN` @16, `UI_DISPLAY_NAME_MAX_LEN` @17, `UI_JSON_REAL_FORMAT` @18, `UI_ENVIRONMENT_KEYS` @23 |
| [simple_ui_visibility.f90:6](../../../src/main/ui/simple_ui_visibility.f90#L6) | 3 | `UI_VIS_STANDARD` @6, `UI_VIS_ADVANCED` @7, `UI_VIS_DEVELOPER` @8 |
| [single_ui_atom.f90:6](../../../src/main/ui/single/single_ui_atom.f90#L6) | 1 | `UI_CATEGORY` @6 |
| [single_ui_map.f90:6](../../../src/main/ui/single/single_ui_map.f90#L6) | 1 | `UI_CATEGORY` @6 |
| [single_ui_nano2D.f90:6](../../../src/main/ui/single/single_ui_nano2D.f90#L6) | 1 | `UI_CATEGORY` @6 |
| [single_ui_nano3D.f90:6](../../../src/main/ui/single/single_ui_nano3D.f90#L6) | 1 | `UI_CATEGORY` @6 |
| [single_ui_trajectory.f90:6](../../../src/main/ui/single/single_ui_trajectory.f90#L6) | 1 | `UI_CATEGORY` @6 |
| [single_ui_tseries.f90:6](../../../src/main/ui/single/single_ui_tseries.f90#L6) | 1 | `UI_CATEGORY` @6 |
| [single_ui_validate.f90:6](../../../src/main/ui/single/single_ui_validate.f90#L6) | 1 | `UI_CATEGORY` @6 |

### src/main/volume

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_dock_vols.f90:40](../../../src/main/volume/simple_dock_vols.f90#L40) | 5 | `PGRP` @40, `NSPACE` @41, `NSPACE_SUB` @42, `NBEST` @43, `LP2SMPD_TARGET` @44 |
| [simple_frozen_accum.f90:24](../../../src/main/volume/simple_frozen_accum.f90#L24) | 7 | `FROZEN_CONTEXT_SCHEMA` @24, `FROZEN_SET_SCHEMA` @25, `FROZEN_SCHEMA_VERSION` @26, `FROZEN_BACKEND` @27, `SMPD_RELTOL` @28, `BUNDLE_SCHEMA` @29, `BUNDLE_SCHEMA_VERSION` @30 |
| [simple_halfmap_diagnostics.f90:18](../../../src/main/volume/simple_halfmap_diagnostics.f90#L18) | 5 | `CFAR_NDIRS` @18, `CFAR_CONE_HALF_ANGLE_DEG` @19, `CFAR_FSC_THRESHOLD` @20, `CFAR_MIN_COUNT` @21, `FSC_MODE_UNKNOWN` @22 |
| [simple_halfmap_enhancer.f90:17](../../../src/main/volume/simple_halfmap_enhancer.f90#L17) | 8 | `FSC_CUTOFF` @17, `ENVELOPE_FULL_SUPPORT` @18, `MAX_SHARPEN_B` @19, `MAX_SPECTRAL_GAIN` @20, `TARGET_WINDOW_A` @21, `CORRELATED_NOISE_RATIO` @22, `MIN_WINDOW` @23, `MAX_WINDOW` @24 |
| [simple_order_prior_offline.f90:23](../../../src/main/volume/simple_order_prior_offline.f90#L23) | 7 | `DATA_SCALE_SHELLS` @23, `SUPPORT_MIN` @24, `EMPTY_SHELL_REL` @25, `EMPTY_SHELL_FAC` @26, `CG_RTOL` @27, `CG_XTOL` @28, `CG_MAXITS` @29 |
| [simple_pcg_atomic_prior.f90:16](../../../src/main/volume/simple_pcg_atomic_prior.f90#L16) | 14 | `SIGMA_MIN_VOX` @16, `FOOT_NSIGMA` @17, `FOOT_RAD_MIN` @18, `FOOT_RAD_MAX` @19, `MIN_SEP_A` @20, `GRAM_RIDGE_REL` @21, `GRAM_TOL` @22, `GRAM_MAXITS` @23, `WIDTH_FIT_CORR` @24, `WIDTH_FIT_MINSHELLS` @25, `FIT_SIGMA_VOX` @26, `FIT_NSIGMA` @27, `C2_MAX` @28, `SIGNAL_FLOOR` @29 |
| [simple_pcg_lattice.f90:22](../../../src/main/volume/simple_pcg_lattice.f90#L22) | 1 | `PCG_FFT_MEASURE_MIN_BOXPD` @22 |
| [simple_pcg_order_prior.f90:35](../../../src/main/volume/simple_pcg_order_prior.f90#L35) | 8 | `ORDER_FRAC_GRID` @35, `ORDER_LAMBDA_GRID` @37, `CORE_RES_FAC` @38, `CORE_SIGNAL_MIN` @39, `CORE_MIN_VOX` @40, `DENSITY_LP_FACTOR` @42, `DENSITY_LP_MIN_A` @43, `SUPPORT_MIN` @44 |
| [simple_reconstructor.f90:19](../../../src/main/volume/simple_reconstructor.f90#L19) | 2 | `RHO_HDR_BYTES` @19, `RHO_MAGIC` @20 |
| [simple_reconstructor_pcg.f90:32](../../../src/main/volume/simple_reconstructor_pcg.f90#L32) | 10 | `PCG_OP_MATRIXFREE` @32, `PCG_OP_KERNEL` @33, `PCG_RAW_ACCUM_VERSION` @36, `PCG_RAW_PROV_LEN` @37, `PCG_SUPPORT_DIV_MIN` @38, `PCG_RHO_SAMPLED_FRAC` @39, `PCG_RAW_ACCUM_MAGIC` @40, `PCG_LAMBDA` @42, `PCG_WEIGHT_THRESHOLD` @44, `PCG_START_MAX_REL_RESID` @47 |
| [simple_refine3D_provenance.f90:28](../../../src/main/volume/simple_refine3D_provenance.f90#L28) | 5 | `VERSION` @28, `MAX_FILES` @28, `REFSET_OK` @29, `REFSET_LEGACY` @29, `REFSET_INVALID` @29 |
| [simple_symanalyzer.f90:14](../../../src/main/volume/simple_symanalyzer.f90#L14) | 7 | `DEBUG_HERE` @14, `WRITE_VOLUMES` @15, `SHSRCH_HWDTH` @16, `ZSCORE_LBOUND_PEAK` @17, `ZSCORE_LBOUND_C1` @18, `MAXKMIT` @19, `NCLS` @20 |
| [simple_trail_chain_manifest.f90:16](../../../src/main/volume/simple_trail_chain_manifest.f90#L16) | 4 | `TRAIL_MANIFEST_OK` @16, `TRAIL_MANIFEST_MISSING` @17, `TRAIL_MANIFEST_UNREADABLE` @18, `MANIFEST_VERSION` @19 |
| [simple_volpair_metrics.f90:19](../../../src/main/volume/simple_volpair_metrics.f90#L19) | 1 | `SMPD_TOL` @19 |
| [simple_volpft_corrcalc.f90:44](../../../src/main/volume/simple_volpft_corrcalc.f90#L44) | 1 | `i_startvec` @44 |
| [simple_volpft_symsrch.f90:13](../../../src/main/volume/simple_volpft_symsrch.f90#L13) | 2 | `NPROJ` @13, `NBEST` @14 |
| [simple_wave_multislice.f90:54](../../../src/main/volume/simple_wave_multislice.f90#L54) | 2 | `ZS_POINT` @54, `ZS_DIR` @54 |
| [simple_wave_projector.f90:96](../../../src/main/volume/simple_wave_projector.f90#L96) | 6 | `PASS_LOSS_GRAD` @96, `PASS_ADJOINT` @96, `PASS_CURVATURE` @96, `PASS_PREDICT` @96, `PASS_LINEAR` @96, `PASS_GAIN` @97 |

### src/utils

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_aff_prop.f90:10](../../../src/utils/clustering/simple_aff_prop.f90#L10) | 1 | `EXEMPLAR_STABLE_CHECKS` @10 |
| [simple_gmm.f90:23](../../../src/utils/clustering/simple_gmm.f90#L23) | 2 | `GMM_MAX_RESPAWN` @23, `GMM_MERGE_D2` @24 |
| [simple_kmeans.f90:15](../../../src/utils/clustering/simple_kmeans.f90#L15) | 1 | `KMEANS_MAXITS` @15 |
| [simple_xd_gmm.f90:23](../../../src/utils/clustering/simple_xd_gmm.f90#L23) | 1 | `XD_RIDGE_REL` @23 |
| [simple_http_post.f90:27](../../../src/utils/comm/simple_http_post.f90#L27) | 1 | `VERIFY_PEER` @27 |
| [simple_ipc_tcp_socket_client.f90:17](../../../src/utils/comm/simple_ipc_tcp_socket_client.f90#L17) | 5 | `TCP_MAX_RETRIES` @17, `TCP_TIMEOUT_MS` @18, `DEBUG` @19, `SO_SNDTIMEO` @38, `SO_RCVTIMEO` @39 |
| [simple_ipc_tcp_socket_helpers.f90:7](../../../src/utils/comm/simple_ipc_tcp_socket_helpers.f90#L7) | 8 | `TCP_MAX_MSG` @7, `POLLIN` @10, `POLLOUT` @11, `POLLERR` @12, `POLLHUP` @13, `POLLNVAL` @14, `MSG_PEEK` @30, `MSG_DONTWAIT` @31 |
| [simple_ipc_tcp_socket_server.f90:22](../../../src/utils/comm/simple_ipc_tcp_socket_server.f90#L22) | 3 | `TCP_PORT_MIN` @22, `TCP_PORT_MAX` @23, `LISTEN_BACKLOG` @46 |
| [simple_butterworth.f90:11](../../../src/utils/filter/simple_butterworth.f90#L11) | 1 | `BW_ORDER` @11 |
| [simple_dose_curve.f90:18](../../../src/utils/filter/simple_dose_curve.f90#L18) | 4 | `VERSION_LINE` @18, `NHEADER` @19, `KV_ATOL` @20, `DOSE_RTOL` @21 |
| [simple_dose_series.f90:20](../../../src/utils/filter/simple_dose_series.f90#L20) | 12 | `DOSE_WEIGHTS_MODEL` @20, `DOSE_WEIGHTS_CONTROL` @20, `DOSE_WEIGHTS_FLAT` @20, `N_SCHEMES` @21, `MIN_BAND_SHELLS` @23, `N_ALT_PASSES` @24, `N_LM_ITER` @25, `N_LM_TRIES` @26, `N_BANDS` @27, `OWN_FSC_FLOOR` @28, `SE_RATIO` @29, `C_DETECT` @30 |
| [simple_estimate_ssnr.f90:23](../../../src/utils/filter/simple_estimate_ssnr.f90#L23) | 1 | `LPSTAGES_STAGE1_LP_FLOOR` @23 |
| [simple_gui_metadata_histogram.f90:17](../../../src/utils/gui/metadata/simple_gui_metadata_histogram.f90#L17) | 1 | `MAX_HISTOGRAM_BINS` @17 |
| [simple_gui_metadata_micrograph.f90:18](../../../src/utils/gui/metadata/simple_gui_metadata_micrograph.f90#L18) | 1 | `MAX_MIC_COORDINATES` @18 |
| [simple_gui_metadata_optics_group.f90:12](../../../src/utils/gui/metadata/simple_gui_metadata_optics_group.f90#L12) | 1 | `MAX_OPTICS_SHIFTS` @12 |
| [simple_gui_metadata_timeplot.f90:17](../../../src/utils/gui/metadata/simple_gui_metadata_timeplot.f90#L17) | 1 | `MAX_TIMEPLOT_POINTS` @17 |
| [simple_gui_metadata_vol3D.f90:21](../../../src/utils/gui/metadata/simple_gui_metadata_vol3D.f90#L21) | 4 | `MAX_FSC_VOL3D` @21, `ORIDIST_BINWIDTH` @22, `ORIDIST_NBINS_X` @23, `ORIDIST_NBINS_Y` @24 |
| [simple_gui_metadata_stream_particle_sieving.f90:19](../../../src/utils/gui/metadata/stream/simple_gui_metadata_stream_particle_sieving.f90#L19) | 1 | `MAX_SIEVE_SELECTION` @19 |
| [simple_gui_metadata_stream_solve3D_multistate.f90:21](../../../src/utils/gui/metadata/stream/simple_gui_metadata_stream_solve3D_multistate.f90#L21) | 1 | `MAX_STATES_SOLVE3D_MULTISTATE` @21 |
| [simple_gui_metadata_stream_update.f90:18](../../../src/utils/gui/metadata/stream/simple_gui_metadata_stream_update.f90#L18) | 3 | `MAX_PICKREFS_SELECTION` @18, `MAX_SNAPSHOT2D_SELECTION` @19, `MAX_SNAPSHOT_FNAME_LEN` @20 |
| [simple_gui_assembler.f90:35](../../../src/utils/gui/simple_gui_assembler.f90#L35) | 15 | `GUI_STAGE_STATUS_UNKNOWN` @35, `GUI_STAGE_STATUS_RUNNING` @36, `GUI_STAGE_STATUS_FAILED` @37, `GUI_STAGE_STATUS_FINISHED` @38, `GUI_STAGE_STATUS_SKIPPED` @39, `SECTION_PREPROCESS` @52, `SECTION_OPTICS_ASSIGNMENT` @53, `SECTION_INITIAL_PICKING` @54, `SECTION_REFERENCE_PICKING` @55, `SECTION_INITIAL_ANALYSIS` @56, `SECTION_PARTICLE_SIEVING` @57, `SECTION_POOL2D` @58, `SECTION_SOLVE3D_MULTISTATE` @59, `SECTION_PROJECT` @60, `NSECTIONS` @61 |
| [simple_gui_project_builder.f90:45](../../../src/utils/gui/simple_gui_project_builder.f90#L45) | 1 | `UNSET` @45 |
| [simple_math.f90:23](../../../src/utils/math/simple_math.f90#L23) | 1 | `warn` @23 |
| [simple_rnd.f90:76](../../../src/utils/math/simple_rnd.f90#L76) | 22 | `RND_KEYED_VERSION` @76, `RND_KEY_NLANES` @77, `RND_MASK32` @78, `RND_MIX_C1` @79, `RND_MIX_C2` @80, `RND_LANE_A0` @81, `RND_LANE_B0` @82, `RND_TWO53` @83, `RND_FAM_SAMPLE` @85, `RND_FAM_INPL` @86, `RND_FAM_DISK` @87, `RND_FAM_CAND` @88, `RND_FAM_ASSIGN` @89, `RND_FAM_DIAG` @90, `RND_FAM_INPLQ` @91, `RND_PROG_REFINE3D` @93, `RND_PROG_REFINE2D` @93, `RND_PASS_ALIGN` @94, `RND_PASS_CONT` @95, `RND_PASS_FINAL` @96, `RND_PASS_DIAG` @97, `RND_SCOPE_STRIDE` @98 |
| [simple_stat.f90:75](../../../src/utils/math/simple_stat.f90#L75) | 1 | `NNET_CONST` @75 |
| [simple_persistent_worker_server.f90:30](../../../src/utils/persistent_worker/simple_persistent_worker_server.f90#L30) | 11 | `WORKER_STATUS_IDLE` @30, `WORKER_STATUS_ERROR` @31, `TASK_QUEUE_SIZE` @32, `TCP_BUFSZ` @33, `KILL_WAIT_TIME_US` @34, `KILL_POLL_TIME_US` @35, `POLL_TIMEOUT_MS` @36, `SWEEP_PERIOD` @37, `SCALE_DOWN_COOLDOWN_S` @38, `SCALE_UP_COOLDOWN_S` @39, `DEBUG` @40 |
| [simple_qsys_async_job.f90:41](../../../src/utils/qsys/simple_qsys_async_job.f90#L41) | 4 | `ASYNC_JOB_IDLE` @41, `ASYNC_JOB_RUNNING` @42, `ASYNC_JOB_DONE` @43, `ASYNC_JOB_FAILED` @44 |
| [simple_qsys_coarray.f90:11](../../../src/utils/qsys/simple_qsys_coarray.f90#L11) | 1 | `MAXENVITEMS` @11 |
| [simple_qsys_ctrl.f90:34](../../../src/utils/qsys/simple_qsys_ctrl.f90#L34) | 2 | `SHORTTIME` @34, `LOCAL_POLL_US` @37 |
| [simple_qsys_funs.f90:12](../../../src/utils/qsys/simple_qsys_funs.f90#L12) | 1 | `SHORTTIME` @12 |
| [simple_qsys_job_record.f90:50](../../../src/utils/qsys/simple_qsys_job_record.f90#L50) | 5 | `JOB_UNKNOWN` @50, `JOB_ALIVE` @51, `JOB_GONE` @52, `JOB_LIVENESS_S` @53, `JOB_LOST_EXIT_CODE` @54 |
| [simple_qsys_local.f90:10](../../../src/utils/qsys/simple_qsys_local.f90#L10) | 1 | `MAXENVITEMS` @10 |
| [simple_qsys_lsf.f90:10](../../../src/utils/qsys/simple_qsys_lsf.f90#L10) | 1 | `MAXENVITEMS` @10 |
| [simple_qsys_pbs.f90:10](../../../src/utils/qsys/simple_qsys_pbs.f90#L10) | 1 | `MAXENVITEMS` @10 |
| [simple_qsys_persistent_worker.f90:17](../../../src/utils/qsys/simple_qsys_persistent_worker.f90#L17) | 1 | `MAXENVITEMS` @17 |
| [simple_qsys_sge.f90:10](../../../src/utils/qsys/simple_qsys_sge.f90#L10) | 1 | `MAXENVITEMS` @10 |
| [simple_qsys_slurm.f90:10](../../../src/utils/qsys/simple_qsys_slurm.f90#L10) | 1 | `MAXENVITEMS` @10 |
| [simple_forked_process.f90:28](../../../src/utils/simple_forked_process.f90#L28) | 5 | `FORK_STATUS_FAILED` @28, `FORK_STATUS_RUNNING` @29, `FORK_STATUS_STOPPED` @30, `FORK_STATUS_SKIPPED` @31, `FORK_POLL_TIME` @32 |
| [simple_gpu_utils.f90:9](../../../src/utils/simple_gpu_utils.f90#L9) | 9 | `CUFFT_R2C` @9, `CUFFT_R2C_MANY` @10, `CUFFT_C2R` @11, `CUDA_SUCCESS` @12, `CUFFT_SUCCESS` @13, `CUBLAS_STATUS_SUCCESS` @14, `CUBLAS_OP_N` @15, `CUBLAS_OP_T` @16, `CUBLAS_OP_C` @17 |
| [simple_magic_boxes.f90:13](../../../src/utils/simple_magic_boxes.f90#L13) | 3 | `NSZS` @13, `boxsizes` @14, `MIN_PFTSZ` @25 |
| [simple_memory_monitor.f90:15](../../../src/utils/simple_memory_monitor.f90#L15) | 1 | `DEFAULT_REPORT_INTERVAL_SECONDS` @15 |
| [simple_timer.f90:36](../../../src/utils/simple_timer.f90#L36) | 3 | `timer_int_kind` @36,38, `MAX_TOKENS` @45, `MAX_TOKEN_CHARSIZE` @46 |
| [simple_chash.f90:21](../../../src/utils/structs/simple_chash.f90#L21) | 1 | `BUFFSZ_DEFAULT` @21 |
| [simple_hash.f90:13](../../../src/utils/structs/simple_hash.f90#L13) | 1 | `BUFFSZ_DEFAULT` @13 |

## Appendix B. All production commander procedure/block constants

These 161 declaration occurrences are additional to Appendix A's 51 commander constants. Repeated names are distinct local declarations; no recommendation to merge them is implied.

| File | Count | Constants and declaration lines |
| --- | ---: | --- |
| [simple_commanders_atoms.f90:412](../../../src/main/commanders/simple/simple_commanders_atoms.f90#L412) | 1 | `SPECIES_KEYS` @412 |
| [simple_commanders_cavgs.f90:101](../../../src/main/commanders/simple/simple_commanders_cavgs.f90#L101) | 16 | `DEBUG` @101, `SCORE_THRES_INCL` @102, `NCLUST_MAX` @103, `QUALITY_MODE_APPLY` @430, `QUALITY_MODE_ANALYZE` @431, `QUALITY_MODE_LEARN` @432, `QUALITY_MODE_EVALUATE` @433, `QUALITY_MODE_PROMOTE` @434, `STATE_SCORE_FNAME` @626, `LP_BIN` @859, `DEBUG` @860, `NCLUST_MAX` @1040, `LP_DEFAULT` @1305, `TRS_DEFAULT` @1306, `NSPACE_DEFAULT` @1307, `LP_BIN` @1368 |
| [simple_commanders_class_student.f90:88](../../../src/main/commanders/simple/simple_commanders_class_student.f90#L88) | 5 | `NPHASES` @88, `NARM` @89, `PHASE_NAMES` @90, `PHASE_ORDER` @96, `RULE` @420 |
| [simple_commanders_flex_pca.f90:145](../../../src/main/commanders/simple/simple_commanders_flex_pca.f90#L145) | 1 | `NATIVE_MEAN_FNAME` @145 |
| [simple_commanders_imgops.f90:122](../../../src/main/commanders/simple/simple_commanders_imgops.f90#L122) | 1 | `SIGMA_DEFAULT` @122 |
| [simple_commanders_imgproc.f90:274](../../../src/main/commanders/simple/simple_commanders_imgproc.f90#L274) | 3 | `FILT` @274, `BIN` @275, `MASKED` @276 |
| [simple_commanders_mic_subset.f90:142](../../../src/main/commanders/simple/simple_commanders_mic_subset.f90#L142) | 1 | `keys` @142 |
| [simple_commanders_ori.f90:184](../../../src/main/commanders/simple/simple_commanders_ori.f90#L184) | 2 | `hlen` @184, `NSPACE_REDUCED` @184 |
| [simple_commanders_pose_student.f90:83](../../../src/main/commanders/simple/simple_commanders_pose_student.f90#L83) | 5 | `NPHASES` @83, `NRECS` @84, `PHASE_NAMES` @85, `REC_NAMES` @89, `REC_HELP` @92 |
| [simple_commanders_project_core.f90:267](../../../src/main/commanders/simple/simple_commanders_project_core.f90#L267) | 8 | `NREMAP_SCOPES` @267, `REMAP_SCOPES` @268, `REMAP_OLD_KEYS` @270, `REMAP_NEW_KEYS` @272, `CONSENSUS_REPORT` @413, `MAX_EXACT_NSTATES` @414, `LP_BIN` @733, `DEFTAB` @773 |
| [simple_commanders_project_ptcl.f90:805](../../../src/main/commanders/simple/simple_commanders_project_ptcl.f90#L805) | 5 | `STKFBODY` @805, `BUFBYTES` @806, `REPORT_SECS` @807, `RATE_WEIGHT` @808, `MAX_NCUNITS` @977 |
| [simple_commanders_refine3D.f90:330](../../../src/main/commanders/simple/simple_commanders_refine3D.f90#L330) | 32 | `DEBUG` @330, `MINBOX` @331, `MAXITS_REFINE3D_CAP` @512, `TARGET_UPDATES_PER_PARTICLE_REFINE3D` @513, `REGPASS_NSPACE` @649, `STATE_OVERLAP_NEIGH_REFINE3D_STATES` @773, `LPSTART_REFINE3D_STATES` @774, `LPSTOP_REFINE3D_STATES` @775, `WORKFLOW_LABEL` @776, `STATE_OVERLAP_TARGET` @860, `LPSTART_DEFAULT` @861, `LPSTOP_DEFAULT` @862, `MIN_STATE_FRAC` @863, `SMPD_TOL` @864, `MIN_SUBAVGS` @865, `FLEX_MIN_NEFF` @866, `FLEX_MIN_FIT_ROWS` @867, `WORKFLOW_LABEL` @868, `WORK_PROJFILE` @869, `WORK_SIGMA2_STATE` @870, `MEAN_VOL_FNAME` @871, `NSAMPLE_PER_STATE_CLASSIFY3D_REFS` @1305, `NSAMPLE_CLASSIFY3D_REFS_CAP` @1306, `NSPACE` @1307, `NSPACE_SUB` @1308, `MINITS_CLASSIFY3D_REFS` @1309, `MAXITS_CLASSIFY3D_REFS` @1310, `TARGET_UPDATES_PER_PARTICLE` @1311, `LPSTART_CLASSIFY3D_REFS` @1312, `LPSTOP_CLASSIFY3D_REFS` @1313, `WORKFLOW_LABEL` @1314, `STEP_IT` @1571 |
| [simple_commanders_resolest.f90:81](../../../src/main/commanders/simple/simple_commanders_resolest.f90#L81) | 6 | `PLOT_RES` @81, `TAB_CHAR` @377, `LPSTART_LB` @791, `LPSTART_DEFAULT` @791, `LPSTOP_BOUNDS` @792, `LPSTART_BOUNDS` @793 |
| [simple_commanders_sim.f90:261](../../../src/main/commanders/simple/simple_commanders_sim.f90#L261) | 3 | `BATCH` @261, `NPIX_CUTOFF` @643, `BATCH` @1162 |
| [simple_commanders_solve2D.f90:889](../../../src/main/commanders/simple/simple_commanders_solve2D.f90#L889) | 3 | `DIR_PROJS` @889, `WAITTIME` @890, `TARGET_NCLS_PER_CHUNK` @891 |
| [simple_commanders_solve3D.f90:924](../../../src/main/commanders/simple/simple_commanders_solve3D.f90#L924) | 3 | `ADDON_DIAG_DIR` @924, `FROZEN_COPY_DIR` @925, `INI3D_DIR` @2307 |
| [simple_commanders_validate.f90:32](../../../src/main/commanders/simple/simple_commanders_validate.f90#L32) | 9 | `PROJFILE_MINI_STREAM` @32, `NCLS_MIN` @33, `NCLS_MAX` @33, `LPSTOP` @34, `DIR_THUMBS` @134, `PROJFILE_CHECK_REFPICK` @135, `NCLS_MIN` @136, `NCLS_MAX` @136, `LPSTOP` @137 |
| [simple_commanders_volops.f90:173](../../../src/main/commanders/simple/simple_commanders_volops.f90#L173) | 6 | `DOCK_HP` @173, `EDGE` @175, `SYMSHTAB` @831, `SYMAXISTAB` @832, `MAXBOX` @916, `MAXBOX` @1005 |
| [single_commanders_experimental.f90:294](../../../src/main/commanders/single/single_commanders_experimental.f90#L294) | 1 | `RAD_LB` @294 |
| [single_commanders_nano2D.f90:34](../../../src/main/commanders/single/single_commanders_nano2D.f90#L34) | 3 | `NCLS_CEN_NANO` @34, `STARTVOL` @162, `LP_EST_DIAM` @163 |
| [single_commanders_nano3D.f90:67](../../../src/main/commanders/single/single_commanders_nano3D.f90#L67) | 5 | `SPLITTED` @67, `FINAL_MAPS` @68, `TAG` @69, `LP_LIST` @425, `HP_LIM` @426 |
| [single_commanders_tomo.f90:77](../../../src/main/commanders/single/single_commanders_tomo.f90#L77) | 12 | `TXT_READ` @77, `MASK_FNAME` @130, `FINAL_FBODY` @131, `LPSTART_DEFAULT` @132, `LPSTOP_DEFAULT` @133, `AMSKLP_DEFAULT` @134, `TRS_DEFAULT` @135, `ATHRES_DEFAULT` @136, `MASK_GROW_A` @137, `NSTAGES_DEFAULT` @138, `MAXITS_DEFAULT` @139, `MAXITS_PCG_DEFAULT` @140 |
| [single_commanders_trajectory.f90:383](../../../src/main/commanders/single/single_commanders_trajectory.f90#L383) | 2 | `SEED` @383, `BATCH` @384 |
| [single_commanders_trajectory_solve3D.f90:40](../../../src/main/commanders/single/single_commanders_trajectory_solve3D.f90#L40) | 25 | `NSTAGES` @40, `S_IMPORT` @41, `S_CENTER2D` @41, `S_DIAM` @41, `S_STARTVOL` @41, `S_REFINE2D` @41, `S_GATE1` @42, `S_GATE2` @42, `S_EXTRACT` @42, `S_DENOISE` @42, `S_SWAP` @42, `S_AUTOREFINE3D` @43, `STAGE_PRGS` @44, `RECORDS` @48, `KEPT` @49, `DENOISED` @50, `STARTVOL` @51, `GATE_TABLE` @52, `GATE_REJ` @53, `GATE_KEPT` @54, `GATE_REPROJ` @55, `NCHUNKS` @56, `NKEPT_PAIRS` @57, `LP_EST_DIAM` @58, `MSK_FAC` @59 |
| [single_commanders_tseries.f90:230](../../../src/main/commanders/single/single_commanders_tseries.f90#L230) | 1 | `LAM_TV` @230 |
| [simple_commanders_stream_p00_master.f90:312](../../../src/main/commanders/stream/simple_commanders_stream_p00_master.f90#L312) | 2 | `POLL_US` @312, `POLL_US` @424 |
