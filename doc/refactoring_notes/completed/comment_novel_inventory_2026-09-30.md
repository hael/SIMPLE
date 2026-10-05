# Comment novels: inventory, September 2026

Line numbers refer to `master` at commit 714a15a99. Scope: `src/` (excluding `src/extlibs`) and
`production/`: 670 Fortran files, about 437,000 lines. This is a read-only audit; no source file
was changed. Replacement drafts for every audited block are in the companion note
`comment_novel_inventory_2026-09-30_drafts.md`.

## Status (2026-10-01)

Two passes are done, both as an uncommitted working-tree diff.

What they covered:
- Pass 1:
  - all 85 audited records (table E1 and section D);
  - the C1 switch mentions and the C2 paths;
  - the section G docs;
  - the dead code named in H.4.
- Pass 2:
  - the remaining blocks of 100 words or more, each checked against the code before it was rewritten;
  - all 33 remaining templated headers;
  - every dated comment line;
  - the further flex and NU dead code the pass-1 agents had found.

Results:
- no comment block of 100 words or more remains;
- no templated header remains;
- words in blocks of 80 words or more fell from 38,750 to 6,277;
- dated comment lines fell from 174 to 5 (4 author attributions and a seed literal);
- C1 and C2 are both empty.

Still open:
- the flex_pca default ruling (B.2);
- the refine3D UI offering `snhc`/`shc_neigh` modes that the matcher rejects;
- blocks under 100 words (not in scope).

## Verdict

Most comments in SIMPLE are fine. Of 14,570 prose comment blocks, 98% are under 80 words. The
problem is concentrated. 298 blocks of 80 words or more hold 38,750 words, a fifth of all comment
prose, and 65 files carry at least one block of 150 words or more. flex_pca alone accounts for
11,400 of those words, more than commanders, strategies and volume combined.

Long comments are usually also inaccurate. Each of the 78 blocks of 150 words or more was checked
against the code:

- 29 make a claim that the code contradicts (WRONG);
- 23 more are out of date (STALE): moved paths, renamed tests, removed switches, wrong dates;
- 26 are accurate (OK).

Section A gives four reasons why long comments rot.

Proposed target, in the spirit of `scripts/check_descr.py` (whose `!@descr:` docstring already
calls the tag "the subject line of a commit message"):

- File header: the `!@descr:` line plus at most 5 lines.
- Procedure header: at most 3 lines.
- In-body comments: 1-2 lines that say *why*. Leave out *what* (the code shows it) and *when* (git
  records it).
- Dates, names, measurements, "used to", and rejected alternatives go in `git log` or a
  `doc/**/..._decision_log.md`, not in code.
- A comment that restates a policy note should point to it instead.

## A. Why the long comments are wrong

1. **The lab notebook lives in the code.** 174 comment lines in 80 files carry a date or a name
   ("2026-09-19, Hans", "measured 2026-09-11", "until 2026-09-12"), and 58 of the 298 long blocks
   contain one. A dated sentence records what was true on that day. Every later change can make it
   wrong without anyone noticing.
   - Example: `simple_nu_filter.f90:10` dates the removal of the shell walk to 2026-09-18. Git
     shows the deletion in b9cc8e31b on 2026-09-16. Line 14 cites `ed36eb4c`, which touches only
     NICE; it names the build that was run, which a reader of the code cannot tell.
   - Densest files (dated or attributed comment lines): `simple_nu_filter.f90` (15), `simple_commanders_refine3D.f90` (9), `simple_nu_state_filter.f90` (9), `simple_commanders_volops.f90` (8), `simple_nu_filter_bank.f90` (8), `simple_abinitio_controller.f90` (6), `simple_image_tester.f90` (5), `simple_final_rec.f90` (4), `simple_commanders_postprocess_nu.f90` (4), `simple_reconstructor_pcg.f90` (4), `simple_commanders_test_class.f90` (3), `simple_commanders_flex_pca.f90` (3), `simple_flex_pca_strategy.f90` (3), `simple_nu_filter_evidence.f90` (3), `simple_nu_filter_stats.f90` (3), `simple_flex_pca_model.f90` (3), `simple_pcg_solvent_sidecar.f90` (3), `simple_pcg_halfset_tester.f90` (3).

2. **Some headers list what the compiler already knows.** 51 files use the templated
   `MODULE / PURPOSE / ENTRY POINT / INTERNAL SUBROUTINES / DEPENDENCIES` header (stream, GUI
   metadata, persistent worker, comm).
   - In 22 of them the DEPENDENCIES list no longer matches the `use` statements (table C3).
     `simple_stream_p03_initial_analysis.f90` omits 27 modules it uses.
   - The routine and field lists have drifted the same way. The p03 header documents
     `micimporter`, shape ranking, `process_selected_refs` and `restart_requested`, none of which
     runs. The `simple_gui_metadata_vol3D.f90` field list misses `cfar`, `i` and `i_max`.

3. **Switches and paths are still described after their removal.**
   - 26 `SIMPLE_COV_*` environment switches are named in comments (39
     mentions) but read nowhere in code, scripts or NICE (table C1). Some of these comments tell
     the user to set the switch. For example, `simple_flex_pca_planes.f90:16` says
     "SIMPLE_COV_RESIDENT_GB overrides it", and the same name appears in a runtime log message at
     line 62.
   - 14 citations of `doc/` files point to paths that do not exist. 7 of these documents have
     moved; the other 7 were never committed (table C2).
   - Two comments cite personal paths: `~/ribo_local/rank_criterion.py` and
     `/Users/elmlundho/model_cavgs_rejection/...`.
   - Six comments cite `file.f90:NNN` line numbers. At least two are already wrong.

4. **Some headers copy a policy note, and the two then drift apart.**
   - `simple_ptcl_cache.f90:1-53` (520 words, the largest non-vendored block) still documents
     its 3D consumers: refine3D, prob_tab and calc_3Drec. Both 3D entry points reject
     `cache=yes` with `THROW_HARD` (`simple_refine3D_strategy.f90:633`,
     `simple_commanders_abinitio.f90:823`).
   - `simple_commanders_flex_pca.f90:16-29` says its PCG defaults "MATCH the flex_pca UI
     declaration ... the cap is 20". The constant below it is 4. The UI declares 20, 1e-3 and
     gridding (`simple_ui_heterogeneity.f90:163-176`).
   - `simple_commanders_refine3D.f90:64` and `:285` say the registration pass is `refine=prob`.
     The code sets `greedy` (:354).

## B. Fix first

The 29 WRONG blocks are listed in table E1. The ones below are most likely to mislead someone who
changes the code. Each fix also makes the comment shorter.

1. `strategies/search/simple_ptcl_cache.f90:1-53`: the header documents 3D consumers, but the cache
   is 2D-only. The policy note has the same leftover (`doc/policies/2D/particle_cache_policy.md`
   :83, :133, :162).
2. `commanders/simple/simple_commanders_flex_pca.f90:16-29`: the defaults described in the
   comment and the defaults in the code have diverged from the UI. **This needs a ruling.** Should
   a UI launch and a command-line launch run the same estimator?
3. `flex/simple_flex_pca_pcg.f90:1-19` and `solve` (:1816-1906): the header says "every solve
   starts from ZERO", and line 1906 says "Opt-in until measured: SIMPLE_COV_PCG_WARM=1". Line 1907
   hard-wires `iwarm = 1`.
4. `commanders/simple/simple_commanders_refine3D.f90:64-74, 211-232, 285`: the registration pass is
   described as `prob` (the code uses `greedy`). The startup is described as "reconstruct ->
   build masks -> re-reconstruct", but the code runs one reconstruction.
5. `volume/simple_pcg_solvent_sidecar.f90:16-17`: "the FSC is solvent-flattened" contradicts lines
   21-22 of the same header and the code. The FSC comes from the prior-free pair.
6. `utils/persistent_worker/message/simple_persistent_worker_message_terminate.f90` and
   `production/simple_persistent_worker.f90`: both say the worker "finishes its current tasks".
   It actually cancels them (`simple_persistent_worker.f90:281-299`). The production header also
   says scripts are validated, but validation is bypassed (`safe_path = .true.`, :235).
7. `stream/simple_stream_p03_initial_analysis.f90:1-35`: the header describes a design that no
   longer runs (see A.2). The 2026-09-30 stream review (A2) already flagged this header.
8. `nu_filt/simple_nu_filter_envmask.f90:1-27` gives the wrong margin formula: the margin comes
   from the raw unaries, not `dmats_mask`. `commanders/simple/simple_commanders_postprocess_nu.f90:1-21`
   describes the retired v1 pipeline.
9. `ori/simple_oris_sampling.f90:78-95`: "mnew = N whenever M > 0 or n = N" fails for M = 0,
   n = N, u < 1 (:113-118). This is a contract comment on the fractional-update population rule.
10. `pftc/simple_pftc_shsrch_grad.f90:416-436`: "without irot_in, a discrete all-angle scan". In
    fact `irot_in` is required and no scan path exists.
11. `flex/simple_flex_pca_planes.f90:1-17` and `flex/simple_flex_pca_plane_cache.f90:1-19`
    describe an environment override and a selection-hash check that do not exist.

## C. Library-wide mechanical checks

### C1. Environment switches named in comments but read nowhere

| switch named in comments | mentions | first mention |
|---|---:|---|
| `SIMPLE_COV_PAIRED_MERGE` | 5 | `main/flex/simple_flex_pca_em_iter.f90:788` |
| `SIMPLE_COV_EM_DEFLATE` | 4 | `main/flex/simple_flex_pca_em_mstep.f90:259` |
| `SIMPLE_COV_XFSC_REG` | 4 | `main/flex/simple_flex_pca_em_mstep.f90:91` |
| `SIMPLE_COV_BASIS_MAX` | 3 | `main/flex/simple_flex_pca_em_fit.f90:211` |
| `SIMPLE_COV_GMM` | 2 | `main/flex/simple_flex_pca_weights.f90:301` |
| `SIMPLE_COV_CONTRAST` | 1 | `main/flex/simple_flex_pca_em_mean.f90:12` |
| `SIMPLE_COV_DECONV` | 1 | `main/flex/simple_flex_pca_model.f90:838` |
| `SIMPLE_COV_EM` | 1 | `main/flex/simple_flex_pca_em_fit.f90:126` |
| `SIMPLE_COV_EM_MIX` | 1 | `main/flex/simple_flex_pca_em_iter.f90:456` |
| `SIMPLE_COV_GMM_AUTO` | 1 | `main/flex/simple_flex_pca_weights.f90:316` |
| `SIMPLE_COV_HALF_CONTRAST` | 1 | `main/flex/simple_flex_pca_em_embed.f90:87` |
| `SIMPLE_COV_MASTER_NTHR` | 1 | `main/strategies/parallelization/simple_flex_pca_strategy.f90:229` |
| `SIMPLE_COV_MIN_NEFF` | 1 | `main/ui/simple/simple_ui_heterogeneity.f90:116` |
| `SIMPLE_COV_MIN_STATE` | 1 | `main/flex/simple_flex_pca_gmm.f90:328` |
| `SIMPLE_COV_PAIRED` | 1 | `main/flex/simple_flex_pca_em_iter.f90:630` |
| `SIMPLE_COV_PCG_WARM` | 1 | `main/flex/simple_flex_pca_pcg.f90:1906` |
| `SIMPLE_COV_POLAR_ESTEP` | 1 | `main/flex/simple_flex_pca_em_iter.f90:55` |
| `SIMPLE_COV_POLAR_EXACT` | 1 | `main/flex/simple_flex_pca_em_polar.f90:25` |
| `SIMPLE_COV_POLAR_NDIR` | 1 | `main/flex/simple_flex_pca_em_polar.f90:34` |
| `SIMPLE_COV_POLAR_OSAMP` | 1 | `main/flex/simple_flex_pca_em_iter.f90:78` |
| `SIMPLE_COV_POLAR_RHYB` | 1 | `main/flex/simple_flex_pca_em_iter.f90:77` |
| `SIMPLE_COV_PROBE_CONTRAST` | 1 | `main/flex/simple_flex_pca_em_iter.f90:504` |
| `SIMPLE_COV_PROBE_CONV` | 1 | `main/flex/simple_flex_pca_em_iter.f90:52` |
| `SIMPLE_COV_PROBE_MLSCALE` | 1 | `main/flex/simple_flex_pca_em_iter.f90:505` |
| `SIMPLE_COV_RESIDENT` | 1 | `main/flex/simple_flex_pca_planes.f90:17` |
| `SIMPLE_COV_XFSC_` | 1 | `main/flex/simple_flex_pca_em_iter.f90:115` |

### C2. Doc paths cited in comments that do not exist

| comment | cites | status |
|---|---|---|
| `main/commanders/test/simple_commanders_test_highlevel.f90:1160` | `doc/policies/reconstruct3D_pcg_policy.md` | moved to `doc/policies/3D/reconstruct3D_pcg_policy.md` |
| `main/interp/simple_gridding.f90:71` | `doc/implementation_notes/drop_legacy_box_division.md` | moved to `doc/implementation_notes/completed/drop_legacy_box_division.md` |
| `main/cavg_quality/simple_cavg_quality_feats.f90:52` | `doc/microchunk_and_rejection/model_cavgs_rejection.md` | moved to `doc/policies/sieving_and_rejection/model_cavgs_rejection.md` |
| `main/image/simple_projector.f90:47` | `doc/implementation_notes/drop_legacy_box_division.md` | moved to `doc/implementation_notes/completed/drop_legacy_box_division.md` |
| `main/image/simple_image_tester.f90:7` | `doc/refactoring_notes/planned/image_rmat_padding_encapsulation.md` | moved to `doc/refactoring_notes/completed/image_rmat_padding_encapsulation.md` |
| `main/pftc/simple_polarft_corr.f90:1900` | `doc/implementation_notes/drop_legacy_box_division.md` | moved to `doc/implementation_notes/completed/drop_legacy_box_division.md` |
| `main/sigma2/simple_euclid_sigma2.f90:21` | `doc/implementation_notes/drop_legacy_box_division.md` | moved to `doc/implementation_notes/completed/drop_legacy_box_division.md` |
| `main/flex/simple_flex_pca_crossfsc.f90:2` | `doc/for_developers/ideas/flex_pca_crossfsc_shrinkage_marching_spec.md` | not in repo |
| `main/flex/simple_flex_pca_rec3D.f90:81` | `doc/refactoring_notes/flex_pca_branch_reconciliation_2026_09_15.md` | not in repo |
| `main/flex/simple_flex_pca_em.f90:53` | `doc/policies/flex_pca_policy.md` | not in repo |
| `main/flex/simple_flex_pca_em.f90:273` | `doc/for_developers/ideas/flex_pca_crossfsc_shrinkage_marching_spec.md` | not in repo |
| `main/flex/simple_flex_pca_targets.f90:24` | `doc/implementation_notes/flex_pca_state_placement_measurements.md` | not in repo |
| `main/flex/simple_flex_pca_gmm.f90:23` | `doc/implementation_notes/flex_pca_state_placement_measurements.md` | not in repo |
| `main/volume/simple_reconstructor_pcg.f90:3` | `doc/policies/reconstruct3D_pcg_policy.md` | moved to `doc/policies/3D/reconstruct3D_pcg_policy.md` |

`src/main/flex/README.md:20` also cites `doc/policies/flex_pca_policy.md`.

### C3. Templated headers whose DEPENDENCIES list has drifted

| header | drift |
|---|---|
| `main/stream/simple_stream_p01_preprocess_new.f90` | omits used modules: simple_gui_metadata_utils, simple_image, simple_motion_gain_analysis, simple_motion_gain_helpers |
| `main/stream/simple_stream_p03_initial_analysis.f90` | omits used modules: 27 modules |
| `main/stream/simple_stream_p05_sieve_cavgs_new.f90` | omits used modules: simple_fileio, simple_gui_metadata_api, simple_stream_state |
| `main/stream/simple_stream_p06_pool2D_new.f90` | omits used modules: simple_gui_metadata_utils |
| `main/stream/simple_stream_p07_abinitio3D_multistate.f90` | omits used modules: simple_commanders_cavgs, simple_gui_utils, simple_imghead, simple_qsys_env |
| `utils/comm/simple_http_post.f90` | omits used modules: simple_error, simple_string |
| `utils/gui/metadata/simple_gui_metadata_optics_group.f90` | lists modules not used: simple_core_module_api; omits used modules: simple_error |
| `utils/gui/metadata/simple_gui_metadata_project.f90` | lists modules not used: simple_eer_factory; omits used modules: 8 modules |
| `utils/gui/metadata/simple_gui_metadata_timeplot.f90` | lists modules not used: simple_core_module_api; omits used modules: simple_defs, simple_error, simple_string |
| `utils/gui/metadata/stream/simple_gui_metadata_stream_abinitio3D_multistate.f90` | omits used modules: simple_error |
| `utils/gui/metadata/stream/simple_gui_metadata_stream_opening2D.f90` | omits used modules: simple_error |
| `utils/gui/metadata/stream/simple_gui_metadata_stream_particle_sieving.f90` | omits used modules: simple_error |
| `utils/gui/metadata/stream/simple_gui_metadata_stream_picking.f90` | omits used modules: simple_error |
| `utils/gui/metadata/stream/simple_gui_metadata_stream_pool2D.f90` | omits used modules: simple_error |
| `utils/gui/metadata/stream/simple_gui_metadata_stream_pool2D_snapshot.f90` | omits used modules: simple_error, simple_string |
| `utils/gui/metadata/stream/simple_gui_metadata_stream_update.f90` | omits used modules: simple_defs, simple_error, simple_string |
| `utils/gui/simple_gui_assembler_tester.f90` | lists modules not used: simple_syslib |
| `utils/persistent_worker/message/simple_persistent_worker_message_status.f90` | omits used modules: simple_defs |
| `utils/persistent_worker/message/simple_persistent_worker_message_task.f90` | omits used modules: simple_defs |
| `utils/persistent_worker/message/simple_persistent_worker_message_terminate.f90` | omits used modules: simple_defs |
| `utils/persistent_worker/simple_persistent_worker_server.f90` | omits used modules: simple_ipc_tcp_socket_client, simple_ipc_tcp_socket_helpers |
| `utils/simple_forked_process.f90` | omits used modules: simple_defs, simple_error, simple_fileio, simple_memory_monitor |

The DEPENDENCIES section repeats the `use` block a few lines below it. Removing it costs nothing.

## D. Procedures where comments rival the code

These procedures have at least 200 comment words and at least 2 comment words per code line.
Declarations count as code. Examples:

- `fit_stage_config` has 462 comment words over 47 code lines. It documents environment switches
  it does not read.
- `build_covariance_eigenbasis` has 285 words over 20 lines. It documents rounds that no longer
  exist, and the routine looks unreachable: every non-compose master run goes to
  `run_flex_pca_paired` (`simple_flex_pca_model.f90:228`).

| comment words per code line | comment words | code lines | procedure |
|---:|---:|---:|---|
| 14.2 | 285 | 20 | `main/flex/simple_flex_pca_em_fit.f90:13` `build_covariance_eigenbasis` |
| 9.8 | 462 | 47 | `main/flex/simple_flex_pca_em_iter.f90:435` `fit_stage_config` |
| 5.2 | 235 | 45 | `main/strategies/parallelization/simple_flex_pca_strategy.f90:183` `master_initialize` |
| 5.1 | 345 | 68 | `main/flex/simple_flex_pca_pcg.f90:1076` `finalize` |
| 5.0 | 347 | 69 | `main/flex/simple_flex_pca_em_iter.f90:782` `probe_subspace_paired` |
| 5.0 | 281 | 56 | `main/nu_filt/simple_nu_filter_bank.f90:274` `optimize_nu_cutoff_finds` |
| 4.9 | 1408 | 285 | `main/flex/simple_flex_pca_em_iter.f90:23` `probe_subspace_iteration` |
| 4.6 | 404 | 88 | `main/flex/simple_flex_pca_pcg.f90:1824` `solve` |
| 4.2 | 416 | 100 | `main/params/simple_parameters_core.f90:8` `init_dynamic_defaults` |
| 4.2 | 299 | 72 | `main/commanders/simple/simple_commanders_rec_distr.f90:140` `blend_trailing_accumulators` |
| 4.1 | 566 | 138 | `main/abinitio/simple_abinitio_controller.f90:522` `emit_refine3D_stage_cfg` |
| 4.0 | 250 | 63 | `main/commanders/simple/simple_commanders_rec_distr.f90:44` `restore_state_from_parts` |
| 3.7 | 451 | 121 | `main/commanders/simple/simple_commanders_refine3D.f90:2096` `exec_bootstrap_rec3D` |
| 3.5 | 619 | 175 | `main/commanders/simple/simple_commanders_volops.f90:292` `postprocess_volume_from_files` |
| 3.5 | 235 | 67 | `main/image/simple_image_calc.f90:1919` `nu_objective` |
| 3.4 | 201 | 59 | `main/class/simple_classaverager_restore.f90:190` `cavger_update_sums` |
| 3.3 | 216 | 66 | `main/nano/simple_nanoparticle.f90:1988` `write_np_stats` |
| 3.2 | 452 | 141 | `utils/math/simple_testfuns.f90:26` `get_testfun` |
| 3.2 | 538 | 168 | `main/strategies/search/simple_ptcl_cache.f90:597` `ptcl_cache_ensure` |
| 3.0 | 1238 | 419 | `main/flex/simple_flex_pca_model.f90:82` `run_flex_pca` |
| 2.9 | 299 | 102 | `main/flex/simple_flex_pca_pcg.f90:1300` `apply_operator` |
| 2.9 | 821 | 282 | `main/commanders/simple/simple_commanders_refine3D.f90:49` `exec_refine3D_auto` |
| 2.9 | 1043 | 364 | `main/flex/simple_flex_pca_em_mstep.f90:28` `fit_iter_finish` |
| 2.7 | 201 | 74 | `fileio/simple_imghead.f90:724` `setMinimal` |
| 2.7 | 244 | 91 | `main/ui/simple/simple_ui_refine3D_pose_cont.f90:16` `construct_refine3D_pose_cont_program` |
| 2.7 | 923 | 347 | `main/flex/simple_flex_pca_em_embed.f90:15` `embed_latents_with_contrast` |
| 2.7 | 313 | 118 | `main/pftc/simple_pftc_shsrch_grad.f90:437` `grad_shsrch_minimize_joint` |
| 2.6 | 246 | 93 | `main/flex/simple_flex_pca_em_iter.f90:638` `run_flex_pca_paired` |
| 2.5 | 830 | 331 | `main/flex/simple_flex_pca_merge.f90:270` `two_gate_state_merge` |
| 2.5 | 209 | 84 | `main/flex/simple_flex_pca_em_basis.f90:100` `orthonormalize_representatives` |
| 2.4 | 679 | 285 | `main/flex/simple_flex_pca_weights.f90:27` `build_covariance_state_weights` |
| 2.1 | 303 | 141 | `main/strategies/search/simple_matcher_refvol_utils.f90:238` `read_mask_filter_refvols` |
| 2.0 | 294 | 147 | `main/flex/simple_flex_pca_pcg.f90:1641` `cg_core` |

These procedures have the most comment words in absolute terms:

| comment words | code lines | procedure |
|---:|---:|---|
| 1582 | 964 | `main/commanders/test/simple_commanders_test_highlevel.f90:1203` `exec_test_pcg_recon` |
| 1408 | 285 | `main/flex/simple_flex_pca_em_iter.f90:23` `probe_subspace_iteration` |
| 1238 | 419 | `main/flex/simple_flex_pca_model.f90:82` `run_flex_pca` |
| 1043 | 364 | `main/flex/simple_flex_pca_em_mstep.f90:28` `fit_iter_finish` |
| 923 | 347 | `main/flex/simple_flex_pca_em_embed.f90:15` `embed_latents_with_contrast` |
| 830 | 331 | `main/flex/simple_flex_pca_merge.f90:270` `two_gate_state_merge` |
| 821 | 282 | `main/commanders/simple/simple_commanders_refine3D.f90:49` `exec_refine3D_auto` |
| 768 | 399 | `main/stream/simple_stream_p03_initial_analysis.f90:83` `exec_stream_p03_initial_analysis` |
| 764 | 457 | `main/commanders/simple/simple_commanders_abinitio.f90:767` `exec_abinitio3D` |
| 679 | 285 | `main/flex/simple_flex_pca_weights.f90:27` `build_covariance_state_weights` |
| 673 | 341 | `main/nu_filt/simple_nu_filter_evidence.f90:16` `build_nu_evidence_state` |
| 619 | 175 | `main/commanders/simple/simple_commanders_volops.f90:292` `postprocess_volume_from_files` |

## E. Inventory

### Roll-up by subsystem

| subsystem | >=150 w | 100-149 w | 80-99 w | 50-79 w | words in blocks >=80 w | audited novels stale/wrong |
|---|---:|---:|---:|---:|---:|---:|
| main/flex | 23 | 27 | 36 | 102 | 11,437 | 18 of 23 |
| main/commanders | 6 | 12 | 9 | 46 | 3,437 | 4 of 6 |
| main/strategies | 3 | 10 | 14 | 25 | 3,222 | 2 of 3 |
| main/volume | 6 | 8 | 5 | 21 | 2,853 | 2 of 6 |
| main/nu_filt | 9 | 2 | 9 | 20 | 2,762 | 5 of 9 |
| utils/gui | 3 | 13 | 5 | 4 | 2,687 | 3 of 3 |
| main/stream | 5 | 3 | 1 | 11 | 1,406 | 5 of 5 |
| utils/persistent_worker | 6 | 1 | 0 | 5 | 1,267 | 6 of 6 |
| main/abinitio | 3 | 4 | 1 | 8 | 1,075 | 0 of 3 |
| main/image | 1 | 3 | 5 | 23 | 960 | 0 of 1 |
| main/pftc | 2 | 3 | 1 | 4 | 790 | 2 of 2 |
| fileio | 1 | 2 | 3 | 4 | 659 | 0 of 1 |
| main/project | 2 | 0 | 2 | 5 | 600 | 0 of 2 |
| utils/(root) | 1 | 1 | 3 | 0 | 532 | 1 of 1 |
| main/ori | 1 | 2 | 0 | 5 | 495 | 1 of 1 |
| main/nano | 1 | 1 | 2 | 4 | 467 | 1 of 1 |
| utils/math | 0 | 3 | 1 | 9 | 459 | - |
| utils/filter | 1 | 1 | 2 | 6 | 446 | 0 of 1 |
| main/interp | 0 | 2 | 2 | 1 | 398 | - |
| utils/comm | 0 | 2 | 2 | 3 | 380 | - |
| main/class | 1 | 1 | 1 | 4 | 376 | 0 of 1 |
| main/pca | 0 | 1 | 2 | 5 | 336 | - |
| main/params | 1 | 0 | 1 | 2 | 308 | 0 of 1 |
| utils/qsys | 0 | 2 | 1 | 0 | 296 | - |
| main/image_processing | 0 | 1 | 1 | 3 | 189 | - |
| main/sigma2 | 0 | 0 | 2 | 2 | 185 | - |
| main/sieve | 1 | 0 | 0 | 7 | 175 | 1 of 1 |
| production | 1 | 0 | 0 | 0 | 163 | 1 of 1 |
| main/opt | 0 | 1 | 0 | 1 | 109 | - |
| main/(root) | 0 | 1 | 0 | 3 | 107 | - |
| utils/clustering | 0 | 0 | 1 | 0 | 92 | - |
| main/ctf | 0 | 0 | 1 | 4 | 82 | - |
| main/ui | 0 | 0 | 0 | 4 | 0 | - |
| main/star | 0 | 0 | 0 | 4 | 0 | - |
| defs | 0 | 0 | 0 | 3 | 0 | - |
| **total** | **78** | **107** | **113** | **354** | **38,750** | **52 of 78** |

### E1. Every block of 150 words or more, with audit verdict

Verdicts:

- **WRONG**: a claim that the current code contradicts.
- **STALE**: an outdated path, test name, date or switch, or a claim that is incomplete in a way
  that matters.
- **OK**: accurate, but still long.

Some OK blocks repeat a policy note almost word for word, e.g. `simple_project_superset.f90` and
`abinitio3D_addon_policy.md` §4. Those can shrink to a pointer.


**fileio**

| words | lines | where | verdict | finding |
|---:|---:|---|---|---|
| 179 | 17 | `fileio/simple_flex_weights_file.f90:1-17` (file header) | OK | accurate; carries a date |

**main/abinitio**

| words | lines | where | verdict | finding |
|---:|---:|---|---|---|
| 194 | 16 | `main/abinitio/simple_abinitio3D_manifest_tester.f90:1-16` (file header) | OK | "optics/CTF" overstated: only dfx is varied |
| 183 | 14 | `main/abinitio/simple_abinitio3D_addon_report.f90:1-14` (file header) | OK | "base mask" is really the mskdiam sphere |
| 159 | 15 | `main/abinitio/simple_abinitio_controller.f90:549-563` (in `emit_refine3D_stage_cfg`) | OK | accurate; repeats abinitio3D_policy.md §4 |

**main/class**

| words | lines | where | verdict | finding |
|---:|---:|---|---|---|
| 188 | 15 | `main/class/simple_cavg_sums.f90:1-15` (file header) | OK | accurate |

**main/commanders**

| words | lines | where | verdict | finding |
|---:|---:|---|---|---|
| 327 | 44 | `main/commanders/test/simple_commanders_test_highlevel.f90:1159-1202` (`exec_test_pcg_recon` header) | STALE | lists 12 stages, code runs 14; cites moved policy path |
| 214 | 19 | `main/commanders/test/simple_commanders_test_highlevel.f90:3482-3500` (`run_abinitio3D_addon_gate` header) | OK | accurate |
| 202 | 14 | `main/commanders/simple/simple_commanders_flex_pca.f90:16-29` (module spec) | WRONG | "these MATCH the UI; the cap is 20": constant is 4, UI says 20/1e-3/gridding |
| 177 | 21 | `main/commanders/simple/simple_commanders_postprocess_nu.f90:1-21` (file header) | WRONG | describes retired v1; "envelope uses the compact state" false; standard postprocess was rewritten |
| 162 | 15 | `main/commanders/test/simple_commanders_test_highlevel.f90:2817-2831` (`exec_test_rec3D_backends` header) | STALE | describes retired /box convention and running in CWD |
| 158 | 17 | `main/commanders/simple/simple_commanders_imgops.f90:505-521` (in `exec_ppca_denoise`) | OK | accurate; "cascade blurrier" unrecorded |

**main/flex**

| words | lines | where | verdict | finding |
|---:|---:|---|---|---|
| 423 | 41 | `main/flex/simple_flex_pca_polar.f90:1-41` (file header) | WRONG | "the ONLY approximation is the direction snap": \|T\|^2 ring mean too; "bank costs no memory" false |
| 389 | 34 | `main/flex/simple_flex_pca_em_pairmerge.f90:1-34` (file header) | STALE | SIMPLE_COV_PAIRED_MERGE has no reader; mode-2 step never runs; caller flow changed |
| 280 | 23 | `main/flex/simple_flex_pca_em_iter.f90:132-154` (in `probe_subspace_iteration`) | STALE | cap is off by default (unstated); STRIDE vs cap headings contradict; duplicates cov_stage_subsample |
| 257 | 19 | `main/flex/simple_flex_pca_pcg.f90:1-19` (file header) | WRONG | "every solve starts from ZERO": solve hard-wires iwarm=1 (pcg:1907) |
| 254 | 20 | `main/flex/simple_flex_pca_em_fit.f90:34-53` (in `build_covariance_eigenbasis`) | WRONG | describes removed COLS/SNR/SOLVE rounds and probe_worker_pass; routine is a wrapper of init_basis_datafree |
| 247 | 19 | `main/flex/simple_flex_pca_plane_cache.f90:1-19` (file header) | WRONG | "selection hash mismatch rebuilds": hash written, never checked (:141-145) |
| 245 | 20 | `main/flex/simple_flex_pca_em_compose.f90:1-20` (file header) | STALE | "embedded ONCE with the union basis": compose_cut_reembed re-embeds by default |
| 219 | 17 | `main/flex/simple_flex_pca_deconv.f90:1-17` (file header) | OK | accurate; omits XD_CV_MAX subsample |
| 209 | 17 | `main/flex/simple_flex_pca_planes.f90:1-17` (file header) | WRONG | SIMPLE_COV_RESIDENT(_GB) read nowhere; budget always 0.25*MemAvailable |
| 208 | 15 | `main/flex/simple_flex_pca_em_estep.f90:933-947` (`box_reset` header) | WRONG | "band boxing" documents box_* helpers that nothing calls; v5 ships index lists |
| 203 | 14 | `main/flex/simple_flex_pca_em_solve.f90:191-204` (`probe_solve_ecm` header) | WRONG | first 6 lines document spd_inv_dp (line 479), not probe_solve_ecm |
| 198 | 16 | `main/flex/simple_flex_pca_em_pairmerge.f90:108-123` (`init_deflated_matchcos` header) | WRONG | first lines are probe_paired_merge's header, misplaced; "rank gate consumes sub_cos" outdated |
| 194 | 22 | `main/flex/simple_flex_pca_merge.f90:1-22` (file header) | STALE | "npreimages=0 enables it implicitly": preimage_auto=yes does (model:150-154) |
| 187 | 18 | `main/flex/simple_flex_reconstructor_latent_ops.f90:799-816` (`add_invtausq2rho_coupled` header) | STALE | cites simple_reconstructor.f90:1151; the line is now :1160 |
| 187 | 14 | `main/flex/simple_flex_pca_em_crossfsc.f90:118-131` (`xfsc_build_invtau2` header) | STALE | "arm 2 is the recommended first A/B": arm hard-wired to 1; "spec par.4.x" doc absent |
| 175 | 16 | `main/flex/simple_flex_pca_em.f90:129-144` (module spec) | OK | accurate; crash story and "proposal §4" cite no document |
| 172 | 18 | `main/flex/simple_flex_weights_state.f90:1-18` (file header) | OK | accurate |
| 169 | 15 | `main/flex/simple_flex_pca_em_fit.f90:126-140` (`init_basis_datafree` header) | STALE | SIMPLE_COV_EM=1 has no reader; "the column path" was removed |
| 167 | 15 | `main/flex/simple_flex_pca_em_polar.f90:21-35` (module spec) | WRONG | argues from the removed moment estimator; "ndir costs no memory" false (estep:139) |
| 165 | 13 | `main/flex/simple_flex_pca_em_iter.f90:494-506` (in `fit_stage_config`) | WRONG | fit_stage_config "reads env switches": reads none; "contrast pinned at 1" contradicts estep:221 |
| 164 | 14 | `main/flex/simple_flex_pca_em_crossfsc.f90:209-222` (`xfsc_paired_record` header) | STALE | Wiener call site moved to em_mstep:151; cites ~/ribo_local path and "impl-map step-3" |
| 156 | 13 | `main/flex/simple_flex_pca_em_polar.f90:142-154` (`project_fplane_mean_banded` header) | OK | accurate; ~80% speed figure is an unrecorded measurement |
| 151 | 11 | `main/flex/simple_flex_pca_tester.f90:1-11` (file header) | OK | accurate |

**main/image**

| words | lines | where | verdict | finding |
|---:|---:|---|---|---|
| 156 | 13 | `main/image/simple_image_tester.f90:690-702` (`test_padding_query` header) | OK | accurate |

**main/nano**

| words | lines | where | verdict | finding |
|---:|---:|---|---|---|
| 199 | 25 | `main/nano/simple_atoms.f90:12-36` (module spec) | WRONG | PDB table disagrees with the formats below it (chainID, segID, iCode) |

**main/nu_filt**

| words | lines | where | verdict | finding |
|---:|---:|---|---|---|
| 341 | 39 | `main/nu_filt/simple_nu_filter_sharpen.f90:1-39` (file header) | OK | accurate; B is from the evidence-pair average |
| 216 | 19 | `main/nu_filt/simple_nu_filter.f90:60-78` (module spec) | WRONG | "cut bounds the lead to 1.25-1.5x": only the upper bound holds; nu_refine control is gone |
| 215 | 27 | `main/nu_filt/simple_nu_filter_envmask.f90:1-27` (file header) | WRONG | margin formula: computed from raw unaries, not dmats_mask |
| 192 | 26 | `main/nu_filt/simple_nu_filter.f90:1-26` (file header) | STALE | shell walk removed 2026-09-16 (b9cc8e31b), not 09-18; cites build ed36eb4c, a NICE-only commit |
| 173 | 18 | `main/nu_filt/simple_nu_filter_bank.f90:65-82` (in `setup_nu_dmats`) | OK | accurate; "finest rung" should read "finest retained rung" |
| 167 | 18 | `main/nu_filt/simple_nu_filter_evidence.f90:104-121` (in `build_nu_evidence_state`) | OK | accurate (the policy doc is the stale side) |
| 164 | 17 | `main/nu_filt/simple_nu_filter.f90:261-277` (module spec) | OK | accurate |
| 155 | 14 | `main/nu_filt/simple_nu_filter_bank.f90:295-308` (in `optimize_nu_cutoff_finds`) | WRONG | "n(n+1)/2-1 instead of n" smoothing passes: the n passes are still paid |
| 151 | 12 | `main/nu_filt/simple_nu_filter.f90:179-190` (module spec) | WRONG | "flex_pca report quotes this floor": flex_pca gets the 5% default (stats:156) |

**main/ori**

| words | lines | where | verdict | finding |
|---:|---:|---|---|---|
| 213 | 18 | `main/ori/simple_oris_sampling.f90:78-95` (`population_blend_weights` header) | WRONG | "mnew = N whenever M>0 or n=N": fails for M=0, n=N, u<1 |

**main/params**

| words | lines | where | verdict | finding |
|---:|---:|---|---|---|
| 212 | 31 | `main/params/simple_parameters.f90:1-31` (file header) | OK | accurate recipe; omits UI layer (inline !< at 222/303/356 stale) |

**main/pftc**

| words | lines | where | verdict | finding |
|---:|---:|---|---|---|
| 186 | 21 | `main/pftc/simple_pftc_shsrch_grad.f90:416-436` (`grad_shsrch_minimize_joint` header) | WRONG | "without irot_in, an all-angle scan": irot_in is required, no scan path |
| 153 | 19 | `main/pftc/simple_polarft_corr.f90:1548-1566` (`gen_corr_grad_at_angle` header) | STALE | cites simple_test_continuous_inplane_cc_grad, merged into simple_pftc_inplane_tester |

**main/project**

| words | lines | where | verdict | finding |
|---:|---:|---|---|---|
| 251 | 23 | `main/project/simple_project_superset.f90:1-23` (file header) | OK | accurate; near-verbatim copy of addon policy §4 |
| 183 | 15 | `main/project/simple_project_superset_tester.f90:1-15` (file header) | OK | accurate test plan |

**main/sieve**

| words | lines | where | verdict | finding |
|---:|---:|---|---|---|
| 175 | 37 | `main/sieve/simple_ptcl_sieve.f90:1-37` (file header) | STALE | workflow order differs from cycle(); REJECTION_FAILED only written by tester |

**main/strategies**

| words | lines | where | verdict | finding |
|---:|---:|---|---|---|
| 520 | 53 | `main/strategies/search/simple_ptcl_cache.f90:1-53` (file header) | WRONG | documents 3D consumers; 3D rejects cache=yes (refine3D_strategy:633, abinitio:823) |
| 177 | 18 | `main/strategies/search/simple_ptcl_cache.f90:322-339` (`cache_key_stkfp` header) | OK | accurate |
| 168 | 14 | `main/strategies/parallelization/simple_flex_pca_strategy.f90:221-234` (in `master_initialize`) | WRONG | "SIMPLE_COV_MASTER_NTHR overrides": read nowhere; block's last line says no override |

**main/stream**

| words | lines | where | verdict | finding |
|---:|---:|---|---|---|
| 243 | 41 | `main/stream/simple_stream_p06_pool2D_new.f90:1-41` (file header) | STALE | omits stage-7 export, pause policy, sieverefs; stepwise cap is cumulative |
| 205 | 35 | `main/stream/simple_stream_p03_initial_analysis.f90:1-35` (file header) | WRONG | documents micimporter, shape ranking, process_selected_refs, restart_requested: none run |
| 196 | 40 | `main/stream/simple_stream_p05_sieve_cavgs_new.f90:1-40` (file header) | STALE | sieve created lazily, not at start-up; flag is single_pass not coarse_only |
| 178 | 36 | `main/stream/simple_stream_p07_abinitio3D_multistate.f90:1-36` (file header) | STALE | "refine3D not yet wired": abinitio3D_addon is wired |
| 164 | 27 | `main/stream/simple_stream_p04_refpick_extract_new.f90:1-27` (file header) | WRONG | input is preprocess, not "stage 3 CTF"; thresholds only with reject_mics=yes; no particle STAR |

**main/volume**

| words | lines | where | verdict | finding |
|---:|---:|---|---|---|
| 447 | 44 | `main/volume/simple_pcg_solvent_sidecar.f90:1-44` (file header) | WRONG | "FSC is solvent-flattened": FSC comes from the prior-free pair (contradicts own lines 21-22) |
| 293 | 25 | `main/volume/simple_frozen_accum.f90:1-25` (file header) | OK | accurate; duplicates abinitio3D_addon_policy.md §5, §7 |
| 177 | 20 | `main/volume/simple_nu_state_filter.f90:1-20` (file header) | STALE | "PCG carries no prior of its own": pcg_solvent=yes passes an apply pair |
| 176 | 18 | `main/volume/simple_reconstructor.f90:698-715` (`floor_rho_shellwise` header) | OK | accurate; measurements and RELION remarks unrecorded |
| 176 | 22 | `main/volume/simple_nu_state_filter.f90:52-73` (`nonuniform_filter_state` header) | OK | accurate; argument catalogue + "2026-09-21, Hans" |
| 160 | 14 | `main/volume/simple_frozen_accum_tester.f90:1-14` (file header) | OK | accurate test plan |

**production**

| words | lines | where | verdict | finding |
|---:|---:|---|---|---|
| 163 | 30 | `production/simple_persistent_worker.f90:1-30` (file header) | WRONG | worker_id is set by qsys_env, not server; script-path validation is bypassed |

**utils/(root)**

| words | lines | where | verdict | finding |
|---:|---:|---|---|---|
| 159 | 30 | `utils/simple_forked_process.f90:1-30` (file header) | STALE | restart count off by one; status() can re-fork (side effect) |

**utils/filter**

| words | lines | where | verdict | finding |
|---:|---:|---|---|---|
| 152 | 10 | `utils/filter/simple_lpstages_tester.f90:1-10` (file header) | OK | accurate; lpstages_ref.py is not in the repo |

**utils/gui**

| words | lines | where | verdict | finding |
|---:|---:|---|---|---|
| 263 | 38 | `utils/gui/metadata/simple_gui_metadata_vol3D.f90:1-38` (file header) | STALE | field list omits cfar, i, i_max; "set/get for all fields" false |
| 227 | 42 | `utils/gui/simple_gui_assembler.f90:1-42` (file header) | STALE | does not send; heartbeats never hash-suppressed |
| 197 | 36 | `utils/gui/metadata/stream/simple_gui_metadata_stream_update.f90:1-36` (file header) | WRONG | "only changed fields are sent": whole object is sent; omits pickrefs_cycle |

**utils/persistent_worker**

| words | lines | where | verdict | finding |
|---:|---:|---|---|---|
| 244 | 39 | `utils/persistent_worker/simple_persistent_worker_server.f90:1-39` (file header) | STALE | "workers connect periodically": one persistent connection, server polls |
| 234 | 50 | `utils/persistent_worker/message/simple_persistent_worker_message_task.f90:1-50` (file header) | STALE | dispatch goes out as WORKER_NEW_TASK_MSG; timing fields never set |
| 182 | 40 | `utils/persistent_worker/message/simple_persistent_worker_message_heartbeat.f90:1-40` (file header) | WRONG | HEARTBEAT_TIMEOUT_MS unused; "four fields" (six); reply can be STATUS |
| 175 | 35 | `utils/persistent_worker/message/simple_persistent_worker_message_status.f90:1-35` (file header) | STALE | omits its role as queue_task acknowledgement |
| 166 | 34 | `utils/persistent_worker/message/simple_persistent_worker_message_base.f90:1-34` (file header) | WRONG | USAGE redeclares msg_type in an extension (illegal); contract text contradicts itself |
| 157 | 34 | `utils/persistent_worker/message/simple_persistent_worker_message_terminate.f90:1-34` (file header) | WRONG | "worker finishes its current tasks": worker cancels them (worker:281-299) |

### E2. Blocks of 100-149 words (not audited)

"dated" marks a block that carries a date or a name.

| words | lines | where | dated |
|---:|---:|---|---|
| 116 | 12 | `fileio/simple_syslib.f90:421-432` (`io_read_nstreams` header) |  |
| 113 | 7 | `fileio/simple_imghead_tester.f90:1-7` (file header) | yes |
| 107 | 12 | `main/simple_final_rec.f90:18-29` (`calc_final_rec` header) |  |
| 123 | 10 | `main/abinitio/simple_abinitio3D_manifest.f90:1-10` (file header) |  |
| 122 | 11 | `main/abinitio/simple_abinitio_controller.f90:567-577` (in `emit_refine3D_stage_cfg`) | yes |
| 112 | 10 | `main/abinitio/simple_abinitio_controller.f90:536-545` (in `emit_refine3D_stage_cfg`) | yes |
| 100 | 9 | `main/abinitio/simple_abinitio_controller.f90:683-691` (in `emit_refine3D_stage_cfg`) | yes |
| 103 | 7 | `main/class/simple_classaverager_restore.f90:536-542` (`commit_carryover` header) |  |
| 144 | 12 | `main/commanders/test/simple_commanders_test_highlevel.f90:3542-3553` (in `run_abinitio3D_addon_gate`) | yes |
| 140 | 15 | `main/commanders/simple/simple_commanders_imgops.f90:542-556` (in `exec_ppca_denoise`) |  |
| 138 | 12 | `main/commanders/test/simple_test_truth_metrics.f90:1-12` (file header) |  |
| 118 | 12 | `main/commanders/simple/simple_commanders_refine3D.f90:64-75` (in `exec_refine3D_auto`) | yes |
| 117 | 14 | `main/commanders/simple/simple_commanders_refine3D.f90:2082-2095` (`exec_bootstrap_rec3D` header) |  |
| 112 | 13 | `main/commanders/simple/simple_commanders_refine3D.f90:211-223` (in `exec_refine3D_auto`) |  |
| 110 | 23 | `main/commanders/test/simple_commanders_test_class.f90:108-130` (module spec) |  |
| 110 | 10 | `main/commanders/simple/simple_commanders_abinitio.f90:1360-1369` (`addon_prologue` header) |  |
| 106 | 8 | `main/commanders/simple/simple_commanders_ori.f90:624-631` (`exec_measure_projspace_angres` header) |  |
| 104 | 10 | `main/commanders/test/simple_commanders_test_highlevel.f90:1857-1866` (in `exec_test_pcg_recon`) |  |
| 101 | 9 | `main/commanders/simple/simple_commanders_rec_distr.f90:194-202` (in `blend_trailing_accumulators`) |  |
| 101 | 11 | `main/commanders/simple/simple_commanders_refine3D.f90:285-295` (`run_registration_pass` header) |  |
| 145 | 11 | `main/flex/simple_flex_pca_em_mean.f90:341-351` (`cov_half_parity` header) | yes |
| 142 | 10 | `main/flex/simple_flex_pca_em_estep.f90:1008-1017` (`pk_mask_r2` header) |  |
| 142 | 13 | `main/flex/simple_flex_pca_em_iter.f90:55-67` (in `probe_subspace_iteration`) |  |
| 137 | 11 | `main/flex/simple_flex_pca_polar.f90:385-395` (`polar_sample_particle_fused` header) |  |
| 129 | 12 | `main/flex/simple_flex_pca_polar.f90:314-325` (`polar_sample_particle` header) |  |
| 129 | 10 | `main/flex/simple_flex_pca_em.f90:245-254` (module spec) |  |
| 127 | 10 | `main/flex/simple_flex_pca_pcg.f90:2202-2211` (`test_flex_pcg_operator` header) |  |
| 126 | 10 | `main/flex/simple_flex_pcg_tester.f90:1-10` (file header) |  |
| 124 | 9 | `main/flex/simple_flex_pca_pcg.f90:1898-1906` (in `solve`) |  |
| 123 | 10 | `main/flex/simple_flex_pca_merge.f90:427-436` (in `two_gate_state_merge`) | yes |
| 120 | 9 | `main/flex/simple_flex_pca_em_iter.f90:456-464` (in `fit_stage_config`) |  |
| 119 | 10 | `main/flex/simple_flex_pca_pcg.f90:1848-1857` (in `solve`) |  |
| 119 | 10 | `main/flex/simple_flex_pca_em_iter.f90:69-78` (in `probe_subspace_iteration`) |  |
| 119 | 10 | `main/flex/simple_flex_pca_em.f90:272-281` (module spec) | yes |
| 117 | 9 | `main/flex/simple_flex_pca_pcg.f90:1689-1697` (in `cg_core`) |  |
| 117 | 10 | `main/flex/simple_flex_pca_em_fit.f90:229-238` (`em_calibrate_noise_prior` header) |  |
| 116 | 8 | `main/flex/simple_flex_pca_pcg.f90:1816-1823` (`solve` header) |  |
| 114 | 8 | `main/flex/simple_flex_pca_em_iter.f90:484-491` (in `fit_stage_config`) |  |
| 113 | 9 | `main/flex/simple_flex_pca_rounds.f90:1-9` (file header) |  |
| 111 | 8 | `main/flex/simple_flex_pca_pcg.f90:1308-1315` (in `apply_operator`) |  |
| 111 | 10 | `main/flex/simple_flex_pca_crossfsc.f90:1-10` (file header) | yes |
| 111 | 9 | `main/flex/simple_flex_pca_em.f90:230-238` (module spec) |  |
| 109 | 9 | `main/flex/simple_flex_pca_em_polar.f90:229-237` (`polar_hybrid_exact_accum` header) |  |
| 109 | 11 | `main/flex/simple_flex_pca_em_embed.f90:436-446` (`write_embed_stats_part` header) |  |
| 108 | 7 | `main/flex/simple_flex_pca_pcg.f90:1097-1103` (in `finalize`) |  |
| 105 | 8 | `main/flex/simple_flex_pca_polar.f90:94-101` (`polar_grid_build` header) |  |
| 100 | 8 | `main/flex/simple_flex_pca_em_basis.f90:206-213` (`align_basis_to_reference` header) |  |
| 134 | 16 | `main/image/simple_image_calc.f90:1925-1940` (in `nu_objective`) |  |
| 103 | 10 | `main/image/simple_image_calc.f90:1777-1786` (`nu_objective_noise_profile` header) |  |
| 103 | 8 | `main/image/simple_image_tester.f90:1-8` (file header) | yes |
| 108 | 15 | `main/image_processing/simple_segmentation.f90:333-347` (`canny` header) |  |
| 124 | 16 | `main/interp/simple_kbinterpol.f90:1-16` (file header) |  |
| 106 | 16 | `main/interp/simple_gridding.f90:63-78` (`kb_stencil_envelope_1d` header) |  |
| 106 | 7 | `main/nano/simple_atoms_tester.f90:1-7` (file header) | yes |
| 112 | 10 | `main/nu_filt/simple_nu_filter.f90:102-111` (module spec) | yes |
| 110 | 12 | `main/nu_filt/simple_nu_filter_stats.f90:120-131` (`get_nu_filtmap_finest_selected_lp` header) | yes |
| 109 | 12 | `main/opt/simple_opt_subs.f90:17-28` (`amoeba` header) |  |
| 142 | 15 | `main/ori/simple_oris_dists.f90:38-52` (`find_angres` header) |  |
| 140 | 14 | `main/ori/simple_oris_reshape.f90:148-161` (`reseed_classes` header) |  |
| 145 | 8 | `main/pca/simple_ppca.f90:185-192` (`calc_bic_ppca` header) | yes |
| 140 | 9 | `main/pftc/simple_polarft_corr_tester.f90:1-9` (file header) |  |
| 125 | 10 | `main/pftc/simple_pftc_inplane_tester.f90:1-10` (file header) |  |
| 105 | 11 | `main/pftc/simple_polarft_corr.f90:1417-1427` (`gen_raw_euclid_grad_at_angle` header) |  |
| 143 | 13 | `main/strategies/search/simple_ptcl_cache.f90:141-153` (`cache_fname` header) |  |
| 134 | 9 | `main/strategies/search/simple_cavg_registration_tester.f90:1-9` (file header) |  |
| 109 | 9 | `main/strategies/search/simple_pose_cont_1jyx_tester.f90:1-9` (file header) |  |
| 106 | 9 | `main/strategies/search/simple_ptcl_cache.f90:123-131` (`cache_run_token` header) |  |
| 106 | 12 | `main/strategies/search/simple_matcher_refvol_utils.f90:423-434` (`mask_matching_reference` header) | yes |
| 105 | 10 | `main/strategies/search/simple_ptcl_cache.f90:587-596` (`ptcl_cache_ensure` header) |  |
| 105 | 12 | `main/strategies/parallelization/simple_calc_pspec_strategy.f90:259-270` (`compute_pspec_partitions` header) |  |
| 104 | 11 | `main/strategies/search/simple_matcher_3Drec.f90:431-441` (`prep_imgs4rec` header) |  |
| 103 | 10 | `main/strategies/search/simple_matcher_ptcl_io.f90:452-461` (`prep_rec_observation` header) | yes |
| 100 | 11 | `main/strategies/search/simple_ptcl_cache.f90:840-850` (in `ptcl_cache_read_batch`) |  |
| 117 | 9 | `main/stream/simple_stream_tester.f90:1-9` (file header) | yes |
| 110 | 18 | `main/stream/simple_stream_abinitio2D_chunks.f90:1-18` (file header) |  |
| 108 | 18 | `main/stream/simple_stream_p00_master.f90:1-18` (file header) |  |
| 147 | 12 | `main/volume/simple_pcg_halfset_tester.f90:1-12` (file header) | yes |
| 137 | 13 | `main/volume/simple_reconstructor_pcg.f90:2920-2932` (`shrink_by_ml_prior` header) | yes |
| 132 | 12 | `main/volume/simple_cartesian_pose_refiner_tester.f90:1-12` (file header) |  |
| 130 | 12 | `main/volume/simple_reconstructor_pcg.f90:563-574` (`fold_solvent_ridge_into_precond` header) | yes |
| 123 | 12 | `main/volume/simple_nu_state_filter.f90:242-253` (in `record_nu_alignment_lowpass_limit`) | yes |
| 114 | 13 | `main/volume/simple_nu_state_filter.f90:139-151` (in `nonuniform_filter_state`) | yes |
| 107 | 14 | `main/volume/simple_halfmap_diagnostics.f90:160-173` (`support_provenance_fname` header) | yes |
| 102 | 9 | `main/volume/simple_pcg_solvent_sidecar.f90:68-76` (module spec) | yes |
| 114 | 25 | `utils/simple_forked_process_tester.f90:1-25` (file header) |  |
| 108 | 8 | `utils/comm/simple_http_post_tester.f90:1-8` (file header) |  |
| 104 | 11 | `utils/comm/simple_ipc_tcp_socket_helpers.f90:104-114` (`poll_fds` header) |  |
| 121 | 8 | `utils/filter/simple_bspline_smoother_tester.f90:1-8` (file header) |  |
| 144 | 27 | `utils/gui/metadata/stream/simple_gui_metadata_stream_abinitio3D_multistate.f90:1-27` (file header) |  |
| 139 | 25 | `utils/gui/metadata/stream/simple_gui_metadata_stream_pool2D.f90:1-25` (file header) |  |
| 132 | 22 | `utils/gui/metadata/simple_gui_metadata_ptcl.f90:1-22` (file header) |  |
| 126 | 25 | `utils/gui/metadata/simple_gui_metadata_project.f90:1-25` (file header) |  |
| 126 | 24 | `utils/gui/metadata/stream/simple_gui_metadata_stream_opening2D.f90:1-24` (file header) |  |
| 124 | 23 | `utils/gui/metadata/stream/simple_gui_metadata_stream_particle_sieving.f90:1-23` (file header) |  |
| 123 | 24 | `utils/gui/simple_gui_assembler_tester.f90:1-24` (file header) |  |
| 116 | 21 | `utils/gui/metadata/simple_gui_metadata_cavg2D.f90:1-21` (file header) |  |
| 104 | 26 | `utils/gui/metadata/simple_gui_metadata_tester.f90:1-26` (file header) |  |
| 103 | 21 | `utils/gui/metadata/stream/simple_gui_metadata_stream_preprocess.f90:1-21` (file header) |  |
| 102 | 21 | `utils/gui/metadata/stream/simple_gui_metadata_stream_optics_assignment.f90:1-21` (file header) |  |
| 102 | 21 | `utils/gui/metadata/stream/simple_gui_metadata_stream_picking.f90:1-21` (file header) |  |
| 101 | 19 | `utils/gui/metadata/simple_gui_metadata_base.f90:1-19` (file header) |  |
| 141 | 14 | `utils/math/simple_rnd.f90:514-527` (`r8po_fa` header) |  |
| 119 | 12 | `utils/math/simple_math.f90:1195-1206` (`SavitzkyGolay_filter` header) |  |
| 110 | 7 | `utils/math/simple_linalg_tester.f90:277-283` (`test_fit_straight_line_recovery` header) |  |
| 109 | 22 | `utils/persistent_worker/message/simple_persistent_worker_message_types.f90:1-22` (file header) |  |
| 107 | 8 | `utils/qsys/simple_qsys_env_tester.f90:1-8` (file header) |  |
| 100 | 22 | `utils/qsys/simple_qsys_persistent_worker.f90:1-22` (file header) |  |

Blocks of 80-99 words (113) and 50-79 words (354) are counted in the roll-up and not listed.
Re-running the scan regenerates the full list.

## F. Left alone

- **`src/main/opt/simple_opt_lbfgsb.f90`**: 15 blocks, 8,700 words. This is the netlib
  L-BFGS-B 3.0 documentation, ported with the code. Keeping it verbatim makes it easy to diff
  against upstream.
- **`src/extlibs`**: third-party code.
- **Inline `!<` descriptors in `simple_parameters.f90`**: 31 inline comments in the library are
  over 100 characters, mostly here. They are documentation only (`simple_args_generator.pl`
  strips them), and the UI metadata is the help text users see. Three are stale:
  - :222 says "requires rec_backend=pcg", but both backends validate;
  - :303 gives the default as `ptcl`, while the UI says `rand`;
  - :356 lists `prob_snhc` and omits `snhc_smpl`/`snhc_smpl_many`.
- **Commented-out code**: 165 blocks, 416 lines. These are not novels and belong in a separate
  cleanup.
- **UI template stubs** ("INPUT PARAMETER SPECIFICATIONS ... `<empty>`"): many lines but few words,
  so they are not counted.

## G. Destination docs to fix before content moves into them

The same audit found the target notes stale on the same points:

- `doc/policies/2D/particle_cache_policy.md` (:83, :133, :162) still describes cache-enabled
  abinitio3D.
- `doc/policies/3D/reconstruct3D_pcg_policy.md:508` describes the gate as having "nine" stages;
  the test now runs 14.
- `doc/policies/3D/refine3D_auto_policy.md:98-100` says the startup `reconstruct3D` runs only
  without a compatible volume. It always runs (`simple_commanders_refine3D.f90:235-253`).
- `doc/code_overview/test_inventory.md:47` says pcg_recon uses box 24; the code uses
  `BOX = 32` (:1211).
- `doc/algorithms/heterogeneity_analysis/flex_pca.md:63` gives the contrast bracket as
  [0.2, 3]; the code clamps to [0.1, 5]. Lines 78-81 describe an even/odd stop rule that is
  disabled (`simple_flex_pca_em_mstep.f90:470`).
- `.github/skills/simple-main-nu-filt/SKILL.md` lists the deleted `simple_nu_filter_extend.f90`
  and `nu_refine`.
- Not re-checked here: the audit also reported `doc/policies/NU/nonuniform_filtering_policy.md`
  (lines 13, 343, 368, 532) and `doc/policies/persistent_worker_policy.md` §5-7 (message-type
  table, mutex claims, queue removal) as stale.

## H. Code issues found in passing

These are not comment problems. They surfaced during the audit.

1. `simple_flex_pca_em_mstep.f90:428-429` relies on `fit%z` being zero-initialised.
   `simple_flex_pca_em_iter.f90:159` allocates it without initialising it, and only the
   in-process path zeroes it (:252). A distributed master in the polish stage therefore tests
   uninitialised memory. The effect is limited to one diagnostic log line.
2. `simple_stream_p07_abinitio3D_multistate.f90:179` vs :668-674: for `addon_it >= 2` the frozen
   input and the stage's output project both resolve to `abinitio3D_addon/it_<n>/...`. It
   probably should be `it_<n-1>`. Not traced end to end, because `CWD_GLOB` changes inside the
   loop.
3. `production/simple_persistent_worker.f90:235`: worker script-path validation is disabled
   (`safe_path = .true. ! for testing`).
4. The long comments are keeping dead code alive:
   - the `box_*` helpers in `simple_flex_pca_em_estep.f90` (no callers);
   - the mode-2 branch of the paired merge and the Wiener-off branch in `fit_iter_finish`;
   - the NU shell-walk helpers (`count_nu_walked_label_voxels`,
     `get_nu_filtmap_highres_shell_depth`, `nu_effective_base_label_for_candidate`,
     `retain_nu_filter_setup`), none with external callers.

## Suggested order of work

1. Correct the 29 WRONG blocks (table E1). The drafts note has an accurate short replacement for
   each.
2. Move the lab notebook into decision logs. `doc/implementation_notes/completed/pcg_decision_log.md`
   already exists. flex_pca and NU need one each, or a section in the existing notes.
3. Collapse the 51 templated headers to `!@descr:` plus at most 4 lines, dropping the
   DEPENDENCIES, routine and field lists.
4. Fix the C2 citations and delete the C1 switch mentions, together with any dead branches they
   describe.
5. Work through the 107 unaudited blocks of 100-149 words (E2), one subsystem at a time.

To keep the problem from growing back, a `scripts/check_comments.py` modelled on `check_descr.py`
could flag:

- comment blocks longer than N lines;
- dates or names in comments;
- `SIMPLE_*` switches that nothing reads;
- `doc/` paths that do not exist.
