# Policy-set audit and gap register, 2026-10-21

Planned work. Source: a read-only audit of every document under `doc/policies/`
against the code at `fd5b5bc8e` (working tree, 2026-10-21), by six auditors, one
per cluster (refine3D, solve3D/PCG/NU, 2D and particle I/O, stream and sieving,
infrastructure, preprocessing). Each concrete claim (names, ownership, defaults,
file names, tests, links) was checked with `rg` and by reading the code; nothing
was built or run. The raw reports, with `file:line` evidence for every finding
below, are in [policy_set_audit_2026-10-21/](policy_set_audit_2026-10-21/). They
are the record before the corrections; the register below says what was done.

The audit followed the restructuring of 2026-10-20
([structure2D_and_project_owner.md](../completed/structure2D_and_project_owner.md)).
Result: the names and files the policies cite almost all resolve, and the
refine3D and solve2D ownership text of that restructuring is accurate. The
drift is in the claims around it (defaults, stage tables, who calls what, test
names), in absolute statements the code no longer keeps, and in what no policy
says. Section 1 lists what was corrected, sections 2 and 3 what needs the
maintainer, section 4 the gaps, section 5 the structural debt.

## 1. Corrections made (documents only)

Confirmed STALE, DRIFT, OLDNAME and CONTRADICTS findings were corrected in place,
fixing the side the code supports: 51 existing documents edited (+1137, -759
lines; `git diff --stat -- doc` shows each). No code, skill or instruction file was
changed. Gap and structure findings, and every decision in section 2, were left
alone except where an existing sentence was false.
Larger examples: `refine3D_policy.md` said matching references are never
multiplied by an envelope (they are) and that `automsk=yes` needs an NU filter
mode (only `automsk=nu` does); `solve3D_policy.md` section 11 described the
add-on entering at stage 3 with a full bootstrap (it enters at the base run's last
stage after a registration pass, with an incremental final reconstruction);
`segdiampick_mics_policy.md` described a version of the code that was reverted the
day after it was committed and is rewritten; `class_average_bootstrap_policy.md`
had not followed the command since its rework of 2026-06-12; the persistent-worker
policy said batch scheduling reads only `JOB_FINISHED*`; false "implemented, not
yet compiled or run" status lines were removed. The three "tests are performed by
the user" sentences now name registered tests.

New files: `doc/policies/README.md` (index of every policy) and
`layering_and_project_ownership_policy.md` (the developer guidance; its section 8
is the ledger of today's deviations).

## 2. Decisions needed

Each is a place where the policy and the code disagree and either could be right.
The documents state the code's behaviour for now.

| # | Decision | Where |
|---|---|---|
| D1 | The PCG master writes the particle field unconditionally; add the `write` flag the gridding path has, or keep the exception | `simple_rec3D_pcg_strategy.f90:931`, `refine3D_policy.md` section 4 |
| D2 | The distributed `refine2D` master writes `ptcl2D` in `merge_algndocs` (no `write=.false.`), after assembly and in `distr_finalize_iteration`: three times per iteration, four in probabilistic mode | `simple_refine2D_strategy.f90:404,771,449`, `simple_prob_align_strategy.f90:314` |
| D3 | Done 2026-10-21 (C4): the distributed `make_cavgs` master now makes the class assignment once and writes it before the workers start; no worker writes | `simple_make_cavgs_strategy.f90` |
| D4 | The default phase-shift window (0 to 180) triggers the `>= 180` warning that the policy says it should not: change the threshold or the text | `simple_parameters_phases.f90:761-765`, `phase_shift_ctf_policy.md` |
| D5 | In `filt_mode=nonuniform_lpset` an explicit `lp` is replaced by the NU handoff, not kept as a hard override (read from one routine; verify), and the open-question paragraph of the NU policy needs closing | `refresh_matching_lp_from_project`, `nonuniform_filtering_policy.md` |
| D6 | The default multi-state `lpstop=6.0` is injected as a ceiling and acts in NU stage 5; the policy says only an explicit `lpstop` is a ceiling | `simple_commanders_solve3D.f90:1003,1163` |
| D7 | `res`/`res05` are not cleared on the `cavg_ini_ext` and `state=` routes, so an old `res05` can advance their entry stage | `simple_commanders_solve3D.f90:1294-1324,2192` |
| D8 | `ptcl_buffer_budget` branches on `qsys=='coarray'`, against "must not leak into strategy or domain branching": allow the read or move it behind `qsys_env` | `simple_matcher_ptcl_io.f90:343` |
| D9 | Extend the umbrella rule to "do not import an umbrella in new code" (58 of 172 files added since 2026-10-01 import `simple_core_module_api`) | `compile_time_policy.md` |
| D10 | Keep `segdiampick_mics` as the `mini_stream` validation program only, or retire it; the stream uses `segdiam_bin_picker`, which has no policy | `simple_mini_stream_utils.f90` |
| D11 | `refine2D` (public since 2026-10-20) does not expose `center_type`, `extr_lim`, `oritype` and leaves `cenlp`, `lpstart`, `lpstop` unresolved in the default audit | `scripts/default_audit.py` output |
| D12 | `ptcl_read_session` lives in `strategies/search` but `class/`, `volume/` and `flex/` use it: move it beside `stack_read_session` or declare the exception | `simple_matcher_ptcl_io.f90` |
| D13 | The pool's `finalise` writes the job protocol's marker files into the folder of the job it has cancelled: declare it the one allowed cross-write, or rely on `cancel` | `simple_stream_pool2D.f90:1156-1163` |
| D14 | Detector-defect curation exists twice with different semantics for EER gain defects (integrated micrograph versus particles extracted from movies): name the contract or unify | `simple_motion_correct.f90:569`, `simple_particle_extractor.f90:647` |
| D15 | `nparts` larger than the item count is an error in `motion_correct` and `preprocess`, a clamp in the other preprocessing strategies | preprocessing strategies |
| D16 | `bootstrap_cavgs` with `mkdir=no` overwrites the project file with the synthetic project | `simple_commanders_mkcavgs.f90:209` |
| D17 | Three `_completed` notes link to `_unfinished` notes that are not in this checkout (qsys completion, prob-align peak memory, `pose_cont` aftermath), as does `for_developers/index.md`: recreate them or repoint | `1bfc75404` split |

## 3. Code defects and hygiene (not changed)

Found while auditing; none is a documentation error. "Read" means the code was
re-read on 2026-10-21 after the audit and the claim holds; "auditor" means the
auditor's evidence was not re-checked. Nothing was run. Severity is the potential
effect, not a measured one: H wrong results or a lost end state, M wrong
intermediate state or a race, L extra I/O, dead code or wrong text.

**Confirmed by reading**

| # | Defect | Where | Sev |
|---|---|---|---|
| C1 | **Fixed 2026-10-21.** `solve3D` `symmetrize` with `l_symran` (start group has more symmetry than the target) randomizes and writes `spproj%os_ptcl3D` that was last read before the symmetry-search stage's refinement; the stage's poses are replaced by the previous stage's randomized poses, and the symmetric reconstruction reads them | `simple_solve3D_utils.f90:621-622`; called from `simple_commanders_solve3D.f90:1397` (`solve3D`) and `:266` (`solve3D_cavgs`), both before the project is re-read | M |
| C2 | **Fixed 2026-10-21.** A pool paused at iteration 25 or later is not woken by the sieve's empty final set (`runs_to_final` is false from iteration 25; only an import with particles or a GUI mask change unpauses), so it publishes no `pool_final=yes` and p07 never runs its final refinement | `simple_stream_stage_pool2D.f90:533-538,736-738,1100-1113` | H for the stream |
| C3 | **Fixed 2026-10-21.** The `nu_stage` key was in neither manifest key list, so an add-on of a base run with `nu_stage` > 5 and `nstages` 5 or 6 decides NU for its last stage differently from its base | `simple_solve3D_manifest.f90:33-60`, `simple_solve3D_controller.f90:427-433` | M, rare |
| C4 | **Fixed 2026-10-21.** The distributed `make_cavgs` worker part 1 rewrote the particle segment in place while the other parts read the file | `simple_make_cavgs_strategy.f90:246-255` | M (race) |
| C5 | `segdiampick_mics` reads only the first line of `mic_diam` and appends its five columns as diameters; no guard when nothing is picked | `simple_mini_stream_utils.f90:122-160,166` | L (validation program) |
| C6 | The PCG master writes the particle field unconditionally, with no `write` flag | `simple_rec3D_pcg_strategy.f90:931` | L |
| C7 | The distributed `refine2D` master writes `ptcl2D` in `merge_algndocs` (default `write`), after assembly and at the iteration's end: three times, four in probabilistic mode | `simple_refine2D_strategy.f90:404,449,771`, `simple_prob_align_strategy.f90:314` | L |
| C8 | `segdiam_bin_picker%pick` writes the caller's project twice (and on its early return); `extract`'s master writes the full project before and after its workers; `make_pickrefs_impl` calls `params%new` again inside the pick master and runs a commander | `simple_segdiam_bin_picker.f90:225,322-323`, `simple_extract_strategy.f90:588,719`, `simple_pick_strategy.f90:213,536` | L |
| C9 | `bootstrap_cavgs` writes the synthetic project to `params%projfile`; with `mkdir=no` it overwrites the user's project | `simple_commanders_mkcavgs.f90:210` | M, only with `mkdir=no` |
| C10 | `cavg_state.bin` is deleted only by `qsys_cleanup(keep2D=.false.)` and no caller passes it, so a finished run leaves it | `simple_qsys_funs.f90:43-48` | L |
| C11 | Wrong text: UI help for `pcg_solvent` says stage 7 (controller: 8); comment says "before NU" for the independent multi-state default (its fifth stage is the first NU stage); `exec_prob_tab` labels itself `simple_commanders_refine3D`, the `make_cavgs` worker `simple_commanders_refine2D :: exec_make_cavgs`; p06 stops with `SIMPLE_STREAM_SOLVE2D NORMAL STOP` | `simple_ui_solve3D.f90:77`, `simple_commanders_solve3D.f90:1001`, `simple_commanders_prob.f90:87,140`, `simple_commanders_stream_p06_pool2D.f90:60` | L |
| C12 | Dead code: `hac_1d_fast` (no caller, exported in the core API), `calc_3Drec_gpu` (a stub that stops), the members `cline_prob_align_distr` (set, never read) and `cline_volassemble` (never used) of the distributed `refine3D` strategy | `simple_math.f90:539`, `simple_reconstructor_openmpoffload.f90`, `simple_refine3D_strategy.f90:100-102,1110-1121` | L |
| C13 | The master forwards `optics_dir` to p03 and nothing in p03's stage or commander reads it | `simple_commanders_stream_p00_master.f90:494` | L |
| C14 | `ptcl_buffer_budget` branches on `qsys=='coarray'` against the coarray policy's "must not leak" rule | `simple_matcher_ptcl_io.f90:343` | L |
| C15 | The executable is recognised by a substring of the full `argv[0]` path (`private` in four places, `simple_private`, `simple_stream`). A binary run from a path that contains `private` (macOS `/private/tmp`, `/private/var`) is taken for `simple_private_exec`: `sp_required` becomes false, `mkdir=yes` makes the run directory but does not copy the project, the project and the sigma2 state are written at the original path, and the log redirect is skipped. Found 2026-10-21 when `solve3D_addon` run from `/private/tmp/...` failed to find `1_solve3D/addon_gate.simple`. Compare the basename instead | `simple_parameters_phases.f90:87,126,130,159`, `simple_exec_helpers.f90:221,242` | L, M for an install under such a path |

**The four fixes of 2026-10-21** (committed with the policy-set audit). Built (24 objects, no
warnings); fast gate 15 of 15. Under CTest, `lib_stream`, `lib_reconstruction`,
`simulated_workflow_1jxy` and `solve3D_addon` passed, but CTest puts `build/bin` on `PATH`
and sets `SIMPLE_PATH=build`, and that install was older than the fixes: the in-process
driver ran the new code, the worker processes (`simple_private_exec`) the old. The re-run
under `test=nightly` on a snapshot of the current binaries passed `solve3D_addon` (1171 s)
and `multistate_assignment` (331 s). It was then invalidated: a clean rebuild of `build/`
removed the snapshot, which sat inside it, so `cont_refine3D_1jxy` (48 s), `flex_pca_blobs`
and the rest did not run, and the failure of `multistate_workflows` (957 s, exit 2, no log
kept) cannot be told from the same loss. Not yet verified on current binaries: the other 16
gate entries, including `simulated_workflow_1jxy` (the end-to-end check of C4 with the new
worker), and the cause of the `multistate_workflows` failure. `make install` was not run.

- C1: `symmetry_randomize_poses` re-reads `ptcl3D` from the file before it
  randomizes and writes; both call sites go through `symmetrize`. Test: `solve3D_cavgs
  stage rules` (`unit_reconstruction`). With the read removed the test fails: none of
  the 12 file poses survives and 7 end as the stale poses.
- C2: `owes_final_run` (a pure rule) and a `l_final_published` flag. Past `FINAL_ITER`
  with the final set in the pool, the pool unpauses until a publication has carried
  `pool_final=yes`; a set after the final one takes the owed publication back. Test:
  `unit_stream` "pool 2D" (the rule and the flag's life cycle). The stage step itself
  is not driven by a test: it needs a pool at a chosen iteration.
- C3: `nu_stage` joins the recorded and the replayed manifest keys (33 and 20). Test:
  `solve3D manifest` (`unit_project`); with the replay key removed it fails with
  `expected=7 actual=0`. The manifest version is unchanged: a manifest without the key
  replays without it.
- C4: the master of `make_cavgs` runs `generate_cluster_centers` and writes the segment
  after the sigma2 rebuild and before the schedule; the worker no longer generates the
  assignment and never writes; `write_oritype_segment_for_compute` is shared-memory
  only. Class expansion draws random splits, so workers that each made it would have
  disagreed with one another and with the file. Covered end to end only by
  `simulated_workflow_1jxy` (two workers) on current binaries, pending in the nightly
  re-run; no unit test drives the strategy.

**Auditor evidence, not re-read**

- The add-on final reconstruction rewrites the working project in full three times
  (`switch_addon_membership`, `calc_final_rec_addon`).
- `distr_initialize` clears `updatecnt`/`sampled` in memory only, so a stand-alone
  `refine2D nparts=N refine=snhc_smpl startit=1` starts from the file's flags
  (reading supports it; `solve2D` masks it by writing the cleaned field first).
- Preprocessing: the masters of `motion_correct`, `ctf_estimate`, `pick`,
  `gen_pspecs_and_thumbs` and `preprocess` read the project twice and their workers
  read it whole; `extract`, `reextract` and `refine_motion_model` assemble part
  projects in three near-identical copies, with count checks only, and `extract`
  copies particle fields one by one so a new field is dropped silently; `nparts`
  larger than the item count is an error in two strategies and a clamp in the rest.
- Detector-defect curation exists twice with different EER gain-defect semantics.
- `default_audit.py` (auditor ran it): 8 errors on the stream programs, 20 warnings
  (six on `refine2D`), scans `commanders`, `exec`, `stream`, `defs`, `params` only
  (108 default-setting sites in `strategies` are invisible), and rewrites
  `simple_ui_default_values.f90` on every run.
- The commit hash is a `parameter` of `simple_defs` (mechanism seen in the `.mod`);
  the rebuild cascade is not measured.
- Vestigial: the `gen_pspecs_and_thumbs` non-worker branch (writes a file nothing
  merges), the `stream` branches of `ctf_estimate` and `pick`, `ctf=flip|no`
  validation in `extract`/`reextract`, two public constants no gate uses
  (`FUZZY_BALL_SIGNAL_HARD_REJECT_MIN`, `SIEVE_BP_CENTER_EDGE_VAR_HARD_REJECT_MIN`),
  unused constants in `simple_defs_stream.f90`, a per-particle `ind_in_stk`
  resolution in `cavger_transf_oridat` that nothing reads, a `stop` and an
  `error stop` in two testers, and a CMake comment saying 5 000 particles for
  `lib_cart_align3D` (code: 1 000).
- Withdrawn: the claim that `set_shmem_flag` has no caller. `solve2D` and
  `make_cavgs` call it.

**Policy and code disagree (the decisions of section 2)** D5 to D7 and D9 are
not defects until the maintainer decides; they are listed there.

**Compile-time rules, counted** (`compile_time_policy.md`; nothing checks them):
8 plain `type(cmdline)` and 2 `type(builder)` components (one added after the
policy, `simple_refine3D_strategy.f90:75`); 4 whole-`parameters` copies (none new);
191 plain `parameters`/`builder`/`cmdline` locals in testers (51 added since
2026-10-04); 768 module-level `use` without `only:` in 481 of 785 files (63 in 59
of the 172 files added since 2026-10-01); `simple_core_module_api` importers 278 to
291.

## 4. Gap register

Priority 1 means a developer could plausibly break a contract the code keeps.
"Home" is where the rule belongs.

| # | Gap | Evidence | Home |
|---|---|---|---|
| G1 | **Project ownership and write budget per workflow.** Only refine3D (section 4) and the general layering policy state it. Missing: the 2D segment table, the solve3D family (solve3D owns the file between stages, each `refine3D` master inside one; `solve3D_cavgs` writes the original once; chunks and pyramid never write it; the add-on rewrites its working project three times), reconstruct3D, the preprocessing steps, the stream | auditors A2, B, C, E | a "Project reads and writes" block in each workflow policy; rules in the layering policy |
| G2 | **Preprocessing-step family.** Four project-I/O patterns, `nparts` required, part-file shapes, keys each step writes on the mic row, failure rules; no CTF-estimation policy (`ctfdoc` format, the `+1+box/2` pixel-origin convention extraction relies on, `forctf` lifecycle, stream thresholds); no pick/extract/reextract contract (box-file format, `state=0` for empty micrographs, the integer-shift convention of reextract); no `refine_motion_model` policy | audit E section 0 and cross-cutting 1-5 | `preprocessing_policy.md` (or a layering-policy section) plus `ctf_estimation_policy.md` |
| G3 | **Queue-system contract.** Completion by sentinel, exit-code file and liveness (`JOB_FINISHED_*`, `EXIT_CODE_JOB_*`, `JOB_LOST_EXIT_CODE`), `stop_on_failed_parts`, cancellation, `qsys_declare_part_finished`, `numlen`, the job record; the coarray and persistent-worker policies both assume it; no skill covers `src/utils/qsys` | audit D cross-cutting 1 | `queue_system_policy.md`, the two backends as deltas |
| G4 | **State weight set and `flex_pca`.** Producer, consumers (`m_estimator=flex`, assembly floor, PCG, trailing chains, `refine3D_auto state=X`, probability table ESS), files, kinds, thresholds, spread over four policies; `flex` has no policy though `check_flex_dag.py` runs in the fast gate | audit A1 cross-cutting 1 | `heterogeneity/state_weight_set_policy.md` |
| G5 | **3D convergence contract.** Single-state, independent multi-state, `input_oris_refine`; `fracsrch`, `overlap`, `OVERLAP_STATE_JOINT`, `scan_exit`, `minits`, `lpexit`, the distributed `iter >= startit+2` floor, `FRAC_SH_LIM`; the 2D counterpart exists | `simple_convergence.f90:455-560` | `refine3D_policy.md` new section |
| G6 | **Who draws the sample, per mode,** and the sampling predicates: `master_draws_sample3D`, the polish/eval/sigma/`update_missing` exceptions; `update_frac >= 0.99` means no sampling; `refine2D_blends_carryover`; `is_fresh_2D_start` | `simple_matcher_smpl_and_lplims.f90:273-355` | `importance_sampling_fractional_update_policy.md` section 1 |
| G7 | **solve3D stage plan as one table** (stages 1-8: backend, `ml_reg`, `trs`, `overlap`/`fracsrch`, `frac_best`, `greedy_sampling`, `trail_rec`, `fillin`, `nspace_sub`), the `state=` continuation route, `nrestarts`/`screen`, the sampling rule with the `nstates` factor, the manifest key tables (32 recorded, 19 replayed), the strip/copy duality of `strip_pcg_backend_keys` | controller `:299-523` | `solve3D_policy.md` |
| G8 | **Optional PCG priors** (soft solvent, atomicity: 744 lines, `pcg_mskfile`) have no contract in any policy; the PCG policy says it has no regularization beyond `P_tau` | `simple_pcg_atomic_prior.f90`, `simple_pcg_solvent_sidecar.f90` | `reconstruct3D_pcg_policy.md` section "Optional priors" |
| G9 | **Stream hand-off.** "Who writes what" (the master owns no project; stage folders; the three exceptions), "Programs the stages run" (program, binary, project in/out, owner during the call), the watcher contract, `.tmp` plus rename for every hand-off, one process per stage folder (no lock exists), the TERM_STREAM marker | audit C cross-cutting 1, 2, 5 | `stream/README.md` |
| G10 | **Particle stack identity.** About 15 sites write `stk`, `nptcls_stk`, ranges, `stkind`, `indstk` by hand; `ptcl_layout_digest` and the cache key hash them with no rule that a change invalidates `sigma2_state` and caches; three project mergers, the policy covers one; readers outside the read session | audit B I2-I4, R2 | `project_stack_indexing_policy.md` |
| G11 | **Compile-time rules have no check.** Proposed `scripts/check_compile_rules.py`, report-only with a baseline, in the fast gate: plain heavy components, whole-`parameters` copies, bare module-level `use` in changed files, umbrella imports in new files, plain heavy tester locals. Plus a rule that no volatile value (commit hash, date) is a `parameter` of an imported module | audit D table of 16 rules | `compile_time_policy.md`, `scripts/` |
| G12 | **Development workflow** (never commit or push, the incremental build loop and parallel install, which tests for which change, master-based working now in test policy 3.4), **Fortran conventions that are enforced** (`!@descr:`, one module per file), and **generated sources and the key vocabulary** (where a new key must be declared) have no policy | audit D cross-cutting 2-4 | `development_workflow_policy.md`; a "Generated sources" section in `ui_layer_policy.md` |
| G13 | **Artifact and lifetime table** of one refine3D iteration (reprojection model, `dist_part*`, `assignment_part.dat`, partial volumes, `pcg_raw_*`, trailing chains, `algndoc_*`, sigma2 files, sentinels, `BENCH/`) and its 2D counterpart; dropped and carried-forward states; `sigma_commit_deferred`; the `update_missing` pass | audit A1 | `refine3D_policy.md` section 11 (replacing the bench table) |
| G14 | **`segdiam_bin_picker`**, the picker the stream uses: first batch decides the diameter bins and box, later batches reuse them, transient keys, deletion of micrographs | audit E section 5 | `segdiam_bin_picker_policy.md` |
| G15 | **Optics-group assignment as one contract** (batch program and stream p02 differ); the meaning of mic `smpd` changes after motion correction; the solve3D manifest digests CTF model and `ogid`; `ori%get_ctfvars` silently returns 0 for absent `kv/cs/fraca` while `get_ctfparams` stops | audit E microscope | `microscope_parameters_policy.md` |
| G16 | **Final-map recipe** (`calc_final_rec`: which route when, sigma bootstrap criteria, the cc-map route without ML regularization, cold PCG floor, sharpening) and the **multi-state controller contract** (what each caller sets, the order of the controller's steps, project hand-over between blocks); the wrapper settings of `classify3D_refs`, `refine3D_states` | audit A1 cross-cutting 5, 6 | a final-map policy; `heterogeneity/` |
| G17 | The relational feature `corr_knn_signal_v1` that every quality-model artifact requires is named and never defined; the callers table of the backend | `simple_cavg_quality_relations.f90` | `model_cavgs_rejection.md` |
| G18 | **Tests.** No policy maps its rules to tests. No gating test for `classify3D_refs`, the particle cache, cache versus no-cache identity, `restore_cavgs_inmem` versus distributed assembly parity, cFAR numerics, `bootstrap_cavgs`, dose weighting and `calc_eer_fraction`, or project reads and writes per iteration. No document maps a code area to the high-level entries a change affects (proposal: a `covers:` field per registry row) | all clusters | each policy gets a "Tests" line |
| G19 | `solve3D_cavgs` restart driver (`nrestarts_collapse`, `exit_collapse`, `SOLVE3D_CAVGS_FINAL_DIR` read by stream p03) and the route's stage constants | audit A2 | `solve3D_cavgs_policy.md` |
| G20 | Re-entrancy of private commanders: under the coarray driver one image runs its partitions in one process, so module state, threads and the random state carry over | audit D coarray | `coarray_execution_policy.md` section 7 |
| G21 | **SINGLE workflows** (`nano3D`, `refine_tomo`, `tseries`, species discovery; 14 modules in `src/main/nano`, six `single_commanders_*` files) have no policy; five policies mention them only in passing | `rg` over `doc/policies` | `single/` |
| G22 | The **sigma2 state** contract (canonical state, commit, deferred commit, project-owned path, bootstrap, add-on copy-on-write, stream pool store) is spread over 14 policies with no owner | `rg sigma2 doc/policies` | `sigma2_policy.md` |
| G23 | Constants and defaults tables that exist only in code: the `fsc_area_score` refinement constants; KB constants and the cropped-box padded grid; `refine3D_states` and `classify3D_refs` wrapper defaults | audits B, A1 | the owning policy |

## 5. Structural debt

Policies describe the software as it is (`for_developers/index.md`); these carry
history, benchmark logs, plans or protocols. Splits are proposed, not made.

| Document | Lines | Problem | Proposal |
|---|---|---|---|
| `test_environment_policy.md` | 1328 | test cards with dated medians (196 lines), history (158), planned work (43), an external workflow (49) | keep about 600; runner and gates policy about 330; cards, history and plans to `code_overview`, `test_review_record.md`, the testing plan |
| `3D/refine3D_policy.md` | 1229 | sections 8 and 11 are mostly implementation notes, bench formats, a changelog | keep about 450; matcher note to `implementation_notes/`; bench reference to `how2s/`; sharpening recipe to the PCG policy (it is stated twice) |
| `3D/reconstruct3D_pcg_policy.md` | 882 | regression stories, review sections, a seed protocol | keep intro, sections 1, 2, 8, 10 and the invariants; the rest to `pcg_decision_log.md` |
| `3D/solve3D_policy.md` | 730 | about 200 lines of dated decisions and trajectories; section 3 duplicates the sampling policy | `implementation_notes/completed/solve3D_stage_decisions.md` |
| `NU/nonuniform_filtering_policy.md` | 675 | about 150 lines of records; section 10 documents a removed feature; two orphan figures | `pcg_decision_log.md`; link or delete the figures |
| `phase_shift_ctf_policy.md` | 625 | about 330 lines of revision log, a future-estimator proposal, a validation protocol with accession tables | keep about 250 |
| `importance_sampling_fractional_update_policy.md` | 573 | sections 5 and 7 duplicate `solve3D_policy.md` and the trailing skill | contracts only |
| `stream/stream_3D_ingestion_policy.md` | 577 | the state-assignment method (sections 9 and 10) is shared with the batch `solve3D_pyramid` | `3D/state_assignment_policy.md` |
| `2D/solve2D_policy.md` | 566 | part-file layouts and algorithms (262-337); carried sums stated in three policies | implementation note; one home for carried sums |
| `ui_layer_policy.md`, `compile_time_policy.md`, `coarray_execution_policy.md`, `solve3D_chunks_policy.md`, `refine3D_auto_policy.md`, `sieving_and_rejection/model_cavgs_rejection.md` | | status sections, an evidence log, design wishes, plan narrative, benchmark narrative, a retired-rule-engine comparison (75 lines), a personal path | move to completed notes |

Duplicates to give one owner: carried class sums (three policies), final-map
sharpening (two), the benchmark inventory (two), the NU ladder (two), the merge
contract (two), first-set triggers (two), the final-ingestion trigger (two).

## 6. Instruction and skill fixes (applied 2026-10-21, with the maintainer's approval)

- A "Layering and project I/O" section was added to `CLAUDE.md` and
  `.github/copilot-instructions.md` (the proposed text is kept in
  [proposed_instruction_edits.md](policy_set_audit_2026-10-21/proposed_instruction_edits.md));
  `CLAUDE.md` now lists the sixth executable, `simple_persistent_worker`.
  `CLAUDE.md`, `AGENTS.md` and `copilot-instructions.md` still repeat the same text
  and can drift; `AGENTS.md` only points to `copilot-instructions.md` and was not changed.
- `simple-architecture`, `simple-main-commanders` and `simple-main-project` point to
  the layering policy and state the worker-frame rule, the missing file lock and the
  single-maker rule.
- `simple-refine3d`, `simple-main-strategies` and `simple-main-project` no longer say
  that assembly and the PCG master write nothing or that every master writes once per
  iteration: the PCG master and the refine2D master are named as the exceptions.
- `simple-main-pick` no longer says `segdiampick_mics` estimates the stream's
  diameter and points to the policy that describes the stream's picker.
- `simple-solve2d` states that the 2D master writes `ptcl2D` more than once per
  iteration. Its `simple_exception` was not wrong: that is the generic interface of
  module `simple_error`, which the skill now names.
- `multi_microscope_parameters.md` lists `phshift` where it listed `phaseplate`
  (the field does not exist; the merge validator requires `phshift`).

## 7. What the audit did not check

Scientific and benchmark numbers (PfCRT, aldolase, beta-gal, msp1, 5.6x, resolution
values); dated rulings individually; bit-identity and byte-identity claims; the
numerics of the continuous in-plane, particle-frame shift and PCG kernels (names
and call sites only); the NICE side of the IPC policy; coarray runtime behaviour;
the external `SIMPLE_data_testing` workflow; the algorithm notes under
`doc/algorithms/`. Counts in section 3 are static greps at `fd5b5bc8e`.
