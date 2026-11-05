# Audit C: stream and sieving_and_rejection policies

Scope: doc/policies/stream/ (README, p00-p07 overviews, ipc, pool2D, reference_generation, restart,
stream_3D_ingestion) and doc/policies/sieving_and_rejection/ (index, class_compatibility,
model_cavgs_rejection, ptcl_sieve). Checked against src/main/stream/{,stages,master,shared,pool2D},
src/main/commanders/stream, src/main/sieve, src/main/cavg_quality, src/main/class, src/main/solve,
the skills simple-main-stream, simple-cluster-cavgs-quality, simple-main-class, and the working tree
as of 2026-10-21. Read-only; nothing built or run.

Baseline result: names, files, routines, constants, test names and relative links in these 16
documents are almost all correct (every `src/...` and `doc/...` path resolves; every routine named
in stream_3D_ingestion_policy exists; the README resource table, the IPC key table and the sieve
sentinel list match the code). No document uses `cluster2D` or `abinitio*` as a living name. The
problems are (a) the pool2D policy does not state how the pool calls refine2D and who owns the
pool project during an iteration, (b) there is no statement anywhere of who owns which stream
project file, (c) model_cavgs_rejection.md is stale on who calls it, (d) a lot of dated history is
mixed into the policies.

Two files in the working tree are new and untracked (not audited): doc/policies/README.md and
doc/policies/layering_and_project_ownership_policy.md. Its section 5.1 ("A stream stage and the 2D
pool each own the files they write. Two processes never write one file.") is the only ownership rule
in the doc set; see the pool findings below for where the code does not follow it.

---

## doc/policies/stream/

### doc/policies/stream/README.md  (lines 88; last commit 2026-10-14; verdict: minor drift)
- [M][GAP] L3-6 claims scope "the watcher and utilities in src/main/stream", but no policy states the
  watcher contract (a file is reported once untouched for `report_time`; history by basename, so two
  movies with one name in different folders collide; `spproj=.true.` mode watches `*.simple`;
  `add2history`; `nretries`). Evidence: src/main/stream/simple_stream_watcher.f90:21-24, 70-80. The only
  place the rule is written is the skill (simple-main-stream Guardrails). Proposal: a "Hand-off between
  stages" section here (see Cross-cutting gap 1).
- [M][GAP] No page lists the stage-to-stage contract in one place: input folder (`dir_target`) -> output
  folder, what is handed on (`spprojs_completed/*.simple`, renamed in from `.tmp`), markers
  (`STREAM_IDLE`, `STREAM_FINISHED`), `moldiam.txt`, `selected_references.mrcs`, optics maps. Pieces
  are spread over ptcl_sieve_policy 6.3, reference_generation 4, restart 4. Proposal: table in README.
- [L][STRUCTURE] L44-45, L63: "decision 14 of the 5 October fix plan", "(decision 10)": plan labels in a
  living policy. Fix: delete the parentheses.
- [L][STRUCTURE] No link to doc/algorithms/streaming_pipeline.md or to the new
  doc/policies/layering_and_project_ownership_policy.md (both exist). Add to "Related policies".
- Verified correct: resource table L50-58 against simple_stream_master_resources.f90:36-52, partition
  variables L60-62 (simple_defs_environment.f90:14-30), worker server L66-74 (p00_master.f90:57-58),
  optics paragraph L76-88 (simple_optics_groups.f90:10-24).

### doc/policies/stream/p00_overview.md  (107; 2026-10-14; current)
- [L][STRUCTURE] L7 "as of 14 October 2026": dated status line. Same on p01 L7, p02 L6, p03 L8, p04 L7,
  p05 L8, p06 L8, p07 L8. Fix: drop the date (the README already says the policies and code rule).
- Content verified against p00_master.f90: fork order p01,p02,p04,p05,p06,p07 then p03 (L168-181),
  5 s heartbeat, 10 ms listener, 60 s optics stop, 600 s kill, 60 s start-up abort, STOP_ONE_BY_ONE.

### doc/policies/stream/p01_overview.md, p02_overview.md, p03_overview.md, p05_overview.md, p07_overview.md  (66/54/89/78/86; 2026-10-14; current)
- Numbers checked against named constants and found correct: 5 movies per set, 16x4 threads, 60 s
  settle, 180 s partial set, 900 s idle, 1000 initial picks, thresholds 10/10/1.0, 1.3 A; optics 50 per
  pass, 5 maps kept, tilt_thres 0.05; p03 100/500 micrographs, 501 rows, 0.8/0.1/1.5 vetoes; p05
  5000/10000 particles, 100 classes, 128 px, 15/10 A, 29 min, 20 sets per pass, prob 0.25; p07
  EO_ANGDIST_MAX 50, EO_OVERLAP_MIN 70, REFS_PAIR_CC_MAX 0.97, FIT_NMAD 3, FIT_REJECT_MAX 0.3,
  MARGIN_AMBIGUOUS 0.02, NSPACE_ASSIGN 1000. No findings other than the dated line.

### doc/policies/stream/p04_overview.md  (79; 2026-10-14; current)
- [L][GAP] L76 "The `lp_pick` setting currently has no effect" is true (no strategy reads it for
  pick_extract) but `lp_pick` is still offered in the pick_extract UI (src/main/ui/simple/simple_ui_stream.f90:394)
  and defaulted by the p04 commander (simple_commanders_stream_p04_refpick_extract.f90:76). An open
  defect sits in an overview. Proposal: list it under a Known gaps bullet of reference_generation_policy
  section 7 (p04's use of references), or remove the parameter.

### doc/policies/stream/p06_overview.md  (75; 2026-10-14; minor drift)
- [M][GAP] L73-75 "Possible gap, not yet tested: a pool already paused past iteration 25 when the sieve's
  empty final set arrives may never run the iteration that carries the final flag" -> confirmed by the code,
  not "possible". `transfer_sets` returns at the empty final set with `nimported = 0`
  (simple_stream_stage_pool2D.f90:597-613), `import_sets` unpauses only inside `nimported > 0`
  (L534-538), and `runs_to_final` is `iter < FINAL_ITER` (L1100-1104), so a pool paused at iteration >= 25
  stays paused and `publishes_final` (L1109) is never reached; p07 then never gets `pool_final=yes`.
  Fix: move to pool2D_policy section 12 as a Known gap with this evidence (or fix: unpause when the final
  set is noted).

### doc/policies/stream/ipc_policy.md  (196; 2026-10-13; minor drift)
- [M][GAP] Section 7 documents SIGTERM as the stop channel. There is a second, file-based one: every
  stage ends its loop when `SIMPLE_TERM_STREAM` exists in its folder (`finished = ... file_exists(TERM_STREAM)`,
  stages/simple_stream_stage_pool2D.f90:350, optics:196, preprocess:334, refpick:352, sieve:293, solve3D:378;
  not p03), and `terminate_stream` makes the refine2D master and cavgassemble exit when `async=yes`
  (simple_stream_utils.f90:32-38, simple_refine2D_strategy.f90:394, 758, simple_cavgassemble_strategy.f90:99).
  Nothing in the master writes it any more; the only writer is the pool's `finalise`
  (simple_stream_pool2D.f90:1159). restart_policy mentions removing a leftover file but nobody says what
  creates it. Proposal: section 7 item 6, "TERM_STREAM is a legacy stop marker: stages poll it, only the pool's
  final write creates it, to bring its refine2D master to a hard stop".
- [L][STRUCTURE] L159 "Added 13 October 2026." Delete.
- Verified: pipe framing constants, frame bound, `clear_stage`, locking, `l_updates` stages (p01,p03,p06,p07),
  the GUI answer table and its drop rules (simple_stream_master_gui_commands.f90:55-170), capacity
  constants, the GUI-metadata import rule (no import reaches src/main), fork-safe memory monitor
  (simple_posix.c:1462), SIGINT ignore (simple_forked_process.f90:118-123), test suite names.

### doc/policies/stream/pool2D_policy.md  (244; 2026-10-14; minor drift, but incomplete on the refine2D contract)
- [H][GAP] The policy never says how the pool runs 2D. The iteration is `prg=refine2D_distr` through
  simple_private_exec (the batch `commander_refine2D`, `refine2D_distr_strategy` master; workers `prg=refine2D`,
  `commander_refine2D_distr_worker`), with the command line built in simple_stream_pool2D.f90:235-281 (`async=yes`,
  `stream2d=yes` "the only place this flag should be turned on", `cls_init=rand`, `objfun=euclid`, `ml_reg=no`,
  `refine=snhc_smpl` only, `nparts`, `nthr`, `frcs`, `refs`, `startit=maxits=iter`, `smpd_crop/box_crop/msk_crop`,
  `lp/gauref/gaufreq/extr_iter/extr_lim/trs` per iteration, `update_frac` only when the user gave it) and, under
  `objfun=euclid` (always, p06 commander sets it), a two-step script `calc_pspec` then `refine2D_distr`
  (L567-575). `stream2d=yes` changes the batch workflow in seven places, none described in any policy:
  FRC-driven low-pass with lpstop floor (simple_matcher_smpl_and_lplims.f90:177-185), carried-sum blend rule
  (same file L350-353), snhc -> snhc_smpl and greedy for unsearched particles
  (simple_strategy2D_matcher.f90:261-269, 369-372), extremal power (simple_strategy2D_alloc.f90:105),
  convergence limits (simple_convergence.f90:180, 195), no progress file (simple_refine2D_strategy.f90:560),
  no cropped-particle sums (simple_matcher_ptcl_batch.f90:318). A developer cleaning "legacy" stream2d code in
  batch refine2D would break the pool silently. Proposal: new section "The refine2D call" after section 6, one
  paragraph per bullet above, plus "a change of a refine2D default or of the `stream2d` branches is a pool change".
- [H][GAP] The pool project hand-over is unstated. The pool writes a sampled project to `refine2D.simple`
  (`POOL_PROJFILE`, simple_defs_stream.f90:21; stacks of the sample, ptcl2D, ptcl3D, cls2D; simple_stream_pool2D.f90:511)
  and keeps a copy `refine2D_input.simple` (L513); the refine2D master then owns the file for the iteration (class-average
  assembly writes cls2D, cls3D, out and ptcl2D, simple_refine2D_strategy.f90:768-771; `distr_finalize_iteration` writes
  ptcl2D again, L449; `distr_finalize_run` writes out, L484) and `sigma2_state.bin`; the pool touches neither until
  `REFINE2D_FINISHED` exists (`update_status` L917, gated by `l_available`), then reads back only `cls2D`, `stk`, `ptcl2D`
  (`update`, L726-735). The only guards are the `l_available` flag, the job's exit status and the liveness check. Section 6.7
  says only "the iteration is a local queued job". Proposal: state the sequence and "the pool never writes `refine2D.simple`
  while `l_available` is false".
- [M][GAP] Section 10 / 6.7 omit what `finalise` does when a stop comes mid-iteration: it touches `TERM_STREAM`,
  `JOB_FINISHED_<part>` and `CAVGASSEMBLE_FINISHED` in the pool folder "tricking the asynchronous master process
  to come to a hard stop" and finally deletes `refine2D.simple` (simple_stream_pool2D.f90:1158-1163, 1174;
  stage finalize cancels first: stage_pool2D.f90:365-366). The stage is writing into the file protocol of the
  job it just cancelled. Contradicts layering_and_project_ownership_policy 5.1 (own files only). Proposal: state
  it as the one allowed cross-write, or remove the hack and rely on `cancel`.
- [M][GAP] Half-sets: the iteration project's `ptcl3D` is a copy of `ptcl2D` with the 2D clustering removed so
  that it carries the halves, because the refine2D master calls `ensure_field_halves` and finds a blank ptcl3D
  inconsistent with odd stack rows (simple_stream_pool2D.f90:487-495, commit af6c578dc; simple_refine2D_strategy.f90:314).
  `eo` returns to the pool with `transfer_2Dparams` (simple_ori.f90:321). The policy mentions `ptcl3D` only for
  the publication and final project (L197-198). Proposal: one sentence in section 6 item 3.
- [M][GAP] Sigma2 is absent: with `objfun=euclid` the pool has a project-owned `sigma2_state.bin`
  (`set_sigma2_state_path`, simple_stream_pool2D.f90:226-228) in its folder, rebuilt by the `calc_pspec` step, and
  deleted on restart with every `*.bin` (`cleanup_root_folder`, simple_stream_refine2D_utils.f90:37). Proposal:
  mention in the refine2D-call section and in restart section 3 (p06 row).
- [L][GAP] The pool takes no part in the particle cache and no cropped sums (stream2d), as solve2D_policy L131,
  L366 and particle_cache_policy L55 say; this policy is silent, so a reader of pool2D_policy would not know.
- [L][DRIFT] L211-215: "The schedule constants are named in the stage (... `OPTICS_ID_DELTA`) and in the pool module
  (`ITERLIM`, `ITERSHIFT`)". `OPTICS_ID_DELTA` and `NPTCLS_FIRST3D` are in src/defs/simple_defs_stream.f90; `NMICS_FIRST3D`
  (stage L96) is missing from the list; `ITERLIM`/`ITERSHIFT` are locals of `update_aln_params`
  (simple_stream_pool2D.f90:786-787). Fix: name the real homes.
- [L][DRIFT] L118-123: the ramp `lp = lpstop + (lpstart-lpstop)*(20-iter)/20`: `update_aln_params` runs before
  `iterate`, so `iter` is the last dispatched iteration and the formula applies to iteration `iter+1`
  (simple_stream_stage_pool2D.f90:776-777). State "for the iteration after `iter`", as the p06 overview's
  "shifts from iteration 6" already does.
- [L][STRUCTURE] History and plan labels in a policy: L62-64 ("until 5 October 2026 the publications started after
  iteration 25..."), L75 ("follow-up plan, decision 7"), L110 "(decision 5)", L141 "(D40)", L146 "(decision 6)",
  L156 "A guard ... was removed on 5 October 2026.", L184-185 ("decided 12 October 2026, before which only the master's
  former `sieve_ini3D=yes`"), L234 ("review G1 closed by ..."). Fix: delete; keep the rule.
- [L][STRUCTURE] Sections 3.2/3.4 repeat stream_3D_ingestion_policy 3.1 (first-set triggers, constants) in full;
  keep one owner (3D ingestion is the contract; pool2D should link).
- Verified: pause table vs `pause_rate_factor` (stage L1082-1090), `runs_to_final`, `exports_after`, `first_set_due`,
  `default_mskdiam`, ring history, sampling, class draws, snapshot rules, constants `MAX_SMPD` 2.67, `STREAM_NPTCLS_MAX`.

### doc/policies/stream/reference_generation_policy.md  (218; 2026-10-13; minor drift)
- [M][STALE] L208-210 "`STATE_CC_MIN_FRAC` and `STATE_POP_FLOOR` (10% each) wait for a validation run; a speck larger
  than 10% of the particle still makes a second object" -> `STATE_CC_MIN_FRAC` no longer exists (only this document and
  stream_followup_plan mention it); the veto is `STATE_DOMINANT_FRAC` = 0.8 (the largest component's share of the
  foreground) plus `STATE_POP_FLOOR` = 0.1 and `STATE_RES_FACTOR` = 1.5 (stage_initial_analysis.f90:156-158,
  state_vetoes L1571-1581), as the table at L76 already says. Fix: "`STATE_DOMINANT_FRAC` (0.8), `STATE_POP_FLOOR`
  (10%) and `STATE_RES_FACTOR` (1.5) wait for a validation run; the fractions are guesses."
- [L][DRIFT] L72 table row `solve2D` omits `cache=yes` (stage_initial_analysis.f90:1346); the p07 route's solve2D
  states it (3D ingestion 9.4). Add.
- [L][STRUCTURE] Plan/decision labels: L40 (follow-up plan, decision 4), L62 (decision 20), L76 (decisions 11-14,
  "the 3 October order ... is gone"), L82 ("decided 8 October 2026"), L196-201 (review ids M1, M2, M5, R2, R4 "done on
  5 October 2026"). Fix: delete, keep the rule.
- [L][STRUCTURE] L51-59 (final-ingestion trigger and the `STREAM_IDLE`/`STREAM_FINISHED` wording) duplicates
  ptcl_sieve_policy 6.3 item 1; keep it in ptcl_sieve_policy and link.
- Verified: routes, plan sizes, state vetoes, balancing, mask-diameter rules, output contract names, `make_pickrefs` settings.

### doc/policies/stream/restart_policy.md  (141; 2026-10-14; minor drift)
- [M][GAP] There is no single-writer rule per stage folder. Section 1 item 3 allows "starting one stage on its own, in
  its folder or with `dir_exec`" while the master's copy may be running; nothing (no lock or pid file) stops two
  processes from writing one stage folder. The master only checks `is_running` before a GUI restart
  (p00_master.f90:351). pool2D_policy 12 records "one pool per folder" for the pool alone. Proposal: state "one
  process per stage folder at a time; the master's `is_running` is the only guard".
- [L][GAP] Section 4 / invariant 1 ("hand-offs are atomic") should record the optics map's two-file order: `optics_map_<id>.simple`
  is written first (temporary + rename), the `.txt` renamed last, and readers find a map by its `.txt`
  (simple_sp_project_io.f90:1298-1311, simple_stream_utils.f90:130); reordering breaks readers.
- [L][DRIFT] L21 "`dir_exec=<folder>` (p04, p06)": p01 (stage_preprocess.f90:229) and p07 (stage_solve3D.f90:306) take it too.
- [L][STRUCTURE] L109, L128 "(follow-up plan, decisions 28-30)", "decision 28".
- Verified: every file name, constant and function named in sections 3-6 exists (`gui_thresholds.txt`, `state_refs.txt`,
  `rejected_publications.txt`, `JOB_LIVENESS_S` 300, `JOB_LOST_EXIT_CODE` 254, `cleanup_root_folder`, `restore_export_id`,
  `fresh_job_dir`, `cancel_unfinished_jobs`, RESTART headers in all seven stages); TERM_STREAM removal in p01,p02,p04,p05,p06,p07
  and none in p03; markers removed at stage start; test names.

### doc/policies/stream/stream_3D_ingestion_policy.md  (577; 2026-10-20; minor drift, history-heavy)
- [M][STRUCTURE] Changelog mixed into the contract: L25 ("since 15 October 2026"), L34-36 (inputs that "went with them"), L76-78,
  L124-127 (Until 7 October / on 8 October / from 12 to 13 October), L200-203 (item 11 "Retired on 13 October 2026", wholly
  history), L220 ("plan, decision 5"), L226, L345-352, L392, L396, L422 and L427 ("50 degrees since 14 October 2026, 20 before";
  "the run of 13 October 2026 that collapsed"), L474-477. Fix: delete 4.11; state values without "since/until".
- [M][STRUCTURE] Sections 9 (items 4-8) and 10 (items 3-6) and their constants are the shared method of
  `simple_state_assignment` and of the batch program `solve3D_pyramid`; solve3D_pyramid_policy.md L4-7 points back to these
  sections, so a 3D policy lives in the stream folder. Propose doc/policies/3D/state_assignment_policy.md for the method
  (check, matching, fit rules, local refinement seeds) and keep here the stage's phases, retries, publications and GUI. Also
  brings the file from 577 dense lines to ~350.
- [L][GAP] 9.2 says the first publication keeps "the sieve's 2D parameters and selection"; `build_sieve_publication` also
  deselects (state 0, class 0) every particle whose class is outside its set's classes
  (simple_stream_refine2D_utils.f90:171-176). Add one clause.
- Verified: ~70 routine names, all constants (`NCLS2D_FIRST` 150, `NCLS_EXPANDED` 500, `ROUTE_TRIES` 2, `LP_ASSIGN` 15, `MIN_PTCLS_PER_STATE` 5,
  `NQUALITY_KEPT` 3, `NPUBLICATIONS_KEPT` 2, minpop 25, athres_cont 10), commander settings (nstates 2..20, nparts3D/nthr3D 8), job programs
  and exec binaries (simple_exec for solve2D, cls_expansion, solve3D_cavgs, refine3D; simple_private_exec for assign_cavgs_states),
  `write_stage_project` (tempfile), `eo_transaction` seeds, all test names.

---

## doc/policies/sieving_and_rejection/

### doc/policies/sieving_and_rejection/index.md  (19; last commit 2026-10-05, working tree modified; verdict: minor drift)
- [L][STRUCTURE] L11 links a plan. The committed version points at implementation_notes/planned/distance_transform_shape_rejection_plan.md,
  which does not exist; the working tree already repoints it to implementation_notes/rejected/ (exists). Keep that fix; better, drop the
  bullet (plans are not policy) or label it "rejected".
- [L][STRUCTURE] No "read before changing" column (the stream README has one) and no link to the stream README; the index lists three
  documents but not that p03/p05 call the sieve and that p03/p07 call the quality backend.

### doc/policies/sieving_and_rejection/ptcl_sieve_policy.md  (432; 2026-10-06; minor drift)
- [M][DRIFT] L141-151 "`cycle(project_list)` must execute in this order: collect_and_reject, generate_chunks_coarse, generate_chunks_fine
  when not coarse-only, submit ... must not be rearranged" -> the code has a fifth step `if( final_ingestion ) hand_off_final_set(project_list)`
  after `submit` (simple_ptcl_sieve.f90:603-607), and 6.3 item 5 depends on it ("collect_and_reject runs before it"). Fix: add step 5.
- [L][DRIFT] L21-22 "not responsible for stream file watching/import logic (owned by stream commanders)" -> owned by the stage types
  `stream_stage_sieve` and `stream_stage_initial_analysis` (src/main/stream/stages); the commanders are thin.
- [L][DRIFT] L43 "the standalone `sieve` commander" -> the CLI program is `particle_sieving` (`commander_sieve_ptcls`,
  simple_commanders_sieve.f90:14; exec_sieve:21; pre_chunked mode, `generate_sieve_projects` in simple_ptcl_sieve_utils.f90). Neither is described anywhere.
- [L][DRIFT] L254, L270-271: `evaluate_cavg_quality_hard_reject` / `evaluate_cavg_quality`; the sieve calls `score_project_cavgs`
  (simple_cavg_quality_selection.f90:35), which calls both (simple_ptcl_sieve.f90:1616-1618). Fix the names.
- [L][DRIFT] L242-245: partition variable. Only `SIMPLE_STREAM_CHUNK_PARTITION` is named; p03's sieve chunks use `SIMPLE_STREAM_REFGEN_PARTITION`
  (stage_initial_analysis.f90:714-718), as the stream README says.
- [L][GAP] Section 10 omits `combine_completed_chunks(..., with_sigma2)`; p03 combines with `with_sigma2=.false.` (stage_initial_analysis.f90:838).
- [L][STRUCTURE] L38 "decision 12", L58-62 ("Before 5 October 2026 scoring always used the box's disc ... need a validation run before release"),
  L265-266 ("stream fix plan, decision 8, confirmed 5 October 2026"), L398-399 ("follow-up plan, decision 29"). Fix: delete.
- Verified: public methods, `new` signature, sentinel and list files, tier defaults, `get_finished` logic, hand-off by temp + rename,
  retry/drop rules, reasons `<tier>_reject:<reason>` with underscores, counters, tester name.

### doc/policies/sieving_and_rejection/class_compatibility_policy.md  (178; 2026-10-05; minor drift)
- [L][DRIFT] L113-114 "rejection reason `class_compatibility: size_incompatible_subset`" -> the code string has no blank:
  `class_compatibility:size_incompatible_subset` (simple_class_compatibility.f90:219), consistent with the "underscores for blanks"
  rule in ptcl_sieve_policy 8. Fix the string.
- [L][GAP] Scope lists no callers. Two exist with different lifetimes: the sieve (one persistent model per tier, trained until converged,
  simple_ptcl_sieve.f90:1633, 1654) and p03 (a fresh model per class selection, trained and applied on the same selection,
  stage_initial_analysis.f90:1438). Add a "Callers" bullet.
- Verified: API, metrics fields, convergence constants (`CONV_ABS_EPS` 0.5, `CONV_REL_EPS` 0.01), `nset < 3` rule, preprocessing steps.

### doc/policies/sieving_and_rejection/model_cavgs_rejection.md  (475; 2026-10-05; stale in parts)
- [M][STALE] L22 "stream initial analysis in the `chunk` context (`simple_stream_p03_initial_analysis`), and the multistate stream stage
  runs `model_cavgs_rejection`" -> p03 is `simple_stream_stage_initial_analysis` (chunk100mics via `score_project_cavgs`,
  stage L1459-1460); p07 does not run the program: it calls `score_cavgs_with_model` (simple_state_assignment.f90:384) with the `pool`
  preset, for the route's class selection (applied) and once per publication (logged only, stage_solve3D.f90:563), and the batch
  `solve3D_pyramid` does the same. CONTRADICTS the skill simple-cluster-cavgs-quality "Other Callers", which is right. Fix: replace the sentence
  with a callers table (Cross-cutting gap 4).
- [M][DRIFT] L58-64 the three "workflow phases" (sieve; chunk = "sieve-cleaned particles processed in larger chunks, typically 10-30k"; pool =
  "final 2D phase before 3D") describe a pipeline the stream no longer has. Today: `sieve` = both tiers of the particle sieve (p03 and p05);
  `chunk` = p03's class selection and the CLI default; `pool` = p07's route and per-publication scoring and solve3D_pyramid. There is no 10-30k
  chunk stage. The code comment in simple_cavg_quality_types.f90:49-52 has the same drift.
- [M][GAP] The integration point is not in the document: `score_project_cavgs` / `write_cavg_selection_stacks` in
  simple_cavg_quality_selection.f90 (shared by sieve, p03, p07, the commander) and its rule that `smpd` must be passed so the relational feature's mask
  radius is right ("harmless only with mskdiam = 0", L35-41). L11-18 also lists 6 of the 9 modules; `_helpers`, `_relations`, `_selection` are missing.
- [M][GAP] The relational feature required by every artifact (`corr_knn_signal_v1`, feature `signal_stats_anchor_topk_mean`, knn 5, hp 100 A, lp 15 A,
  trs 10 px; simple_cavg_quality_relations.f90, types.f90:35-39) is named in L37, L137, L169-188, L352 but never defined: what is correlated, what
  is aggregated, how the coefficient enters. Add a "Relational feature" section.
- [M][STRUCTURE] L20-96 and L87-96: about 75 lines compare the model with a retired rule engine (`simple_cluster2D_rejector`, `simple_microchunked2D`, both
  removed), with a "Retired microchunk tier thresholds" table that no code uses. History, not contract. Fix: delete or move to
  doc/refactoring_notes/completed.
- [M][STRUCTURE] Training-set and benchmark content: L163, L165, L289, L293 name `/Users/elmlundho/model_cavgs_rejection/chunk_training5` and `pool_training4`
  (a personal path); L295-301 a stray config block and "macro_evaluate_score=0.50847 ... tp=297, fp=95, tn=151, fn=33". Fix: delete; record the training
  set in the model spec if needed.
- [L][OLDNAME] L22, L27 `simple_cluster2D_rejector`, `simple_microchunked2D`, `simple_stream_p03_initial_analysis`: historical or renamed module names.
- [L][DRIFT] L83, L253 "non-positive" local variance: the gate is `<= EPS` with EPS 1e-6 for both (simple_cavg_quality_feats.f90:299-300). L225-239
  omits `ANALYSIS_BOXSIZE` = 128 (features are measured at a 128-px working box, feats.f90:68, 451). `FUZZY_BALL_SIGNAL_HARD_REJECT_MIN` and
  `SIEVE_BP_CENTER_EDGE_VAR_HARD_REJECT_MIN` are public constants that no gate uses (dead; delete in code, the doc is right to omit them).
- [L][GAP] No Tests section; the backend is covered by `unit_numerics` sub-suites "cavg quality relations" and "cavg quality selection"
  (simple_commanders_test_class.f90:358-359). Add, as the other two sieving policies have.
- Verified: modes and file names written by `apply`/`analyze`/`learn`/`evaluate`/`promote`, `score_states` rules, model_version 11, built-in thresholds
  (chunk 0.35, pool 0.30, sieve 0.25), hard gate tables per context, defaults `mkdir=yes`/`prune=no`/`oritype=cls2D`.

---

## Cross-cutting gaps in this cluster

1. **Who owns which stream project file, and the single-writer rule: stated nowhere in the stream policies.** As built (verified):
   the master p00 owns no `sp_project` at all (p00_master.f90 has no project I/O); it forks the stages with `projfile=<job>.simple`,
   `outdir=<job>` and `dir_target=<upstream>` (make_stage_clines L447-589). Each stage owns its folder (names in simple_defs_stream.f90:11-39) and
   the stage project in it (written at finalize: p01 stage L1020-1021, p02 L202/321, p04 L846-853, p07 `write_stage_project` tempfile, p06 final project
   in `finalise`). Stages hand on only by `spprojs_completed/*.simple` written as `.tmp` and renamed (p01/p04 job sets, p05 `hand_off`, p06 publications)
   and by single files `selected_references.mrcs`, `moldiam.txt`, optics maps. Downstream stages never write upstream folders. Exceptions: (a) the master
   symlinks the user's earlier preprocessing folder and touches `STREAM_FINISHED` in it (p00_master.f90:166-181), (b) p01/p04 rewrite the `mic` segment of a
   finished job's project in the job folder before moving it (p01 stage L803), (c) p02 writes a stray project in the launch folder before `params%new` copies it
   (stage_optics.f90:126-132), (d) the pool and its refine2D master share `refine2D.simple` in turn, and the pool's final write puts the job-protocol markers in
   that folder (pool2D findings above). Protection against two writers is by sequencing only (exit status files, `REFINE2D_FINISHED`, the liveness check, set-aside
   `_unfinished<k>` folders on restart, `is_running` before a GUI restart); there is no lock. Proposal: a README section "Who writes what" with this table, link from
   layering_and_project_ownership_policy 5.1 and name the pool and the three exceptions there.
2. **The stream's calls into the 2D/3D workflows have no home.** The facts exist only in code and the skill: p03 `solve2D` (cache=yes, nparts=1 default) and
   `solve3D_cavgs`, `reproject`, `extract` through simple_exec; p05 sieve chunks `solve2D` through simple_exec; p06 `calc_pspec` + `refine2D_distr` through
   simple_private_exec; p07 `solve2D`, `cls_expansion`, `solve3D_cavgs`, `refine3D refine=cont maxits=1 nparts` through simple_exec and `assign_cavgs_states` through
   simple_private_exec. For each: program, exec binary, project it reads, project/segments it returns, who owns the project during the call, and the
   keys the stage must set. Proposal: README table "Programs the stages run", and the long form for the pool in pool2D_policy. None of the stream policies names a moved
   function (`volassemble_exec`, `cavgassemble_exec`, `prob_align2D`, `commit_sigma2_exec`) or `execute_from_project`/`spproj_src`; the restructure left them untouched,
   because the stream only calls programs.
3. **Dated history in living policies.** About 45 occurrences (month-day dates, "decision N", "follow-up plan", "fix plan", "review Gn/Mn", "until/since ...").
   Per doc/for_developers/index.md policies describe the software as it is. Proposal: one cleanup pass; history stays in doc/refactoring_notes/completed.
4. **Quality-model contexts and callers table** (belongs in model_cavgs_rejection.md, linked from ptcl_sieve_policy 8 and reference_generation 3):
   sieve preset+context: ptcl_sieve both tiers (p03 coarse-only, p05 two tier; `use_model=no` -> sieve hard gates only); chunk100mics+chunk: p03 cycles 1 and 2;
   pool: p07 route class selection (applied to the first set), p07 per-publication scoring (logged), solve3D_pyramid; plus the program. State once that the fine tier
   keeps `sieve`.
5. **Stage hand-off and watcher contract** (see README findings): untouched-for-`report_time`, basename history, `.tmp` + rename for everything handed on, marker
   semantics, optics map two-file order. Currently split between the skill, ptcl_sieve 6.3, restart 5.1.
6. **Particle cache in stream jobs.** p03's two solve2D and the p07 route solve2D run with `cache=yes`; the sieve chunks and the pool do not (pool: `stream2d` excludes it
   and the cropped sums). Neither the stream policies nor particle_cache_policy tie these together (the latter only says "streaming 2D variants" are not consumers).
7. **Code smells found on the way (not doc errors):** the master forwards `optics_dir` to p03 (p00_master.f90:494) which never reads it; p06 stops with the banner
   "SIMPLE_STREAM_SOLVE2D NORMAL STOP" (commanders/stream/simple_commanders_stream_p06_pool2D.f90:60); simple_defs_stream.f90 keeps many constants no stage reads (`POOL_FREQ_REJECTION`,
   `CLS_REJECT_STD`, `FLUSH_TIMELIMIT`, `POOL_DISTR_EXEC_FNAME`, `SIEVING_*`).

## Unverified

- Algorithm internals quoted in the p01, p03, p04 overviews (patch size ~200 A, defocus range, 20-500 A diameter bins, 3 MAD, coarse 4 A/px and 2 A/px
  refinement, Otsu refinement, 1% usable area, balancing arithmetic, the "astigmatism is relative" and "CTF fit worse than 30 A recorded at the threshold" remarks):
  they live in picking, motion and CTF code outside this cluster and were not traced.
- NICE side of ipc_policy (keeps `restart_<key>` until the stage runs; marks a stream failed after a minute without a heartbeat): not read (nice/).
- Runtime behaviour of liveness checks and scheduler commands (`squeue`, `bjobs`, `qstat`): checked by name and constants only.
- Statistics and training-table statements in model_cavgs_rejection.md (AUC, macro scores) are not checkable from code.
- Tests: names exist and are registered (simple_commanders_test_class.f90 suites), none was run; the registry file
  simple_test_registry.f90 holds only `mini_stream_*` and `stream_preproc` entries, the stream suites are `suites_stream`/`suites_lib_stream` in
  simple_commanders_test_class.f90.
- Claims about the brief (suspect): (a) for the 2D master "writes the particle field once before the workers and once per iteration" does not hold.
  `refine2D_distr_strategy` writes nothing before scheduling (distr_initialize / distr_execute_iteration), and per iteration it writes ptcl2D twice when
  `restore_cavgs=yes` (the pool's case): in `run_distributed_cavg_assembly` (simple_refine2D_strategy.f90:771) and again in `distr_finalize_iteration`
  (L449). The structure2D report assigns the sentence to the 3D master (Phase 5). The pool relies on its own write before the call. This also sits
  against layering_and_project_ownership_policy 5.4 ("one write of a segment per iteration (ruled)"). (b) "workers prg=refine2D": correct for
  simple_private_exec (`refine2D` = `commander_refine2D_distr_worker`, `refine2D_distr` = the master), the reverse of the public simple_exec, where
  `refine2D` is the master.
