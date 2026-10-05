# Stream follow-up plan, 5 October 2026

Follow-up to the 5 October stream review
(`doc/refactoring_notes/completed/stream_area_review_2026-10-05.md`) after its fix plan
(`stream_fix_plan_2026-10-05.md`) landed: D41 ("smaller defects, by area") and the review items
that plan left out. Each item was checked against the current code. Short
names (M1, P1, ...) are this plan's; file references are to the current tree. The review's
methodology items M6 (p03's state choice) and M8 (p07's addon route), its defect D16 (p02's
regrouping) and its architecture items G12 to G15 are here too, in their own sections and
workstreams G to L; they are the review's items of those names, not D41 items.

## 1. Status of each D41 item

| Id | Area | Defect (review wording, shortened) | Status now | Where |
|---|---|---|---|---|
| M1 | master | Ctrl-C kills every stage at once: stages share the master's process group and reset SIGINT to the default | **applies** | `simple_forked_process.f90`:92-96 |
| M2 | master | nothing is flushed before `fork`, so start-up lines can appear once per stage | **applies** (suspected) | `simple_forked_process.f90`:84-90 |
| M3 | master | the store copies a frame into its tag's type with no length check | fixed (WS10, `frame_fits`) | |
| M4a | master | a start-up failure after the forks leaves the stages and the worker server running | **applies** | p00:193, 197 |
| M4b | master | `failtime` rewritten on every poll | **applies** (`stoptime` is now set once) | `simple_forked_process.f90`:219-221 |
| M4c | master | no heartbeat during the 60 s wait for optics assignment to stop | **applies** | p00 `wait_for_stop` |
| M4d | master | `force_stop` kills only the stage process; its queued jobs keep running | **applies** | `simple_stream_master_stage.f90`:160-165 |
| P1 | p01/p06 | the raw-project path reads the never-updated root stub `optics_assignment.simple` and overwrites the groups just applied; p01's call does nothing | **applies, latent**: p06 never reaches the raw path today (O2) | `simple_stream_refine2D_utils.f90`:247-249; `simple_starproject_stream.f90`:419-423; p01 stage:961 |
| P2 | p01 | the movie rate spikes after a restart: the restored history counts as new movies | **applies** | `simple_stream_watcher.f90`:171-181 |
| P3a | p01 | a failed job counts as one movie, not five | **applies** | p01 stage:707, 869-870 |
| P3b | p02 | reports accepted micrographs as imported | **applies** | p02 stage:262-271 |
| P3c | p01 | the last one to four movies of a session are never submitted; preprocessing then goes idle with them unprocessed | **applies** | p01 stage:626 |
| P4 | p02 | writes its project to the bare basename of `projinfo`, and makes the root stub before `params%new` | **applies** | p02 stage:110-113, 171, 239 |
| R1 | p03 | no accepted diameter bin leaves the box at 0, and extraction runs with it | **applies** (suspected) | p03 stage:715-724 |
| R2 | p04 | a forced `box_extract` below 128 is raised to 128 silently | **applies** | p04 stage:271-272 |
| R3 | picking | the template low-pass clamp is inverted: `min(max(30, x), 15)` is always 15 Å | **applies** | `simple_pick_strategy.f90`:199, 333 |
| R4 | p04 | `make_pickrefs` rewrites `moldiam.txt` and the templates in place | fixed (WS5: local job, renamed into place) | |
| R5a | p03 | `params_sieve%nmics = 100` is never read | fixed (WS6: sieve settings) | |
| R5b | p03 | the 500-micrograph cap is checked once per pass | **applies** | p03 stage:536 |
| R5c | p03 | cycle 2's GUI status shows cycle 1's counts | **applies** | p03 `send_opening2D_status`, `send_picking_status` |
| S1a | sieve | coarse chunks are marked complete before the fine project is written | **applies** | `simple_ptcl_sieve.f90` `generate_chunks_fine` |
| S1b | sieve | a coarse chunk's project is written before `chunked_mics.txt`: a crash between them imports its micrographs twice after a restart | **applies** (now with `chunked_mics.txt`) | `generate_chunks_coarse` |
| S2 | sieve | new chunk ids come from the array size, not the highest folder: a skipped folder makes the next chunk overwrite one | **applies** | `simple_ptcl_sieve.f90`:760, 806, 849, 918, 983 |
| S3a | sieve | its queue ignores the chunk partition | fixed (WS6/WS7) | |
| S3b | sieve | reusing a worker server turns its autoscaling (warm-up cooldown) off | **applies** | `simple_qsys_env.f90` reuse branch |
| S4a | sieve | counters restart at 0 on restart (the policy calls them cumulative) | **applies** | `import_existing_chunks_*` |
| S4b | sieve | chunks with nothing selected are never counted | **applies** | `collect_and_reject` |
| S4c | sieve | `kill` leaves tuning defaults and previews in place | **partly**: defaults and mask are reset; previews are not | `ptcl_sieve%kill` |
| S4d | sieve | rejection reasons with blanks are truncated on write | **applies** (as read) | `reject_cavgs` (`'coarse_reject: '//...`) |
| O1 | p06 | `stepwise=yes` imports one set per iteration past the threshold | fixed (WS6) | |
| O2 | p06 | no final project if p06 stops before its first complete iteration, against `pool2D_policy.md` §10 | **applies** | p06 `finalize`; `terminate_stream2D` gets no project list |
| O3 | p06 | a snapshot of the current iteration registers `frcs.bin` before checking it exists, which stops p06 | **applies** | `simple_stream_refine2D_utils.f90`:393-396 (`add_frcs2os_out` requires the file) |
| O4 | p06 | publications carry no optics table and a raw copy of `ptcl2D` as `ptcl3D` | **applies**: the optics table is copied from the pool, which never has one | `build_pool_publication`:535, 561 |
| O5a | p06 | snapshot requests before the pool starts go unanswered | **applies** | p06 `apply_gui_updates` |
| O5b | p06 | a command-line `mskdiam` is replaced by the sieve's at iteration 10 | **applies** | p06 `apply_final_mskdiam` |
| O5c | p06 | resolution limits stay on the first mask diameter | fixed (WS3) | |
| O5d | p06 | stopping mid-iteration pairs iteration N-1's class averages with a live `frcs.bin` | **applies** (suspected) | `terminate_stream2D` |
| T1 | p07 | `nstates=1` or above 20 stops the stage | **applies** | p07 commander; metadata `set_state_stats` |
| T2a | p07 | an addon is reported as "running refine3D" | fixed (status "running solve3D_addon") | |
| T2b | p07 | `particles_imported` counts deselected rows | **applies** | p07 `send_status` |
| T2c | master | old vol3D entries stay in the master's store after a p07 restart | **applies** | `simple_stream_master_meta_store.f90` |
| T3 | p07 | `finish_run` adopts the job's project wholesale (its `projinfo`) and carries no optics table | **applies** | p07 `finish_run` |
| T4a | p07 | a malformed publication stops p07, and a restart takes the same one again | **applies** | p07 `merge_publication` |
| T4b | p07 | the mask diameter is read once, from the first export taken | **applies** | p07 `read_mskdiam` |

Fixed: 8 (M3, R4, R5a, S3a, O1, O5c, T2a, and S4c in part). Still to fix: 36.

## 1a. Review M6: p03's state choice

The review: the 3 October state choice can prefer a compact junk volume. WS2 of the fix plan
(decision 4) did what the review proposed. Components are counted inside the mask and above
`STATE_CC_MIN_FRAC` (10%) of the largest; one component is a veto, not the first key; there is a
population floor (`STATE_POP_FLOOR`, 10%); then distinct projection directions, then population
decide (`choose_state`, p03 stage:1417; `vol_shape_descr`, `simple_image_bin.f90`:564). Four
gaps remain, and with them the failure the review describes:

| Id | Gap | Where |
|---|---|---|
| G1 | The veto needs exactly one component (`nccs == 1`). An elongated or multi-domain particle that the strict second threshold splits in two fails it, while a compact junk state (one blob, at least 10% of the population) passes and wins outright. The tester pins this: a one-component state with one direction beats two states with nine (`test_choose_state`, first case). | p03 stage:1423 |
| G2 | When no state passes, the fallback is the 3 October order, fewest components first: the order that prefers compact junk. | p03 stage:1437-1450 |
| G3 | Both Otsu thresholds are computed over the whole box, solvent included; the mask is applied after thresholding. | `simple_image_bin.f90`:586-603 |
| G4 | Nothing tests for structure: a featureless state that passes is ranked on view coverage alone, which junk can have too. The evidence exists: `solve3D_cavgs` writes a per-state FSC (`refine3D_fsc_fname(state)`) with its final reconstruction. | p03 `finish_solve3D` |
| G5 | The 10% and 10% thresholds were never validated on known datasets. | p03 stage:144-145 |

## 1b. Review D16: p02 regroups every micrograph on every pass

Still applies in full. Each pass that imports anything calls `assign_optics_groups` over every
micrograph (p02 stage:235), which clusters each tilt group with `h_clust`
(`simple_starproject_utils.f90`:200-303). That allocates an N×N distance matrix (1.6 GB at 20,000
micrographs, 10 GB at 50,000). Without `dir_meta`, p01 sets every shift to 0 (p01 stage:665-666),
all centroids coincide, and the merge step costs on the order of N³. p02 checks for no stop
request while it clusters. Two further effects found in this check:

| Id | Defect | Where |
|---|---|---|
| H1 | Regrouping from scratch can move an old micrograph to another group, and renumber groups, between two maps. The readers (sieve, pool, p04) always apply the newest map. | `simple_optics_groups.f90` |
| H2 | After a restart, p02 re-imports at most 50 projects per pass (`MAX_PROJECTS_IMPORT`) but publishes a map every pass. A reader moves every micrograph the map does not list into group 1 (`import_optics_map` fills unlisted ones with 1), so during the catch-up the later stages regroup most of the session into group 1. | p02 stage:57, 242; `simple_sp_project_optics.f90`:39-48 |

## 1c. Review M8: p07's addon route

The review: the addon route has no cadence, rebase or verdict rule. WS2 of the fix plan added two
of them: the cadence (decision 1: an addon run starts when the cohort reaches
max(`MIN_PTCLS_PER_STATE` × nstates, `ADDON_COHORT_FRAC` (10%) of the frozen particles), which
bounds the total cost) and the verdict (decision 2: `solve3D_addon_report.txt` is read, logged
and shown in the GUI status). Rebasing was left a known gap (decision 3). Two parts still apply:

| Id | Gap | Where |
|---|---|---|
| I1 | A REGRESSED addon still becomes the base of every later run. The ingestion policy calls rejecting it "the user's call", but p07 takes no GUI updates (the master forwards them to p01, p03 and p06 only), so the user cannot make it. | p07 `finish_run` (stage:707-736), `read_addon_verdict` |
| I2 | Frozen particles keep the poses and states they got against the maps of their time, and later runs inherit the first `solve3D`'s ladder; nothing ever realigns them, not even at the end of a session. | `stream_3D_ingestion_policy.md` §4 item 9, §7 |

## 1d. Review G12 to G15: compile time, fork, the NICE protocol, queued jobs

| Id | Review item | Status now |
|---|---|---|
| G12 | Compile-time policy: heavy components plain, `parameters` snapshots, testers with plain locals, p00's module-level imports | **applies**, except the snapshots (p03's sieve copy went with WS6's settings type; `estimate_mskdiam` no longer declares a `parameters`). Plain components: `cmdline` in p01 (`cline_exec`), p04 (`cline_exec`, `cline_pickrefs`), `forked_process`, `stream_master_stage` (which holds a second copy in its fork) and the sieve's `chunk2D_state`; `sp_project` in every stage (three in p03); `qsys_env` in p01, p03 (two), p04 (two), p05, p07 and `ptcl_sieve`. Plain test locals: the master tester's `stream_master_stage`, p03's `test_choose_state` stage, the sieve tester's `parameters` and `cmdline` (about 15). p00 imports the seven stage commanders at module level for `stage_commander`. |
| G13 | Fork without exec | **partly**: WS0 makes a child forget the inherited worker server; workstream A of this plan will flush before fork and have stages ignore SIGINT. A restarted stage is still forked from a process with running threads (`restart_policy.md` §7). **Deferred** (decision 26). |
| G14 | The update protocol with NICE is undefined | **mostly fixed** by WS0: the master forwards only changes per stage, nothing to a stopping stage, and acts on a restart flag once (`ipc_policy.md` §5, §6). Left: no written contract on NICE's side, and NICE still sends two keys SIMPLE ignores (`ref_selection`, `increase_nmics`: `nice/nice_lite/data_structures/streamjob.py`:159-160, 307-308, 447). |
| G15 | Queued jobs have no common lifecycle | **partly**: WS8 gave every job a record, a cancel, an exit code and a fresh folder on restart. Left: three ways of running a job (`qsys_async_job` in p03, p07 and p04's `make_pickrefs`; job sets with marker files in p01 and p04, `simple_stream_job_sets`; the sieve's chunks and the pool iteration, markers with exit codes); no liveness (a job killed without an exit code is waited for for good, fix plan decision 26); a local cancel signals the script and its direct children only (`pkill -P`, `simple_qsys_job_record.f90`:67), so a distributed program's local part jobs outlive it. |

## 2. Workstreams

### A. Master and fork (M1, M2, M4a-d). Small to medium

- **M1** (decision 1): each stage ignores SIGINT after the fork (`simple_forked_process`), and
  the master stops the stages in order on Ctrl-C, as `ipc_policy.md` promises.
- **M2:** flush `logfhandle` (and `OUTPUT_UNIT`) before every `fork`.
- **M4a:** a start-up failure after the forks stops the forked stages (request, wait, kill) and
  the worker server before the master stops.
- **M4b:** `failtime` is set once, when the failure is first seen.
- **M4c:** `wait_for_stop` sends the heartbeat while it waits.
- **M4d** (decision 2): after `force_stop`, the master cancels the jobs in the killed stage's
  job records.
- Tests: the fork tester (SIGINT ignored by a child; flush), the master tester (start-up failure
  cleanup with a no-op stage).

### B. Preprocessing and optics counts and files (P2, P3a-c, P4). Small

- **P2:** the watcher's rate starts after the restored history (`raten` set to `n_history` on
  the first watch, or by a `reset_rate` the stage calls after restoring).
- **P3a:** a failed job counts the movies of its set.
- **P3b:** p02 reports every imported micrograph as imported, and the accepted ones as assigned.
- **P3c** (decision 3): the movies short of a set are submitted as a partial set once no new
  movie has arrived for a quiet period.
- **P4:** p02 writes to `params%projfile`, and makes no project in the master's folder before
  `params%new`.
- Tests: p01 tester (partial set, failed-job count, rate after restore), p02 tester (counts,
  project path).

### C. Picking (R1, R2, R3, R5b, R5c). Small

- **R1** (decision 4): with no accepted bin, p03 warns, skips extraction and picks again with
  more micrographs.
- **R2** (decision 5): a forced `box_extract` is used as given, with a warning below 128.
- **R3** (decision 6): the template low-pass is clamped to [15, 30] Å (both places in
  `simple_pick_strategy.f90`).
- **R5b:** picking of the "all" set stops within the pass once `NMICS_PLAN(2)` micrographs are
  picked.
- **R5c:** cycle 2's status reports cycle 2's project.
- Tests: p03 tester (zero box, the cap within a pass, cycle 2 status), the pick strategy's clamp.

### D. Sieve integrity (S1a, S1b, S2, S3b, S4a-d). Medium

- **S1** (decision 10): the writes are reordered so the completion marker comes last, and
  `chunked_mics.txt` is rebuilt from the complete chunks' projects on restart.
- **S2:** a new chunk id is one past the highest chunk folder on disk.
- **S3b:** reusing a worker server never turns its warm-up cooldown off (only a streaming owner
  turns it on).
- **S4a:** the counters are rebuilt on restart from the complete chunks' projects.
- **S4b:** a chunk with nothing selected counts its particles as rejected.
- **S4c:** `kill` resets the previews.
- **S4d:** rejection reasons are written without blanks (`coarse_reject:<reason>`), and read
  back the same.
- Tests: ptcl_sieve tester (ids past a skipped folder, counters after a restore, the reason
  round trip, the crash windows simulated by leaving the second write out).

### E. Pool outputs (O2, O3, O4, O5a, O5b, O5d, P1). Medium

- **O2:** before the first complete iteration, p06's final project is the raw project of every
  imported particle, as `pool2D_policy.md` §10 says (`terminate_stream2D` gets the pool's imported
  sets).
- **P1:** that raw project takes its groups from the newest optics map only: the copy from the
  root stub (`copy_micrographs_optics`) goes, and so does p01's call that does nothing.
- **O3:** a snapshot of the current iteration registers `frcs.bin` only when it exists; otherwise
  it is reported unservable, as for an iteration the history no longer keeps.
- **O4:** a publication carries the newest optics map's groups and table, and its `ptcl3D` is
  prepared as the final project's (2D clustering removed, shifts kept, `updatecnt` and `sampled`
  removed).
- **O5a:** a snapshot request before the pool starts is answered at once: no particles, no file.
- **O5b** (decision 7): a `mskdiam` given on the command line is kept; the sieve's applies only
  without one.
- **O5d:** a stop mid-iteration writes iteration N-1's class averages with its kept FRCs
  (`frcs_iterNNN.bin`), not the live `frcs.bin`.
- Tests: pool stage tester (raw final project, early snapshot answer, publication's optics and
  `ptcl3D`), refine2D utils (snapshot without `frcs.bin`).

### F. Multistate 3D (T1, T2b, T2c, T3, T4a, T4b). Small to medium

- **T1** (decision 8): the commander checks `2 <= nstates <= 20` up front.
- **T2b:** `particles_imported` counts the selected rows (the total stays in the log).
- **T2c:** a restart of a stage clears its list entries in the master's store (micrographs,
  class averages, volumes), so the GUI shows only the new process's.
- **T3:** `finish_run` takes the job project's data segments and keeps the stage's `projinfo`,
  `compenv` and optics table.
- **T4a** (decision 9): a malformed publication is skipped with a warning and recorded as
  rejected, so a restart skips it too.
- **T4b:** the mask diameter is read from every publication taken; a change applies from the
  next run, and is logged.
- Tests: p07 tester (state checks, counts, adoption, malformed publication), master tester (store
  cleared for a restarted stage).

### G. p03's state choice (review M6: G1-G5). Medium, needs datasets

- **G3** (decision 14): `vol_shape_descr` computes both Otsu thresholds from the voxels inside
  the mask only, and keeps the double threshold.
- **G1** (decision 11): `vol_shape_descr` also returns the dominant fraction, the largest
  component's share of the foreground voxels inside the mask (0 for an empty binarisation). The
  veto passes a state whose dominant fraction is at least `STATE_DOMINANT_FRAC` (a named
  constant, 0.8 to start) and whose population is at least `STATE_POP_FLOOR`. The component count
  stays in the log; `STATE_CC_MIN_FRAC` goes once nothing reads the count.
- **G4** (decision 12): `finish_solve3D` reads each candidate's FSC 0.143 resolution from its
  FSC file in the final folder (`get_resolution`, as `solve3D_state_fsc_lowpass` does). A state
  passes only when its resolution is within `STATE_RES_FACTOR` (a named constant, 1.5 to start)
  of the best candidate's. A state without an FSC file is not vetoed on resolution; the log says
  so. Among the states that pass both vetoes: distinct projection directions, then population,
  then the lowest state, as now.
- **G2** (decision 13): when no state passes, the fallback ranks on directions, then population;
  the component key goes. p03 logs a warning and sends it to the GUI with its status.
- `choose_state` stays pure and stateless (bound `nopass`), and takes the dominant fractions and
  resolutions in place of the component counts. p03 logs one line per state with every key, so
  a validation run can be read from the log.
- **G5:** the three constants (`STATE_DOMINANT_FRAC`, `STATE_RES_FACTOR`, `STATE_POP_FLOOR`) are
  set by a validation run on datasets with a known answer, before release. The datasets are the
  user's to name; until then the constants are marked provisional in the code and in
  `reference_generation_policy.md`.
- To confirm while implementing: that the FSC files are in the folder `find_final_solve3D_cavgs_dir`
  names, for both the first run and a restart (`mkdir=yes`).
- Tests: `simple_volume_shape_tester` (the dominant fraction of one blob, two equal blobs, a blob
  and a speck; thresholds unaffected by solvent outside the mask); `test_choose_state` rewritten
  for the new keys (a split particle with a dominant component beats a compact low-resolution
  state; the resolution veto is relative; a missing FSC does not veto; the fallback ranks on
  directions). The first case of today's test, which asserts the junk-like choice, goes.
- `reference_generation_policy.md` states the rule, its constants and why.

### H. p02's optics grouping (review D16, H1, H2). Medium

- **Algorithm** (decisions 15, 19): every pass still regroups every micrograph, but by single
  linkage without a matrix. Within a tilt group, two micrographs are linked when their shifts lie
  within `tilt_thres`, and the groups are the connected sets (union-find). Shifts are binned into
  grid cells `tilt_thres/sqrt(2)` wide, so a cell's members are linked without comparison and a
  shift is compared only with the cells up to two away; a pair of cells already in one set is
  skipped, and the first link between two cells ends their comparison. Equal shifts (no
  `dir_meta`) fall into one cell: O(N). Memory is O(N). The remaining worst case, two dense
  neighbouring cells that never link, is n_a × n_b comparisons, still without a matrix. The
  result does not depend on import order.
- **Centroid** (decision 16): a group's `opcx`/`opcy` is the mean shift of its members, as now.
- **Ids** (decision 20): after each regroup, groups are matched to the previous assignment. Pairs
  (new group, previous id) are taken in order of shared members, largest first; each new group
  and each previous id is used once. A group without a match gets the next id never used.
  Ids are therefore stable but can have gaps once groups merge.
- **Restart** (decisions 17, 20): p02 reads the newest map (`import_latest_optics_map` on its own
  folder) as the previous assignment, and publishes no map until it has re-imported every
  micrograph that map lists. During the catch-up, the readers keep applying the newest map, which
  is complete (H2).
- **Scope** (decision 18): stream only. `assign_optics_groups` in `simple_optics_groups` gets the
  new clustering and takes the previous assignment; `h_clust` stays for the batch paths
  (`simple_starproject_utils`, `simple_relion`).
- **Gaps in ids:** the consumers that treat `ogid` as 1..n are fixed: the GUI optics plot in
  `simple_sp_project_io.f90`:598-603 loops over row numbers. Others look ids up by value or offset
  them by the largest (`simple_stream_meta_plots`, `simple_projfile_utils`,
  `sp_project%append_project`). Check the rest with `rg ogid` while implementing.
- Tests: `simple_optics_groups_tester` (the existing two clusters, with and without beam tilt;
  50,000 equal shifts give one group quickly; a chain within `tilt_thres` is one group; ids kept
  when a micrograph is added to a group, the larger group keeps its id when two merge, a new
  group takes the next unused id); `simple_stream_stage_optics_tester` (a restart publishes no
  map until the newest map's micrographs are back; ids continue across a restart).
- `reference_generation_policy.md`, or a short optics section in `restart_policy.md`, states the
  grouping rule, id stability and the restart rule.

### I. p07's addon verdict and final run (review M8: I1, I2). Medium

- **No rebase during a session** (decision 21). Frozen particles keep their poses until the final
  run; the ingestion policy records this as decided, no longer as a gap.
- **A REGRESSED verdict rolls the run back** (decision 22). `finish_run` reads the addon report
  before it adopts anything. A run with any REGRESSED state (a run's states are one result) is
  not adopted: the frozen project, `frozen_active`, the stage's project and the GUI's volumes stay
  the previous base's. The next attempt waits until the cohort has grown by the cadence step
  again, max(`MIN_PTCLS_PER_STATE` × nstates, `ADDON_COHORT_FRAC` of the frozen particles) beyond
  the rolled-back cohort, so a retry is not started for a handful of new particles. A failed run
  keeps today's rule (any larger cohort). The log and the GUI status say "rolled back" with the
  verdict and the count of rollbacks in a row. A run that keeps regressing never replaces the
  base; the final run realigns everything at the end.
- **The final run** (decisions 23, 24). The pool marks its final publication: the one after
  iteration `FINAL_ITER` while the sieve's final set is in the pool (`l_sieve_final`), with
  `pool_final=yes` in the publication's out segment, as the sieve marks `sieve_final`. Once p07 has
  merged a final publication and no job runs, it starts one multistate `refine3D` of all selected
  particles, started from the base's state volumes at the base's resolution, in its own folder
  (`refine3D_final`). Its result becomes the base like any other run's (`finish_run`, with the
  project handling of T3), and its per-state resolution is logged beside the base's. It runs once
  per final publication: not again for the same selected particles (after a restart, or a repeated
  final publication). If the sieve takes finality back (a later set, `pool2D_policy.md`), the pool
  publishes again unmarked, p07 goes on with addon runs, and the next final publication starts a
  new final run. A stop request cancels it like any running job (`finalize`). The GUI status says
  "running final refine3D".
- To confirm while implementing: the `refine3D` command line for several states from given
  volumes (one per state), its low-pass start from the base's FSC, and its project and volume
  outputs for `finish_run` and `send_volumes`.
- Tests: p07 tester (`next_job` waits for the cadence step after a rollback; a REGRESSED report
  leaves the base and the stage project unchanged; an IMPROVED one is adopted; a final
  publication starts the final run once; a non-final one does not); pool stage tester (the
  publication after `FINAL_ITER` with the sieve's final set carries `pool_final`, an earlier one
  does not).
- `stream_3D_ingestion_policy.md` (items 8 and 9, known gaps), `pool2D_policy.md` (the final
  publication's flag) and a pointer from `solve3D_addon_policy.md` §11 state the rules.

### J. Compile-time policy in the stream (G12). Medium, measured

- Every plain component listed in section 1d becomes `allocatable`, allocated in `new` (or before
  the constructor that fills it) and released in `kill`: the `cmdline`s, the stages'
  `sp_project`s and `qsys_env`s, `ptcl_sieve%qenv`, `chunk2D_state%cline`.
- `stream_master_stage` keeps one command line: the fork's, which the stage record reads through
  the fork, instead of a second copy made on every start.
- Testers declare stages, the sieve, `parameters` and `cmdline` as `class(T), allocatable` (the
  master tester, p03's tester, the sieve tester).
- p00's seven commander imports move into `stage_commander`.
- No behaviour change; the existing suites are the test. A compile-time gain is claimed only from
  a like-for-like `scripts/profile_build.sh` comparison before and after, which the user runs.

### K. The NICE contract (G14). Small

- `ipc_policy.md` §5 and §6 are the contract: NICE answers every heartbeat with its whole state;
  the master forwards to each stage only what changed since its last update, nothing to a
  stopping stage, and acts on a `restart_<key>` once, again only after the key has left an answer.
  The section says so in those words, lists every key per stage, and is marked as the contract.
- NICE stops sending `ref_selection` and `increase_nmics` (`streamjob.py`), with the stream view
  that offers the sieve-reference selection (`stream_views.py`:1028) checked with NICE's owners;
  `ipc_policy.md` §9 drops the two entries.
- The NICE side points to `ipc_policy.md` (a short note beside `streamjob.py`'s update code).
- Tests: the master's GUI-command tests already cover the keys; NICE's own tests if it has them.

### L. One job type (G15). Large

- **One stream job type** (decision 28): `qsys_async_job`, extended, runs every queued job of the
  stream: submit with its record, poll its exit code, cancel, a fresh folder on restart, and the
  caller's retry policy. p01's and p04's job sets (`simple_stream_job_sets`), the sieve's chunks
  and the pool iteration move onto it; their own polling for a job's end goes. Domain markers
  (`SOLVE2D_FINISHED`, `REJECTION_FINISHED`, `COMPLETE`, `REFINE2D_FINISHED` and the like) stay
  where they mean a result rather than a job's end. `solve2D_chunks` stays as it is (fix plan
  decision 10).
- **Liveness by asking the scheduler** (decision 29, replacing fix plan decision 26). While a job
  has no exit code, the job type asks every `JOB_LIVENESS_S` (a named constant, a few minutes)
  whether its recorded id still exists: `squeue`, `bjobs` or `qstat` for a scheduler job, its
  pid for a local one on the recorded host. A job seen gone twice in a row with no exit code (one
  miss can be a scheduler's completing state) has failed and takes the caller's failure path. A
  local job recorded on another host than the stage's is not checked, and the log says so.
- **A local job in its own process group** (decision 30). The job script starts under `setsid`;
  the record holds its process group, and the cancel is `kill -TERM -<pgid>`, which reaches the
  script, its children and theirs. Scheduler jobs are unchanged.
- To confirm while implementing: persistent-worker tasks (what their record holds, and whether a
  lost worker leaves its task's script running or gone); `setsid` on the platforms SIMPLE runs on.
- Tests: `unit_parallel` "qsys control" (a local job's grandchild ends with the cancel; a job
  removed from a fake scheduler without an exit code is failed after two checks); per stage, a job
  that vanishes takes the stage's failure path (p01 and p04 sets, a sieve chunk, the pool
  iteration, p03's and p07's jobs).
- `restart_policy.md` (the job lifecycle, its §7 gaps on liveness and grandchildren removed), the
  stage policies and `ptcl_sieve_policy.md` say so.

### G13, deferred (decision 26)

Fork+exec, or `posix_spawn`, of `simple_stream prg=<stage>` with the pipe descriptors on its
command line is not in this plan. The child keeps resetting what it inherits (the server, WS0;
the flush and SIGINT, workstream A), and `restart_policy.md` §7 keeps the fork from a process
with running threads as a known gap.

## 3. Decisions

Taken 5 October 2026 (1-10 for D41, 11-14 for M6, 15-20 for D16, 21-24 for M8, 25-30 for
G12 to G15).

1. **Ctrl-C (M1).** Each stage ignores SIGINT after the fork. The master catches Ctrl-C and
   stops the stages in order (request, wait, kill). There is no separate process group.
2. **The jobs of a killed stage (M4d).** After `force_stop` kills a stage, the master reads the
   stage's job records (`simple_qsys_job_record`) from its folder and cancels each job, local or
   cluster.
3. **The movies short of a set (P3c).** p01 submits the leftover movies as a partial set once no
   new movie has arrived for a quiet period (a few watcher polls, as a named constant). This covers
   sessions that end without a stop.
4. **No accepted diameter bin (R1).** p03 logs a warning, skips extraction for the pass and picks
   again once more micrographs have arrived. It never extracts with box 0.
5. **A forced `box_extract` below 128 (R2).** The user's value is used, with a warning that it is
   below 128.
6. **The template low-pass (R3).** The clamp is fixed to `max(15, min(30, x))` Å, as intended.
   This changes picking: the templates' low-pass now depends on the diameter.
7. **A given mask diameter against the sieve's (O5b).** A `mskdiam` given on the command line is
   kept. The sieve's diameter replaces the mask at iteration 10 only when none was given.
   `pool2D_policy.md` says so.
8. **p07's state count (T1).** The p07 commander checks `2 <= nstates <= 20`
   (`MAX_STATES_SOLVE3D_MULTISTATE`) before `params%new` and stops with a clear error otherwise.
9. **A malformed publication (T4a).** p07 logs a warning, records the publication as rejected in
   its folder so that a restart skips it too, and waits for the pool's next one.
10. **The sieve's crash windows (S1).** The writes are reordered so the completion marker comes
    last (fine project before the coarse chunks' completion; coarse project and
    `chunked_mics.txt` before the chunk counts as complete). On restart, `chunked_mics.txt` is
    also rebuilt from the complete chunks' projects, so a crash between the writes can neither
    duplicate nor drop micrographs.
11. **p03's shape veto (G1).** A dominant component: a state passes when its largest component
    holds at least a set fraction (0.8 to start) of the foreground inside the mask. It replaces
    "exactly one component", so specks and a minor split no longer fail a state.
12. **Evidence of structure (G4).** A relative FSC veto: a state passes only when its FSC 0.143
    resolution is within a set factor (1.5 to start) of the best candidate's. View coverage, then
    population, rank the states that pass.
13. **The fallback when no state passes (G2).** Directions, then population; the component key
    goes. A warning goes to the log and the GUI.
14. **The binarisation (G3).** Both Otsu thresholds come from the voxels inside the mask; the
    double threshold stays.
15. **p02's grouping (D16).** A full regroup on every pass, as now, with a matrix-free clustering.
    Incremental assignment, which keeps groups as they were formed, was not chosen: the groups stay
    global and independent of import order.
16. **A group's centroid (D16).** The running mean of its members.
17. **p02's restart (D16, H2).** The newest optics map gives the previous assignment.
18. **Scope (D16).** Stream only; `h_clust` stays for the batch paths.
19. **The clustering (D16).** Single linkage on a grid: shifts within `tilt_thres` are linked,
    and groups are the connected sets. Closely spaced shift positions can chain into one group;
    that is the rule `h_clust` approximates.
20. **Ids and the restart's maps (D16, H1, H2).** A regrouped group keeps the id of the previous
    group it shares most members with; a new group takes the next unused id. After a restart, p02
    publishes no map until every micrograph of the newest map is re-imported.
21. **Rebasing during a session (M8, I2).** None. Frozen particles keep their poses until the
    final run.
22. **A REGRESSED addon (M8, I1).** Rolled back: the previous base stays, and the next attempt
    waits for a larger cohort (the cadence step beyond the rolled-back one). No rebase follows.
23. **A final run (M8, I2).** On the pool's final publication, p07 realigns every selected
    particle once.
24. **The final run's kind (M8).** Multistate `refine3D` from the base's state volumes, at its
    resolution; states keep their identity.
25. **The compile-time policy (G12).** Applied in full to the stream code (heavy components
    allocatable, testers' objects polymorphic allocatables, p00's imports local), with a
    `profile_build.sh` comparison before and after.
26. **Fork without exec (G13).** Not yet: no fork+exec in this plan.
27. **The protocol with NICE (G14).** State, as now, written down as the contract; NICE stops
    sending the keys SIMPLE ignores.
28. **Queued jobs (G15).** One job type for every queued job of the stream.
29. **Liveness (G15).** The job type asks the scheduler (or the local host) whether a job without
    an exit code still exists; gone twice in a row means failed. This replaces fix plan decision
    26.
30. **A local job's cancel (G15).** Local jobs run in their own process group, which the cancel
    signals as a whole.

## 4. Order

A, then D (both guard data and processes), then H (D16 caps a session's length), then E and
F, then I (after F: both change p07's `finish_run`), then B and C, then G. K can land at any
time. L follows the functional workstreams, since it rebuilds the job paths several of them
touch; J comes last, once the code it reshapes has settled, with its build comparison. G changes which
references p04 picks with, so its constants need the validation run before its default ships.
Each workstream updates the policy it touches (`ipc_policy.md`, `restart_policy.md`,
`ptcl_sieve_policy.md`, `pool2D_policy.md`, `stream_3D_ingestion_policy.md`,
`reference_generation_policy.md`).
