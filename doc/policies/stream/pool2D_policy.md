# Pool 2D Policy

The 2D pool's schedule: when it runs and pauses, what each iteration samples, how its resolution
and dimensions evolve, its snapshots and its final project. What a change must preserve.

## 1. Scope

- The stage (p06): `src/main/stream/stages/simple_stream_stage_pool2D.f90`, driven by
  `src/main/commanders/stream/simple_commanders_stream_p06_pool2D.f90`.
- The pool: `stream_pool2D` in `src/main/stream/pool2D/simple_stream_pool2D.f90`, its state
  (private components) and its iterations, sampling, resolution, dimensions, snapshots,
  publications for 3D and final project. p06 holds one, empty from its start, hands it the sets
  it imports (`append_sets`) and starts it at the first import; what the GUI shows of it is one
  record (`stats`). Its public procedures are those p06 calls, and `init_state` for its tester.
- Stateless helpers (folder clean-up, the rows a set adds, class draws, publication building and
  naming, iteration files): `src/main/stream/pool2D/simple_stream_refine2D_utils.f90`.

The publications for 3D are governed by `doc/policies/stream/stream_3D_ingestion_policy.md`; restarts
by `doc/policies/stream/restart_policy.md`.

## 2. Inputs

1. The particle sets the sieve (p05) hands off to its completed folder; each is a fine chunk's
   project. The set marked `sieve_final=yes` in `os_out` is the sieve's last. The sieve marks
   only its empty final set, handed off once every chunk has ended, which ends the intake and
   transfers nothing; the pool also takes the flag on a set with particles. A later set with
   particles and without the flag takes the note back: the sieve's upstream had more movies
   (`ptcl_sieve_policy.md`, section 6.3).
2. The mask diameter of the sieve's 2D, read from the first set with class averages
   (`read_final_mskdiam`).
3. The master's settings: `ncls` 150, `nparts` 6, `nthr` 8 (p00's table), `optics_dir`.
4. GUI updates: a new 2D mask diameter (`mskdiam2D`), a snapshot request (`snapshot2D`).

## 3. A pass

Each pass of p06 runs, in this order:

1. take the sets handed off since the last pass, in the order the sieve handed them off: by the
   number ending their name (a chunk's id, so `chunk_fine_10` after `chunk_fine_2`), the sieve's
   final sets (`sieve_final_c<n>_f<m>`) after every other (`sort_sets`), whatever order the folder
   lists them in. A final set is therefore never taken before a set the sieve handed off earlier,
   which would otherwise read as a set after the final one and take the final note back. Until
   the first publication, each set's class averages and FRCs are kept in `sieve_cavgs/` and its
   accepted particles and micrographs counted, as it is handed off;
2. publish the first set for 3D as soon as it is due, without waiting for a pool iteration: the
   sieve's own 2D of the sets handed off, once they hold `NPTCLS_FIRST3D` (100,000) accepted
   particles, cover `NMICS_FIRST3D` (500) micrographs, or the sieve's final set has come
   (`publish_first_set`; the 3D ingestion policy,
   sections 3 and 9);
3. refresh the pool: a paused pool refreshes its statistics; a running one is checked for a
   completed iteration, whose particle parameters, classes and resolution come back into the
   pool, with a dimension update (section 7);
4. publish for 3D when the refresh brought back iteration `EXPORT_START_ITER` (25) or a later
   one: the pool's classified state, before anything new is imported or dispatched (the 3D
   ingestion policy). In a pool that has published nothing yet, iteration `FIRST_EXPORT_ITER`
   (10, the `MSKDIAM_SWITCH_ITER`) is the first set's fallback trigger. When the first set cannot
   be built from the sieve's 2D, the pool's classified state goes instead, after the first
   iteration whose publication would select `NPTCLS_FIRST3D` particles or iteration 10, whichever
   comes first (`exports_after`); at iteration 10 it carries the sieve's mask diameter, which the
   pool took when it dispatched iteration 10, and earlier the pool's starting one. A restarted
   pool with publications on disk skips the first and resumes at 25. The final
   run stops at `FINAL_ITER` (25), so a short session still publishes its last iteration (until 5
   October 2026 the publications started after iteration 25, and a session whose final set came
   before it never reached 3D);
5. import the new sets when the pool is free (the first import starts it, section 4): p06 chooses
   the sets and the pool appends them (`append_sets`: their micrographs, their stacks renumbered
   after the pool's, their particles as new rows with no 2D parameters but their shifts). With
   `stepwise=yes` (p06's default, set by its commander; a registered parameter) an import takes
   sets in order until its own particles reach the starting threshold (or `ncls` * 20 before it
   is known); the rest wait for the next import. Only the particles of the import count, so a
   pool past its threshold still takes as many sets as one threshold's worth;
6. apply the pause rule (section 5);
7. start the next iteration unless paused or below the starting threshold;
8. from iteration `MSKDIAM_SWITCH_ITER` (10), switch once to the sieve's mask diameter, unless
   one was given on the command line, which is kept;
9. send the class averages to the GUI;
10. apply the GUI's updates.

## 4. Start

1. The first import gives the pool its pixel size and box (from the data) and, when none was
   given, a mask diameter of (box/2 - soft edge - 1 px) * 2 * smpd (`default_mskdiam`).
2. The pool runs downscaled when that is possible (`set_dimensions`): to a pixel size of up to
   `MAX_SMPD` (2.67 Å), and never below a 128-pixel box (`CHUNK_MINBOXSZ`). Everything written for
   use downstream is rescaled to the native sampling or labelled with the pool's.
3. The first iteration waits for max(`ncls` * 20, rate * 500) particles, where rate is the
   particles per micrograph of the sets so far; the sieve's final set lifts the wait.

## 5. Pause rule

| Iteration | Pauses when | Waits for |
|---|---|---|
| 2 to `LATE_ITER` (20) | more than one iteration has passed since the last import | max(`ncls` * 20, rate * `EARLY_RATE_FACTOR` (50)) more particles |
| after 20 | an iteration has passed since the last import | max(`ncls` * 20, rate * `LATE_RATE_FACTOR` (500)) more particles |

- An import that brings the awaited particles resumes the pool, as does a new mask diameter from
  the GUI.
- **Final run:** from the sieve's final set until `FINAL_ITER` (25), the pool runs without
  pausing (`runs_to_final`). Past `FINAL_ITER`, with the final set in the pool, it runs until a
  publication has carried `pool_final=yes` (`owes_final_run`): a pool that paused past iteration 25
  and published before the final set came runs one more iteration, because the empty final set
  brings no particles to wake it. A set after the final one takes the owed publication back.

## 6. An iteration

1. **History:** before dispatch, the completed iteration's pool project, with its class averages
   and a copy of its FRCs (`frcs_iterNNN.bin`), goes into the history. The history is a ring of
   `POOL_NHISTORY` (5) iterations (`simple_stream_pool2D`, slot `history_slot(iter)`). The
   new entry replaces the iteration `POOL_NHISTORY` before it, so memory is bounded at five copies
   and no copy is made beyond the one per iteration.
2. **Sampling:** stacks are shuffled and taken until more than `STREAM_NPTCLS_MAX` (500,000)
   selected particles, or `nsample_max` when given. There is no memory of earlier iterations.
3. **The sample is the update set**: every particle of the sampled stacks is
   updated, and the particles outside the sample keep their parameters; the class averages come
   from the sample. Only the user's `update_frac` thins the sample again, and then the class
   averages are not centered (`center=no`); otherwise the pool's own `center` applies, so a
   return to a full update centers again. The iteration's project carries `ptcl3D` as a copy of
   `ptcl2D` with the 2D clustering removed, so that it holds the particles' halves (`eo`), which the
   refine2D master checks (`ensure_field_halves`); `eo` returns to the pool with the 2D parameters.
4. **New particles:** from iteration 2, particles never updated get a populated class before
   their first alignment, drawn in one thread from a generator seeded with the iteration
   (`draw_new_classes`), so a run is reproducible; the process's generator is left as it was.
5. **Resolution:** `lpstart` and `lpstop` come from the mask diameter (`mskdiam2lplimits`), with
   `lpstop` at least twice the pool's pixel size, or the user's `lpstop`. A new mask diameter (the
   sieve's at iteration 10, or the GUI's) recomputes `lpstart`, the centering limit and the ramp. Until iteration 20 the
   low-pass limit of the next iteration follows lp = lpstop + (lpstart - lpstop) * (20 - iter) / 20,
   `iter` being the latest iteration dispatched (the limit is set just before the next is dispatched),
   with a Gaussian filter at lp. The shift search is 0 until iteration 5, then `MINSHIFT`. After
   iteration 20 the limit is free.
6. **Tidying:** each iteration deletes the files of the iteration that has just left the history
   (`tidy_2Dstream_iter`): its class averages (and even and odd halves), JPEG, class STAR file and
   `frcs_iterNNN.bin`. History and files always hold the same iterations.
7. **Job:** the iteration is a local queued job with an exit status (`EXIT_CODE_refine2D_pool`)
   and a job record; its project as made is kept (`refine2D_input.simple`). A job that writes a
   status without `REFINE2D_FINISHED` has failed: its log is kept as
   `simple_log_refine2D_pool_failed_iterNNN_attempt<k>`, the part files are cleared, and the
   iteration is submitted again once from its project as made. A second failure stops the stage,
   which reports it and writes its final project from the last complete iteration. The iteration
   runs as the pool's `qsys_async_job` (label `refine2D_pool`), whose liveness check finds a job
   killed before it writes a status (`restart_policy.md`, job lifecycle). On stop, the running
   iteration is cancelled.

## 6a. The refine2D call

The iteration runs `prg=refine2D_distr` through `simple_private_exec`: the batch refine2D master
(`commander_refine2D`, distributed strategy), whose workers run `prg=refine2D`
(`commander_refine2D_distr_worker`).

1. **Command line** (`stream_pool2D`: `init_state`, `iterate`, `update_aln_params`, `set_dimensions`,
   `set_resolution_limits`). Fixed: `projfile=refine2D.simple`, `oritype=ptcl2D`, `objfun=euclid`,
   `ml_reg=no` and `refine=snhc_smpl` (set by p06's commander; the pool accepts `snhc` and
   `snhc_smpl`), `cls_init=rand` (dropped from the second iteration), `center`, `center_type`,
   `async=yes`, `stream2d=yes`, `mkdir=no`, `nparts`, `nthr` and the worker server's address.
   Per iteration: `ncls`, `startit` and `maxits` set to the iteration, `frcs`, `refs` (the latest
   class averages), `update_frac` only when the user gave one, `mskdiam`, the working `smpd_crop`,
   `box_crop` and `msk_crop`, `lpstart`, `lpstop`, `cenlp`, and, until iteration 20, `lp`,
   `gauref`, `gaufreq`, `extr_iter`, `extr_lim` and `trs`.
2. **`stream2d=yes`** is set here and nowhere else. It changes the batch workflow in seven places:
   - with no `lp` on the command line (the pool gives one until iteration 20), the search low-pass
     limit is the class FRCs' estimate, or `lplims2D(3)` without FRCs, and not below `lpstop`
     (`simple_matcher_smpl_and_lplims.f90`);
   - carried class sums are blended from the second iteration on whenever `update_frac` is below
     0.99 (`refine2D_blends_carryover`, same file);
   - in the matcher a sampled update set is used only when the carried sums are blended (otherwise
     `update_frac` is 1), `snhc` becomes `snhc_smpl`, and a particle that is new or not yet
     searched is aligned greedily (`simple_strategy2D_matcher.f90`);
   - the sampling power is `EXTR_POWER` (`simple_strategy2D_alloc.f90`);
   - no cropped particles are made in memory for the class sums (`crop_particles2D`,
     `simple_matcher_ptcl_batch.f90`);
   - convergence uses the limits of a fractional update, `OVERLAP_2D_FRAC` and
     `FRACSRCHSPACE_FRAC` (`simple_convergence.f90`);
   - no progress file is written (`check_refine2D_convergence`).
3. **The script.** The pool's objective is `euclid`, so the job is one queued script (`start_seq`):
   `calc_pspec` on the command line without `stream2d` and `update_frac`, then `refine2D_distr`.
   Pool membership changes between iterations, so the canonical sigma2 bootstrap is made again for
   the exact pool layout in the script that clusters. The sigma2 state is the pool project's own:
   `sigma2_state.bin` in the pool folder, registered in the project's `projinfo` at `init_state`,
   and removed with the other `*.bin` files on a restart.
4. **Project hand-over.**
   1. Before dispatch the pool writes the sampled project (the stacks of the sample with their
      `ptcl2D`, `ptcl3D` and `cls2D`) as `refine2D.simple`, and keeps a copy as
      `refine2D_input.simple`.
   2. The refine2D master owns `refine2D.simple` and `sigma2_state.bin` until it ends. It merges the
      workers' documents and writes `ptcl2D`, `cls2D`, `cls3D` and `out` into the file
      (`run_distributed_cavg_assembly`, `distr_finalize_iteration`, `distr_finalize_run`), and
      touches `REFINE2D_FINISHED` as it ends. The pool does not write the file while an iteration
      runs (`l_available` is false).
   3. When `REFINE2D_FINISHED` exists (`update_status`), the pool takes the iteration back: it
      reads `cls2D`, `stk` and `ptcl2D` of `refine2D.simple` and the iteration's statistics file,
      and transfers the 2D parameters of the sampled particles that are still selected to its own
      project (`update`).
   4. A first attempt that failed is repeated once: the pool copies `refine2D_input.simple` over
      `refine2D.simple` and submits again (`update_status`).

## 7. Dimensions

The pool's native box and pixel size are those of its first import, and its mask diameter is
the given one or the default (pool state, not parameters); its working dimensions are downscaled
from them. The mask radius is clamped to (box - `COSMSKHALFWIDTH`)/2 pixels, and the clamp is
logged, so a diameter beyond the box never reaches the workers. The pool's command line
carries the working dimensions only.

With `dynreslim=yes` (p06's default) and `autoscale=yes`, from iteration 10, the pool grows its
box when the resolution has sat at Nyquist for `POOL_NPREV_RES` (5) iterations:

- the box goes to the next magic size, and the pixel size follows, when the new pixel size is
  more than 5/4 of the native one;
- the class averages (and their even and odd halves) and FRCs are Fourier-padded to it, and
  registered again in the pool's project, as are the carried class sums. Padding adds no
  information: it only lets the next iterations reach a finer resolution with the new Nyquist;
- the hard low-pass limit follows the new Nyquist (or the user's `lpstop`);
- the pixel size never goes below `POOL_SMPD_HARD_LIMIT`.

## 8. Classes

The pool never rejects a class: `ncls_rejected_glob` stays 0. Class selection happens once per publication
in 3D, and in snapshots by the user's selection.

## 9. Snapshots

1. A snapshot request from the GUI (id, iteration, selected classes, file name) is written once
   per id. A request before the pool has started is answered at once, as not written (item 6).
2. It holds the selected classes of the requested iteration (the pool's `write_snapshot`, the request
   passed as arguments):
   - for the current iteration, the pool project;
   - for an earlier one, its history entry, when the history still holds it and its files exist.
3. It is written in `snapshots/<name>/` with its class averages, FRCs, micrograph and particle
   STAR files and the newest optics map.
4. The STAR files' optics group ids are offset by (`nicedispid` - 1) * 500, so the GUI's displays
   do not collide.
5. The snapshot's class averages go back to the GUI, with its particle count.
6. A request the pool cannot serve is not fatal: no iteration yet, an iteration no longer kept,
   or files missing (also the current iteration's `frcs.bin`, which is registered only when it
   exists). Nothing is written, a warning is logged, and the GUI is told the snapshot has 0
   particles and no file.

## 10. The final project

The first publication for 3D is the sieve's own 2D, not the pool's
(`stream_3D_ingestion_policy.md`, sections 3 and 9; decided 12 October 2026, before which only
the master's former `sieve_ini3D=yes` did so): p06 keeps the class averages and FRCs of the sets
handed off in `sieve_cavgs/` until it, and removes them after.

When p06 stops, the history of completed iterations (`POOL_NHISTORY` full copies of the pool,
kept for snapshots) is freed first, so the final write does not hold it too, and the resident
memory is logged at each step (`log_rss`). The pool's last complete iteration is then written as
the stage's project:
- the class averages and FRCs at the native sampling; when the stop came mid-iteration, the
  previous iteration's class averages go with its own FRCs (`frcs_iterNNN.bin`), not with
  `frcs.bin`, which the cancelled iteration was rewriting;
- the newest optics map's groups;
- `ptcl3D` prepared for 3D: 2D clustering removed (class, in-plane angle, corr, frac), shifts
  kept, `updatecnt` and `sampled` removed;
- the micrograph and particle STAR files.

Its class averages are then ranked. Before any complete iteration, the pool's imported
micrographs, stacks and particles are written instead, as they came (no 2D parameters but their
shifts), with the newest optics map's groups and the STAR files (the pool's `finalise`). The groups
are applied to the project as written; the root folder's optics project is not read.

Publications for 3D (`stream_3D_ingestion_policy.md`): the one after iteration `FINAL_ITER` or a
later one, while the sieve's final set is in the pool, carries `pool_final=yes` in its out
segment (`publishes_final`); multistate 3D then runs its final refine3D.

## 11. Change rules

- The schedule constants are named in the stage (`EARLY_RATE_FACTOR`, `LATE_RATE_FACTOR`,
  `LATE_ITER`, `FINAL_ITER`, `MSKDIAM_SWITCH_ITER`, `FIRST_EXPORT_ITER`, `EXPORT_START_ITER`,
  `NMICS_FIRST3D`, `NPUBLICATIONS_KEPT`, `NPTCLS_PER_CLS_MIN`), in `simple_defs_stream`
  (`NPTCLS_FIRST3D`, `OPTICS_ID_DELTA`) and in the pool module (`ITERLIM`, `ITERSHIFT`, locals of
  `update_aln_params`). A change of any of them updates this policy.
- The pool command line sets `msk_crop` explicitly; changing `mskdiam` recomputes `msk_crop`,
  clamped to the box (the pool's `set_mask`), and the low-pass ramp (`set_mskdiam`).
- Anything written for downstream use is rescaled to the native sampling or labelled with the
  pool's.
- `POOL_NHISTORY` (`simple_stream_pool2D`) sets both the history and the iterations whose files
  are kept; a change updates sections 6 and 9.
- Tests: `unit_stream` "pool 2D" covers the rules (pause, final run, default mask), the set
  transfer as the stage counts it, the final-set flag, the publication's contents and numbering,
  and the snapshot's report to the GUI (written, and not written). "pool 2D object"
  (`simple_stream_pool2D_tester`) covers the pool: an empty pool and its kill, the rows a set
  adds and the pool's counts after an import, the reproducible class draws, the mask clamp and
  the mask in `stats`, a restart after kill, an unclassified pool publishing nothing, and a
  snapshot of an iteration without FRCs. Iterations, the update set, the history, the writing of
  snapshots and dimension changes need a queue and have no unit test.

## 12. Known gaps

- **One pool per folder:** the pool is a type, but its files have fixed names (the pool folder,
  its exit status, its iteration files).
- **The user's `update_frac`** thins the sample again inside refine2D (about frac of the sample).
- **Late particles:** the low-pass ramp follows the global iteration, so particles arriving after
  iteration 20 never see the coarse limits.
- **Dimension growth** pads class averages and FRCs, which adds no information.
- **The mask diameter changes mid-run** at iteration 10, when none was given.
- **Optics group ids** of the STAR files depend on the GUI display id.
- **History memory:** five full copies of the pool project in memory. A pool of millions of
  particles may want the history on disk instead (writing one project per iteration).
