# Pool 2D Policy

The 2D pool's schedule: when it runs and pauses, what each iteration samples, how its resolution
and dimensions evolve, its snapshots and its final project. What a change must preserve.

## 1. Scope

- The stage (p06): `src/main/stream/stages/simple_stream_stage_pool2D.f90`, driven by
  `src/main/commanders/stream/simple_commanders_stream_p06_pool2D.f90`.
- The pool's iterations, sampling, resolution and dimensions: `src/main/stream/pool2D/simple_stream_pool2D_utils.f90`.
- Snapshots, publications for 3D, the final project: `src/main/stream/pool2D/simple_stream_refine2D_utils.f90`.
- The pool's state: `src/main/stream/pool2D/simple_stream2D_state.f90` (module variables).

The publications for 3D are governed by `doc/policies/stream/stream_3D_ingestion_policy.md`; restarts
by `doc/policies/stream/restart_policy.md`.

## 2. Inputs

1. The particle sets the sieve (p05) hands off to its completed folder; each is a fine chunk's
   project. The set marked `sieve_final=yes` in `os_out` is the sieve's last.
2. The mask diameter of the sieve's 2D, read from the first set's class averages
   (`read_final_mskdiam`).
3. The master's settings: `ncls` 150, `nparts` 6, `nthr` 8 (p00's table), `optics_dir`.
4. GUI updates: a new 2D mask diameter (`mskdiam2D`), a snapshot request (`snapshot2D`).

## 3. A pass

Each pass of p06 runs, in this order:

1. take the sets handed off since the last pass;
2. refresh the pool: a paused pool refreshes its statistics; a running one is checked for a
   completed iteration, whose particle parameters, classes and resolution come back into the
   pool, with a dimension update (section 7);
3. publish the pool's classified state for 3D when the refresh brought back an iteration past
   `EXPORT_START_ITER` (25), before anything new is imported or dispatched (the 3D ingestion
   policy);
4. import the new sets when the pool is free (the first import starts it, section 4);
5. apply the pause rule (section 5);
6. start the next iteration unless paused or below the starting threshold;
7. from iteration `MSKDIAM_SWITCH_ITER` (10), switch once to the sieve's mask diameter;
8. send the class averages to the GUI;
9. apply the GUI's updates.

## 4. Start

1. The first import gives the pool its pixel size and box (from the data) and, when none was
   given, a mask diameter of (box/2 - soft edge - 1 px) * 2 * smpd (`default_mskdiam`).
2. The pool runs downscaled when that is possible (`setup_downscaling`): to a pixel size of up to
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
  pausing (`runs_to_final`).

## 6. An iteration

1. **History:** before dispatch, the completed iteration's pool project, with its class averages
   and a copy of its FRCs (`frcs_iterNNN.bin`), goes into the history. The history is a ring of
   `POOL_NHISTORY` (5) iterations (`simple_stream2D_state`, slot `pool_history_slot(iter)`). The
   new entry replaces the iteration `POOL_NHISTORY` before it, so memory is bounded at five copies
   and no copy is made beyond the one per iteration.
2. **Sampling:** stacks are shuffled and taken until more than `STREAM_NPTCLS_MAX` (500,000)
   selected particles, or `nsample_max` when given. There is no memory of earlier iterations.
3. **Fractional update:** beyond the cap, the iteration runs with
   `update_frac` = (selected particles - particles outside the sample) / (selected particles),
   or the user's `update_frac`. The out-of-sample class populations are recorded in `cls2D`
   (`prev_pop_even`, `prev_pop_odd`), and `center=no` is set on the pool's command line.
4. **New particles:** from iteration 2, particles never updated get a random populated class
   before their first alignment.
5. **Resolution:** `lpstart` and `lpstop` come from the mask diameter (`mskdiam2lplimits`), with
   `lpstop` at least twice the pool's pixel size, or the user's `lpstop`. Until iteration 20 the
   low-pass limit follows lp = lpstop + (lpstart - lpstop) * (20 - iter) / 20 with a Gaussian
   filter at lp. The shift search is 0 until iteration 5, then `MINSHIFT`. After iteration 20 the
   limit is free.
6. **Tidying:** each iteration deletes the files of the iteration that has just left the history
   (`tidy_2Dstream_iter`): its class averages (and even and odd halves), JPEG, class STAR file and
   `frcs_iterNNN.bin`. History and files always hold the same iterations.

## 7. Dimensions

With `dynreslim=yes` (p06's default) and `autoscale=yes`, from iteration 10, the pool grows its
box when the resolution has sat at Nyquist for `POOL_NPREV_RES` (5) iterations:

- the box goes to the next magic size, and the pixel size follows;
- the class averages and FRCs are Fourier-padded to it, and registered again in the pool's
  project;
- the pixel size never goes below `POOL_SMPD_HARD_LIMIT`.

## 8. Classes

The pool never rejects a class: `ncls_rejected_glob` stays 0. Class selection happens once per publication
in 3D, and in snapshots by the user's selection.

## 9. Snapshots

1. A snapshot request from the GUI (id, iteration, selected classes, file name) is written once
   per id, once the pool has started.
2. It holds the selected classes of the requested iteration (`write_pool_snapshot`, the request
   passed as arguments):
   - for the current iteration, the pool project;
   - for an earlier one, its history entry, when the history still holds it and its files exist.
3. It is written in `snapshots/<name>/` with its class averages, FRCs, micrograph and particle
   STAR files and the newest optics map.
4. The STAR files' optics group ids are offset by (`nicedispid` - 1) * 500, so the GUI's displays
   do not collide.
5. The snapshot's class averages go back to the GUI, with its particle count.
6. A request the pool cannot serve is not fatal: an iteration no longer kept, or files missing.
   Nothing is written, a warning is logged, and the GUI is told the snapshot has 0 particles and
   no file.

## 10. The final project

When p06 stops, the pool's last complete iteration is written as the stage's project:
- the class averages and FRCs at the native sampling;
- the newest optics map's groups;
- `ptcl3D` prepared for 3D: 2D clustering removed (class, in-plane angle, corr, frac), shifts
  kept, `updatecnt` and `sampled` removed;
- the micrograph and particle STAR files.

Its class averages are then ranked. Before any iteration, the raw imported particles are written
instead.

## 11. Change rules

- The schedule constants are named in the stage (`EARLY_RATE_FACTOR`, `LATE_RATE_FACTOR`,
  `LATE_ITER`, `FINAL_ITER`, `MSKDIAM_SWITCH_ITER`, `EXPORT_START_ITER`, `NPUBLICATIONS_KEPT`,
  `NPTCLS_PER_CLS_MIN`,
  `OPTICS_ID_DELTA`) and in the pool module (`ITERLIM`, `ITERSHIFT`). A change of any of them
  updates this policy.
- The pool command line sets `msk_crop` explicitly; changing `mskdiam` recomputes `msk_crop`
  (`update_mskdiam`).
- Anything written for downstream use is rescaled to the native sampling or labelled with the
  pool's.
- `POOL_NHISTORY` (`simple_stream2D_state`) sets both the history and the iterations whose files
  are kept; a change updates sections 6 and 9.
- Tests: `unit_stream` "pool 2D" covers the rules (pause, final run, default mask), the set
  transfer, the final-set flag, the publication's contents and numbering, and the snapshot's
  report to the GUI (written, and not written). Iterations, the history, the writing of
  snapshots and dimension changes have no unit test.

## 12. Known gaps

- **Module state** (review G1). The pool is about thirty-five public module variables, read and
  written by three modules; p06 holds a pointer the pool reads as `master_cline`. One process
  runs one pool, and the pool cannot be unit-tested. Proposal R6: a `stream_pool2D` type.
- **Fractional update** (M4):
  - suspected double sub-sampling (about frac²);
  - the out-of-sample populations are written and never read;
  - `center=no` stays set for good.
- **Not reproducible:** new particles' classes are drawn inside an OpenMP loop.
- **Late particles:** the low-pass ramp follows the global iteration, so particles arriving after
  iteration 20 never see the coarse limits.
- **Dimension growth** pads class averages and FRCs, which adds no information. The guard
  `ncls_glob < ncls_max` can never be true.
- **The mask diameter changes mid-run** at iteration 10.
- **Optics group ids** of the STAR files depend on the GUI display id.
- **History memory:** five full copies of the pool project in memory. A pool of millions of
  particles may want the history on disk instead (writing one project per iteration).
