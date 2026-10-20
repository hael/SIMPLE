# Offline stream emulation: emulate_solve3D_stream and solve3D_chunks

Status (2026-10-14): `emulate_solve3D_stream` is removed from the code base
(its commander, UI entry, the parameters `nptcls_addon` and `rollback`, the
explicit split `partition_selected` and its report). `solve3D_chunks`
remains; this note stays as its design record. The implementation is
done. The open work is four validation runs, listed in
[emulate_solve3D_stream_unfinished.md](../planned/emulate_solve3D_stream_unfinished.md).
Its leaf module
`simple_solve3D_stream_emulation` is now `simple_solve3D_chunks_plan`, with
`chunks_step`, `chunks_report`, `prepare_chunks_source` and
`select_chunks_step` for the `emulation_*` names below.

Implementation note. `emulate_solve3D_stream` was committed on 2026-10-10
(`77fb358a0`, `9be7dd841`), `solve3D_chunks` on 2026-10-11 (`da5a17e89`). The
whole-stack run cuts, the `solve3D_chunks` defaults and UI, and
`SOLVE3D_BENCH.txt` followed in `76c37002b`. Validation is still open
(section 9); it moved to the unfinished note above.

Policies: `doc/policies/3D/solve3D_chunks_policy.md` is the contract of
`solve3D_chunks`. `doc/policies/3D/solve3D_policy.md` and
`doc/policies/3D/solve3D_addon_policy.md` are the contracts of the runs both
workflows execute.

## 1. Purpose

Stream task 7 (`simple_stream_stage_solve3D`) runs one `solve3D` on the first
particles the pool publishes. It then runs one `solve3D_addon` for each growth
of the pool. Its run time and map quality depend on a few settings: base size,
cohort size, `nstates`, `nstages`, `lpstop`, the sample size, `nparts` and
`nthr`. A live session cannot be replayed, and its timing is entangled with
the upstream tasks. So these settings cannot be tuned there.

Two `simple_exec` programs replay the stream's `solve3D` and `solve3D_addon`
chain offline, on an existing project:

- **`emulate_solve3D_stream`** is the explicit control. The caller sets the
  base and add-on sizes and chooses whether a regressed add-on is rolled
  back.

  ```text
  simple_exec prg=emulate_solve3D_stream projfile=<project> nptcls_base=<n> nptcls_addon=<m> <solve3D inputs>
  ```

- **`solve3D_chunks`** sizes the runs itself. It is meant mainly for stream
  emulation and for finding stream parameters. Its runs are cut as the stream
  cuts them: consecutive whole stacks in arrival order.

  ```text
  simple_exec prg=solve3D_chunks projfile=<project> <solve3D inputs>
  ```

Each invocation is one setting; a parameter sweep is several invocations.
Both produce the chain of projects the stream would publish, and a report that
puts the steps' timing and quality side by side.

## 2. The stream and the two emulations

Sources: `simple_stream_stage_solve3D.f90` (`start_solve3D`, `start_addon`,
`cap_rows`, `roll_back_addon`) and
`simple_commanders_stream_p07_solve3D_multistate.f90`.

| | Stream | `emulate_solve3D_stream` | `solve3D_chunks` |
| --- | --- | --- | --- |
| Rows | one row set that grows; rows not in a run are deselected in the job's project | all rows from the start; rows of later runs deselected | same as the emulation |
| Base | the first published particles, whole stacks up to `nptcls3D_max` | the first `nptcls_base` selected particles | the first block of whole stacks, sized automatically |
| Cohort | whole stacks in order, up to `nptcls_addon_max` (`cap_rows`); cadence gate | `nptcls_addon` particles each; no cadence gate | the next block of whole stacks; no cadence gate |
| Add-on | `solve3D_addon`, compute keys only, the rest replayed from the manifest | same commander | same commander |
| Regressed add-on | rolled back, cohort retried with more particles | `rollback=yes` (default): previous result stays frozen, the cohort joins the next | always adopted, logged |
| Sample size | from the stream's settings | caller's `nsample` | fixed total of 5000 per iteration |

Neither workflow reproduces the upstream tasks: the pool's 2D, its selection
model, or the sieve's `cavg_ini_ext` start. Nor does either run the final
multi-state `refine3D` the stream runs on `pool_final`.

## 3. Run plans

Both workflows count the rows whose input `ptcl2D` state is above 0, in row
order. That is the stream's order: pool stacks in import order. Inactive rows
belong to no run and stay inactive in every step. Let `N` be the number of
active rows.

### 3.1 emulate_solve3D_stream

```text
base       = active particles 1 .. nptcls_base
rem        = N - nptcls_base
nchunks    = max(1, rem / nptcls_addon)          (integer division)
chunk k    = nptcls_addon particles, k = 1 .. nchunks-1
chunk last = rem - (nchunks-1)*nptcls_addon
```

Example: `N=10000, nptcls_base=3000, nptcls_addon=2000` gives
3000 | 2000 | 2000 | 3000. A remainder below one chunk becomes a single
add-on; the base is never enlarged, since its size is the control variable.
The cuts fall at exact particle counts, not at stack boundaries.

Refused before any run: `nptcls_base < 5*nstates`, `nptcls_addon <
5*nstates`, `rem < 5*nstates`, and `nptcls_base >= N`.

### 3.2 solve3D_chunks

The full rule is in the policy (`solve3D_chunks_policy.md`, sections 3 and 4).
In short:

- **Sample size.** Every run samples 5000 particles per iteration, all states
  together. `solve3D` samples `nsample * nstates`, a rule that stays
  unchanged, so the workflow sets `nsample = 5000 / nstates`, rounded down.
  The add-ons replay it from the base manifest.
- **Run sizes.** Every run searches 10000 to 30000 new particles, two to six
  times the sample size. The run count is `R = ceiling(N / 30000)`.
- **Cuts.** The active particles are cut in row order into `R` blocks of
  consecutive whole stacks. Each cut goes to the stack boundary nearest its
  balanced target that keeps every run within the limits. Failing that, it
  falls inside a stack; the log and the report count those cuts.
- **Refusal.** Fewer than 10000 active particles is the only one. With 10000
  to 30000, only the base runs.

Run 0 is the base, the start of the session as in the stream. Uneven views
within a run are left to `balance=cavg`.

## 4. Steps

### 4.1 Directory layout

The workflow's `params%new` creates the execution directory `E` and copies the
input project into it.

```text
E/<project>.simple                  input copy; solve3D_chunks publishes the final result over it
E/source.simple                     the 3D-free source (4.2)
E/base.simple                       base step project
E/1_solve3D/base.simple             base result, registers its run manifest
E/addon_01.simple                   add-on 1 step project, replaced by its published union
E/2_solve3D_addon/                  add-on 1 run directory: frozen copy, maps, manifest, add-on report
E/addon_02.simple, E/3_solve3D_addon/, ...
E/emulate_solve3D_stream_report.txt or E/solve3D_chunks_report.txt
```

Every run directory also holds `SOLVE3D_BENCH.txt` (section 6).

### 4.2 The source copy

Once, on `source.simple`, everything 3D and every run registration is
removed:

- `cls3D`;
- the 3D entries of `out`;
- the `ptcl3D` alignment;
- the `projinfo` entries `sigma2_state`, `solve3D_manifest` and
  `solve3D_run_id`.

The 2D solution stays: `cls2D`, and the class-average and `frc2D` entries of
`out`. That keeps `balance=cavg` and the FRC-planned ladder available. These
are the classes of the whole dataset, a mild look-ahead the stream does not
have. Without a usable 2D solution the base runs under `balance=none`, with
the `mskdiam` ladder.

### 4.3 Base and add-on k

The step project of run `k` is the source with runs `k+1` and later
deselected in `ptcl2D` and `ptcl3D`. It is rebuilt from the source each time.
The frozen particles' poses and labels come from the frozen project, not the
step project.

- **Base:** `solve3D projfile=E/base.simple mkdir=yes`, with the base command
  line (section 5). Its result is the first frozen project.
- **Add-on k:** `solve3D_addon projfile=E/addon_kk.simple
  projfile_frozen=<current frozen> mkdir=yes`, with the forwarded keys. It
  publishes the union over `E/addon_kk.simple`, which becomes the frozen
  project of add-on `k+1`. Each add-on runs a registration pass and the base
  run's last stage.

Both commanders run in-process, unchanged. The workflow returns to `E` after
each run.

### 4.4 Verdicts and adoption

After every add-on, both workflows read its `solve3D_addon_report.txt` and
keep the per-state verdicts.

- **`emulate_solve3D_stream`:** with `rollback=yes` (the default, as in the
  stream), a REGRESSED state leaves the frozen project unchanged. The rejected
  cohort stays selected for the next add-on, the closest deterministic
  analogue of the stream's retry with a grown cohort. `rollback=no` adopts
  every add-on.
- **`solve3D_chunks`:** every add-on is adopted, so each block is searched
  exactly once. A regression is logged with the states it affects and shown in
  the report.

### 4.5 Publication

At the end, `solve3D_chunks` publishes the final result over
`E/<project>.simple`, where `simple_exec` records the job. It does this the
way `solve3D_addon` publishes its own result: temporary file renamed over the
target, manifest registered by absolute path. The emulation does not publish.

## 5. Command line and routing

The UI entries share `add_solve3D_base_inputs` (`simple_ui_solve3D.f90`).
Both workflows ask for the stream-compatible inputs, without `vol1`,
`cavg_ini`, `cavg_ini_ext` and `state`. `solve3D_chunks` also drops `nsample`,
`rec_backend` and the PCG controls.

### 5.1 emulate_solve3D_stream

| Key | Status |
| --- | --- |
| `nptcls_base`, `nptcls_addon` | required |
| `rollback` | optional, default yes |
| `addon_diag` | optional, default no; goes to every add-on |
| other stream-compatible `solve3D` inputs, including `nsample` and `rec_backend` | as in `solve3D`; go to the base |
| `vol1`, `cavg_ini`, `cavg_ini_ext`, `state`, `projfile_frozen` | refused |

### 5.2 solve3D_chunks

| Key | Status |
| --- | --- |
| `projfile` | required |
| `addon_diag` | optional, default yes; goes to every add-on |
| `filt_mode` | optional, default `nonuniform_lpset` |
| other stream-compatible `solve3D` inputs | as in `solve3D`; go to the base; `nstages` at least 3 |
| `nsample` | not accepted; set to `5000 / nstates` |
| `rec_backend` and the PCG controls | not accepted; `rec_backend=gridding` is set |
| `nptcls_base`, `nptcls_addon`, `rollback` | not accepted |
| `vol1`, `cavg_ini`, `cavg_ini_ext`, `state`, `projfile_frozen` | not accepted |

Every key is checked against the program's UI entry, so nothing is ignored
silently. The defaults are set before `params%new`, in the form the UI
default audit (`scripts/default_audit.py`) reads, so the UI shows them.
`balance` keeps `solve3D`'s default: `cavg` with a usable 2D solution, `none`
without.

### 5.3 Routing to the add-ons

An entry key goes to an add-on when the `solve3D_addon` UI accepts it and the
base manifest does not record it as an input (`manifest_records_input`). That
covers `addon_diag`, `nparts`, `nthr` and `overlap`. Recorded inputs
(`balance`, `nclust`, `mskdiam`, `filt_mode`, `rec_backend`, ...) reach the
add-ons through the base manifest. `nsample` is not a recorded input, but the
manifest's `replay` sets it, with `nstates`, from the base run's sampling
record. `solve3D_chunks` therefore leaves its own `nsample` off the add-on
command line. The workflows' own keys never reach the child runs.

## 6. Reports and benchmarks

The report is rewritten after every run, so a run that stops part-way keeps
what it measured. It has one row per run:

- kind, run index, frozen and added particles, union size;
- last stage, with its low-pass limit and box;
- wall time, and seconds per 1000 added particles;
- adopted or rolled back.

Then come per-state lines:

- **Every run:** the FSC=0.143 and FSC=0.5 resolutions of the result.
- **Add-ons, in addition:** the verdict, FSC shell shift, map correlation
  with the frozen map, and the cohort-only resolution. That last figure is
  0.00 when `addon_diag=no`, which means "not computed".

The headers differ:

- **`emulate_solve3D_stream`** records the requested sizes and the rollback
  policy.
- **`solve3D_chunks`** records the run count, `nstates`, `nsample` with the
  total it samples, every run's size, and the number of cuts inside a stack.

Every `solve3D` and `solve3D_addon` run also writes `SOLVE3D_BENCH.txt` in its
run directory, under `L_BENCH_GLOB`, in the format of the other SIMPLE
benchmark files. It times six phases:

- setup;
- the starting state (for an add-on, the frozen accumulation);
- the registration pass;
- the refine3D stages;
- the final reconstruction;
- the epilogue (for an add-on, its report and cohort-only map).

Read across the chain, these separate the capped search from the passes over
the whole union.

## 7. Code map

| File | Owns |
| --- | --- |
| `src/main/commanders/simple/simple_commanders_solve3D_chunks.f90` | `commander_emulate_solve3D_stream` and `commander_solve3D_chunks`, each a complete executor: validation, the run plan, the child runs, routing, statistics, the report, adoption, and for `solve3D_chunks` publication; shared leaf helpers for step projects, key forwarding and statistics |
| `src/main/solve/simple_solve3D_stream_emulation.f90` | `partition_selected` (explicit), `chunk_nsample`, `derive_chunk_run_count`, `partition_selected_balanced` (whole-stack cuts), the constants, source preparation, step selection, the report |
| `src/main/solve/simple_solve3D_stream_emulation_tester.f90` | the unit sub-suite `solve3D stream emulation` (`unit_project`) |
| `src/main/commanders/simple/simple_commanders_solve3D.f90` | the child runs, unchanged except `SOLVE3D_BENCH.txt` |
| `src/main/ui/simple/simple_ui_solve3D.f90`, `src/main/exec/simple_exec_solve3D.f90` | the UI entries and dispatch |
| `simple_parameters.f90`, `simple_parameters_parse.f90` | `nptcls_addon` and `rollback`; `solve3D_chunks` adds no keys |

## 8. Evidence so far

All three runs below used the version committed in `da5a17e89`, which still
dealt particles to the runs class by class. It did not yet have the
whole-stack cuts.

| Project | Particles | Runs | States | Storage | Result |
| --- | --- | --- | --- | --- | --- |
| msp1 | 141838 | 5 of about 28368 | 2 | `/usr/local/cache` | base 2238 s; add-ons 1461, 1510, 1560 s; resolution improves at every add-on, to 5.6 and 5.9 Å after add-on 3; cohort-only maps 4.6 to 6.6 Å |
| embb | 57970 | 2 of 28985 | 3 | `/usr/local/data` | base 2320 s, add-on 1179 s; states 1 and 3 improve to 6.0 and 6.2 Å; state 2 stays at 19.7 Å (collapsed); run without `addon_diag` |
| ribosome | 436455 of 751729 rows | 15 of 29097 | 3 | spinning disk | base at iteration 65 of at most 84 after almost two hours; a late iteration took 204 s, 140 s of it the probabilistic pre-step, with the CPUs often idle; one state collapsed to 1979 particles at 19.7 Å |

What the runs show:

- **The add-on chain works.** It improves the maps at every step, and with
  growing union size each add-on takes slightly longer: about 50 s per step in
  msp1.
- **Reading dominated the ribosome run.** Dealt by class, each iteration's
  sample touched about one particle per stack file across the whole dataset.
  That motivated the whole-stack cuts.
- **One state collapsed in both three-state runs.** Each add-on starts its
  cohort with states assigned uniformly at random, so a collapsed state still
  receives a third of every cohort.

## 9. Open validation

- The ribosome rerun with whole-stack cuts, with `SOLVE3D_BENCH.txt` from
  every run directory.
- A base-only run, a split near the 30000 boundary, and a cut that falls
  inside a stack.
- An `emulate_solve3D_stream` rollback case.
- The published project: its solution and its manifest registration.

These four checks are now tracked in
[emulate_solve3D_stream_unfinished.md](../planned/emulate_solve3D_stream_unfinished.md).
Note that `emulate_solve3D_stream` has since been removed, so its rollback
case (third item) no longer applies.

## 10. Known limits

- **Row growth is not exercised.** All rows exist from the start, so the
  superset relation never sees appended rows. The stream's append path is
  covered by the snapshot generator test.
- **The explicit emulation cuts at exact counts,** not at stack boundaries.
- **Look-ahead in the 2D solution** (4.2).
- **The cap bounds the search, not the union.** Each add-on also
  accumulates the frozen term over every frozen particle and reconstructs the
  whole union. That cost grows with every step, in the stream as well.
- **Fewer sampled particles per state.** With a fixed total, each state gets
  about `5000 / nstates` freshly searched particles per iteration.

## 11. Decisions

For `emulate_solve3D_stream` (reviewed 2026-10-10):

1. Masked steps on one row set, not physical snapshots.
2. Counts are of selected particles, not rows.
3. A remainder below one chunk becomes a single add-on, never merged into the
   base.
4. `rollback=yes` by default; the rejected cohort joins the next.
5. `balance` follows `solve3D`'s default.
6. The source's 2D classes and FRCs are kept, despite their look-ahead.
7. Add-ons get only forwarded keys.
8. The runs execute the commanders in-process; no subprocesses.
9. Later, optional: a batch reference `solve3D` on all particles, with the
   final union compared against it.

For `solve3D_chunks` (Hans, 2026-10-11):

10. All sizing derives from a fixed sample of 5000 per iteration. Runs hold
    10000 to 30000 particles. `nsample = 5000 / nstates`, while `solve3D`'s own
    rule (`nsample * nstates` per iteration) stays unchanged.
11. It is meant mainly for stream emulation. Runs are the fewest that fit,
    cut in row order at stack boundaries, as the stream cuts them. A
    class-by-class deal was tried and dropped: it hid the stream's base
    condition and made reading bound by disk seeks.
12. A base-only plan is allowed. Fewer than 10000 particles is refused, and so
    is `nstages` below 3, before anything is written.
13. Every add-on is adopted. Regressions are logged.
14. The 30000 cap is on the newly searched particles of a run, not on the
    union.
15. The final result is published over the execution directory's project.
16. Defaults: `addon_diag=yes`, `filt_mode=nonuniform_lpset`, and gridding
    reconstruction always, with no PCG controls; `balance` follows
    `solve3D`.
